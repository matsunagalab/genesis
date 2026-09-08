"""End-to-end checks on the chignolin trajectory (the demo notebook workflow).

Requires the chignolin data: ``python -m genepie.tests.download_test_data``.
"""
import os
import tempfile

import numpy as np
import pytest

from .. import genesis_exe
from ..exceptions import GenesisError, GenesisFortranError, GenesisValidationError
from ..s_molecule import SMolecule
from .conftest import (
    CHIGNOLIN_DCD, CHIGNOLIN_PDB, CHIGNOLIN_PSF,
    has_module, requires_chignolin,
)

pytestmark = requires_chignolin


def _atom_names(mol):
    return [''.join(n).strip() for n in mol.atom_name]


@pytest.fixture(scope="module")
def allatom():
    return SMolecule.from_file(pdb=CHIGNOLIN_PDB, psf=CHIGNOLIN_PSF)


@pytest.fixture(scope="module")
def proa(allatom):
    """Fitted PROA trajectory and its molecule (the workhorse of the demo)."""
    trajs, mol = genesis_exe.crd_convert(
        molecule=allatom,
        trj_files=[str(CHIGNOLIN_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="segid:PROA",
        fitting_selection="segid:PROA",
        fitting_method="TR+ROT",
        pbc_correct="NO",
    )
    return trajs[0], mol


def test_smolecule_loading(allatom):
    assert allatom.num_atoms > 0
    assert len(allatom.atom_name) == allatom.num_atoms
    assert allatom.atom_coord.shape == (allatom.num_atoms, 3)
    assert allatom.subset_atoms(np.array([0, 1, 2])).num_atoms == 3


def test_crd_convert_selections(allatom, proa):
    traj, mol = proa
    assert traj.nframe > 0
    assert traj.natom == mol.num_atoms

    trajs_ca, mol_ca = genesis_exe.crd_convert(
        molecule=allatom, trj_files=[str(CHIGNOLIN_DCD)], trj_format="DCD",
        trj_type="COOR+BOX", selection="an:CA", fitting_selection="an:CA",
        fitting_method="TR+ROT", pbc_correct="NO")
    assert trajs_ca[0].nframe > 0
    assert trajs_ca[0].natom == mol_ca.num_atoms
    assert all(n == "CA" for n in _atom_names(mol_ca))

    _trajs_heavy, mol_heavy = genesis_exe.crd_convert(
        molecule=allatom, trj_files=[str(CHIGNOLIN_DCD)], trj_format="DCD",
        trj_type="COOR+BOX", selection="heavy", fitting_selection="heavy",
        fitting_method="TR+ROT", pbc_correct="NO")
    assert not any(n.startswith("H") for n in _atom_names(mol_heavy))


def test_trj_analysis(proa):
    traj, mol = proa
    ca = genesis_exe.selection(mol, "an:CA")
    assert len(ca) >= 4
    result = genesis_exe.trj_analysis(
        trajs=traj,
        distance_pairs=np.array([[ca[0], ca[1]]], dtype=np.int32),
        angle_triplets=np.array([[ca[0], ca[1], ca[2]]], dtype=np.int32),
        torsion_quadruplets=np.array([[ca[0], ca[1], ca[2], ca[3]]], dtype=np.int32),
    )
    assert result.distance.shape[0] == traj.nframe
    assert result.angle.shape[0] == traj.nframe
    assert result.torsion.shape[0] == traj.nframe
    assert 2.0 < result.distance.mean() < 6.0


def test_rg_analysis(proa):
    traj, mol = proa
    rg = genesis_exe.rg_analysis(molecule=mol, trajs=traj,
                                 analysis_selection="an:CA", mass_weighted=True)
    assert rg.rg.shape[0] == traj.nframe
    assert 4.0 < rg.rg.mean() < 15.0


def test_rmsd_analysis(proa):
    traj, mol = proa
    rmsd = genesis_exe.rmsd_analysis(
        molecule=mol, trajs=traj, analysis_selection="an:CA",
        fitting_selection="an:CA", fitting_method="TR+ROT")
    assert rmsd.rmsd.shape[0] == traj.nframe
    assert rmsd.rmsd[0] < 3.0
    assert np.all(rmsd.rmsd >= 0)


def test_avecrd_analysis(proa):
    traj, mol = proa
    ave = genesis_exe.avecrd_analysis(
        mol, traj, selection_group=["segid:PROA and heavy"],
        fitting_method="TR+ROT", fitting_atom=1, check_only=False,
        num_iterations=5, analysis_atom=1)
    assert ave.pdb is not None
    assert len(ave.pdb) > 100
    assert "ATOM" in ave.pdb


def test_msd_and_diffusion_analysis(proa):
    traj, mol = proa
    msd = genesis_exe.msd_analysis(
        molecule=mol, trajs=traj, selection_group=["an:CA"],
        oversample=True, delta=min(traj.nframe - 1, 1000))
    assert msd.msd is not None
    assert np.all(msd.msd >= 0)

    # diffusion_analysis expects column 0 to be the time
    nsteps = msd.msd.shape[0]
    time_col = np.arange(nsteps, dtype=np.float64).reshape(-1, 1)
    diffusion = genesis_exe.diffusion_analysis(
        msd_data=np.hstack([time_col, msd.msd]), time_step=1.0,
        start_step=max(1, int(nsteps * 0.2)))
    assert diffusion.diffusion_coefficients is not None


def test_wham_analysis_on_synthetic_windows():
    with tempfile.TemporaryDirectory() as tmpdir:
        rows = np.arange(1, 201)
        rng = np.random.RandomState(42)
        for i, loc in enumerate([1.0, 2.0, 3.0], 1):
            values = rng.normal(loc=loc, scale=0.3, size=200)
            np.savetxt(os.path.join(tmpdir, f"window_{i}.dat"),
                       np.column_stack((rows, values)), fmt="%d %.5f")
        pmf = genesis_exe.wham_analysis(
            cvfile=os.path.join(tmpdir, "window_{:d}.dat"),
            dimension=1, nblocks=1, temperature=300.0, tolerance=1.0e-7,
            rest_function=[1], grids=[(0.0, 4.0, 81)],
            constant=[(30.0, 30.0, 30.0)], reference=[(1.0, 2.0, 3.0)],
            is_periodic=[False])
    assert pmf is not None
    assert pmf.shape[0] > 0


@pytest.mark.skipif(not has_module("mdtraj"), reason="mdtraj is not installed")
def test_mdtraj_interface(proa):
    traj, mol = proa
    assert mol.to_mdtraj_topology().n_atoms == mol.num_atoms
    assert traj.to_mdtraj_trajectory(mol).n_frames == traj.nframe


@pytest.mark.skipif(not has_module("MDAnalysis"), reason="MDAnalysis is not installed")
def test_mdanalysis_interface(proa):
    traj, mol = proa
    assert mol.to_mdanalysis_universe().atoms.n_atoms == mol.num_atoms
    assert traj.to_mdanalysis_universe(mol).trajectory.n_frames == traj.nframe


def _ca_features(traj, mol):
    ca = [i for i, n in enumerate(_atom_names(mol)) if n == "CA"]
    return traj.coords[:, ca, :].reshape(traj.nframe, -1).astype(np.float32)


@pytest.mark.skipif(not has_module("sklearn"), reason="scikit-learn is not installed")
def test_sklearn_integration(proa):
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler

    traj, mol = proa
    features = _ca_features(traj, mol)
    scaled = StandardScaler().fit_transform(features)
    assert scaled.shape == features.shape
    n_samples = min(500, traj.nframe)
    tsne = TSNE(n_components=2, random_state=42,
                perplexity=min(30.0, n_samples - 1), max_iter=250)
    assert tsne.fit_transform(scaled[:n_samples]).shape == (n_samples, 2)


@pytest.mark.skipif(not has_module("torch"), reason="PyTorch is not installed")
def test_torch_integration(proa):
    import torch

    traj, mol = proa
    features = _ca_features(traj, mol)
    assert torch.tensor(features).shape == (traj.nframe, features.shape[1])


def test_error_types():
    with pytest.raises(GenesisError):
        raise GenesisValidationError("test error")
    err = GenesisFortranError("test error", code=1, stderr_output="test stderr")
    assert isinstance(err, GenesisError)
    assert err.code == 1
    assert err.stderr_output == "test stderr"
