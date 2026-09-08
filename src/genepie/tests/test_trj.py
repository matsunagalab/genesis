"""trj_analysis: distances, angles, torsions and their COM variants."""
import numpy as np
import pytest

from .. import genesis_exe
from ..s_molecule import SMolecule
from .conftest import (
    ANALYSIS_REF_ROOT, BPTI_DCD, BPTI_PDB, BPTI_PSF,
    load_trajectories, requires_regression_data,
)


def _ca_measurements(mol):
    """Distance pairs, angle triplets and torsion quads on the first CA atoms."""
    ca = [genesis_exe.selection(mol, f"rno:{i} and an:CA")[0] for i in (1, 2, 3, 4)]
    dist_pairs = np.array([[ca[0], ca[1]], [ca[1], ca[2]]], dtype=np.int32)
    angle_triplets = np.array([[ca[0], ca[1], ca[2]]], dtype=np.int32)
    torsion_quads = np.array([[ca[0], ca[1], ca[2], ca[3]]], dtype=np.int32)
    return dist_pairs, angle_triplets, torsion_quads


@requires_regression_data
def test_trj_analysis_atoms():
    """Atom-based measurements must match the CLI reference output."""
    trajs, mol = load_trajectories(BPTI_DCD, pdb=BPTI_PDB, psf=BPTI_PSF)
    dist_pairs, angle_triplets, torsion_quads = _ca_measurements(mol)
    ref_root = ANALYSIS_REF_ROOT / "test_trj_analysis"

    for t in trajs:
        result = genesis_exe.trj_analysis(
            t,
            distance_pairs=dist_pairs,
            angle_triplets=angle_triplets,
            torsion_quadruplets=torsion_quads,
        )

        assert result.distance.shape == (t.nframe, 2)
        assert result.angle.shape == (t.nframe, 1)
        assert result.torsion.shape == (t.nframe, 1)
        assert np.all(result.distance > 0)
        assert np.all((result.angle >= 0) & (result.angle <= 180))

        ref_dist = np.loadtxt(ref_root / "Distance/ref")
        ref_ang = np.loadtxt(ref_root / "Angle/ref")
        ref_tor = np.loadtxt(ref_root / "Dihedral/ref")
        np.testing.assert_allclose(result.distance, ref_dist[:, 1:], rtol=0, atol=1e-3)
        np.testing.assert_allclose(result.angle, ref_ang[:, 1:], rtol=0, atol=1e-3)
        np.testing.assert_allclose(result.torsion, ref_tor[:, 1:], rtol=0, atol=1e-3)


@pytest.mark.parametrize("ana_period", [1, 2])
def test_trj_lazy_vs_memory(ana_period):
    """Lazy trj_analysis must match the memory-based result."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    dist_pairs, angle_triplets, torsion_quads = _ca_measurements(mol)
    common = dict(trj_files=[str(BPTI_DCD)], trj_format="DCD",
                  trj_type="COOR+BOX", selection="all")

    mem_trajs, _ = genesis_exe.crd_convert(mol, lazy=False, **common)
    lazy_trajs, _ = genesis_exe.crd_convert(mol, lazy=True, **common)
    assert lazy_trajs[0].is_lazy

    kwargs = dict(distance_pairs=dist_pairs, angle_triplets=angle_triplets,
                  torsion_quadruplets=torsion_quads, ana_period=ana_period)
    mem = genesis_exe.trj_analysis(mem_trajs[0], **kwargs)
    lazy = genesis_exe.trj_analysis(lazy_trajs[0], **kwargs)

    assert mem.distance.shape == lazy.distance.shape
    np.testing.assert_allclose(lazy.distance, mem.distance, rtol=1e-4, atol=1e-5)
    np.testing.assert_allclose(lazy.angle, mem.angle, rtol=1e-4, atol=1e-5)
    np.testing.assert_allclose(lazy.torsion, mem.torsion, rtol=1e-4, atol=1e-5)


def _residue_groups(mol):
    res = [list(genesis_exe.selection(mol, f"rno:{i}")) for i in (1, 2, 3, 4)]
    cdis_groups = [(res[0], res[1])]
    cang_groups = [(res[0], res[1], res[2])]
    ctor_groups = [(res[0], res[1], res[2], res[3])]
    return cdis_groups, cang_groups, ctor_groups


@pytest.mark.parametrize("ana_period", [1, 2])
def test_trj_com_lazy_vs_memory(ana_period):
    """Lazy COM trj_analysis must match the memory-based result."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    ca1 = genesis_exe.selection(mol, "rno:1 and an:CA")[0]
    ca2 = genesis_exe.selection(mol, "rno:2 and an:CA")[0]
    cdis_groups, cang_groups, ctor_groups = _residue_groups(mol)
    common = dict(trj_files=[str(BPTI_DCD)], trj_format="DCD",
                  trj_type="COOR+BOX", selection="all")

    mem_trajs, _ = genesis_exe.crd_convert(mol, lazy=False, **common)
    lazy_trajs, _ = genesis_exe.crd_convert(mol, lazy=True, **common)
    assert lazy_trajs[0].is_lazy

    # Mix atom-based and COM-based measurements to exercise every stream.
    kwargs = dict(
        distance_pairs=np.array([[ca1, ca2]], dtype=np.int32),
        cdis_groups=cdis_groups, cang_groups=cang_groups,
        ctor_groups=ctor_groups, molecule=mol, ana_period=ana_period,
    )
    mem = genesis_exe.trj_analysis(mem_trajs[0], **kwargs)
    lazy = genesis_exe.trj_analysis(lazy_trajs[0], **kwargs)

    assert mem.distance.shape == lazy.distance.shape
    np.testing.assert_allclose(lazy.distance, mem.distance, rtol=1e-4, atol=1e-5)
    np.testing.assert_allclose(lazy.com_distance, mem.com_distance, rtol=1e-4, atol=1e-5)
    np.testing.assert_allclose(lazy.com_angle, mem.com_angle, rtol=1e-4, atol=1e-5)
    np.testing.assert_allclose(lazy.com_torsion, mem.com_torsion, rtol=1e-4, atol=1e-5)


def test_trj_analysis_com():
    """COM-based measurements return finite values of the expected shape."""
    trajs, mol = load_trajectories(BPTI_DCD, pdb=BPTI_PDB, psf=BPTI_PSF)
    cdis_groups, cang_groups, _ctor_groups = _residue_groups(mol)

    for t in trajs:
        result = genesis_exe.trj_analysis(
            t, cdis_groups=cdis_groups, cang_groups=cang_groups, molecule=mol)
        assert result.com_distance.shape == (t.nframe, 1)
        assert result.com_angle.shape == (t.nframe, 1)
        assert np.all(result.com_distance > 0)
        assert np.all((result.com_angle >= 0) & (result.com_angle <= 180))
        assert np.all(np.isfinite(result.com_distance))
        assert np.all(np.isfinite(result.com_angle))


def test_trj_analysis_mixed():
    """Atom-based and COM-based measurements in one call."""
    trajs, mol = load_trajectories(BPTI_DCD, pdb=BPTI_PDB, psf=BPTI_PSF)
    ca1 = genesis_exe.selection(mol, "rno:1 and an:CA")[0]
    ca2 = genesis_exe.selection(mol, "rno:2 and an:CA")[0]
    cdis_groups, _cang_groups, _ctor_groups = _residue_groups(mol)

    for t in trajs:
        result = genesis_exe.trj_analysis(
            t,
            distance_pairs=np.array([[ca1, ca2]], dtype=np.int32),
            cdis_groups=cdis_groups,
            molecule=mol,
        )
        assert result.distance.shape == (t.nframe, 1)
        assert result.com_distance.shape == (t.nframe, 1)
        assert np.all(result.distance > 0)
        assert np.all(result.com_distance > 0)
        assert np.all(np.isfinite(result.distance))
        assert np.all(np.isfinite(result.com_distance))
