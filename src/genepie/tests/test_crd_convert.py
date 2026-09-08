import numpy as np
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_crd_convert_basic():
    """Test basic trajectory loading without fitting."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
    )

    assert len(trajs) == 1
    traj = trajs[0]
    assert traj.nframe > 0
    assert traj.natom == mol.num_atoms
    assert traj.coords.shape == (traj.nframe, traj.natom, 3)
    print(f"  Loaded {traj.nframe} frames, {traj.natom} atoms")


def test_crd_convert_with_selection():
    """Test trajectory loading with atom selection."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="an:CA",
    )

    assert len(trajs) == 1
    traj = trajs[0]
    assert traj.nframe > 0
    # Subset should have fewer atoms
    assert traj.natom == subset_mol.num_atoms
    assert traj.natom < mol.num_atoms
    print(f"  Selected {traj.natom} CA atoms from {mol.num_atoms} total atoms")


def test_crd_convert_with_fitting():
    """Test trajectory loading with TR+ROT fitting."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="an:CA",
        fitting_selection="an:CA",
        fitting_method="TR+ROT",
    )

    assert len(trajs) == 1
    traj = trajs[0]
    assert traj.nframe > 0
    assert traj.natom == subset_mol.num_atoms
    # Check coordinates are reasonable (not NaN/Inf)
    assert np.all(np.isfinite(traj.coords))
    print(f"  Fitted trajectory with {traj.nframe} frames")


def test_crd_convert_with_centering():
    """Test trajectory loading with centering."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="an:CA",
        centering=True,
        center_coord=(0.0, 0.0, 0.0),
    )

    assert len(trajs) == 1
    traj = trajs[0]
    # Check that coordinates are centered (COM near origin)
    com = np.mean(traj.coords, axis=1)  # COM per frame
    # COM should be near (0, 0, 0) for each frame
    assert np.allclose(com, 0.0, atol=5.0)  # Allow 5 Angstrom tolerance
    print(f"  Centered trajectory with {traj.nframe} frames")


def test_crd_convert_info():
    """Test crd_convert_info function."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    info = genesis_exe.crd_convert_info(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
    )

    assert len(info.frame_counts) == 1
    assert info.frame_counts[0] > 0
    print(f"  Trajectory info: {info.frame_counts[0]} frames")
