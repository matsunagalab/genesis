import numpy as np
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_crd_convert_lazy_basic():
    """Test basic crd_convert with lazy=True."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)

    # Create lazy trajectory
    lazy_trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )

    assert len(lazy_trajs) == 1, "Should return one lazy trajectory"
    lazy_traj = lazy_trajs[0]

    # Check lazy trajectory properties
    assert lazy_traj.is_lazy, "Trajectory should be lazy"
    assert lazy_traj.lazy_dcd_file is not None, "lazy_dcd_file should be set"
    assert lazy_traj.lazy_trj_type is not None, "lazy_trj_type should be set"
    assert lazy_traj.nframe > 0, "nframe should be positive"
    assert lazy_traj.natom > 0, "natom should be positive"

    # Lazy trajectory should NOT have coords loaded
    assert lazy_traj.coords is None, "Lazy trajectory should not have coords loaded"

    print(f"Lazy trajectory: nframe={lazy_traj.nframe}, natom={lazy_traj.natom}")
    print(f"  lazy_dcd_file: {lazy_traj.lazy_dcd_file}")
    print(f"  lazy_trj_type: {lazy_traj.lazy_trj_type}")


def test_crd_convert_lazy_with_selection():
    """Test crd_convert lazy with atom selection."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)

    # Create lazy trajectory with CA selection
    lazy_trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="an:CA",
        lazy=True,
    )

    lazy_traj = lazy_trajs[0]
    assert lazy_traj.is_lazy, "Trajectory should be lazy"
    assert lazy_traj.selection_indices is not None, "selection_indices should be set"

    # selection_indices should contain CA atoms
    n_ca = len(lazy_traj.selection_indices)
    assert n_ca > 0, "Should have CA atoms"
    assert n_ca < mol.num_atoms, "CA atoms should be subset of all atoms"
    assert lazy_traj.natom == n_ca, "Lazy trajectory should expose selected atoms"
    assert lazy_traj.lazy_dcd_natom == mol.num_atoms, \
        "Lazy trajectory should retain the physical DCD atom count"

    # Check subset molecule
    assert subset_mol.num_atoms == n_ca, f"Subset molecule should have {n_ca} atoms"

    print(f"Lazy trajectory with selection: natom={lazy_traj.natom}, selected={n_ca}")
    print(f"  selection_indices (first 5): {lazy_traj.selection_indices[:5]}")


def test_crd_convert_lazy_selection_and_period_parity():
    """Selected lazy views and two-stage periods must match memory mode."""
    mol = SMolecule.from_file(
        pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB
    )
    common = dict(
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="an:CA",
        ana_period=2,
    )
    mem_trajs, mem_mol = genesis_exe.crd_convert(
        mol, lazy=False, **common
    )
    lazy_trajs, lazy_mol = genesis_exe.crd_convert(
        mol, lazy=True, **common
    )

    mem_traj = mem_trajs[0]
    lazy_traj = lazy_trajs[0]
    assert lazy_traj.natom == lazy_mol.num_atoms == mem_traj.natom
    assert lazy_traj.nframe == mem_traj.nframe
    assert lazy_traj.lazy_ana_period == 2

    mem_rmsd = genesis_exe.rmsd_analysis(
        mem_mol, mem_traj, analysis_selection="all", ana_period=2
    ).rmsd
    lazy_rmsd = genesis_exe.rmsd_analysis(
        lazy_mol, lazy_traj, analysis_selection="all", ana_period=2
    ).rmsd
    np.testing.assert_allclose(lazy_rmsd, mem_rmsd, rtol=1e-4, atol=1e-6)

    mem_rg = genesis_exe.rg_analysis(
        mem_mol, mem_traj, analysis_selection="all", ana_period=2
    ).rg
    lazy_rg = genesis_exe.rg_analysis(
        lazy_mol, lazy_traj, analysis_selection="all", ana_period=2
    ).rg
    np.testing.assert_allclose(lazy_rg, mem_rg, rtol=1e-4, atol=1e-6)


def test_crd_convert_lazy_vs_memory_nframe():
    """Test that lazy and memory crd_convert report same nframe."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)

    # Memory mode
    mem_trajs, _ = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=False,
    )

    # Lazy mode
    lazy_trajs, _ = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )

    # Compare frame counts
    assert mem_trajs[0].nframe == lazy_trajs[0].nframe, \
        f"Frame count mismatch: memory={mem_trajs[0].nframe}, lazy={lazy_trajs[0].nframe}"
    assert mem_trajs[0].natom == lazy_trajs[0].natom, \
        f"Atom count mismatch: memory={mem_trajs[0].natom}, lazy={lazy_trajs[0].natom}"

    print(f"Memory vs Lazy frame count matches: {mem_trajs[0].nframe} frames, {mem_trajs[0].natom} atoms")


def test_crd_convert_lazy_single_file_required():
    """Test that lazy=True requires single DCD file."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)

    # Lazy mode with multiple files should raise error
    try:
        genesis_exe.crd_convert(
            mol,
            trj_files=[str(BPTI_DCD), str(BPTI_DCD)],  # Two files
            trj_format="DCD",
            trj_type="COOR+BOX",
            selection="all",
            lazy=True,
        )
        assert False, "Should have raised error for multiple files with lazy=True"
    except Exception as e:
        # Should raise an error about single file requirement
        assert "lazy" in str(e).lower() or "single" in str(e).lower(), \
            f"Error should mention lazy/single file: {e}"
        print(f"Correctly raised error for multiple files: {e}")


def test_crd_convert_lazy_trj_type_coor():
    """Test crd_convert lazy with TRJ_TYPE_COOR (no box)."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)

    # Create lazy trajectory with COOR only
    lazy_trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR",  # No box
        selection="all",
        lazy=True,
    )

    lazy_traj = lazy_trajs[0]
    assert lazy_traj.is_lazy, "Trajectory should be lazy"
    # TRJ_TYPE_COOR = 1
    assert lazy_traj.lazy_trj_type == 1, f"lazy_trj_type should be 1, got {lazy_traj.lazy_trj_type}"

    print(f"Lazy trajectory with COOR: lazy_trj_type={lazy_traj.lazy_trj_type}")
