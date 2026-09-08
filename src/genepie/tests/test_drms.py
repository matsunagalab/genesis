import tempfile
from pathlib import Path
import numpy as np
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..s_molecule import SMolecule
from .. import genesis_exe
from ._dcd_test_utils import rewrite_dcd


def compute_contact_list_from_refcoord(mol, selection_indices, min_dist=1.0, max_dist=6.0, exclude_residues=4):
    """Compute contact list and reference distances from molecule reference coordinates.

    This is a simplified version of the Fortran setup_contact_list for testing purposes.
    It computes contacts between atoms within min_dist-max_dist and excludes nearby residues.
    """
    ref_coord = mol.atom_refcoord  # (num_atoms, 3), 1-indexed in mol
    residue_no = mol.residue_no

    contact_list = []
    contact_dist = []

    for i, idx_i in enumerate(selection_indices):
        for j, idx_j in enumerate(selection_indices):
            if j <= i:
                continue  # Skip duplicates and self

            # Check residue exclusion
            res_i = residue_no[idx_i - 1]  # Convert to 0-indexed
            res_j = residue_no[idx_j - 1]
            if abs(res_i - res_j) < exclude_residues:
                continue

            # Calculate distance from reference coordinates
            coord_i = ref_coord[idx_i - 1, :]  # Convert to 0-indexed
            coord_j = ref_coord[idx_j - 1, :]
            d = np.sqrt(np.sum((coord_i - coord_j) ** 2))

            if min_dist <= d < max_dist:
                contact_list.append([min(idx_i, idx_j), max(idx_i, idx_j)])
                contact_dist.append(d)

    return np.array(contact_list, dtype=np.int32).T, np.array(contact_dist, dtype=np.float64)


def test_drms_analysis():
    """Test DRMS analysis with zerocopy interface."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
    )
    _ = subset_mol

    for t in trajs:
        # Get CA atom indices for contact computation
        ca_indices = genesis_exe.selection(mol, "an: CA")
        print(f"Number of CA atoms: {len(ca_indices)}")

        # Compute contact list using Python
        contact_list, contact_dist = compute_contact_list_from_refcoord(
            mol, ca_indices,
            min_dist=1.0, max_dist=6.0, exclude_residues=4
        )
        print(f"Number of contacts: {contact_list.shape[1]}")
        print(f"contact_list shape: {contact_list.shape}")
        print(f"contact_dist shape: {contact_dist.shape}")

        # Run DRMS analysis
        result = genesis_exe.drms_analysis(
            t,
            contact_list=contact_list,
            contact_dist=contact_dist,
            ana_period=1,
            pbc_correct=False,
        )

        # Validate results
        assert result.drms is not None, "DRMS result should not be None"
        assert len(result.drms) > 0, "DRMS result should have values"
        assert all(d >= 0 for d in result.drms), "DRMS values should be non-negative"

        print(f"DRMS (n={len(result.drms)}): "
              f"min={min(result.drms):.5f}, max={max(result.drms):.5f}")


def test_drms_lazy():
    """Test DRMS analysis with lazy loading using unified API."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)

    # Create lazy trajectory via crd_convert(lazy=True)
    lazy_trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )

    lazy_traj = lazy_trajs[0]
    assert lazy_traj.is_lazy, "Trajectory should be lazy"

    # Get CA atom indices for contact computation
    ca_indices = genesis_exe.selection(mol, "an: CA")
    print(f"Number of CA atoms: {len(ca_indices)}")

    # Compute contact list using Python
    contact_list, contact_dist = compute_contact_list_from_refcoord(
        mol, ca_indices,
        min_dist=1.0, max_dist=6.0, exclude_residues=4
    )
    print(f"Number of contacts: {contact_list.shape[1]}")

    # Run DRMS analysis with lazy trajectory
    result = genesis_exe.drms_analysis(
        lazy_traj,
        contact_list=contact_list,
        contact_dist=contact_dist,
        ana_period=1,
        pbc_correct=False,
    )

    # Validate results
    assert result.drms is not None, "DRMS result should not be None"
    assert len(result.drms) > 0, "DRMS result should have values"
    assert all(d >= 0 for d in result.drms), "DRMS values should be non-negative"

    print(f"Lazy DRMS (n={len(result.drms)}): "
          f"min={min(result.drms):.5f}, max={max(result.drms):.5f}")


def test_drms_lazy_vs_memory():
    """Compare lazy DRMS analysis with memory-based DRMS analysis."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)

    # Get CA atom indices for contact computation
    ca_indices = genesis_exe.selection(mol, "an: CA")

    # Compute contact list using Python
    contact_list, contact_dist = compute_contact_list_from_refcoord(
        mol, ca_indices,
        min_dist=1.0, max_dist=6.0, exclude_residues=4
    )

    # Memory-based DRMS
    trajs_mem, mol_mem = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=False,
    )
    result_mem = genesis_exe.drms_analysis(
        trajs_mem[0],
        contact_list=contact_list,
        contact_dist=contact_dist,
        ana_period=1,
        pbc_correct=False,
    )

    # Lazy-based DRMS
    trajs_lazy, mol_lazy = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )
    result_lazy = genesis_exe.drms_analysis(
        trajs_lazy[0],
        contact_list=contact_list,
        contact_dist=contact_dist,
        ana_period=1,
        pbc_correct=False,
    )

    # Compare results
    assert len(result_mem.drms) == len(result_lazy.drms), \
        f"Frame count mismatch: memory={len(result_mem.drms)}, lazy={len(result_lazy.drms)}"

    for i, (mem_val, lazy_val) in enumerate(zip(result_mem.drms, result_lazy.drms)):
        assert np.isclose(mem_val, lazy_val, rtol=1e-4, atol=1e-6), \
            f"Frame {i}: memory={mem_val}, lazy={lazy_val}"

    print(f"Memory vs Lazy DRMS comparison passed: {len(result_mem.drms)} frames")
    print(f"  Memory: min={min(result_mem.drms):.5f}, max={max(result_mem.drms):.5f}")
    print(f"  Lazy:   min={min(result_lazy.drms):.5f}, max={max(result_lazy.drms):.5f}")


def test_drms_lazy_selected_view():
    """DRMS contact indices use the selected lazy trajectory index space."""
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
    assert lazy_trajs[0].natom == lazy_mol.num_atoms == mem_mol.num_atoms

    selected_indices = np.arange(1, lazy_mol.num_atoms + 1, dtype=np.int32)
    contact_list, contact_dist = compute_contact_list_from_refcoord(
        lazy_mol, selected_indices,
        min_dist=1.0, max_dist=8.0, exclude_residues=4,
    )
    assert contact_list.shape[1] > 0

    mem_result = genesis_exe.drms_analysis(
        mem_trajs[0], contact_list, contact_dist, ana_period=2
    ).drms
    lazy_result = genesis_exe.drms_analysis(
        lazy_trajs[0], contact_list, contact_dist, ana_period=2
    ).drms
    np.testing.assert_allclose(
        lazy_result, mem_result, rtol=1e-4, atol=1e-6
    )


def test_drms_lazy_triclinic_box_matches_memory():
    """Lazy DCD unit-cell conversion must match the standard reader."""
    mol = SMolecule.from_file(
        pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB
    )
    ca_indices = genesis_exe.selection(mol, "an:CA")
    contact_list, contact_dist = compute_contact_list_from_refcoord(
        mol, ca_indices, min_dist=1.0, max_dist=6.0, exclude_residues=4
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        triclinic_dcd = Path(tmpdir) / "triclinic.dcd"
        rewrite_dcd(
            BPTI_DCD,
            triclinic_dcd,
            box_values=(35.0, 2.0, 36.0, 1.0, 3.0, 37.0),
        )
        common = dict(
            trj_files=[str(triclinic_dcd)],
            trj_type="COOR+BOX",
            selection="all",
        )
        mem_trajs, _ = genesis_exe.crd_convert(
            mol, lazy=False, **common
        )
        lazy_trajs, _ = genesis_exe.crd_convert(
            mol, lazy=True, **common
        )
        mem_result = genesis_exe.drms_analysis(
            mem_trajs[0], contact_list, contact_dist, pbc_correct=True
        ).drms
        lazy_result = genesis_exe.drms_analysis(
            lazy_trajs[0], contact_list, contact_dist, pbc_correct=True
        ).drms
        np.testing.assert_allclose(
            lazy_result, mem_result, rtol=1e-4, atol=1e-6
        )
