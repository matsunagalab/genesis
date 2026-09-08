import numpy as np
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_rg_analysis():
    """Test RG analysis with zerocopy interface."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    trajs, subset_mol = genesis_exe.crd_convert(
            mol,
            trj_files=[str(BPTI_DCD)],
            trj_format="DCD",
            trj_type="COOR+BOX",
            selection="all",
    )
    _ = subset_mol

    for t in trajs:
        result = genesis_exe.rg_analysis(
            mol, t,
            analysis_selection="all",
            ana_period=1,
            mass_weighted=True,
        )

        # Validate Rg results
        assert result.rg is not None, "Rg result should not be None"
        assert len(result.rg) > 0, "Rg result should have at least one frame"
        assert all(r > 0 for r in result.rg), "Rg values should be positive"
        # Rg should be in reasonable range for proteins (5-100 Angstroms)
        assert all(r < 100.0 for r in result.rg), "Rg values should be reasonable (< 100 A)"

        print(f"Rg values (n={len(result.rg)}): min={min(result.rg):.3f}, max={max(result.rg):.3f}")


def test_rg_lazy():
    """Test RG analysis with lazy loading using unified API."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)

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

    # Run RG analysis with lazy trajectory
    result = genesis_exe.rg_analysis(
        subset_mol,
        lazy_traj,
        analysis_selection="all",
        ana_period=1,
        mass_weighted=True,
    )

    # Validate Rg results
    assert result.rg is not None, "Rg result should not be None"
    assert len(result.rg) > 0, "Rg result should have at least one frame"
    assert all(r > 0 for r in result.rg), "Rg values should be positive"
    assert all(r < 100.0 for r in result.rg), "Rg values should be reasonable (< 100 A)"

    print(f"Lazy Rg values (n={len(result.rg)}): min={min(result.rg):.3f}, max={max(result.rg):.3f}")


def test_rg_lazy_vs_memory():
    """Compare lazy RG analysis with memory-based RG analysis."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)

    # Memory-based RG
    trajs_mem, mol_mem = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=False,
    )
    result_mem = genesis_exe.rg_analysis(
        mol_mem,
        trajs_mem[0],
        analysis_selection="an:CA",
        ana_period=1,
        mass_weighted=True,
    )

    # Lazy-based RG
    trajs_lazy, mol_lazy = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )
    result_lazy = genesis_exe.rg_analysis(
        mol_lazy,
        trajs_lazy[0],
        analysis_selection="an:CA",
        ana_period=1,
        mass_weighted=True,
    )

    # Compare results
    assert len(result_mem.rg) == len(result_lazy.rg), \
        f"Frame count mismatch: memory={len(result_mem.rg)}, lazy={len(result_lazy.rg)}"

    for i, (mem_val, lazy_val) in enumerate(zip(result_mem.rg, result_lazy.rg)):
        assert np.isclose(mem_val, lazy_val, rtol=1e-4, atol=1e-6), \
            f"Frame {i}: memory={mem_val}, lazy={lazy_val}"

    print(f"Memory vs Lazy RG comparison passed: {len(result_mem.rg)} frames")
    print(f"  Memory: min={min(result_mem.rg):.5f}, max={max(result_mem.rg):.5f}")
    print(f"  Lazy:   min={min(result_lazy.rg):.5f}, max={max(result_lazy.rg):.5f}")
