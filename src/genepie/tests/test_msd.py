import pytest
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..s_molecule import SMolecule
from .. import genesis_exe


@pytest.mark.slow
def test_msd_analysis():
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
        d = genesis_exe.msd_analysis(
            mol, t,
            selection_group=["rnam:TIP3"],
            selection=[1],
            mode=["ALL"],
            oversample=True,
            delta=9,
        )
        # Validate MSD results
        assert d.msd is not None, "MSD result should not be None"
        assert d.msd.shape[0] > 0, "MSD result should have at least one time point"
        assert d.msd.shape[1] > 0, "MSD result should have at least one molecule"
        # MSD should be non-negative
        assert (d.msd >= 0).all(), "MSD values should be non-negative"
        print(f"MSD shape: {d.msd.shape}, min={d.msd.min():.3f}, max={d.msd.max():.3f}", flush=True)
