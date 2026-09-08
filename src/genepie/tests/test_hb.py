"""hb_analysis: hydrogen-bond counts per atom and per snapshot."""
import pytest

from .. import genesis_exe
from ..s_molecule import SMolecule
from .conftest import RALP_DCD, RALP_PDB, RALP_PSF


def _ralp_trajectories():
    mol = SMolecule.from_file(pdb=RALP_PDB, psf=RALP_PSF)
    trajs, _subset = genesis_exe.crd_convert(
        mol,
        trj_files=[str(RALP_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        centering=True,
        centering_selection="all",
        center_coord=(0.0, 0.0, 0.0),
        rename_res=["HSE HIS", "HSD HIS"],
    )
    return mol, trajs


@pytest.mark.parametrize("output_type", ["Count_atom", "Count_Snap"])
def test_hb_analysis(output_type):
    mol, trajs = _ralp_trajectories()
    for t in trajs:
        d = genesis_exe.hb_analysis(
            mol, t,
            selection_group=["sid:PROA",
                             "resname:DPPC & (an:O11 | an:O12 | an:O13 | an:O14)"],
            check_only=False,
            output_type=output_type,
            solvent_list="DPPC",
            analysis_atom=1,
            target_atom=2,
            boundary_type="PBC",
            hb_distance=3.4,
            dha_angle=120.0,
            hda_angle=30.0,
        )
        assert d is not None
        print(d, flush=True)
