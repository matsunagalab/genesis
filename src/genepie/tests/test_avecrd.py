from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_avecrd_analysis():
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
        d = genesis_exe.avecrd_analysis(
            mol, t,
            selection_group=["an:CA"],
            fitting_method="TR+ROT",
            fitting_atom=1,
            check_only=False,
            num_iterations=5,
            analysis_atom=1,
        )
        print(d.pdb)
