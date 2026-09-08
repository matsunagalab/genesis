"""Round trips between STrajectories and MDAnalysis universes."""
from ..s_trajectories import STrajectories
from .conftest import (
    BPTI_DCD, BPTI_PDB, BPTI_PSF, assert_trajectories_close,
    load_trajectories, requires_mdanalysis,
)

pytestmark = requires_mdanalysis


def test_from_mdanalysis_universe():
    import MDAnalysis as mda

    uni = mda.Universe(BPTI_PDB, BPTI_DCD)
    uni.atoms.guess_bonds()
    trj, _mol = STrajectories.from_mdanalysis_universe(uni)
    gtrajs, _gmol = load_trajectories(BPTI_DCD, pdb=BPTI_PDB, psf=BPTI_PSF)
    assert_trajectories_close(gtrajs[0], trj)


def test_to_mdanalysis_universe():
    gtrajs, gmol = load_trajectories(BPTI_DCD, pdb=BPTI_PDB, psf=BPTI_PSF)
    for t in gtrajs:
        uni = t.to_mdanalysis_universe(gmol)
        gt, _gm = STrajectories.from_mdanalysis_universe(uni)
        assert_trajectories_close(t, gt)
