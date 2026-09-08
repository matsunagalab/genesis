"""Round trips between STrajectories and mdtraj trajectories."""
from ..s_trajectories import STrajectories
from .conftest import (
    BPTI_DCD, BPTI_PDB, BPTI_PSF, assert_trajectories_close,
    load_trajectories, requires_mdtraj,
)

pytestmark = requires_mdtraj


def test_from_mdtraj_trajectory():
    import mdtraj

    mdt = mdtraj.load(str(BPTI_DCD), top=str(BPTI_PDB))
    trj, _mol = STrajectories.from_mdtraj_trajectory(mdt)
    gtrajs, _gmol = load_trajectories(BPTI_DCD, pdb=BPTI_PDB)
    assert_trajectories_close(gtrajs[0], trj)


def test_to_mdtraj_trajectory():
    strajs, smol = load_trajectories(BPTI_DCD, pdb=BPTI_PDB, psf=BPTI_PSF)
    for t in strajs:
        mdt = t.to_mdtraj_trajectory(smol)
        gt, _gm = STrajectories.from_mdtraj_trajectory(mdt)
        assert_trajectories_close(t, gt)
