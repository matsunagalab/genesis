"""Shared paths, helpers and pytest configuration for the genepie test suite.

Run the suite with ``pytest`` from the repository root (or
``pytest src/genepie/tests/test_rmsd.py`` for one module). Tests are plain
``test_*`` functions; the helpers below replace the former CustomTestCase.

The shared library must be discoverable: either ``make install`` it, or point
``GENEPIE_LIB_DIR`` at ``src/analysis/interface/python_interface/.libs``.
"""
import importlib.util
import os
import pathlib
from typing import List, Tuple, Union

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------------
TEST_DIR = pathlib.Path(__file__).parent
DATA_DIR = TEST_DIR / "data"
# Reference data shared with the CLI regression harness (source checkout only)
REPO_ROOT = TEST_DIR.parent.parent.parent
REGRESSION_ROOT = REPO_ROOT / "tests" / "regression_test"
ANALYSIS_REF_ROOT = REGRESSION_ROOT / "test_analysis"
ANALYSIS_TRJ_ROOT = ANALYSIS_REF_ROOT / "trajectories"
ATDYN_BUILD_ROOT = REGRESSION_ROOT / "build"
ATDYN_PARAM_ROOT = REGRESSION_ROOT / "param"

# BPTI system paths
BPTI_PDB = DATA_DIR / "bpti" / "BPTI_ionize.pdb"
BPTI_PSF = DATA_DIR / "bpti" / "BPTI_ionize.psf"
BPTI_DCD = DATA_DIR / "bpti" / "BPTI_run.dcd"

# RALP-DPPC system paths
RALP_PDB = DATA_DIR / "ralp_dppc" / "RALP_DPPC_run.pdb"
RALP_PSF = DATA_DIR / "ralp_dppc" / "RALP_DPPC.psf"
RALP_DCD = DATA_DIR / "ralp_dppc" / "RALP_DPPC_run.dcd"

# Chignolin system paths (for integration tests - downloaded from Google Drive)
CHIGNOLIN_PDB = DATA_DIR / "chignolin" / "chignolin.pdb"
CHIGNOLIN_PSF = DATA_DIR / "chignolin" / "chignolin.psf"
CHIGNOLIN_DCD = DATA_DIR / "chignolin" / "chignolin.dcd"

# T-REMD trialanine (Ala3) tutorial data (downloaded from Google Drive).
# Populated by ``python -m genepie.tests.download_tremd_data``; see that module.
TREMD_DIR = DATA_DIR / "remd_alat_tutorial"
TREMD_NREPLICA = 20
TREMD_PDB = TREMD_DIR / "trialanine.pdb"
TREMD_PSF = TREMD_DIR / "ala3.psf"
# Per-state (parameter-sorted) inputs; {} expands to the 1-based state index.
TREMD_POT_PATTERN = str(TREMD_DIR / "remd_paramID{}.pot")
TREMD_DCD_PATTERN = str(TREMD_DIR / "remd_paramID{}_trialanine.dcd")
# Reference MBAR results produced in the original paper.
TREMD_REF_DIR = TREMD_DIR / "reference"
TREMD_FENE_REF = TREMD_REF_DIR / "fene.dat"
TREMD_WEIGHT_PATTERN = str(TREMD_REF_DIR / "weight{}.dat")
TREMD_TOR_PATTERN = str(TREMD_REF_DIR / "remd_paramID{}.tor")
# Temperature ladder (K) of the 20 parameter-sorted states.
TREMD_TEMPERATURES = [
    300.00, 302.53, 305.09, 307.65, 310.24,
    312.85, 315.47, 318.12, 320.78, 323.46,
    326.16, 328.87, 331.61, 334.37, 337.14,
    339.94, 342.75, 345.59, 348.44, 351.26,
]
TREMD_TARGET_TEMPERATURE = 300.0

# Other test data
MOLECULE_PDB = DATA_DIR / "molecule.pdb"
MSD_DATA = DATA_DIR / "msd.data"


# ---------------------------------------------------------------------------
# Skip markers for optional data and optional dependencies
# ---------------------------------------------------------------------------
def has_module(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


requires_regression_data = pytest.mark.skipif(
    not ANALYSIS_REF_ROOT.is_dir(),
    reason="tests/regression_test is only available in a source checkout")
requires_atdyn_data = pytest.mark.skipif(
    not ATDYN_BUILD_ROOT.is_dir(),
    reason="tests/regression_test/build is only available in a source checkout")
requires_chignolin = pytest.mark.skipif(
    not CHIGNOLIN_DCD.is_file(),
    reason="chignolin data missing (python -m genepie.tests.download_test_data)")
requires_tremd = pytest.mark.skipif(
    not (TREMD_DIR / "README.md").is_file(),
    reason="T-REMD data missing (python -m genepie.tests.download_tremd_data)")
requires_mdtraj = pytest.mark.skipif(
    not has_module("mdtraj"), reason="mdtraj is not installed")
requires_mdanalysis = pytest.mark.skipif(
    not has_module("MDAnalysis"), reason="MDAnalysis is not installed")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
PathLike = Union[str, bytes, os.PathLike]


def load_trajectories(dcd: PathLike, **molecule_files) -> Tuple[List, object]:
    """Load ``dcd`` with crd_convert (all atoms, no fitting) and return
    ``(trajectories, molecule)``. ``molecule_files`` are passed to
    ``SMolecule.from_file`` (pdb=..., psf=..., ...)."""
    from .. import genesis_exe
    from ..s_molecule import SMolecule

    mol = SMolecule.from_file(**molecule_files)
    trajs, _subset = genesis_exe.crd_convert(
        mol,
        trj_files=[str(dcd)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        fitting_method="NO",
        pbc_correct="NO",
    )
    return trajs, mol


def assert_objects_close(expected, actual, atol: float = 1e-4,
                         exclude: frozenset = frozenset()):
    """Compare two objects attribute by attribute: float arrays and floats
    within ``atol``, everything else exactly."""
    ekeys = set(vars(expected)) - exclude
    akeys = set(vars(actual)) - exclude
    assert ekeys == akeys, f"attribute names differ: {ekeys ^ akeys}"
    for key in sorted(ekeys):
        ev, av = vars(expected)[key], vars(actual)[key]
        if isinstance(ev, np.ndarray) and isinstance(av, np.ndarray):
            if ev.dtype.kind == "f" and av.dtype.kind == "f":
                np.testing.assert_allclose(av, ev, rtol=0.0, atol=atol,
                                           err_msg=f"{key} differs")
            else:
                assert np.array_equal(ev, av), f"{key} differs"
        elif isinstance(ev, float) and isinstance(av, float):
            assert abs(ev - av) <= atol, f"{key} differs: {ev} != {av}"
        else:
            assert ev == av, f"{key} differs: {ev!r} != {av!r}"


def assert_trajectories_close(expected, actual, atol: float = 1e-4):
    """Compare two STrajectories objects (coordinates within ``atol``)."""
    from ..s_trajectories import STrajectories

    assert isinstance(expected, STrajectories)
    assert isinstance(actual, STrajectories)
    assert_objects_close(
        expected, actual, atol=atol,
        exclude=frozenset({"c_obj", "src_c_obj", "_mem_owner",
                           "_numpy_coords", "_numpy_pbc_boxes"}))


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def bpti_molecule():
    """BPTI molecule loaded from PDB + PSF (ref = PDB)."""
    from ..s_molecule import SMolecule
    return SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)


@pytest.fixture
def bpti_trajectory(bpti_molecule):
    """One in-memory BPTI trajectory (all atoms, no fitting)."""
    from .. import genesis_exe
    trajs, _subset = genesis_exe.crd_convert(
        bpti_molecule,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
    )
    assert len(trajs) == 1
    return trajs[0]

