"""mbar_analysis regression tests (1-D umbrella dataset, with and without blocks)."""
import numpy as np
import pytest

from .. import genesis_exe
from ..exceptions import GenesisFortranNotSupportedError, GenesisValidationError
from .conftest import ANALYSIS_REF_ROOT, ANALYSIS_TRJ_ROOT, requires_regression_data

pytestmark = requires_regression_data

# 61 umbrella windows spaced 3 degrees apart along a periodic dihedral.
NREPLICA = 61
CONSTANT = (0.06092,) * NREPLICA
REFERENCE = tuple(3.0 * i for i in range(NREPLICA))
CVFILE = str(ANALYSIS_TRJ_ROOT / "umbrella_1d" / "{}.dat")
# The MBAR iteration converges to 1e-8 and the CLI regression harness accepts
# 0.01 across platforms, so 1e-6 catches regressions while leaving room for
# BLAS-dependent differences in the solver.
TOLERANCE = 1.0e-6

# (nblocks, reference directory): the plain 1-D case and the same dataset
# split into blocks for error estimation.
CASES = [(None, "umbrella_1d"), (5, "umbrella_block")]


def load_fene_reference(path):
    """Load a fene reference file as a 2-D (n_replica, n_blocks) array."""
    ref = np.loadtxt(path)
    if ref.ndim == 1:
        ref = ref[:, np.newaxis]
    return ref


def run_mbar(nblocks=None, **overrides):
    kwargs = dict(
        cvfile=CVFILE,
        nreplica=NREPLICA,
        input_type="US",
        dimension=1,
        temperature=300.0,
        target_temperature=300.0,
        tolerance=1E-08,
        rest_function=(1,),
        grids=((-1.0, 181.0, 81),),
        constant=(CONSTANT,),
        reference=(REFERENCE,),
        is_periodic=(True,),
        box_size=(360.0,),
    )
    if nblocks is not None:
        kwargs["nblocks"] = nblocks
    kwargs.update(overrides)
    return genesis_exe.mbar_analysis(**kwargs)


@pytest.mark.parametrize("nblocks, ref_dir", CASES, ids=[c[1] for c in CASES])
def test_mbar_analysis_matches_cli_reference(nblocks, ref_dir):
    """Free energies must reproduce what the CLI mbar_analysis writes."""
    fene = run_mbar(nblocks)
    ref = load_fene_reference(
        ANALYSIS_REF_ROOT / "test_mbar_analysis" / ref_dir / "fene.dat.ref")
    assert ref.shape == fene.shape
    assert np.all(np.isfinite(fene))
    np.testing.assert_allclose(fene, ref, rtol=0, atol=TOLERANCE)


def test_mbar_analysis_leaves_cwd_clean(tmp_path, monkeypatch):
    """The Fortran writer must not drop fene.dat into the cwd."""
    monkeypatch.chdir(tmp_path)
    run_mbar()
    assert not (tmp_path / "fene.dat").exists()


def test_mbar_analysis_rejects_dcd_input():
    """DCD input is unsupported because molecules are never defined."""
    with pytest.raises(GenesisFortranNotSupportedError):
        run_mbar(dcdfile="whatever.dcd")


def test_mbar_analysis_requires_cvfile():
    with pytest.raises(GenesisValidationError):
        run_mbar(cvfile=None)


def test_mbar_analysis_rejects_missing_cvfile():
    """A nonexistent cvfile must raise, not abort the whole process.

    The Fortran open_file() calls exit(1) on a missing input file, so the
    wrapper has to check existence up front.
    """
    with pytest.raises(GenesisValidationError):
        run_mbar(cvfile="/no/such/dir/{}.dat")
