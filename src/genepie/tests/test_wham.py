"""wham_analysis regression tests against the CLI reference output."""
from unittest import mock

import numpy as np
import pytest

from .. import genesis_exe
from ..exceptions import (
    GenesisFortranFileError,
    GenesisFortranNotSupportedError,
    GenesisValidationError,
)
from .conftest import ANALYSIS_REF_ROOT, ANALYSIS_TRJ_ROOT, requires_regression_data

pytestmark = requires_regression_data

# Same tolerance the CLI regression harness uses for this test case, see
# tests/regression_test/test_analysis/test_wham_analysis/triala_cv/config.ini
TOLERANCE = 1.0e-8
# 14 umbrella windows along an end-to-end distance of trialanine.
CONSTANT = (1.2,) * 14
REFERENCE = (1.80, 2.72, 3.64, 4.56, 5.48, 6.40, 7.32,
             8.24, 9.16, 10.08, 11.00, 11.92, 12.84, 13.76)
CVFILE = str(ANALYSIS_TRJ_ROOT / "triala_cv" / "{}.dis")


def run_wham(**overrides):
    kwargs = dict(
        cvfile=CVFILE,
        dimension=1,
        nblocks=1,
        temperature=300.0,
        tolerance=10E-08,
        rest_function=(1,),
        grids=((0.0, 15.0, 301),),
        constant=(CONSTANT,),
        reference=(REFERENCE,),
        is_periodic=(False,),
    )
    kwargs.update(overrides)
    return genesis_exe.wham_analysis(**kwargs)


def test_wham_analysis_matches_cli_reference():
    """PMF must reproduce the value the CLI wham_analysis writes."""
    pmf = run_wham()
    ref = np.loadtxt(ANALYSIS_REF_ROOT / "test_wham_analysis" / "triala_cv" / "ref")
    assert ref.shape == pmf.shape
    assert np.all(np.isfinite(pmf))
    # Column 0 holds the bin centers, column 1 the free energy.
    np.testing.assert_allclose(pmf, ref, rtol=0, atol=TOLERANCE)


def test_wham_analysis_rejects_dcd_input():
    """DCD input is unsupported because molecules are never defined."""
    with pytest.raises(GenesisFortranNotSupportedError):
        run_wham(dcdfile="whatever.dcd")


def test_wham_analysis_requires_cvfile():
    with pytest.raises(GenesisValidationError):
        run_wham(cvfile=None)


def test_wham_analysis_rejects_missing_cvfile():
    """A nonexistent cvfile must raise, not abort the whole process.

    The Fortran open_file() calls exit(1) on a missing input file, so the
    wrapper has to check existence up front.
    """
    with pytest.raises(GenesisValidationError):
        run_wham(cvfile="/no/such/dir/{}.dis")


def test_wham_fortran_catches_missing_cvfile():
    """Even without the Python pre-check, Fortran must not exit(1).

    The Python guard (_validate_cvfiles_exist) is bypassed here so the
    missing path reaches Fortran's open_file(). The library-mode error
    guard converts the would-be exit(1) into a catchable
    GenesisFortranFileError (ERROR_FILE_NOT_FOUND=201) instead of killing
    the interpreter.
    """
    with mock.patch(
        "genepie.analysis.free_energy._validate_cvfiles_exist",
        lambda *a, **k: None,
    ):
        with pytest.raises(GenesisFortranFileError) as excinfo:
            run_wham(cvfile="/no/such/dir/{}.dis")
    assert excinfo.value.code == 201
