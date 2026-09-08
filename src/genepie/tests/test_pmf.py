"""pmf_analysis regression tests against the CLI reference output."""
import numpy as np
import pytest

from .. import genesis_exe
from ..exceptions import GenesisValidationError
from .conftest import ANALYSIS_REF_ROOT, ANALYSIS_TRJ_ROOT, requires_regression_data

pytestmark = requires_regression_data

# Same tolerance the CLI regression harness uses for this test case, see
# tests/regression_test/test_analysis/test_pmf_analysis/pathcv/config.ini
TOLERANCE = 0.1


def _traj(name):
    return str(ANALYSIS_TRJ_ROOT / "umbrella_path" / name)


def run_pmf(**overrides):
    kwargs = dict(
        cvfile=_traj("{}.pathcv"),
        weightfile=_traj("{}.weight"),
        distfile=_traj("{}.pathdist"),
        nreplica=16,
        dimension=1,
        temperature=300.0,
        cutoff=0.04,
        grids=((1.0, 16.0, 50),),
        band_width=(0.3,),
        is_periodic=(False,),
    )
    kwargs.update(overrides)
    return genesis_exe.pmf_analysis(**kwargs)


def test_pmf_analysis_matches_cli_reference():
    """PMF must reproduce the values the CLI pmf_analysis writes."""
    res = run_pmf()
    ref = np.loadtxt(ANALYSIS_REF_ROOT / "test_pmf_analysis" / "pathcv" / "pmf.dat.ref")
    # ref columns: bin center, standard PMF, Gaussian-kernel PMF
    assert ref.shape[0] == res.cv.shape[0]
    assert np.all(np.isfinite(res.cv))
    assert np.all(np.isfinite(res.pmf))
    assert np.all(np.isfinite(res.pmf_gaussian))
    np.testing.assert_allclose(res.cv, ref[:, 0], rtol=0, atol=TOLERANCE)
    np.testing.assert_allclose(res.pmf, ref[:, 1], rtol=0, atol=TOLERANCE)
    np.testing.assert_allclose(res.pmf_gaussian, ref[:, 2], rtol=0, atol=TOLERANCE)


def test_pmf_analysis_memory_matches_file():
    """In-memory arrays must give the same PMF as the equivalent files.

    Uses a single replica (no dist cutoff) so the memory path and the
    file path see identical samples.
    """
    cv = np.loadtxt(_traj("1.pathcv"))[:, 1]
    weight = np.loadtxt(_traj("1.weight"))[:, 1]
    file_res = genesis_exe.pmf_analysis(
        cvfile=_traj("{}.pathcv"),
        weightfile=_traj("{}.weight"),
        nreplica=1,
        dimension=1,
        temperature=300.0,
        grids=((1.0, 16.0, 50),),
        band_width=(0.3,),
        is_periodic=(False,),
    )
    mem_res = genesis_exe.pmf_analysis(
        cv=cv,
        weight=weight,
        temperature=300.0,
        grids=((1.0, 16.0, 50),),
        band_width=(0.3,),
        is_periodic=(False,),
    )
    np.testing.assert_allclose(mem_res.cv, file_res.cv, rtol=0, atol=1e-8)
    np.testing.assert_allclose(mem_res.pmf, file_res.pmf, rtol=0, atol=1e-6)
    np.testing.assert_allclose(mem_res.pmf_gaussian, file_res.pmf_gaussian,
                               rtol=0, atol=1e-6)


def test_pmf_analysis_2d_runs():
    """A 2-D reweighted PMF returns a finite matrix with matching axes."""
    rng = np.random.default_rng(0)
    phi = rng.uniform(-170, 170, size=2000)
    psi = rng.uniform(-170, 170, size=2000)
    weight = rng.uniform(0.5, 1.5, size=2000)
    res = genesis_exe.pmf_analysis(
        cv=np.column_stack([phi, psi]),
        weight=weight,
        temperature=300.0,
        grids=((-180.0, 180.0, 37), (-180.0, 180.0, 37)),
        band_width=(15.0, 15.0),
        is_periodic=(True, True),
        box_size=(360.0, 360.0),
    )
    assert res.pmf.shape == (36, 36)
    assert res.cv1.shape[0] == 36
    assert res.cv2.shape[0] == 36
    assert np.all(np.isfinite(res.pmf))
    assert abs(float(res.pmf.min())) <= 1e-9


def test_pmf_analysis_requires_cv():
    with pytest.raises(GenesisValidationError):
        genesis_exe.pmf_analysis(grids=((1.0, 16.0, 50),), band_width=(0.3,))


def test_pmf_analysis_rejects_both_cv_and_cvfile():
    with pytest.raises(GenesisValidationError):
        genesis_exe.pmf_analysis(
            cv=np.zeros(10), cvfile="x{}.cv",
            grids=((1.0, 16.0, 50),), band_width=(0.3,))
