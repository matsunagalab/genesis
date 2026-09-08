import pytest
from .conftest import MSD_DATA
from .. import genesis_exe
import numpy as np


@pytest.mark.slow
def test_diffusion_analysis():
    """Test diffusion analysis with zerocopy interface."""
    msd = np.loadtxt(MSD_DATA, dtype=np.float64)
    assert msd.ndim == 2, "MSD data should be 2D array"
    print(f"MSD data shape: {msd.shape}")

    ndata = msd.shape[0]
    start_step = int(ndata * 0.2)  # 20% start

    result = genesis_exe.diffusion_analysis(
            msd,
            time_step=2.0,
            start_step=start_step
            )

    # Validate result structure
    assert result is not None, "Result should not be None"
    assert result.out_data is not None, "out_data should not be None"
    assert result.diffusion_coefficients is not None, "diffusion_coefficients should not be None"

    # Validate result values
    assert (result.diffusion_coefficients >= 0).all(), "Diffusion coefficients should be non-negative"

    # Validate output data shape
    n_sets = msd.shape[1] - 1  # number of MSD columns (excluding time)
    expected_out_cols = 2 * n_sets + 1  # time + (msd + fit) * n_sets
    assert result.out_data.shape == (ndata, expected_out_cols), \
        f"out_data shape mismatch: {result.out_data.shape} vs {(ndata, expected_out_cols)}"

    print(f"Diffusion coefficients: {result.diffusion_coefficients}")
    print(f"Output data shape: {result.out_data.shape}")
