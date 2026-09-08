"""Mean squared displacement and diffusion coefficients."""

import ctypes
import io
from collections import namedtuple
import numpy as np
import numpy.typing as npt
from typing import (
    Iterable,
    Optional,
)

from ..libgenesis import LibGenesis
from ..s_molecule import SMolecule
from ..s_trajectories import STrajectories
from .. import ctrl_files
from .. import c2py_util
from ..exceptions import GenesisValidationError
from .._fortran import (
    ctrl_to_bytes,
    fortran_status,
    molecule_c,
)


MsdAnalysisResult = namedtuple(
        'MsdAnalysisResult',
        ['msd'])


def msd_analysis(
        molecule: SMolecule, trajs: STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        selection: Optional[Iterable[str]] = None,
        mode: Optional[Iterable[int]] = None,
        check_only: Optional[bool] = None,
        oversample: Optional[bool] = None,
        delta: Optional[int] = None,
        ) -> MsdAnalysisResult:
    """
    Executes msd_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        msd
    """
    ana_period_c = ctypes.c_int(ana_period)
    num_analysis_mols_c = ctypes.c_int(0)
    num_delta_c = ctypes.c_int(0)
    result_msd_c = ctypes.c_void_p()
    try:
        with molecule_c(molecule) as mol_c:
            ctrl = io.BytesIO()
            ctrl_files.write_ctrl_output(
                    ctrl,
                    msdfile="dummy.msd")
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl_files.write_ctrl_molecule_selection(
                    ctrl, selection, mode)
            ctrl.write(b"[OPTION]\n")
            ctrl_files.write_kwargs(
                    ctrl,
                    check_only=check_only,
                    oversample=oversample,
                    delta=delta,
                    )

            ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)
            with fortran_status() as (status, msgbuf, msglen):
                LibGenesis().lib.ma_analysis_c(
                        ctypes.byref(mol_c),
                        ctypes.byref(trajs.get_c_obj()),
                        ctypes.byref(ana_period_c),
                        ctrl_bytes,
                        ctrl_len,
                        ctypes.byref(result_msd_c),
                        ctypes.byref(num_analysis_mols_c),
                        ctypes.byref(num_delta_c),
                        ctypes.byref(status),
                        msgbuf,
                        ctypes.c_int(msglen),
                        )
        result_msd = c2py_util.conv_double_ndarray(
            result_msd_c, [num_delta_c.value, num_analysis_mols_c.value])
        return MsdAnalysisResult(
                result_msd)
    finally:
        if result_msd_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_msd_c),
                    ctypes.byref(num_delta_c),
                    ctypes.byref(num_analysis_mols_c))


DiffusionAnalysisResult = namedtuple(
        'DiffusionAnalysisResult',
        ['out_data', 'diffusion_coefficients'])


def diffusion_analysis(
        msd_data: npt.NDArray[np.float64],
        time_step: float = 1.0,
        distance_unit: float = 1.0,
        ndofs: int = 3,
        start_step: int = 1,
        stop_step: Optional[int] = None,
        ) -> DiffusionAnalysisResult:
    """
    Executes diffusion analysis on mean square displacement data.

    This function analyzes MSD data to calculate diffusion coefficients
    using linear fitting. It uses zero-copy memory sharing with pre-allocated
    result arrays for optimal performance.

    Args:
        msd_data: 2D numpy array of shape (ndata, ncols) where:
                  - column 0: time steps (integers or floats)
                  - columns 1+: MSD values for each set
        time_step: Time per step in ps (default: 1.0)
        distance_unit: Distance unit factor (default: 1.0 for Angstrom)
        ndofs: Degrees of freedom (default: 3 for 3D diffusion)
        start_step: Start step for linear fitting (1-indexed, default: 1)
        stop_step: Stop step for linear fitting (1-indexed, default: ndata)

    Returns:
        DiffusionAnalysisResult containing:
        - out_data: 2D array with time, MSD, and fitted values
        - diffusion_coefficients: 1D array of diffusion coefficients (cm^2/s)

    Example:
        >>> # MSD data with time column and one MSD set
        >>> msd = np.array([[0, 0.0], [1, 0.1], [2, 0.3], ...])
        >>> result = diffusion_analysis(msd, time_step=0.01)
        >>> print(f"D = {result.diffusion_coefficients[0]:.2e} cm^2/s")
    """
    lib = LibGenesis().lib

    # Validate input
    if msd_data.ndim != 2:
        raise GenesisValidationError(f"msd_data must be 2D, got {msd_data.ndim}D")
    if msd_data.shape[1] < 2:
        raise GenesisValidationError("msd_data must have at least 2 columns")

    ndata = msd_data.shape[0]
    ncols = msd_data.shape[1]
    n_sets = ncols - 1

    # Set default stop_step
    if stop_step is None:
        stop_step = ndata

    # Ensure array is Fortran order (column-major)
    msd_f = np.asfortranarray(msd_data.T, dtype=np.float64)

    # Get pointer (zero-copy)
    msd_ptr = msd_f.ctypes.data_as(ctypes.c_void_p)

    # Pre-allocate output arrays (full zero-copy)
    out_ncols = 2 * n_sets + 1  # time + (msd + fit) * n_sets
    out_data_f = np.zeros((out_ncols, ndata), dtype=np.float64, order='F')
    diff_coeff = np.zeros(n_sets, dtype=np.float64)

    out_data_ptr = out_data_f.ctypes.data_as(ctypes.c_void_p)
    diff_coeff_ptr = diff_coeff.ctypes.data_as(ctypes.c_void_p)


    with fortran_status() as (status, msg, msglen):
        lib.diffusion_analysis_c(
            msd_ptr,
            ctypes.c_int(ndata),
            ctypes.c_int(ncols),
            ctypes.c_double(time_step),
            ctypes.c_double(distance_unit),
            ctypes.c_int(ndofs),
            ctypes.c_int(start_step),
            ctypes.c_int(stop_step),
            out_data_ptr,
            ctypes.c_int(out_ncols * ndata),
            diff_coeff_ptr,
            ctypes.c_int(n_sets),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )


    # Transpose to (ndata, ncols) for Python convention
    out_data = out_data_f.T

    return DiffusionAnalysisResult(out_data, diff_coeff)
