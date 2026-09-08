"""Running the ATDYN MD engine and energy minimiser."""

import ctypes
import io
import os
import numpy as np
import numpy.typing as npt
from typing import (
    NamedTuple,
    Optional,
)

from ..libgenesis import LibGenesis
from .. import ctrl_files
from .. import c2py_util
from .._fortran import (
    ctrl_to_bytes,
    fortran_status,
)


class AtdynMDResult(NamedTuple):
    """Result from atdyn MD simulation."""
    energies: npt.NDArray[np.float64]  # Shape: (nterms, nframes)
    final_coords: npt.NDArray[np.float64]  # Shape: (3, natom)
    energy_labels: tuple  # ('total', 'bond', 'angle', ...)


class AtdynMinResult(NamedTuple):
    """Result from atdyn energy minimization."""
    energies: npt.NDArray[np.float64]  # Shape: (nterms, nsteps)
    final_coords: npt.NDArray[np.float64]  # Shape: (3, natom)
    converged: bool
    final_gradient: float
    energy_labels: tuple  # ('total', 'bond', 'angle', ...)


_ENERGY_LABELS = ('total', 'bond', 'angle', 'urey_bradley', 'dihedral',
                  'improper', 'electrostatic', 'van_der_waals')


def run_atdyn_md(
        # Input files
        psffile: Optional[str] = None,
        prmtopfile: Optional[str] = None,
        ambcrdfile: Optional[str] = None,
        pdbfile: Optional[str] = None,
        rstfile: Optional[str] = None,
        topfile: Optional[str] = None,
        parfile: Optional[str] = None,
        strfile: Optional[str] = None,
        grotopfile: Optional[str] = None,
        grocrdfile: Optional[str] = None,
        # Energy parameters
        forcefield: Optional[str] = None,
        electrostatic: Optional[str] = None,
        switchdist: Optional[float] = None,
        cutoffdist: Optional[float] = None,
        pairlistdist: Optional[float] = None,
        vdw_force_switch: Optional[bool] = None,
        implicit_solvent: Optional[str] = None,
        output_style: Optional[str] = None,
        pme_alpha: Optional[float] = None,
        pme_ngrid_x: Optional[int] = None,
        pme_ngrid_y: Optional[int] = None,
        pme_ngrid_z: Optional[int] = None,
        pme_nspline: Optional[int] = None,
        dispersion_corr: Optional[str] = None,
        # Dynamics parameters
        integrator: Optional[str] = None,
        nsteps: Optional[int] = None,
        timestep: Optional[float] = None,
        eneout_period: Optional[int] = None,
        crdout_period: Optional[int] = None,
        rstout_period: Optional[int] = None,
        nbupdate_period: Optional[int] = None,
        iseed: Optional[int] = None,
        verbose: Optional[bool] = None,
        # Boundary parameters
        boundary_type: Optional[str] = None,
        box_size_x: Optional[float] = None,
        box_size_y: Optional[float] = None,
        box_size_z: Optional[float] = None,
        # Ensemble parameters
        ensemble: Optional[str] = None,
        tpcontrol: Optional[str] = None,
        temperature: Optional[float] = None,
        pressure: Optional[float] = None,
        gamma_t: Optional[float] = None,
        # Constraints parameters
        rigid_bond: Optional[bool] = None,
        shake_iteration: Optional[int] = None,
        shake_tolerance: Optional[float] = None,
        water_model: Optional[str] = None,
        # Output files (optional)
        dcdfile: Optional[str] = None,
        ) -> AtdynMDResult:
    """
    Run atdyn molecular dynamics simulation.

    Args:
        psffile: PSF topology file (CHARMM format)
        prmtopfile: AMBER parameter/topology file
        ambcrdfile: AMBER coordinate file
        pdbfile: PDB coordinate file
        rstfile: GENESIS restart file
        topfile: CHARMM topology RTF file
        parfile: CHARMM parameter file
        strfile: CHARMM stream file
        grotopfile: GROMACS topology file
        grocrdfile: GROMACS coordinate file (.gro)
        forcefield: Force field type (CHARMM, AMBER, GROAMBER, etc.)
        electrostatic: Electrostatic method (PME, CUTOFF, etc.)
        switchdist: Switch distance for vdW (angstrom)
        cutoffdist: Cutoff distance (angstrom)
        pairlistdist: Pairlist distance (angstrom)
        vdw_force_switch: Use force switching for vdW
        implicit_solvent: Implicit solvent model (NONE, GBSA, EEF1, etc.)
        output_style: Output style (GENESIS, CHARMM, etc.)
        pme_alpha: PME alpha parameter
        pme_ngrid_x: PME grid points in X
        pme_ngrid_y: PME grid points in Y
        pme_ngrid_z: PME grid points in Z
        pme_nspline: PME spline order
        dispersion_corr: Dispersion correction (NO, epress, etc.)
        integrator: Integrator (VVER, LEAP, etc.)
        nsteps: Number of MD steps
        timestep: Time step (ps)
        eneout_period: Energy output period
        crdout_period: Coordinate output period
        rstout_period: Restart output period
        nbupdate_period: Nonbond list update period
        iseed: Random seed
        verbose: Verbose output
        boundary_type: Boundary type (PBC, NOBC)
        box_size_x: Box size X (angstrom)
        box_size_y: Box size Y (angstrom)
        box_size_z: Box size Z (angstrom)
        ensemble: Ensemble type (NVE, NVT, NPT, etc.)
        tpcontrol: Thermostat/barostat (LANGEVIN, BUSSI, etc.)
        temperature: Temperature (K)
        pressure: Pressure (atm)
        gamma_t: Langevin friction coefficient
        rigid_bond: Use rigid bonds (SHAKE/SETTLE)
        shake_iteration: Maximum SHAKE iterations
        shake_tolerance: SHAKE convergence tolerance
        water_model: Water model for SETTLE (TIP3, WAT, SOL, etc.)
        dcdfile: Output DCD trajectory file

    Returns:
        AtdynMDResult containing:
            - energies: Array of energy terms (nterms x nframes)
            - final_coords: Final coordinates (3 x natom)
            - energy_labels: Tuple of energy term names

    Raises:
        GenesisValidationError: If input validation fails
        GenesisFortranError: If Fortran code returns an error
    """
    # === INPUT VALIDATION ===
    from ..file_validators import (
        validate_file_exists,
        validate_topology_combination,
    )
    from ..param_validators import (
        validate_enum,
        validate_positive,
        validate_distance_ordering,
        validate_pme_params,
        validate_shake_params,
        validate_ensemble_params,
        validate_pbc_params,
        FORCEFIELDS,
        ELECTROSTATICS,
        INTEGRATORS,
        BOUNDARY_TYPES,
        IMPLICIT_SOLVENTS,
    )

    # Validate topology - at least one format required
    validate_topology_combination(psffile, prmtopfile, grotopfile)

    # Validate input files exist
    validate_file_exists(psffile, "psffile", required=False)
    validate_file_exists(prmtopfile, "prmtopfile", required=False)
    validate_file_exists(ambcrdfile, "ambcrdfile", required=False)
    validate_file_exists(pdbfile, "pdbfile", required=False)
    validate_file_exists(rstfile, "rstfile", required=False)
    validate_file_exists(topfile, "topfile", required=False)
    validate_file_exists(parfile, "parfile", required=False)
    validate_file_exists(strfile, "strfile", required=False)
    validate_file_exists(grotopfile, "grotopfile", required=False)
    validate_file_exists(grocrdfile, "grocrdfile", required=False)

    # Validate enum parameters
    validate_enum(forcefield, FORCEFIELDS, "forcefield")
    validate_enum(electrostatic, ELECTROSTATICS, "electrostatic")
    validate_enum(integrator, INTEGRATORS, "integrator")
    validate_enum(boundary_type, BOUNDARY_TYPES, "boundary_type")
    validate_enum(implicit_solvent, IMPLICIT_SOLVENTS, "implicit_solvent")

    # Validate numeric parameters
    validate_positive(nsteps, "nsteps", allow_none=False)
    validate_positive(timestep, "timestep")
    validate_positive(eneout_period, "eneout_period")
    validate_positive(crdout_period, "crdout_period")
    validate_positive(rstout_period, "rstout_period")
    validate_positive(nbupdate_period, "nbupdate_period")

    # Validate distance ordering
    validate_distance_ordering(switchdist, cutoffdist, pairlistdist)

    # Validate conditional parameters
    validate_pme_params(
        electrostatic, pme_alpha, pme_ngrid_x, pme_ngrid_y, pme_ngrid_z, pme_nspline
    )
    validate_shake_params(rigid_bond, shake_iteration, shake_tolerance)
    validate_ensemble_params(ensemble, temperature, pressure, tpcontrol)
    validate_pbc_params(boundary_type, box_size_x, box_size_y, box_size_z)
    # === END VALIDATION ===

    result_energies_c = ctypes.c_void_p(None)
    result_nframes_c = ctypes.c_int(0)
    result_nterms_c = ctypes.c_int(0)
    result_final_coords_c = ctypes.c_void_p(None)
    result_natom_c = ctypes.c_int(0)

    try:
        # Build control string
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_input(
            ctrl,
            psffile=psffile,
            prmtopfile=prmtopfile,
            ambcrdfile=ambcrdfile,
            pdbfile=pdbfile,
            rstfile=rstfile,
            topfile=topfile,
            parfile=parfile,
            strfile=strfile,
            grotopfile=grotopfile,
            grocrdfile=grocrdfile,
        )
        ctrl_files.write_ctrl_output(
            ctrl,
            dcdfile=dcdfile,
        )
        ctrl_files.write_ctrl_energy(
            ctrl,
            forcefield=forcefield,
            electrostatic=electrostatic,
            switchdist=switchdist,
            cutoffdist=cutoffdist,
            pairlistdist=pairlistdist,
            vdw_force_switch=vdw_force_switch,
            implicit_solvent=implicit_solvent,
            output_style=output_style,
            pme_alpha=pme_alpha,
            pme_ngrid_x=pme_ngrid_x,
            pme_ngrid_y=pme_ngrid_y,
            pme_ngrid_z=pme_ngrid_z,
            pme_nspline=pme_nspline,
            dispersion_corr=dispersion_corr,
        )
        ctrl_files.write_ctrl_dynamics(
            ctrl,
            integrator=integrator,
            nsteps=nsteps,
            timestep=timestep,
            eneout_period=eneout_period,
            crdout_period=crdout_period,
            rstout_period=rstout_period,
            nbupdate_period=nbupdate_period,
            iseed=iseed,
            verbose=verbose,
        )
        ctrl_files.write_ctrl_boundary(
            ctrl,
            type=boundary_type,
            box_size_x=box_size_x,
            box_size_y=box_size_y,
            box_size_z=box_size_z,
        )
        ctrl_files.write_ctrl_ensemble(
            ctrl,
            ensemble=ensemble,
            tpcontrol=tpcontrol,
            temperature=temperature,
            pressure=pressure,
            gamma_t=gamma_t,
        )
        ctrl_files.write_ctrl_constraints(
            ctrl,
            rigid_bond=rigid_bond,
            shake_iteration=shake_iteration,
            shake_tolerance=shake_tolerance,
            water_model=water_model,
        )

        ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)

        with fortran_status() as (status, msg, msglen):
            LibGenesis().lib.atdyn_md_c(
                ctrl_bytes,
                ctrl_len,
                ctypes.byref(result_energies_c),
                ctypes.byref(result_nframes_c),
                ctypes.byref(result_nterms_c),
                ctypes.byref(result_final_coords_c),
                ctypes.byref(result_natom_c),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )

        # Convert results to numpy arrays
        nframes = result_nframes_c.value
        nterms = result_nterms_c.value
        natom = result_natom_c.value

        energies = c2py_util.conv_double_ndarray(
            result_energies_c, [nterms, nframes])
        final_coords = c2py_util.conv_double_ndarray(
            result_final_coords_c, [3, natom])

        return AtdynMDResult(
            energies=energies,
            final_coords=final_coords,
            energy_labels=_ENERGY_LABELS,
        )

    finally:
        # Deallocate results
        LibGenesis().lib.deallocate_atdyn_results_c()


def run_atdyn_min(
        # Input files
        psffile: Optional[str] = None,
        prmtopfile: Optional[str] = None,
        ambcrdfile: Optional[str] = None,
        pdbfile: Optional[str] = None,
        rstfile: Optional[str] = None,
        topfile: Optional[str] = None,
        parfile: Optional[str] = None,
        strfile: Optional[str] = None,
        grotopfile: Optional[str] = None,
        grocrdfile: Optional[str] = None,
        # Energy parameters
        forcefield: Optional[str] = None,
        electrostatic: Optional[str] = None,
        switchdist: Optional[float] = None,
        cutoffdist: Optional[float] = None,
        pairlistdist: Optional[float] = None,
        vdw_force_switch: Optional[bool] = None,
        implicit_solvent: Optional[str] = None,
        output_style: Optional[str] = None,
        pme_alpha: Optional[float] = None,
        pme_ngrid_x: Optional[int] = None,
        pme_ngrid_y: Optional[int] = None,
        pme_ngrid_z: Optional[int] = None,
        pme_nspline: Optional[int] = None,
        dispersion_corr: Optional[str] = None,
        # Minimize parameters
        method: Optional[str] = None,
        nsteps: Optional[int] = None,
        eneout_period: Optional[int] = None,
        crdout_period: Optional[int] = None,
        rstout_period: Optional[int] = None,
        nbupdate_period: Optional[int] = None,
        force_scale_init: Optional[float] = None,
        force_scale_max: Optional[float] = None,
        verbose: Optional[bool] = None,
        tol_rmsg: Optional[float] = None,
        tol_maxg: Optional[float] = None,
        # Constraints parameters
        rigid_bond: Optional[bool] = None,
        # Boundary parameters
        boundary_type: Optional[str] = None,
        box_size_x: Optional[float] = None,
        box_size_y: Optional[float] = None,
        box_size_z: Optional[float] = None,
        # Output files (optional)
        dcdfile: Optional[str] = None,
        ) -> AtdynMinResult:
    """
    Run atdyn energy minimization.

    Args:
        psffile: PSF topology file (CHARMM format)
        prmtopfile: AMBER parameter/topology file
        ambcrdfile: AMBER coordinate file
        pdbfile: PDB coordinate file
        rstfile: GENESIS restart file
        topfile: CHARMM topology RTF file
        parfile: CHARMM parameter file
        strfile: CHARMM stream file
        grotopfile: GROMACS topology file
        grocrdfile: GROMACS coordinate file (.gro)
        forcefield: Force field type (CHARMM, AMBER, GROAMBER, etc.)
        electrostatic: Electrostatic method (PME, CUTOFF, etc.)
        switchdist: Switch distance for vdW (angstrom)
        cutoffdist: Cutoff distance (angstrom)
        pairlistdist: Pairlist distance (angstrom)
        vdw_force_switch: Use force switching for vdW
        implicit_solvent: Implicit solvent model (NONE, GBSA, EEF1, etc.)
        output_style: Output style (GENESIS, CHARMM, etc.)
        pme_alpha: PME alpha parameter
        pme_ngrid_x: PME grid points in X
        pme_ngrid_y: PME grid points in Y
        pme_ngrid_z: PME grid points in Z
        pme_nspline: PME spline order
        dispersion_corr: Dispersion correction (NO, epress, etc.)
        method: Minimization method (SD, LBFGS)
        nsteps: Maximum number of minimization steps
        eneout_period: Energy output period
        crdout_period: Coordinate output period
        rstout_period: Restart output period
        nbupdate_period: Nonbond list update period
        force_scale_init: Initial force scale (SD)
        force_scale_max: Maximum force scale (SD)
        verbose: Verbose output
        tol_rmsg: RMS gradient tolerance for convergence
        tol_maxg: Max gradient tolerance for convergence
        boundary_type: Boundary type (PBC, NOBC)
        box_size_x: Box size X (angstrom)
        box_size_y: Box size Y (angstrom)
        box_size_z: Box size Z (angstrom)
        dcdfile: Output DCD trajectory file

    Returns:
        AtdynMinResult containing:
            - energies: Array of energy terms (nterms x nsteps)
            - final_coords: Final coordinates (3 x natom)
            - converged: Whether minimization converged
            - final_gradient: Final RMS gradient
            - energy_labels: Tuple of energy term names

    Raises:
        GenesisValidationError: If input validation fails
        GenesisFortranError: If Fortran code returns an error
    """
    # === INPUT VALIDATION ===
    from ..file_validators import (
        validate_file_exists,
        validate_topology_combination,
    )
    from ..param_validators import (
        validate_enum,
        validate_positive,
        validate_distance_ordering,
        validate_pme_params,
        validate_pbc_params,
        FORCEFIELDS,
        ELECTROSTATICS,
        MINIMIZERS,
        BOUNDARY_TYPES,
        IMPLICIT_SOLVENTS,
    )

    # Validate topology - at least one format required
    validate_topology_combination(psffile, prmtopfile, grotopfile)

    # Validate input files exist
    validate_file_exists(psffile, "psffile", required=False)
    validate_file_exists(prmtopfile, "prmtopfile", required=False)
    validate_file_exists(ambcrdfile, "ambcrdfile", required=False)
    validate_file_exists(pdbfile, "pdbfile", required=False)
    validate_file_exists(rstfile, "rstfile", required=False)
    validate_file_exists(topfile, "topfile", required=False)
    validate_file_exists(parfile, "parfile", required=False)
    validate_file_exists(strfile, "strfile", required=False)
    validate_file_exists(grotopfile, "grotopfile", required=False)
    validate_file_exists(grocrdfile, "grocrdfile", required=False)

    # Validate enum parameters
    validate_enum(forcefield, FORCEFIELDS, "forcefield")
    validate_enum(electrostatic, ELECTROSTATICS, "electrostatic")
    validate_enum(method, MINIMIZERS, "method")
    validate_enum(boundary_type, BOUNDARY_TYPES, "boundary_type")
    validate_enum(implicit_solvent, IMPLICIT_SOLVENTS, "implicit_solvent")

    # Validate numeric parameters
    validate_positive(nsteps, "nsteps", allow_none=False)
    validate_positive(eneout_period, "eneout_period")
    validate_positive(crdout_period, "crdout_period")
    validate_positive(rstout_period, "rstout_period")
    validate_positive(nbupdate_period, "nbupdate_period")
    validate_positive(force_scale_init, "force_scale_init")
    validate_positive(force_scale_max, "force_scale_max")
    validate_positive(tol_rmsg, "tol_rmsg")
    validate_positive(tol_maxg, "tol_maxg")

    # Validate distance ordering
    validate_distance_ordering(switchdist, cutoffdist, pairlistdist)

    # Validate conditional parameters
    validate_pme_params(
        electrostatic, pme_alpha, pme_ngrid_x, pme_ngrid_y, pme_ngrid_z, pme_nspline
    )
    validate_pbc_params(boundary_type, box_size_x, box_size_y, box_size_z)
    # === END VALIDATION ===

    result_energies_c = ctypes.c_void_p(None)
    result_nsteps_c = ctypes.c_int(0)
    result_nterms_c = ctypes.c_int(0)
    result_final_coords_c = ctypes.c_void_p(None)
    result_natom_c = ctypes.c_int(0)
    result_converged_c = ctypes.c_int(0)
    result_final_gradient_c = ctypes.c_double(0.0)

    try:
        # Build control string
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_input(
            ctrl,
            psffile=psffile,
            prmtopfile=prmtopfile,
            ambcrdfile=ambcrdfile,
            pdbfile=pdbfile,
            rstfile=rstfile,
            topfile=topfile,
            parfile=parfile,
            strfile=strfile,
            grotopfile=grotopfile,
            grocrdfile=grocrdfile,
        )
        ctrl_files.write_ctrl_output(
            ctrl,
            dcdfile=dcdfile,
        )
        ctrl_files.write_ctrl_energy(
            ctrl,
            forcefield=forcefield,
            electrostatic=electrostatic,
            switchdist=switchdist,
            cutoffdist=cutoffdist,
            pairlistdist=pairlistdist,
            vdw_force_switch=vdw_force_switch,
            implicit_solvent=implicit_solvent,
            output_style=output_style,
            pme_alpha=pme_alpha,
            pme_ngrid_x=pme_ngrid_x,
            pme_ngrid_y=pme_ngrid_y,
            pme_ngrid_z=pme_ngrid_z,
            pme_nspline=pme_nspline,
            dispersion_corr=dispersion_corr,
        )
        ctrl_files.write_ctrl_minimize(
            ctrl,
            method=method,
            nsteps=nsteps,
            eneout_period=eneout_period,
            crdout_period=crdout_period,
            rstout_period=rstout_period,
            nbupdate_period=nbupdate_period,
            force_scale_init=force_scale_init,
            force_scale_max=force_scale_max,
            verbose=verbose,
            tol_rmsg=tol_rmsg,
            tol_maxg=tol_maxg,
        )
        ctrl_files.write_ctrl_constraints(
            ctrl,
            rigid_bond=rigid_bond,
        )
        ctrl_files.write_ctrl_boundary(
            ctrl,
            type=boundary_type,
            box_size_x=box_size_x,
            box_size_y=box_size_y,
            box_size_z=box_size_z,
        )

        ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)

        with fortran_status() as (status, msg, msglen):
            LibGenesis().lib.atdyn_min_c(
                ctrl_bytes,
                ctrl_len,
                ctypes.byref(result_energies_c),
                ctypes.byref(result_nsteps_c),
                ctypes.byref(result_nterms_c),
                ctypes.byref(result_final_coords_c),
                ctypes.byref(result_natom_c),
                ctypes.byref(result_converged_c),
                ctypes.byref(result_final_gradient_c),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )

        # Convert results to numpy arrays
        nsteps_out = result_nsteps_c.value
        nterms = result_nterms_c.value
        natom = result_natom_c.value

        energies = c2py_util.conv_double_ndarray(
            result_energies_c, [nterms, nsteps_out])
        final_coords = c2py_util.conv_double_ndarray(
            result_final_coords_c, [3, natom])

        return AtdynMinResult(
            energies=energies,
            final_coords=final_coords,
            converged=bool(result_converged_c.value),
            final_gradient=result_final_gradient_c.value,
            energy_labels=_ENERGY_LABELS,
        )

    finally:
        # Deallocate results
        LibGenesis().lib.deallocate_atdyn_results_c()


def _import_root() -> str:
    """Return the directory the ``genepie`` package can be imported from.

    Derived from the package itself rather than from this module's path so that
    it stays correct regardless of how deeply nested this module is.
    """
    from .. import __file__ as genepie_init
    return os.path.dirname(os.path.dirname(os.path.abspath(genepie_init)))


def _run_isolated_subprocess(
    func_name: str,
    result_fields: list,
    timeout: Optional[float],
    task_description: str,
    **kwargs
) -> dict:
    """
    Common helper for running atdyn functions in isolated subprocess.

    Args:
        func_name: Name of the genesis_exe function to call
        result_fields: List of field names to extract from result
        timeout: Maximum time in seconds to wait
        task_description: Description for error messages (e.g., "MD simulation")
        **kwargs: Arguments to pass to the function

    Returns:
        Dictionary with the result fields
    """
    import subprocess
    import sys
    import pickle
    import base64

    kwargs_bytes = base64.b64encode(pickle.dumps(kwargs)).decode('ascii')
    result_fields_str = repr(result_fields)

    script = f'''
import sys
import pickle
import base64

try:
    from genepie import genesis_exe
except ImportError:
    sys.path.insert(0, "{_import_root()}")
    from genepie import genesis_exe

kwargs = pickle.loads(base64.b64decode("{kwargs_bytes}"))

try:
    result = genesis_exe.{func_name}(**kwargs)
    output = {{"success": True}}
    for field in {result_fields_str}:
        val = getattr(result, field)
        output[field] = val.tolist() if hasattr(val, 'tolist') else val
except Exception as e:
    import traceback
    output = {{
        "success": False,
        "error": str(e),
        "error_type": type(e).__name__,
        "traceback": traceback.format_exc(),
    }}

sys.stdout.buffer.write(base64.b64encode(pickle.dumps(output)))
'''

    try:
        proc = subprocess.run(
            [sys.executable, '-c', script],
            capture_output=True,
            timeout=timeout,
            env={**os.environ, 'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS', '1')},
        )
    except subprocess.TimeoutExpired as e:
        raise TimeoutError(f"atdyn {task_description} timed out after {timeout} seconds") from e

    if proc.returncode != 0:
        # GENESIS errors normally come back as a pickled GenesisFortranError
        # (the engine runs under the library error guard), so a bare non-zero
        # exit means the engine died outside it: a Fortran runtime error, a
        # crash, or a signal. Its own output was discarded by the child, so
        # point at the way to see it.
        stderr_text = proc.stderr.decode('utf-8', errors='replace').strip()
        if proc.returncode < 0:
            import signal
            try:
                how = f"was killed by {signal.Signals(-proc.returncode).name}"
            except ValueError:
                how = f"was killed by signal {-proc.returncode}"
        else:
            how = f"exited with code {proc.returncode}"
        raise RuntimeError(
            f"atdyn {task_description} subprocess {how}.\n"
            + (f"stderr:\n{stderr_text}\n" if stderr_text else "")
            + "The engine's own output is not captured; run the same input "
            "with the atdyn command-line program to see it."
        )

    try:
        output = pickle.loads(base64.b64decode(proc.stdout))
    except Exception as e:
        raise RuntimeError(
            f"Failed to decode subprocess output: {e}\n"
            f"stdout: {proc.stdout[:500]}\n"
            f"stderr: {proc.stderr.decode('utf-8', errors='replace')}"
        )

    if not output["success"]:
        from .. import exceptions as _exc
        error_type = output.get("error_type", "")
        error_msg = output.get("error", "Unknown error")

        # Re-raise the same typed error the child saw (e.g. GenesisFortranFileError),
        # falling back to the base class for a subclass this module does not know.
        exc_cls = getattr(_exc, error_type, None)
        if isinstance(exc_cls, type) and issubclass(exc_cls, _exc.GenesisFortranError):
            raise exc_cls(error_msg)
        elif "GenesisFortran" in error_type:
            raise _exc.GenesisFortranError(error_msg)
        elif error_type == "GenesisValidationError":
            raise _exc.GenesisValidationError(error_msg)
        else:
            raise RuntimeError(f"{error_type}: {error_msg}\n{output.get('traceback', '')}")

    return output


def run_atdyn_md_isolated(
    timeout: Optional[float] = None,
    **kwargs
) -> AtdynMDResult:
    """
    Run atdyn MD simulation in an isolated subprocess.

    This function runs the MD simulation in a separate Python subprocess,
    providing maximum isolation from Fortran errors or crashes. Use this
    when running multiple sequential simulations or when reliability is
    more important than performance.

    Args:
        timeout: Maximum time in seconds to wait for completion (None = no limit)
        **kwargs: All arguments passed to run_atdyn_md()

    Returns:
        AtdynMDResult containing energies, final coordinates, and labels

    Raises:
        GenesisValidationError: If input validation fails
        GenesisFortranError: If Fortran code returns an error
        TimeoutError: If simulation exceeds timeout
        RuntimeError: If subprocess fails unexpectedly

    Example:
        >>> # Run multiple simulations safely
        >>> for i in range(10):
        ...     result = run_atdyn_md_isolated(
        ...         prmtopfile="system.prmtop",
        ...         ambcrdfile=f"frame_{i}.rst",
        ...         nsteps=1000,
        ...         timeout=300,  # 5 minute timeout
        ...     )
    """
    output = _run_isolated_subprocess(
        func_name="run_atdyn_md",
        result_fields=["energies", "final_coords", "energy_labels"],
        timeout=timeout,
        task_description="MD simulation",
        **kwargs
    )
    return AtdynMDResult(
        energies=np.array(output["energies"]),
        final_coords=np.array(output["final_coords"]),
        energy_labels=output["energy_labels"],
    )


def run_atdyn_min_isolated(
    timeout: Optional[float] = None,
    **kwargs
) -> AtdynMinResult:
    """
    Run atdyn energy minimization in an isolated subprocess.

    This function runs the minimization in a separate Python subprocess,
    providing maximum isolation from Fortran errors or crashes.

    Args:
        timeout: Maximum time in seconds to wait for completion (None = no limit)
        **kwargs: All arguments passed to run_atdyn_min()

    Returns:
        AtdynMinResult containing energies, final coordinates, convergence info

    Raises:
        GenesisValidationError: If input validation fails
        GenesisFortranError: If Fortran code returns an error
        TimeoutError: If minimization exceeds timeout
        RuntimeError: If subprocess fails unexpectedly
    """
    output = _run_isolated_subprocess(
        func_name="run_atdyn_min",
        result_fields=["energies", "final_coords", "converged", "final_gradient", "energy_labels"],
        timeout=timeout,
        task_description="minimization",
        **kwargs
    )
    return AtdynMinResult(
        energies=np.array(output["energies"]),
        final_coords=np.array(output["final_coords"]),
        converged=output["converged"],
        final_gradient=output["final_gradient"],
        energy_labels=output["energy_labels"],
    )
