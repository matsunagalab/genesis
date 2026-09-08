"""ATDYN (in-process MD engine) smoke tests on the regression-test systems.

The first test drives ``run_atdyn_min`` directly. The remaining ones use the
``*_isolated`` variants, which run ATDYN in a fresh subprocess: that is the
recommended way to run several simulations from one interpreter because the
engine keeps global Fortran state between calls.
"""
import os

import pytest

from .. import genesis_exe
from .conftest import ATDYN_BUILD_ROOT, ATDYN_PARAM_ROOT, requires_atdyn_data

pytestmark = requires_atdyn_data


@pytest.fixture(autouse=True)
def _single_thread(monkeypatch):
    monkeypatch.setenv("OMP_NUM_THREADS", "1")


def _system_dir(name):
    d = ATDYN_BUILD_ROOT / name
    if not d.is_dir():
        pytest.skip(f"regression system {name} not found under {ATDYN_BUILD_ROOT}")
    return str(d)


def _param(name):
    return str(ATDYN_PARAM_ROOT / name)


def test_atdyn_min_glycam():
    """Energy minimization with an AMBER (glycam) system."""
    test_dir = _system_dir("glycam")
    result = genesis_exe.run_atdyn_min(
        prmtopfile=os.path.join(test_dir, "glycam.top"),
        ambcrdfile=os.path.join(test_dir, "glycam.rst"),
        rstfile=os.path.join(test_dir, "rst"),
        forcefield="AMBER",
        electrostatic="PME",
        switchdist=12.0,
        cutoffdist=12.0,
        pairlistdist=14.0,
        pme_alpha=0.34,
        pme_ngrid_x=64,
        pme_ngrid_y=64,
        pme_ngrid_z=64,
        pme_nspline=4,
        dispersion_corr="epress",
        method="SD",
        nsteps=20,
        eneout_period=2,
        nbupdate_period=4,
        rigid_bond=False,
        boundary_type="PBC",
        box_size_x=69.5294360,
        box_size_y=68.0597930,
        box_size_z=56.2256950,
    )
    assert result.final_coords.shape[1] > 0
    assert result.energies[0, 0] < 0, "total energy should be negative"


def test_atdyn_md_glycam():
    """MD with an AMBER (glycam) system, PME."""
    test_dir = _system_dir("glycam")
    result = genesis_exe.run_atdyn_md_isolated(
        prmtopfile=os.path.join(test_dir, "glycam.top"),
        ambcrdfile=os.path.join(test_dir, "glycam.rst"),
        rstfile=os.path.join(test_dir, "rst"),
        forcefield="AMBER",
        electrostatic="PME",
        switchdist=12.0,
        cutoffdist=12.0,
        pairlistdist=14.0,
        pme_alpha=0.34,
        pme_ngrid_x=64,
        pme_ngrid_y=64,
        pme_ngrid_z=64,
        pme_nspline=4,
        dispersion_corr="epress",
        integrator="VVER",
        nsteps=20,
        timestep=0.001,
        eneout_period=2,
        nbupdate_period=5,
        iseed=314159,
        verbose=True,
        rigid_bond=True,
        shake_iteration=500,
        shake_tolerance=1.0e-10,
        water_model="WAT",
        ensemble="NVE",
        tpcontrol="NO",
        temperature=0,
        boundary_type="PBC",
        box_size_x=69.5294360,
        box_size_y=68.0597930,
        box_size_z=56.2256950,
    )
    assert result.energies[0, 0] < 0, "total energy should be negative"


def test_atdyn_md_bpti():
    """MD with a GROMACS (bpti) system, PME."""
    test_dir = _system_dir("bpti")
    result = genesis_exe.run_atdyn_md_isolated(
        grotopfile=os.path.join(test_dir, "bpti.top"),
        grocrdfile=os.path.join(test_dir, "bpti.gro"),
        rstfile=os.path.join(test_dir, "rst"),
        forcefield="GROAMBER",
        electrostatic="PME",
        switchdist=12.0,
        cutoffdist=12.0,
        pairlistdist=14.0,
        pme_alpha=0.34,
        pme_ngrid_x=64,
        pme_ngrid_y=64,
        pme_ngrid_z=64,
        pme_nspline=4,
        output_style="GENESIS",
        integrator="VVER",
        nsteps=20,
        timestep=0.001,
        eneout_period=2,
        nbupdate_period=5,
        iseed=314159,
        verbose=True,
        rigid_bond=True,
        shake_iteration=500,
        shake_tolerance=1.0e-10,
        water_model="SOL",
        ensemble="NVE",
        tpcontrol="NO",
        temperature=0,
        boundary_type="PBC",
        box_size_x=65.3318,
        box_size_y=65.3318,
        box_size_z=65.3318,
    )
    assert result.energies[0, 0] < 0, "total energy should be negative"


def test_atdyn_md_jac_param27():
    """MD with a CHARMM (jac_param27) system, PME."""
    test_dir = _system_dir("jac_param27")
    result = genesis_exe.run_atdyn_md_isolated(
        topfile=_param("top_all27_prot_lipid.rtf"),
        parfile=_param("par_all27_prot_lipid.prm"),
        psffile=os.path.join(test_dir, "jac_param27.psf"),
        pdbfile=os.path.join(test_dir, "jac_param27.pdb"),
        rstfile=os.path.join(test_dir, "rst"),
        forcefield="CHARMM",
        electrostatic="PME",
        switchdist=8.0,
        cutoffdist=10.0,
        pairlistdist=12.0,
        pme_alpha=0.34,
        pme_ngrid_x=64,
        pme_ngrid_y=64,
        pme_ngrid_z=64,
        pme_nspline=4,
        vdw_force_switch=False,
        output_style="GENESIS",
        integrator="LEAP",
        nsteps=20,
        timestep=0.001,
        eneout_period=2,
        nbupdate_period=5,
        iseed=314159,
        verbose=True,
        rigid_bond=True,
        shake_iteration=500,
        shake_tolerance=1.0e-10,
        water_model="TIP3",
        ensemble="NVE",
        tpcontrol="NO",
        temperature=0,
        boundary_type="PBC",
        box_size_x=65.5,
        box_size_y=65.5,
        box_size_z=65.5,
    )
    assert result.energies[0, 0] < 0, "total energy should be negative"


def test_atdyn_md_dppc_nvt():
    """MD with a CHARMM (dppc) system in the NVT ensemble (Langevin)."""
    test_dir = _system_dir("dppc")
    result = genesis_exe.run_atdyn_md_isolated(
        topfile=_param("top_all36_lipid.rtf"),
        parfile=_param("par_all36_lipid.prm"),
        strfile=_param("toppar_water_ions.str"),
        psffile=os.path.join(test_dir, "dppc.psf"),
        pdbfile=os.path.join(test_dir, "dppc.pdb"),
        rstfile=os.path.join(test_dir, "rst"),
        forcefield="CHARMM",
        electrostatic="PME",
        switchdist=10.0,
        cutoffdist=12.0,
        pairlistdist=13.5,
        pme_alpha=0.34,
        pme_ngrid_x=72,
        pme_ngrid_y=72,
        pme_ngrid_z=72,
        pme_nspline=4,
        output_style="GENESIS",
        integrator="VVER",
        nsteps=20,
        timestep=0.001,
        eneout_period=2,
        nbupdate_period=5,
        iseed=314159,
        verbose=True,
        rigid_bond=True,
        shake_iteration=500,
        shake_tolerance=1.0e-10,
        ensemble="NVT",
        tpcontrol="LANGEVIN",
        temperature=300.0,
        pressure=1.0,
        boundary_type="PBC",
        box_size_x=69.4792,
        box_size_y=69.4792,
        box_size_z=71.6508,
    )
    assert result.energies[0, 0] < 0, "total energy should be negative"


def test_atdyn_min_dppc():
    """Minimization with a CHARMM (dppc) system."""
    test_dir = _system_dir("dppc")
    result = genesis_exe.run_atdyn_min_isolated(
        topfile=_param("top_all36_lipid.rtf"),
        parfile=_param("par_all36_lipid.prm"),
        strfile=_param("toppar_water_ions.str"),
        psffile=os.path.join(test_dir, "dppc.psf"),
        pdbfile=os.path.join(test_dir, "dppc.pdb"),
        rstfile=os.path.join(test_dir, "rst"),
        forcefield="CHARMM",
        electrostatic="PME",
        switchdist=10.0,
        cutoffdist=12.0,
        pairlistdist=13.5,
        pme_alpha=0.34,
        pme_ngrid_x=72,
        pme_ngrid_y=72,
        pme_ngrid_z=72,
        pme_nspline=4,
        output_style="GENESIS",
        method="SD",
        nsteps=20,
        eneout_period=2,
        nbupdate_period=4,
        rigid_bond=False,
        boundary_type="PBC",
        box_size_x=69.4792,
        box_size_y=69.4792,
        box_size_z=71.6508,
    )
    assert result.energies[0, 0] < 0, "total energy should be negative"
