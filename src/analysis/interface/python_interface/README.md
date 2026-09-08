### Python interface for GENESIS analysis tools (June, 2025)

# Contents

* Required environment
* Download and compile source code
* Download chignolin data
* Jupyter notebook
* Regression test
* Python script editing
* Folder and file description

# Required environment

```
gfortran 10 (GCC 10) or later is supported
uv
autoconf
automake
libtool
wget
```

# Download and compile source code

Download GENESIS code from GitHub and switch to the working branch `python_interface`.

```
$ git clone https://github.com/matsunagalab/genepie.git
$ cd genepie/
$ git checkout python_interface
```

Constrct a virtual environment for python by using `uv`. In the following, please use the following virtual environment.

```
$ cd /path/to/genesis/

# For Mac users
$ brew install uv
# For Ubuntu users, see https://docs.astral.sh/uv/getting-started/installation/

$ uv venv --python=python3.11
$ source .venv/bin/activate
(genesis) $ uv pip install torch torchvision torchaudio nglview numpy mdtraj MDAnalysis plotly jupyterlab py3Dmol scikit-learn gdown

# if you want to deactivate the virtual environment, use the following command
# (genesis) $ deactivate
```

To avoid module name conflicts between mbar_analysis and msd_analysis, generate files with replaced module names related to mbar_analysis in `src/analysis/interface/mbar_analysis` (required before automake etc.)

```
$ cd src/analysis/interface/python_interface
$ python mbar_rename.py
```

Compile GENESIS. The Python interface is also compiled simultaneously (LAPACK is required for GENESIS).

```
# autoconf, automake, and libtool are required.
# In case of Mac, you can install them by the following command.
$ brew install autoconf automake libtool

(genesis) $ cd /path/to/genesis/
(genesis) $ autoscan
(genesis) $ autoheader
(genesis) $ mkdir m4
(genesis) $ aclocal
(genesis) $ autoconf
(genesis) $ libtoolize #In case of Mac, please use `glibtoolize` instead of `libtoolize`
(genesis) $ automake -a
(genesis) $ ./configure LAPACK_LIBS="-L/usr/local/lib -llapack -lblas"
# In case of Mac, please use GNU gcc by specifying CC=gcc-14 or CC=gcc-15 other verions
# e.g, $ CC=gcc-14 ./configure LAPACK_LIBS="-L/usr/local/lib -llapack -lblas"
(genesis) $ make
(genesis) $ make install

# Check the compiled binaries in `bin/`
(genesis) $ ls bin/
atdyn*                 drms_analysis*         kmeans_clustering*     pmf_analysis*          rmsd_analysis*
avecrd_analysis*       dssp_interface*        lipidthick_analysis*   prjcrd_analysis*       rpath_generator*
cg_convert*            eigmat_analysis*       mbar_analysis*         qmmm_generator*        rst_convert*
comcrd_analysis*       emmap_generator*       meanforce_analysis*    qval_analysis*         rst_upgrade*
contact_analysis*      energy_analysis*       morph_generator*       qval_residcg_analysis* sasa_analysis*
crd_convert*           flccrd_analysis*       msd_analysis*          rdf_analysis*          spdyn*
density_analysis*      fret_analysis*         pathcv_analysis*       remd_convert*          tilt_analysis*
diffusion_analysis*    hb_analysis*           pcavec_drawer*         rg_analysis*           trj_analysis*
distmat_analysis*      hbond_analysis*        pcrd_convert*          ring_analysis*         wham_analysis*

# Check the compiled libraries in `bin/`
(genesis) $ ls lib/
libpython_interface.la* libpython_interface.so*

# Set environment variables
(genesis) $ cd /path/to/genesis/
(genesis) $  pip install -e .

```

# Download chignolin data

```
(genesis) $ cd /path/to/genesis/demo/
(genesis) $ gdown --id 1WyFzvhuMjlwp2pNjga9B8RvTKoygBh-a -O chignolin.pdb
(genesis) $ gdown --id 1L1Y7YdSz46sTI1lQ7PoQJIqqbzM4F9Vh -O chignolin.psf
(genesis) $ gdown --id 1DZFUbCBVdCsfKzzrroIslre0eSctMaY- -O chignolin.dcd
```

# Jupyter notebook

```
(genesis) $ cd /path/to/genesis/demo/
(genesis) $ jupyter-lab
# Open demo.ipynb in JupyterLab and let's execte cells!
```

# Regression test

The Python tests live in `src/genepie/tests` and run with pytest:

```
(genesis) $ cd /path/to/genesis
(genesis) $ pytest -m "not slow"        # what CI runs
(genesis) $ pytest                      # everything
```

Test for an individual analysis tool can be run, for example

```
(genesis) $ pytest src/genepie/tests/test_rmsd.py
```

## Python script editing

For example, in the case of rmsd_analysis.py
* Write the path to PDB/PSF files for generating s_molecule: pdb_path / psf_path
* Write the keywords (from inp file) for crd_convert execution to generate s_trajectory in the arguments of crd_convert
* Write the keywords (from inp file) for analysis in the arguments of trj_analysis

```
import os
import pathlib
from ctrl_files import TrajectoryParameters
from s_molecule import SMolecule
import genesis_exe


def test_rmsd_analysis():
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")

    mol = SMolecule.from_file(pdb=pdb_path, psf=psf_path)
    with genesis_exe.crd_convert(
            mol,
            traj_params = [
                TrajectoryParameters(
                    trjfile = "BPTI_run.dcd",
                    md_step = 10,
                    mdout_period = 1,
                    ana_period = 1,
                    repeat = 1,
                    ),
                ],
            trj_format = "DCD",
            trj_type = "COOR+BOX",
            trj_natom = 0,
            selection_group = ["all", ],
            fitting_method = "NO",
            fitting_atom = 1,
            check_only = False,
            pbc_correct = "NO",
            ) as trajs:
        for t in trajs:
            d = genesis_exe.rmsd_analysis(
                    mol, t,
                    selection_group = ["sid:BPTI and an:CA", ],
                    fitting_method = "TR+ROT",
                    fitting_atom = 1,
                    check_only = False,
                    analysis_atom  = 1,
                    )
            print(d.rmsd, flush=True)


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_rmsd_analysis()


if __name__ == "__main__":
    main()
```

The `trj_analysis.py` script implements regression tests using the CustomTestCase class, so the keywords for crd_convert are written in `custom_test_case.py`.

## Adding a new analysis wrapper (library-mode error guard)

GENESIS reports fatal errors with `error_msg`, which calls `exit(1)`. Inside
Python that would kill the interpreter, so every `bind(C)` entry point runs its
work under the error guard in `error_mod` (`run_guarded`, setjmp/longjmp based).
The body that the guard calls must be a **module procedure** that receives its
data through a context variable. Never use an internal (`contains`) procedure
as the callback: taking `c_funloc()` of it makes gfortran generate a
trampoline, which requires an executable stack, which `dlopen()` refuses on
glibc >= 2.41 (and which crashes on non-executable stacks). gfortran's
`-Wtrampolines` reports any such callback.

Skeleton for a new tool `foo` (see `rmsd_analysis_lazy_c` in `rmsd_c_mod.fpp`
for a complete, tested example):

```fortran
  ! 1. Context: everything the body reads, writes or acquires.
  type :: t_foo_ctx
    type(c_ptr) :: coords_ptr = c_null_ptr   ! inputs (copies of the C arguments)
    integer     :: natom = 0
    integer     :: nframe = 0
    type(c_ptr) :: result_ptr = c_null_ptr
    integer     :: nstru = 0                 ! outputs (copied back by the wrapper)
    type(s_error)       :: err               ! error state
    type(s_trj_source)  :: source            ! resources released by the wrapper
    type(s_result_sink) :: sink
  end type t_foo_ctx

  ! 2. Entry point: pack -> run_guarded -> unpack -> release.
  subroutine foo_analysis_c(coords_ptr, natom, nframe, result_ptr, nstru, &
                            status, msg, msglen) bind(C, name="foo_analysis_c")
    type(c_ptr),    value       :: coords_ptr, result_ptr
    integer(c_int), value       :: natom, nframe, msglen
    integer(c_int), intent(out) :: nstru, status
    character(kind=c_char), intent(out) :: msg(*)
    type(t_foo_ctx), target :: c

    c%coords_ptr = coords_ptr
    c%natom      = natom
    c%nframe     = nframe
    c%result_ptr = result_ptr

    call run_guarded(foo_analysis_body, c, c%err, status, msg, msglen)

    nstru = c%nstru
    call finalize_sink(c%sink)
    call finalize_source(c%source)
  end subroutine foo_analysis_c

  ! 3. Body: a bind(C) module procedure; failures go into c%err only.
  subroutine foo_analysis_body(ctx) bind(C, name="foo_analysis_body")
    type(c_ptr), value :: ctx
    type(t_foo_ctx), pointer :: c
    real(wp), pointer :: coords(:,:,:), result(:)

    call c_f_pointer(ctx, c)
    if (c%natom <= 0) then
      call error_set(c%err, ERROR_INVALID_PARAM, "foo_analysis_c: natom must be positive")
      return
    end if
    call c_f_pointer(c%coords_ptr, coords, [3, c%natom, c%nframe])
    call c_f_pointer(c%result_ptr, result, [c%nframe])
    ! ... init_source_*(c%source, ...), init_sink_array(c%sink, result, c%nframe),
    !     then call the shared analysis core. If it calls error_msg, the guard
    !     turns that into c%err and the wrapper reports it to Python.
  end subroutine foo_analysis_body
```

Rules of thumb:

- Give the body and the context type a name unique to the module: `bind(C)`
  names are global symbols of the shared library.
- Inputs are plain copies. Reference arguments (`character(*)` strings,
  `s_molecule_c`) are either converted first, as `c_filename_to_fortran` does
  for file names, or stored as `c_ptr` via `c_loc` (add `target` to the dummy).
- Allocatable components of the context are freed automatically when the entry
  point returns. Only units/files need explicit `finalize_*` calls, and those
  are safe to call even when the body never initialised the object.
- Never call `error_to_c` from the body: `run_guarded` reports `c%err`.
- The interface is not thread-safe, and the guard is not either.

# Folder and file description

## Folders

All files exist in `genesis/src/analysis/interface/python_interface`.

Files with changed module names related to mbar_analysis created by `mbar_rename.py` exist in `genesis/src/analysis/interface/mbar_analysis` (content description omitted as only module names were changed)

## Files

### s_molecule structure related

|File name            |Content                                      |
|:--------------------|:--------------------------------------------|
|s_molecule.py        |Python s_molecule class definition, Python s_molecule/s_molecule_c mutual conversion, Python s_molecule object creation/deallocation from input files|
|s_molecule_c.py      |s_molecule_c class definition                |
|s_molecule_c_mod.fpp |s_molecule_c/s_molecule mutual conversion, s_molecule_c memory allocation/deallocation|
|define_molecule.fpp  |s_molecule_c creation from input file       |

### s_trajectories structure related

|File name                |Content                                        |
|:------------------------|:----------------------------------------------|
|s_trajectories.py        |Python s_trajectory class definition          |
|s_trajectories_c.py      |s_trajectory_c class definition               |
|s_trajectories_c_mod.fpp |s_trajectory_c initialization/copy/frame acquisition/memory allocation/deallocation etc.|

### Common program related

|File name                |Content                                        |
|:------------------------|:----------------------------------------------|
|c2py_util.py             |C data -> numpy NDArray conversion            |
|py2c_util.py             |numpy NDArray -> C data conversion            |
|conv_f_c_util.fpp        |C data <-> FORTRAN data conversion            |
|genesis_exe.py           |Call FORTRAN functions from Python            |
|libgenesis.py            |C type definition of FORTRAN subroutines      |
|mbar_rename.py           |mbar_analysis related module name replacement file generation tool|
|ctrl_c_mod.fpp           |iso_c_binding version control data            |
|ctrl_files.py            |Temporary control file output                  |

### crd_convert related

|File name                |Content                                        |
|:------------------------|:----------------------------------------------|
|crd_convert.py           |crd_convert regression test execution program |
|crd_convert_c_mod.fpp    |GENESIS cc_main.fpp and cc_setup.fpp functionality|
|crd_convert_convert.fpp  |GENESIS cc_convert.fpp functionality          |

### trj_analysis related

|File name                 |Content                                       |
|:-------------------------|:---------------------------------------------|
|trj_analysis.py           |trj_analysis regression test execution program (Distance/Angle/Dihedral)|
|trj_analysis_c_mod.fpp    |GENESIS ta_main.fpp and ta_setup.fpp functionality|
|trj_analysis_analysis.fpp |GENESIS ta_analyze.fpp functionality         |

### wham_analysis related

|File name                 |Content                                       |
|:-------------------------|:---------------------------------------------|
|wham_analysis.py          |wham_analysis regression test execution program|
|wa_analysis_c_mod.fpp     |GENESIS wa_main.fpp and wa_setup.fpp functionality|
|wa_analysis_analysis.fpp  |GENESIS wa_analyze.fpp functionality         |

### mbar_analysis related

|File name                        |Content                                                    |
|:-------------------------------|:------------------------------------------------------|
|mbar_analysis.py                |mbar_analysis regression test execution program        |
|mbar_analysis_umbrella_1d.py    |mbar_analysis regression test execution program (Umbrella 1D)|
|mbar_analysis_umbrella_block.py |mbar_analysis regression test execution program (Umbrella Block)|
|mbar_analysis_c_mod.fpp         |GENESIS ma_main.fpp and ma_setup.fpp functionality    |
|mbar_analysis_analysis.fpp      |GENESIS ma_analyze.fpp functionality                  |

* After executing mbar_rename.py, files with changed module names related to mbar_analysis exist in genesis/src/analysis/interface/mbar_analysis

### pmf_analysis related

|File name        |Content                                                                 |
|:----------------|:-----------------------------------------------------------------------|
|pmf_c_mod.fpp    |bind(C) wrapper (`pmf_analysis_c`); receives the control string, calls `analyze_pmf_unified` of the CLI module `pm_analyze_mod`, and returns the PMF array pointer + dimensions|

* Reuses `control_from_string` in `../../free_energy/pmf_analysis/pm_control.fpp` and links the static `libpmf_analysis.a` built from the pmf_analysis CLI sources.

### avecrd_analysis related

|File name                 |Content                                       |
|:-------------------------|:---------------------------------------------|
|avecrd_analysis.py        |avecrd_analysis regression test execution program|
|aa_analysis_c_mod.fpp     |GENESIS aa_main.fpp and aa_setup.fpp functionality|
|aa_analysis_analysis.fpp  |GENESIS aa_analyze.fpp functionality         |

### kmeans_clustering related

|File name                 |Content                                       |
|:-------------------------|:---------------------------------------------|
|kmeans_clustering.py      |kmeans_clustering regression test execution program|
|kc_analysis_c_mod.fpp     |GENESIS kc_main.fpp and kc_setup.fpp functionality|
|kc_analysis_analysis.fpp  |GENESIS kc_analyze.fpp functionality         |

### hb_analysis related

|File name                 |Content                                       |
|:-------------------------|:---------------------------------------------|
|hb_analysis_count_atom.py |hb_analysis regression test execution program|
|hb_analysis_count_snap.py |hb_analysis regression test execution program|
|hb_analysis_c_mod.fpp     |GENESIS hb_main.fpp and hb_setup.fpp functionality|
|hb_analysis_analysis.fpp  |GENESIS hb_analyze.fpp functionality         |

### rmsd_analysis related

|File name                 |Content                                       |
|:-------------------------|:---------------------------------------------|
|rmsd_analysis.py          |rmsd_analysis regression test execution program|
|ra_analysis_c_mod.fpp     |GENESIS ra_main.fpp and ra_setup.fpp functionality|
|ra_analysis_analysis.fpp  |GENESIS ra_analyze.fpp functionality         |

### drms_analysis related

|File name                 |Content                                       |
|:-------------------------|:---------------------------------------------|
|drms_analysis.py          |drms_analysis regression test execution program|
|dr_analysis_c_mod.fpp     |GENESIS dr_main.fpp and dr_setup.fpp functionality|
|dr_analysis_analysis.fpp  |GENESIS dr_analyze.fpp functionality         |

### rg_analysis related

|File name                 |Content                                       |
|:-------------------------|:---------------------------------------------|
|rg_analysis.py            |rg_analysis regression test execution program|
|rg_analysis_c_mod.fpp     |GENESIS rg_main.fpp and rg_setup.fpp functionality|
|rg_analysis_analysis.fpp  |GENESIS rg_analyze.fpp functionality         |

### msd_analysis related

|File name                   |Content                                       |
|:---------------------------|:---------------------------------------------|
|msd_analysis.py             |msd_analysis regression test execution program|
|ma_analysis_c_mod.fpp       |GENESIS ma_main.fpp and ma_setup.fpp functionality|
|ma_analysis_analysis.fpp    |GENESIS ma_analyze.fpp functionality         |

### diffusion_analysis related

|File name                        |Content                                       |
|:---------------------------------|:---------------------------------------------|
|diffusion_analysis.py             |diffusion_analysis regression test execution program|
|diffusion_analysis_main_c_mod.fpp |GENESIS da_main.fpp and da_setup.fpp functionality|
|diffusion_analysis_analyze.fpp    |GENESIS da_analyze.fpp functionality         |

### MDTraj, MDAnalysis related

|File name                        |Content                              |
|:---------------------------------|:------------------------------------|
|test_mdanalysis.py                |MdAnalysis regression test execution program|
|test_mdtraj.py                    |MdTraj regression test execution program|
