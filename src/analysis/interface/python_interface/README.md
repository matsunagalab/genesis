# python_interface: the Fortran side of the Python interface

This directory builds `libpython_interface`, the shared library that Python
loads through `ctypes`. It contains only the `bind(C)` wrappers and the small
runtime they share. The analysis code itself lives in the CLI modules under
`src/analysis/` (and `src/atdyn/` for the MD engine) and is linked in as
static libraries; see `Makefile.am`.

Everything user- and developer-facing is documented in the repository README:

- Build, install and tests:
  [README.md, For Developers](../../../../README.md#for-developers)
- How to add a tool, the wrapper skeleton (context type + `run_guarded` +
  `bind(C)` body) and the rules of thumb:
  [README.md, Adding a New Analysis Tool](../../../../README.md#adding-a-new-analysis-tool)
- Python side: `src/genepie/` (ctypes signatures in `libgenesis.py`, wrappers
  in `analysis/`, pytest suite in `tests/`)

When you add or remove a `.fpp` file here, update `Makefile.am` and
`Makefile.depends` (`make depend` regenerates the latter).

## Shared runtime

| File | Content |
|------|---------|
| `error_mod.fpp` | Library-mode error guard: `s_error`, `error_set`, `run_guarded`, status code / message hand-over to Python (`fi_msg_len`) |
| `conv_f_c_util.fpp` | C/Fortran conversion helpers (`c_filename_to_fortran`, `cptr_to_fstring`) and the `allocate_c_*` / `deallocate_*` routines Python uses to free Fortran-allocated results |
| `ctrl_c_mod.fpp` | C structs for control-file sections and their conversion to GENESIS `s_inp_info` / `s_out_info` / `s_trj_info` / selection info |
| `s_molecule_c_mod.fpp` | `s_molecule_c` (C mirror of `s_molecule`), `allocate_s_molecule_c` / `deallocate_s_molecule_c` |
| `define_molecule.fpp` | `define_molecule_from_file`: build a molecule from PDB/PSF/PRMTOP/GRO... for Python |
| `export_molecule_c_mod.fpp` | `export_pdb_to_string_c`: PDB text of a molecule returned to Python |
| `s_trajectories_c_mod.fpp` | `s_trajectories_c` buffers shared with NumPy: init, deep copy, join, arrays of trajectories, and their deallocation (respect the ownership rules) |
| `selection_c_mod.fpp` | `selection_c`: atom selection strings evaluated from Python |
| `dynamic_string_mod.fpp`, `internal_file_type_mod.fpp` | Growing strings and PDB-to-string output used by wrappers that return text (avecrd, kmeans, export) |

## Tool wrappers

One `<tool>_c_mod.fpp` per Python API. Each packs its arguments into a context
type, runs a module-procedure body under `run_guarded`, and calls the unified
core of the CLI module listed here (`analyze_<tool>_unified`, or the shared
conversion routines for `crd_convert`).

| Wrapper | CLI module (under `src/analysis/`) |
|---------|------------------------------------|
| `crd_convert_c_mod.fpp` | `converter/crd_convert` (`cc_convert.fpp`) |
| `trj_c_mod.fpp` | `trj_analysis/trj_analysis` (`ta_analyze.fpp`) |
| `rmsd_c_mod.fpp` | `trj_analysis/rmsd_analysis` (`ra_analyze.fpp`), in-memory and lazy DCD |
| `rg_c_mod.fpp` | `trj_analysis/rg_analysis` (`rg_analyze.fpp`) |
| `drms_c_mod.fpp` | `trj_analysis/drms_analysis` (`dr_analyze.fpp`) |
| `msd_c_mod.fpp` | `trj_analysis/msd_analysis` (`ma_analyze.fpp`) |
| `hb_c_mod.fpp` | `trj_analysis/hb_analysis` (`hb_analyze.fpp`), text sink |
| `diffusion_c_mod.fpp` | `trj_analysis/diffusion_analysis` (`da_analyze.fpp`) |
| `avecrd_c_mod.fpp` | `mode_analysis/avecrd_analysis` (`aa_analyze.fpp`) |
| `kmeans_c_mod.fpp` | `clustering/kmeans_clustering` (`kc_analyze.fpp`) |
| `wham_c_mod.fpp` | `free_energy/wham_analysis` (`wa_analyze.fpp`) |
| `mbar_c_mod.fpp` | `free_energy/mbar_analysis` (`mbar_analyze.fpp`) |
| `pmf_c_mod.fpp` | `free_energy/pmf_analysis` (`pm_analyze.fpp`) |
| `atdyn_c_mod.fpp` | `src/atdyn` (`atdyn_md_c`, `atdyn_min_c`, `reset_atdyn_state_c`) |

`rmsd_c_mod.fpp` is the reference example for a trajectory tool and
`hb_c_mod.fpp` for a control-text based tool with a text sink.
