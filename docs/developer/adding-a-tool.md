# Adding a New Analysis Tool

There is one pattern: **the science lives in the CLI module as
`analyze_<tool>_unified(...)`, and the Python interface only marshals
arguments** through a `bind(C)` wrapper that runs under the library-mode error
guard. The test is a `src/genepie/tests/test_<tool>.py` with plain pytest
`test_*` functions.

The step-by-step guide, the Fortran wrapper skeleton (context type +
`run_guarded` + `bind(C)` body), the Python wrapper conventions and the
rebuild matrix are kept in one place, the repository README:

[README.md → Adding a New Analysis Tool](https://github.com/matsunagalab/genepie/blob/main/README.md#adding-a-new-analysis-tool)

Complete examples: `rmsd_c_mod.fpp` (in-memory and lazy DCD trajectory,
scalar sink) and `hb_c_mod.fpp` (control-text based tool, text sink) in
`src/analysis/interface/python_interface/`. Errors raised inside Fortran
surface as the typed exceptions described in [the error model](../reference/errors.md).
