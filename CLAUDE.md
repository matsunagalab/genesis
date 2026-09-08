# GENESIS / genepie Repository Guide

## Scope

GENESIS is a Fortran molecular-dynamics and analysis codebase. This repository
also ships `genepie`, a Python 3.9+ interface to selected analysis routines and
the ATDYN engine.

Keep this file short and stable. For details, use the code and configuration as
the source of truth:

1. `pyproject.toml` — package metadata, dependencies, and installed CLI commands
2. `src/genepie/analysis/` — Python API implementations
3. `src/analysis/interface/python_interface/` — C-callable Fortran wrappers
4. `.github/workflows/` — supported CI, wheel, docs, and release procedures
5. `README.md` and `docs/` — user-facing behavior and tutorials

## Repository Map

- `src/genepie/`: Python package, ctypes bindings, validation, and tests
- `src/genepie/analysis/`: analysis wrappers grouped by tool family
- `src/genepie/genesis_exe.py`: compatibility re-export; do not add substantial
  implementations here
- `src/analysis/`: Fortran analysis programs and libraries
- `src/analysis/libana/`: shared trajectory-source/result-sink abstractions
- `src/analysis/interface/python_interface/`: Python-facing Fortran interface
- `src/atdyn/`: atom-decomposition MD engine
- `src/spdyn/`: MPI domain-decomposition engine; not included in PyPI wheels
- `src/lib/`: shared Fortran infrastructure
- `src/genepie/tests/`: executable test modules and packaged test data
- `tests/regression_test/`: source-only reference data
- `docs/`: Jupyter Book and runnable tutorials

## Architecture and Invariants

- Python calls `libpython_interface` through `ctypes`; Fortran exports entry
  points with `ISO_C_BINDING`/`bind(C)`.
- The Fortran library is not thread-safe. Do not run Fortran-backed APIs
  concurrently in threads; use process isolation when concurrency is required.
- Public Python implementations belong in `src/genepie/analysis/*.py`. Re-export
  public names from `analysis/__init__.py`; preserve the historical
  `genesis_exe.<name>` API when practical.
- ctypes signatures belong in `src/genepie/libgenesis.py`. Keep Python
  `argtypes` exactly aligned with the Fortran `bind(C)` signature.
- Prefer the shared CLI/Python analysis core when one exists. `trj_source_mod`
  abstracts file, memory, and lazy-DCD input; `result_sink_mod` abstracts file
  and NumPy-array output.
- NumPy trajectory/result buffers may be shared with Fortran. They must be
  `float64`, C-contiguous, correctly shaped, and alive for the whole call.
  Python `(nframe, natom, 3)` C-order aliases Fortran `(3, natom, nframe)`.
- `SMolecule.to_SMoleculeC()` is copy-based. Always deallocate its C/Fortran
  object in `finally`; use the helpers in `src/genepie/_fortran.py` where
  possible.
- Respect ownership in `STrajectories`: never free NumPy-owned buffers through
  Fortran or leak Fortran-owned buffers.
- `bind(C)` wrappers run their work through `run_guarded` in `error_mod` with
  a context derived type; the guarded body must be a module procedure. Never
  take `c_funloc` of an internal procedure: the trampoline needs an executable
  stack, which `dlopen` rejects on glibc 2.41+. See the wrapper skeleton in
  `src/analysis/interface/python_interface/README.md`.
- Validate paths, enums, dimensions, sizes, and contiguity before crossing the
  language boundary. Map Fortran status/message outputs through
  `exceptions.py`; do not replace typed errors with generic exceptions.
- Stateful ATDYN and free-energy calls have isolated subprocess variants.
  Prefer those for repeated runs or crash containment.
- Lazy trajectories contain metadata instead of all coordinates. New analysis
  code must either support lazy input explicitly or reject it clearly.

## Build and Environment

Typical Python setup:

```bash
uv venv --python=python3.11
source .venv/bin/activate
uv pip install -e ".[dev]"
```

Build the non-MPI library and tools:

```bash
autoreconf -fi

# Linux
./configure --disable-mpi CC=gcc FC=gfortran \
  LAPACK_LIBS="-llapack -lblas"

# macOS (Homebrew)
GCC="$(ls "$(brew --prefix)"/bin/gcc-* | sort -V | tail -1)"
./configure --disable-mpi CC="$GCC" FC=gfortran \
  LAPACK_LIBS="-L$(brew --prefix lapack)/lib -llapack -lblas"

make -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
uv pip install -e ".[dev]"
```

Rebuild rules:

- Python-only change: no native rebuild
- Existing `.fpp` change: `make`
- New/removed Fortran source: update the relevant `Makefile.am`, regenerate
  dependencies if required, then use a clean rebuild
- `configure.ac`/`Makefile.am` change: `autoreconf -fi`, configure, and rebuild

Edit Autotools source files, not generated `Makefile`, `Makefile.in`, or helper
scripts. Regenerate generated files only when the task requires it.

## Testing

Run the narrowest relevant test first:

```bash
python -m genepie.tests.test_rmsd
python -m genepie.tests.test_error_handling
python -m genepie.tests.test_atdyn
```

Run the local suite with fail-fast shell behavior:

```bash
bash -e src/genepie/tests/all_run.sh
```

`all_run.sh` normally continues after a failed command, so do not rely on a
plain invocation's final exit status. It skips integration tests when the
downloaded chignolin data is absent. Optional datasets:

```bash
python -m genepie.tests.download_test_data
python -m genepie.tests.download_tremd_data
```

Testing expectations:

- Python API change: focused test plus relevant integration coverage
- Fortran/interface change: focused test plus regression/error tests
- Bug fix: add a regression that fails before the fix
- New analysis: add `test_<name>.py` and include it in `all_run.sh`
- `test_msd` and `test_diffusion` are memory-heavy and intentionally omitted
  from CI

Tests are mostly directly executable modules, not a conventional pytest-only
suite. Match the style of the neighboring test file.

## Common Change Paths

For a Python analysis API change:

1. Update the appropriate `src/genepie/analysis/` module.
2. Update exports in `analysis/__init__.py` and compatibility exports if needed.
3. Add validation and typed error behavior.
4. Update focused tests and user docs for externally visible changes.

For a new Fortran-backed operation:

1. Reuse an existing Fortran analysis core instead of duplicating science code.
2. Add a `bind(C)` wrapper with explicit status/message outputs.
3. Register the source in the interface `Makefile.am`.
4. Add the exact ctypes signature in `libgenesis.py`.
5. Add a thin Python wrapper with deterministic cleanup.
6. Test normal, invalid, empty, and lazy/ownership cases as applicable.

## Packaging, CI, and Releases

- Wheels bundle `libpython_interface`, ATDYN, and CLI binaries listed in
  `pyproject.toml`; `setup.py` forces platform-specific `py3-none-*` wheels.
- Supported wheel targets are Linux x86_64 (`manylinux_2_28`) and macOS arm64
  and x86_64. Linux wheels require glibc 2.28+.
- `tests.yml` builds and tests Linux and macOS on pushes/PRs to `main`.
- `docs.yml` executes the Jupyter Book and deploys Pages from `main`.
- `build-wheels.yml` builds wheels on PRs, publishes to TestPyPI on manual
  dispatch, and publishes to PyPI for `v*` tags.
- `pyproject.toml` is the package-version source of truth. If changing a
  release version, check any separately exposed `__version__` value too.

## Working Principles

- Inspect nearby code before editing; follow existing ownership and error
  patterns.
- Keep changes scoped. Do not modify unrelated generated files or user work.
- Preserve public API compatibility unless the task explicitly permits a
  breaking change.
- Do not claim CLI/Python numerical parity without a regression comparison.
- Verify commands and API details against current files instead of expanding
  this guide with volatile inventories.
