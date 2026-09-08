# GENESIS

**GENeralized-Ensemble SImulation System** - A molecular dynamics simulation software for biomolecular systems.

## Python Interface (genepie)

The `genepie` package provides a Python interface to GENESIS analysis tools and the ATDYN MD engine.

### Why genepie?

An MD study is not a single command — it is a workflow of *setup → simulation → analysis*, with many
analysis methods to choose from, and increasingly a final step that feeds the results into AI/ML tooling.
Doing this by stitching together CLI tools and intermediate files is tedious and hard to reproduce.

`genepie` turns GENESIS into a **single, programmable Python workflow**: run simulations, run analyses, and
hand the numerical results straight to NumPy / PyTorch / scikit-learn without leaving Python. It bridges the
GENESIS Fortran engine with researchers and the AI/ML ecosystem.

### Architecture

`genepie` is a thin interface layer over the existing GENESIS Fortran engine — it does not reimplement any
science. Three layers are connected by standard interop mechanisms:

```
┌────────────────────────────────────────────────────┐
│ User layer       Python script / Jupyter / NumPy   │
├────────────────────────────────────────────────────┤   ctypes  (C ↔ Python)
│ Interface layer  genepie (this package)            │
├────────────────────────────────────────────────────┤   ISO_C_BINDING (Fortran ↔ C)
│ Engine layer     GENESIS (Fortran: MD & analysis)  │
└────────────────────────────────────────────────────┘
```

The interface layer is built around four design ideas (see [Design Highlights](#design-highlights) for details):

- **Zerocopy memory** — NumPy arrays and Fortran share the same memory (no copy).
- **Lazy DCD loading** — read one frame at a time for large trajectories.
- **CLI/Python unified core** — the same Fortran routine serves both the CLI and Python (identical results).
- **Structured error handling** — a typed exception hierarchy with numeric error codes.

---

## Table of Contents

- [Python Interface (genepie)](#python-interface-genepie)
  - [Why genepie?](#why-genepie)
  - [Architecture](#architecture)
- [For Users](#for-users)
  - [Installation](#installation)
  - [Testing Your Installation](#testing-your-installation)
  - [Quick Start](#quick-start)
  - [Available Analysis Functions](#available-analysis-functions)
  - [MD Engine Functions](#md-engine-functions)
  - [Supported File Formats](#supported-file-formats)
  - [Tutorials](#tutorials)
- [For Developers](#for-developers)
  - [Installation from Source](#installation-from-source)
  - [Running Tests](#running-tests)
  - [Developer Workflow](#developer-workflow)
  - [Project Structure](#project-structure)
  - [Design Highlights](#design-highlights)
  - [Adding a New Analysis Tool](#adding-a-new-analysis-tool)
- [Documentation](#documentation)
- [License](#license)

---

## For Users

### Installation

```bash
# Create virtual environment with uv
uv venv --python=python3.11
source .venv/bin/activate

# From PyPI (coming soon)
uv pip install genepie

# Currently available from TestPyPI:
uv pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ genepie
```

**Requirements:**
- Python 3.9+
- Linux (x86_64) or macOS (arm64, x86_64)
- glibc 2.28+ for Linux (see table below)

| Ubuntu Version | glibc | Status |
|---------------|-------|--------|
| 24.04 LTS | 2.39 | Supported |
| 22.04 LTS | 2.35 | Supported |
| 20.04 LTS | 2.31 | Supported |
| 18.04 LTS | 2.27 | Not supported (build from source) |

### Testing Your Installation

```bash
uv pip install pytest

# Basic analysis tests (no additional data required)
pytest --pyargs genepie.tests.test_rmsd genepie.tests.test_crd_convert \
       genepie.tests.test_rg genepie.tests.test_drms genepie.tests.test_avecrd

# Integration tests (requires ~500 MB download)
uv pip install gdown mdtraj MDAnalysis
python -m genepie.tests.download_test_data
pytest --pyargs genepie.tests.test_integration
```

Tests that compare against the CLI reference data (test_trj, test_wham,
test_pmf, test_mbar, test_atdyn) need the full source repository and skip
themselves otherwise.

### Quick Start

```python
from genepie import genesis_exe, SMolecule

# Load molecular structure (pdb/psf topology, ref = reference structure for RMSD)
mol = SMolecule.from_file(pdb="protein.pdb", psf="protein.psf", ref="protein.pdb")
print(f"Loaded {mol.num_atoms} atoms")

# Load trajectory. crd_convert returns (list_of_trajectories, subset_molecule)
trajs, subset_mol = genesis_exe.crd_convert(
    mol,
    trj_files=["trajectory.dcd"],
    trj_format="DCD",
    trj_type="COOR+BOX",
    selection="all",
)

# Calculate RMSD of C-alpha atoms with translational + rotational fitting
result = genesis_exe.rmsd_analysis(
    mol,
    trajs[0],
    analysis_selection="an:CA",
    fitting_selection="an:CA",
    fitting_method="TR+ROT",
)
print(f"RMSD: {result.rmsd.mean():.2f} Å")  # result.rmsd is a NumPy array

# Run MD simulation (returns a namedtuple: energies, final_coords, energy_labels)
md = genesis_exe.run_atdyn_md(
    prmtopfile="protein.prmtop",
    ambcrdfile="protein.inpcrd",
    nsteps=1000,
    ensemble="NVT",
    temperature=300.0,
)
print(md.energies.shape, md.final_coords.shape)
```

For memory-efficient processing of large trajectories, pass `lazy=True` to `crd_convert()` and feed the
resulting lazy trajectory to `rmsd_analysis()` — the same call reads frames from disk on demand.

### Available Analysis Functions

`Zerocopy` = results/coordinates share memory with NumPy (no copy).
`Lazy` = supports frame-by-frame DCD loading for large trajectories.

| Function | Description | Zerocopy | Lazy |
|----------|-------------|:--------:|:----:|
| `crd_convert()` | Coordinate/trajectory conversion | ✓ | – |
| `trj_analysis()` | Distance, angle, dihedral analysis | ✓ | ✓ |
| `rmsd_analysis()` | RMSD calculation | ✓ | ✓ |
| `drms_analysis()` | Distance RMSD calculation | ✓ | ✓ |
| `rg_analysis()` | Radius of gyration | ✓ | ✓ |
| `diffusion_analysis()` | Diffusion coefficient calculation | ✓ | – |
| `msd_analysis()` | Mean squared displacement | – | – |
| `hb_analysis()` | Hydrogen bond analysis | – | – |
| `avecrd_analysis()` | Average coordinate calculation | – | – |
| `wham_analysis()` | WHAM free energy analysis | – | – |
| `mbar_analysis()` | MBAR free energy analysis | – | – |
| `pmf_analysis()` | PMF from CV samples + optional weights (histogram / Gaussian kernel) | – | – |
| `mbar_resample_trajectory()` | MBAR-weighted trajectory resampling at a target temperature or bias-free condition | – | – |
| `kmeans_clustering()` | K-means trajectory clustering | – | – |

The full GENESIS analysis suite (43 tools) is also installed as CLI commands; the table above lists the
tools currently wrapped as native Python functions.

`wham_analysis`, `mbar_analysis`, and `pmf_analysis` each also provide an
`*_isolated(..., timeout=...)` variant that runs the solve in a throwaway subprocess, giving every call
clean Fortran state so many estimates can be computed back to back in one session. `pmf_analysis` accepts
either in-memory NumPy arrays (`cv`, `weight`) or CLI-style filename patterns (`cvfile`, `weightfile`).

### MD Engine Functions

- `run_atdyn_md()` - Run MD simulation
- `run_atdyn_min()` - Run energy minimization
- `run_atdyn_md_isolated()` - Run MD in subprocess (crash-safe)
- `run_atdyn_min_isolated()` - Run minimization in subprocess

### Supported File Formats

| Format | Topology | Coordinates | Parameters |
|--------|----------|-------------|------------|
| AMBER | `prmtopfile` | `ambcrdfile` | (in prmtop) |
| GROMACS | `grotopfile` | `grocrdfile` | (in grotop) |
| CHARMM | `psffile` | `pdbfile`/`crdfile` | `parfile`, `strfile` |

### Tutorials

The best way to learn genepie is the **hands-on tutorial series** in [`docs/`](docs/).
It is a [Jupyter Book](https://jupyterbook.org): every chapter is a runnable notebook,
and (except for the ML chapter) they all run on the small trajectories bundled with
the package — **no downloads required**. New to genepie? Start with
**[Getting Started](docs/getting-started.ipynb)**: a complete load → analyze → plot
workflow on the BPTI system in about ten lines.

| # | Chapter | What you will do |
|---|---------|------------------|
| — | [Getting Started](docs/getting-started.ipynb) | Full RMSD workflow on BPTI, end to end |
| 1 | [Molecules](docs/tutorials/01_molecule.ipynb) | Load and visualize a structure with `SMolecule` |
| 2 | [Trajectories](docs/tutorials/02_trajectory.ipynb) | Load DCDs, apply selections, fitting, centering |
| 2b | [Lazy Loading](docs/tutorials/02b_lazy_loading.ipynb) | Analyze trajectories too big for RAM |
| 3 | [Structural Analysis](docs/tutorials/03_structure.ipynb) | RMSD, Rg, DRMS, distances/angles, average structure |
| 4 | [Dynamics](docs/tutorials/04_dynamics.ipynb) | Mean squared displacement and diffusion |
| 5 | [Free Energy](docs/tutorials/05_free_energy.ipynb) | Recover a PMF with WHAM & MBAR |
| 5b | [MBAR Weighted Resampling](docs/tutorials/05b_mbar_resampling.ipynb) | Reweight T-REMD to 300 K, resample frames, and save a uniform-weight DCD |
| 6 | [Clustering](docs/tutorials/06_clustering.ipynb) | k-means conformational clustering (GENESIS + scikit-learn) |
| 7 | [MD Engine](docs/tutorials/07_md_engine.ipynb) | Run minimization and MD with ATDYN |
| 8 | [ML Integration](docs/tutorials/08_ml_integration.ipynb) | Feed results into MDTraj, scikit-learn, PyTorch |

**Run them locally:**

```bash
# Open any chapter interactively
jupyter lab docs/getting-started.ipynb

# ...or build the whole book into a browsable HTML site
uv pip install -r docs/requirements.txt
jupyter-book build docs
open docs/_build/html/index.html          # macOS (use xdg-open on Linux)
```

#### MBAR Weighted Resampling (T-REMD)

[`docs/tutorials/05b_mbar_resampling.ipynb`](docs/tutorials/05b_mbar_resampling.ipynb) demonstrates
the complete T-REMD workflow: calculate unbiased per-sample MBAR weights at a
target temperature, verify the reweighted Ramachandran PMF, resample frames with
replacement, and save the resulting uniform-weight ensemble as a DCD file.

The tutorial dataset is downloaded on demand:

```bash
uv pip install gdown mdtraj matplotlib
python -m genepie.tests.download_tremd_data
jupyter lab docs/tutorials/05b_mbar_resampling.ipynb
```

---

## For Developers

> **genepie is an open-source public repository.** External contributors do **not** have push access
> to `matsunagalab/genepie`. To contribute, **fork** the repository, push your changes to a branch on
> your fork, and open a **pull request** back to `matsunagalab/genepie`. See
> [Contributing (Fork & Pull Request)](#4-contributing-fork--pull-request) below.

### Installation from Source

```bash
# 1. Fork matsunagalab/genepie on GitHub (click "Fork"), then clone YOUR fork
git clone https://github.com/<your-username>/genepie.git
cd genepie

# 2. Add the upstream repository so you can keep your fork in sync
git remote add upstream https://github.com/matsunagalab/genepie.git

# Set up Python environment with uv
uv venv --python=python3.11
source .venv/bin/activate
uv pip install numpy

# Install build dependencies
# Linux (Ubuntu/Debian):
#   sudo apt install gfortran liblapack-dev libblas-dev autoconf automake libtool
# macOS:
#   brew install gcc lapack autoconf automake libtool

# Build GENESIS
autoreconf -fi

# Linux:
./configure --disable-mpi CC=gcc FC=gfortran LAPACK_LIBS="-llapack -lblas"

# macOS:
./configure --disable-mpi CC=gcc-14 FC=gfortran \
    LAPACK_LIBS="-L$(brew --prefix lapack)/lib -llapack -lblas"

make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)

# Install in editable mode
uv pip install -e .
```

### Running Tests

The test suite is plain **pytest**: every test is a `test_*` function in
`src/genepie/tests/test_<tool>.py` (no `unittest` classes, no custom base
class). Shared paths, fixtures and skip markers (`requires_regression_data`,
`requires_chignolin`, `requires_tremd`, ...) live in
`src/genepie/tests/conftest.py`; `testpaths` and the `slow` marker are
configured in `pyproject.toml`. `uv pip install -e ".[dev]"` installs pytest,
and tests find the freshly built library in the build tree, so no
`make install` is needed.

```bash
# Everything except the memory-heavy tests (what CI runs)
pytest -m "not slow"

# One module, one test
pytest src/genepie/tests/test_rmsd.py
pytest src/genepie/tests/test_rmsd.py -k lazy

# Everything, including test_msd / test_diffusion
pytest

# Optional data; the tests that need it skip themselves until it is present
uv pip install gdown mdtraj MDAnalysis
python -m genepie.tests.download_test_data    # chignolin -> test_integration
python -m genepie.tests.download_tremd_data   # T-REMD   -> test_mbar_resample
```

Tests that compare against the CLI reference values (test_trj, test_wham,
test_pmf, test_mbar, test_atdyn) read `tests/regression_test/` and therefore
need the full source repository.

### Developer Workflow

#### 1. Local Development Cycle

| Step | Action | Notes |
|------|--------|-------|
| 1 | Edit code | Python: instant; Fortran: requires `make` |
| 2 | Run tests | See testing strategy below |
| 3 | Commit changes | Use descriptive commit messages |
| 4 | Push to your fork | Create a feature branch, push to `origin` (your fork) |
| 5 | Open a PR | From your fork's branch → `matsunagalab/genepie` `main` |

#### 2. Rebuilding After Changes

| Change Type | Required Action |
|-------------|-----------------|
| Python files (.py) | No rebuild needed |
| Fortran files (.fpp) | `make` |
| New Fortran files | `make clean && make` |
| configure.ac changes | `autoreconf -fi && ./configure ... && make` |

#### 3. Testing Strategy

| Change Type | Tests to Run |
|-------------|--------------|
| Python API changes | Basic tests + `test_integration` |
| Fortran interface | All tests including regression tests |
| Bug fixes | Relevant test + add new regression test |
| New analysis function | Create `src/genepie/tests/test_<name>.py` (pytest picks it up) |

#### 4. Contributing (Fork & Pull Request)

`genepie` is a public repository, so contributions go through the standard GitHub **fork & pull request**
flow. You do not need write access to the upstream repository.

```bash
# 1. Fork matsunagalab/genepie on GitHub, then clone your fork
git clone https://github.com/<your-username>/genepie.git
cd genepie
git remote add upstream https://github.com/matsunagalab/genepie.git

# 2. Sync your fork with upstream before starting
git checkout main
git fetch upstream
git merge upstream/main
git push origin main

# 3. Create a feature branch
git checkout -b feature/my-change

# 4. Make changes, run tests, and commit
git add -A
git commit -m "Describe your change"

# 5. Push the branch to YOUR fork
git push origin feature/my-change

# 6. Open a pull request on GitHub:
#    base: matsunagalab/genepie  main  ←  compare: <your-username>/genepie  feature/my-change
```

Then:

1. GitHub Actions automatically builds wheels for all platforms on the PR
2. Review CI build results (Linux x86_64, macOS arm64, macOS x86_64)
3. Address review feedback by pushing more commits to the same branch
4. A maintainer merges after approval and CI passes

#### 5. Release Workflow

> **Maintainers only.** The steps below require push access to `matsunagalab/genepie` and are performed
> by maintainers after a contributor's pull request has been merged.

```
Fork → PR → main → TestPyPI (optional) → Tag → PyPI
```

**Step 1: Test on TestPyPI (optional but recommended)**

```bash
# Trigger manually via GitHub Actions UI:
# Actions → "Build and publish wheels" → Run workflow (workflow_dispatch)

# Then test the package:
pip install --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ genepie
python -c "import genepie; print(genepie.__version__)"
```

**Step 2: Release to PyPI**

```bash
# 1. Update version in pyproject.toml
# 2. Commit the version bump
git add pyproject.toml
git commit -m "Bump version to 0.1.x"
git push origin main

# 3. Create and push version tag
git tag v0.1.x
git push origin v0.1.x
# GitHub Actions automatically builds and publishes to PyPI
```

#### 6. CI/CD Pipeline Overview

| Workflow | Trigger | Action |
|----------|---------|--------|
| `tests.yml` | Push/PR to main | Run test suite on Linux & macOS |
| `build-wheels.yml` | Push/PR to main | Build wheels for all platforms |
| `build-wheels.yml` | `workflow_dispatch` | Build + publish to TestPyPI |
| `build-wheels.yml` | Push tag `v*` | Build + publish to PyPI |

**Test Platforms (`tests.yml`):**
- Linux x86_64 (manylinux_2_28 container)
- macOS arm64 (Apple Silicon)

**Build Platforms (`build-wheels.yml`):**
- Linux x86_64 (manylinux_2_28)
- macOS arm64 (Apple Silicon)
- macOS x86_64 (Intel)

### Project Structure

```
genepie/
├── src/
│   ├── genepie/           # Python interface (main package)
│   │   ├── genesis_exe.py # Analysis function wrappers
│   │   ├── libloader.py   # Shared library loader
│   │   └── tests/         # Test files and data
│   ├── atdyn/             # MD engine
│   └── analysis/          # Analysis tools
├── CLAUDE.md              # Developer guide for Claude Code
└── pyproject.toml         # Package configuration
```

### Design Highlights

These four ideas define how the interface layer connects Python to the GENESIS Fortran engine.

#### 1. Zerocopy memory management

Coordinate arrays and analysis results are shared between Python and Fortran instead of being copied.
Python (NumPy) owns and allocates the memory; Fortran creates an alias to it with `C_F_POINTER` and writes
results directly into the NumPy array.

```
Copy-based (before):  Python → copy → Fortran → copy → Python   (2× memory + overhead)
Zerocopy (after):     Python numpy array  ⇄  Fortran C_F_POINTER (1× memory, in-place)
```

Applied to: `crd_convert`, `rmsd_analysis`, `rg_analysis`, `drms_analysis`, `trj_analysis`, `diffusion_analysis`.

#### 2. Lazy DCD loading

Instead of loading an entire trajectory into memory, lazy mode reads one frame at a time directly from disk.
Because DCD files have a fixed frame size, any frame is reachable in O(1) via a computed byte offset
(`byte_offset = header_size + (N - 1) * frame_size`), using Fortran stream I/O. This keeps memory usage
minimal for large trajectories. Enabled via `crd_convert(..., lazy=True)` and currently supported by
`rmsd_analysis`, `rg_analysis`, `drms_analysis`, and `trj_analysis` (atom- and COM-based measurements).

#### 3. CLI / Python unified core

The CLI tool and the Python interface call **the same** analysis routine (`analyze_*_unified()`), so results
are guaranteed identical and bugs are fixed in one place. The routine is agnostic to its caller; only the I/O
is abstracted:

- `trj_source` — where frames are read from: a file on disk (CLI) or a NumPy array in memory (Python).
  Implemented as `TRJ_SOURCE_FILE`, `TRJ_SOURCE_MEMORY`, and `TRJ_SOURCE_LAZY_DCD` in `trj_source_mod`.
- `result_sink` — where results are written to: an output file (CLI) or a zerocopy NumPy array (Python).

Every tool with a Python API follows this pattern (RMSD, RG, DRMS, trj, avecrd, HB, MSD, k-means,
diffusion, WHAM, MBAR, PMF); see [Adding a New Analysis Tool](#adding-a-new-analysis-tool).

#### 4. Structured error handling

Fortran return codes are mapped to a typed Python exception hierarchy, and messages written to Fortran's
stderr are captured and attached to the exception. Numeric codes (`e.code`) allow programmatic handling.

| Code range | Category | Exception |
|-----------|----------|-----------|
| 100–199 | Memory (alloc/dealloc) | `GenesisFortranMemoryError` |
| 200–299 | File I/O (not found, format) | `GenesisFortranFileError` |
| 300–399 | Validation (invalid/missing params) | `GenesisFortranValidationError` |
| 400–499 | Data (mismatch, no data) | `GenesisFortranDataError` |
| 500–599 | Not supported (feature, dimension) | `GenesisFortranNotSupportedError` |
| 600–699 | Internal (internal, syntax) | `GenesisFortranInternalError` |

All inherit from `GenesisFortranError` → `GenesisError`. See [CLAUDE.md](CLAUDE.md) for usage examples.

### Adding a New Analysis Tool

There is one pattern: **the science lives in the CLI module, and the Python interface only marshals
arguments**. Every tool that has both a command-line program and a Python API shares a single implementation
(`analyze_<tool>_unified`), so CLI and Python results are identical by construction and there is exactly one
place to fix a bug. The former `*_impl.fpp` copies of CLI code in `src/analysis/interface/python_interface/`
are gone; do not add new ones.

#### Step 1. Give the CLI module a unified core

In `src/analysis/<family>/<tool>_analysis/*_analyze.fpp`, split `analyze()` into two routines:

- `analyze_<tool>_unified(...)` — the computation. It reads frames through `trj_source_mod`
  (`s_trj_source`, `get_next_frame`, `reset_source` when the trajectory is read more than once) and returns
  its results through `result_sink_mod` sinks or through array arguments. It never opens a file by name and
  knows nothing about the CLI control file or NumPy.
- `analyze(...)` — the CLI driver. It keeps the old signature, builds the inputs from `option`/`input`,
  calls `init_source_file(...)`, opens file sinks (`init_sink_file`, `init_sink_file_rows`,
  `init_sink_text`), calls the core and writes whatever the core returned in memory.

`result_sink_mod` has three flavours, each with a file and an in-memory form: one scalar per frame
(`write_result`), one row of values per frame (`write_result_row`) and text lines (`write_result_line`).
Free-energy tools (WHAM, MBAR, PMF) and file-based tools (diffusion) have no trajectory source; their core
simply returns arrays and the driver writes them.

Examples: `ra_analyze.fpp` (one scalar per frame), `ta_analyze.fpp` (rows per frame, six streams),
`aa_analyze.fpp` (scalar sink + a result written back into `molecule`), `hb_analyze.fpp` (text lines),
`ma_analyze.fpp` / `kc_analyze.fpp` (arrays, several passes over the trajectory),
`wa_analyze.fpp` / `mbar_analyze.fpp` / `pm_analyze.fpp` (arrays).

Check that the CLI output is unchanged with the regression harness:

```bash
cd tests/regression_test/test_analysis
python test_analysis.py /path/to/src/analysis/<family>/<tool>_analysis -- Test_<tool>_analysis
```

#### Step 2. Write the `bind(C)` wrapper (`<tool>_c_mod.fpp`)

GENESIS reports fatal errors with `error_msg`, which calls `exit(1)`. Inside Python that would kill the
interpreter, so every `bind(C)` entry point runs its work under the library-mode error guard in `error_mod`
(`run_guarded`, setjmp/longjmp based): an `error_msg` inside the analysis core becomes a Python exception
instead of a process exit. The body that the guard calls must be a **module procedure** that receives its
data through a context variable. Never use an internal (`contains`) procedure as the callback: taking
`c_funloc()` of it makes gfortran generate a trampoline, which needs an executable stack, which `dlopen()`
refuses on glibc >= 2.41. The build fails with `-Werror=trampolines` on purpose if this happens.

Declare every `bind(C)` entry point that Python looks up **`public`** in its module, even though the
binding label alone gives it external linkage: gfortran 16.2.0 gives `private` `bind(C)` procedures hidden
visibility (GCC PR fortran/126872), and `dlsym()` then cannot find them.

Skeleton for a new tool `foo` in `src/analysis/interface/python_interface/foo_c_mod.fpp`
(`rmsd_c_mod.fpp` is a complete, tested example; `hb_c_mod.fpp` shows a control-text based tool with a
text sink):

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
    !     then call analyze_foo_unified(...). If it calls error_msg, the guard
    !     turns that into c%err and the wrapper reports it to Python.
  end subroutine foo_analysis_body
```

Rules of thumb:

- Give the body and the context type a name unique to the module: `bind(C)` names are global symbols of the
  shared library.
- Inputs are plain copies. Reference arguments (`character(*)` strings, `s_molecule_c`) are either converted
  first, as `c_filename_to_fortran` does for file names, or stored as `c_ptr` via `c_loc` (add `target` to
  the dummy argument).
- In-memory trajectories: `init_source_memory(...)`. Lazy DCD input: `init_source_lazy_dcd(...)`, with the
  Python wrapper dispatching on `trajs.is_lazy` (see `rmsd.py`). NumPy result buffers: `init_sink_array` /
  `init_sink_array_rows` on views made with `c_f_pointer`.
- Allocatable components of the context are freed automatically when the entry point returns. Only
  units/files need explicit `finalize_*` calls, and those are safe to call even when the body never
  initialised the object.
- Results handed to Python (`c_ptr` to a C string or array) are freed by Python only
  (`deallocate_c_string`, `deallocate_double2`). Never keep them in a Fortran `save` pointer that a later
  call frees again.
- Never call `error_to_c` from the body: `run_guarded` reports `c%err`.
- The interface is not thread-safe, and the guard is not either.

Register the file in `src/analysis/interface/python_interface/Makefile.am` and `Makefile.depends`, then
`make clean && make`.

#### Step 3. Python side

1. Add the exact C signature (`argtypes` / `restype`) to `src/genepie/libgenesis.py`.
2. Add the wrapper in `src/genepie/analysis/<tool>.py` and re-export it from
   `src/genepie/analysis/__init__.py` (which flows through to `genesis_exe`). Conventions:
   - Convert the molecule once with `molecule.to_SMoleculeC()` and free it in `finally:` with
     `deallocate_s_molecule_c`.
   - Pre-allocate result arrays in NumPy (`float64`, C-contiguous) and pass `.ctypes.data_as(c_void_p)`
     for zerocopy output. Python `(nframe, natom, 3)` C-order aliases Fortran `(3, natom, nframe)`.
   - Wrap the call in `with fortran_status() as (status, msg, msglen):` so Fortran errors surface as the
     typed exceptions above.
   - Return a small `namedtuple` so results are self-describing.
3. Add `src/genepie/tests/test_<tool>.py` with plain `test_*` functions; pytest collects it automatically
   (see [Running Tests](#running-tests)). Compare against the CLI reference data in
   `tests/regression_test/` when it exists (`test_trj.py` is an example).

---

## Documentation

- **[Tutorials & docs site](docs/)** — runnable Jupyter Book (see [Tutorials](#tutorials) above)
- [Getting Started](docs/getting-started.ipynb) — a complete workflow in ~10 lines
- [GENESIS Website](https://www.r-ccs.riken.jp/labs/cbrt/)
- [CLAUDE.md](CLAUDE.md) - Developer guide

## License

LGPL-3.0-or-later. See LICENSE file for details.
