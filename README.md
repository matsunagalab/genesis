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

The suite is a plain pytest suite (`uv pip install -e ".[dev]"` installs
pytest). Tests find the freshly built library in the build tree, so no
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

Unified tools: RMSD (`ra_analyze.fpp`), RG (`rg_analyze.fpp`), DRMS (`dr_analyze.fpp`).

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

**Recommended: Unified Architecture** (for RMSD, RG, DRMS pattern)

When CLI already exists with `trj_source_mod` support:
1. Export `analyze_*_unified()` from CLI's `*_analyze.fpp` with primitive arguments
2. Create `*_c_mod.fpp` that calls unified function via `init_source_memory()` + `init_sink_array()`
3. Add function signature to `libgenesis.py` and wrapper to `genesis_exe.py`
4. Add test as `tests/test_<name>.py`

**Alternative: Separate Implementation** (for HB, WHAM, MBAR pattern)

When unified pattern doesn't fit:
1. Create `*_c_mod.fpp` with `bind(C)` interface
2. Create `*_impl.fpp` for analysis implementation
3. Update `Makefile.am` to include new `.fpp` files
4. Add function signature to `libgenesis.py` and wrapper to `genesis_exe.py`
5. Add test as `tests/test_<name>.py`

See [CLAUDE.md](CLAUDE.md) for detailed instructions.

---

## Documentation

- **[Tutorials & docs site](docs/)** — runnable Jupyter Book (see [Tutorials](#tutorials) above)
- [Getting Started](docs/getting-started.ipynb) — a complete workflow in ~10 lines
- [GENESIS Website](https://www.r-ccs.riken.jp/labs/cbrt/)
- [CLAUDE.md](CLAUDE.md) - Developer guide

## License

LGPL-3.0-or-later. See LICENSE file for details.
