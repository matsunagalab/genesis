# Installation

genepie ships as a platform-specific wheel that bundles the compiled GENESIS
binaries (the `atdyn` MD engine, 43 analysis CLI tools, and the
`libpython_interface` shared library used by the Python API).

## Requirements

- Python 3.9+
- Linux (x86_64) or macOS (arm64 / x86_64)
- Linux additionally needs glibc 2.28+ (the wheels are `manylinux_2_28`)

| Ubuntu version | glibc | Status |
|----------------|-------|--------|
| 24.04 LTS | 2.39 | Supported |
| 22.04 LTS | 2.35 | Supported |
| 20.04 LTS | 2.31 | Supported |
| 18.04 LTS | 2.27 | Not supported (build from source) |

## Install from PyPI

```bash
# Create a virtual environment (uv shown here; venv/conda also work)
uv venv --python=python3.11
source .venv/bin/activate

# From PyPI (coming soon)
uv pip install genepie

# Currently available from TestPyPI:
uv pip install --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ genepie
```

## Verify the installation

The basic analysis tests run on the small trajectories bundled inside the
package, so they need no extra downloads:

```bash
uv pip install pytest
pytest --pyargs genepie.tests.test_rmsd genepie.tests.test_crd_convert \
       genepie.tests.test_rg genepie.tests.test_drms genepie.tests.test_avecrd
```

If pytest reports every test as passed, your installation is working.

## Install from source (developers)

Building from source is required to run the regression tests (WHAM/MBAR, the MD
engine, etc.) and to work on the Fortran layer.

```bash
git clone https://github.com/matsunagalab/genepie.git
cd genepie

uv venv --python=python3.11
source .venv/bin/activate
uv pip install numpy

# Build dependencies:
#   Linux:  sudo apt install gfortran liblapack-dev libblas-dev autoconf automake libtool
#   macOS:  brew install gcc lapack autoconf automake libtool

autoreconf -fi

# Linux:
./configure --disable-mpi CC=gcc FC=gfortran LAPACK_LIBS="-llapack -lblas"
# macOS:
./configure --disable-mpi CC=gcc-14 FC=gfortran \
    LAPACK_LIBS="-L$(brew --prefix lapack)/lib -llapack -lblas"

make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)
make install
uv pip install -e .
```

```{admonition} Locating the shared library in a source build
:class: note
`make install` copies `libpython_interface` next to the `genepie` package so it
is found automatically. If you build without installing, point genepie at the
libtool output directory instead:

    export GENEPIE_LIB_DIR="$PWD/src/analysis/interface/python_interface/.libs"

This is the same variable the documentation build uses locally.
```

## Building this documentation

```bash
uv pip install -r docs/requirements.txt
uv pip install -e .                 # genepie itself must be importable
cd docs && jupyter-book build .
open _build/html/index.html         # macOS (use xdg-open on Linux)
```

The notebooks are executed during the build, so a working genepie installation
(with the shared library discoverable) is required.
