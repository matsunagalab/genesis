# Adding a New Analysis Tool

There is one pattern: **the science lives in the CLI module, and the Python
interface only marshals arguments**. Every tool that has both a command-line
program and a Python API shares a single implementation, so CLI and Python
results are identical by construction and there is exactly one place to fix a
bug. (The former `*_impl.fpp` copies of the CLI code in
`src/analysis/interface/python_interface/` are gone; do not add new ones.)

## 1. Give the CLI module a unified core

In the tool's `src/analysis/<family>/<tool>_analysis/*_analyze.fpp`, split
`analyze()` into

- `analyze_<tool>_unified(...)`: the computation. It reads frames through
  `trj_source_mod` (`s_trj_source`, `get_next_frame`, `reset_source` when the
  trajectory is read several times) and returns its results either through
  `result_sink_mod` sinks or through allocatable/array arguments. It never opens
  files by name; it knows nothing about the CLI control file or NumPy.
- `analyze(...)`: the CLI driver. It keeps the old signature, builds the inputs
  from `option`/`input`, calls `init_source_file(...)`, opens file sinks
  (`init_sink_file`, `init_sink_file_rows`, `init_sink_text`), calls the core and
  writes whatever the core returned in memory.

`result_sink_mod` offers three flavours, each with a file and an in-memory
form: one scalar per frame (`write_result`), one row of values per frame
(`write_result_row`) and text lines (`write_result_line`). Free-energy tools
(WHAM, MBAR, PMF) and file-based tools (diffusion) have no trajectory source;
their core simply returns arrays and the driver writes them.

Examples: `ta_analyze.fpp` (rows per frame, six streams), `aa_analyze.fpp`
(scalar sink + a result written back into `molecule`), `hb_analyze.fpp` (text
lines), `ma_analyze.fpp` / `kc_analyze.fpp` (arrays, several passes over the
trajectory), `wa_analyze.fpp` / `mbar_analyze.fpp` / `pm_analyze.fpp` (arrays).

Verify the CLI still produces identical output with the regression harness:

```bash
cd tests/regression_test/test_analysis
python test_analysis.py /path/to/src/analysis/<family>/<tool>_analysis -- Test_<tool>_analysis
```

## 2. Write the bind(C) entry point

Create `src/analysis/interface/python_interface/<tool>_c_mod.fpp`. The entry
point packs its arguments into a context type, runs a module-procedure body under
the library-mode error guard, copies the outputs back and releases resources.
The skeleton and the rules are in
`src/analysis/interface/python_interface/README.md` ("Adding a new analysis
wrapper"); `rmsd_c_mod.fpp` is a complete example, `hb_c_mod.fpp` shows a
control-text based tool with a text sink.

- In-memory trajectories: `init_source_memory(...)`; lazy DCD input:
  `init_source_lazy_dcd(...)` and let the Python wrapper dispatch on
  `trajs.is_lazy` (see `rmsd.py`).
- NumPy result buffers: `init_sink_array` / `init_sink_array_rows` on views made
  with `c_f_pointer`.
- Never take `c_funloc` of an internal procedure: the build fails with
  `-Werror=trampolines` on purpose (an executable stack is not loadable on
  glibc 2.41).

Add the new `.fpp` to `Makefile.am` and to `Makefile.depends` in that directory.

## 3. Python side

1. Add the C signature to `src/genepie/libgenesis.py` (`argtypes` / `restype`).
2. Add the wrapper in `src/genepie/analysis/<tool>.py` and re-export it from
   `genepie/analysis/__init__.py` (which flows through to `genesis_exe`).
3. Add `src/genepie/tests/test_<tool>.py` with plain `test_*` functions; pytest
   collects it automatically, shared paths and helpers live in `conftest.py`.
   Compare against the CLI reference data in `tests/regression_test/` when it
   exists (see `test_trj.py`).

## Python wrapper conventions

- Convert the molecule once: `mol_c = molecule.to_SMoleculeC()` and free it in a
  `finally:` with `deallocate_s_molecule_c`.
- Pre-allocate result arrays in NumPy and pass `.ctypes.data_as(c_void_p)` for
  zerocopy output.
- Wrap the call in `with fortran_status() as (status, msg, msglen):` so Fortran
  errors surface as the typed exceptions in [the error model](../reference/errors.md).
- Return a small `namedtuple` so results are self-describing.
- Resolve user-facing enums with the case-insensitive helpers in `_common.py`.

## Rebuild matrix

| Change | Command |
|--------|---------|
| Python only (`.py`) | none — instant |
| Fortran (`.fpp`) | `make` |
| New Fortran files | update `Makefile.am`, then `make clean && make` |
| `configure.ac` | `autoreconf -fi && ./configure ... && make` |

See the [CLAUDE.md](https://github.com/matsunagalab/genepie/blob/main/CLAUDE.md)
developer guide for the fully worked signatures.
