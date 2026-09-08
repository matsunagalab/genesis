# SPDX-License-Identifier: LGPL-3.0-or-later
# Minimal, robust loader for libgenesis.* without LD_LIBRARY_PATH.
from __future__ import annotations

import ctypes
import os
import platform
from pathlib import Path
from importlib.resources import files
from typing import Iterable, List, Optional, Tuple

CANDIDATES: tuple[str, ...] = (
        "libgenesis.so",
        "libpython_interface.so",
        "libgenesis.dylib",
        "genesis.dll")


def _get_glibc_version() -> Optional[Tuple[int, int]]:
    """Get glibc version on Linux systems.

    Uses multiple detection methods for robustness:
    1. platform.libc_ver() - standard library method
    2. ctypes CDLL + gnu_get_libc_version() - direct C library call
    3. /lib libc.so.6 symlink parsing - fallback for edge cases

    Returns:
        Tuple of (major, minor) version numbers, or None if not on Linux/glibc.
    """
    if platform.system() != 'Linux':
        return None

    # Method 1: platform.libc_ver()
    try:
        libc_name, version = platform.libc_ver()
        if libc_name.lower() == 'glibc' and version:
            parts = version.split('.')
            if len(parts) >= 2:
                return (int(parts[0]), int(parts[1]))
    except Exception:
        pass

    # Method 2: ctypes direct call to gnu_get_libc_version()
    try:
        import ctypes.util
        libc_path = ctypes.util.find_library('c')
        if libc_path:
            libc = ctypes.CDLL(libc_path, mode=ctypes.RTLD_LOCAL)
            gnu_get_libc_version = libc.gnu_get_libc_version
            gnu_get_libc_version.restype = ctypes.c_char_p
            version_bytes = gnu_get_libc_version()
            if version_bytes:
                version_str = version_bytes.decode('utf-8')
                parts = version_str.split('.')
                if len(parts) >= 2:
                    return (int(parts[0]), int(parts[1]))
    except Exception:
        pass

    # Method 3: Parse /lib/libc.so.6 or similar symlink
    try:
        import re
        for lib_dir in ['/lib', '/lib64', '/lib/x86_64-linux-gnu',
                        '/lib/aarch64-linux-gnu']:  # ARM64 support
            libc_path = Path(lib_dir) / 'libc.so.6'
            if libc_path.exists():
                # Try to read the symlink target or library itself
                target = libc_path.resolve().name if libc_path.is_symlink() else libc_path.name
                # Match patterns like "libc-2.31.so", "libc.so.6", or "libc-2.31-0ubuntu9.so"
                patterns = [
                    r'libc[.-](\d+)\.(\d+)',  # libc-2.31.so
                    r'(\d+)\.(\d+)',           # Version numbers anywhere
                ]
                for pattern in patterns:
                    match = re.search(pattern, target)
                    if match:
                        major, minor = int(match.group(1)), int(match.group(2))
                        # Sanity check: glibc versions are 2.x or potentially 3.x
                        if 2 <= major <= 3:
                            return (major, minor)
    except Exception:
        pass

    # Method 4: Parse ldd --version output (last resort)
    try:
        import re
        import subprocess
        result = subprocess.run(['ldd', '--version'],
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            # Output example: "ldd (Ubuntu GLIBC 2.31-0ubuntu9.7) 2.31"
            first_line = result.stdout.split('\n')[0]
            match = re.search(r'(\d+)\.(\d+)\s*$', first_line)
            if match:
                return (int(match.group(1)), int(match.group(2)))
    except Exception:
        pass

    return None


def _check_glibc_version(min_major: int = 2, min_minor: int = 28) -> None:
    """Check if glibc version meets minimum requirements.

    Raises:
        OSError: If glibc version is too old.
    """
    version = _get_glibc_version()
    if version is None:
        return  # Not Linux or can't detect - let ctypes handle it
    major, minor = version
    if (major, minor) < (min_major, min_minor):
        raise OSError(
            f"genepie requires glibc {min_major}.{min_minor}+ "
            f"(found glibc {major}.{minor}).\n"
            f"Ubuntu 20.04+ is required. Ubuntu 18.04 is not supported.\n"
            f"Options:\n"
            f"  1. Upgrade to Ubuntu 20.04 or later\n"
            f"  2. Build genepie from source"
        )


def _validate_platform_compatibility() -> None:
    """Validate platform is compatible with GENESIS.

    Checks:
    1. Endianness: GENESIS requires little-endian (x86_64, ARM64)
    2. Pointer size: GENESIS requires 64-bit pointers

    Raises:
        OSError: If platform is not supported
    """
    import sys
    import struct

    # Check endianness
    if sys.byteorder != 'little':
        raise OSError(
            f"GENESIS requires a little-endian platform. "
            f"Detected: {sys.byteorder}-endian. "
            f"Big-endian systems (SPARC, some IBM Power) are not supported."
        )

    # Check pointer size (must be 64-bit)
    pointer_size = struct.calcsize('P') * 8
    if pointer_size != 64:
        raise OSError(
            f"GENESIS requires a 64-bit platform. "
            f"Detected: {pointer_size}-bit pointers. "
            f"32-bit systems are not supported."
        )


def _first_existing(paths: Iterable[str]) -> Optional[str]:
    for p in paths:
        if os.path.exists(p):
            return p
    return None

def _candidate_paths_in_dir(d: Path) -> List[Path]:
    return [d / name for name in CANDIDATES]

def _find_in_package_dotlib(pkg: str) -> Optional[Path]:
    try:
        libdir = Path(files(pkg).joinpath(".lib"))
    except Exception:
        return None
    return _first_existing(_candidate_paths_in_dir(libdir))


def _find_repo_root_lib_near(mod_file: str) -> Optional[Path]:
    here = Path(mod_file).resolve()
    for parent in list(here.parents)[:8]:
        cand_dir = parent / "lib"
        if cand_dir.is_dir():
            hit = _first_existing(_candidate_paths_in_dir(cand_dir))
            if hit:
                return hit
    return None

def _find_build_tree_lib(mod_file: str) -> Optional[Path]:
    """Library freshly built in a source checkout (before ``make install``)."""
    here = Path(mod_file).resolve()
    for parent in list(here.parents)[:8]:
        cand_dir = (parent / "src" / "analysis" / "interface"
                    / "python_interface" / ".libs")
        if cand_dir.is_dir():
            hit = _first_existing(_candidate_paths_in_dir(cand_dir))
            if hit:
                return hit
    return None

def _find_pkg_lib_python_interface_dotlib(mod_file: str) -> Optional[Path]:
    here = Path(mod_file).resolve()
    for parent in list(here.parents)[:8]:
        cand_dir = parent / "pkg" / "lib" / "python_interface" / ".lib"
        if cand_dir.is_dir():
            hit = _first_existing(_candidate_paths_in_dir(cand_dir))
            if hit:
                return hit
    return None

def _diagnose_load_error(path: Path, error: OSError) -> str:
    """Diagnose library load failure and provide helpful suggestions.

    Args:
        path: Path to the library that failed to load
        error: The OSError that was raised

    Returns:
        Diagnostic message with suggestions
    """
    error_msg = str(error)
    suggestions = []

    # Check for libmvec symbol errors (legacy wheels or -ffast-math builds)
    if '_ZGVdN4v_' in error_msg or '_ZGVbN' in error_msg or 'libmvec' in error_msg:
        suggestions.append(
            "libmvec symbol error detected. This may occur with older wheel versions.\n"
            "  Workaround:\n"
            "    ARCH_DIR=$(gcc -print-multiarch 2>/dev/null || echo \"x86_64-linux-gnu\")\n"
            "    export LD_PRELOAD=/lib/${ARCH_DIR}/libmvec.so.1:/lib/${ARCH_DIR}/libm.so.6\n"
            "    python your_script.py\n"
            "  Or upgrade to the latest genepie version."
        )

    # Check for libgfortran errors
    if 'libgfortran' in error_msg:
        suggestions.append(
            "libgfortran not found. Install gfortran runtime:\n"
            "  Ubuntu/Debian: sudo apt install libgfortran5\n"
            "  RHEL/CentOS:   sudo yum install libgfortran\n"
            "  Fedora:        sudo dnf install libgfortran"
        )

    # Check for OpenMP runtime errors (libgomp)
    if 'libgomp' in error_msg or 'GOMP_' in error_msg or 'omp_' in error_msg.lower():
        suggestions.append(
            "OpenMP runtime library (libgomp) not found.\n"
            "  This is required for parallel execution in GENESIS.\n"
            "  Install OpenMP runtime:\n"
            "    Ubuntu/Debian: sudo apt install libgomp1\n"
            "    RHEL/CentOS:   sudo yum install libgomp\n"
            "    Fedora:        sudo dnf install libgomp\n"
            "    macOS:         brew install gcc  # includes OpenMP support"
        )

    # Check for LAPACK/BLAS errors
    if 'lapack' in error_msg.lower() or 'blas' in error_msg.lower():
        suggestions.append(
            "LAPACK/BLAS not found. Install linear algebra libraries:\n"
            "  Ubuntu/Debian: sudo apt install liblapack3 libblas3\n"
            "  RHEL/CentOS:   sudo yum install lapack blas\n"
            "  Fedora:        sudo dnf install lapack blas"
        )

    # Check for general undefined symbol errors
    if 'undefined symbol' in error_msg:
        suggestions.append(
            "Undefined symbol error. This may indicate:\n"
            "  - Missing system library dependency\n"
            "  - glibc version incompatibility\n"
            "  - Library built with incompatible compiler"
        )

    # Check for GLIBC version errors
    if 'GLIBC_' in error_msg:
        import re
        match = re.search(r'GLIBC_(\d+\.\d+)', error_msg)
        if match:
            required = match.group(1)
            suggestions.append(
                f"glibc {required} or newer is required.\n"
                "  Check your glibc version: ldd --version\n"
                "  Ubuntu 20.04+ provides glibc 2.31+\n"
                "  Options:\n"
                "    1. Upgrade your Linux distribution\n"
                "    2. Build genepie from source"
            )

    if suggestions:
        return "\n\n".join(suggestions)
    else:
        return (
            "Library load failed. Possible causes:\n"
            "  - Missing runtime dependencies\n"
            "  - Architecture mismatch (x86_64 vs arm64)\n"
            "  - Corrupted library file\n"
            "Check dependencies: ldd " + str(path)
        )


def _load(path: Path) -> ctypes.CDLL:
    """Load shared library with detailed error diagnostics."""
    try:
        return ctypes.CDLL(os.fspath(path), mode=ctypes.RTLD_GLOBAL)
    except OSError as e:
        diagnostic = _diagnose_load_error(path, e)
        raise OSError(
            f"Failed to load {path}:\n"
            f"  Original error: {e}\n\n"
            f"Diagnosis:\n{diagnostic}"
        ) from e


def _find_in_same_dir(mod_file: str) -> Optional[Path]:
    """Search in the same directory as this module (for meson-python builds)."""
    here = Path(mod_file).resolve().parent
    return _first_existing(_candidate_paths_in_dir(here))


def _find_in_installed_genepie() -> Optional[Path]:
    """Search in the installed genepie package location.

    Searches sys.path for genepie packages that contain .so files,
    skipping the local source directory.
    """
    import sys
    local_dir = Path(__file__).resolve().parent
    for path_entry in sys.path:
        try:
            candidate_dir = Path(path_entry) / "genepie"
            if not candidate_dir.is_dir():
                continue
            # Skip local source directory
            if candidate_dir.resolve() == local_dir:
                continue
            hit = _first_existing(_candidate_paths_in_dir(candidate_dir))
            if hit:
                return hit
        except Exception:
            continue
    return None


def load_genesis_lib() -> ctypes.CDLL:
    # Validate platform compatibility (endianness, 64-bit)
    _validate_platform_compatibility()

    # Check glibc version first on Linux (manylinux_2_28 wheels require glibc 2.28+)
    _check_glibc_version()

    tried: list[str] = []

    env_dir = os.environ.get("GENEPIE_LIB_DIR")
    if env_dir:
        d = Path(env_dir)
        tried.append(f"{d}/(env)")
        hit = _first_existing(_candidate_paths_in_dir(d))
        if hit:
            return _load(hit)

    # in same directory as this module (meson-python builds)
    p = _find_in_same_dir(__file__)
    tried.append("<package_dir>/")
    if p:
        return _load(p)

    # in installed genepie package (for running local tests with installed package)
    p = _find_in_installed_genepie()
    tried.append("<installed_genepie>/")
    if p:
        return _load(p)

    # in python_interface
    p = _find_in_package_dotlib("python_interface")
    tried.append("python_interface/.lib/")
    if p:
        return _load(p)

    # in genepie/.lib
    p = _find_in_package_dotlib("genepie")
    tried.append("genepie/.lib/")
    if p:
        return _load(p)

    # in <repo>/lib
    p = _find_repo_root_lib_near(__file__)
    tried.append("<repo>/lib/")
    if p:
        return _load(p)

    # in <repo>/pkg/lib/python_interface/.lib/
    p = _find_pkg_lib_python_interface_dotlib(__file__)
    tried.append("<repo>/pkg/lib/python_interface/.lib/")
    if p:
        return _load(p)

    # in the build tree of a source checkout (after make, before make install)
    p = _find_build_tree_lib(__file__)
    tried.append("<repo>/src/analysis/interface/python_interface/.libs/")
    if p:
        return _load(p)

    raise OSError(
        "genesis shared library not found.\nSearched:\n  - "
        + "\n  - ".join(tried)
        + "\nHint: set GENEPIE_LIB_DIR=/absolute/path/to/dir containing one of: "
        + ", ".join(CANDIDATES)
    )

