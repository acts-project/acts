#!/usr/bin/env python3

import os
import argparse
import json
import urllib.request
import urllib.error
import re
import subprocess
import hashlib
import tempfile
from pathlib import Path
from typing import Tuple, Dict, Optional, Callable, TypeVar, Any
import contextlib
import time
import functools
from contextvars import ContextVar

# Modify the default cache dir to use a temporary directory
DEFAULT_CACHE_SIZE_LIMIT = 1 * 1024 * 1024  # 1MB

T = TypeVar("T")
remaining_retries: ContextVar[int] = ContextVar("remaining_retries", default=0)


def retry_on_http_error(max_retries: int = 3, base_delay: float = 1.0):
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> T:
            for attempt in range(max_retries):
                remaining_retries.set(max_retries - attempt - 1)
                try:
                    return func(*args, **kwargs)
                except urllib.error.HTTPError as e:
                    if attempt < max_retries - 1:
                        delay = base_delay * (2**attempt)
                        print(
                            f"Got HTTP error {e.code}, retrying in {delay} seconds..."
                        )
                        time.sleep(delay)
                        continue
                    raise
                except urllib.error.URLError as e:
                    if attempt < max_retries - 1:
                        delay = base_delay * (2**attempt)
                        print(f"Got URL error {e}, retrying in {delay} seconds...")
                        time.sleep(delay)
                        continue
                    raise
            return func(*args, **kwargs)  # Final attempt

        return wrapper

    return decorator


def compute_cache_key(url: str) -> str:
    """Compute a cache key for a URL"""
    return hashlib.sha256(url.encode()).hexdigest()


def compute_cache_digest(cache_dir: Path) -> str:
    """Compute a digest of all cache files except digest.txt"""
    files = sorted(
        f
        for f in os.listdir(cache_dir)
        if (cache_dir / f).is_file() and f != "digest.txt"
    )

    digest = hashlib.sha256()
    for fname in files:
        fpath = cache_dir / fname
        digest.update(fname.encode())
        digest.update(str(fpath.stat().st_size).encode())
        digest.update(fpath.read_bytes())
    return digest.hexdigest()


def update_cache_digest(cache_dir: Path):
    """Update the cache digest file"""
    digest = compute_cache_digest(cache_dir)
    (cache_dir / "digest.txt").write_text(digest)


def prune_cache(cache_dir: Optional[Path], size_limit: int):
    """Prune the cache to keep it under the size limit"""
    if cache_dir is None or not cache_dir.exists():
        return

    # Get all cache files with their modification times
    cache_files = [
        (cache_dir / f, (cache_dir / f).stat().st_mtime)
        for f in os.listdir(cache_dir)
        if (cache_dir / f).is_file()
        and f != "digest.txt"  # Exclude digest from pruning
    ]
    total_size = sum(f.stat().st_size for f, _ in cache_files)

    if total_size <= size_limit:
        return

    # Sort by modification time (oldest first)
    cache_files.sort(key=lambda x: x[1])

    # Remove files until we're under the limit
    for file_path, _ in cache_files:
        if total_size <= size_limit:
            break
        total_size -= file_path.stat().st_size
        file_path.unlink()

    # Update digest after pruning
    update_cache_digest(cache_dir)


@retry_on_http_error()
def fetch_github(base_url: str, cache_dir: Optional[Path], cache_limit: int) -> bytes:
    headers = {}
    token = os.environ.get("GITHUB_TOKEN")

    # Always authenticate when a token is available. Anonymous requests are
    # subject to GitHub's much lower unauthenticated rate limit, so dropping the
    # token on any attempt makes rate-limit failures *more* likely, not less.
    if token:
        headers["Authorization"] = f"Bearer {token}"

    print(f"Remaining retries: {remaining_retries.get()} for {base_url}")
    # Avoid printing the header values, which would leak the token into CI logs.
    print(f"Authenticated: {'Authorization' in headers}")

    with contextlib.ExitStack() as stack:
        if cache_dir is not None:
            cache_dir.mkdir(parents=True, exist_ok=True)
        else:
            cache_dir = Path(stack.enter_context(tempfile.TemporaryDirectory()))

        # Check cache first
        cache_key = compute_cache_key(base_url)
        cache_file = cache_dir / cache_key

        if cache_file.exists():
            print("Cache hit on", base_url)
            return cache_file.read_bytes()
        else:
            print("Cache miss on", base_url)

        try:
            req = urllib.request.Request(base_url, headers=headers)
            with urllib.request.urlopen(req) as response:
                content = response.read()

                # Write to cache
                cache_file.write_bytes(content)

                # Update digest after adding new file
                update_cache_digest(cache_dir)

                # Prune cache if necessary (this will update digest again if pruning occurs)
                prune_cache(cache_dir, cache_limit)

                return content
        except urllib.error.URLError as e:
            print(f"Failed to fetch from {base_url}: {e}")
            raise e
        except json.JSONDecodeError as e:
            print(f"Failed to parse JSON response: {e}")
            raise e


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", type=str, required=True, help="Tag to use")
    parser.add_argument("--arch", type=str, required=True, help="Architecture to use")
    parser.add_argument(
        "--compiler-binary",
        type=str,
        default=os.environ.get("CXX"),
        help="Compiler to use (defaults to CXX environment variable if set)",
    )
    parser.add_argument(
        "--compiler",
        type=str,
        default=None,
        help="Compiler to use (defaults to compiler binary if set)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file to write lockfile to",
    )
    parser.add_argument(
        "--cache-dir",
        type=lambda x: Path(x).expanduser() if x else None,
        default=os.environ.get("LOCKFILE_CACHE_DIR"),
        help="Directory to use for caching (defaults to LOCKFILE_CACHE_DIR env var)",
    )
    parser.add_argument(
        "--cache-limit",
        type=int,
        default=int(os.environ.get("LOCKFILE_CACHE_LIMIT", DEFAULT_CACHE_SIZE_LIMIT)),
        help="Cache size limit in bytes (defaults to LOCKFILE_CACHE_LIMIT env var)",
    )
    parser.add_argument(
        "--cxx",
        type=str,
        default=os.environ.get("CXXSTD", "20"),
        help="C++ standard (e.g. 20, 23). Defaults to CXXSTD env var or 20.",
    )
    parser.add_argument(
        "--flavor",
        type=str,
        default=os.environ.get("FLAVOR"),
        help=(
            "Accelerator flavor (e.g. cuda13, rocm-gfx90a). 'host' or unset "
            "selects the plain CPU stack. Defaults to the FLAVOR env var."
        ),
    )
    args = parser.parse_args()

    # Normalize to cxxNN format used in lockfile names
    cxx = args.cxx.strip()
    if not re.match(r"^\d+$", cxx):
        print(f"Invalid C++ standard: {args.cxx} (expected a number, e.g. 20 or 23)")
        exit(1)
    args.cxx = f"cxx{cxx}"

    # ci-dependencies spells the plain CPU stack as the flavor `host`, which
    # produces a triplet with no flavor token at all. Normalise both spellings
    # to None so everything downstream compares against the parsed token.
    args.flavor = normalize_flavor(args.flavor)

    print("Fetching lockfiles for tag:", args.tag)
    print("Architecture:", args.arch)
    print("C++ standard variant:", args.cxx)
    print("Flavor:", args.flavor or "host")

    base_url = f"https://api.github.com/repos/acts-project/ci-dependencies/releases/tags/{args.tag}"

    data = json.loads(fetch_github(base_url, args.cache_dir, args.cache_limit))

    lockfiles = parse_assets(data)

    print("Available lockfiles:")
    for arch, compilers in lockfiles.items():
        print(f"> {arch}:")
        for c, (n, _) in compilers.items():
            print(f"  - {c} [flavor: {extract_flavor(c) or 'host'}]: {n}")

    if args.arch not in lockfiles:
        print(f"No lockfile found for architecture {args.arch}")
        exit(1)

    if args.compiler_binary is not None:
        compiler = determine_compiler_version(args.compiler_binary)
        print("Compiler:", args.compiler_binary, f"{compiler}")
    elif args.compiler is not None:
        if not re.match(r"^([\w-]+)@(\d+\.\d+\.\d+)$", args.compiler):
            print(f"Invalid compiler format: {args.compiler}")
            exit(1)
        compiler = args.compiler
        print("Compiler:", f"{compiler}")
    else:
        compiler = None

    lockfile = select_lockfile(lockfiles, args.arch, compiler, args.cxx, args.flavor)

    print("Selected lockfile:", lockfile)

    if args.output:
        with open(args.output, "wb") as f:
            f.write(fetch_github(lockfile, args.cache_dir, args.cache_limit))


def parse_assets(data: Dict) -> Dict[str, Dict[str, Tuple[str, str]]]:
    lockfiles: Dict[str, Dict[str, Tuple[str, str]]] = {}

    for asset in data["assets"]:
        url = asset["browser_download_url"]
        name = asset["name"]

        if not name.startswith("spack_"):
            continue

        # New format: spack_<arch>_<compiler>_<cxx>-locks[.ext]
        # e.g. spack_linux-ubuntu24.04-aarch64_gcc@13.3.0_cxx20-locks
        # e.g. spack_linux-ubuntu24.04-x86_64_llvm@22.1.1_cxx23-locks.tar.gz
        m_new = re.match(
            r"spack_(.*(?:aarch64|x86_64))_(.+)-locks(?:\.[a-zA-Z0-9]+)*$",
            name,
        )
        if m_new:
            arch, compiler_spec = m_new.groups()
            lockfiles.setdefault(arch, {})[compiler_spec] = (name, url)
            continue

        # Legacy format: spack_<arch>[_(compiler)].lock
        # e.g. spack_linux-ubuntu24.04-x86_64.lock
        # e.g. spack_linux-ubuntu24.04-x86_64_gcc@13.3.0.lock
        if name.endswith(".lock"):
            m = re.match(r"spack_(.*(?:aarch64|x86_64))(?:_(.*))?\.lock", name)
            if m:
                arch, compiler = m.groups()
                compiler = compiler if compiler else "default"
                lockfiles.setdefault(arch, {})[compiler] = (name, url)

    return lockfiles


def normalize_flavor(flavor: Optional[str]) -> Optional[str]:
    """`host`, "" and None all mean the plain CPU stack, which has no token."""
    flavor = flavor.strip() if flavor else ""
    return None if flavor in ("", "host") else flavor


def extract_flavor(spec: str) -> Optional[str]:
    """Flavor token of a compiler spec, or None for the host stack.

    Lockfiles are named `<arch>_<compiler>_cxx<std>[_<flavor>]`
    (ci-dependencies flavors/README.md) and `parse_assets` strips the arch, so
    the flavor is whatever follows `cxx<std>` -- rejoined, since a flavor name
    may contain underscores. Pre-cxx specs (`gcc@13.3.0`, `default`) have none.
    """
    parts = spec.split("_")
    for i, part in enumerate(parts):
        if re.fullmatch(r"cxx\d+", part):
            rest = "_".join(parts[i + 1 :])
            return rest or None
    return None


def extract_cxx(spec: str) -> Optional[str]:
    """`cxx<std>` token of a compiler spec, or None.

    Not a suffix test: on a flavored spec it sits in the middle
    (`gcc@15.2.0_cxx23_cuda13`), where `endswith("_cxx23")` silently misses.
    """
    for part in spec.split("_"):
        if re.fullmatch(r"cxx\d+", part):
            return part
    return None


def select_lockfile(
    lockfiles: Dict[str, Dict[str, Tuple[str, str]]],
    arch: str,
    compiler: Optional[str],
    cxx: str = "cxx20",
    flavor: Optional[str] = None,
) -> str:
    # Filter first, so every rule below (family, version, cxx) stays inside one
    # accelerator stack rather than falling back across them.
    arch_lockfiles = {
        comp: v for comp, v in lockfiles[arch].items() if extract_flavor(comp) == flavor
    }

    if not arch_lockfiles:
        # Fail hard: a CUDA job handed the host lockfile fails much later, in a
        # configure step, with an error pointing nowhere near here.
        available = sorted({extract_flavor(c) or "host" for c in lockfiles[arch]})
        print(
            f"No lockfile found for flavor '{flavor or 'host'}' on {arch}; "
            f"available flavors: {', '.join(available)}"
        )
        exit(1)

    # Resolve default lockfile (legacy format has "default", new format may not)
    if "default" in arch_lockfiles:
        _, lockfile = arch_lockfiles["default"]
    else:
        # New format: no default, use first matching cxx variant when compiler
        # unspecified
        cxx_defaults = [
            (comp, url)
            for comp, (_, url) in arch_lockfiles.items()
            if extract_cxx(comp) == cxx
        ]
        if cxx_defaults:
            _, lockfile = cxx_defaults[0]
        else:
            _, lockfile = next(iter(arch_lockfiles.values()))

    if compiler is None:
        return lockfile

    def extract_version(spec: str) -> list:
        """Extract version tuple for comparison (handles gcc@13.3.0 or gcc@13.3.0_cxx20)."""
        base = spec.split("_")[0]
        if "@" not in base:
            return []
        try:
            return [int(v) for v in base.split("@")[1].split(".")]
        except (ValueError, IndexError):
            return []

    # Spack uses "llvm" for clang; treat them as aliases
    compiler_family = compiler.split("@")[0]
    llvm_families = ["clang", "llvm"]
    compiler_families = (
        llvm_families if compiler_family in llvm_families else [compiler_family]
    )

    # Find all specs of matching compiler family (clang and llvm are aliases)
    family_matches = {
        comp: (name, url)
        for comp, (name, url) in arch_lockfiles.items()
        if comp != "default" and comp.split("@")[0] in compiler_families
    }

    if not family_matches:
        print(
            f"No lockfile found for compiler family '{compiler_family}', "
            f"falling back to default: {lockfile}"
        )
        return lockfile

    # Filter to specs matching our exact compiler version
    compiler_version = extract_version(compiler)
    compiler_matches = {
        comp: v
        for comp, v in family_matches.items()
        if compiler_version and extract_version(comp) == compiler_version
    }

    if not compiler_matches:
        # No exact version match: take the highest of the family, but keep the
        # requested C++ standard if it publishes one. Without that second
        # filter `max` breaks the tie in dict order, so an unknown compiler
        # version (clang 22.1.7 -> 22.1.8) silently loses a cxx23 request too.
        highest_version = max(map(extract_version, family_matches.keys()))
        candidates = {
            comp: v
            for comp, v in family_matches.items()
            if extract_version(comp) == highest_version
        }
        cxx_candidates = {
            comp: v for comp, v in candidates.items() if extract_cxx(comp) == cxx
        }
        highest = next(iter(cxx_candidates or candidates))
        _, lockfile = family_matches[highest]
        print(
            f"No lockfile found for compiler '{compiler}', "
            f"falling back to highest version '{highest}': {lockfile}"
        )
        return lockfile

    # Prefer cxx-specific variant; exact family name takes priority over alias
    # (compiler_matches already contains both clang and llvm aliases)
    cxx_matches = {
        comp: v for comp, v in compiler_matches.items() if extract_cxx(comp) == cxx
    }
    if cxx_matches:
        exact_family = next(
            (
                v
                for comp, v in cxx_matches.items()
                if comp.split("@")[0] == compiler_family
            ),
            None,
        )
        _, lockfile = exact_family or next(iter(cxx_matches.values()))
        return lockfile

    # Old format (no cxx suffix)
    if compiler in compiler_matches:
        _, lockfile = compiler_matches[compiler]
        return lockfile

    # Fallback: first available. The normal path for the flavors, published at
    # cxx23 only, so a cxx20 job asking for `cuda13` gets the cxx23 stack --
    # deliberate; see the `flavor` input in .github/actions/dependencies.
    first_comp, (_, lockfile) = next(iter(compiler_matches.items()))
    print(
        f"No lockfile found for compiler '{compiler}' with {cxx} "
        f"(flavor '{flavor or 'host'}'), falling back to '{first_comp}': {lockfile}"
    )
    return lockfile


def determine_compiler_version(binary: str):
    try:
        result = subprocess.run([binary, "--version"], capture_output=True, text=True)

        line = result.stdout.split("\n", 1)[0]
        print(line)
        if "clang" in line:
            compiler = "clang"
            if "Apple" in line:
                compiler = "apple-clang"
        elif "gcc" in line or "GCC" in line or "g++" in line:
            compiler = "gcc"
        else:
            print(f"Unknown compiler: {binary}")
            exit(1)

        m = re.search(r"(\d+\.\d+\.\d+)", line)
        if m is None:
            print(f"Failed to determine version for compiler: {binary}")
            exit(1)
        (version,) = m.groups()
        return f"{compiler}@{version}"

    except (subprocess.SubprocessError, FileNotFoundError):
        print(f"Failed to determine version for compiler: {binary}")
        exit(1)


if __name__ == "__main__":
    main()
