"""Standard-library-only manifests and integrity seals for benchmark hand-offs.

A seal detects changed, missing and additional files. It is a checksum manifest,
not a signature or a claim that a host administrator cannot rewrite a run.
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath


SCORERS = ("potts-metropolis", "typicality-gate", "softmax-repaired", "dphi-product")
SCHEMA_VERSION = "1.0"


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_digest(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _files(run_dir):
    run_dir = Path(run_dir)
    if run_dir.is_symlink() or not run_dir.is_dir():
        raise ValueError("Run must be a real directory, not a symlink")
    files = {}
    for path in sorted(run_dir.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Symlinks are not allowed in a sealed run: {path}")
        if path.is_file() and path.name != "seal.json":
            files[path.relative_to(run_dir).as_posix()] = sha256_file(path)
        elif not path.is_dir() and not path.is_file():
            raise ValueError(f"Unsupported artifact type: {path}")
    # Only the root seal is excluded. A nested seal is an ordinary artifact.
    for path in run_dir.rglob("seal.json"):
        if path.parent != run_dir:
            files[path.relative_to(run_dir).as_posix()] = sha256_file(path)
    return dict(sorted(files.items()))


def seal_run(run_dir):
    run_dir = Path(run_dir)
    if (run_dir / "seal.json").exists():
        raise ValueError("Run is already sealed; create a new output directory")
    files = _files(run_dir)
    if not files:
        raise ValueError("Cannot seal an empty run")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "sealed_at": utc_now(),
        "algorithm": "sha256",
        "files": files,
        "content_digest": canonical_digest(files),
        "scope": "Every regular run artifact except this root seal; no symlinks",
    }
    write_json(run_dir / "seal.json", manifest)
    return manifest


def verify_seal(run_dir):
    run_dir = Path(run_dir)
    if run_dir.is_symlink() or (run_dir / "seal.json").is_symlink():
        raise ValueError("Sealed run and seal must not be symlinks")
    try:
        manifest = json.loads((run_dir / "seal.json").read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read benchmark seal: {exc}") from exc
    if (not isinstance(manifest, dict) or manifest.get("schema_version") != SCHEMA_VERSION
            or manifest.get("algorithm") != "sha256"):
        raise ValueError("Unsupported benchmark seal schema or hash algorithm")
    expected = manifest.get("files")
    if not isinstance(expected, dict) or not expected:
        raise ValueError("Invalid benchmark seal file manifest")
    for name, digest in expected.items():
        path = PurePosixPath(name)
        if (path.is_absolute() or ".." in path.parts or str(path) != name
                or name == "seal.json" or not isinstance(digest, str)
                or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest)):
            raise ValueError(f"Invalid sealed artifact path or digest: {name}")
    if canonical_digest(expected) != manifest.get("content_digest"):
        raise ValueError("Seal manifest checksum mismatch")
    actual = _files(run_dir)
    if actual != expected:
        missing = sorted(set(expected) - set(actual))
        added = sorted(set(actual) - set(expected))
        changed = sorted(k for k in set(actual) & set(expected) if actual[k] != expected[k])
        raise ValueError(f"Run seal verification failed: missing={missing}, added={added}, changed={changed}")
    return manifest
