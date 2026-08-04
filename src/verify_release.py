#!/usr/bin/env python3
"""Verify the canonical core certificates shipped with this repository."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("REPRODUCIBILITY_MANIFEST.json"),
        help="Repository-relative manifest path.",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    manifest_path = args.manifest if args.manifest.is_absolute() else root / args.manifest
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    failures: list[str] = []
    checked: list[dict[str, str]] = []
    for item in manifest["canonical_outputs"]:
        path = root / item["path"]
        if not path.is_file():
            failures.append(f"missing: {item['path']}")
            continue
        actual = digest(path)
        checked.append({"path": item["path"], "sha256": actual})
        if actual != item["sha256"]:
            failures.append(
                f"hash mismatch: {item['path']} expected={item['sha256']} actual={actual}"
            )

    print(json.dumps({"checked": checked, "failures": failures}, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
