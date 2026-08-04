#!/usr/bin/env python3
"""Publication correction for the affine off-axis character-sum bound.

The complete projective sum has the 4*sqrt(q) bound.  In the off-axis case one
point at infinity is regular (after flex-tangent cancellation) and contributes
1.  The affine sum therefore satisfies |S_aff| <= 4*sqrt(q)+1.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    old_path = root / "results/current/manuscript_uniform_theorem.json"
    scan_path = root / "results/current/trace_switch_dimension_scan.json"
    out_path = root / "results/current/publication_uniform_theorem_corrected.json"

    old = json.loads(old_path.read_text(encoding="utf-8"))
    scan = json.loads(scan_path.read_text(encoding="utf-8"))
    rows = {row["n"]: row for row in scan["dimensions"]}

    threshold = []
    for n in range(10, 42, 2):
        q = 1 << n
        sqrt_q = 1 << (n // 2)
        numerator = q - 11 - 74 * sqrt_q
        threshold.append(
            {
                "n": n,
                "q": q,
                "sqrt_q": sqrt_q,
                "numerator": numerator,
                "lower_bound": numerator / 16,
                "positive": numerator > 0,
            }
        )

    data = {
        "schema": "publication-uniform-theorem-corrected-v1",
        "date": "2026-08-02",
        "theorem": old["theorem"],
        "correction": {
            "issue": (
                "The earlier draft applied the complete-curve 4*sqrt(q) bound "
                "directly to the affine off-axis sum."
            ),
            "geometry": (
                "At the cancelled infinity point the phase is regular with value 0, "
                "so the complete sum contains a contribution +1 that is absent from "
                "the affine sum."
            ),
            "off_axis_affine_bound": "|S_aff| <= 4*sqrt(q)+1",
            "old_lower_bound": old["lower_bound"],
            "corrected_lower_bound": "h_n(theta) >= (q-11-74 sqrt(q))/16",
            "theorem_unchanged": True,
        },
        "error_accounting": {
            "hasse": 2,
            "axis_pairs": 6,
            "axis_each": 6,
            "off_axis_pairs": 9,
            "off_axis_each_complete": 4,
            "off_axis_affine_additive_constant_each": 1,
            "sqrt_constant": 74,
            "constant_term": 11,
        },
        "threshold": threshold,
        "finite_bridge": [
            {"n": 10, "minimum_h": rows[10]["minimum_positive_representation_count"]},
            {"n": 12, "minimum_h": rows[12]["minimum_positive_representation_count"]},
        ],
        "validation": {
            "positive_from_n14": all(r["positive"] for r in threshold if r["n"] >= 14),
            "n10_bridge": rows[10]["minimum_positive_representation_count"] == 46,
            "n12_bridge": rows[12]["minimum_positive_representation_count"] == 208,
            "old_json_sha256": sha256(old_path),
            "dimension_scan_sha256": sha256(scan_path),
        },
    }
    out_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    print(sha256(out_path))


if __name__ == "__main__":
    main()
