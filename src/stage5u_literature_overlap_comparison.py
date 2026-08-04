#!/usr/bin/env python3
"""Stage 5U: exact local comparison with a published hyperplane modification.

Reconstruct Taniguchi--Polujan--Pott--Arshad, Example 4.2, in the
project's GF(256) model and compare its portable orthoderivative signature
with the exact historical correspondence already stored in the checkpoint.

This script does not establish novelty by itself. It records an auditable
identity of computed invariants and links it to an existing exact CCZ witness
from results/evidence/historical_correspondence.csv.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code/current"))

from stage5_canonical_low_rank_locus import (  # noqa: E402
    absolute_trace,
    gf_mult,
    gf_pow,
    portable_ortho_signature,
    sbox_sha256,
)

LAMBDA = 0xBD  # primitive element of the unique GF(4) subfield in this model


def relative_trace_to_f4(x: int) -> int:
    return x ^ gf_pow(x, 4) ^ gf_pow(x, 16) ^ gf_pow(x, 64)


def is_apn(function: list[int]) -> bool:
    for a in range(1, 256):
        counts = [0] * 256
        for x in range(256):
            counts[function[x ^ a] ^ function[x]] += 1
        if max(counts) > 2:
            return False
    return True


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    f4 = sorted(x for x in range(256) if gf_pow(x, 4) == x)
    assert f4 == [0x00, 0x01, 0xBC, 0xBD]
    assert gf_pow(LAMBDA, 3) == 1 and LAMBDA not in (0, 1)

    # Published Example 4.2, after choosing the project's representation of GF(4):
    # G(x)=x^3 + lambda * Tr_1^8(x) * Tr_2^8(x).
    function = []
    for x in range(256):
        modification = (
            gf_mult(LAMBDA, relative_trace_to_f4(x))
            if absolute_trace(x)
            else 0
        )
        function.append(gf_pow(x, 3) ^ modification)

    assert is_apn(function)
    signature = portable_ortho_signature(function)

    rows = list(csv.DictReader(
        (ROOT / "results/evidence/historical_correspondence.csv").open(
            newline="", encoding="utf-8"
        )
    ))
    matches = [r for r in rows if r["ortho_signature_sha256"] == signature["sha256"]]
    assert len(matches) == 1
    historical = matches[0]
    assert historical["internal_id"] == "CLASS-B"
    assert historical["exact_ccz"] == "True"

    canonical = json.loads(
        (ROOT / "results/current/canonical_low_rank_locus.json").read_text(
            encoding="utf-8"
        )
    )
    class_b_points = []

    def walk(obj):
        if isinstance(obj, dict):
            if (
                obj.get("ortho_signature_sha256") == signature["sha256"]
                and "CLASS-B" in obj.get("signature_compatible_classes", [])
            ):
                class_b_points.append(obj.get("label", obj.get("orbit_index")))
            for value in obj.values():
                walk(value)
        elif isinstance(obj, list):
            for value in obj:
                walk(value)

    walk(canonical)
    assert class_b_points

    result = {
        "schema": "literature-overlap-comparison-v1",
        "date": "2026-08-02",
        "published_example": {
            "source": (
                "Taniguchi--Polujan--Pott--Arshad, Changing Almost Perfect "
                "Nonlinear Functions on Affine Subspaces of Small Codimensions, "
                "Example 4.2"
            ),
            "normalized_formula": (
                "G(x)=x^3+lambda*Tr_{GF(256)/GF(2)}(x)*"
                "Tr_{GF(256)/GF(4)}(x)"
            ),
            "lambda_hex": "0xbd",
            "field_model": "GF(2)[z]/(z^8+z^4+z^3+z+1)",
            "apn_verified": True,
            "function_sha256": sbox_sha256(function),
            "orthoderivative_signature": signature,
        },
        "local_match": {
            "internal_id": historical["internal_id"],
            "historical_family": historical["literature_family"],
            "historical_origin": historical["literature_origin"],
            "representative_formula_latex": historical["representative_formula_latex"],
            "exact_ccz_previously_verified": historical["exact_ccz"] == "True",
            "known_representative_number": int(historical["ep09_representative_no"]),
            "matched_canonical_low_rank_entries": class_b_points,
        },
        "interpretation": (
            "The published dimension-eight hyperplane modification has the "
            "same portable orthoderivative signature as CLASS-B. The checkpoint "
            "already contains an independent exact CCZ witness identifying "
            "CLASS-B with Edel--Pott catalogue representative No. 6, "
            "x^9+Tr(x^3). Therefore this example is not evidence for a new "
            "CCZ class in the present project; it is a literature-overlap anchor."
        ),
    }

    output = ROOT / "results/current/literature_overlap_comparison.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "sha256": sha256(output),
        "signature_sha256": signature["sha256"],
        "matched_class": historical["internal_id"],
    }, indent=2))


if __name__ == "__main__":
    main()
