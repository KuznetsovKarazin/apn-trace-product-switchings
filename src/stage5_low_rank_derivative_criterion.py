#!/usr/bin/env python3
"""Stage 5E pilot: exact low-rank derivative criterion for APN switching.

Let F: F_2^n -> F_2^n be quadratic APN and, for a != 0,

    L_a^F(x) = F(x+a) + F(x) + F(a).

Then ker L_a^F = <a>.  For a perturbation

    Delta(x) = sum_i v_i q_i(x),

with Boolean quadratic q_i, write

    ell_{a,i}(x) = q_i(x+a)+q_i(x)+q_i(a).

If U(y)=sum_i y_i v_i and R_a(x)=(ell_{a,i}(x))_i, then
F+Delta fails APN at a exactly when there exists y != 0 such that

    U(y) is in im L_a^F,
    R_a((L_a^F)^(-1) U(y)) = y.

The inverse is well-defined modulo <a>, and every ell_{a,i}(a)=0, so the
criterion does not depend on the chosen preimage.

This script validates the criterion on the five canonical scalar families
already scanned by brute-force APN testing.  It performs no Gröbner search.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

GF_MOD = 0x11B
N = 256


def gf_mult(a: int, b: int) -> int:
    a &= 0xFF
    b &= 0xFF
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 0x100:
            a ^= GF_MOD
    return out & 0xFF


def gf_pow(a: int, exponent: int) -> int:
    result = 1
    base = a & 0xFF
    e = int(exponent)
    while e:
        if e & 1:
            result = gf_mult(result, base)
        base = gf_mult(base, base)
        e >>= 1
    return result


def gf_order(a: int) -> int:
    if a == 0:
        return 0
    value = 1
    for order in range(1, 256):
        value = gf_mult(value, a)
        if value == 1:
            return order
    raise RuntimeError("invalid nonzero GF(256) element")


def gold_sbox(exponent: int) -> list[int]:
    return [gf_pow(x, exponent) for x in range(N)]


def q_value(edges: list[tuple[int, int]], x: int) -> int:
    value = 0
    for i, j in edges:
        value ^= ((x >> i) & 1) & ((x >> j) & 1)
    return value


def derivative_lut(function: list[int], a: int) -> list[int]:
    return [function[x ^ a] ^ function[x] ^ function[a] for x in range(N)]


def preimage_map_of_derivative(function: list[int], a: int) -> dict[int, int]:
    derivative = derivative_lut(function, a)
    preimage: dict[int, int] = {}
    for x, value in enumerate(derivative):
        preimage.setdefault(value, x)
    if len(preimage) != 128:
        raise RuntimeError(
            f"base function is not APN at a=0x{a:02x}: image size {len(preimage)}"
        )
    # Every image value must have exactly two preimages separated by a.
    counts = Counter(derivative)
    if set(counts.values()) != {2}:
        raise RuntimeError(f"unexpected derivative fibre sizes at a=0x{a:02x}")
    for value, x in preimage.items():
        if derivative[x ^ a] != value:
            raise RuntimeError("derivative kernel is not <a>")
    return preimage


def output_rank(vectors: list[int]) -> int:
    rows = [value for value in vectors if value]
    rank = 0
    for bit in range(8):
        pivot = next(
            (i for i in range(rank, len(rows)) if (rows[i] >> bit) & 1),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and ((rows[i] >> bit) & 1):
                rows[i] ^= rows[rank]
        rank += 1
        if rank == len(rows):
            break
    return rank


def precompute_base(exponent: int) -> dict[str, Any]:
    function = gold_sbox(exponent)
    preimages = {
        a: preimage_map_of_derivative(function, a)
        for a in range(1, N)
    }
    return {
        "exponent": exponent,
        "function": function,
        "preimages": preimages,
    }


def precompute_boolean_derivatives(
    kernels: dict[str, list[tuple[int, int]]],
) -> dict[str, dict[int, list[int]]]:
    result: dict[str, dict[int, list[int]]] = {}
    for name, edges in kernels.items():
        truth = [q_value(edges, x) for x in range(N)]
        result[name] = {
            a: [truth[x ^ a] ^ truth[x] ^ truth[a] for x in range(N)]
            for a in range(1, N)
        }
    return result


def criterion_witness(
    base: dict[str, Any],
    boolean_derivatives: dict[str, dict[int, list[int]]],
    components: list[tuple[int, str]],
) -> dict[str, int] | None:
    """Return first APN-destroying witness, or None if APN is preserved."""
    rank = len(components)
    if rank == 0:
        return None
    vectors = [value for value, _ in components]
    if output_rank(vectors) != rank:
        raise ValueError("components must use an F_2-independent output basis")

    for a in range(1, N):
        preimage = base["preimages"][a]
        for y in range(1, 1 << rank):
            output = 0
            for i, vector in enumerate(vectors):
                if (y >> i) & 1:
                    output ^= vector
            x = preimage.get(output)
            if x is None:
                continue
            recovered = 0
            for i, (_, kernel_name) in enumerate(components):
                if boolean_derivatives[kernel_name][a][x]:
                    recovered |= 1 << i
            if recovered == y:
                return {
                    "a": a,
                    "a_hex": f"0x{a:02x}",
                    "y": y,
                    "preimage_x": x,
                    "preimage_x_hex": f"0x{x:02x}",
                    "output_combination": output,
                    "output_combination_hex": f"0x{output:02x}",
                }
    return None


def parameter_components(family: dict[str, Any], parameter: int):
    mode = family["parameter_mode"]
    if mode == "direct_single_coefficient":
        return [(parameter, family["components"][0][1])]
    if mode == "global_field_multiplier":
        return [
            (gf_mult(parameter, vector), kernel)
            for vector, kernel in family["components"]
        ]
    raise ValueError(mode)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kernels",
        type=Path,
        default=Path("results/evidence/canonical_switching_kernels.json"),
    )
    parser.add_argument(
        "--archived-scan",
        type=Path,
        default=Path("results/evidence/scalar_family_completion/scalar_family_apn_scan.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/low_rank_derivative_criterion.json"),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]

    def resolve(path: Path) -> Path:
        return path if path.is_absolute() else root / path

    kernels_path = resolve(args.kernels)
    archived_path = resolve(args.archived_scan)
    output_path = resolve(args.output)

    kernel_payload = json.loads(kernels_path.read_text(encoding="utf-8"))
    archived = json.loads(archived_path.read_text(encoding="utf-8"))
    kernels = {
        name: [tuple(edge) for edge in record["expanded_terms"]]
        for name, record in kernel_payload["kernels"].items()
    }
    boolean_derivatives = precompute_boolean_derivatives(kernels)
    bases = {3: precompute_base(3), 9: precompute_base(9)}

    families = [
        {
            "id": "X3-Q0",
            "gold_exponent": 3,
            "parameter_mode": "direct_single_coefficient",
            "components": [(1, "q0")],
        },
        {
            "id": "X3-R2",
            "gold_exponent": 3,
            "parameter_mode": "global_field_multiplier",
            "components": [(0x3D, "q1"), (0x24, "q0")],
        },
        {
            "id": "X9-Q0",
            "gold_exponent": 9,
            "parameter_mode": "global_field_multiplier",
            "components": [(0xF2, "q0")],
        },
        {
            "id": "X9-Q2",
            "gold_exponent": 9,
            "parameter_mode": "global_field_multiplier",
            "components": [(0x0C, "q2")],
        },
        {
            "id": "X9-R2",
            "gold_exponent": 9,
            "parameter_mode": "global_field_multiplier",
            "components": [(0x6F, "q1"), (0xF2, "q0")],
        },
    ]

    archived_by_id = {
        record["id"]: record for record in archived["families"]
    }
    family_records = []
    all_match = True
    for family in families:
        accepted = []
        rejected_witness_histogram = Counter()
        first_rejected_examples = []
        component_rank_distribution = Counter()
        for parameter in range(1, 256):
            components = parameter_components(family, parameter)
            component_rank = output_rank([value for value, _ in components])
            component_rank_distribution[component_rank] += 1
            witness = criterion_witness(
                bases[family["gold_exponent"]],
                boolean_derivatives,
                components,
            )
            if witness is None:
                accepted.append(parameter)
            else:
                rejected_witness_histogram[witness["a_hex"]] += 1
                if len(first_rejected_examples) < 5:
                    first_rejected_examples.append({
                        "parameter": parameter,
                        "parameter_hex": f"0x{parameter:02x}",
                        "witness": witness,
                    })

        expected = [
            int(value, 16)
            for value in archived_by_id[family["id"]]["apn_parameters"]
        ]
        match = accepted == expected
        all_match &= match
        family_records.append({
            **family,
            "base_component_output_rank": output_rank(
                [value for value, _ in family["components"]]
            ),
            "parameter_count_tested": 255,
            "accepted_parameters": accepted,
            "accepted_parameters_hex": [f"0x{x:02x}" for x in accepted],
            "accepted_parameter_orders": {
                f"0x{x:02x}": gf_order(x) for x in accepted
            },
            "archived_bruteforce_parameters_hex": [
                f"0x{x:02x}" for x in expected
            ],
            "exact_match_with_archived_bruteforce_scan": match,
            "component_rank_distribution": {
                str(rank): count
                for rank, count in sorted(component_rank_distribution.items())
            },
            "rejected_first_witness_a_histogram": dict(
                rejected_witness_histogram.most_common()
            ),
            "sample_rejected_witnesses": first_rejected_examples,
        })

    x3_q0 = next(record for record in family_records if record["id"] == "X3-Q0")
    q0_values = x3_q0["accepted_parameters"]
    lambda_value = gf_mult(q0_values[1], gf_pow(q0_values[0], 254))
    q0_subspace_checks = {
        "accepted_plus_zero_closed_under_xor": all(
            (a ^ b) in ({0} | set(q0_values))
            for a in ({0} | set(q0_values))
            for b in ({0} | set(q0_values))
        ),
        "dimension_over_F2": 2,
        "base_alpha_hex": f"0x{q0_values[0]:02x}",
        "lambda_hex": f"0x{lambda_value:02x}",
        "lambda_order": gf_order(lambda_value),
        "alpha_times_lambda_hex": f"0x{gf_mult(q0_values[0], lambda_value):02x}",
        "alpha_times_lambda_squared_hex": (
            f"0x{gf_mult(q0_values[0], gf_mult(lambda_value, lambda_value)):02x}"
        ),
        "interpretation": (
            "The three nonzero X3-Q0 coefficients are exactly alpha*GF(4)^*, "
            "with alpha=0x24 and lambda=0xbd of order three."
        ),
    }

    report = {
        "schema": "apn-low-rank-derivative-criterion-v1",
        "date": "2026-08-02",
        "scope": (
            "Exact rank-one/rank-two derivative update criterion, validated "
            "against five archived scalar-family APN scans"
        ),
        "provenance": {
            "canonical_kernels_file": str(kernels_path.relative_to(root)),
            "canonical_kernels_sha256": sha256_file(kernels_path),
            "archived_scalar_scan_file": str(archived_path.relative_to(root)),
            "archived_scalar_scan_sha256": sha256_file(archived_path),
        },
        "criterion": {
            "base_derivative": "L_a^F(x)=F(x+a)+F(x)+F(a)",
            "APN_base_condition": "ker L_a^F=<a> for every a!=0",
            "perturbation": "Delta=sum_i v_i q_i",
            "boolean_derivatives": "ell_{a,i}(x)=q_i(x+a)+q_i(x)+q_i(a)",
            "failure_condition": (
                "There exist a!=0 and y!=0 such that U(y) is in im L_a^F "
                "and R_a((L_a^F)^(-1)U(y))=y."
            ),
            "well_defined_reason": (
                "Two preimages differ by a and ell_{a,i}(a)=0, so the recovered "
                "Boolean derivative vector is independent of the preimage."
            ),
            "rank_one_specialization": (
                "For Delta=v q, failure at a occurs exactly when v is in im "
                "L_a^F and ell_{a,q}((L_a^F)^(-1)v)=1."
            ),
            "complexity": (
                "After precomputing the 255 derivative image/preimage tables, "
                "a rank-r candidate needs at most 255*(2^r-1) small lookups."
            ),
        },
        "families": family_records,
        "all_five_families_match_archived_bruteforce_scan": all_match,
        "x3_q0_GF4_coset_structure": q0_subspace_checks,
        "main_conclusions": [
            "The low-rank derivative criterion reproduces exactly every accepted scalar parameter in all five archived families.",
            "No Groebner basis and no exhaustive differential table are needed once derivative preimages are cached.",
            "The X3-Q0 accepted set {0x24,0x88,0xac} is the multiplicative coset 0x24*GF(4)^*; together with zero it is a two-dimensional F_2 subspace.",
            "The criterion is suitable for symbolic analysis of q0,q1,q2 and for symmetry-guided synthesis of new low-rank APN switchings.",
        ],
    }

    if not all_match:
        raise RuntimeError("low-rank criterion does not match archived scan")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "all_five_families_match": all_match,
        "accepted": {
            record["id"]: record["accepted_parameters_hex"]
            for record in family_records
        },
        "x3_q0_structure": q0_subspace_checks,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
