#!/usr/bin/env python3
"""Stage 5E/F: exact rank<=2 APN locus in the canonical q0,q1,q2 kernel space.

The script works entirely over GF(2^8) with modulus 0x11B and requires no
Sage, Magma, sboxU, or Gröbner basis computation.

For a quadratic APN base F and Boolean kernels Q=(q_1,...,q_r), define

    L_a(x) = F(x+a)+F(x)+F(a),
    R_a(x) = (D_a q_1(x),...,D_a q_r(x)).

For each nonzero y in F_2^r, precompute the forbidden output set

    B_y = {u : exists a != 0, u in im L_a and
                 R_a(L_a^{-1}(u)) = y}.

Then Delta=sum_i v_i q_i preserves APN exactly when

    U(y)=sum_i y_i v_i not in B_y  for every y != 0.

This gives a closed finite-set criterion.  The script enumerates every
nonzero coefficient-rank-one and coefficient-rank-two perturbation in

    GF(256) tensor_F2 span(q0,q1,q2),

computes portable orthoderivative signatures, partitions the accepted
perturbations under the standard AΓL(1,256) stabilizer of Gold x^3/x^9,
and records the cubic/trace bridge between the q0,q1 families.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, deque
from pathlib import Path
from typing import Any

GF_MOD = 0x11B
N = 256
GROUP_ORDER = 256 * 255 * 8
KERNEL_NAMES = ("q0", "q1", "q2")


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


def gf_inv(a: int) -> int:
    if a == 0:
        raise ZeroDivisionError("zero has no inverse")
    return gf_pow(a, 254)


def gf_order(a: int) -> int:
    if a == 0:
        return 0
    value = 1
    for order in range(1, 256):
        value = gf_mult(value, a)
        if value == 1:
            return order
    raise RuntimeError("invalid field element")


def absolute_trace(a: int) -> int:
    value = 0
    term = a & 0xFF
    for _ in range(8):
        value ^= term
        term = gf_mult(term, term)
    if value not in (0, 1):
        raise AssertionError("absolute trace did not land in GF(2)")
    return value


def parity(a: int) -> int:
    return a.bit_count() & 1


def q_value(edges: list[tuple[int, int]], x: int) -> int:
    value = 0
    for i, j in edges:
        value ^= ((x >> i) & 1) & ((x >> j) & 1)
    return value


def gold_sbox(exponent: int) -> list[int]:
    return [gf_pow(x, exponent) for x in range(N)]


def mask_label(mask: int) -> str:
    terms = [KERNEL_NAMES[i] for i in range(3) if (mask >> i) & 1]
    return "+".join(terms)


def mask_truth(base_truths: list[list[int]], mask: int) -> list[int]:
    return [
        ((base_truths[0][x] if mask & 1 else 0)
         ^ (base_truths[1][x] if mask & 2 else 0)
         ^ (base_truths[2][x] if mask & 4 else 0))
        for x in range(N)
    ]


def derivative_preimages(function: list[int]) -> dict[int, dict[int, int]]:
    result: dict[int, dict[int, int]] = {}
    for a in range(1, N):
        preimage: dict[int, int] = {}
        derivative = [
            function[x ^ a] ^ function[x] ^ function[a]
            for x in range(N)
        ]
        for x, value in enumerate(derivative):
            preimage.setdefault(value, x)
        if len(preimage) != 128:
            raise RuntimeError(
                f"base is not APN at a=0x{a:02x}: image={len(preimage)}"
            )
        result[a] = preimage
    return result


def boolean_derivatives(truth: list[int]) -> dict[int, list[int]]:
    return {
        a: [truth[x ^ a] ^ truth[x] ^ truth[a] for x in range(N)]
        for a in range(1, N)
    }


def rank_one_forbidden_set(
    preimages: dict[int, dict[int, int]],
    derivative: dict[int, list[int]],
) -> set[int]:
    forbidden: set[int] = set()
    for a, preimage in preimages.items():
        d = derivative[a]
        for output, x in preimage.items():
            if d[x]:
                forbidden.add(output)
    return forbidden


def rank_two_exact_sets(
    preimages: dict[int, dict[int, int]],
    derivative_1: dict[int, list[int]],
    derivative_2: dict[int, list[int]],
) -> list[set[int]]:
    exact = [set() for _ in range(4)]
    for a, preimage in preimages.items():
        d1 = derivative_1[a]
        d2 = derivative_2[a]
        for output, x in preimage.items():
            recovered = d1[x] | (d2[x] << 1)
            exact[recovered].add(output)
    return exact


def dot3(a: int, b: int) -> int:
    return parity(a & b)


def output_rank_two(a: int, b: int) -> bool:
    return a != 0 and b != 0 and a != b


def perturbation_rank_one(
    kernel_truths: dict[int, list[int]], mask: int, coefficient: int
) -> bytes:
    truth = kernel_truths[mask]
    return bytes(coefficient if truth[x] else 0 for x in range(N))


def perturbation_rank_two(
    kernel_truths: dict[int, list[int]],
    mask_1: int,
    coefficient_1: int,
    mask_2: int,
    coefficient_2: int,
) -> bytes:
    t1 = kernel_truths[mask_1]
    t2 = kernel_truths[mask_2]
    return bytes(
        (coefficient_1 if t1[x] else 0)
        ^ (coefficient_2 if t2[x] else 0)
        for x in range(N)
    )


def add_perturbation(base: list[int], perturbation: bytes) -> list[int]:
    return [base[x] ^ perturbation[x] for x in range(N)]


def sbox_sha256(values: list[int] | bytes) -> str:
    return hashlib.sha256(bytes(values)).hexdigest()


def orthoderivative(function: list[int]) -> list[int]:
    result = [0] * N
    for a in range(1, N):
        image = {
            function[x ^ a] ^ function[x] ^ function[a]
            for x in range(N)
        }
        if len(image) != 128:
            raise RuntimeError("orthoderivative requested for non-APN function")
        candidates = [
            w for w in range(1, N)
            if all(parity(w & value) == 0 for value in image)
        ]
        if len(candidates) != 1:
            raise RuntimeError(
                f"left kernel at a=0x{a:02x} has {len(candidates)} generators"
            )
        result[a] = candidates[0]
    return result


def sage_counter_repr(counter: Counter[int]) -> str:
    return "{" + ", ".join(
        f"{key}:{counter[key]}" for key in sorted(counter)
    ) + "}"


def portable_ortho_signature(function: list[int]) -> dict[str, Any]:
    ortho = orthoderivative(function)

    differential = Counter()
    for a in range(1, N):
        row = [0] * N
        for x in range(N):
            row[ortho[x ^ a] ^ ortho[x]] += 1
        differential.update(row)

    walsh = Counter()
    for output_mask in range(1, N):
        values = [
            1 if parity(output_mask & ortho[x]) == 0 else -1
            for x in range(N)
        ]
        step = 1
        while step < N:
            for offset in range(0, N, 2 * step):
                for index in range(offset, offset + step):
                    left = values[index]
                    right = values[index + step]
                    values[index] = left + right
                    values[index + step] = left - right
            step *= 2
        walsh.update(abs(value) for value in values)

    signature = {
        "differential_spectrum": sage_counter_repr(differential),
        "absolute_walsh_spectrum": sage_counter_repr(walsh),
    }
    serialization = json.dumps(
        signature, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    )
    return {
        **signature,
        "serialization": serialization,
        "sha256": hashlib.sha256(serialization.encode("utf-8")).hexdigest(),
    }


def primitive_element() -> int:
    for value in range(2, N):
        if gf_order(value) == 255:
            return value
    raise RuntimeError("primitive element not found")


def build_gold_action(a: int, frobenius: int, translation: int, exponent: int):
    inverse_a = gf_inv(a)
    inverse_frobenius = (8 - frobenius) % 8
    input_permutation = [
        gf_pow(gf_mult(x ^ translation, inverse_a), 1 << inverse_frobenius)
        for x in range(N)
    ]
    output_factor = gf_pow(a, exponent)
    output_lut = bytes(
        gf_mult(output_factor, gf_pow(value, 1 << frobenius))
        for value in range(N)
    )
    return input_permutation, output_lut


def gold_stabilizer_generators(exponent: int):
    primitive = primitive_element()
    inverse_primitive = gf_inv(primitive)
    return [
        build_gold_action(1, 0, 1, exponent),
        build_gold_action(primitive, 0, 0, exponent),
        build_gold_action(inverse_primitive, 0, 0, exponent),
        build_gold_action(1, 1, 0, exponent),
        build_gold_action(1, 7, 0, exponent),
    ]


def apply_gold_action(state: bytes, action) -> bytes:
    input_permutation, output_lut = action
    return bytes(
        output_lut[state[input_permutation[x]]]
        for x in range(N)
    )


def orbit_summary(
    start: bytes,
    generators,
    labelled_states: dict[bytes, list[str]],
) -> dict[str, Any]:
    seen = {start}
    queue = deque([start])
    canonical = start
    hits = set(labelled_states.get(start, []))
    while queue:
        state = queue.popleft()
        for generator in generators:
            image = apply_gold_action(state, generator)
            if image in seen:
                continue
            seen.add(image)
            queue.append(image)
            if image < canonical:
                canonical = image
            hits.update(labelled_states.get(image, []))
    return {
        "orbit_size": len(seen),
        "stabilizer_size": GROUP_ORDER // len(seen),
        "canonical_sha256": sbox_sha256(canonical),
        "accepted_labels_hit": sorted(hits),
    }


def load_signature_classes(path: Path) -> dict[str, list[str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result: dict[str, list[str]] = {}
    for record in payload["results"]:
        result[record["ortho_signature_sha256"]] = record["signature_matches"]
    return result


def enumerate_locus(
    exponent: int,
    kernel_truths: dict[int, list[int]],
    derivatives: dict[int, dict[int, list[int]]],
    signature_classes: dict[str, list[str]],
) -> dict[str, Any]:
    base = gold_sbox(exponent)
    preimages = derivative_preimages(base)

    accepted_states: dict[bytes, list[str]] = {}
    accepted_records: list[dict[str, Any]] = []
    rank_one_records = []

    for mask in range(1, 8):
        forbidden = rank_one_forbidden_set(preimages, derivatives[mask])
        accepted = [value for value in range(1, N) if value not in forbidden]
        mask_record = {
            "kernel_mask": mask,
            "kernel_mask_binary": f"{mask:03b}",
            "kernel": mask_label(mask),
            "forbidden_size": len(forbidden),
            "accepted_coefficients": accepted,
            "accepted_coefficients_hex": [f"0x{x:02x}" for x in accepted],
        }
        rank_one_records.append(mask_record)
        for coefficient in accepted:
            label = f"R1:{mask_label(mask)}:0x{coefficient:02x}"
            perturbation = perturbation_rank_one(
                kernel_truths, mask, coefficient
            )
            function = add_perturbation(base, perturbation)
            signature = portable_ortho_signature(function)
            record = {
                "label": label,
                "coefficient_rank": 1,
                "kernel_mask": mask,
                "kernel": mask_label(mask),
                "coefficient": coefficient,
                "coefficient_hex": f"0x{coefficient:02x}",
                "coefficient_cube_hex": f"0x{gf_pow(coefficient, 3):02x}",
                "coefficient_trace": absolute_trace(coefficient),
                "perturbation_sha256": sbox_sha256(perturbation),
                "function_sha256": sbox_sha256(function),
                "ortho_signature_sha256": signature["sha256"],
                "signature_compatible_classes": signature_classes.get(
                    signature["sha256"], []
                ),
            }
            accepted_records.append(record)
            accepted_states.setdefault(perturbation, []).append(label)

    rank_two_records = []
    for kernel_vector in range(1, 8):
        perpendicular = [
            mask for mask in range(1, 8)
            if dot3(mask, kernel_vector) == 0
        ]
        if len(perpendicular) != 3:
            raise AssertionError("unexpected projective plane")
        mask_1, mask_2 = perpendicular[0], perpendicular[1]
        if mask_1 ^ mask_2 != perpendicular[2]:
            raise AssertionError("chosen masks do not span the perpendicular")

        exact = rank_two_exact_sets(
            preimages, derivatives[mask_1], derivatives[mask_2]
        )
        accepted_pairs = []
        for coefficient_1 in range(1, N):
            if coefficient_1 in exact[1]:
                continue
            for coefficient_2 in range(1, N):
                if not output_rank_two(coefficient_1, coefficient_2):
                    continue
                if coefficient_2 in exact[2]:
                    continue
                if (coefficient_1 ^ coefficient_2) in exact[3]:
                    continue
                accepted_pairs.append((coefficient_1, coefficient_2))

        rank_two_records.append({
            "coefficient_map_kernel": kernel_vector,
            "coefficient_map_kernel_binary": f"{kernel_vector:03b}",
            "basis_masks": [mask_1, mask_2],
            "basis_kernels": [mask_label(mask_1), mask_label(mask_2)],
            "forbidden_exact_sizes": {
                "01": len(exact[1]),
                "10": len(exact[2]),
                "11": len(exact[3]),
            },
            "accepted_pair_count": len(accepted_pairs),
            "accepted_pairs": [
                {
                    "coefficient_1": a,
                    "coefficient_1_hex": f"0x{a:02x}",
                    "coefficient_2": b,
                    "coefficient_2_hex": f"0x{b:02x}",
                    "projective_slope_hex": f"0x{gf_mult(b, gf_inv(a)):02x}",
                    "projective_slope_order": gf_order(
                        gf_mult(b, gf_inv(a))
                    ),
                }
                for a, b in accepted_pairs
            ],
        })

        for coefficient_1, coefficient_2 in accepted_pairs:
            label = (
                f"R2:{mask_label(mask_1)}=0x{coefficient_1:02x}:"
                f"{mask_label(mask_2)}=0x{coefficient_2:02x}:"
                f"ker={kernel_vector:03b}"
            )
            perturbation = perturbation_rank_two(
                kernel_truths,
                mask_1,
                coefficient_1,
                mask_2,
                coefficient_2,
            )
            function = add_perturbation(base, perturbation)
            signature = portable_ortho_signature(function)
            record = {
                "label": label,
                "coefficient_rank": 2,
                "coefficient_map_kernel": kernel_vector,
                "coefficient_map_kernel_binary": f"{kernel_vector:03b}",
                "basis_masks": [mask_1, mask_2],
                "basis_kernels": [mask_label(mask_1), mask_label(mask_2)],
                "coefficients": [coefficient_1, coefficient_2],
                "coefficients_hex": [
                    f"0x{coefficient_1:02x}",
                    f"0x{coefficient_2:02x}",
                ],
                "projective_slope_hex": (
                    f"0x{gf_mult(coefficient_2, gf_inv(coefficient_1)):02x}"
                ),
                "perturbation_sha256": sbox_sha256(perturbation),
                "function_sha256": sbox_sha256(function),
                "ortho_signature_sha256": signature["sha256"],
                "signature_compatible_classes": signature_classes.get(
                    signature["sha256"], []
                ),
            }
            accepted_records.append(record)
            accepted_states.setdefault(perturbation, []).append(label)

    generators = gold_stabilizer_generators(exponent)
    unassigned = set(accepted_states)
    orbit_records = []
    while unassigned:
        start = min(unassigned)
        orbit = orbit_summary(start, generators, accepted_states)
        hit_states = {
            state for state, labels in accepted_states.items()
            if any(label in orbit["accepted_labels_hit"] for label in labels)
        }
        labels = orbit["accepted_labels_hit"]
        compatible = sorted({
            target
            for record in accepted_records
            if record["label"] in labels
            for target in record["signature_compatible_classes"]
        })
        orbit_records.append({
            "orbit_index": len(orbit_records) + 1,
            **orbit,
            "signature_compatible_classes": compatible,
        })
        unassigned.difference_update(hit_states)

    orbit_records.sort(key=lambda rec: rec["canonical_sha256"])
    for index, record in enumerate(orbit_records, 1):
        record["orbit_index"] = index

    return {
        "gold_exponent": exponent,
        "ambient_kernel_space": "GF(256) tensor_F2 span(q0,q1,q2)",
        "ambient_binary_dimension": 24,
        "zero_perturbation_excluded": True,
        "rank_one": {
            "accepted_count": sum(
                len(record["accepted_coefficients"])
                for record in rank_one_records
            ),
            "by_kernel": rank_one_records,
        },
        "rank_two": {
            "accepted_count": sum(
                record["accepted_pair_count"]
                for record in rank_two_records
            ),
            "by_coefficient_map_kernel": rank_two_records,
        },
        "accepted_total": len(accepted_records),
        "accepted_perturbations": accepted_records,
        "gold_stabilizer_orbits": orbit_records,
        "gold_stabilizer_orbit_count": len(orbit_records),
        "all_orbit_sizes_divide_group_order": all(
            GROUP_ORDER % record["orbit_size"] == 0
            for record in orbit_records
        ),
        "all_signatures_match_archived_classes": all(
            bool(record["signature_compatible_classes"])
            for record in accepted_records
        ),
        "classification_warning": (
            "Orthoderivative-signature agreement is a necessary CCZ filter, "
            "not an exact CCZ/EA proof for candidates not previously exact-tested."
        ),
    }


def find_rank_one_record(locus: dict[str, Any], mask: int) -> dict[str, Any]:
    return next(
        record for record in locus["rank_one"]["by_kernel"]
        if record["kernel_mask"] == mask
    )


def q0_q1_bridge(x3: dict[str, Any], x9: dict[str, Any]) -> dict[str, Any]:
    x3_q0 = find_rank_one_record(x3, 1)["accepted_coefficients"]
    x3_q1 = find_rank_one_record(x3, 2)["accepted_coefficients"]
    x9_q0 = find_rank_one_record(x9, 1)["accepted_coefficients"]
    x9_q1 = find_rank_one_record(x9, 2)["accepted_coefficients"]

    x3_plane = next(
        record for record in x3["rank_two"]["by_coefficient_map_kernel"]
        if record["coefficient_map_kernel"] == 4
    )
    x9_plane = next(
        record for record in x9["rank_two"]["by_coefficient_map_kernel"]
        if record["coefficient_map_kernel"] == 4
    )
    pairs3 = [
        (record["coefficient_1"], record["coefficient_2"])
        for record in x3_plane["accepted_pairs"]
    ]
    pairs9 = [
        (record["coefficient_1"], record["coefficient_2"])
        for record in x9_plane["accepted_pairs"]
    ]

    alpha = 0x24
    beta = 0x3D
    lambda_value = 0xBD
    lambda_squared = gf_mult(lambda_value, lambda_value)
    slope3 = gf_mult(beta, gf_inv(alpha))
    coefficient9_q0 = gf_pow(alpha, 3)
    coefficient9_q1 = gf_pow(beta, 3)
    slope9 = gf_pow(slope3, 3)

    all_cube_roots_q0 = sorted(
        value for value in range(1, N)
        if gf_pow(value, 3) == coefficient9_q0
    )
    all_cube_roots_q1 = sorted(
        value for value in range(1, N)
        if gf_pow(value, 3) == coefficient9_q1
    )

    expected_pairs3 = sorted([
        (value, gf_mult(slope3, value))
        for value in all_cube_roots_q0
        if absolute_trace(value) == 1
    ])

    return {
        "field_equations": {
            "x3_q0": (
                "mu is accepted iff mu^3=0xf2 (mu nonzero); "
                "the three roots are 0x24*GF(4)^*."
            ),
            "x3_q1": (
                "nu is accepted iff nu^3=0x6f (nu nonzero); "
                "the three roots are 0x3d*GF(4)^*."
            ),
            "x3_rank_two_q0_q1": (
                "(mu,nu) is accepted iff mu^3=0xf2, "
                "nu=0xed*mu, and Tr(mu)=1."
            ),
            "x9_rank_two_q0_q1": (
                "The unique accepted pair is (0xf2,0x6f), "
                "with slope 0xb0."
            ),
        },
        "lambda": lambda_value,
        "lambda_hex": f"0x{lambda_value:02x}",
        "lambda_order": gf_order(lambda_value),
        "lambda_squared_hex": f"0x{lambda_squared:02x}",
        "x3_q0_coefficients_hex": [f"0x{x:02x}" for x in x3_q0],
        "x3_q1_coefficients_hex": [f"0x{x:02x}" for x in x3_q1],
        "x9_q0_coefficients_hex": [f"0x{x:02x}" for x in x9_q0],
        "x9_q1_coefficients_hex": [f"0x{x:02x}" for x in x9_q1],
        "q0_cube_target_hex": f"0x{coefficient9_q0:02x}",
        "q1_cube_target_hex": f"0x{coefficient9_q1:02x}",
        "x3_q0_is_full_cube_fibre": sorted(x3_q0) == all_cube_roots_q0,
        "x3_q1_is_full_cube_fibre": sorted(x3_q1) == all_cube_roots_q1,
        "x3_rank_two_pairs_hex": [
            [f"0x{a:02x}", f"0x{b:02x}"] for a, b in pairs3
        ],
        "x3_rank_two_expected_by_cube_slope_trace_hex": [
            [f"0x{a:02x}", f"0x{b:02x}"]
            for a, b in expected_pairs3
        ],
        "x3_rank_two_exact_trace_characterisation": (
            sorted(pairs3) == expected_pairs3
        ),
        "excluded_third_cube_lift": {
            "q0_coefficient_hex": "0xac",
            "q1_coefficient_hex": "0xb2",
            "common_multiplier_hex": f"0x{lambda_squared:02x}",
            "absolute_trace_q0_coefficient": absolute_trace(0xAC),
            "reason_in_closed_characterisation": "Tr(0xac)=0",
        },
        "x3_slope_hex": f"0x{slope3:02x}",
        "x3_slope_order": gf_order(slope3),
        "x9_slope_hex": f"0x{slope9:02x}",
        "x9_slope_order": gf_order(slope9),
        "slope_cube_relation": slope9 == gf_pow(slope3, 3),
        "x9_rank_two_pairs_hex": [
            [f"0x{a:02x}", f"0x{b:02x}"] for a, b in pairs9
        ],
        "componentwise_cube_images_of_x3_pairs_hex": sorted({
            (f"0x{gf_pow(a, 3):02x}", f"0x{gf_pow(b, 3):02x}")
            for a, b in pairs3
        }),
        "componentwise_cube_collapses_to_unique_x9_pair": (
            { (gf_pow(a, 3), gf_pow(b, 3)) for a, b in pairs3 }
            == set(pairs9)
        ),
        "interpretation": (
            "The x^3 q0/q1 APN locus is a cubic lift of the x^9 locus. "
            "The cube map has kernel GF(4)^*, explaining the three rank-one "
            "lifts.  In rank two, a trace-one condition retains two of the "
            "three synchronized lifts; componentwise cubing identifies both "
            "with the unique x^9 pair and sends the fifth-root slope "
            "0xed to 0xb0."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kernels",
        type=Path,
        default=Path("results/evidence/canonical_switching_kernels.json"),
    )
    parser.add_argument(
        "--signatures",
        type=Path,
        default=Path(
            "results/evidence/scalar_family_completion/"
            "scalar_family_signature_classification.json"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/canonical_low_rank_locus.json"),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    resolve = lambda path: path if path.is_absolute() else root / path
    kernels_path = resolve(args.kernels)
    signatures_path = resolve(args.signatures)
    output_path = resolve(args.output)

    payload = json.loads(kernels_path.read_text(encoding="utf-8"))
    base_truths = [
        [
            q_value(
                [tuple(edge) for edge in payload["kernels"][name]["expanded_terms"]],
                x,
            )
            for x in range(N)
        ]
        for name in KERNEL_NAMES
    ]
    kernel_truths = {
        mask: mask_truth(base_truths, mask)
        for mask in range(1, 8)
    }
    derivatives = {
        mask: boolean_derivatives(truth)
        for mask, truth in kernel_truths.items()
    }
    signature_classes = load_signature_classes(signatures_path)

    x3 = enumerate_locus(3, kernel_truths, derivatives, signature_classes)
    x9 = enumerate_locus(9, kernel_truths, derivatives, signature_classes)
    bridge = q0_q1_bridge(x3, x9)

    report = {
        "schema": "canonical-low-rank-apn-locus-v1",
        "date": "2026-08-02",
        "scope": (
            "Complete coefficient-rank <=2 APN enumeration in the canonical "
            "24-dimensional q0,q1,q2 tensor space around Gold x^3 and x^9"
        ),
        "method": {
            "forbidden_set_criterion": (
                "For fixed Boolean kernels, Delta=sum v_i q_i is APN iff "
                "U(y) avoids the exact forbidden set B_y for every y!=0."
            ),
            "rank_one_test": "v not in B_1",
            "rank_two_test": (
                "For Delta=v1 Q1+v2 Q2: v1 not in B_01, "
                "v2 not in B_10, and v1+v2 not in B_11."
            ),
            "search_status": (
                "Exhaustive only inside the rank<=2 determinantal subset of "
                "GF(256) tensor span(q0,q1,q2); no broad APN/Gröbner search."
            ),
        },
        "provenance": {
            "kernels_file": str(kernels_path.relative_to(root)),
            "kernels_sha256": sbox_sha256(kernels_path.read_bytes()),
            "archived_signature_file": str(signatures_path.relative_to(root)),
            "archived_signature_sha256": sbox_sha256(
                signatures_path.read_bytes()
            ),
        },
        "loci": {"x3": x3, "x9": x9},
        "q0_q1_cubic_trace_bridge": bridge,
        "main_results": [
            "Gold x^3 has exactly 7 rank-one and 2 rank-two nonzero APN perturbations in the full q0,q1,q2 kernel tensor space.",
            "Gold x^9 has exactly 7 rank-one and 1 rank-two nonzero APN perturbations in the same space.",
            "Every rank-two solution lies in the q0,q1 plane; all other rank-two coefficient-map kernels are empty.",
            "The q0,q1 x^3 locus admits an exact cube-fibre plus trace-one description, and componentwise cubing maps it to the unique x^9 rank-two point.",
            "Accepted perturbations split into 5 Gold-stabilizer orbits for x^3 and 6 for x^9, showing that one ordinary signature/class can have several distinct marked local secant orbits.",
        ],
        "validation": {
            "x3_counts": {
                "rank_one": x3["rank_one"]["accepted_count"],
                "rank_two": x3["rank_two"]["accepted_count"],
                "total": x3["accepted_total"],
            },
            "x9_counts": {
                "rank_one": x9["rank_one"]["accepted_count"],
                "rank_two": x9["rank_two"]["accepted_count"],
                "total": x9["accepted_total"],
            },
            "all_rank_two_solutions_are_q0_q1": all(
                record["accepted_pair_count"] == 0
                or record["coefficient_map_kernel"] == 4
                for locus in (x3, x9)
                for record in locus["rank_two"]["by_coefficient_map_kernel"]
            ),
            "cube_trace_bridge_exact": all([
                bridge["x3_q0_is_full_cube_fibre"],
                bridge["x3_q1_is_full_cube_fibre"],
                bridge["x3_rank_two_exact_trace_characterisation"],
                bridge["componentwise_cube_collapses_to_unique_x9_pair"],
                bridge["slope_cube_relation"],
            ]),
            "all_orbit_sizes_valid": (
                x3["all_orbit_sizes_divide_group_order"]
                and x9["all_orbit_sizes_divide_group_order"]
            ),
            "all_signatures_covered_by_archived_classes": (
                x3["all_signatures_match_archived_classes"]
                and x9["all_signatures_match_archived_classes"]
            ),
        },
    }

    expected = {
        "x3": (7, 2, 9, 5),
        "x9": (7, 1, 8, 6),
    }
    observed = {
        "x3": (
            x3["rank_one"]["accepted_count"],
            x3["rank_two"]["accepted_count"],
            x3["accepted_total"],
            x3["gold_stabilizer_orbit_count"],
        ),
        "x9": (
            x9["rank_one"]["accepted_count"],
            x9["rank_two"]["accepted_count"],
            x9["accepted_total"],
            x9["gold_stabilizer_orbit_count"],
        ),
    }
    if observed != expected:
        raise RuntimeError(f"unexpected canonical locus counts: {observed}")
    if not all(report["validation"].values()):
        raise RuntimeError("one or more validation checks failed")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "x3": observed["x3"],
        "x9": observed["x9"],
        "bridge": report["validation"]["cube_trace_bridge_exact"],
        "all_rank_two_q0_q1": report["validation"]["all_rank_two_solutions_are_q0_q1"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
