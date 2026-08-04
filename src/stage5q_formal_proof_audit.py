#!/usr/bin/env python3
r"""Stage 5Q: formal proof audit and corrected permutation-condition screen.

This audit does not introduce a new APN search.  It checks the logical and
computational dependencies of Stages 5E--5P, corrects the previously too narrow
permutation necessary-condition screen, and records which statements are pure,
computer-assisted, externally classified, or still require manuscript work.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Iterable, Any

from stage5g_trace_coordinate_theory import (
    absolute_trace,
    gf_mult,
    gf_pow,
)
from stage5o_small_dimension_structure import (
    BinaryField,
    build_q_and_lut,
    component_data,
    first_irreducible_polynomial,
    is_apn_quadratic,
    absolute_trace as small_absolute_trace,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binary_rank(rows: Iterable[int], n: int) -> int:
    basis = [0] * n
    rank = 0
    for value in rows:
        x = value
        while x:
            p = x.bit_length() - 1
            if basis[p]:
                x ^= basis[p]
            else:
                basis[p] = x
                rank += 1
                break
    return rank


def rref_subspaces(k: int, n: int):
    """Yield a canonical RREF basis and a bit-mask of nonzero vectors."""
    total = 0
    for pivots in itertools.combinations(range(n), k):
        nonpivots = [c for c in range(n) if c not in pivots]
        allowed = [(i, c) for c in nonpivots for i, p in enumerate(pivots) if p < c]
        for assignment in range(1 << len(allowed)):
            rows = [1 << p for p in pivots]
            for bit, (i, c) in enumerate(allowed):
                if (assignment >> bit) & 1:
                    rows[i] |= 1 << c
            assert binary_rank(rows, n) == k
            span_mask = 0
            for subset in range(1, 1 << k):
                x = 0
                for i, row in enumerate(rows):
                    if (subset >> i) & 1:
                        x ^= row
                span_mask |= 1 << x
            total += 1
            yield tuple(rows), span_mask
    expected = {
        (1, 6): 63, (2, 6): 651, (3, 6): 1395, (4, 6): 651, (5, 6): 63,
        (1, 8): 255, (2, 8): 10795, (3, 8): 97155,
        (4, 8): 200787, (5, 8): 97155, (6, 8): 10795, (7, 8): 255,
    }.get((k, n))
    if expected is not None:
        assert total == expected


def nb_mask(nb: set[int]) -> int:
    result = 0
    for x in nb:
        result |= 1 << x
    return result


def full_complementary_screen(names: list[str], nb_sets: list[set[int]], n: int) -> dict[str, Any]:
    """Test the actual theorem: complementary dimensions may be unequal."""
    # At least one member of a complementary pair has dimension >= ceil(n/2).
    dims = list(range(1, n))
    masks = [nb_mask(s) for s in nb_sets]
    contained: dict[int, list[list[tuple[tuple[int, ...], int]]]] = {
        d: [[] for _ in names] for d in dims
    }
    total_by_dim: dict[int, int] = {}
    for d in dims:
        total = 0
        for basis, span in rref_subspaces(d, n):
            total += 1
            for i, mask in enumerate(masks):
                if span & ~mask == 0:
                    contained[d][i].append((basis, span))
        total_by_dim[d] = total

    records = []
    for idx, name in enumerate(names):
        split_records = []
        total_pairs = 0
        for d in range(1, n // 2 + 1):
            e = n - d
            left = contained[d][idx]
            right = contained[e][idx]
            count = 0
            sample = None
            if d == e:
                for i in range(len(left)):
                    for j in range(i + 1, len(right)):
                        if left[i][1] & right[j][1] == 0:
                            count += 1
                            if sample is None:
                                sample = (left[i][0], right[j][0])
            else:
                for a in left:
                    for b in right:
                        if a[1] & b[1] == 0:
                            count += 1
                            if sample is None:
                                sample = (a[0], b[0])
            total_pairs += count
            split_records.append({
                "dimension_split": f"{d}+{e}",
                "contained_first": len(left),
                "contained_second": len(right),
                "complementary_pairs": count,
                "sample_pair_bases_hex": None if sample is None else [
                    [f"0x{x:0{(n+3)//4}x}" for x in sample[0]],
                    [f"0x{x:0{(n+3)//4}x}" for x in sample[1]],
                ],
            })
        records.append({
            "id": name,
            "non_bent_component_count": len(nb_sets[idx]),
            "contained_subspaces_by_dimension": {
                str(d): len(contained[d][idx]) for d in dims
            },
            "dimension_splits": split_records,
            "total_complementary_pairs": total_pairs,
            "passes_actual_necessary_condition": total_pairs > 0,
        })
    return {
        "ambient_dimension": n,
        "theorem_checked": (
            "NB(F) union {0} contains subspaces V,W of arbitrary complementary "
            "dimensions with V direct_sum W = F_2^n"
        ),
        "subspaces_enumerated_by_dimension": {str(k): v for k, v in total_by_dim.items()},
        "records": records,
        "all_fail": all(not r["passes_actual_necessary_condition"] for r in records),
    }


def q_c(c: int, x: int, lam: int) -> int:
    return absolute_trace(gf_mult(c, x)) & absolute_trace(gf_mult(gf_mult(lam, c), x))


def is_boolean_linear(truth: list[int], n: int) -> bool:
    if truth[0] != 0:
        return False
    columns = [truth[1 << i] for i in range(n)]
    for x in range(1 << n):
        value = 0
        for i, bit in enumerate(columns):
            if (x >> i) & 1:
                value ^= bit
        if value != truth[x]:
            return False
    return True


def projective_well_definedness() -> dict[str, Any]:
    E = sorted(x for x in range(256) if gf_pow(x, 4) == x)
    E_star = [x for x in E if x]
    L = sorted(x for x in range(256) if gf_pow(x, 16) == x)
    lam = next(x for x in E if x not in (0, 1))
    checked = 0
    for c in L:
        if c == 0:
            continue
        for alpha in E_star:
            d = [q_c(gf_mult(alpha, c), x, lam) ^ q_c(c, x, lam) for x in range(256)]
            assert is_boolean_linear(d, 8)
            assert gf_pow(gf_mult(alpha, c), 252) == gf_pow(c, 252)  # inverse cube
            assert gf_pow(gf_mult(alpha, c), 6) == gf_pow(c, 6)
            checked += 1
    return {
        "pairs_checked": checked,
        "q_alpha_c_plus_q_c_is_linear": True,
        "rho_inverse_cube_representative_independent": True,
        "rho_sixth_power_representative_independent": True,
        "formal_correction": (
            "q_[c] is an equivalence class modulo linear Boolean functions; the "
            "resulting vectorial switching is representative-independent only up "
            "to addition of a linear map, which preserves APN."
        ),
    }


def exact_small_dimension_checks() -> dict[str, Any]:
    f4 = BinaryField(4, first_irreducible_polynomial(4))
    n4 = []
    for theta in range(1, 16):
        _, lut = build_q_and_lut(f4, theta)
        actual = is_apn_quadratic(lut)
        predicted = small_absolute_trace(f4, theta) == 0
        assert actual == predicted
        n4.append({"theta_hex": f"0x{theta:x}", "APN": actual, "predicted": predicted})

    f6 = BinaryField(6, first_irreducible_polynomial(6))
    n6 = []
    for theta in range(1, 64):
        _, lut = build_q_and_lut(f6, theta)
        actual = is_apn_quadratic(lut)
        predicted = f6.pow(theta, 6) ^ f6.pow(theta, 3) ^ 1 == 0
        assert actual == predicted
        n6.append({"theta_hex": f"0x{theta:02x}", "APN": actual, "predicted": predicted})
    return {
        "n4_all_15_nonzero_coefficients_checked": True,
        "n4_safe_count": sum(r["APN"] for r in n4),
        "n6_all_63_nonzero_coefficients_checked": True,
        "n6_safe_count": sum(r["APN"] for r in n6),
        "all_predictions_match_direct_APN": True,
    }


def build_n6_nb_sets(root: Path) -> tuple[list[str], list[set[int]], dict[str, Any]]:
    small = json.loads((root / "results/current/small_dimension_trace_switch_structure.json").read_text())
    f6 = BinaryField(6, first_irreducible_polynomial(6))
    names: list[str] = []
    sets: list[set[int]] = []
    for rec in small["dimension_6"]["direct_records"]:
        theta = int(rec["theta_hex"], 16)
        _, lut = build_q_and_lut(f6, theta)
        _, nb, _ = component_data(lut)
        names.append(rec["theta_hex"])
        sets.append(nb)

    kim_field = BinaryField(6, 0x5B)
    alpha = 0x02
    kim = [
        kim_field.pow(x, 3)
        ^ kim_field.pow(x, 10)
        ^ kim_field.mul(alpha, kim_field.pow(x, 24))
        for x in range(64)
    ]
    _, kim_nb, _ = component_data(kim)
    names.append("Kim-positive-control")
    sets.append(kim_nb)
    return names, sets, {"family_count": 6, "positive_control_included": True}


def load_n8_nb_sets(root: Path) -> tuple[list[str], list[set[int]], dict[str, Any]]:
    projective = json.loads((root / "results/current/projective_switch_permutation_screen.json").read_text())
    archived = json.loads((root / "results/current/novel_class_permutation_screen.json").read_text())
    names: list[str] = []
    sets: list[set[int]] = []
    for rec in projective["records"]:
        names.append(f"projective-{rec['rho0_hex']}-{rec['eta_hex']}")
        sets.append({int(x, 16) for x in rec["non_bent_components_hex"]})
    for rec in archived["records"]:
        names.append(rec["id"])
        sets.append({int(x, 16) for x in rec["non_bent_components_hex"]})
    return names, sets, {"projective_count": 8, "archived_count": 8}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/current/formal_proof_audit.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.output if args.output.is_absolute() else root / args.output

    small_checks = exact_small_dimension_checks()
    projective_checks = projective_well_definedness()

    names6, nb6, meta6 = build_n6_nb_sets(root)
    screen6 = full_complementary_screen(names6, nb6, 6)
    # Six trace switchings fail, Kim positive control passes.
    assert all(not r["passes_actual_necessary_condition"] for r in screen6["records"][:-1])
    assert screen6["records"][-1]["passes_actual_necessary_condition"]
    assert screen6["records"][-1]["total_complementary_pairs"] == 3

    names8, nb8, meta8 = load_n8_nb_sets(root)
    screen8 = full_complementary_screen(names8, nb8, 8)
    assert screen8["all_fail"]

    scan = json.loads((root / "results/current/trace_switch_dimension_scan.json").read_text())
    rows = {row["n"]: row for row in scan["dimensions"]}
    assert rows[10]["minimum_positive_representation_count"] == 46
    assert rows[12]["minimum_positive_representation_count"] == 208
    total_sqrt_constant = 2 + 6 * 6 + 9 * 4
    assert total_sqrt_constant == 74
    assert (2**14 - 2 - 74 * 2**7) > 0

    claim_status = [
        {
            "claim": "General low-rank derivative-incidence criterion",
            "status": "proved after normalization correction",
            "required_wording": (
                "Use L_aF(x)=F(x+a)+F(x)+F(a)+F(0), or explicitly assume "
                "F(0)=0 and q_i(0)=0.  The criterion works for every normalized "
                "quadratic APN centre, not only Gold."
            ),
        },
        {
            "claim": "n=4 coefficient classification and Gold EA reduction",
            "status": "proved with exhaustive finite necessity check",
            "evidence": "all 15 nonzero coefficients directly checked; closed-form output map proves the safe EA reduction",
        },
        {
            "claim": "n=6 order-nine coefficient classification",
            "status": "computer-assisted exact theorem",
            "evidence": "256-case Kummer certificate over GF(4), independently matched by all 63 direct APN checks",
        },
        {
            "claim": "n=8 rank-one coefficient set GF(4)^*",
            "status": "computer-assisted exact theorem",
            "evidence": "symbolic Bezout certificate for sufficiency; 4096-section exact certificate for necessity",
        },
        {
            "claim": "n=8 projective rank-two selector",
            "status": "computer-assisted exact theorem",
            "evidence": "linear cases excluded symbolically; involutions reduced to a 16-pair/four-row table; projective transport checked",
        },
        {
            "claim": "nonexistence for every even n>=10",
            "status": "proof architecture valid; manuscript expansion required",
            "evidence": "Artin-Schreier character bound gives constant 74 from n>=14; exact convolution closes n=10,12",
        },
        {
            "claim": "permutation-class exclusion",
            "status": "corrected and revalidated",
            "required_wording": (
                "The theorem permits arbitrary complementary dimensions.  The "
                "new full scan checks every split; conclusions remain negative "
                "for all six n=6 trace switchings and all sixteen n=8 representatives."
            ),
        },
        {
            "claim": "n=6 identification with Banff class 2",
            "status": "external-classification dependent",
            "required_wording": (
                "State that the orthoderivative differential spectrum equals the "
                "unique class-2 row in the published exhaustive Banff table; cite "
                "that classification."
            ),
        },
        {
            "claim": "eight n=8 projective switchings are eight distinct classes",
            "status": "not proved and must not be claimed",
            "required_wording": "They are eight marked parameter points in two Frobenius orbits; exact EA/CCZ partition remains separate.",
        },
    ]

    tracked = [
        "results/current/low_rank_derivative_criterion.json",
        "results/current/trace_coordinate_theory.json",
        "results/current/rank_two_selector_proof.json",
        "results/current/trace_switch_dimension_scan.json",
        "results/current/projective_displacement_theorem.json",
        "results/current/uniform_trace_switch_nonexistence.json",
        "results/current/coordinate_free_projective_theory.json",
        "results/current/small_dimension_trace_switch_structure.json",
        "results/current/kummer_n6_trace_switch_certificate.json",
    ]
    provenance = {p: sha256(root / p) for p in tracked}

    result = {
        "schema": "formal-proof-audit-stage5-v1",
        "date": "2026-08-02",
        "scope": "Stages 5E--5P, including the corrected CCZ-permutation necessary-condition screen",
        "normalization_and_projective_checks": projective_checks,
        "small_dimension_exact_checks": small_checks,
        "permutation_condition_correction": {
            "old_incorrect_restriction": "only n/2+n/2 complementary subspaces were searched",
            "actual_condition": "two subspaces of arbitrary dimensions summing to n",
            "n6_metadata": meta6,
            "n6_full_screen": screen6,
            "n8_metadata": meta8,
            "n8_full_screen": screen8,
            "conclusions_survive": True,
        },
        "uniform_bound_audit": {
            "Hasse_error_constant": 2,
            "axis_character_pairs": 6,
            "axis_bound_constant_each": 6,
            "off_axis_character_pairs": 9,
            "off_axis_bound_constant_each": 4,
            "total_sqrt_constant": total_sqrt_constant,
            "first_even_dimension_with_positive_uniform_bound": 14,
            "finite_bridge_minima": {"10": 46, "12": 208},
            "status": "numerics verified; final paper must state and cite the Artin-Schreier sum theorem and prove the pole counts",
        },
        "claim_status": claim_status,
        "publication_readiness": {
            "ready_as_pure_theorem": ["general low-rank criterion after wording correction"],
            "ready_as_computer_assisted_theorem": [
                "n=4 exact classification", "n=6 order-nine theorem",
                "n=8 rank-one theorem", "n=8 projective rank-two theorem",
            ],
            "requires_full_written_proof": ["uniform n>=10 Artin-Schreier/Weil bound"],
            "requires_external_citation_or_exact_equivalence": ["Banff class-2 identification"],
            "must_be_removed_or_weakened": [
                "equal-half formulation of the permutation condition",
                "claim that eight marked switchings are eight inequivalent classes",
                "unqualified use of P(L/E) instead of P_E(L)",
            ],
        },
        "provenance_sha256": provenance,
        "validation": {
            "all_assertions_passed": True,
            "n6_negative_conclusion_survives_corrected_screen": True,
            "Kim_positive_control_passes_corrected_screen": True,
            "n8_negative_conclusions_survive_corrected_screen": True,
            "all_small_dimension_coefficients_checked": True,
            "projective_representative_independence_checked": True,
        },
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "sha256": sha256(output),
        "n6_trace_switches_pass_actual_condition": sum(
            r["passes_actual_necessary_condition"] for r in screen6["records"][:-1]
        ),
        "kim_complementary_pairs": screen6["records"][-1]["total_complementary_pairs"],
        "n8_representatives_pass_actual_condition": sum(
            r["passes_actual_necessary_condition"] for r in screen8["records"]
        ),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
