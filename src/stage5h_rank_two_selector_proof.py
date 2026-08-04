#!/usr/bin/env python3
"""Stage 5H: compact analytic proof of the rank-two trace selector.

This script replaces the exhaustive 10,795-subspace obstruction certificate
from Stage 5G by a calculation in the tower

    GF(2) subset GF(4) subset GF(16) subset GF(256).

The six transition maps in GL(2,2) split as

    y -> alpha*y       (identity and two order-three maps),
    y -> alpha*y^2     (the three involutions),

with alpha in GF(4)^*.  The first family cannot produce a section product in
GF(16).  The second family reduces to a 16-pair calculation in GF(16),
independent of alpha, and meets the target line only for eta=lambda^2, with
four sections per involution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from stage5g_trace_coordinate_theory import (
    gf_inv,
    gf_mult,
    gf_order,
    gf_pow,
    relative_trace_to_f4,
)

N = 256
R = 0xB0
LAMBDA = 0xBD


def tr_f16_f2(x: int) -> int:
    value = x ^ gf_pow(x, 2) ^ gf_pow(x, 4) ^ gf_pow(x, 8)
    if value not in (0, 1):
        raise AssertionError(f"not an F16/F2 trace value: 0x{value:02x}")
    return value


def f16_elements() -> list[int]:
    return [x for x in range(N) if gf_pow(x, 16) == x]


def f4_elements() -> list[int]:
    return [x for x in range(N) if gf_pow(x, 4) == x]


def expr_map(r: int) -> dict[int, str]:
    """Unique expressions in the F2-basis (1,r,r^2,r^3) of GF(16)."""
    basis = [1, r, gf_pow(r, 2), gf_pow(r, 3)]
    names = ["1", "r", "r^2", "r^3"]
    result: dict[int, str] = {}
    for mask in range(16):
        value = 0
        terms = []
        for i, item in enumerate(basis):
            if (mask >> i) & 1:
                value ^= item
                terms.append(names[i])
        if value in result:
            raise AssertionError("basis is not independent")
        result[value] = "+".join(terms) if terms else "0"
    return result


def map_order(kind: str, alpha: int) -> int:
    if kind == "linear":
        return 1 if alpha == 1 else 3
    if kind == "semilinear":
        return 2
    raise ValueError(kind)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage5g",
        type=Path,
        default=Path("results/current/trace_coordinate_theory.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/rank_two_selector_proof.json"),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    resolve = lambda p: p if p.is_absolute() else root / p
    stage5g_path = resolve(args.stage5g)
    output_path = resolve(args.output)
    stage5g = json.loads(stage5g_path.read_text(encoding="utf-8"))

    r = R
    lam = LAMBDA
    r2, r3, r4 = (gf_pow(r, i) for i in (2, 3, 4))
    lam2 = gf_pow(lam, 2)
    assert gf_order(r) == 5
    assert gf_order(lam) == 3
    assert lam == (r2 ^ r3)
    assert lam2 == (r ^ r4)
    assert r4 ^ r == lam2
    assert gf_pow(r, 5) == 1

    f16 = f16_elements()
    f4 = f4_elements()
    f4_star = [x for x in f4 if x]
    assert set(f4_star) == {1, lam, lam2}
    expressions = expr_map(r)

    # Choose an Artin-Schreier generator w for GF(256)/GF(16).
    w_candidates = [
        x for x in range(N)
        if (gf_pow(x, 2) ^ x) == r and (gf_pow(x, 16) ^ x) == 1
    ]
    assert len(w_candidates) == 2
    w = min(w_candidates)

    # For x=A+s*w, Tr_K/F16(x)=s.  The pair
    # (T(x),T(r*x)) determines s uniquely.
    d = r ^ r4
    assert d == lam2
    d_inv = gf_inv(d)
    assert d_inv == lam

    def relative_trace_pair_from_s(s: int) -> tuple[int, int]:
        return (
            s ^ gf_pow(s, 4),
            gf_mult(r, s) ^ gf_mult(r4, gf_pow(s, 4)),
        )

    pair_to_s: dict[tuple[int, int], int] = {}
    for s in f16:
        pair = relative_trace_pair_from_s(s)
        assert pair not in pair_to_s
        pair_to_s[pair] = s
    assert len(pair_to_s) == 16

    # GL(2,2) = {alpha*y} union {alpha*y^2}, alpha in GF(4)^*.
    transition_maps: list[dict[str, Any]] = []
    for kind in ("linear", "semilinear"):
        for alpha in f4_star:
            if kind == "linear":
                image_1 = alpha
                image_lambda = gf_mult(alpha, lam)
                u_coeff, v_coeff = alpha, 0
            else:
                image_1 = alpha
                image_lambda = gf_mult(alpha, lam2)
                u_coeff, v_coeff = 0, alpha
            s_a = pair_to_s[(1, u_coeff)]
            s_b = pair_to_s[(0, v_coeff)]
            # Closed forms obtained from the two trace equations.
            assert s_a == gf_mult(u_coeff ^ r4, d_inv)
            assert s_b == gf_mult(v_coeff, d_inv)
            transition_maps.append({
                "kind": kind,
                "alpha_hex": f"0x{alpha:02x}",
                "images_of_1_lambda_hex": [
                    f"0x{image_1:02x}", f"0x{image_lambda:02x}"
                ],
                "order": map_order(kind, alpha),
                "Tr_K_F16_a_hex": f"0x{s_a:02x}",
                "Tr_K_F16_b_hex": f"0x{s_b:02x}",
            })

    assert Counter(row["order"] for row in transition_maps) == {1: 1, 2: 3, 3: 2}

    # Linear transitions: b is in GF(16).  Write a=s(u+w).  A product
    # a^3+b^3 can lie in GF(16) only if g(u)=0, where
    # g(z)=z^2+z+1+r.  But Tr_F16/F2(1+r)=1, so g has no root.
    assert tr_f16_f2(1 ^ r) == 1

    def g(z: int) -> int:
        return gf_pow(z, 2) ^ z ^ 1 ^ r

    assert all(g(z) != 0 for z in f16)

    # Semilinear transitions: s=r^4/d, h=alpha/d.  Because alpha^3=1,
    # kappa=(h/s)^3=r^3 is independent of alpha and s^3=r^2.
    s = gf_mult(r4, d_inv)
    assert gf_pow(s, 3) == r2
    kappa = r3
    semilinear_reductions = []
    for alpha in f4_star:
        h = gf_mult(alpha, d_inv)
        assert gf_pow(gf_mult(h, gf_inv(s)), 3) == kappa
        semilinear_reductions.append({
            "alpha_hex": f"0x{alpha:02x}",
            "s_hex": f"0x{s:02x}",
            "h_hex": f"0x{h:02x}",
            "h_over_s_hex": f"0x{gf_mult(h, gf_inv(s)):02x}",
            "kappa_hex": f"0x{kappa:02x}",
        })

    # If a=s(u+w), b=h(v+w), then
    #   a^3+b^3 in GF(16) <=> g(u)=r^3*g(v),
    # and in that case
    #   (a^3+b^3)/r^2 = (u+v)g(u)+(1+r^3).
    reduced_pairs = []
    for u in f16:
        for v in f16:
            if g(u) != gf_mult(kappa, g(v)):
                continue
            normalized_product = gf_mult(u ^ v, g(u)) ^ (1 ^ kappa)
            reduced_pairs.append((u, v, g(v), g(u), normalized_product))
    assert len(reduced_pairs) == 16

    y_values = sorted({row[2] for row in reduced_pairs})
    expected_y = {r4, 1 ^ r2, r3, r}
    assert set(y_values) == expected_y

    table = []
    normalized_product_counter: Counter[int] = Counter()
    for y in y_values:
        rows = [row for row in reduced_pairs if row[2] == y]
        assert len(rows) == 4
        u_roots = sorted({row[0] for row in rows})
        v_roots = sorted({row[1] for row in rows})
        products = Counter(row[4] for row in rows)
        assert len(u_roots) == len(v_roots) == 2
        assert u_roots[0] ^ u_roots[1] == 1
        assert v_roots[0] ^ v_roots[1] == 1
        assert sorted(products.values()) == [2, 2]
        normalized_product_counter.update(products)
        table.append({
            "y_equals_g_v_hex": f"0x{y:02x}",
            "y_expression": expressions[y],
            "r3_y_equals_g_u_hex": f"0x{gf_mult(r3, y):02x}",
            "r3_y_expression": expressions[gf_mult(r3, y)],
            "u_roots_hex": [f"0x{x:02x}" for x in u_roots],
            "u_roots_expressions": [expressions[x] for x in u_roots],
            "v_roots_hex": [f"0x{x:02x}" for x in v_roots],
            "v_roots_expressions": [expressions[x] for x in v_roots],
            "normalized_products": [
                {
                    "hex": f"0x{x:02x}",
                    "expression": expressions[x],
                    "multiplicity": products[x],
                }
                for x in sorted(products)
            ],
        })

    # The target line after division by s^3=r^2 is eta*(1+r^3).
    targets = {}
    target_counts = {}
    for eta in f4_star:
        target = gf_mult(eta, 1 ^ r3)
        targets[eta] = target
        target_counts[eta] = normalized_product_counter[target]
    assert target_counts == {1: 0, lam: 0, lam2: 4}

    # Match the Stage 5G full-subspace certificate exactly.
    full_rows = stage5g["rank_two_q0_q1_theory"][
        "trace_transition_map_classification"
    ]
    expected_full = sorted(
        (
            tuple(row["images_of_1_lambda_hex"]),
            row["order"],
            row["target_line_witness_counts"],
        )
        for row in full_rows
    )
    compact_full = []
    for row in transition_maps:
        compact_full.append((
            tuple(row["images_of_1_lambda_hex"]),
            row["order"],
            {
                f"0x{eta:02x}": (
                    target_counts[eta] if row["order"] == 2 else 0
                )
                for eta in sorted(f4_star)
            },
        ))
    assert sorted(compact_full) == expected_full

    result = {
        "schema": "rank-two-selector-proof-v1",
        "date": "2026-08-02",
        "statement": {
            "family": "x^3 + eta*(Q_1+r^2*Q_r), eta in GF(4)^*",
            "theorem": "The family is APN iff Tr_GF4/F2(lambda*eta)=1.",
            "accepted_eta_hex": [f"0x{x:02x}" for x in (1, lam)],
            "excluded_eta_hex": f"0x{lam2:02x}",
            "excluded_witness_count": 12,
            "witness_distribution": "four sections for each of the three involutions in GL(2,2)",
        },
        "field_identities": {
            "r_hex": f"0x{r:02x}",
            "r_order": gf_order(r),
            "lambda_hex": f"0x{lam:02x}",
            "lambda_expression": "r^2+r^3",
            "lambda_squared_hex": f"0x{lam2:02x}",
            "lambda_squared_expression": "r+r^4",
            "artin_schreier_w_hex": f"0x{w:02x}",
            "w_equations": "w^2+w=r and w^16+w=1",
            "GF16_basis": ["1", "r", "r^2", "r^3"],
        },
        "GL2_F2_decomposition": {
            "description": "GL(2,2)={y->alpha*y} disjoint_union {y->alpha*y^2}, alpha in GF(4)^*.",
            "maps": transition_maps,
        },
        "linear_transition_proof": {
            "maps": "y->alpha*y (orders 1,3,3)",
            "argument": "For b in GF(16) and a=s(u+w), a^3+b^3 in GF(16) would require g(u)=0, g(u)=u^2+u+1+r. Since Tr_GF16/F2(1+r)=1, g has no GF(16) root.",
            "trace_1_plus_r": tr_f16_f2(1 ^ r),
            "conclusion": "No linear transition can meet any target eta*(1+r^2) in GF(16).",
        },
        "semilinear_transition_proof": {
            "maps": "y->alpha*y^2 (the three involutions)",
            "reductions": semilinear_reductions,
            "equations": {
                "membership_in_GF16": "g(u)=r^3*g(v), g(z)=z^2+z+1+r",
                "normalized_product": "(a^3+b^3)/r^2=(u+v)g(u)+(1+r^3)",
                "target": "eta*(1+r^3)",
            },
            "candidate_y_set": [
                {"hex": f"0x{x:02x}", "expression": expressions[x]}
                for x in y_values
            ],
            "four_row_table": table,
            "normalized_product_multiset": [
                {
                    "hex": f"0x{x:02x}",
                    "expression": expressions[x],
                    "multiplicity": normalized_product_counter[x],
                }
                for x in sorted(normalized_product_counter)
            ],
            "target_intersections": [
                {
                    "eta_hex": f"0x{eta:02x}",
                    "eta_expression": expressions[eta],
                    "normalized_target_hex": f"0x{targets[eta]:02x}",
                    "normalized_target_expression": expressions[targets[eta]],
                    "sections_per_involution": target_counts[eta],
                }
                for eta in (1, lam, lam2)
            ],
            "conclusion": "Only eta=lambda^2 occurs, exactly four times for each semilinear involution; alpha disappears because alpha^3=1.",
        },
        "proof_status": {
            "replacement": "The 10,795-subspace exhaustive certificate is replaced by a two-family proof and a four-row GF(16) table.",
            "finite_calculation_size": "16 reduced (u,v) pairs, organized by four y-values",
            "agreement_with_stage5g": True,
        },
        "provenance": {
            "stage5g_file": str(stage5g_path.relative_to(root)),
            "stage5g_sha256": hashlib.sha256(stage5g_path.read_bytes()).hexdigest(),
        },
        "validation": {
            "field_tower_identities": True,
            "six_GL2_maps_recovered": True,
            "linear_maps_excluded_symbolically": True,
            "semilinear_reduction_independent_of_alpha": True,
            "four_row_table_exact": True,
            "stage5g_transition_counts_reproduced": True,
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
        "linear_transition_hits": 0,
        "semilinear_hits_per_involution": target_counts[lam2],
        "total_excluded_witnesses": 3 * target_counts[lam2],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
