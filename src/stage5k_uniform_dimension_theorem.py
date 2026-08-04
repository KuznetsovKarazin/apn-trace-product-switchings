#!/usr/bin/env python3
r"""Stage 5K: uniform nonexistence theorem for the relative-trace rank-one switch.

Let K=GF(2^n), n even, GF(4) subset K, lambda in GF(4)\GF(2),

    Q(x)=Tr_K/F2(x) Tr_K/F2(lambda*x),
    F_theta(x)=x^3+theta Q(x).

Stages 5G/5I reduced APN to theta not belonging to

    S_n={a^3+b^3 : Tr_K/GF4(a)=1, Tr_K/GF4(b)=0}.

This stage proves S_n=K^* for every even n>=10.  The proof counts

    h_n(theta)=#{(a,b): T(a)=1,T(b)=0,a^3+b^3=theta}

on the Fermat cubic C_theta: a^3+b^3=theta.  Additive-character
orthogonality over GF(4) gives 16 curve sums.  The projective cubic is a
smooth genus-one curve with three rational points at infinity.  For a
nonzero linear phase u*a+v*b, the Artin--Schreier bound is 6*sqrt(q) on
an axis (three simple poles) and 4*sqrt(q) off the axes (two simple poles).
Hence

    h_n(theta) >= (q-2-74*sqrt(q))/16.

This is positive for q>=2^14.  Exact Stage-5I convolution certificates
cover n=10 and n=12, closing all even n>=10.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dimension-scan",
        type=Path,
        default=Path("results/current/trace_switch_dimension_scan.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/uniform_trace_switch_nonexistence.json"),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    resolve = lambda p: p if p.is_absolute() else root / p
    scan_path = resolve(args.dimension_scan)
    output_path = resolve(args.output)
    scan = json.loads(scan_path.read_text(encoding="utf-8"))
    rows = {row["n"]: row for row in scan["dimensions"]}

    # Exact finite certificates required below the uniform Weil threshold.
    finite_cases = []
    for n in (10, 12):
        row = rows[n]
        assert row["represented_nonzero_count"] == (1 << n) - 1
        assert row["safe_nonzero_count"] == 0
        assert row["minimum_positive_representation_count"] > 0
        finite_cases.append({
            "n": n,
            "q": 1 << n,
            "all_nonzero_theta_represented": True,
            "minimum_h_n_theta": row["minimum_positive_representation_count"],
            "certificate": "exact XOR convolution from Stage 5I",
        })

    # Curve-bound constants.
    axis_pairs = 6       # (u,0),(0,v), u/v nonzero in GF(4)
    off_axis_pairs = 9   # u,v both nonzero
    axis_bound = 6       # 3 simple poles, genus 1: 2*3 sqrt(q)
    off_axis_bound = 4   # 2 simple poles, genus 1: 2*2 sqrt(q)
    hasse_bound = 2      # projective genus-one point-count error
    total_sqrt_constant = (
        hasse_bound + axis_pairs * axis_bound + off_axis_pairs * off_axis_bound
    )
    assert total_sqrt_constant == 74

    threshold_rows = []
    for n in range(10, 31, 2):
        q = 1 << n
        sqrt_q = 1 << (n // 2)
        numerator = q - 2 - total_sqrt_constant * sqrt_q
        lower_bound = numerator / 16
        threshold_rows.append({
            "n": n,
            "q": q,
            "sqrt_q": sqrt_q,
            "lower_bound_numerator": numerator,
            "lower_bound_h": lower_bound,
            "strictly_positive": numerator > 0,
        })
    assert not next(row for row in threshold_rows if row["n"] == 12)["strictly_positive"]
    assert next(row for row in threshold_rows if row["n"] == 14)["strictly_positive"]
    assert all(row["strictly_positive"] for row in threshold_rows if row["n"] >= 14)

    result: dict[str, Any] = {
        "schema": "uniform-relative-trace-switch-nonexistence-v1",
        "date": "2026-08-02",
        "family": {
            "field": "K=GF(2^n), n even, with GF(4) embedded",
            "kernel": "Q(x)=Tr_K/F2(x)*Tr_K/F2(lambda*x), lambda in GF(4)\\GF(2)",
            "function": "F_theta(x)=x^3+theta*Q(x)",
        },
        "previous_exact_reduction": {
            "criterion": "F_theta is APN iff theta is not in S_n",
            "S_n": "{a^3+b^3 : Tr_K/GF4(a)=1, Tr_K/GF4(b)=0}",
            "multiplicity": "h_n(theta)=#{(a,b) in the two trace fibres with a^3+b^3=theta}",
        },
        "character_expansion": {
            "formula": "h_n(theta)=1/16 * sum_{u,v in GF(4)} (-1)^Tr_GF4/F2(u) * sum_{a^3+b^3=theta} (-1)^Tr_K/F2(u*a+v*b)",
            "base_curve": "C_theta: X^3+Y^3=theta*Z^3 is smooth of genus 1 for theta!=0",
            "points_at_infinity": 3,
            "base_term": "#C_theta(affine)=q-2+e_theta, |e_theta|<=2*sqrt(q)",
        },
        "nontrivial_sum_bounds": {
            "justification": "The phase u*x+v*y has simple poles at infinity and cannot be an Artin-Schreier coboundary; the standard Artin-Schreier/Weil bound is (2g-2+sum(m_P+1))*sqrt(q).",
            "axis_pairs": axis_pairs,
            "axis_simple_poles": 3,
            "axis_bound": "6*sqrt(q)",
            "off_axis_pairs": off_axis_pairs,
            "off_axis_simple_poles": 2,
            "off_axis_bound": "4*sqrt(q)",
            "total_error_bound": "74*sqrt(q)",
        },
        "uniform_lower_bound": {
            "formula": "h_n(theta) >= (q-2-74*sqrt(q))/16 for every theta!=0",
            "first_even_dimension_proved_positive_by_bound": 14,
            "threshold_table": threshold_rows,
        },
        "finite_bridge": finite_cases,
        "theorem": {
            "statement": "For every even n>=10 and every theta!=0, h_n(theta)>0; equivalently S_n=GF(2^n)^*.",
            "APN_consequence": "For every even n>=10, x^3+theta*Q(x) is not APN for any theta!=0.",
            "exceptional_dimension": "n=8 is the last nontrivial dimension of this fixed relative-trace rank-one construction; there the APN coefficients are exactly GF(4)^*.",
            "status": "uniform theorem, using the curve bound for n>=14 and exact finite certificates for n=10,12",
        },
        "provenance": {
            "dimension_scan_file": str(scan_path.relative_to(root)),
            "dimension_scan_sha256": sha256(scan_path),
        },
        "validation": {
            "error_constant_is_74": total_sqrt_constant == 74,
            "bound_positive_for_all_even_n_ge_14_in_table": all(
                row["strictly_positive"] for row in threshold_rows if row["n"] >= 14
            ),
            "n10_n12_exactly_closed": all(
                row["all_nonzero_theta_represented"] for row in finite_cases
            ),
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "sha256": sha256(output_path),
        "uniform_theorem": True,
        "finite_cases": [10, 12],
        "curve_bound_from_n": 14,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
