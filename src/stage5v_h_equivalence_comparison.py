#!/usr/bin/env python3
"""Stage 5V: exact comparison with the H-equivalence criterion.

The script exhaustively enumerates all GF(2)-linear maps L on GF(16) satisfying
L(e0)=0, and verifies equivalence of:
  (i) direct APN-ness of G(x)=x^3+Tr(x)L(x),
 (ii) the all-direction polar-kernel criterion,
(iii) Taniguchi et al.'s hyperplane injectivity criterion.

It also reproduces the published count 448 for n=4.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

MOD = 0x13  # x^4+x+1
N = 16


def mul(a: int, b: int) -> int:
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 0x10:
            a ^= MOD
    return out & 0xF


def power(a: int, e: int) -> int:
    out = 1
    while e:
        if e & 1:
            out = mul(out, a)
        a = mul(a, a)
        e >>= 1
    return out


def trace(a: int) -> int:
    out = 0
    x = a
    for _ in range(4):
        out ^= x
        x = mul(x, x)
    assert out in (0, 1)
    return out


def linear(columns: tuple[int, int, int, int], x: int) -> int:
    out = 0
    for i in range(4):
        if (x >> i) & 1:
            out ^= columns[i]
    return out


def base(x: int) -> int:
    return power(x, 3)


def polar_f(a: int, x: int) -> int:
    return base(x ^ a) ^ base(x) ^ base(a)


def switched(columns, x: int) -> int:
    return base(x) ^ (linear(columns, x) if trace(x) else 0)


def direct_apn(columns) -> bool:
    g = [switched(columns, x) for x in range(N)]
    for a in range(1, N):
        counts = [0] * N
        for x in range(N):
            counts[g[x ^ a] ^ g[x]] += 1
        if max(counts) > 2:
            return False
    return True


def all_direction_kernel(columns) -> bool:
    for a in range(1, N):
        kernel = []
        for x in range(N):
            delta = (
                polar_f(a, x)
                ^ (linear(columns, x) if trace(a) else 0)
                ^ (linear(columns, a) if trace(x) else 0)
            )
            if delta == 0:
                kernel.append(x)
        if set(kernel) != {0, a}:
            return False
    return True


def h_injectivity(columns, e0: int, t0: list[int]) -> bool:
    for a in t0:
        direction = a ^ e0  # every trace-one direction exactly once
        values = [linear(columns, x) ^ polar_f(x, direction) for x in t0]
        if len(set(values)) != len(t0):
            return False
    return True


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    t0 = [x for x in range(N) if trace(x) == 0]
    t1 = [x for x in range(N) if trace(x) == 1]
    e0 = min(t1)

    checked = 0
    accepted = 0
    disagreements = []
    examples = []

    # Enumerate all 4x4 binary matrices as four output columns.
    for encoded in range(1 << 16):
        columns = tuple((encoded >> (4 * i)) & 0xF for i in range(4))
        if linear(columns, e0) != 0:
            continue
        checked += 1
        h = h_injectivity(columns, e0, t0)
        k = all_direction_kernel(columns)
        d = direct_apn(columns)
        if not (h == k == d):
            disagreements.append({
                "encoded": encoded,
                "columns": list(columns),
                "h_injectivity": h,
                "all_direction_kernel": k,
                "direct_apn": d,
            })
            if len(disagreements) >= 10:
                break
        if d:
            accepted += 1
            if len(examples) < 8:
                examples.append({"encoded": encoded, "columns": list(columns)})

    assert checked == 4096
    assert not disagreements
    assert accepted == 448

    result = {
        "schema": "h-equivalence-comparison-v1",
        "date": "2026-08-02",
        "field": "GF(16)=GF(2)[z]/(z^4+z+1)",
        "centre": "F(x)=x^3",
        "switch": "G(x)=F(x)+Tr(x)L(x)",
        "normalization": {
            "e0": e0,
            "e0_hex": f"0x{e0:x}",
            "trace_zero_size": len(t0),
            "linear_maps_with_L_e0_zero": checked,
        },
        "criteria": {
            "direct_apn": "DDT multiplicity at most 2",
            "all_direction_kernel": (
                "ker(B_F(a,.)+Tr(a)L(.)+Tr(.)L(a))={0,a} for every a!=0"
            ),
            "h_injectivity": (
                "x in T0 -> L(x)+B_F(x,a+e0) injective for every a in T0"
            ),
        },
        "result": {
            "all_three_criteria_agree": True,
            "disagreement_count": 0,
            "accepted_map_count": accepted,
            "published_count_reproduced": 448,
            "sample_accepted_maps": examples,
        },
        "formal_relation": {
            "trace_one_direction": (
                "For d with Tr(d)=1, every coset modulo <d> has a unique "
                "trace-zero representative x; the perturbed polar kernel "
                "condition becomes B_F(d,x)+L(x)=0."
            ),
            "trace_zero_direction": (
                "For d with Tr(d)=0, an extra kernel vector must have trace "
                "one and satisfies B_F(d,x)+L(d)=0. By symmetry B_F(d,x)="
                "B_F(x,d), this is the H-injectivity test for the trace-one "
                "direction x evaluated at d in T0."
            ),
            "conclusion": (
                "The H-equivalence criterion is the common-trace-factor "
                "specialization of the symmetric all-direction polar-kernel "
                "criterion. The general low-rank update lemma extends the "
                "setup to arbitrary quadratic Boolean kernels and coefficient "
                "rank r, while no priority is claimed for the hyperplane case."
            ),
        },
    }

    root = Path(__file__).resolve().parents[1]
    output = root / "results/current/h_equivalence_criterion_comparison.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "sha256": sha256(output),
        "checked_maps": checked,
        "accepted_maps": accepted,
        "criteria_agree": True,
    }, indent=2))


if __name__ == "__main__":
    main()
