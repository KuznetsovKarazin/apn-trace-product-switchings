#!/usr/bin/env python3
r"""Stage 5P: field-independent Kummer certificate for the n=6 coefficient theorem.

Let E=GF(4)=GF(2)(lambda), lambda^2+lambda+1=0, and
K=E(u) with u^3=lambda.  Then K=GF(64), Tr_{K/E}(a+b*u+c*u^2)=a.
For an F2-linear section phi(y)=a*y+b*y^2 of the relative trace, write

a=1+p*u+q*u^2,  b=r*u+s*u^2,  p,q,r,s in E.

The bad-coefficient set is the set of section products a^3+b^3.  This script
performs the exact 4^4 calculation entirely over E, producing a portable
certificate that its complement is E^*u union E^*u^2, equivalently the six
roots of X^6+X^3+1, i.e. the elements of multiplicative order 9.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def e_add(a: int, b: int) -> int:
    return a ^ b


def e_mul(a: int, b: int) -> int:
    out = 0
    aa, bb = a, b
    while bb:
        if bb & 1:
            out ^= aa
        bb >>= 1
        aa <<= 1
        if aa & 0b100:
            aa ^= 0b111  # t^2+t+1
    return out & 0b11


def e_square(a: int) -> int:
    return e_mul(a, a)


def e_cube(a: int) -> int:
    return e_mul(e_square(a), a)


LAMBDA = 0b10
LAMBDA2 = 0b11


def k_add(x: tuple[int, int, int], y: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(a ^ b for a, b in zip(x, y))


def k_mul(x: tuple[int, int, int], y: tuple[int, int, int]) -> tuple[int, int, int]:
    # Multiply polynomials in u and reduce u^3=lambda.
    raw = [0] * 5
    for i, a in enumerate(x):
        for j, b in enumerate(y):
            raw[i + j] ^= e_mul(a, b)
    # u^4=lambda*u, u^3=lambda.
    raw[1] ^= e_mul(raw[4], LAMBDA)
    raw[0] ^= e_mul(raw[3], LAMBDA)
    return raw[0], raw[1], raw[2]


def k_pow(x: tuple[int, int, int], exponent: int) -> tuple[int, int, int]:
    out = (1, 0, 0)
    base = x
    e = exponent
    while e:
        if e & 1:
            out = k_mul(out, base)
        base = k_mul(base, base)
        e >>= 1
    return out


def cube_formula(x: tuple[int, int, int]) -> tuple[int, int, int]:
    a, b, c = x
    constant = e_cube(a) ^ e_mul(LAMBDA, e_cube(b)) ^ e_mul(LAMBDA2, e_cube(c))
    u_coeff = (
        e_mul(e_square(a), b)
        ^ e_mul(LAMBDA, e_mul(a, e_square(c)))
        ^ e_mul(LAMBDA, e_mul(e_square(b), c))
    )
    u2_coeff = (
        e_mul(a, e_square(b))
        ^ e_mul(e_square(a), c)
        ^ e_mul(LAMBDA, e_mul(b, e_square(c)))
    )
    result = constant, u_coeff, u2_coeff
    assert result == k_pow(x, 3)
    return result


def order(x: tuple[int, int, int]) -> int:
    if x == (0, 0, 0):
        return 0
    for d in (1, 3, 7, 9, 21, 63):
        if 63 % d == 0 and k_pow(x, d) == (1, 0, 0):
            return d
    raise AssertionError("order not found")


def triple_str(x: tuple[int, int, int]) -> str:
    names = {0: "0", 1: "1", 2: "lambda", 3: "lambda^2"}
    return f"({names[x[0]]},{names[x[1]]},{names[x[2]]})"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/kummer_n6_trace_switch_certificate.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.output if args.output.is_absolute() else root / args.output

    # Irreducibility of X^3-lambda over E follows because every nonzero
    # element of E has cube 1, so lambda is not a cube.
    assert all(e_cube(x) != LAMBDA for x in range(4))

    counts: Counter[tuple[int, int, int]] = Counter()
    records = []
    for p in range(4):
        for q in range(4):
            a = (1, p, q)
            a3 = cube_formula(a)
            for r in range(4):
                for s in range(4):
                    b = (0, r, s)
                    theta = k_add(a3, cube_formula(b))
                    counts[theta] += 1
                    records.append({
                        "p": p,
                        "q": q,
                        "r": r,
                        "s": s,
                        "theta": list(theta),
                    })
    assert sum(counts.values()) == 256
    assert counts[(0, 0, 0)] == 0

    all_nonzero = [
        (a, b, c)
        for a in range(4) for b in range(4) for c in range(4)
        if (a, b, c) != (0, 0, 0)
    ]
    missing = sorted(x for x in all_nonzero if counts[x] == 0)
    expected_missing = sorted(
        [(0, alpha, 0) for alpha in (1, LAMBDA, LAMBDA2)]
        + [(0, 0, alpha) for alpha in (1, LAMBDA, LAMBDA2)]
    )
    assert missing == expected_missing
    assert all(order(x) == 9 for x in missing)
    assert all(k_add(k_add(k_pow(x, 6), k_pow(x, 3)), (1, 0, 0)) == (0, 0, 0) for x in missing)

    represented = [x for x in all_nonzero if counts[x] > 0]
    assert len(represented) == 57
    multiplicity_distribution = Counter(counts[x] for x in represented)
    assert multiplicity_distribution == {3: 26, 4: 9, 6: 15, 7: 6, 10: 1}

    # Aggregate by Kummer coordinate support; this is a compact auditable
    # certificate of surjectivity outside the two missing E-lines.
    support_multiplicity = Counter()
    for x in represented:
        support = "".join(str(i) for i, value in enumerate(x) if value)
        support_multiplicity[(support, counts[x])] += 1

    result = {
        "schema": "kummer-n6-trace-switch-certificate-v1",
        "date": "2026-08-02",
        "field_model": {
            "E": "GF(4)=GF(2)[lambda]/(lambda^2+lambda+1)",
            "K": "E[u]/(u^3-lambda)=GF(64)",
            "relative_trace": "Tr_K/E(a+b*u+c*u^2)=a",
        },
        "section_parameterization": {
            "trace_one_element": "a=1+p*u+q*u^2",
            "trace_zero_element": "b=r*u+s*u^2",
            "parameters": "p,q,r,s in E",
            "cases": 256,
        },
        "cube_formula": {
            "constant": "a0^3+lambda*a1^3+lambda^2*a2^3",
            "u": "a0^2*a1+lambda*a0*a2^2+lambda*a1^2*a2",
            "u2": "a0*a1^2+a0^2*a2+lambda*a1*a2^2",
        },
        "theorem_certificate": {
            "represented_nonzero_count": len(represented),
            "missing_nonzero_count": len(missing),
            "missing_set": [triple_str(x) for x in missing],
            "missing_set_coordinate_free": "E^*u union E^*u^2",
            "equivalent_polynomial_condition": "theta^6+theta^3+1=0",
            "equivalent_group_condition": "ord(theta)=9",
            "multiplicity_distribution_on_represented_values": {
                str(k): v for k, v in sorted(multiplicity_distribution.items())
            },
            "support_and_multiplicity_summary": {
                f"support_{support}_count_{mult}": count
                for (support, mult), count in sorted(support_multiplicity.items())
            },
        },
        "validation": {
            "all_4_power_4_cases_enumerated_over_base_field_E": len(records) == 256,
            "cube_formula_checked_against_field_multiplication_every_case": True,
            "zero_never_represented": counts[(0, 0, 0)] == 0,
            "exactly_two_projective_E_lines_missing": True,
            "all_missing_elements_have_order_9": True,
        },
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "sha256": sha256(output),
        "represented": len(represented),
        "missing": len(missing),
        "multiplicity_distribution": dict(sorted(multiplicity_distribution.items())),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
