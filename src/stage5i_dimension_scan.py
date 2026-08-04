#!/usr/bin/env python3
r"""Stage 5I: exact dimension scan for rank-one relative-trace switching.

For even n, K=GF(2^n), GF(4) subset K, lambda in GF(4)\GF(2), and

    Q(x)=Tr_K/F2(x) Tr_K/F2(lambda*x),
    F_theta(x)=x^3+theta*Q(x),

a bad derivative is equivalent to a transverse two-dimensional F2-subspace.
Writing the inverse relative-trace section as phi(y)=a*y+b*y^2 gives the
exact criterion

    F_theta APN <=> theta not in
        S_n={a^3+b^3 : Tr_K/GF4(a)=1, Tr_K/GF4(b)=0}.

The representation counts of S_n are computed exactly as an XOR convolution
using the Walsh-Hadamard transform.  No random APN testing is used for the
reported dimension table.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


def poly_degree(p: int) -> int:
    return p.bit_length() - 1


def poly_mul(a: int, b: int) -> int:
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
    return out


def poly_mod(a: int, modulus: int) -> int:
    dm = poly_degree(modulus)
    while poly_degree(a) >= dm:
        a ^= modulus << (poly_degree(a) - dm)
    return a


def poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, poly_mod(a, b)
    return a


def poly_square_mod(a: int, modulus: int) -> int:
    return poly_mod(poly_mul(a, a), modulus)


def is_irreducible_binary(polynomial: int, degree: int) -> bool:
    x = 2
    power = x
    for _ in range(1, degree // 2 + 1):
        power = poly_square_mod(power, polynomial)
        if poly_gcd(power ^ x, polynomial) != 1:
            return False
    power = x
    for _ in range(degree):
        power = poly_square_mod(power, polynomial)
    return power == x


def first_irreducible_polynomial(degree: int) -> int:
    for low in range(1, 1 << degree, 2):
        candidate = (1 << degree) | low
        if is_irreducible_binary(candidate, degree):
            return candidate
    raise RuntimeError(f"no irreducible polynomial found for degree {degree}")


class BinaryField:
    def __init__(self, degree: int, modulus: int):
        self.degree = degree
        self.modulus = modulus
        self.size = 1 << degree
        self.mask = self.size - 1

    def mul(self, a: int, b: int) -> int:
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a >> self.degree:
                a ^= self.modulus
        return out & self.mask

    def pow(self, a: int, exponent: int) -> int:
        out = 1
        base = a
        e = exponent
        while e:
            if e & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            e >>= 1
        return out


def fwht(values: np.ndarray) -> np.ndarray:
    result = values.astype(np.int64, copy=True)
    width = 1
    while width < result.size:
        blocks = result.reshape(-1, 2 * width)
        left = blocks[:, :width].copy()
        right = blocks[:, width:].copy()
        blocks[:, :width] = left + right
        blocks[:, width:] = left - right
        width *= 2
    return result


def field_tables(field: BinaryField) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = field.degree
    N = field.size
    square = np.empty(N, dtype=np.uint32)
    for x in range(N):
        square[x] = field.mul(x, x)
    fourth = square[square]
    relative_trace = np.arange(N, dtype=np.uint32)
    term = np.arange(N, dtype=np.uint32)
    for _ in range(1, n // 2):
        term = fourth[term]
        relative_trace ^= term
    cube = np.empty(N, dtype=np.uint32)
    for x in range(N):
        cube[x] = field.mul(x, int(square[x]))
    return square, relative_trace, cube


def trace_f4_to_f2(field: BinaryField, value: int) -> int:
    result = value ^ field.pow(value, 2)
    if result not in (0, 1):
        raise AssertionError("value is not in GF(4)")
    return result


def direct_apn_check(
    field: BinaryField,
    cube: np.ndarray,
    relative_trace: np.ndarray,
    lambda_value: int,
    theta: int,
) -> bool:
    N = field.size
    q = np.zeros(N, dtype=np.uint8)
    for x in range(N):
        tx = int(relative_trace[x])
        q[x] = trace_f4_to_f2(field, tx) & trace_f4_to_f2(
            field, field.mul(lambda_value, tx)
        )
    lut = cube.copy()
    lut[q == 1] ^= theta
    # F is quadratic and F(0)=0; APN iff every nonzero polar derivative has
    # kernel exactly {0,a}.
    for a in range(1, N):
        fa = int(lut[a])
        zeros = 0
        for x in range(N):
            if (int(lut[x ^ a]) ^ int(lut[x]) ^ fa) == 0:
                zeros += 1
                if zeros > 2:
                    return False
        if zeros != 2:
            return False
    return True


def scan_dimension(n: int, direct_validate: bool) -> dict[str, Any]:
    modulus = first_irreducible_polynomial(n)
    field = BinaryField(n, modulus)
    N = field.size
    square, relative_trace, cube = field_tables(field)

    trace_values = sorted(int(x) for x in np.unique(relative_trace))
    # The image of the relative trace is exactly the embedded GF(4).
    # Reusing the four observed trace values avoids scanning the whole field
    # with a separate x^4=x test in the large dimensions.
    f4 = trace_values
    f4_star = [x for x in f4 if x]
    lambda_value = next(x for x in f4 if x not in (0, 1))
    assert len(f4) == 4

    trace_one_mask = relative_trace == 1
    trace_zero_mask = relative_trace == 0
    fibre_size = 1 << (n - 2)
    assert int(trace_one_mask.sum()) == fibre_size
    assert int(trace_zero_mask.sum()) == fibre_size

    first_frequency = np.bincount(
        cube[trace_one_mask], minlength=N
    ).astype(np.int64)
    second_frequency = np.bincount(
        cube[trace_zero_mask], minlength=N
    ).astype(np.int64)
    counts = fwht(fwht(first_frequency) * fwht(second_frequency)) // N
    assert np.all(counts >= 0)
    assert int(counts.sum()) == fibre_size * fibre_size
    # Zero cannot be represented: a^3=b^3 would give b=omega*a with
    # omega in GF(4)^*, contradicting T(a)=1 and T(b)=0.
    assert int(counts[0]) == 0

    safe = (np.flatnonzero(counts[1:] == 0) + 1).astype(int).tolist()
    represented_nonzero = int(np.count_nonzero(counts[1:]))
    positive = counts[1:][counts[1:] > 0]
    multiplicity_distribution = (
        Counter(int(x) for x in positive.tolist()) if n <= 14 else Counter()
    )

    direct_record: dict[str, Any] | None = None
    if direct_validate:
        direct_safe = []
        for theta in range(1, N):
            if direct_apn_check(field, cube, relative_trace, lambda_value, theta):
                direct_safe.append(theta)
        assert direct_safe == safe
        direct_record = {
            "all_nonzero_coefficients_checked": N - 1,
            "safe_coefficients_match_convolution": True,
        }

    return {
        "n": n,
        "extension_degree_over_GF4": n // 2,
        "field_modulus_hex": f"0x{modulus:x}",
        "field_size": N,
        "GF4_hex": [f"0x{x:x}" for x in f4],
        "lambda_hex": f"0x{lambda_value:x}",
        "relative_trace_fibre_size": fibre_size,
        "section_pair_count": fibre_size * fibre_size,
        "represented_nonzero_count": represented_nonzero,
        "safe_nonzero_count": len(safe),
        "safe_equals_GF4_star": set(safe) == set(f4_star),
        "safe_coefficients_hex": [f"0x{x:x}" for x in safe] if n <= 10 else None,
        "minimum_positive_representation_count": int(positive.min()) if positive.size else None,
        "maximum_representation_count": int(positive.max()) if positive.size else None,
        "expected_average_nonzero_count": (fibre_size * fibre_size) / (N - 1),
        "multiplicity_distribution": {
            str(k): v for k, v in sorted(multiplicity_distribution.items())
        } if n <= 14 else None,
        "direct_APN_validation": direct_record,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dimensions",
        default="2,4,6,8,10,12,14,16,18,20,22",
        help="comma-separated even dimensions",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/trace_switch_dimension_scan.json"),
    )
    args = parser.parse_args()

    dimensions = [int(x.strip()) for x in args.dimensions.split(",") if x.strip()]
    if any(n <= 0 or n % 2 for n in dimensions):
        raise ValueError("all dimensions must be positive and even")

    root = Path(__file__).resolve().parents[2]
    output_path = args.output if args.output.is_absolute() else root / args.output

    rows = []
    for n in dimensions:
        rows.append(scan_dimension(n, direct_validate=n <= 8))

    row_by_n = {row["n"]: row for row in rows}
    if 8 in row_by_n:
        assert row_by_n[8]["safe_equals_GF4_star"]
    tested_large = [row for row in rows if row["n"] >= 10]
    assert all(row["safe_nonzero_count"] == 0 for row in tested_large)

    result = {
        "schema": "trace-switch-dimension-scan-v1",
        "date": "2026-08-02",
        "general_exact_criterion": {
            "setting": "K=GF(2^n), n even, T=Tr_K/GF4, Q(x)=Tr_K/F2(x)Tr_K/F2(lambda*x), lambda in GF(4)\\GF(2).",
            "switch": "F_theta(x)=x^3+theta*Q(x)",
            "bad_set": "S_n={a^3+b^3 : T(a)=1, T(b)=0}",
            "theorem": "For theta!=0, F_theta is APN iff theta is not in S_n.",
            "derivation": "A nontrivial derivative-kernel element defines a 2D F2-subspace V transverse to ker(T), with product theta. The inverse section phi(y)=a*y+b*y^2 satisfies T(a)=1,T(b)=0 and prod(V)=a^3+b^3.",
            "convolution": "If f(u)=#{a:T(a)=1,a^3=u} and g(v)=#{b:T(b)=0,b^3=v}, then the witness count is the XOR convolution h(theta)=sum_u f(u)g(theta+u).",
        },
        "dimensions": rows,
        "observations": {
            "dimension_8": "The safe set is exactly GF(4)^*, agreeing with the symbolic Stage 5G theorem.",
            "dimensions_10_through_22": "Every nonzero coefficient has at least one bad section; no nontrivial switch in this fixed relative-trace family is APN.",
            "small_dimensions": "Dimensions 2,4,6 have exceptional safe sets different from the dimension-8 GF(4)^* fibre.",
            "conjecture": "For every even n>=10, S_n=GF(2^n)^*. This is exact for n=10,12,...,22 but is not yet proved uniformly.",
        },
        "method": {
            "type": "exact integer XOR convolution via Walsh-Hadamard transform",
            "random_testing": False,
            "direct_validation": "All nonzero coefficients were independently APN-tested for n=2,4,6,8.",
        },
        "validation": {
            "all_convolution_counts_nonnegative": True,
            "all_pair_counts_conserved": True,
            "zero_never_represented": True,
            "direct_checks_n_le_8_match": True,
            "stage5g_n8_reproduced": True,
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
        "summary": {str(row["n"]): row["safe_nonzero_count"] for row in rows},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
