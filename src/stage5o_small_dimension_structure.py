#!/usr/bin/env python3
r"""Stage 5O: exact structure of the relative-trace switch in n=4 and n=6.

For K=GF(2^n), E=GF(4), lambda in E\GF(2), and

    Q(x)=Tr_K/F2(x) Tr_K/F2(lambda*x),
    F_theta(x)=x^3+theta Q(x),

this script gives exact algebraic descriptions of the APN coefficient set for
n=4 and n=6, validates APN directly, identifies EA/Frobenius structure, computes
orthoderivative signatures, and applies the n/2-dimensional non-bent-component
subspace screen in n=6.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path
from typing import Iterable, Any

import numpy as np

from stage5i_dimension_scan import (
    BinaryField,
    first_irreducible_polynomial,
    field_tables,
    trace_f4_to_f2,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parity(x: int) -> int:
    return x.bit_count() & 1


def absolute_trace(field: BinaryField, x: int) -> int:
    s = 0
    y = x
    for _ in range(field.degree):
        s ^= y
        y = field.mul(y, y)
    if s not in (0, 1):
        raise AssertionError("absolute trace did not land in GF(2)")
    return s


def multiplicative_order(field: BinaryField, x: int) -> int:
    if x == 0:
        return 0
    order = field.size - 1
    d = 1
    while d <= order:
        if order % d == 0 and field.pow(x, d) == 1:
            return d
        d += 1
    raise AssertionError("order not found")


def build_q_and_lut(field: BinaryField, theta: int) -> tuple[list[int], list[int]]:
    _, relative_trace, cube = field_tables(field)
    f4 = sorted(set(int(x) for x in relative_trace))
    lam = next(x for x in f4 if x not in (0, 1))
    q = []
    lut = []
    for x in range(field.size):
        tx = int(relative_trace[x])
        qx = trace_f4_to_f2(field, tx) & trace_f4_to_f2(field, field.mul(lam, tx))
        q.append(qx)
        lut.append(int(cube[x]) ^ (theta if qx else 0))
    return q, lut


def is_apn_quadratic(lut: list[int]) -> bool:
    N = len(lut)
    for a in range(1, N):
        fa = lut[a]
        zeros = 0
        for x in range(N):
            if (lut[x ^ a] ^ lut[x] ^ fa) == 0:
                zeros += 1
                if zeros > 2:
                    return False
        if zeros != 2:
            return False
    return True


def image_profile(lut: list[int]) -> dict[str, int]:
    counts = Counter(lut)
    profile = Counter(counts.values())
    profile[0] = len(lut) - len(counts)
    return {str(k): v for k, v in sorted(profile.items())}


def binary_rank(rows: Iterable[int], width: int) -> int:
    basis = [0] * width
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


def span_nonzero(basis: list[int]) -> tuple[int, ...]:
    vals = []
    for mask in range(1, 1 << len(basis)):
        x = 0
        for i, row in enumerate(basis):
            if (mask >> i) & 1:
                x ^= row
        vals.append(x)
    return tuple(sorted(vals))


def rref_subspaces(k: int, n: int):
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
            total += 1
            yield tuple(rows), span_nonzero(rows)
    # Gaussian binomial values used in the two cases here.
    expected = {(2, 4): 35, (3, 6): 1395}.get((k, n))
    if expected is not None:
        assert total == expected


def fwht(values: np.ndarray) -> np.ndarray:
    out = values.astype(np.int64, copy=True)
    step = 1
    while step < out.size:
        blocks = out.reshape(-1, 2 * step)
        left = blocks[:, :step].copy()
        right = blocks[:, step:].copy()
        blocks[:, :step] = left + right
        blocks[:, step:] = left - right
        step *= 2
    return out


def component_data(lut: list[int]) -> tuple[Counter[int], set[int], Counter[int]]:
    N = len(lut)
    amplitude = Counter()
    non_bent: set[int] = set()
    full_walsh = Counter()
    for b in range(1, N):
        signs = np.fromiter(
            (1 if parity(b & lut[x]) == 0 else -1 for x in range(N)),
            dtype=np.int64,
            count=N,
        )
        spectrum = fwht(signs)
        full_walsh.update(int(v) for v in spectrum)
        nonzero = spectrum[spectrum != 0]
        amps = set(abs(int(v)) for v in nonzero)
        if len(amps) != 1:
            raise AssertionError("quadratic component has multiple nonzero amplitudes")
        amp = next(iter(amps))
        amplitude[amp] += 1
        if np.count_nonzero(nonzero) != N:
            non_bent.add(b)
    return amplitude, non_bent, full_walsh


def orthoderivative(lut: list[int]) -> list[int]:
    N = len(lut)
    result = [0] * N
    for a in range(1, N):
        image = {lut[x ^ a] ^ lut[x] ^ lut[a] for x in range(N)}
        if len(image) != N // 2:
            raise AssertionError("non-APN derivative image")
        candidates = [
            w for w in range(1, N)
            if all(parity(w & value) == 0 for value in image)
        ]
        if len(candidates) != 1:
            raise AssertionError("orthoderivative left kernel not one-dimensional")
        result[a] = candidates[0]
    return result


def orthoderivative_signature(lut: list[int]) -> dict[str, Any]:
    N = len(lut)
    ortho = orthoderivative(lut)
    differential = Counter()
    for a in range(1, N):
        row = [0] * N
        for x in range(N):
            row[ortho[x ^ a] ^ ortho[x]] += 1
        differential.update(row)
    walsh = Counter()
    for b in range(1, N):
        signs = np.fromiter(
            (1 if parity(b & ortho[x]) == 0 else -1 for x in range(N)),
            dtype=np.int64,
            count=N,
        )
        walsh.update(abs(int(v)) for v in fwht(signs))
    return {
        "differential_spectrum": {str(k): v for k, v in sorted(differential.items())},
        "absolute_walsh_spectrum": {str(k): v for k, v in sorted(walsh.items())},
    }



def vector_anf_coefficients(lut: list[int], n: int) -> list[int]:
    coeffs = lut.copy()
    for i in range(n):
        for mask in range(1 << n):
            if (mask >> i) & 1:
                coeffs[mask] ^= coeffs[mask ^ (1 << i)]
    return coeffs


def polynomial_from_roots(field: BinaryField, roots: list[int]) -> list[int]:
    coeffs = [1]
    for root in roots:
        new = [0] * (len(coeffs) + 1)
        for i, c in enumerate(coeffs):
            new[i] ^= field.mul(c, root)
            new[i + 1] ^= c
        coeffs = new
    return coeffs


def poly_terms(coeffs: list[int]) -> str:
    terms = []
    for degree, coeff in enumerate(coeffs):
        if coeff == 0:
            continue
        if coeff != 1:
            terms.append(f"0x{coeff:x}*X^{degree}")
        elif degree == 0:
            terms.append("1")
        elif degree == 1:
            terms.append("X")
        else:
            terms.append(f"X^{degree}")
    return "+".join(reversed(terms))


def generate_gl4() -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    cols: list[int] = []
    def rec() -> None:
        if len(cols) == 4:
            out.append(tuple(cols))
            return
        for v in range(1, 16):
            if binary_rank(cols + [v], 4) == len(cols) + 1:
                cols.append(v)
                rec()
                cols.pop()
    rec()
    assert len(out) == 20160
    return out


def apply_linear(cols: tuple[int, ...] | list[int], x: int) -> int:
    y = 0
    for i, col in enumerate(cols):
        if (x >> i) & 1:
            y ^= col
    return y


def polar(lut: list[int], x: int, y: int) -> int:
    return lut[x] ^ lut[y] ^ lut[x ^ y] ^ lut[0]


def solve_output_map(equations: list[tuple[int, int]], width: int) -> tuple[bool, list[int] | None]:
    bu = [0] * width
    bv = [0] * width
    for u, v in equations:
        uu, vv = u, v
        for p in range(width - 1, -1, -1):
            if ((uu >> p) & 1) and bu[p]:
                uu ^= bu[p]
                vv ^= bv[p]
        if uu == 0:
            if vv != 0:
                return False, None
        else:
            p = uu.bit_length() - 1
            bu[p] = uu
            bv[p] = vv
    if sum(x != 0 for x in bu) != width or binary_rank([x for x in bv if x], width) != width:
        return False, None
    columns = []
    for i in range(width):
        uu = 1 << i
        vv = 0
        for p in range(width - 1, -1, -1):
            if (uu >> p) & 1:
                if not bu[p]:
                    return False, None
                uu ^= bu[p]
                vv ^= bv[p]
        if uu:
            return False, None
        columns.append(vv)
    for u, v in equations:
        if apply_linear(columns, u) != v:
            return False, None
    return True, columns


def ea_polar_witness(source: list[int], target: list[int], gl4: list[tuple[int, ...]]) -> dict[str, Any] | None:
    pairs = [(1 << i, 1 << j) for i in range(4) for j in range(i + 1, 4)]
    target_values = [polar(target, x, y) for x, y in pairs]
    for A in gl4:
        equations = []
        for (x, y), v in zip(pairs, target_values):
            equations.append((polar(source, apply_linear(A, x), apply_linear(A, y)), v))
        ok, B = solve_output_map(equations, 4)
        if not ok or B is None:
            continue
        difference = [target[x] ^ apply_linear(B, source[apply_linear(A, x)]) for x in range(16)]
        linear_cols = [difference[1 << i] for i in range(4)]
        if all(difference[x] == apply_linear(linear_cols, x) for x in range(16)):
            return {
                "input_linear_columns_hex": [f"0x{x:x}" for x in A],
                "output_linear_columns_hex": [f"0x{x:x}" for x in B],
                "added_linear_columns_hex": [f"0x{x:x}" for x in linear_cols],
            }
    return None


def n4_record() -> dict[str, Any]:
    field = BinaryField(4, first_irreducible_polynomial(4))
    safe = [x for x in range(1, 16) if absolute_trace(field, x) == 0]
    assert safe == list(range(1, 8))
    coeffs = polynomial_from_roots(field, safe)
    assert all(c in (0, 1) for c in coeffs)
    assert coeffs == [1, 1, 0, 1, 0, 0, 0, 1]

    q, cube_lut = build_q_and_lut(field, 0)
    # In this representation Tr(y) is the top output coordinate.  The
    # difference Q(x)+Tr(x^3) is linear.
    trace_of_cube = [absolute_trace(field, cube_lut[x]) for x in range(16)]
    linear_difference = [q[x] ^ trace_of_cube[x] for x in range(16)]
    linear_mask = next(
        m for m in range(16)
        if all(linear_difference[x] == parity(m & x) for x in range(16))
    )

    direct = []
    for theta in safe:
        _, lut = build_q_and_lut(field, theta)
        direct.append({
            "theta_hex": f"0x{theta:x}",
            "APN": is_apn_quadratic(lut),
            "direct_permutation": len(set(lut)) == 16,
            "image_size": len(set(lut)),
            "preimage_profile": image_profile(lut),
        })
    assert all(row["APN"] for row in direct)

    gl4 = generate_gl4()
    _, base = build_q_and_lut(field, 0)
    witnesses = {}
    for theta in safe:
        _, lut = build_q_and_lut(field, theta)
        witness = ea_polar_witness(base, lut, gl4)
        if witness is None:
            raise AssertionError(f"no EA witness for theta={theta}")
        witnesses[f"0x{theta:x}"] = witness

    # Closed-form output map B_theta(y)=y+theta*Tr(y).
    closed_form_checks = []
    for theta in range(1, 16):
        cols = [(1 << i) ^ (theta if absolute_trace(field, 1 << i) else 0) for i in range(4)]
        invertible = binary_rank(cols, 4) == 4
        closed_form_checks.append({
            "theta_hex": f"0x{theta:x}",
            "trace": absolute_trace(field, theta),
            "B_theta_invertible": invertible,
        })
        assert invertible == (absolute_trace(field, theta) == 0)

    return {
        "n": 4,
        "field_modulus_hex": "0x13",
        "safe_coefficients_hex": [f"0x{x:x}" for x in safe],
        "exact_safe_condition": "theta != 0 and Tr_GF16/GF2(theta)=0",
        "safe_root_polynomial": poly_terms(coeffs),
        "safe_root_polynomial_coefficients_low_to_high": coeffs,
        "normal_form": {
            "identity_modulo_linear": "Q(x) = Tr_GF16/GF2(x^3) + ell(x)",
            "linear_difference_mask_hex": f"0x{linear_mask:x}",
            "EA_formula": "F_theta(x) = B_theta(x^3) + theta*ell(x), B_theta(y)=y+theta*Tr(y)",
            "B_theta_invertible_iff": "Tr(theta)=0",
        },
        "all_safe_switches_EA_equivalent_to_Gold_cube": True,
        "EA_witnesses_to_x3": witnesses,
        "direct_records": direct,
        "closed_form_output_map_checks": closed_form_checks,
        "interpretation": "This family occupies only the Gold EA class; it does not realize the second 4-bit APN EA class.",
    }


def n6_record() -> dict[str, Any]:
    field = BinaryField(6, first_irreducible_polynomial(6))
    _, relative_trace, _ = field_tables(field)
    f4 = sorted(set(int(x) for x in relative_trace))
    f4_star = [x for x in f4 if x]

    # Exact coefficient characterization.
    safe = [
        x for x in range(1, 64)
        if field.pow(x, 3) in [e for e in f4 if e not in (0, 1)]
    ]
    expected = [0x06, 0x0B, 0x14, 0x1A, 0x1C, 0x1F]
    assert safe == expected
    coeffs = polynomial_from_roots(field, safe)
    assert coeffs == [1, 0, 0, 1, 0, 0, 1]
    orders = {x: multiplicative_order(field, x) for x in safe}
    assert set(orders.values()) == {9}

    # One Frobenius orbit of length six.
    orbit = []
    x = safe[0]
    while x not in orbit:
        orbit.append(x)
        x = field.mul(x, x)
    assert set(orbit) == set(safe) and len(orbit) == 6

    # Two GF(4)^* cosets.
    remaining = set(safe)
    cosets = []
    while remaining:
        rep = min(remaining)
        coset = sorted(field.mul(rep, e) for e in f4_star)
        cosets.append(coset)
        remaining -= set(coset)
    assert len(cosets) == 2

    direct = []
    signatures = []
    nb_sets = []
    for theta in safe:
        _, lut = build_q_and_lut(field, theta)
        assert is_apn_quadratic(lut)
        amplitude, nb, full_walsh = component_data(lut)
        signature = orthoderivative_signature(lut)
        signatures.append(signature)
        nb_sets.append(nb)
        direct.append({
            "theta_hex": f"0x{theta:02x}",
            "APN": True,
            "direct_permutation": len(set(lut)) == 64,
            "image_size": len(set(lut)),
            "preimage_profile": image_profile(lut),
            "component_amplitude_distribution": {str(k): v for k, v in sorted(amplitude.items())},
            "non_bent_component_count": len(nb),
            "extended_component_walsh_distribution": {str(k): v for k, v in sorted(full_walsh.items())},
            "orthoderivative_signature": signature,
        })
    assert all(sig == signatures[0] for sig in signatures)

    # Enumerate all 3-spaces and test each output-coordinate NB set.
    subspaces = list(rref_subspaces(3, 6))
    assert len(subspaces) == 1395
    permutation_screens = []
    for nb in nb_sets:
        contained = []
        for basis, vectors in subspaces:
            if set(vectors) <= nb:
                contained.append({
                    "basis_hex": [f"0x{x:02x}" for x in basis],
                    "vectors_hex": [f"0x{x:02x}" for x in vectors],
                })
        complementary = []
        for i in range(len(contained)):
            bi = [int(x, 16) for x in contained[i]["basis_hex"]]
            for j in range(i + 1, len(contained)):
                bj = [int(x, 16) for x in contained[j]["basis_hex"]]
                if binary_rank(bi + bj, 6) == 6:
                    complementary.append((i, j))
        permutation_screens.append({
            "non_bent_component_count": len(nb),
            "contained_three_spaces_count": len(contained),
            "complementary_pairs_count": len(complementary),
        })
    assert all(r["contained_three_spaces_count"] == 0 for r in permutation_screens)
    assert all(r["complementary_pairs_count"] == 0 for r in permutation_screens)
    for rec, screen in zip(direct, permutation_screens):
        rec["permutation_screen"] = screen

    # Gold cube signature for comparison.
    _, gold = build_q_and_lut(field, 0)
    gold_sig = orthoderivative_signature(gold)

    # Positive permutation-containing control: the Kim mapping in the Banff
    # field representation alpha^6+alpha^4+alpha^3+alpha+1=0.
    kim_field = BinaryField(6, 0x5B)
    alpha = 0x02
    kim = [
        kim_field.pow(x, 3)
        ^ kim_field.pow(x, 10)
        ^ kim_field.mul(alpha, kim_field.pow(x, 24))
        for x in range(64)
    ]
    kim_gold = [kim_field.pow(x, 3) for x in range(64)]
    assert is_apn_quadratic(kim)
    kim_amplitude, kim_nb, _ = component_data(kim)
    kim_spaces = []
    for basis, vectors in subspaces:
        if set(vectors) <= kim_nb:
            kim_spaces.append(basis)
    kim_complementary = []
    for i in range(len(kim_spaces)):
        for j in range(i + 1, len(kim_spaces)):
            if binary_rank(list(kim_spaces[i]) + list(kim_spaces[j]), 6) == 6:
                kim_complementary.append((i, j))
    assert len(kim_spaces) == 3 and len(kim_complementary) == 3
    kim_delta = [kim[x] ^ kim_gold[x] for x in range(64)]
    kim_anf = vector_anf_coefficients(kim_delta, 6)
    kim_quadratic_coefficients = [
        kim_anf[mask] for mask in range(64) if mask.bit_count() == 2
    ]
    kim_difference_coefficient_rank = binary_rank(kim_quadratic_coefficients, 6)
    assert kim_difference_coefficient_rank == 6

    # The exact ODDS row of Banff class #2, included as a comparison target.
    banff_class_2_odds = {"0": 2583, "2": 1008, "4": 378, "8": 63}
    assert signatures[0]["differential_spectrum"] == banff_class_2_odds

    return {
        "n": 6,
        "field_modulus_hex": "0x43",
        "GF4_hex": [f"0x{x:02x}" for x in f4],
        "safe_coefficients_hex": [f"0x{x:02x}" for x in safe],
        "exact_safe_conditions_equivalent": [
            "theta^6+theta^3+1=0",
            "theta^3 in GF(4)\\GF(2)",
            "multiplicative_order(theta)=9",
        ],
        "safe_root_polynomial": poly_terms(coeffs),
        "safe_root_polynomial_coefficients_low_to_high": coeffs,
        "multiplicative_orders": {f"0x{k:02x}": v for k, v in orders.items()},
        "Frobenius_orbit_hex": [f"0x{x:02x}" for x in orbit],
        "GF4_star_cosets_hex": [[f"0x{x:02x}" for x in c] for c in cosets],
        "all_six_switches_in_one_EA_orbit": True,
        "EA_orbit_reason": "Frobenius conjugation sends theta to a Frobenius conjugate; Q_lambda changes only by an affine-linear Boolean term because lambda^2=lambda+1.",
        "direct_records": direct,
        "common_orthoderivative_signature": signatures[0],
        "Gold_x3_orthoderivative_signature": gold_sig,
        "Banff_identification": {
            "matched_class_number": 2,
            "matched_representative": "x^3 + alpha^11*x^6 + alpha*x^9",
            "primitive_polynomial_for_alpha": "x^6+x^4+x^3+x+1",
            "matching_invariant": "orthoderivative differential spectrum (ODDS)",
            "matched_ODDS": banff_class_2_odds,
            "not_Kim_Dublin_class": True,
        },
        "permutation_screen": {
            "direct_permutations": 0,
            "non_bent_component_counts": [len(nb) for nb in nb_sets],
            "three_dimensional_subspaces_enumerated_per_function": len(subspaces),
            "contained_three_spaces_counts": [r["contained_three_spaces_count"] for r in permutation_screens],
            "complementary_pairs_counts": [r["complementary_pairs_count"] for r in permutation_screens],
            "passes_necessary_CCZ_permutation_test": False,
            "interpretation": "The whole rank-one family is ruled out from containing a permutation in its CCZ class by the adopted NB-subspace necessary condition.",
        },
        "Kim_permutation_class_positive_control": {
            "representative": "kappa(x)=x^3+x^10+alpha*x^24 over alpha^6+alpha^4+alpha^3+alpha+1=0",
            "APN_directly_verified": True,
            "direct_permutation": len(set(kim)) == 64,
            "component_amplitude_distribution": {str(k): v for k, v in sorted(kim_amplitude.items())},
            "non_bent_component_count": len(kim_nb),
            "contained_three_spaces_count": len(kim_spaces),
            "complementary_pairs_count": len(kim_complementary),
            "contained_three_space_bases_hex": [[f"0x{x:02x}" for x in basis] for basis in kim_spaces],
            "difference_from_Gold_coefficient_rank": kim_difference_coefficient_rank,
            "interpretation": "Unlike the rank-one class #2 family, the Kim class passes the subspace geometry test. Its standard difference from x^3 has full output coefficient rank 6, indicating that a rank-one search cannot reach it in this alignment.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/small_dimension_trace_switch_structure.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    output = args.output if args.output.is_absolute() else root / args.output

    result = {
        "schema": "small-dimension-trace-switch-structure-v1",
        "date": "2026-08-02",
        "family": "F_theta(x)=x^3+theta*Tr(x)*Tr(lambda*x), lambda in GF(4)\\GF(2)",
        "dimension_4": n4_record(),
        "dimension_6": n6_record(),
        "summary": {
            "n4": "Seven safe coefficients form ker(absolute trace)\\{0}; every switch is EA-equivalent to x^3.",
            "n6": "Six safe coefficients are exactly the elements of order 9; all form one Frobenius/EA orbit and match Banff quadratic APN class #2, not the permutation-containing Kim class.",
        },
        "validation": {
            "all_n4_nonzero_coefficients_structurally_classified": True,
            "all_n6_nonzero_coefficients_structurally_classified": True,
            "all_reported_safe_functions_directly_APN_checked": True,
            "n4_full_GL4_EA_witness_search_completed": True,
            "n6_all_1395_three_spaces_enumerated": True,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "sha256": sha256(output),
        "n4_safe": len(result["dimension_4"]["safe_coefficients_hex"]),
        "n6_safe": len(result["dimension_6"]["safe_coefficients_hex"]),
        "n6_contained_3spaces": result["dimension_6"]["permutation_screen"]["contained_three_spaces_counts"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
