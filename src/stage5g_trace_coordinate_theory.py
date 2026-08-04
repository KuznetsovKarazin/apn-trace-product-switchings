#!/usr/bin/env python3
"""Stage 5G: trace-coordinate normal form and exact certificates.

This script identifies q0 and q1 as rank-two trace switches

    Q_c(x) = Tr_K/F2(c x) Tr_K/F2(lambda c x),

with lambda in GF(4)^*, normalizes the x^3 and x^9 switchings, and builds
machine-checkable certificates for the cube-fibre and trace-selector laws.

The rank-one sufficiency certificate is symbolic: a univariate resultant and
Bezout identity over GF(2) prove that the product of a two-dimensional F2
subspace transverse to ker Tr_K/GF(4) cannot lie in GF(4).

The rank-one necessity and rank-two trace selector are recorded as exact
finite-geometry certificates over GF(256), not as heuristic scans.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import sympy as sp

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


def gf_inv(a: int) -> int:
    if a == 0:
        raise ZeroDivisionError
    return gf_pow(a, 254)


def gf_order(a: int) -> int:
    if a == 0:
        return 0
    value = 1
    for order in range(1, 256):
        value = gf_mult(value, a)
        if value == 1:
            return order
    raise AssertionError


def absolute_trace(a: int) -> int:
    value = 0
    term = a & 0xFF
    for _ in range(8):
        value ^= term
        term = gf_mult(term, term)
    if value not in (0, 1):
        raise AssertionError
    return value


def relative_trace_to_f4(a: int) -> int:
    return a ^ gf_pow(a, 4) ^ gf_pow(a, 16) ^ gf_pow(a, 64)


def trace_f4_to_f2(a: int) -> int:
    value = a ^ gf_pow(a, 2)
    if value not in (0, 1):
        raise AssertionError(f"not an F4 trace value: 0x{value:02x}")
    return value


def parity(a: int) -> int:
    return a.bit_count() & 1


def polynomial_basis_dual() -> list[int]:
    # M[j][k] = Tr(z^j z^k); solve M c_i=e_i over GF(2).
    matrix = [
        [absolute_trace(gf_mult(1 << j, 1 << k)) for k in range(8)]
        for j in range(8)
    ]
    aug = [row[:] + [1 if i == j else 0 for j in range(8)]
           for i, row in enumerate(matrix)]
    row = 0
    for col in range(8):
        pivot = next(i for i in range(row, 8) if aug[i][col])
        aug[row], aug[pivot] = aug[pivot], aug[row]
        for i in range(8):
            if i != row and aug[i][col]:
                aug[i] = [x ^ y for x, y in zip(aug[i], aug[row])]
        row += 1
    inverse = [r[8:] for r in aug]
    dual = []
    for i in range(8):
        # column i of M^{-1}
        value = sum(inverse[k][i] << k for k in range(8))
        dual.append(value)
    for i, d in enumerate(dual):
        assert all(
            absolute_trace(gf_mult(d, x)) == ((x >> i) & 1)
            for x in range(N)
        )
    return dual


def xor_elements(values: Iterable[int]) -> int:
    out = 0
    for value in values:
        out ^= value
    return out


def linear_form_coefficient(dual: list[int], bit_indices: list[int]) -> int:
    return xor_elements(dual[i] for i in bit_indices)


def q_product_truth(a: int, b: int) -> list[int]:
    return [
        absolute_trace(gf_mult(a, x))
        & absolute_trace(gf_mult(b, x))
        for x in range(N)
    ]


def q_edge_truth(edges: list[tuple[int, int]]) -> list[int]:
    return [
        xor_elements((((x >> i) & 1) & ((x >> j) & 1)) for i, j in edges)
        for x in range(N)
    ]


def polar_trace_switch(c: int, a: int, x: int) -> int:
    ta = relative_trace_to_f4(gf_mult(c, a))
    tx = relative_trace_to_f4(gf_mult(c, x))
    return trace_f4_to_f2(gf_mult(ta, gf_pow(tx, 2)))


def polar_from_truth(truth: list[int], a: int, x: int) -> int:
    return truth[x ^ a] ^ truth[x] ^ truth[a]


def subspaces_2d() -> list[tuple[int, int, int, int]]:
    spaces: set[tuple[int, int, int, int]] = set()
    for a in range(1, N):
        for x in range(1, N):
            if x == a:
                continue
            spaces.add(tuple(sorted((0, a, x, a ^ x))))
    result = sorted(spaces)
    assert len(result) == 10795
    return result


def subspace_product(space: tuple[int, int, int, int]) -> int:
    product = 1
    for value in space[1:]:
        product = gf_mult(product, value)
    return product


def scale_space(space: tuple[int, int, int, int], scalar: int) -> tuple[int, int, int, int]:
    return tuple(sorted(gf_mult(scalar, value) for value in space))


def roots_polynomial(roots: Iterable[int]) -> list[int]:
    # Coefficients low-to-high over GF(256).
    coeffs = [1]
    for root in sorted(roots):
        nxt = [0] * (len(coeffs) + 1)
        for i, coeff in enumerate(coeffs):
            nxt[i] ^= gf_mult(coeff, root)
            nxt[i + 1] ^= coeff
        coeffs = nxt
    return coeffs


def poly_exponents(poly: sp.Poly) -> list[int]:
    return [int(monomial[0]) for monomial, coeff in poly.terms() if int(coeff) & 1]


def symbolic_resultant_certificate() -> dict[str, Any]:
    r, s = sp.symbols("r s")
    trace_poly = lambda z: z + z**4 + z**16 + z**64
    # If T(r)=1, T(s)=0 and r^3+s^3 lies in F4, these polynomials have a common root.
    resultant_expr = sp.resultant(
        trace_poly(s),
        s**12 + s**3 + r**12 + r**3,
        s,
    )
    resultant = sp.Poly(resultant_expr, r, modulus=2)
    trace_one = sp.Poly(trace_poly(r) + 1, r, modulus=2)
    bezout_a, bezout_b, gcd_poly = sp.gcdex(trace_one, resultant)
    assert gcd_poly == sp.Poly(1, r, modulus=2)
    assert bezout_a * trace_one + bezout_b * resultant == gcd_poly
    payload = {
        "variable": "r",
        "trace_one_polynomial_exponents": poly_exponents(trace_one),
        "resultant_degree": resultant.degree(),
        "resultant_term_count": len(resultant.terms()),
        "resultant_exponents": poly_exponents(resultant),
        "bezout_a_exponents": poly_exponents(bezout_a),
        "bezout_b_exponents": poly_exponents(bezout_b),
        "gcd": 1,
        "identity": "A(r)*(T(r)+1)+B(r)*Res_s(T(s),s^12+s^3+r^12+r^3)=1 over GF(2)",
        "consequence": "No section product r^3+s^3 with T(r)=1,T(s)=0 lies in GF(4).",
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["certificate_sha256"] = hashlib.sha256(encoded).hexdigest()
    return payload



def f4_coordinates(value: int, lambda_value: int) -> tuple[int, int]:
    for a in (0, 1):
        for b in (0, 1):
            if (a ^ (lambda_value if b else 0)) == value:
                return a, b
    raise ValueError(f"0x{value:02x} is not in GF(4)")


def apply_f4_linear_map(
    signature: tuple[int, int], value: int, lambda_value: int
) -> int:
    a, b = f4_coordinates(value, lambda_value)
    return (signature[0] if a else 0) ^ (signature[1] if b else 0)


def compose_f4_linear_maps(
    left: tuple[int, int], right: tuple[int, int], lambda_value: int
) -> tuple[int, int]:
    return (
        apply_f4_linear_map(left, right[0], lambda_value),
        apply_f4_linear_map(left, right[1], lambda_value),
    )


def f4_linear_map_order(signature: tuple[int, int], lambda_value: int) -> int:
    identity = (1, lambda_value)
    power = identity
    for order in range(1, 7):
        power = compose_f4_linear_maps(signature, power, lambda_value)
        if power == identity:
            return order
    raise AssertionError("invalid GL(2,2) map")

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kernels",
        type=Path,
        default=Path("results/evidence/canonical_switching_kernels.json"),
    )
    parser.add_argument(
        "--locus",
        type=Path,
        default=Path("results/current/canonical_low_rank_locus.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/trace_coordinate_theory.json"),
    )
    parser.add_argument(
        "--certificate-output",
        type=Path,
        default=Path("results/current/trace_coordinate_bezout_certificate.json"),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    resolve = lambda p: p if p.is_absolute() else root / p
    kernels_path = resolve(args.kernels)
    locus_path = resolve(args.locus)
    output_path = resolve(args.output)
    certificate_path = resolve(args.certificate_output)

    kernels = json.loads(kernels_path.read_text(encoding="utf-8"))
    locus = json.loads(locus_path.read_text(encoding="utf-8"))
    dual = polynomial_basis_dual()

    # Factor choices are read from the formulas in canonical_switching_kernels.json.
    factor_bits = {
        "q0": ([0, 5], [3]),
        "q1": ([0, 2], [1, 6, 7]),
    }
    trace_forms: dict[str, dict[str, Any]] = {}
    for name, (bits_c, bits_other) in factor_bits.items():
        c = linear_form_coefficient(dual, bits_c)
        other = linear_form_coefficient(dual, bits_other)
        ratio = gf_mult(other, gf_inv(c))
        edges = [tuple(edge) for edge in kernels["kernels"][name]["expanded_terms"]]
        edge_truth = q_edge_truth(edges)
        trace_truth = q_product_truth(c, other)
        assert edge_truth == trace_truth
        assert ratio not in (0, 1)
        assert gf_pow(ratio, 4) == ratio
        assert gf_order(ratio) == 3
        trace_forms[name] = {
            "c": c,
            "c_hex": f"0x{c:02x}",
            "lambda_c": other,
            "lambda_c_hex": f"0x{other:02x}",
            "lambda": ratio,
            "lambda_hex": f"0x{ratio:02x}",
            "factor_bits": [bits_c, bits_other],
            "truth_table_verified": True,
        }

    lambda_value = trace_forms["q0"]["lambda"]
    assert trace_forms["q1"]["lambda"] == lambda_value
    assert gf_pow(lambda_value, 2) == (lambda_value ^ 1)

    c0 = trace_forms["q0"]["c"]
    c1 = trace_forms["q1"]["c"]
    fifth_root = gf_mult(c1, gf_inv(c0))
    assert gf_order(fifth_root) == 5
    assert gf_pow(fifth_root, 16) == fifth_root
    fifth_root_sq = gf_pow(fifth_root, 2)

    # Verify the relative-trace polar formula for every a,x for q0 and q1.
    polar_checks = 0
    for name, record in trace_forms.items():
        c = record["c"]
        edges = [tuple(edge) for edge in kernels["kernels"][name]["expanded_terms"]]
        truth = q_edge_truth(edges)
        for a in range(N):
            for x in range(N):
                assert polar_from_truth(truth, a, x) == polar_trace_switch(c, a, x)
                polar_checks += 1

    f4 = sorted(value for value in range(N) if gf_pow(value, 4) == value)
    f4_nonzero = [value for value in f4 if value]
    assert set(f4_nonzero) == {1, lambda_value, gf_pow(lambda_value, 2)}

    # Rank-one normalization: x -> c*x, output -> c^3*output.
    x3_rank_one = {
        row["kernel"]: row["accepted_coefficients"]
        for row in locus["loci"]["x3"]["rank_one"]["by_kernel"]
    }
    normalized_rank_one: dict[str, Any] = {}
    for name, c in (("q0", c0), ("q1", c1)):
        accepted = x3_rank_one[name]
        normalized = sorted(gf_mult(mu, gf_pow(c, 3)) for mu in accepted)
        assert normalized == f4_nonzero
        cube_constant = gf_pow(gf_inv(c), 9)
        assert all(gf_pow(mu, 3) == cube_constant for mu in accepted)
        normalized_rank_one[name] = {
            "c_hex": f"0x{c:02x}",
            "accepted_original_hex": [f"0x{x:02x}" for x in accepted],
            "accepted_normalized_hex": [f"0x{x:02x}" for x in normalized],
            "condition": f"mu*c^3 in GF(4)^*, equivalently mu^3=0x{cube_constant:02x}",
            "cube_constant_hex": f"0x{cube_constant:02x}",
        }

    # Exact section-product characterization of the forbidden set for Q_1.
    trace_one = [x for x in range(N) if relative_trace_to_f4(x) == 1]
    trace_zero = [x for x in range(N) if relative_trace_to_f4(x) == 0]
    assert len(trace_one) == len(trace_zero) == 64
    section_product_multiplicity: Counter[int] = Counter()
    for r_value in trace_one:
        for s_value in trace_zero:
            section_product_multiplicity[gf_pow(r_value, 3) ^ gf_pow(s_value, 3)] += 1
    section_products = set(section_product_multiplicity)
    assert section_products == set(range(N)) - set(f4)
    product_poly = roots_polynomial(section_products)
    nonbinary = [i for i, coeff in enumerate(product_poly) if coeff not in (0, 1)]
    assert not nonbinary
    product_poly_exponents = [i for i, coeff in enumerate(product_poly) if coeff]
    assert product_poly_exponents == list(range(0, 253, 3))

    symbolic_certificate = symbolic_resultant_certificate()
    certificate_path.parent.mkdir(parents=True, exist_ok=True)
    certificate_path.write_text(
        json.dumps(symbolic_certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    # Normalize rank-two x^3 data.
    x3_rank2 = next(
        row for row in locus["loci"]["x3"]["rank_two"]["by_coefficient_map_kernel"]
        if row["coefficient_map_kernel"] == 4
    )["accepted_pairs"]
    normalized_pairs = []
    for row in x3_rank2:
        eta = gf_mult(row["coefficient_1"], gf_pow(c0, 3))
        second = gf_mult(row["coefficient_2"], gf_pow(c0, 3))
        assert second == gf_mult(eta, fifth_root_sq)
        normalized_pairs.append((eta, second))
    normalized_etas = sorted(eta for eta, _ in normalized_pairs)
    expected_etas = sorted(
        eta for eta in f4_nonzero
        if trace_f4_to_f2(gf_mult(lambda_value, eta)) == 1
    )
    assert normalized_etas == expected_etas

    # Finite-geometry certificate for the excluded third lift.
    spaces = subspaces_2d()
    rank2_bad_counts: dict[int, int] = {}
    rank2_bad_spaces: dict[int, list[tuple[int, int, int, int]]] = {}
    for eta in f4_nonzero:
        bad: list[tuple[int, int, int, int]] = []
        for space in spaces:
            a, x = space[1], space[2]
            b1 = polar_trace_switch(1, a, x)
            br = polar_trace_switch(fifth_root, a, x)
            update = (eta if b1 else 0) ^ (
                gf_mult(eta, fifth_root_sq) if br else 0
            )
            if update and subspace_product(space) == update:
                bad.append(space)
        rank2_bad_counts[eta] = len(bad)
        rank2_bad_spaces[eta] = bad
    safe_etas = sorted(eta for eta, count in rank2_bad_counts.items() if count == 0)
    excluded_etas = [eta for eta, count in rank2_bad_counts.items() if count]
    assert safe_etas == expected_etas
    assert excluded_etas == [gf_pow(lambda_value, 2)]
    excluded_eta = excluded_etas[0]
    excluded_spaces = rank2_bad_spaces[excluded_eta]
    assert len(excluded_spaces) == 12
    excluded_updates = {
        subspace_product(space) for space in excluded_spaces
    }
    assert len(excluded_updates) == 1
    excluded_update = next(iter(excluded_updates))
    assert excluded_update == gf_mult(
        excluded_eta, 1 ^ fifth_root_sq
    )

    # Partition the 12 witnesses under GF(4)^* scaling.
    remaining = set(excluded_spaces)
    scaling_orbits = []
    while remaining:
        representative = min(remaining)
        orbit = sorted({scale_space(representative, eta) for eta in f4_nonzero})
        assert set(orbit) <= set(excluded_spaces)
        scaling_orbits.append(orbit)
        remaining -= set(orbit)
    assert len(scaling_orbits) == 4
    assert all(len(orbit) == 3 for orbit in scaling_orbits)

    # Refine the y=(1,1) obstruction by the transition map
    # L_V = T_r o (T|_V)^(-1) in GL(2,2).
    transition_records: dict[tuple[int, int], dict[str, Any]] = {}
    target_line = {
        eta: gf_mult(eta, 1 ^ fifth_root_sq) for eta in f4_nonzero
    }
    for space in spaces:
        a, x = space[1], space[2]
        if polar_trace_switch(1, a, x) != 1:
            continue
        if polar_trace_switch(fifth_root, a, x) != 1:
            continue
        inverse_section = {
            relative_trace_to_f4(value): value for value in space
        }
        assert len(inverse_section) == 4
        signature = (
            relative_trace_to_f4(gf_mult(fifth_root, inverse_section[1])),
            relative_trace_to_f4(
                gf_mult(fifth_root, inverse_section[lambda_value])
            ),
        )
        record = transition_records.setdefault(signature, {
            "subspace_count": 0,
            "products": Counter(),
        })
        record["subspace_count"] += 1
        record["products"][subspace_product(space)] += 1
    assert len(transition_records) == 6
    transition_summary = []
    involution_witness_total = 0
    for signature, record in sorted(transition_records.items()):
        order = f4_linear_map_order(signature, lambda_value)
        intersections = {
            f"0x{eta:02x}": record["products"][target]
            for eta, target in sorted(target_line.items())
        }
        if order == 2:
            involution_witness_total += intersections[f"0x{excluded_eta:02x}"]
        transition_summary.append({
            "images_of_1_lambda_hex": [
                f"0x{signature[0]:02x}", f"0x{signature[1]:02x}"
            ],
            "order": order,
            "subspace_count": record["subspace_count"],
            "distinct_product_count": len(record["products"]),
            "target_line_witness_counts": intersections,
        })
    assert all(row["subspace_count"] == 256 for row in transition_summary)
    assert involution_witness_total == 12
    assert all(
        row["target_line_witness_counts"][f"0x{excluded_eta:02x}"] == (4 if row["order"] == 2 else 0)
        for row in transition_summary
    )

    # Normalize the unique x^9 rank-two pair and verify componentwise cubing.
    x9_rank2 = next(
        row for row in locus["loci"]["x9"]["rank_two"]["by_coefficient_map_kernel"]
        if row["coefficient_map_kernel"] == 4
    )["accepted_pairs"]
    assert len(x9_rank2) == 1
    x9_first = gf_mult(x9_rank2[0]["coefficient_1"], gf_pow(c0, 9))
    x9_second = gf_mult(x9_rank2[0]["coefficient_2"], gf_pow(c0, 9))
    assert (x9_first, x9_second) == (1, fifth_root)
    for eta, second in normalized_pairs:
        assert gf_pow(eta, 3) == 1
        assert gf_pow(second, 3) == fifth_root

    report = {
        "schema": "trace-coordinate-theory-v1",
        "date": "2026-08-02",
        "field": {
            "representation": "GF(2)[z]/(0x11B)",
            "polynomial_basis_dual_hex": [f"0x{x:02x}" for x in dual],
            "GF4_hex": [f"0x{x:02x}" for x in f4],
            "lambda_hex": f"0x{lambda_value:02x}",
            "lambda_order": gf_order(lambda_value),
        },
        "trace_switch_normal_form": {
            "definition": "Q_c(x)=Tr_K/F2(c*x)*Tr_K/F2(lambda*c*x)",
            "polar": "B_c(a,x)=Tr_GF4/F2(T(ca)*T(cx)^2), T=Tr_K/GF4",
            "q0": trace_forms["q0"],
            "q1": trace_forms["q1"],
            "c1_over_c0_hex": f"0x{fifth_root:02x}",
            "c1_over_c0_order": gf_order(fifth_root),
            "c1_over_c0_in_GF16": gf_pow(fifth_root, 16) == fifth_root,
            "polar_identity_checks": polar_checks,
        },
        "rank_one_theorem_dimension_8": {
            "normalized_statement": "x^3+theta*Q_1 is APN iff theta in GF(4)^*.",
            "coordinate_statement": "x^3+mu*Q_c is APN iff mu*c^3 in GF(4)^*, equivalently mu^3=c^(-9).",
            "families": normalized_rank_one,
            "derivative_geometry": {
                "bad_subspace_condition": "A bad derivative corresponds to a 2D F2-subspace V with product prod(V)=theta and rank(T|V)=2.",
                "section_parameterization": "If rank(T|V)=2, the inverse section has phi(y)=r*y+s*y^2 with T(r)=1,T(s)=0 and prod(V)=r^3+s^3.",
                "section_count": len(trace_one) * len(trace_zero),
                "distinct_section_products": len(section_products),
                "section_product_set": "GF(256) \\ GF(4)",
                "multiplicity_distribution": {
                    str(k): v for k, v in sorted(Counter(section_product_multiplicity.values()).items())
                },
                "squarefree_product_polynomial": "product_{u notin GF(4)}(X+u)=sum_{j=0}^{84} X^(3j)",
                "squarefree_product_exponents": product_poly_exponents,
            },
            "symbolic_sufficiency_certificate_file": str(certificate_path.relative_to(root)),
            "symbolic_sufficiency_certificate_sha256": hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
            "proof_status": {
                "sufficiency_theta_in_GF4_star": "symbolic resultant/Bezout certificate",
                "necessity_theta_notin_GF4": "exact 4096-section finite-geometry image certificate",
            },
        },
        "rank_two_q0_q1_theory": {
            "normalized_x3_family": "x^3 + eta*(Q_1 + r^2*Q_r), r=0xb0, eta in GF(4)^*",
            "r_hex": f"0x{fifth_root:02x}",
            "r_order": gf_order(fifth_root),
            "r_squared_hex": f"0x{fifth_root_sq:02x}",
            "accepted_eta_condition": "Tr_GF4/F2(lambda*eta)=1",
            "accepted_eta_hex": [f"0x{x:02x}" for x in safe_etas],
            "excluded_eta_hex": f"0x{excluded_eta:02x}",
            "bad_subspace_counts": {
                f"0x{eta:02x}": count for eta, count in sorted(rank2_bad_counts.items())
            },
            "excluded_update_hex": f"0x{excluded_update:02x}",
            "excluded_witness_count": len(excluded_spaces),
            "excluded_witness_GF4_scaling_orbits": [
                [[f"0x{x:02x}" for x in space] for space in orbit]
                for orbit in scaling_orbits
            ],
            "trace_transition_map_classification": transition_summary,
            "trace_transition_interpretation": "All 12 excluded-lift witnesses occur for the three order-2 elements of GL(2,2), four witnesses per involution; identity and order-3 transition maps contribute none.",
            "normalized_x9_point": "x^9 + Q_1 + r*Q_r",
            "componentwise_cube": "(eta,eta*r^2)^3=(1,r) for every eta in GF(4)^*; APN x^3 lifts are selected by Tr(lambda*eta)=1.",
            "proof_status": "exact 10795-subspace incidence certificate; a compact symbolic proof of the trace selector remains open",
        },
        "interpretation": [
            "The cube fibres are forced by trace-coordinate normalization, not accidental field constants.",
            "The order-5 slopes arise because q1 is Q_r after normalizing q0, with r in GF(16)^* of order 5.",
            "The excluded third cube lift is obstructed by exactly 12 bad derivative subspaces, forming four GF(4)^*-scaling orbits.",
            "The unique x^9 rank-two point is the componentwise cube image of the three formal x^3 lifts; the trace selector removes one lift before cubing.",
        ],
        "provenance": {
            "kernels_file": str(kernels_path.relative_to(root)),
            "kernels_sha256": hashlib.sha256(kernels_path.read_bytes()).hexdigest(),
            "locus_file": str(locus_path.relative_to(root)),
            "locus_sha256": hashlib.sha256(locus_path.read_bytes()).hexdigest(),
        },
        "validation": {
            "trace_truth_tables_verified": True,
            "relative_trace_polar_verified_all_pairs": True,
            "rank_one_section_image_exact": True,
            "symbolic_resultant_gcd_one": True,
            "rank_two_trace_selector_exact": True,
            "x9_componentwise_cube_exact": True,
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": str(output_path),
        "sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
        "certificate": str(certificate_path),
        "certificate_sha256": hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
        "q0_c": f"0x{c0:02x}",
        "q1_c": f"0x{c1:02x}",
        "lambda": f"0x{lambda_value:02x}",
        "r": f"0x{fifth_root:02x}",
        "section_products": len(section_products),
        "excluded_rank2_witnesses": len(excluded_spaces),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
