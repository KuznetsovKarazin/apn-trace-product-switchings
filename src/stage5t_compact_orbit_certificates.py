#!/usr/bin/env python3
"""Stage 5T: compact Frobenius-orbit certificates for n=6 and n=8.

This replaces the raw 4^4 and 64^2 section tables by one explicit section
witness for each Frobenius orbit of represented coefficients.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

# ---------- GF(4) and GF(64)=GF(4)[u]/(u^3-lambda) ----------

def e_mul(a: int, b: int) -> int:
    out = 0
    aa, bb = a, b
    while bb:
        if bb & 1:
            out ^= aa
        bb >>= 1
        aa <<= 1
        if aa & 0b100:
            aa ^= 0b111
    return out & 0b11

LAMBDA = 0b10
LAMBDA2 = 0b11


def k6_add(x, y):
    return tuple(a ^ b for a, b in zip(x, y))


def k6_mul(x, y):
    raw = [0] * 5
    for i, a in enumerate(x):
        for j, b in enumerate(y):
            raw[i + j] ^= e_mul(a, b)
    raw[1] ^= e_mul(raw[4], LAMBDA)
    raw[0] ^= e_mul(raw[3], LAMBDA)
    return tuple(raw[:3])


def k6_pow(x, exponent: int):
    out = (1, 0, 0)
    base = x
    e = exponent
    while e:
        if e & 1:
            out = k6_mul(out, base)
        base = k6_mul(base, base)
        e >>= 1
    return out


def k6_order(x) -> int:
    if x == (0, 0, 0):
        return 0
    y = (1, 0, 0)
    for i in range(1, 64):
        y = k6_mul(y, x)
        if y == (1, 0, 0):
            return i
    raise AssertionError


def k6_code(x) -> int:
    return x[0] | (x[1] << 2) | (x[2] << 4)


def k6_cube_formula(x):
    a, b, c = x
    sq = lambda z: e_mul(z, z)
    cube = lambda z: e_mul(sq(z), z)
    return (
        cube(a) ^ e_mul(LAMBDA, cube(b)) ^ e_mul(LAMBDA2, cube(c)),
        e_mul(sq(a), b) ^ e_mul(LAMBDA, e_mul(a, sq(c))) ^ e_mul(LAMBDA, e_mul(sq(b), c)),
        e_mul(a, sq(b)) ^ e_mul(sq(a), c) ^ e_mul(LAMBDA, e_mul(b, sq(c))),
    )


def k6_poly_mul(p, q):
    out = [(0, 0, 0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] = k6_add(out[i + j], k6_mul(a, b))
    return out


def k6_minpoly_exponents(x):
    orbit = []
    y = x
    while y not in orbit:
        orbit.append(y)
        y = k6_pow(y, 2)
    poly = [(1, 0, 0)]
    for root in orbit:
        poly = k6_poly_mul(poly, [root, (1, 0, 0)])
    assert all(c in ((0, 0, 0), (1, 0, 0)) for c in poly)
    return [i for i, c in enumerate(poly) if c == (1, 0, 0)]


# ---------- GF(256) AES polynomial basis ----------
GF_MOD = 0x11B


def g8_mul(a: int, b: int) -> int:
    out = 0
    a &= 0xFF
    b &= 0xFF
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 0x100:
            a ^= GF_MOD
    return out & 0xFF


def g8_pow(a: int, exponent: int) -> int:
    out = 1
    base = a & 0xFF
    e = exponent
    while e:
        if e & 1:
            out = g8_mul(out, base)
        base = g8_mul(base, base)
        e >>= 1
    return out


def g8_abs_trace(a: int) -> int:
    out = 0
    t = a
    for _ in range(8):
        out ^= t
        t = g8_mul(t, t)
    assert out in (0, 1)
    return out


def g8_rel_trace(a: int) -> int:
    return a ^ g8_pow(a, 4) ^ g8_pow(a, 16) ^ g8_pow(a, 64)


def g8_poly_mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] ^= g8_mul(a, b)
    return out


def g8_minpoly_exponents(x: int):
    orbit = []
    y = x
    while y not in orbit:
        orbit.append(y)
        y = g8_pow(y, 2)
    poly = [1]
    for root in orbit:
        poly = g8_poly_mul(poly, [root, 1])
    assert all(c in (0, 1) for c in poly)
    return [i for i, c in enumerate(poly) if c]


def choose_simple_pair(pairs):
    return min(pairs, key=lambda z: (z[0].bit_count() + z[1].bit_count(), z))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_n6():
    witnesses = {}
    multiplicities = Counter()
    for p in range(4):
        for q in range(4):
            a = (1, p, q)
            a3 = k6_cube_formula(a)
            for r in range(4):
                for s in range(4):
                    b = (0, r, s)
                    theta = k6_add(a3, k6_cube_formula(b))
                    multiplicities[theta] += 1
                    witnesses.setdefault(theta, []).append((p, q, r, s))
    assert (0, 0, 0) not in witnesses
    safe = []
    represented = []
    for a in range(4):
        for b in range(4):
            for c in range(4):
                x = (a, b, c)
                if x == (0, 0, 0):
                    continue
                (represented if x in witnesses else safe).append(x)
    assert len(safe) == 6 and all(k6_order(x) == 9 for x in safe)

    unseen = set(represented)
    rows = []
    while unseen:
        rep = min(unseen, key=k6_code)
        orbit = []
        y = rep
        while y not in orbit:
            orbit.append(y)
            y = k6_pow(y, 2)
        unseen -= set(orbit)
        witness = min(
            witnesses[rep],
            key=lambda z: (sum(v != 0 for v in z), z),
        )
        p, q, r, s = witness
        a = (1, p, q)
        b = (0, r, s)
        assert k6_add(k6_pow(a, 3), k6_pow(b, 3)) == rep
        rows.append({
            "representative_coordinates": list(rep),
            "representative_code": k6_code(rep),
            "orbit_size": len(orbit),
            "multiplicative_order": k6_order(rep),
            "minimal_polynomial_exponents": k6_minpoly_exponents(rep),
            "section_witness": {"a": list(a), "b": list(b)},
            "representation_multiplicity": multiplicities[rep],
        })
    assert len(rows) == 12
    return {
        "field": "GF(64)=GF(4)[u]/(u^3-lambda)",
        "safe_polynomial_exponents": [0, 3, 6],
        "safe_description": "six roots of X^6+X^3+1, equivalently elements of order 9",
        "represented_orbit_count": len(rows),
        "represented_element_count": len(represented),
        "rows": rows,
    }


def build_n8():
    f4 = {x for x in range(256) if g8_pow(x, 4) == x}
    witnesses = {}
    multiplicities = Counter()
    for r in range(256):
        if g8_rel_trace(r) != 1:
            continue
        for s in range(256):
            if g8_rel_trace(s) != 0:
                continue
            theta = g8_pow(r, 3) ^ g8_pow(s, 3)
            multiplicities[theta] += 1
            witnesses.setdefault(theta, []).append((r, s))
    assert set(witnesses) == set(range(256)) - f4

    unseen = set(witnesses)
    rows = []
    while unseen:
        rep = min(unseen)
        orbit = []
        y = rep
        while y not in orbit:
            orbit.append(y)
            y = g8_pow(y, 2)
        unseen -= set(orbit)
        r, s = choose_simple_pair(witnesses[rep])
        assert g8_rel_trace(r) == 1 and g8_rel_trace(s) == 0
        assert g8_pow(r, 3) ^ g8_pow(s, 3) == rep
        rows.append({
            "representative": rep,
            "representative_hex": f"0x{rep:02x}",
            "orbit_size": len(orbit),
            "minimal_polynomial_exponents": g8_minpoly_exponents(rep),
            "section_witness": {"r_hex": f"0x{r:02x}", "s_hex": f"0x{s:02x}"},
            "representation_multiplicity": multiplicities[rep],
        })
    assert len(rows) == 33
    return {
        "field": "GF(256)=GF(2)[z]/(z^8+z^4+z^3+z+1)",
        "safe_subfield": sorted(f4),
        "safe_polynomial_exponents": [1, 4],
        "safe_description": "GF(4), roots of X^4+X",
        "represented_orbit_count": len(rows),
        "represented_element_count": len(witnesses),
        "rows": rows,
    }


def main():
    root = Path(__file__).resolve().parents[1]
    output = root / "results/current/compact_section_orbit_certificates.json"
    result = {
        "schema": "compact-section-orbit-certificates-v1",
        "date": "2026-08-02",
        "principle": (
            "The section-product set is Frobenius invariant. One explicit "
            "section witness for each Frobenius orbit proves coverage of the "
            "whole orbit."
        ),
        "n6": build_n6(),
        "n8": build_n8(),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "sha256": sha256(output),
        "n6_rows": len(result["n6"]["rows"]),
        "n8_rows": len(result["n8"]["rows"]),
    }, indent=2))


if __name__ == "__main__":
    main()
