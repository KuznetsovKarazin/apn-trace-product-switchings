#!/usr/bin/env python3
r"""Stage 5L: permutation-oriented screening of the eight projective switches.

For the eight marked APN switchings from Stage 5J, compute:
  * direct image/preimage profile;
  * component Walsh-amplitude distribution and non-bent set NB(F);
  * all 4-dimensional F2-subspaces contained in NB(F) union {0};
  * whether two such subspaces are complementary.

The final test is a necessary, not sufficient, condition for CCZ-equivalence
to a permutation: NB(F) union {0} must contain V,W of dimension n/2 with
V direct-sum W = F2^n.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np

from stage5g_trace_coordinate_theory import absolute_trace, gf_inv, gf_mult, gf_pow

N = 256
LAMBDA = 0xBD


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parity(x: int) -> int:
    return x.bit_count() & 1


def q_c(c: int, x: int) -> int:
    return absolute_trace(gf_mult(c, x)) & absolute_trace(
        gf_mult(gf_mult(LAMBDA, c), x)
    )


def fwht(values: np.ndarray) -> np.ndarray:
    out = values.astype(np.int64, copy=True)
    width = 1
    while width < out.size:
        blocks = out.reshape(-1, 2 * width)
        left = blocks[:, :width].copy()
        right = blocks[:, width:].copy()
        blocks[:, :width] = left + right
        blocks[:, width:] = left - right
        width *= 2
    return out


def binary_rank(rows: Iterable[int], width: int = 8) -> int:
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
    values = []
    for mask in range(1, 1 << len(basis)):
        x = 0
        for i, row in enumerate(basis):
            if (mask >> i) & 1:
                x ^= row
        values.append(x)
    assert len(set(values)) == (1 << len(basis)) - 1
    return tuple(sorted(values))


def rref_subspaces_4_in_8():
    """Yield canonical RREF bases and their 15 nonzero vectors."""
    total = 0
    columns = range(8)
    for pivots in itertools.combinations(columns, 4):
        nonpivots = [c for c in columns if c not in pivots]
        allowed = [
            (i, c)
            for c in nonpivots
            for i, p in enumerate(pivots)
            if p < c
        ]
        for assignment in range(1 << len(allowed)):
            rows = [1 << p for p in pivots]
            for bit, (i, c) in enumerate(allowed):
                if (assignment >> bit) & 1:
                    rows[i] |= 1 << c
            assert binary_rank(rows) == 4
            total += 1
            yield tuple(rows), span_nonzero(rows)
    assert total == 200787


def build_lut(rho: int, eta: int) -> list[int]:
    second = gf_mult(eta, gf_pow(gf_inv(rho), 3))
    return [
        gf_pow(x, 3)
        ^ (eta if q_c(1, x) else 0)
        ^ (second if q_c(rho, x) else 0)
        for x in range(N)
    ]


def is_apn_quadratic(lut: list[int]) -> bool:
    for a in range(1, N):
        fa = lut[a]
        zeros = sum(
            (lut[x ^ a] ^ lut[x] ^ fa) == 0
            for x in range(N)
        )
        if zeros != 2:
            return False
    return True


def function_record(rho: int, eta: int, lut: list[int]) -> dict:
    preimages = Counter(lut)
    profile = Counter(preimages.values())
    # Missing outputs have multiplicity zero.
    profile[0] = N - len(preimages)

    amplitude_distribution = Counter()
    non_bent = set()
    component_details = []
    for b in range(1, N):
        signs = np.fromiter(
            (1 if parity(b & lut[x]) == 0 else -1 for x in range(N)),
            dtype=np.int64,
            count=N,
        )
        spectrum = fwht(signs)
        nonzero = spectrum[spectrum != 0]
        amplitudes = sorted(set(abs(int(v)) for v in nonzero))
        assert len(amplitudes) == 1
        amp = amplitudes[0]
        amplitude_distribution[amp] += 1
        bent = len(nonzero) == N
        if not bent:
            non_bent.add(b)
        component_details.append((b, amp, int(np.count_nonzero(spectrum == 0))))

    assert sum(amplitude_distribution.values()) == 255
    return {
        "rho0_hex": f"0x{rho:02x}",
        "eta_hex": f"0x{eta:02x}",
        "second_coefficient_hex": f"0x{gf_mult(eta, gf_pow(gf_inv(rho), 3)):02x}",
        "APN_directly_verified": is_apn_quadratic(lut),
        "direct_permutation": len(preimages) == N,
        "image_size": len(preimages),
        "preimage_multiplicity_distribution": {
            str(k): v for k, v in sorted(profile.items())
        },
        "component_amplitude_distribution": {
            str(k): v for k, v in sorted(amplitude_distribution.items())
        },
        "bent_component_count": 255 - len(non_bent),
        "non_bent_component_count": len(non_bent),
        "non_bent_components_hex": [f"0x{x:02x}" for x in sorted(non_bent)],
        "component_details": component_details,
        "_non_bent_set": non_bent,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage5j",
        type=Path,
        default=Path("results/current/projective_displacement_theorem.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/projective_switch_permutation_screen.json"),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    resolve = lambda p: p if p.is_absolute() else root / p
    stage5j_path = resolve(args.stage5j)
    output_path = resolve(args.output)
    stage5j = json.loads(stage5j_path.read_text(encoding="utf-8"))

    params = []
    for orbit in stage5j["accepted_Frobenius_orbits"]:
        for row in orbit:
            params.append((int(row["rho0_hex"], 16), int(row["eta_equals_zeta_hex"], 16)))
    params = sorted(set(params))
    assert len(params) == 8

    records = []
    for rho, eta in params:
        records.append(function_record(rho, eta, build_lut(rho, eta)))

    # Enumerate the 4-spaces once and test against all eight NB sets.
    contained: list[list[dict]] = [[] for _ in records]
    nb_sets = [rec["_non_bent_set"] for rec in records]
    total_subspaces = 0
    for basis, vectors in rref_subspaces_4_in_8():
        total_subspaces += 1
        vector_set = set(vectors)
        for i, nb in enumerate(nb_sets):
            if vector_set <= nb:
                contained[i].append({
                    "basis_hex": [f"0x{x:02x}" for x in basis],
                    "vectors": vectors,
                })
    assert total_subspaces == 200787

    for rec, spaces in zip(records, contained):
        complementary_pairs = []
        for i in range(len(spaces)):
            basis_i = [int(x, 16) for x in spaces[i]["basis_hex"]]
            for j in range(i + 1, len(spaces)):
                basis_j = [int(x, 16) for x in spaces[j]["basis_hex"]]
                if binary_rank(basis_i + basis_j) == 8:
                    complementary_pairs.append({
                        "first_basis_hex": spaces[i]["basis_hex"],
                        "second_basis_hex": spaces[j]["basis_hex"],
                    })
        rec["permutation_necessary_condition"] = {
            "n_over_2_subspace_dimension": 4,
            "contained_4_spaces_count": len(spaces),
            "contained_4_space_bases_hex": [s["basis_hex"] for s in spaces],
            "complementary_pairs_count": len(complementary_pairs),
            "complementary_pairs": complementary_pairs,
            "passes_two_complementary_subspaces_test": bool(complementary_pairs),
            "interpretation": "Passing is necessary but not sufficient for CCZ-equivalence to a permutation; failing rules it out.",
        }
        rec.pop("_non_bent_set")
        # Keep the canonical output compact.
        rec.pop("component_details")

    result = {
        "schema": "projective-switch-permutation-screen-v1",
        "date": "2026-08-02",
        "scope": "Eight projectively distinct APN switchings of Stage 5J",
        "necessary_condition": {
            "statement": "If F is CCZ-equivalent to a permutation, NB(F) union {0} contains two 4-dimensional subspaces V,W with V direct-sum W=GF(2)^8.",
            "status": "necessary only",
            "subspaces_enumerated": total_subspaces,
        },
        "records": records,
        "summary": {
            "functions": len(records),
            "direct_permutations": sum(rec["direct_permutation"] for rec in records),
            "pass_complementary_subspace_test": sum(
                rec["permutation_necessary_condition"]["passes_two_complementary_subspaces_test"]
                for rec in records
            ),
            "ruled_out_by_test": sum(
                not rec["permutation_necessary_condition"]["passes_two_complementary_subspaces_test"]
                for rec in records
            ),
        },
        "provenance": {
            "stage5j_file": str(stage5j_path.relative_to(root)),
            "stage5j_sha256": sha256(stage5j_path),
        },
        "validation": {
            "all_eight_APN_directly_verified": all(rec["APN_directly_verified"] for rec in records),
            "all_200787_four_spaces_enumerated": total_subspaces == 200787,
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "sha256": sha256(output_path),
        "functions": len(records),
        "direct_permutations": result["summary"]["direct_permutations"],
        "pass_test": result["summary"]["pass_complementary_subspace_test"],
        "ruled_out": result["summary"]["ruled_out_by_test"],
        "contained_4_spaces": [
            rec["permutation_necessary_condition"]["contained_4_spaces_count"]
            for rec in records
        ],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
