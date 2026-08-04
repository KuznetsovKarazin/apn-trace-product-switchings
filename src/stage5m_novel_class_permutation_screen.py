#!/usr/bin/env python3
"""Stage 5M: CCZ-permutation necessary-condition screen for eight archived classes."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from stage5l_permutation_screen import (
    N,
    binary_rank,
    fwht,
    is_apn_quadratic,
    parity,
    rref_subspaces_4_in_8,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyse_sbox(identifier: str, lut: list[int]) -> dict:
    preimages = Counter(lut)
    profile = Counter(preimages.values())
    profile[0] = N - len(preimages)
    amp_dist = Counter()
    non_bent = set()
    for b in range(1, N):
        signs = np.fromiter(
            (1 if parity(b & lut[x]) == 0 else -1 for x in range(N)),
            dtype=np.int64,
            count=N,
        )
        spec = fwht(signs)
        nonzero = spec[spec != 0]
        amps = set(abs(int(v)) for v in nonzero)
        if len(amps) != 1:
            raise AssertionError((identifier, b, sorted(amps)))
        amp = next(iter(amps))
        amp_dist[amp] += 1
        if len(nonzero) != N:
            non_bent.add(b)
    return {
        "id": identifier,
        "APN_directly_verified": is_apn_quadratic(lut),
        "direct_permutation": len(preimages) == N,
        "image_size": len(preimages),
        "preimage_multiplicity_distribution": {str(k): v for k, v in sorted(profile.items())},
        "component_amplitude_distribution": {str(k): v for k, v in sorted(amp_dist.items())},
        "bent_component_count": 255 - len(non_bent),
        "non_bent_component_count": len(non_bent),
        "non_bent_components_hex": [f"0x{x:02x}" for x in sorted(non_bent)],
        "_nb": non_bent,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--representatives",
        type=Path,
        default=Path("results/evidence/all_8_novel_class_representatives.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/novel_class_permutation_screen.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    resolve = lambda p: p if p.is_absolute() else root / p
    reps_path = resolve(args.representatives)
    output_path = resolve(args.output)
    source = json.loads(reps_path.read_text(encoding="utf-8"))
    records = [analyse_sbox(row["id"], row["sbox"]) for row in source["functions"]]
    assert len(records) == 8

    nb_sets = [r["_nb"] for r in records]
    spaces = [[] for _ in records]
    total = 0
    for basis, vectors in rref_subspaces_4_in_8():
        total += 1
        s = set(vectors)
        for i, nb in enumerate(nb_sets):
            if s <= nb:
                spaces[i].append(tuple(basis))
    assert total == 200787

    for rec, contained in zip(records, spaces):
        pairs = []
        for i in range(len(contained)):
            for j in range(i + 1, len(contained)):
                if binary_rank(list(contained[i]) + list(contained[j])) == 8:
                    pairs.append((contained[i], contained[j]))
        rec["permutation_necessary_condition"] = {
            "contained_4_spaces_count": len(contained),
            "contained_4_space_bases_hex": [[f"0x{x:02x}" for x in b] for b in contained],
            "complementary_pairs_count": len(pairs),
            "complementary_pairs": [
                {
                    "first_basis_hex": [f"0x{x:02x}" for x in a],
                    "second_basis_hex": [f"0x{x:02x}" for x in b],
                }
                for a, b in pairs
            ],
            "passes": bool(pairs),
            "meaning": "Necessary, not sufficient, for CCZ-equivalence to a permutation.",
        }
        rec.pop("_nb")

    result = {
        "schema": "novel-class-permutation-screen-v1",
        "date": "2026-08-02",
        "scope": "Eight archived class representatives; no article files modified",
        "condition": "NB(F) union {0} must contain two complementary 4-spaces if F is CCZ-equivalent to an 8-bit permutation.",
        "records": records,
        "summary": {
            "classes": 8,
            "APN_verified": sum(r["APN_directly_verified"] for r in records),
            "direct_permutations": sum(r["direct_permutation"] for r in records),
            "pass_necessary_condition": sum(r["permutation_necessary_condition"]["passes"] for r in records),
            "ruled_out": sum(not r["permutation_necessary_condition"]["passes"] for r in records),
        },
        "provenance": {
            "representatives_file": str(reps_path.relative_to(root)),
            "representatives_sha256": sha256(reps_path),
        },
        "validation": {
            "all_200787_four_spaces_enumerated": total == 200787,
            "all_APN_verified": all(r["APN_directly_verified"] for r in records),
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "sha256": sha256(output_path),
        "summary": result["summary"],
        "rows": [
            {
                "id": r["id"],
                "image": r["image_size"],
                "amp": r["component_amplitude_distribution"],
                "NB": r["non_bent_component_count"],
                "spaces4": r["permutation_necessary_condition"]["contained_4_spaces_count"],
                "pairs": r["permutation_necessary_condition"]["complementary_pairs_count"],
            }
            for r in records
        ],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
