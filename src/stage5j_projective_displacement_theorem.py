#!/usr/bin/env python3
r"""Stage 5J: classify all projective GF(16)/GF(4) trace displacements.

For rho in GF(16)\GF(4), consider the rank-two family

    x^3 + eta*Q_1 + zeta*rho^(-3)*Q_rho,
    eta,zeta in GF(4)^*.

The coefficient normalization makes both rank-one margins individually APN.
The theorem proved from Stage 5H, GF(4)^*-scaling invariance, and Frobenius
transport is:

    APN iff eta=zeta and Tr_GF4/F2(delta(rho)*eta)=1,

where rho0=rho^6 is the unique order-five representative of rho*GF(4)^* and

    delta(rho)=(rho0+rho0^4)^2 in GF(4)^*.

An exhaustive 108-case derivative-subspace check is included only as an
independent certificate.
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
    polar_trace_switch,
    subspace_product,
    subspaces_2d,
    trace_f4_to_f2,
)

N = 256
R = 0xB0
LAMBDA = 0xBD


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage5h",
        type=Path,
        default=Path("results/current/rank_two_selector_proof.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/projective_displacement_theorem.json"),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    resolve = lambda p: p if p.is_absolute() else root / p
    stage5h_path = resolve(args.stage5h)
    output_path = resolve(args.output)
    stage5h = json.loads(stage5h_path.read_text(encoding="utf-8"))

    r = R
    lam = LAMBDA
    lam2 = gf_pow(lam, 2)
    f4 = sorted(x for x in range(N) if gf_pow(x, 4) == x)
    f4_star = [x for x in f4 if x]
    f16 = sorted(x for x in range(N) if gf_pow(x, 16) == x)
    outside = [x for x in f16 if x not in f4]
    assert len(outside) == 12
    assert set(f4_star) == {1, lam, lam2}

    # Recover the complete base-r target multiplicities from the compact
    # 16-pair table of Stage 5H.  The normalized combined target is
    # eta*r^3+zeta, and every occurrence contributes for each of the three
    # semilinear involutions.
    product_multiset: Counter[int] = Counter()
    for row in stage5h["semilinear_transition_proof"][
        "normalized_product_multiset"
    ]:
        product_multiset[int(row["hex"], 16)] = row["multiplicity"]
    base_pair_table = []
    base_counts: dict[tuple[int, int], int] = {}
    for eta in f4_star:
        for zeta in f4_star:
            target = gf_mult(eta, gf_pow(r, 3)) ^ zeta
            per_involution = product_multiset[target]
            total = 3 * per_involution
            base_counts[(eta, zeta)] = total
            base_pair_table.append({
                "eta_hex": f"0x{eta:02x}",
                "zeta_hex": f"0x{zeta:02x}",
                "normalized_target_hex": f"0x{target:02x}",
                "witnesses_per_involution": per_involution,
                "total_bad_sections": total,
                "APN": total == 0,
            })
    assert Counter(base_counts.values()) == {0: 2, 6: 6, 12: 1}
    assert {
        pair for pair, count in base_counts.items() if count == 0
    } == {(1, 1), (lam, lam)}

    # Projective classes: rho0=rho^6 removes the GF(4)^* component and is the
    # unique order-five representative of rho*GF(4)^*.
    classes: dict[int, list[int]] = {}
    for rho in outside:
        rho0 = gf_pow(rho, 6)
        assert gf_order(rho0) == 5
        classes.setdefault(rho0, []).append(rho)
    assert len(classes) == 4
    assert all(len(members) == 3 for members in classes.values())
    assert set(classes) == {gf_pow(r, j) for j in range(1, 5)}

    spaces = subspaces_2d()
    space_products = [subspace_product(space) for space in spaces]
    base_polar = [
        polar_trace_switch(1, space[1], space[2]) for space in spaces
    ]
    rho_polars = {
        rho: [polar_trace_switch(rho, space[1], space[2]) for space in spaces]
        for rho in outside
    }
    projective_records = []
    all_case_records = []
    for rho0, members in sorted(classes.items()):
        d = rho0 ^ gf_pow(rho0, 4)
        delta = gf_pow(d, 2)
        assert d in f4_star and delta in f4_star
        canonical_pair_records = []
        for eta in f4_star:
            for zeta in f4_star:
                predicted_apn = (
                    eta == zeta
                    and trace_f4_to_f2(gf_mult(delta, eta)) == 1
                )
                member_counts = []
                for rho in sorted(members):
                    second_coefficient = gf_mult(
                        zeta, gf_pow(gf_inv(rho), 3)
                    )
                    bad_count = 0
                    for b1, brho, product in zip(
                        base_polar, rho_polars[rho], space_products
                    ):
                        update = (eta if b1 else 0) ^ (
                            second_coefficient if brho else 0
                        )
                        if update and product == update:
                            bad_count += 1
                    member_counts.append(bad_count)
                    all_case_records.append({
                        "rho_hex": f"0x{rho:02x}",
                        "rho0_hex": f"0x{rho0:02x}",
                        "eta_hex": f"0x{eta:02x}",
                        "zeta_hex": f"0x{zeta:02x}",
                        "second_coefficient_hex": f"0x{second_coefficient:02x}",
                        "bad_section_count": bad_count,
                        "APN": bad_count == 0,
                        "predicted_APN": predicted_apn,
                    })
                    assert (bad_count == 0) == predicted_apn
                assert len(set(member_counts)) == 1
                bad_count = member_counts[0]
                # The four order-five representatives are the Frobenius orbit
                # of r, so counts must be the transported base-r counts after
                # the corresponding squaring of eta,zeta.
                k = next(
                    k for k in range(4) if gf_pow(r, 1 << k) == rho0
                )
                inverse_frobenius_power = 1 << ((4 - k) % 4)
                eta_base = gf_pow(eta, inverse_frobenius_power)
                zeta_base = gf_pow(zeta, inverse_frobenius_power)
                assert bad_count == base_counts[(eta_base, zeta_base)]
                canonical_pair_records.append({
                    "eta_hex": f"0x{eta:02x}",
                    "zeta_hex": f"0x{zeta:02x}",
                    "bad_section_count": bad_count,
                    "APN": bad_count == 0,
                    "selector_value": (
                        trace_f4_to_f2(gf_mult(delta, eta))
                        if eta == zeta else None
                    ),
                })
        accepted = [
            row for row in canonical_pair_records if row["APN"]
        ]
        assert len(accepted) == 2
        assert Counter(
            row["bad_section_count"] for row in canonical_pair_records
        ) == {0: 2, 6: 6, 12: 1}
        projective_records.append({
            "canonical_order5_rho_hex": f"0x{rho0:02x}",
            "projective_members_hex": [f"0x{x:02x}" for x in sorted(members)],
            "d_equals_rho0_plus_rho0_fourth_hex": f"0x{d:02x}",
            "delta_equals_d_squared_hex": f"0x{delta:02x}",
            "pair_table": canonical_pair_records,
            "accepted_pairs": accepted,
        })

    assert len(all_case_records) == 108
    assert all(row["APN"] == row["predicted_APN"] for row in all_case_records)

    # Frobenius orbits of the eight accepted projective pairs.
    accepted_pairs = {
        (int(row["canonical_order5_rho_hex"], 16),
         int(pair["eta_hex"], 16))
        for row in projective_records
        for pair in row["accepted_pairs"]
    }
    assert len(accepted_pairs) == 8
    remaining = set(accepted_pairs)
    frobenius_orbits = []
    while remaining:
        start = min(remaining)
        orbit = []
        current = start
        while current not in orbit:
            orbit.append(current)
            current = (gf_pow(current[0], 2), gf_pow(current[1], 2))
        assert current == start
        assert set(orbit) <= accepted_pairs
        frobenius_orbits.append(orbit)
        remaining -= set(orbit)
    assert len(frobenius_orbits) == 2
    assert all(len(orbit) == 4 for orbit in frobenius_orbits)

    result = {
        "schema": "projective-displacement-theorem-v1",
        "date": "2026-08-02",
        "theorem": {
            "family": "x^3+eta*Q_1+zeta*rho^(-3)*Q_rho, rho in GF(16)\\GF(4), eta,zeta in GF(4)^*",
            "canonical_displacement": "rho0=rho^6, the unique order-5 representative of rho*GF(4)^*",
            "selector": "delta(rho)=(rho0+rho0^4)^2 in GF(4)^*",
            "APN_condition": "eta=zeta and Tr_GF4/F2(delta(rho)*eta)=1",
            "bad_count_rule": "0 for the two selected diagonal pairs; 12 for the rejected diagonal pair; 6 for every off-diagonal pair",
        },
        "proof_mechanism": {
            "base_case": "The complete nine-pair table for rho=r is read from the Stage 5H four-row product table.",
            "GF4_scaling": "Replacing rho by alpha*rho, alpha in GF(4)^*, leaves the polar form B_rho and rho^(-3) unchanged.",
            "Frobenius_transport": "The four projective displacements are the orbit r,r^2,r^4,r^3 under squaring; x^3 and APN are preserved under field Frobenius conjugation.",
            "exhaustive_check_role": "The 108-case subspace scan is an independent certificate, not the source of the theorem.",
        },
        "base_r_nine_pair_table": base_pair_table,
        "projective_classes": projective_records,
        "accepted_Frobenius_orbits": [
            [
                {
                    "rho0_hex": f"0x{rho0:02x}",
                    "eta_equals_zeta_hex": f"0x{eta:02x}",
                }
                for rho0, eta in orbit
            ]
            for orbit in frobenius_orbits
        ],
        "exhaustive_certificate": {
            "projective_classes": 4,
            "representatives_per_class": 3,
            "coefficient_pairs_per_representative": 9,
            "total_cases": len(all_case_records),
            "all_predictions_verified": True,
            "case_records": all_case_records,
        },
        "interpretation": [
            "The order-five displacement found in the campaign is not a single accidental pair: every nontrivial projective GF(16)/GF(4) displacement supports exactly two APN synchronized lifts.",
            "The selector is intrinsic after removing the GF(4)^* component by rho0=rho^6.",
            "Off-diagonal choices of the two GF(4)^* fibre coefficients always fail, with exactly six bad derivative sections.",
            "The eight projective APN lifts form two Frobenius orbits of length four.",
        ],
        "provenance": {
            "stage5h_file": str(stage5h_path.relative_to(root)),
            "stage5h_sha256": hashlib.sha256(stage5h_path.read_bytes()).hexdigest(),
        },
        "validation": {
            "base_nine_pair_table_from_compact_products": True,
            "four_projective_classes": True,
            "GF4_scaling_invariance_checked": True,
            "Frobenius_transport_counts_checked": True,
            "all_108_subspace_cases_match_theorem": True,
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
        "projective_classes": len(projective_records),
        "APN_projective_pairs": len(accepted_pairs),
        "Frobenius_orbits": len(frobenius_orbits),
        "exhaustive_cases": len(all_case_records),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
