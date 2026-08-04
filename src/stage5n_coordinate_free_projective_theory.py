#!/usr/bin/env python3
r"""Stage 5N: coordinate-free packaging of the n=8 projective trace-switch theorem.

The byte-level theorem of Stage 5J is recast using the field tower

    F2 < E=F4 < L=F16 < K=F256.

For c in L, let ell_c(x)=Tr_{K/E}(c x) and let q_[c] be the quadratic
Boolean form defined, modulo affine-linear terms, by

    q_[c](x)=Tr_{E/F2}(ell_c(x)) Tr_{E/F2}(lambda ell_c(x)),

where lambda in E\F2.  The class q_[c] depends only on the projective point
[c] in P(L/E)=P^1(E).  Fixing p_infinity=[1], every other point p=[rho]
has a canonical order-five representative s(p)=rho^6.  The selector

    delta(p)=Tr_{L/E}(s(p))^2 in E^*

is representative-independent.  Stage 5J becomes a theorem on marked pairs
of projective points and normalized E^*-coefficients.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from stage5g_trace_coordinate_theory import gf_mult, gf_pow, gf_order, trace_f4_to_f2


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
        default=Path("results/current/coordinate_free_projective_theory.json"),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    stage5j_path = args.stage5j if args.stage5j.is_absolute() else root / args.stage5j
    output_path = args.output if args.output.is_absolute() else root / args.output
    stage5j = json.loads(stage5j_path.read_text(encoding="utf-8"))

    E = sorted(x for x in range(256) if gf_pow(x, 4) == x)
    E_star = [x for x in E if x]
    L = sorted(x for x in range(256) if gf_pow(x, 16) == x)
    outside = [x for x in L if x not in E]
    assert len(E) == 4 and len(L) == 16 and len(outside) == 12

    # Projective points of P(L/E): E^*-orbits on L^*.
    remaining = set(x for x in L if x)
    projective_classes = []
    while remaining:
        representative = min(remaining)
        members = sorted(gf_mult(representative, alpha) for alpha in E_star)
        projective_classes.append(members)
        remaining -= set(members)
    assert len(projective_classes) == 5

    base_class = next(c for c in projective_classes if 1 in c)
    nonbase = [c for c in projective_classes if c != base_class]
    assert len(nonbase) == 4

    stage5j_records = {
        int(row["canonical_order5_rho_hex"], 16): row
        for row in stage5j["projective_classes"]
    }

    records = []
    accepted_marked_points = []
    for members in sorted(nonbase, key=lambda xs: min(xs)):
        # rho^6 is invariant under rho -> alpha*rho, alpha in E^*, since alpha^3=1.
        s_values = {gf_pow(rho, 6) for rho in members}
        assert len(s_values) == 1
        s = next(iter(s_values))
        assert gf_order(s) == 5

        # Trace L/E for the quadratic subextension is z+z^4.
        t = s ^ gf_pow(s, 4)
        assert t in E_star
        delta = gf_pow(t, 2)
        assert delta in E_star

        accepted_eta = [eta for eta in E_star if trace_f4_to_f2(gf_mult(delta, eta)) == 1]
        rejected_eta = [eta for eta in E_star if eta not in accepted_eta]
        assert len(accepted_eta) == 2 and len(rejected_eta) == 1

        source = stage5j_records[s]
        stage5j_accepted = sorted(int(row["eta_hex"], 16) for row in source["accepted_pairs"])
        assert sorted(accepted_eta) == stage5j_accepted

        # rho^{-3} is also representative-independent because alpha^{-3}=1.
        norm_values = {gf_pow(rho, 252) for rho in members}  # rho^{-3}=rho^(255-3)
        assert len(norm_values) == 1

        record = {
            "projective_members_hex": [f"0x{x:02x}" for x in members],
            "canonical_mu5_representative_hex": f"0x{s:02x}",
            "relative_trace_L_over_E_hex": f"0x{t:02x}",
            "selector_delta_hex": f"0x{delta:02x}",
            "representative_independent_rho_inverse_cube_hex": f"0x{next(iter(norm_values)):02x}",
            "accepted_normalized_coefficients_hex": [f"0x{x:02x}" for x in accepted_eta],
            "rejected_normalized_coefficient_hex": f"0x{rejected_eta[0]:02x}",
        }
        records.append(record)
        accepted_marked_points.extend((s, eta) for eta in accepted_eta)

    # Frobenius orbits of the eight marked switchings.
    remaining_pairs = set(accepted_marked_points)
    frobenius_orbits = []
    while remaining_pairs:
        start = min(remaining_pairs)
        orbit = []
        current = start
        while current not in orbit:
            orbit.append(current)
            current = (gf_pow(current[0], 2), gf_pow(current[1], 2))
        assert current == start
        frobenius_orbits.append(orbit)
        remaining_pairs -= set(orbit)
    assert sorted(len(o) for o in frobenius_orbits) == [4, 4]

    result = {
        "schema": "coordinate-free-projective-trace-switch-theory-v1",
        "date": "2026-08-02",
        "field_tower": "GF(2) subset E=GF(4) subset L=GF(16) subset K=GF(256)",
        "objects": {
            "relative_trace_functional": "ell_c(x)=Tr_K/E(c*x)",
            "projective_quadratic_form": "q_[c](x)=Tr_E/F2(ell_c(x))*Tr_E/F2(lambda*ell_c(x)), modulo affine-linear Boolean forms",
            "projective_parameter_space": "P(L/E)=P^1(E), with 5 points",
            "base_point": "p_infinity=[1]",
            "canonical_displacement": "s(p)=rho^6 in mu_5 for p=[rho] != p_infinity",
            "selector": "delta(p)=Tr_L/E(s(p))^2 in E^*",
        },
        "coordinate_free_theorem": {
            "family": "F_{p,eta,zeta}=x^3+eta*q_[1]+zeta*rho^(-3)*q_[rho], eta,zeta in E^*",
            "well_definedness": "Changing rho to alpha*rho with alpha in E^* changes neither rho^(-3), s(p), delta(p), nor the quadratic form modulo affine-linear terms.",
            "APN_criterion": "F_{p,eta,zeta} is APN iff eta=zeta and Tr_E/F2(delta(p)*eta)=1.",
            "local_count": "Each of the 4 non-base projective points supports exactly 2 APN normalized coefficients.",
            "global_count": "There are 8 marked rank-two switchings, split into 2 Frobenius orbits of length 4.",
        },
        "projective_records": records,
        "accepted_Frobenius_orbits": [
            [
                {"s_hex": f"0x{s:02x}", "eta_hex": f"0x{eta:02x}"}
                for s, eta in orbit
            ]
            for orbit in frobenius_orbits
        ],
        "proof_dependencies": [
            "general low-rank derivative-incidence criterion",
            "relative-trace normal form for q_[c]",
            "rank-one coefficient theorem in n=8",
            "rank-two involution selector proof",
            "GF(4)^*-scaling invariance and Frobenius transport",
        ],
        "validation": {
            "all_5_projective_points_constructed": len(projective_classes) == 5,
            "all_4_nonbase_points_constructed": len(records) == 4,
            "representative_independence_checked_for_all_12_nonbase_representatives": True,
            "selector_reproduces_all_stage5j_accepted_pairs": True,
            "eight_marked_switchings_recovered": len(accepted_marked_points) == 8,
            "two_Frobenius_orbits_recovered": len(frobenius_orbits) == 2,
        },
        "provenance": {
            "stage5j_file": str(stage5j_path.relative_to(root)),
            "stage5j_sha256": sha256(stage5j_path),
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "sha256": sha256(output_path),
        "projective_points": len(projective_classes),
        "marked_switchings": len(accepted_marked_points),
        "frobenius_orbits": [len(o) for o in frobenius_orbits],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
