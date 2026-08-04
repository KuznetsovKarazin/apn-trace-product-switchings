#!/usr/bin/env python3
"""Stage 5B/5C: coefficient-support convention and exact slice stabilizer orbits.

A three-pair Gröbner slice does not select a subspace of monomial vectors.
It selects the 24-dimensional coefficient variation space

    E_T = V_out tensor U_T^*  <= Hom(Lambda^2 V, V_out),

where U_T^* is spanned by the three coordinate covectors belonging to the
selected input pairs.  Under the normalized centralizer action

    C |-> B C Lambda^2(B^{-1}),

the output factor remains the whole V_out, while the support covectors act by

    r |-> r Lambda^2(B^{-1}),

or, in column convention, Lambda^2(B^{-1})^T.

This script:
  * freezes and validates that contragredient convention;
  * recomputes the GF(16) module fingerprints of all 3276 coefficient-support
    slices;
  * enumerates their full GL_2(GF(16)) module orbits;
  * computes exact full slice orbits under each cached centre stabilizer;
  * checks invariance of exact solution counts on every coordinate-coordinate
    intersection of a stabilizer orbit.

It does not launch any APN search and does not touch manuscript files.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any, Iterable

GROUP_ORDER = 61200
GL2_GENERATORS = [
    (0x3, 0x0, 0x0, 0x1),
    (0x0, 0x1, 0x1, 0x0),
    (0x1, 0x1, 0x0, 0x1),
]
EXPECTED_CATEGORY_DISTRIBUTION = {
    "all_11_gold": 201,
    "no_gold": 2972,
    "x9_only_all": 100,
    "productive_10_of_11": 3,
}


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gf16_inverse_2x2(stage3, entry):
    a, b, c, d = entry
    determinant = stage3.f16_mul(a, d) ^ stage3.f16_mul(b, c)
    if determinant == 0:
        raise RuntimeError("singular GF(16) matrix")
    scale = stage3.f16_inv(determinant)
    # Minus equals plus in characteristic two.
    return (
        stage3.f16_mul(scale, d),
        stage3.f16_mul(scale, b),
        stage3.f16_mul(scale, c),
        stage3.f16_mul(scale, a),
    )


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def canonical_space16(stage3, rows):
    return tuple(
        tuple(row)
        for row in stage3.rref16([list(row) for row in rows])[0]
        if any(row)
    )


def apply_matrix16_to_space(stage3, space, matrix):
    transformed = []
    for vector in space:
        image = []
        for i in range(len(matrix)):
            value = 0
            for j in range(len(vector)):
                if matrix[i][j] and vector[j]:
                    value ^= stage3.f16_mul(matrix[i][j], vector[j])
            image.append(value)
        transformed.append(image)
    return canonical_space16(stage3, transformed)


def orbit16(stage3, start, generators):
    seen = {start}
    queue = deque([start])
    while queue:
        current = queue.popleft()
        for generator in generators:
            image = apply_matrix16_to_space(stage3, current, generator)
            if image not in seen:
                seen.add(image)
                queue.append(image)
    return seen


def rep_hex(representative) -> str:
    return "_".join("".join(f"{value:x}" for value in row) for row in representative)


def module_orbit_label(rank: int, size: int, representative) -> str:
    return f"CS-R{rank}-O{size}-{rep_hex(representative)}"


def bitrow_from_binary_row(row: Iterable[int]) -> int:
    return sum((value & 1) << index for index, value in enumerate(row))


def binary_rref(rows: Iterable[int], width: int) -> tuple[int, ...]:
    work = list(rows)
    rank = 0
    for column in range(width):
        pivot = next(
            (i for i in range(rank, len(work)) if (work[i] >> column) & 1),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for i in range(len(work)):
            if i != rank and ((work[i] >> column) & 1):
                work[i] ^= work[rank]
        rank += 1
        if rank == len(work):
            break
    return tuple(work[:rank])


def row_times_binary_matrix(row_bits: int, matrix_rows_bits: list[int]) -> int:
    output = 0
    value = row_bits
    while value:
        low = value & -value
        output ^= matrix_rows_bits[low.bit_length() - 1]
        value ^= low
    return output


def transform_binary_rowspace(
    space: tuple[int, ...],
    matrix_rows_bits: list[int],
    width: int,
) -> tuple[int, ...]:
    return binary_rref(
        [row_times_binary_matrix(row, matrix_rows_bits) for row in space],
        width,
    )


def parse_gf16_matrix_json(value) -> tuple[int, int, int, int]:
    return tuple(int(item, 16) for row in value for item in row)


def matrix_entry_json(entry):
    a, b, c, d = entry
    return [[f"0x{a:x}", f"0x{b:x}"], [f"0x{c:x}", f"0x{d:x}"]]


def solution_count_for_seed(row: dict[str, Any], seed: int) -> int:
    if str(seed) in row["x3_solution_counts"]:
        return int(row["x3_solution_counts"][str(seed)])
    if str(seed) in row["x9_solution_counts"]:
        return int(row["x9_solution_counts"][str(seed)])
    raise KeyError(seed)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage3",
        type=Path,
        default=Path("code/current/analyse_c5_module_fingerprints.py"),
    )
    parser.add_argument(
        "--stage4",
        type=Path,
        default=Path("code/current/analyse_centralizer_module_orbits.py"),
    )
    parser.add_argument(
        "--stage5a-script",
        type=Path,
        default=Path("code/current/stage5_center_action.py"),
    )
    parser.add_argument(
        "--stage5a",
        type=Path,
        default=Path("results/current/stage5_center_action.json"),
    )
    parser.add_argument(
        "--correlation",
        type=Path,
        default=Path("results/current/module_productivity_correlation_canonical.json"),
    )
    parser.add_argument(
        "--module-output",
        type=Path,
        default=Path("results/current/coefficient_support_module_orbits_canonical.json"),
    )
    parser.add_argument(
        "--slice-output",
        type=Path,
        default=Path("results/current/stage5_slice_stabilizer_orbits.json"),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]

    def resolve(path: Path) -> Path:
        return path if path.is_absolute() else root / path

    stage3_path = resolve(args.stage3)
    stage4_path = resolve(args.stage4)
    stage5a_script_path = resolve(args.stage5a_script)
    stage5a_path = resolve(args.stage5a)
    correlation_path = resolve(args.correlation)
    module_output_path = resolve(args.module_output)
    slice_output_path = resolve(args.slice_output)

    stage3 = load_module(stage3_path, "stage3_slice_geometry")
    stage4 = load_module(stage4_path, "stage4_slice_geometry")
    stage5a_code = load_module(stage5a_script_path, "stage5a_slice_geometry")

    correlation = json.loads(correlation_path.read_text(encoding="utf-8"))
    stage5a = json.loads(stage5a_path.read_text(encoding="utf-8"))
    rows = correlation["campaign"]["per_triple"]
    if correlation["campaign"]["category_distribution"] != EXPECTED_CATEGORY_DISTRIBUTION:
        raise RuntimeError("canonical Stage 3 category distribution mismatch")
    if len(rows) != 3276:
        raise RuntimeError(f"expected 3276 triples, got {len(rows)}")

    decomposition, diagonalising, inverse, offsets = (
        stage4.reconstruct_diagonalisation(stage3)
    )
    omega_start, omega_end = offsets[stage3.OMEGA]

    # Vector fingerprints from Stage 3 use D^{-1} e_i.  The true coefficient
    # support covectors use e_i^* D.  The latter are rows of the eigenbasis
    # matrix D and are not generally equal to the former.
    coefficient_support_by_pair = {
        pair: diagonalising[pair_index][omega_start:omega_end]
        for pair_index, pair in enumerate(stage3.PAIRS)
    }

    dual_generators = []
    dual_generator_checks = []
    for entry in GL2_GENERATORS:
        inverse_entry = gf16_inverse_2x2(stage3, entry)
        representation_inverse = stage4.representation5(
            stage3,
            diagonalising,
            inverse,
            offsets,
            inverse_entry,
        )
        dual_matrix = transpose(representation_inverse)
        dual_generators.append(dual_matrix)

        binary = stage4.gl2_to_binary(stage3, *entry)
        binary_inverse = stage5a_code.inverse2(binary)
        wedge_inverse = stage3.wedge_action(binary_inverse)

        # Validate on all 28 coordinate covectors: transform in 28D first,
        # then project to the omega block, and compare with the 5D dual matrix.
        all_pairs_agree = True
        for pair_index, pair in enumerate(stage3.PAIRS):
            row = [0] * len(stage3.PAIRS)
            row[pair_index] = 1
            transformed_row = [
                sum(
                    row[k] & wedge_inverse[k][j]
                    for k in range(len(stage3.PAIRS))
                ) & 1
                for j in range(len(stage3.PAIRS))
            ]
            direct_coordinates = []
            for column in range(omega_start, omega_end):
                value = 0
                for k in range(len(stage3.PAIRS)):
                    if transformed_row[k]:
                        value ^= diagonalising[k][column]
                direct_coordinates.append(value)
            source = coefficient_support_by_pair[pair]
            matrix_coordinates = []
            for i in range(5):
                value = 0
                for j in range(5):
                    if dual_matrix[i][j] and source[j]:
                        value ^= stage3.f16_mul(dual_matrix[i][j], source[j])
                matrix_coordinates.append(value)
            if direct_coordinates != matrix_coordinates:
                all_pairs_agree = False
                break
        if not all_pairs_agree:
            raise RuntimeError("28D covector action and 5D dual action disagree")
        dual_generator_checks.append({
            "GL2_GF16_entry": matrix_entry_json(entry),
            "dual_formula": "transpose(rep_omega(B^{-1}))",
            "all_28_coordinate_covectors_agree": True,
        })

    dual_generators_with_inverses = dual_generators + [
        stage3.inverse16(matrix) for matrix in dual_generators
    ]

    triple_spaces = {}
    observed_spaces = set()
    for row in rows:
        space = canonical_space16(
            stage3,
            [coefficient_support_by_pair[tuple(pair)] for pair in row["pairs"]],
        )
        triple_spaces[row["triple_index"]] = space
        observed_spaces.add(space)

    unassigned = set(observed_spaces)
    unsorted_orbits = []
    while unassigned:
        start_space = min(unassigned)
        full_orbit = orbit16(
            stage3,
            start_space,
            dual_generators_with_inverses,
        )
        representative = min(full_orbit)
        observed = observed_spaces.intersection(full_orbit)
        unsorted_orbits.append({
            "module_rank": len(representative),
            "full_GL2_GF16_orbit_size": len(full_orbit),
            "representative": representative,
            "observed_spaces": observed,
        })
        unassigned.difference_update(observed)

    unsorted_orbits.sort(key=lambda record: (
        record["module_rank"],
        record["full_GL2_GF16_orbit_size"],
        record["representative"],
    ))
    space_to_label = {}
    module_orbits = []
    for record in unsorted_orbits:
        label = module_orbit_label(
            record["module_rank"],
            record["full_GL2_GF16_orbit_size"],
            record["representative"],
        )
        for space in record["observed_spaces"]:
            space_to_label[space] = label
        members = [
            row for row in rows
            if triple_spaces[row["triple_index"]] in record["observed_spaces"]
        ]
        module_orbits.append({
            "orbit_label": label,
            "module_rank": record["module_rank"],
            "full_GL2_GF16_orbit_size": record["full_GL2_GF16_orbit_size"],
            "canonical_representative_rref": [
                list(row) for row in record["representative"]
            ],
            "observed_unique_spaces": len(record["observed_spaces"]),
            "observed_triples": len(members),
            "category_distribution": dict(Counter(
                row["productivity_category"] for row in members
            )),
            "graph_type_distribution": dict(Counter(
                row["graph_type"] for row in members
            )),
            "triple_indices": [row["triple_index"] for row in members],
        })

    for row in rows:
        row["coefficient_support_module_orbit_label"] = space_to_label[
            triple_spaces[row["triple_index"]]
        ]
        row["coefficient_support_module_rank"] = len(
            triple_spaces[row["triple_index"]]
        )

    rank_distribution = Counter(
        row["coefficient_support_module_rank"] for row in rows
    )
    rank_by_category = defaultdict(Counter)
    for row in rows:
        rank_by_category[row["productivity_category"]][
            row["coefficient_support_module_rank"]
        ] += 1

    productive_labels = [
        record["orbit_label"] for record in module_orbits
        if any(
            category != "no_gold" and count
            for category, count in record["category_distribution"].items()
        )
    ]
    pure_negative_labels = [
        record["orbit_label"] for record in module_orbits
        if set(record["category_distribution"]) == {"no_gold"}
    ]
    pure_negative_excluded = sum(
        record["observed_triples"] for record in module_orbits
        if record["orbit_label"] in pure_negative_labels
    )

    # Find the two mixed rank-two strata in the corrected coefficient-support
    # convention.
    mixed_rank_two_labels = [
        record["orbit_label"] for record in module_orbits
        if record["module_rank"] == 2
        and len(record["category_distribution"]) > 1
    ]
    if len(mixed_rank_two_labels) != 2:
        raise RuntimeError(
            f"expected two mixed rank-two orbits, got {mixed_rank_two_labels}"
        )

    module_report = {
        "schema": "coefficient-support-c5-module-orbits-v1-canonical",
        "date": "2026-08-02",
        "provenance": {
            "canonical_correlation_file": str(correlation_path.relative_to(root)),
            "canonical_correlation_sha256": sha256_file(correlation_path),
            "stage3_script": str(stage3_path.relative_to(root)),
            "stage3_script_sha256": sha256_file(stage3_path),
            "stage4_script": str(stage4_path.relative_to(root)),
            "stage4_script_sha256": sha256_file(stage4_path),
        },
        "object_convention": {
            "ambient_quadratic_space": "Hom(Lambda^2 V, V_out)",
            "three_pair_slice_direction": "E_T = V_out tensor U_T^*",
            "U_T_dual_definition": (
                "span of the three coordinate coefficient covectors selected "
                "by the pair triple"
            ),
            "coefficient_action": "C maps to B C Lambda^2(B^{-1})",
            "support_covector_action_row": "r maps to r Lambda^2(B^{-1})",
            "support_covector_action_column": (
                "r^T maps to Lambda^2(B^{-1})^T r^T"
            ),
            "omega_fingerprint": (
                "row e_pair^* D restricted to the omega eigenspace columns; "
                "equivalently a row of the diagonalising eigenbasis matrix D"
            ),
            "warning": (
                "The earlier vector fingerprint D^{-1}e_pair classifies "
                "monomial vectors, not the coefficient support of the actual "
                "Groebner slice."
            ),
        },
        "generator_validation": dual_generator_checks,
        "group": {
            "name": "GL_2(GF(16))",
            "order": GROUP_ORDER,
            "dual_omega_action": "transpose(rep_omega(B^{-1}))",
        },
        "summary": {
            "triple_count": len(rows),
            "observed_unique_support_spaces": len(observed_spaces),
            "coefficient_support_orbit_count": len(module_orbits),
            "module_rank_distribution": {
                str(rank): count for rank, count in sorted(rank_distribution.items())
            },
            "module_rank_by_category": {
                category: {
                    str(rank): count for rank, count in sorted(counts.items())
                }
                for category, counts in sorted(rank_by_category.items())
            },
            "productive_orbit_labels": productive_labels,
            "pure_negative_orbit_labels": pure_negative_labels,
            "pure_negative_filter_excludes": pure_negative_excluded,
            "mixed_rank_two_orbit_labels": mixed_rank_two_labels,
        },
        "orbits": module_orbits,
        "main_conclusions": [
            "The actual coefficient-support convention yields 15 observed GL_2(GF(16)) orbit types, not eight vector-support types.",
            "All 40 coefficient-support rank-one triples are productive: 28 all-Gold and 12 x9-only.",
            "All 2452 coefficient-support rank-three triples are nonproductive in the completed campaign.",
            "Only two coefficient-support rank-two orbit types contain productive triples.",
            "The pure-negative module-orbit filter excludes 2552 of 3276 triples, improving the earlier vector-support filter by 192 triples.",
        ],
    }

    # Exact stabilizer action on full 3-dimensional coefficient-support
    # subspaces in the 28-dimensional dual monomial space.
    coordinate_spaces = {}
    triple_pairs_by_index = {}
    for triple_index, pair_indices in enumerate(
        itertools.combinations(range(len(stage3.PAIRS)), 3)
    ):
        space = binary_rref((1 << i for i in pair_indices), len(stage3.PAIRS))
        coordinate_spaces[space] = triple_index
        triple_pairs_by_index[triple_index] = pair_indices
    if len(coordinate_spaces) != 3276:
        raise RuntimeError("coordinate support space construction failed")

    kernel_entries = []
    for exponent in range(5):
        scalar = stage3.f16_pow(stage3.OMEGA, exponent)
        kernel_entries.append((scalar, 0, 0, scalar))

    seed98_record = next(
        record for record in stage5a["centre_orbits"]["centres"]
        if record["seed"] == 98
    )
    seed98_entries = [
        parse_gf16_matrix_json(value)
        for value in seed98_record["stabilizer"]["all_stabilizer_elements"]
    ]

    def support_action_matrix_bits(entry):
        binary = stage4.gl2_to_binary(stage3, *entry)
        binary_inverse = stage5a_code.inverse2(binary)
        wedge_inverse = stage3.wedge_action(binary_inverse)
        return [bitrow_from_binary_row(row) for row in wedge_inverse]

    def classify_coordinate_intersections(entries):
        matrices = [support_action_matrix_bits(entry) for entry in entries]
        triple_to_class = {}
        classes = []
        for start_space, start_index in coordinate_spaces.items():
            if start_index in triple_to_class:
                continue
            full_orbit = {
                transform_binary_rowspace(
                    start_space,
                    matrix,
                    len(stage3.PAIRS),
                )
                for matrix in matrices
            }
            coordinate_intersection = sorted(
                coordinate_spaces[space]
                for space in full_orbit
                if space in coordinate_spaces
            )
            class_index = len(classes)
            for triple_index in coordinate_intersection:
                triple_to_class[triple_index] = class_index
            classes.append({
                "full_stabilizer_orbit_size": len(full_orbit),
                "coordinate_intersection": coordinate_intersection,
                "coordinate_intersection_size": len(coordinate_intersection),
                "canonical_full_orbit_minimum_hex": [
                    f"0x{row:07x}" for row in min(full_orbit)
                ],
            })
        if len(triple_to_class) != 3276:
            raise RuntimeError("not all coordinate triples assigned")
        return classes, triple_to_class

    c5_classes, c5_lookup = classify_coordinate_intersections(kernel_entries)
    c10_classes, c10_lookup = classify_coordinate_intersections(seed98_entries)

    c5_partition = sorted(
        tuple(record["coordinate_intersection"]) for record in c5_classes
    )
    c10_partition = sorted(
        tuple(record["coordinate_intersection"]) for record in c10_classes
    )
    same_coordinate_partition = c5_partition == c10_partition
    if not same_coordinate_partition:
        raise RuntimeError("C5 and seed98 C10 coordinate partitions differ")

    gold_seeds = correlation["campaign"]["x3_seeds"] + correlation["campaign"]["x9_seeds"]
    solution_invariance = []
    for seed in gold_seeds:
        lookup = c10_lookup if seed == 98 else c5_lookup
        values_by_class = defaultdict(set)
        for row in rows:
            values_by_class[lookup[row["triple_index"]]].add(
                solution_count_for_seed(row, seed)
            )
        bad = {
            str(class_index): sorted(values)
            for class_index, values in values_by_class.items()
            if len(values) != 1
        }
        if bad:
            raise RuntimeError(
                f"solution count is not invariant under stabilizer for seed {seed}: {bad}"
            )
        solution_invariance.append({
            "seed": seed,
            "stabilizer": "C10" if seed == 98 else "C5",
            "coordinate_classes_tested": len(values_by_class),
            "exact_solution_count_invariant": True,
        })

    category_by_index = {
        row["triple_index"]: row["productivity_category"] for row in rows
    }
    category_purity_failures = []
    for class_index, record in enumerate(c5_classes):
        categories = {
            category_by_index[index]
            for index in record["coordinate_intersection"]
        }
        if len(categories) != 1:
            category_purity_failures.append({
                "class_index": class_index,
                "triple_indices": record["coordinate_intersection"],
                "categories": sorted(categories),
            })
    if category_purity_failures:
        raise RuntimeError("stabilizer coordinate classes are not category-pure")

    orbit_label_by_triple = {
        row["triple_index"]: row["coefficient_support_module_orbit_label"]
        for row in rows
    }

    mixed_strata = []
    for label in mixed_rank_two_labels:
        triple_indices = sorted(
            index for index, value in orbit_label_by_triple.items()
            if value == label
        )
        class_indices = sorted({c5_lookup[index] for index in triple_indices})
        class_size_distribution = Counter(
            c5_classes[class_index]["coordinate_intersection_size"]
            for class_index in class_indices
        )
        triple_category_distribution = Counter(
            category_by_index[index] for index in triple_indices
        )
        class_category_distribution = Counter(
            category_by_index[
                c5_classes[class_index]["coordinate_intersection"][0]
            ]
            for class_index in class_indices
        )
        mixed_strata.append({
            "coefficient_support_module_orbit_label": label,
            "triple_count": len(triple_indices),
            "stabilizer_coordinate_class_count": len(class_indices),
            "class_size_distribution": {
                str(size): count for size, count in sorted(class_size_distribution.items())
            },
            "triple_category_distribution": dict(triple_category_distribution),
            "class_category_distribution": dict(class_category_distribution),
        })

    multi_classes = [
        record for record in c5_classes
        if record["coordinate_intersection_size"] > 1
    ]
    pair_equivalence_reduction = 3276 - len(c5_classes)

    slice_report = {
        "schema": "apn-stage5-slice-stabilizer-orbits-v1",
        "date": "2026-08-02",
        "scope": (
            "Exact stabilizer orbits of full 24-dimensional coefficient-support "
            "slice directions; no APN search"
        ),
        "provenance": {
            "stage5a_file": str(stage5a_path.relative_to(root)),
            "stage5a_sha256": sha256_file(stage5a_path),
            "coefficient_support_module_file": str(module_output_path.relative_to(root)),
            "canonical_correlation_file": str(correlation_path.relative_to(root)),
            "canonical_correlation_sha256": sha256_file(correlation_path),
        },
        "slice_object": {
            "affine_slice": "F + E_T",
            "direction": "E_T = V_out tensor U_T^*",
            "dimension_over_F2": 24,
            "stabilizer_action": (
                "For B fixing F under normalized conjugation, E_T maps to "
                "V_out tensor (U_T^* Lambda^2(B^{-1}))."
            ),
            "coordinate_intersection_meaning": (
                "Only orbit members that again equal a three-coordinate support "
                "subspace correspond to one of the 3276 executed Groebner slices."
            ),
        },
        "centre_orbit_obstruction": {
            "cached_centres_in_distinct_centralizer_orbits": stage5a[
                "cached_membership_conclusions"
            ]["all_27_cached_centres_are_in_distinct_centralizer_orbits"],
            "consequence": (
                "Centralizer pair-orbits cannot connect two different cached "
                "centres.  Cross-centre patterns such as all-11 and x9-only "
                "must be explained after transport to canonical Gold centres "
                "under their larger EA stabilizers, not by C_GL(A) alone."
            ),
        },
        "stabilizers": {
            "generic_26_centres": {
                "structure": "C5=<A>",
                "entries": [matrix_entry_json(entry) for entry in kernel_entries],
            },
            "Gold_x3_seed98": {
                "structure": "C10",
                "entries": [matrix_entry_json(entry) for entry in seed98_entries],
            },
        },
        "coordinate_partition": {
            "total_coordinate_triples": 3276,
            "C5_coordinate_classes": len(c5_classes),
            "C10_coordinate_classes": len(c10_classes),
            "C5_and_C10_coordinate_partitions_identical": same_coordinate_partition,
            "reduction_from_coordinate_equivalence": pair_equivalence_reduction,
            "class_size_distribution": dict(Counter(
                str(record["coordinate_intersection_size"])
                for record in c5_classes
            )),
            "triples_by_class_size": dict(Counter(
                str(c5_classes[c5_lookup[index]]["coordinate_intersection_size"])
                for index in range(3276)
            )),
            "multi_coordinate_class_count": len(multi_classes),
            "maximum_coordinate_intersection_size": max(
                record["coordinate_intersection_size"] for record in c5_classes
            ),
            "multi_coordinate_classes": multi_classes,
        },
        "validation": {
            "gold_seed_solution_count_invariance": solution_invariance,
            "all_11_gold_seeds_pass": True,
            "productivity_category_pure_on_coordinate_classes": True,
            "category_purity_failures": [],
        },
        "mixed_rank_two_strata": mixed_strata,
        "main_conclusions": [
            "The exact centre stabilizers collapse 3276 executed coordinate slices to 3098 coordinate-equivalence classes, a reduction of only 178.",
            "There are 2924 singleton classes, 172 classes meeting two coordinate slices, and two classes meeting four coordinate slices.",
            "Exact solution counts are invariant on every such class for all 11 Gold centres, independently validating the contragredient action and the campaign data.",
            "The extra C10 symmetry of Gold-x3 seed 98 doubles full non-coordinate orbits but creates no additional coordinate-coordinate identifications beyond C5.",
            "Stabilizer slice orbits do not compactly resolve the two mixed rank-two strata: most coordinate classes remain singletons.",
            "Because all cached centres lie in distinct centralizer orbits, the next explanatory step must transport full slices to canonical Gold coordinates and use the full canonical EA stabilizer.",
        ],
    }

    module_output_path.parent.mkdir(parents=True, exist_ok=True)
    module_output_path.write_text(
        json.dumps(module_report, indent=2), encoding="utf-8"
    )
    slice_output_path.parent.mkdir(parents=True, exist_ok=True)
    slice_output_path.write_text(
        json.dumps(slice_report, indent=2), encoding="utf-8"
    )

    print(json.dumps({
        "module_output": str(module_output_path),
        "slice_output": str(slice_output_path),
        "coefficient_support_orbit_count": len(module_orbits),
        "rank_distribution": module_report["summary"]["module_rank_distribution"],
        "pure_negative_filter_excludes": pure_negative_excluded,
        "mixed_rank_two_labels": mixed_rank_two_labels,
        "stabilizer_coordinate_classes": len(c5_classes),
        "coordinate_equivalence_reduction": pair_equivalence_reduction,
        "all_solution_invariance_checks_pass": True,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
