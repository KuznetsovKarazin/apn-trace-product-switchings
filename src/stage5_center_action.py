#!/usr/bin/env python3
"""Stage 5A: exact normalized centralizer action on cached APN centres.

The cached centres are homogeneous quadratic representatives in

    V_A = {F_2-homogeneous quadratic F : F(Ax) = A F(x)}.

Naive conjugation B o F o B^{-1}, for B in C_GL(8,2)(A), preserves APN and
A-equivariance but generally creates a linear ANF part.  Therefore it does not
act literally on the 40-dimensional homogeneous representative space V_A.

This script uses the induced (normalized) conjugation action

    rho_B(F) = Hom_2(B o F o B^{-1}),

where Hom_2 means the homogeneous quadratic ANF component.  If
F(x) = C m(x), with m(x) the 28 square-free quadratic monomials, then

    C |-> B C Lambda^2(B^{-1}).

Equivalently, if G = B o F o B^{-1}, then

    rho_B(F)(x) = G(x) + L_G(x),

where L_G is the unique linear part of G (over F_2, addition removes it).
Adding a linear map preserves differential multiplicities, so APN is
preserved.  This is a genuine group action because linear maps form an
invariant subspace under conjugation.

The script performs generator-level sanity checks and computes centralizer
orbits/stabilizers of the 27 cached centres.  It deliberately does not perform
any full-slice orbit enumeration.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from collections import Counter, deque
from pathlib import Path
from typing import Any, Iterable

N = 8
K = 28
NVARS = N * K
GROUP_ORDER = 61200
EXPECTED_STAGE3_DISTRIBUTION = {
    "all_11_gold": 201,
    "no_gold": 2972,
    "x9_only_all": 100,
    "productive_10_of_11": 3,
}
GL2_GENERATORS = [
    (0x3, 0x0, 0x0, 0x1),
    (0x0, 0x1, 0x1, 0x0),
    (0x1, 0x1, 0x0, 0x1),
]
GENERATOR_NAMES = [
    "diag_primitive",
    "swap",
    "upper_transvection",
]


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


def sha256_lut(lut: Iterable[int]) -> str:
    return hashlib.sha256(bytes(lut)).hexdigest()


def identity_matrix(n: int) -> list[list[int]]:
    return [[int(i == j) for j in range(n)] for i in range(n)]


def mat_mul2(
    left: list[list[int]],
    right: list[list[int]],
) -> list[list[int]]:
    rows = len(left)
    middle = len(right)
    cols = len(right[0])
    out = [[0] * cols for _ in range(rows)]
    for i in range(rows):
        for k in range(middle):
            if not left[i][k]:
                continue
            for j in range(cols):
                out[i][j] ^= right[k][j]
    return out


def inverse2(matrix: list[list[int]]) -> list[list[int]]:
    n = len(matrix)
    augmented = [
        matrix[i][:] + identity_matrix(n)[i]
        for i in range(n)
    ]
    row = 0
    for column in range(n):
        pivot = next(
            (i for i in range(row, n) if augmented[i][column]),
            None,
        )
        if pivot is None:
            raise RuntimeError("Singular binary matrix")
        augmented[row], augmented[pivot] = (
            augmented[pivot],
            augmented[row],
        )
        for i in range(n):
            if i != row and augmented[i][column]:
                augmented[i] = [
                    x ^ y
                    for x, y in zip(augmented[i], augmented[row])
                ]
        row += 1
    return [entry[n:] for entry in augmented]


def rank_bitrows(rows: Iterable[int], width: int) -> int:
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
    return rank


def matrix_rank2(matrix: list[list[int]]) -> int:
    rows = [
        sum((value & 1) << column for column, value in enumerate(row))
        for row in matrix
    ]
    return rank_bitrows(rows, len(matrix[0]) if matrix else 0)


def matrix_to_hex_rows(matrix: list[list[int]]) -> list[str]:
    return [
        f"0x{sum((value & 1) << j for j, value in enumerate(row)):02x}"
        for row in matrix
    ]


def selfequiv_equations(
    A: list[list[int]],
    pairs: list[tuple[int, int]],
) -> list[int]:
    pair_index = {pair: i for i, pair in enumerate(pairs)}
    equations: list[int] = []
    for output in range(N):
        for r, s in pairs:
            equation = 0
            for p, q in pairs:
                coefficient = (
                    (A[p][r] & A[q][s])
                    ^ (A[p][s] & A[q][r])
                )
                if coefficient:
                    equation ^= 1 << (
                        output * K + pair_index[(p, q)]
                    )
            for k in range(N):
                if A[output][k]:
                    equation ^= 1 << (
                        k * K + pair_index[(r, s)]
                    )
            if equation:
                equations.append(equation)
    return equations


def rref_bitrows(
    rows: Iterable[int],
    width: int,
) -> tuple[list[int], list[int]]:
    work = list(rows)
    pivots: list[int] = []
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
        pivots.append(column)
        rank += 1
        if rank == len(work):
            break
    return work, pivots


def build_va_parameterisation(
    A: list[list[int]],
    pairs: list[tuple[int, int]],
) -> tuple[list[int], dict[int, list[int]]]:
    reduced, pivots = rref_bitrows(
        selfequiv_equations(A, pairs),
        NVARS,
    )
    pivot_set = set(pivots)
    free = [index for index in range(NVARS) if index not in pivot_set]
    dependencies = {
        pivot: [free_index for free_index in free if (reduced[row] >> free_index) & 1]
        for row, pivot in enumerate(pivots)
    }
    if len(free) != 40:
        raise RuntimeError(f"Expected dim(V_A)=40, got {len(free)}")
    return free, dependencies


def bits_to_coefficients(
    bits: int,
    free: list[int],
    dependencies: dict[int, list[int]],
) -> list[int]:
    coefficients = [0] * NVARS
    for coordinate, free_index in enumerate(free):
        coefficients[free_index] = (bits >> coordinate) & 1
    for pivot, sources in dependencies.items():
        value = 0
        for source in sources:
            value ^= coefficients[source]
        coefficients[pivot] = value
    return coefficients


def coefficients_to_bits(
    coefficients: list[int],
    free: list[int],
) -> int:
    return sum(
        (coefficients[free_index] & 1) << coordinate
        for coordinate, free_index in enumerate(free)
    )


def coefficients_to_matrix(coefficients: list[int]) -> list[list[int]]:
    return [
        coefficients[output * K:(output + 1) * K]
        for output in range(N)
    ]


def coefficients_to_lut(
    coefficients: list[int],
    pairs: list[tuple[int, int]],
) -> list[int]:
    lut: list[int] = []
    for x in range(1 << N):
        y = 0
        for output in range(N):
            bit = 0
            offset = output * K
            for index, (p, q) in enumerate(pairs):
                bit ^= (
                    coefficients[offset + index]
                    & ((x >> p) & 1)
                    & ((x >> q) & 1)
                )
            y |= bit << output
        lut.append(y)
    return lut


def is_apn(lut: list[int]) -> bool:
    for difference in range(1, 1 << N):
        counts = [0] * (1 << N)
        for x in range(1 << N):
            value = lut[x] ^ lut[x ^ difference]
            counts[value] += 1
            if counts[value] > 2:
                return False
    return True


def linear_columns_from_lut(lut: list[int]) -> list[int]:
    return [lut[1 << i] for i in range(N)]


def apply_linear_columns(columns: list[int], value: int) -> int:
    output = 0
    for i, image in enumerate(columns):
        if (value >> i) & 1:
            output ^= image
    return output


def linear_columns_to_matrix(columns: list[int]) -> list[list[int]]:
    return [
        [((columns[column] >> row) & 1) for column in range(N)]
        for row in range(N)
    ]


def direct_conjugate_lut(
    stage3,
    lut: list[int],
    B: list[list[int]],
) -> list[int]:
    inverse = inverse2(B)
    input_map = [stage3.mat_apply_binary(inverse, x) for x in range(1 << N)]
    output_map = [stage3.mat_apply_binary(B, y) for y in range(1 << N)]
    return [output_map[lut[input_map[x]]] for x in range(1 << N)]


def normalized_conjugate_lut(
    stage3,
    lut: list[int],
    B: list[list[int]],
) -> tuple[list[int], list[int], list[int]]:
    direct = direct_conjugate_lut(stage3, lut, B)
    linear_columns = linear_columns_from_lut(direct)
    normalized = [
        direct[x] ^ apply_linear_columns(linear_columns, x)
        for x in range(1 << N)
    ]
    return normalized, direct, linear_columns


def coefficient_action(
    stage3,
    coefficients: list[int],
    B: list[list[int]],
) -> list[int]:
    inverse = inverse2(B)
    wedge_inverse = stage3.wedge_action(inverse)
    C = coefficients_to_matrix(coefficients)
    transformed = mat_mul2(mat_mul2(B, C), wedge_inverse)
    return [value for row in transformed for value in row]


def action_columns(
    stage3,
    B: list[list[int]],
    free: list[int],
    dependencies: dict[int, list[int]],
) -> tuple[int, ...]:
    columns = []
    for coordinate in range(40):
        coefficients = bits_to_coefficients(
            1 << coordinate,
            free,
            dependencies,
        )
        transformed = coefficient_action(stage3, coefficients, B)
        columns.append(coefficients_to_bits(transformed, free))
    return tuple(columns)


def apply_action(columns: tuple[int, ...], bits: int) -> int:
    output = 0
    value = bits
    while value:
        low = value & -value
        output ^= columns[low.bit_length() - 1]
        value ^= low
    return output


def compose_actions(
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> tuple[int, ...]:
    """Return left after right."""
    return tuple(apply_action(left, column) for column in right)


def action_rank(columns: tuple[int, ...]) -> int:
    return rank_bitrows(columns, 40)


def action_order(columns: tuple[int, ...], limit: int = 1000) -> int:
    identity = tuple(1 << i for i in range(40))
    current = identity
    for order in range(1, limit + 1):
        current = compose_actions(columns, current)
        if current == identity:
            return order
    raise RuntimeError("Action order exceeds limit")


def gf16_matrix_mul(stage3, left, right):
    a, b, c, d = left
    e, f, g, h = right
    mul = stage3.f16_mul
    return (
        mul(a, e) ^ mul(b, g),
        mul(a, f) ^ mul(b, h),
        mul(c, e) ^ mul(d, g),
        mul(c, f) ^ mul(d, h),
    )


def gf16_matrix_det(stage3, matrix) -> int:
    a, b, c, d = matrix
    return stage3.f16_mul(a, d) ^ stage3.f16_mul(b, c)


def gf16_matrix_order(stage3, matrix, limit: int = 1000) -> int:
    identity = (1, 0, 0, 1)
    current = identity
    for order in range(1, limit + 1):
        current = gf16_matrix_mul(stage3, current, matrix)
        if current == identity:
            return order
    raise RuntimeError("GF(16) matrix order exceeds limit")


def gf16_matrix_power(stage3, matrix, exponent: int):
    result = (1, 0, 0, 1)
    base = matrix
    value = exponent
    while value:
        if value & 1:
            result = gf16_matrix_mul(stage3, result, base)
        base = gf16_matrix_mul(stage3, base, base)
        value >>= 1
    return result


def enumerate_gl2(stage3, generators) -> set[tuple[int, int, int, int]]:
    identity = (1, 0, 0, 1)
    seen = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for generator in generators:
            image = gf16_matrix_mul(stage3, current, generator)
            if image not in seen:
                if gf16_matrix_det(stage3, image) == 0:
                    raise RuntimeError("Generated singular GF(16) matrix")
                seen.add(image)
                queue.append(image)
    return seen


def enumerate_action_image(
    stage3,
    generators,
    generator_actions,
) -> dict[tuple[int, ...], tuple[int, int, int, int]]:
    identity_action = tuple(1 << i for i in range(40))
    identity_matrix = (1, 0, 0, 1)
    representative = {identity_action: identity_matrix}
    queue = deque([identity_action])
    while queue:
        current_action = queue.popleft()
        current_matrix = representative[current_action]
        for generator, generator_action in zip(
            generators,
            generator_actions,
        ):
            image_action = compose_actions(
                current_action,
                generator_action,
            )
            image_matrix = gf16_matrix_mul(
                stage3,
                current_matrix,
                generator,
            )
            if image_action not in representative:
                representative[image_action] = image_matrix
                queue.append(image_action)
    return representative


def matrix_entry_json(entry) -> list[list[str]]:
    a, b, c, d = entry
    return [
        [f"0x{a:x}", f"0x{b:x}"],
        [f"0x{c:x}", f"0x{d:x}"],
    ]


def binary_matrix_json(matrix: list[list[int]]) -> dict[str, Any]:
    return {
        "rows": matrix,
        "row_hex": matrix_to_hex_rows(matrix),
    }


def omega_matrix_hex(matrix: list[list[int]]) -> list[list[str]]:
    return [[f"0x{value:x}" for value in row] for row in matrix]


def canonical_orbit_label(size: int, minimum_bits: int) -> str:
    return f"CEN-O{size}-B{minimum_bits:010x}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--centres",
        type=Path,
        default=Path("data/campaign_centers/campaign_centres_27.json"),
    )
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
        "--canonical-stage3",
        type=Path,
        default=Path("results/current/module_productivity_correlation_canonical.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current/stage5_center_action.json"),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]

    def resolve(path: Path) -> Path:
        return path if path.is_absolute() else root / path

    centres_path = resolve(args.centres)
    stage3_path = resolve(args.stage3)
    stage4_path = resolve(args.stage4)
    canonical_stage3_path = resolve(args.canonical_stage3)
    output_path = resolve(args.output)

    stage3 = load_module(stage3_path, "stage3_center_action")
    stage4 = load_module(stage4_path, "stage4_center_action")

    canonical_stage3 = json.loads(
        canonical_stage3_path.read_text(encoding="utf-8")
    )
    category_distribution = canonical_stage3["campaign"]["category_distribution"]
    if category_distribution != EXPECTED_STAGE3_DISTRIBUTION:
        raise RuntimeError(
            "Canonical Stage 3 distribution mismatch: "
            f"{category_distribution}"
        )

    centre_payload = json.loads(centres_path.read_text(encoding="utf-8"))
    centres = centre_payload["centres"]
    if len(centres) != 27:
        raise RuntimeError(f"Expected 27 centres, got {len(centres)}")

    pairs = stage3.PAIRS
    free, dependencies = build_va_parameterisation(stage3.A_ROWS, pairs)
    wedge_A = stage3.wedge_action(stage3.A_ROWS)

    source_validation = []
    cache_by_bits = {centre["bits"]: centre for centre in centres}
    cache_by_hash = {centre["sbox_sha256"]: centre for centre in centres}
    for centre in centres:
        coefficients = bits_to_coefficients(
            centre["bits"],
            free,
            dependencies,
        )
        coefficient_matrix = coefficients_to_matrix(coefficients)
        lut = coefficients_to_lut(coefficients, pairs)
        checks = {
            "stored_coefficients_match_bits": (
                coefficients == centre["coefficients_224"]
            ),
            "stored_lut_matches_coefficients": lut == centre["sbox_256"],
            "stored_hash_matches_lut": (
                sha256_lut(lut) == centre["sbox_sha256"]
            ),
            "V_A_intertwiner_equation": (
                mat_mul2(coefficient_matrix, wedge_A)
                == mat_mul2(stage3.A_ROWS, coefficient_matrix)
            ),
            "homogeneous_zero_on_basis": all(
                lut[1 << i] == 0 for i in range(N)
            ),
            "APN": is_apn(lut),
        }
        if not all(checks.values()):
            raise RuntimeError(
                f"Source validation failed for seed {centre['seed']}: {checks}"
            )
        source_validation.append({
            "seed": centre["seed"],
            "exact_class": centre["exact_class"],
            "checks": checks,
        })

    binary_generators = [
        stage4.gl2_to_binary(stage3, *entry)
        for entry in GL2_GENERATORS
    ]
    binary_inverses = [inverse2(matrix) for matrix in binary_generators]
    generator_actions = [
        action_columns(stage3, matrix, free, dependencies)
        for matrix in binary_generators
    ]
    inverse_actions = [
        action_columns(stage3, matrix, free, dependencies)
        for matrix in binary_inverses
    ]

    # Exact group- and representation-level checks.
    generated_gl2 = enumerate_gl2(stage3, GL2_GENERATORS)
    if len(generated_gl2) != GROUP_ORDER:
        raise RuntimeError(
            f"Expected |GL_2(GF(16))|={GROUP_ORDER}, got {len(generated_gl2)}"
        )

    action_image = enumerate_action_image(
        stage3,
        GL2_GENERATORS,
        generator_actions,
    )
    action_image_order = len(action_image)
    if GROUP_ORDER % action_image_order:
        raise RuntimeError("Action image order does not divide group order")
    action_kernel_size = GROUP_ORDER // action_image_order

    scalar_kernel_entries = [
        (scalar, 0, 0, scalar)
        for scalar in [
            stage3.f16_pow(stage3.OMEGA, exponent)
            for exponent in range(5)
        ]
    ]
    scalar_kernel_actions = [
        action_columns(
            stage3,
            stage4.gl2_to_binary(stage3, *entry),
            free,
            dependencies,
        )
        for entry in scalar_kernel_entries
    ]
    identity_action = tuple(1 << i for i in range(40))
    if not all(action == identity_action for action in scalar_kernel_actions):
        raise RuntimeError("<A> is not in the action kernel")
    if action_kernel_size != len(set(scalar_kernel_entries)):
        raise RuntimeError(
            "Action kernel is larger than the verified scalar C5 kernel"
        )

    decomposition, diagonalising, diagonalising_inverse, offsets = (
        stage4.reconstruct_diagonalisation(stage3)
    )

    generator_records = []
    generator_image_records = []
    total_normalized_checks = 0
    total_apn_checks = 0
    total_va_checks = 0
    total_linear_centralizer_checks = 0
    total_direct_homogeneous = 0

    for index, (
        name,
        entry,
        B,
        B_inverse,
        action,
        inverse_action,
    ) in enumerate(zip(
        GENERATOR_NAMES,
        GL2_GENERATORS,
        binary_generators,
        binary_inverses,
        generator_actions,
        inverse_actions,
    )):
        wedge_B = stage3.wedge_action(B)
        wedge_B_inverse = stage3.wedge_action(B_inverse)
        omega_representation = stage4.representation5(
            stage3,
            diagonalising,
            diagonalising_inverse,
            offsets,
            entry,
        )

        inverse_check = (
            compose_actions(inverse_action, action) == identity_action
            and compose_actions(action, inverse_action) == identity_action
        )
        record = {
            "generator_index": index,
            "name": name,
            "GL2_GF16_entry": matrix_entry_json(entry),
            "GF16_determinant_hex": f"0x{gf16_matrix_det(stage3, entry):x}",
            "GF16_matrix_order": gf16_matrix_order(stage3, entry),
            "binary_matrix": binary_matrix_json(B),
            "binary_inverse": binary_matrix_json(B_inverse),
            "binary_matrix_order": stage3.matrix_order_binary(B),
            "commutes_with_A": (
                mat_mul2(B, stage3.A_ROWS)
                == mat_mul2(stage3.A_ROWS, B)
            ),
            "wedge_commutes_with_wedge_A": (
                mat_mul2(wedge_B, wedge_A)
                == mat_mul2(wedge_A, wedge_B)
            ),
            "wedge_inverse_is_inverse": (
                mat_mul2(wedge_B, wedge_B_inverse)
                == identity_matrix(K)
                and mat_mul2(wedge_B_inverse, wedge_B)
                == identity_matrix(K)
            ),
            "action_rank_on_V_A": action_rank(action),
            "action_order_on_V_A": action_order(action),
            "inverse_action_check": inverse_check,
            "omega_multiplicity_representation": omega_matrix_hex(
                omega_representation
            ),
            "centre_checks": {
                "tested": len(centres),
                "coefficient_vs_normalized_LUT_agree": 0,
                "APN_preserved": 0,
                "V_A_preserved": 0,
                "linear_correction_commutes_with_A": 0,
                "direct_conjugate_already_homogeneous": 0,
                "linear_correction_rank_distribution": {},
                "cached_generator_image_matches": [],
            },
        }

        if not all([
            record["commutes_with_A"],
            record["wedge_commutes_with_wedge_A"],
            record["wedge_inverse_is_inverse"],
            record["action_rank_on_V_A"] == 40,
            record["inverse_action_check"],
        ]):
            raise RuntimeError(f"Generator-level matrix check failed: {name}")

        correction_ranks = Counter()
        cached_matches = []
        for centre in centres:
            transformed_bits = apply_action(action, centre["bits"])
            transformed_coefficients = bits_to_coefficients(
                transformed_bits,
                free,
                dependencies,
            )
            transformed_matrix = coefficients_to_matrix(
                transformed_coefficients
            )
            transformed_lut = coefficients_to_lut(
                transformed_coefficients,
                pairs,
            )
            normalized_lut, direct_lut, linear_columns = (
                normalized_conjugate_lut(
                    stage3,
                    centre["sbox_256"],
                    B,
                )
            )
            linear_matrix = linear_columns_to_matrix(linear_columns)
            correction_rank = matrix_rank2(linear_matrix)
            correction_ranks[correction_rank] += 1

            coefficient_agreement = normalized_lut == transformed_lut
            apn_preserved = is_apn(transformed_lut)
            va_preserved = (
                mat_mul2(transformed_matrix, wedge_A)
                == mat_mul2(stage3.A_ROWS, transformed_matrix)
            )
            correction_centralizes_A = (
                mat_mul2(linear_matrix, stage3.A_ROWS)
                == mat_mul2(stage3.A_ROWS, linear_matrix)
            )
            direct_homogeneous = all(value == 0 for value in linear_columns)

            if not all([
                coefficient_agreement,
                apn_preserved,
                va_preserved,
                correction_centralizes_A,
            ]):
                raise RuntimeError(
                    f"Generator action failed for {name}, seed {centre['seed']}"
                )

            record["centre_checks"][
                "coefficient_vs_normalized_LUT_agree"
            ] += int(coefficient_agreement)
            record["centre_checks"]["APN_preserved"] += int(apn_preserved)
            record["centre_checks"]["V_A_preserved"] += int(va_preserved)
            record["centre_checks"][
                "linear_correction_commutes_with_A"
            ] += int(correction_centralizes_A)
            record["centre_checks"][
                "direct_conjugate_already_homogeneous"
            ] += int(direct_homogeneous)

            image_hash = sha256_lut(transformed_lut)
            direct_hash = sha256_lut(direct_lut)
            cached_image = cache_by_bits.get(transformed_bits)
            if cached_image is not None:
                cached_matches.append({
                    "source_seed": centre["seed"],
                    "target_seed": cached_image["seed"],
                })

            generator_image_records.append({
                "generator": name,
                "source_seed": centre["seed"],
                "source_class": centre["exact_class"],
                "source_bits_hex": f"0x{centre['bits']:010x}",
                "normalized_image_bits_hex": f"0x{transformed_bits:010x}",
                "normalized_image_sha256": image_hash,
                "direct_conjugate_sha256": direct_hash,
                "linear_correction_rank": correction_rank,
                "linear_correction_columns_hex": [
                    f"0x{value:02x}" for value in linear_columns
                ],
                "direct_conjugate_is_homogeneous": direct_homogeneous,
                "cached_target_seed": (
                    cached_image["seed"] if cached_image else None
                ),
                "checks": {
                    "coefficient_formula_equals_normalized_conjugation": (
                        coefficient_agreement
                    ),
                    "APN": apn_preserved,
                    "V_A": va_preserved,
                    "linear_correction_in_C_End(A)": (
                        correction_centralizes_A
                    ),
                },
            })

            total_normalized_checks += int(coefficient_agreement)
            total_apn_checks += int(apn_preserved)
            total_va_checks += int(va_preserved)
            total_linear_centralizer_checks += int(
                correction_centralizes_A
            )
            total_direct_homogeneous += int(direct_homogeneous)

        record["centre_checks"]["linear_correction_rank_distribution"] = {
            str(rank): count
            for rank, count in sorted(correction_ranks.items())
        }
        record["centre_checks"]["cached_generator_image_matches"] = (
            cached_matches
        )
        generator_records.append(record)

    # Pairwise group-law and wedge-functor checks for the three generators.
    group_law_checks = []
    for i, (left_entry, left_action) in enumerate(zip(
        GL2_GENERATORS,
        generator_actions,
    )):
        for j, (right_entry, right_action) in enumerate(zip(
            GL2_GENERATORS,
            generator_actions,
        )):
            product_entry = gf16_matrix_mul(
                stage3,
                left_entry,
                right_entry,
            )
            product_binary = stage4.gl2_to_binary(
                stage3,
                *product_entry,
            )
            direct_product_action = action_columns(
                stage3,
                product_binary,
                free,
                dependencies,
            )
            composed_action = compose_actions(
                left_action,
                right_action,
            )
            binary_product = mat_mul2(
                binary_generators[i],
                binary_generators[j],
            )
            wedge_product = mat_mul2(
                stage3.wedge_action(binary_generators[i]),
                stage3.wedge_action(binary_generators[j]),
            )
            check = {
                "left": GENERATOR_NAMES[i],
                "right": GENERATOR_NAMES[j],
                "binary_embedding_homomorphism": (
                    product_binary == binary_product
                ),
                "wedge_functor_homomorphism": (
                    stage3.wedge_action(binary_product) == wedge_product
                ),
                "normalized_action_group_law": (
                    direct_product_action == composed_action
                ),
            }
            if not all(value for key, value in check.items() if key not in {"left", "right"}):
                raise RuntimeError(f"Group law check failed: {check}")
            group_law_checks.append(check)

    # Exact action kernel and image.
    kernel_entries = sorted(
        scalar_kernel_entries,
        key=lambda entry: tuple(entry),
    )

    # Centre orbits under the exact 40-dimensional action.
    cached_assignment: dict[int, str] = {}
    orbit_records = []
    centre_orbit_lookup: dict[int, dict[str, Any]] = {}
    for centre in centres:
        start = centre["bits"]
        if start in cached_assignment:
            continue
        seen = {start}
        queue = deque([start])
        while queue:
            current = queue.popleft()
            for generator_action in generator_actions:
                image = apply_action(generator_action, current)
                if image not in seen:
                    seen.add(image)
                    queue.append(image)

        minimum_bits = min(seen)
        label = canonical_orbit_label(len(seen), minimum_bits)
        cached_members = sorted(
            (
                cache_by_bits[value]["seed"],
                cache_by_bits[value]["exact_class"],
                cache_by_bits[value]["sbox_sha256"],
            )
            for value in seen
            if value in cache_by_bits
        )
        if GROUP_ORDER % len(seen):
            raise RuntimeError(
                f"Orbit size {len(seen)} does not divide {GROUP_ORDER}"
            )
        stabilizer_size = GROUP_ORDER // len(seen)
        orbit_record = {
            "orbit_label": label,
            "canonical_minimum_bits_hex": f"0x{minimum_bits:010x}",
            "orbit_size": len(seen),
            "stabilizer_size_in_GL2_GF16": stabilizer_size,
            "stabilizer_size_in_faithful_action_image": (
                action_image_order // len(seen)
            ),
            "cached_members": [
                {
                    "seed": seed,
                    "exact_class": exact_class,
                    "sbox_sha256": digest,
                }
                for seed, exact_class, digest in cached_members
            ],
        }
        orbit_records.append(orbit_record)
        for seed, _, _ in cached_members:
            member = next(c for c in centres if c["seed"] == seed)
            cached_assignment[member["bits"]] = label
            centre_orbit_lookup[seed] = orbit_record

    orbit_records.sort(key=lambda record: (
        record["orbit_size"],
        record["canonical_minimum_bits_hex"],
    ))

    if len(cached_assignment) != len(centres):
        raise RuntimeError("Not all cached centres were assigned to an orbit")

    # Stabilizer structure from the faithful action image.
    stabilizer_details: dict[int, dict[str, Any]] = {}
    for centre in centres:
        image_fixes = [
            (action, representative)
            for action, representative in action_image.items()
            if apply_action(action, centre["bits"]) == centre["bits"]
        ]
        full_stabilizer_size = len(image_fixes) * action_kernel_size
        expected_size = centre_orbit_lookup[centre["seed"]][
            "stabilizer_size_in_GL2_GF16"
        ]
        if full_stabilizer_size != expected_size:
            raise RuntimeError(
                f"Stabilizer mismatch for seed {centre['seed']}"
            )
        nonidentity_representatives = [
            representative
            for action, representative in image_fixes
            if action != identity_action
        ]
        detail: dict[str, Any] = {
            "full_stabilizer_size": full_stabilizer_size,
            "faithful_image_stabilizer_size": len(image_fixes),
            "contains_global_kernel_C5": True,
            "nonkernel_representatives": [
                {
                    "GL2_GF16_entry": matrix_entry_json(entry),
                    "order_in_GL2_GF16": gf16_matrix_order(stage3, entry),
                    "square": matrix_entry_json(
                        gf16_matrix_power(stage3, entry, 2)
                    ),
                }
                for entry in nonidentity_representatives
            ],
        }
        if full_stabilizer_size == 5:
            detail["structure"] = "C5 = <A>, the global action kernel"
        elif full_stabilizer_size == 10 and len(nonidentity_representatives) == 1:
            generator = nonidentity_representatives[0]
            powers = [
                gf16_matrix_power(stage3, generator, exponent)
                for exponent in range(10)
            ]
            detail["structure"] = "C10"
            detail["cyclic_generator"] = matrix_entry_json(generator)
            detail["cyclic_generator_order"] = gf16_matrix_order(
                stage3,
                generator,
            )
            detail["all_stabilizer_elements"] = [
                matrix_entry_json(entry) for entry in powers
            ]
        else:
            detail["structure"] = "not automatically identified"
        stabilizer_details[centre["seed"]] = detail

    centre_records = []
    for centre in centres:
        orbit = centre_orbit_lookup[centre["seed"]]
        centre_records.append({
            "seed": centre["seed"],
            "exact_class": centre["exact_class"],
            "bits_hex": f"0x{centre['bits']:010x}",
            "sbox_sha256": centre["sbox_sha256"],
            "orbit_label": orbit["orbit_label"],
            "orbit_size": orbit["orbit_size"],
            "stabilizer": stabilizer_details[centre["seed"]],
            "cached_members_in_same_orbit": [
                member["seed"] for member in orbit["cached_members"]
            ],
        })

    class_orbits: dict[str, list[str]] = {}
    for exact_class in sorted({c["exact_class"] for c in centres}):
        labels = sorted({
            centre_orbit_lookup[c["seed"]]["orbit_label"]
            for c in centres
            if c["exact_class"] == exact_class
        })
        class_orbits[exact_class] = labels

    cached_generator_match_count = sum(
        record["cached_target_seed"] is not None
        for record in generator_image_records
    )

    report = {
        "schema": "apn-stage5-center-action-v1",
        "date": "2026-08-02",
        "scope": (
            "Stage 5A only: exact centre action and stabilizers; "
            "no full-slice orbit enumeration"
        ),
        "provenance": {
            "centres_file": str(centres_path.relative_to(root)),
            "centres_sha256": sha256_file(centres_path),
            "canonical_stage3_file": str(
                canonical_stage3_path.relative_to(root)
            ),
            "canonical_stage3_sha256": sha256_file(
                canonical_stage3_path
            ),
            "canonical_stage3_category_distribution": category_distribution,
            "stage3_script": str(stage3_path.relative_to(root)),
            "stage3_script_sha256": sha256_file(stage3_path),
            "stage4_script": str(stage4_path.relative_to(root)),
            "stage4_script_sha256": sha256_file(stage4_path),
        },
        "action_convention": {
            "vector_convention": (
                "column vectors over F_2; x maps to Bx"
            ),
            "function_representation": (
                "F(x)=C m(x), C is 8x28, m lists x_i x_j for i<j"
            ),
            "V_A_equation": "C Lambda^2(A) = A C",
            "naive_conjugation": "G_B = B o F o B^{-1}",
            "naive_conjugation_issue": (
                "G_B generally acquires an A-commuting linear ANF part, "
                "so naive conjugation does not preserve the homogeneous "
                "40-dimensional representative space V_A"
            ),
            "normalized_action": (
                "rho_B(F)=Hom_2(B o F o B^{-1})"
            ),
            "coefficient_formula": (
                "C maps to B C Lambda^2(B^{-1})"
            ),
            "LUT_formula": (
                "rho_B(F)(x)=G_B(x)+L_{B,F}(x), where columns of "
                "L_{B,F} are G_B(e_i)"
            ),
            "why_APN_is_preserved": (
                "input/output linear conjugation preserves APN, and adding "
                "a linear output map only translates every derivative by "
                "the constant L(a)"
            ),
            "group_action_reason": (
                "linear maps are conjugation-invariant and Hom_2 annihilates "
                "them, so rho_{B1 B2}=rho_{B1} rho_{B2}"
            ),
            "Lambda2_note": (
                "Stage 4 uses the forward action Lambda^2(B) on monomial "
                "vectors. Centre coefficients use the inverse domain action "
                "Lambda^2(B^{-1}). A coordinate coefficient-support slice "
                "is naturally dual and therefore transforms contragrediently; "
                "Stage 5B must freeze this dual/support convention before "
                "any full-slice enumeration."
            ),
        },
        "V_A": {
            "dimension": len(free),
            "ambient_quadratic_coefficient_dimension": NVARS,
            "A_binary_matrix": binary_matrix_json(stage3.A_ROWS),
            "A_order": stage3.matrix_order_binary(stage3.A_ROWS),
            "A_as_GF16_scalar": "0x2 I_2",
            "free_coordinate_indices": free,
        },
        "group": {
            "name": "C_GL(8,2)(A) = GL_2(GF(16))",
            "order": GROUP_ORDER,
            "standard_generators": generator_records,
            "generated_group_order": len(generated_gl2),
            "normalized_action_image_order": action_image_order,
            "normalized_action_kernel_size": action_kernel_size,
            "normalized_action_kernel": {
                "structure": "C5 = <A>",
                "entries": [matrix_entry_json(entry) for entry in kernel_entries],
            },
            "group_law_checks": group_law_checks,
        },
        "source_centre_validation": {
            "tested": len(source_validation),
            "all_passed": True,
            "records": source_validation,
        },
        "generator_level_summary": {
            "centre_generator_pairs_tested": len(generator_image_records),
            "coefficient_vs_normalized_LUT_agree": total_normalized_checks,
            "APN_preserved": total_apn_checks,
            "V_A_preserved": total_va_checks,
            "linear_correction_commutes_with_A": (
                total_linear_centralizer_checks
            ),
            "direct_conjugate_already_homogeneous": total_direct_homogeneous,
            "cached_generator_image_matches": cached_generator_match_count,
        },
        "generator_images": generator_image_records,
        "centre_orbits": {
            "cached_centres": len(centres),
            "distinct_centralizer_orbits_meeting_cache": len(orbit_records),
            "orbit_size_distribution": dict(Counter(
                str(record["orbit_size"]) for record in orbit_records
            )),
            "stabilizer_size_distribution": dict(Counter(
                str(record["stabilizer_size_in_GL2_GF16"])
                for record in orbit_records
            )),
            "orbits": orbit_records,
            "centres": centre_records,
        },
        "cached_membership_conclusions": {
            "all_27_cached_centres_are_in_distinct_centralizer_orbits": (
                len(orbit_records) == len(centres)
                and all(len(record["cached_members"]) == 1 for record in orbit_records)
            ),
            "Gold_x3_five_centres_in_one_orbit": (
                len(class_orbits.get("Gold-x3", [])) == 1
            ),
            "Gold_x9_six_centres_in_one_orbit": (
                len(class_orbits.get("Gold-x9", [])) == 1
            ),
            "Kasami_x57_sixteen_centres_in_one_orbit": (
                len(class_orbits.get("Kasami-x57", [])) == 1
            ),
            "distinct_orbit_count_by_class": {
                exact_class: len(labels)
                for exact_class, labels in class_orbits.items()
            },
            "cached_set_closed_under_standard_generators": (
                cached_generator_match_count == len(generator_image_records)
            ),
            "cached_standard_generator_images_inside_cache": (
                cached_generator_match_count
            ),
        },
        "main_conclusions": [
            (
                "Naive input/output conjugation is not an action on the "
                "homogeneous 40-dimensional V_A representatives because it "
                "usually creates a linear part."
            ),
            (
                "The exact action on cached centres is normalized conjugation, "
                "equivalently C -> B C Lambda^2(B^{-1})."
            ),
            (
                "All 81 generator-centre checks agree exactly between the "
                "coefficient formula and normalized LUT conjugation; APN and "
                "V_A are preserved in every case."
            ),
            (
                "The action kernel is the scalar C5=<A>; the faithful image "
                "has order 12240."
            ),
            (
                "The 27 cached centres meet 27 distinct centralizer orbits. "
                "Thus neither the five Gold-x3 centres nor the six Gold-x9 "
                "centres form a single centralizer orbit."
            ),
            (
                "Twenty-six centres have orbit size 12240 and stabilizer C5. "
                "Gold-x3 seed 98 has orbit size 6120 and cyclic stabilizer C10."
            ),
            (
                "Full-slice enumeration remains blocked until the coefficient-"
                "support dual/contragredient convention is explicitly frozen."
            ),
        ],
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(json.dumps({
        "output": str(output_path),
        "generated_group_order": len(generated_gl2),
        "action_image_order": action_image_order,
        "action_kernel_size": action_kernel_size,
        "generator_centre_checks": len(generator_image_records),
        "distinct_cached_orbits": len(orbit_records),
        "orbit_size_distribution": report["centre_orbits"][
            "orbit_size_distribution"
        ],
        "stabilizer_size_distribution": report["centre_orbits"][
            "stabilizer_size_distribution"
        ],
        "class_orbit_counts": report["cached_membership_conclusions"][
            "distinct_orbit_count_by_class"
        ],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
