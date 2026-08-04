# APN Trace-Product Switchings

Reproducible code, exact certificates, and manuscript sources for:

**Dimension Rigidity and Projective Geometry of Trace-Product Switchings of the Gold Cube**  
Oleksandr Kuznetsov

## Main results

For the scalar trace-product switching of the Gold cube over `GF(2^n)`:

| even dimension | admissible nonzero coefficients |
|---|---|
| `n = 4` | nonzero trace-zero elements |
| `n = 6` | the six elements of multiplicative order 9 |
| `n = 8` | `GF(4)^*` |
| `n >= 10` | none |

The repository also contains the exact dimension-eight projective rank-two classification: eight marked switchings forming exactly two EA/CCZ classes.

## Scope and relation to the earlier project

This is a **separate companion repository** for the theoretical trace-product switching paper. The earlier repository [`apn-gb-search`](https://github.com/KuznetsovKarazin/apn-gb-search) records the Gröbner-slice discovery pipeline and its historical campaign. The two repositories should not be merged: they have different claims, dependencies, and reproducibility targets.

The raw coefficient lists in dimensions 6 and 8 were previously computed in Razi Arshad's 2018 dissertation. This project supplies intrinsic descriptions, proofs across all even dimensions, the uniform nonexistence theorem, and the projective rank-two geometry.

## Repository layout

- `src/` — exact Python scripts used for all finite certificates and audits.
- `results/` — canonical JSON outputs.
- `manuscript/` — Springer/DCC LaTeX source, compiled manuscript, and supplement.
- `docs/` — proof status, Zenodo instructions, and relationship to the earlier search project.

## Reproduction

Python 3.11+ is recommended.

```bash
python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Core checks:

```bash
python src/stage5p_kummer_n6_certificate.py
python src/stage5s_projective_equivalence_partition.py
python src/stage5t_compact_orbit_certificates.py
python src/stage5v_h_equivalence_comparison.py
python src/stage5y_publication_uniform_bound_correction.py
```

The outputs are deterministic and are written to `results/current/`. Verify the committed or regenerated core certificates with:

```bash
python src/verify_release.py
```

The GitHub Actions workflow reruns the core scripts, checks their SHA-256 values, and requires byte-identical canonical JSON outputs.

## Citation

Until the Zenodo deposit is published, use the repository URL. Replace the placeholder DOI after release:

```bibtex
@software{kuznetsov2026traceproduct,
  author  = {Oleksandr Kuznetsov},
  title   = {APN Trace-Product Switchings: Code and Exact Certificates},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/KuznetsovKarazin/apn-trace-product-switchings},
  note    = {Zenodo DOI to be added after the first archived release}
}
```

## License

Code: MIT. Documentation, manuscript-side data, and exact certificates: CC BY 4.0.
