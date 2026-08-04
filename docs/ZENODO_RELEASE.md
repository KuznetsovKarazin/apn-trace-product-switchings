# GitHub–Zenodo release procedure

1. Create the public GitHub repository `KuznetsovKarazin/apn-trace-product-switchings`.
2. Push this directory to the `main` branch.
3. In Zenodo, enable the GitHub integration for this repository.
4. Create GitHub release `v1.0.0` with title `Manuscript companion release`.
5. Zenodo will mint a version DOI and a concept DOI.
6. Replace `10.5281/zenodo.XXXXXXX` in:
   - `README.md`;
   - `CITATION.cff`;
   - `manuscript/main.tex`;
   - `manuscript/Supplementary_Material.tex` if present.
7. Commit the DOI update and create `v1.0.1`, or reserve the Zenodo DOI before publication and insert it into `v1.0.0` metadata if using a manual Zenodo deposit.

Recommended release archive contents: the complete repository, excluding `.git` and local virtual environments.
