# GitHub Actions path hotfix

The initial repository bundle flattened the canonical JSON files into `results/`, while the scripts retained their original `results/current/` paths. In addition, scripts copied from the checkpoint layout resolved the project root with `parents[2]`; after moving them to the top-level `src/` directory, the correct root is `parents[1]`.

This release fixes both packaging defects:

1. canonical JSON files are stored in `results/current/`;
2. all public scripts resolve the repository root with `Path(__file__).resolve().parents[1]`;
3. the non-portable internal manuscript builder with a hard-coded `/mnt/data/...` path has been removed;
4. GitHub Actions now validates committed hashes, reruns five core certificates, and checks byte-identical outputs;
5. stale release manifests and the invalid placeholder DOI in `CITATION.cff` have been corrected.

The mathematical results and canonical JSON hashes are unchanged.
