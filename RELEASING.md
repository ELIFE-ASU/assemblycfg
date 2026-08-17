# Releasing `assemblycfg`

Releases are published to PyPI from a GitHub Release through Trusted
Publishing. The publisher is tied to `.github/workflows/python-publish.yml` and
the `pypi` GitHub environment; keep both names in sync with the PyPI publisher
configuration.

1. Update `assemblycfg/_version.py`, `CITATION.cff` (`version` and
   `date-released`), and any release notes. Use a version that has not already
   been tagged or uploaded to PyPI. Confirm whether the Zenodo identifier is a
   concept DOI or a version DOI; reserve and update a version DOI for the new
   archive when applicable, and keep the README, `CITATION.cff`, and
   `pyproject.toml` links consistent.
2. Install the development dependencies and run the same checks as CI:

   Start with an empty `dist/` directory so older artifacts cannot be mixed into
   the release.

   ```bash
   python -m pip install --upgrade pip
   python -m pip install --group dev -e .
   python -m pytest
   python -m build
   python -m twine check --strict dist/*
   python scripts/check_dist.py --tag vX.Y.Z
   ```

3. Merge the release commit only after CI is green.
4. Create and push the matching annotated tag, for example `v1.2.4`.
5. Publish a GitHub Release for that tag. The publishing workflow reruns the
   tests, builds the artifacts once, verifies the tag/source/CITATION versions,
   and uploads those exact artifacts to PyPI.
6. Confirm that the GitHub environment deployment and the new PyPI files and
   provenance attestations succeeded.

Do not move or reuse a release tag. PyPI versions and files are immutable.
