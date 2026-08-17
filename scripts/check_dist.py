#!/usr/bin/env python3
"""Validate release artifacts using only the Python standard library."""

from __future__ import annotations

import argparse
from email.parser import Parser
from pathlib import Path, PurePosixPath
import re
import runpy
import tarfile
from zipfile import ZipFile


PROJECT_NAME = "assemblycfg"
DIST_FILENAME_NAME = PROJECT_NAME.replace("-", "_")
EXPECTED_EXTRAS = {"dev", "plot"}
EXPECTED_SDIST_FILES = {
    "CITATION.cff",
    "LICENSE",
    "README.md",
    "RELEASING.md",
    "cfg.bib",
    "examples/example.py",
    "examples/lipids.py",
    "examples/ln_9999.txt",
    "scripts/check_dist.py",
    "tests/test_det.py",
    "tests/test_general.py",
}
EXPECTED_WHEEL_FILES = {
    "assemblycfg/__init__.py",
    "assemblycfg/_version.py",
    "assemblycfg/cfg_ai.py",
    "assemblycfg/det.py",
    "assemblycfg/utils.py",
}


def fail(message: str) -> None:
    raise SystemExit(f"distribution check failed: {message}")


def one_artifact(dist_dir: Path, pattern: str) -> Path:
    matches = sorted(dist_dir.glob(pattern))
    if len(matches) != 1:
        fail(f"expected one {pattern!r} artifact in {dist_dir}, found {len(matches)}")
    return matches[0]


def normalized_name(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def citation_version(path: Path) -> str:
    match = re.search(
        r'^version:\s*["\']?([^"\'\s]+)["\']?\s*$',
        path.read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    if match is None:
        fail(f"could not read a version from {path}")
    return match.group(1)


def inspect_wheel(wheel: Path) -> str:
    with ZipFile(wheel) as archive:
        members = archive.namelist()
        metadata_files = [name for name in members if name.endswith(".dist-info/METADATA")]
        if len(metadata_files) != 1:
            fail(f"expected one METADATA file in {wheel.name}")

        metadata = Parser().parsestr(archive.read(metadata_files[0]).decode("utf-8"))
        if normalized_name(metadata["Name"]) != normalized_name(PROJECT_NAME):
            fail(f"wheel project name is {metadata['Name']!r}, expected {PROJECT_NAME!r}")

        version = metadata["Version"]
        expected_filename = f"{DIST_FILENAME_NAME}-{version}-py3-none-any.whl"
        if wheel.name != expected_filename:
            fail(f"wheel filename is {wheel.name!r}, expected {expected_filename!r}")

        provided_extras = set(metadata.get_all("Provides-Extra", []))
        missing_extras = EXPECTED_EXTRAS - provided_extras
        if missing_extras:
            fail(f"wheel is missing expected extras: {sorted(missing_extras)}")

        metadata_root = PurePosixPath(metadata_files[0]).parts[0]
        unexpected_roots = {
            PurePosixPath(name).parts[0]
            for name in members
            if name
            and PurePosixPath(name).parts[0] not in {PROJECT_NAME, metadata_root}
        }
        if unexpected_roots:
            fail(f"wheel contains unexpected top-level paths: {sorted(unexpected_roots)}")

        missing = EXPECTED_WHEEL_FILES - set(members)
        if missing:
            fail(f"wheel is missing expected files: {sorted(missing)}")

        return version


def inspect_sdist(sdist: Path, expected_version: str) -> str:
    expected_root = f"{DIST_FILENAME_NAME}-{expected_version}"
    expected_filename = f"{expected_root}.tar.gz"
    if sdist.name != expected_filename:
        fail(f"sdist filename is {sdist.name!r}, expected {expected_filename!r}")

    with tarfile.open(sdist, "r:gz") as archive:
        members = archive.getmembers()
        roots = {PurePosixPath(member.name).parts[0] for member in members if member.name}
        if roots != {expected_root}:
            fail(f"sdist root paths are {sorted(roots)}, expected only {expected_root!r}")

        pkg_info_path = f"{expected_root}/PKG-INFO"
        pkg_info_member = archive.getmember(pkg_info_path)
        pkg_info_file = archive.extractfile(pkg_info_member)
        if pkg_info_file is None:
            fail(f"could not read {pkg_info_path} from {sdist.name}")
        metadata = Parser().parsestr(pkg_info_file.read().decode("utf-8"))

        relative_members = {
            name.split("/", 1)[1]
            for member in members
            if (name := member.name).count("/") >= 1
        }

    if normalized_name(metadata["Name"]) != normalized_name(PROJECT_NAME):
        fail(f"sdist project name is {metadata['Name']!r}, expected {PROJECT_NAME!r}")
    if metadata["Version"] != expected_version:
        fail(f"sdist version is {metadata['Version']!r}, expected {expected_version!r}")

    missing = EXPECTED_SDIST_FILES - relative_members
    if missing:
        fail(f"sdist is missing expected files: {sorted(missing)}")

    return metadata["Version"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dist_dir", nargs="?", type=Path, default=Path("dist"))
    parser.add_argument("--tag", help="release tag to compare with the artifact version")
    args = parser.parse_args()

    wheel = one_artifact(args.dist_dir, "*.whl")
    sdist = one_artifact(args.dist_dir, "*.tar.gz")
    built_version = inspect_wheel(wheel)
    sdist_version = inspect_sdist(sdist, built_version)

    source_version = runpy.run_path("assemblycfg/_version.py")["__version__"]
    cited_version = citation_version(Path("CITATION.cff"))
    versions = {
        "wheel": built_version,
        "sdist": sdist_version,
        "source": source_version,
        "citation": cited_version,
    }
    if len(set(versions.values())) != 1:
        fail(f"versions disagree: {versions}")

    if args.tag is not None and args.tag.removeprefix("v") != built_version:
        fail(f"release tag {args.tag!r} does not match version {built_version!r}")

    print(f"validated {wheel.name} and {sdist.name} (version {built_version})")


if __name__ == "__main__":
    main()
