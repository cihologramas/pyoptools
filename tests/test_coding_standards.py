import os
from subprocess import run

import pytest

import pyoptools

EXCLUDE_CYTHON_FILES = []


@pytest.mark.linting
def test_pep8_conformance():
    """Verify that all pure Python files adhere to PEP8/Ruff standards."""
    dirname = os.path.dirname(pyoptools.__file__)
    tests_dir = os.path.join(os.path.dirname(dirname), "tests")

    result = run(
        ["ruff", "check", dirname, tests_dir],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"Ruff found lint errors:\n{result.stdout}{result.stderr}"
    )


@pytest.mark.linting
def test_cython_conformance():
    """Verify that all Cython files adhere to cython-lint standards."""
    dirname = os.path.dirname(pyoptools.__file__)

    for root, dirs, files in os.walk(dirname):
        for file in files:
            if file.endswith((".pyx", ".pxd")) and file not in EXCLUDE_CYTHON_FILES:
                filepath = os.path.join(root, file)
                result = run(
                    ["cython-lint", filepath],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                assert result.returncode == 0, (
                    f"Found Cython code syntax errors in {file}:\n{result.stdout}{result.stderr}"
                )
