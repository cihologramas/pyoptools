import os
import pytest
import pycodestyle  # This is the updated name for the pep8 tool
from subprocess import run, PIPE

import pyoptools

PEP8_ADDITIONAL_IGNORE = ["E501"]

EXCLUDE_FILES = []


EXCLUDE_CYTHON_FILES = []


@pytest.mark.linting
def test_pep8_conformance():
    dirname = os.path.dirname(pyoptools.__file__)
    pep8style = pycodestyle.StyleGuide()

    # Extend the number of PEP8 guidelines which are not checked.
    pep8style.options.ignore = pep8style.options.ignore + tuple(PEP8_ADDITIONAL_IGNORE)
    pep8style.options.exclude.extend(EXCLUDE_FILES)
    pep8style.options.filename = ["*.py"]

    result = pep8style.check_files([dirname])
    msg = "Found Python code syntax errors (and warnings)."
    assert result.total_errors == 0, msg


@pytest.mark.linting
def test_cython_conformance():
    dirname = os.path.dirname(pyoptools.__file__)

    for root, dirs, files in os.walk(dirname):
        for file in files:
            if (
                file.endswith(".pyx") or file.endswith(".pxd")
            ) and file not in EXCLUDE_CYTHON_FILES:
                filepath = os.path.join(root, file)
                result = run(
                    ["cython-lint", filepath], stdout=PIPE, stderr=PIPE, text=True
                )

                assert (
                    result.returncode == 0
                ), f"Found Cython code syntax errors in {file}:\n{result.stdout}{result.stderr}"
