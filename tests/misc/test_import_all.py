import pytest

def test_import_pyoptools():
    try:
        import pyoptools.all
        assert True
    except ImportError:
        assert False