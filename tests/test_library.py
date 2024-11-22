import pytest
from pyoptools.raytrace.library import library
from pyoptools.raytrace._comp_lib.optic_factory import optic_factory
from pyoptools.raytrace._comp_lib.spherical_lens import SphericalLens
import tempfile


@pytest.mark.skip(reason="Skipping this test for now takes too long")
def test_all_optics():
    for i, (part, descriptor) in enumerate(library.items()):
        optic = optic_factory(**descriptor)


def test_direct_access():
    lens = library["LB1862"]
    assert isinstance(lens, SphericalLens)


def test_attribute_access():
    lens = library.Thorlabs["LB1862"]
    assert isinstance(lens, SphericalLens)

    lens = library.Edmund["08068"]
    assert isinstance(lens, SphericalLens)


def test_descriptor():
    assert library.Edmund.descriptor("08068")["radius"] == 12.2
    assert library.descriptor("08068")["radius"] == 12.2

    with pytest.raises(KeyError):
        library.Edmund.descriptor("None")

    with pytest.raises(KeyError):
        library.descriptor("None")


def test_parts():
    assert "LB1862" in library.Thorlabs.parts()
    assert "LB1862" in library.parts()


def test_user_lib():
    ntf = tempfile.NamedTemporaryFile("w", suffix=".json")
    with ntf as fp:
        fp.write(
            """
           {
                "a_test_lens": {
                    "description": "test",
                    "material": "BK7",
                    "glass_catalogs": "SCHOTT",
                    "thickness": 3.5,
                    "radius": 12.7,
                    "curvature_s2": 0.0,
                    "curvature_s1": 0.02,
                    "type": "SphericalLens"
                }
            }
        """
        )
        fp.flush()
        fp.seek(0)

        library.add(ntf.name)

        assert library.descriptor("a_test_lens")["thickness"] == 3.5
        assert library.descriptor("a_test_lens")["description"] == "test"
