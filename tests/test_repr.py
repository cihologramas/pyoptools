from pyoptools.raytrace.surface.cylindrical import Cylindrical
from pyoptools.raytrace.surface.spherical import Spherical


def test_spherical_repr():
    # Create a Spherical object with known parameters
    curvature = 0.01
    reflectivity = 0.5
    aperture_shape = None  # Assuming no aperture shape for simplicity
    spherical_surface = Spherical(
        curvature=curvature, reflectivity=reflectivity, shape=aperture_shape
    )
    # Expected representation string
    expected_repr = (
        f"Spherical(shape={repr(spherical_surface.shape)}, "
        f"reflectivity={reflectivity}, "
        f"curvature={curvature}, "
        f"id={spherical_surface.id})"
    )
    # Check if the actual representation matches the expected one
    assert repr(spherical_surface) == expected_repr


def test_cylindrical_repr():
    # Create a Cylindrical object with known parameters
    curvature = 0.01
    reflectivity = 0.5
    aperture_shape = None  # Assuming no aperture shape for simplicity
    cylindrical_surface = Cylindrical(
        curvature=curvature, reflectivity=reflectivity, shape=aperture_shape
    )
    # Expected representation string
    expected_repr = (
        f"Cylindrical(shape={repr(cylindrical_surface.shape)}, "
        f"reflectivity={reflectivity}, "
        f"curvature={curvature}, "
        f"id={cylindrical_surface.id})"
    )
    # Check if the actual representation matches the expected one
    assert repr(cylindrical_surface) == expected_repr
