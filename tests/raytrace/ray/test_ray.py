import pytest
import numpy as np
import pyoptools.raytrace.ray.ray as ray
from math import nan

def test_ray_equal():
    "Rays with exactly same attributes are equals."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    assert ray1 == ray2


def test_ray_almost_equal_pos():
    "Rays with different pos are different."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0000001]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    assert ray1 != ray2


def test_ray_almost_equal_dir():
    "Rays with different dir are different."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.00000001, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    assert ray1 != ray2


def test_ray_dir_different_norm():
    "Rays with same normalized dir are the same."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.01]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    assert ray1 == ray2


def test_ray_almost_equal_intensity():
    "Rays with different intensity are different."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0000001,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    assert ray1 != ray2


def test_ray_almost_equal_wavelength():
    "Rays with different wavelength are different."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength + 0.00000001,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    assert ray1 != ray2


def test_ray_different_n():
    "Rays with different n are different."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=1,
        label=label,
        orig_surf=None,
        order=0,
    )

    assert ray1 != ray2


def test_ray_different_label():
    "Rays with different label are different."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label + "_different",
        orig_surf=None,
        order=0,
    )

    assert ray1 != ray2


def test_ray_different_orig_surf():
    "Rays with different orig_surf are different."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=[],
        order=0,
    )

    assert ray1 != ray2


def test_ray_different_order():
    "Rays with different order are different."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([5.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=1,
    )

    assert ray1 != ray2


def test_ray_almost_equal_same():
    "Rays with exactly same attributes are equals."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([1.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([1.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    assert ray.Ray.almost_equal(ray1, ray2)


def test_ray_almost_equal_almost_same():
    "Rays with almost same attributes are almost equal."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([1.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([1.0, 0.0, 0.0000001]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    assert ray.Ray.almost_equal(ray1, ray2)


def test_ray_almost_equal_less_almost_same():
    "Rays with less almost same attributes are not almost equal."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        origin=np.array([1.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        origin=np.array([1.0, 0.0, 0.000001]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=nan,
        label=label,
        orig_surf=None,
        order=0,
    )

    assert not ray.Ray.almost_equal(ray1, ray2)


def test_ray_dir():
    "The direction is always a unitary vector. This is done with the setter."
    ray1 = ray.Ray(direction=(0, 0, 5))
    assert np.array_equal(ray1.direction, np.array([0, 0, 1]))

    ray1.direction = (0, 5, -5)
    np.testing.assert_array_almost_equal(ray1.direction, [0, 0.70710678, -0.70710678])


def test_ray_pos():
    ray1 = ray.Ray()

    ray1.origin = (1, 2, 3)
    np.testing.assert_array_almost_equal(ray1.origin, [1, 2, 3])


def test_ch_coord_sys_inv():
    origin = np.array([0, 0, 0])
    direction = np.array([0, 1, 0])

    ray_expected = ray.Ray(
        origin=(0.0, 0.0, 0.0), direction=(0.84147098, 0.0, 0.54030231)
    )

    ray1 = ray.Ray()
    ray_calculated = ray1.ch_coord_sys_inv(origin, direction)

    assert ray.Ray.almost_equal(ray_expected, ray_calculated)


def test_ch_coord_sys():
    origin = np.array([0, 0, 0])
    direction = np.array([0, 1, 0])

    ray_expected = ray.Ray()

    ray1 = ray.Ray(origin=(0.0, 0.0, 0.0), direction=(0.84147098, 0.0, 0.54030231))

    ray_calculated = ray1.ch_coord_sys(origin, direction)

    assert ray.Ray.almost_equal(ray_expected, ray_calculated)


@pytest.mark.skip(reason="Test for get_final_rays is pending.")
def test_get_final_rays():
    pass


def test_copy():
    ray1 = ray.Ray(
        origin=np.array([1.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=0.6328,
        n=nan,
        label="HeNe",
        orig_surf=None,
        order=0,
    )

    ray2 = ray1.copy()

    assert ray2 == ray1


def test_reverse():
    ray_expected = ray.Ray(origin=(0.0, 0.0, 0.0), direction=(0.0, 0.0, -1.0))

    ray1 = ray.Ray()
    ray_calculated = ray1.reverse()

    assert ray.Ray.almost_equal(ray_expected, ray_calculated)


@pytest.mark.skip(reason="Test for add_child is pending.")
def test_add_child():
    pass


@pytest.mark.skip(reason="Test for optical_path_parent is pending.")
def test_optical_path_parent():
    pass


@pytest.mark.skip(reason="Test for optical_path is pending.")
def test_optical_path():
    pass
