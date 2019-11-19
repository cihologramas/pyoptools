# standard imports

# third-party imports
import numpy as np

# local imports
import pyoptools.raytrace.ray.ray as ray


def test_ray_equal():
    "Rays with exactly same attributes are equals."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert ray1 == ray2


def test_ray_almost_equal_pos():
    "Rays with differents pos are differents."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0000001]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert not (ray1 == ray2)


def test_ray_almost_equal_dir():
    "Rays with different dir are differents."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.00000001, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert not (ray1 == ray2)


def test_ray_dir_different_norm():
    "Rays with same normalized dir are same."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.01]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert ray1 == ray2


def test_ray_almost_equal_intensity():
    "Rays with different intensity are differents."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0000001,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert not (ray1 == ray2)


def test_ray_almost_equal_wavelength():
    "Rays with different wavelength are differents."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength + 0.00000001,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert not (ray1 == ray2)


def test_ray_different_n():
    "Rays with differents n are differents."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=1,
        label=label,
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert not (ray1 == ray2)


def test_ray_different_label():
    "Rays with differents label are differents."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label + "_different",
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert not (ray1 == ray2)


def test_ray_different_orig_surf():
    "Rays with differents orig_surf are differents."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label ,
        orig_surf=[],
        order=0,
    )

    print(ray1)
    print(ray2)
    assert not (ray1 == ray2)


def test_ray_different_order():
    "Rays with different order are differents."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([5.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=1,
    )

    print(ray1)
    print(ray2)
    assert not (ray1 == ray2)


def test_ray_almost_equal_same():
    "Rays with exactly same attributes are equals."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([1.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([1.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert ray1.almost_equal(ray2)


def test_ray_almost_equal_almost_same():
    "Rays with exactly same attributes are equals."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([1.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([1.0, 0.0, 0.0000001]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert ray1.almost_equal(ray2)


def test_ray_almost_equal_less_almost_same():
    "Rays with exactly same attributes are equals."
    wavelength = 0.1234
    label = "dummy_label"

    ray1 = ray.Ray(
        pos=np.array([1.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    ray2 = ray.Ray(
        pos=np.array([1.0, 0.0, 0.000001]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=wavelength,
        n=None,
        label=label,
        orig_surf=None,
        order=0,
    )

    print(ray1)
    print(ray2)
    assert not ray1.almost_equal(ray2)