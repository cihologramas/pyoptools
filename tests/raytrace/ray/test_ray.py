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
        label=label,
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

    assert ray.Ray.almost_equal(ray1, ray2)


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
    assert ray.Ray.almost_equal(ray1, ray2)


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
    assert not ray.Ray.almost_equal(ray1, ray2)


def test_ray_dir():
    "The direction is always a unitary vector. This is done with the setter"
    ray1 = ray.Ray(dir=(0, 0, 5))
    assert np.array_equal(ray1.dir, np.array([0, 0, 1]))

    ray1.dir = (0, 5, -5)
    np.testing.assert_array_almost_equal(
        ray1.dir, [0, 0.70710678, -0.70710678])


def test_ray_pos():
    ray1 = ray.Ray()

    ray1.pos = (1, 2, 3)
    np.testing.assert_array_almost_equal(ray1.pos, [1, 2, 3])


def test_ch_coord_sys_inv():
    # ch_coord_sys_inv(self,no,ae,childs=False)
    origin = np.array([0, 0, 0])
    direction = np.array([0, 1, 0])

    ray_expected = ray.Ray(pos=(0.0, 0.0, 0.0),
                           dir=(0.84147098, 0.0, 0.54030231))

    ray1 = ray.Ray()
    ray_calculated = ray1.ch_coord_sys_inv(origin, direction)

    assert ray.Ray.almost_equal(ray_expected, ray_calculated)


# TODO this seems broken as it pass with any values of origin
def test_ch_coord_sys_inv_f():
    # ch_coord_sys_inv_f(self,np.ndarray no ,np.ndarray ae,bool childs)
    origin = np.array([0, 0, 0])
    direction = np.array([0, 1, 0])

    ray_expected = ray.Ray(pos=(0.0, 0.0, 0.0),
                           dir=(0.84147098, 0.0, 0.54030231))

    ray1 = ray.Ray()
    ray_calculated = ray1.ch_coord_sys_inv_f(origin, direction, False)
    print(ray_calculated)

    assert ray.Ray.almost_equal(ray_expected, ray_calculated)


# TODO this seems broken as it pass with many values of origin and direction
def test_ch_coord_sys():
    # ch_coord_sys(self, np.ndarray no, np.ndarray ae)
    origin = np.array([2, 1, 2])
    direction = np.array([1, 2, 1])

    ray_expected = ray.Ray(pos=(0.0, 0.0, 0.0), dir=(0, 0, 1))

    ray1 = ray.Ray()
    ray_calculated = ray1.ch_coord_sys(origin, direction)
    print(ray_calculated)

    assert ray.Ray.almost_equal(ray_expected, ray_calculated)


# TODO Please add a test here
def test_get_final_rays():
    # get_final_rays(self, inc_zeros=True)
    pass


def test_copy():
    ray1 = ray.Ray(
        pos=np.array([1.0, 0.0, 0.0]),
        dir=np.array([0.0, 0.0, 1.0]),
        intensity=1.0,
        wavelength=0.6328,
        n=None,
        label="HeNe",
        orig_surf=None,
        order=0,
    )

    ray2 = ray1.copy()

    assert ray2 == ray1


def test_reverse():
    ray_expected = ray.Ray(pos=(0.0, 0.0, 0.0), dir=(0.0, 0.0, -1.0))

    ray1 = ray.Ray()
    ray_calculated = ray1.reverse()

    assert ray.Ray.almost_equal(ray_expected, ray_calculated)


# TODO Please add a test here
def test_add_child():
    pass


# TODO Please add a test here
def test_optical_path_parent():
    pass


# TODO Please add a test here
def test_optical_path():
    pass
