# standard imports

# third-party imports
import numpy as np

# local imports
import pyoptools.raytrace.ray.ray as ray
import pyoptools.raytrace.ray.ray_source as ray_source
from math import nan


def beam_equal(beam1, beam2):
    "Beams are equals i.e. each Ray in them are equal. Rays are compared in the same order in the two beams."
    are_equal = True
    for ray1, ray2 in zip(beam1, beam2):
        are_equal = are_equal and (ray1 == ray2)
    print(are_equal)
    return are_equal


def beam_almost_equal(beam1, beam2):
    "Beams are almost equals i.e. each Ray in them are equal. Rays are compared in the same order in the two beams."
    are_almost_equal = True
    for ray1, ray2 in zip(beam1, beam2):
        are_almost_equal = are_almost_equal and (ray.Ray.almost_equal(ray1, ray2))
    print(are_almost_equal)
    return are_almost_equal


def test_beam_equal_same():
    "Ray must have the same .origin to be equal"
    wavelength = 0.1234
    label = "dummy_label"

    beam1 = [
        ray.Ray(
            origin=np.array([i, 0.0, 0.0]),
            direction=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-0.5, 0.5)
    ]

    beam2 = [
        ray.Ray(
            origin=np.array([i, 0.0, 0.0]),
            direction=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-0.5, 0.5)
    ]

    print(beam1)
    print(beam2)
    assert beam_equal(beam1, beam2)


def test_beam_equal_different_pos():
    "Ray must have the same .origin to be equal"
    wavelength = 0.1234
    label = "dummy_label"

    beam1 = [
        ray.Ray(
            origin=np.array([i, 0.0, 0.0]),
            direction=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-0.5, 0.5)
    ]

    beam2 = [
        ray.Ray(
            origin=np.array([i, 0.0, 0.0]),
            direction=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-1, 1)
    ]

    assert not beam_equal(beam1, beam2)


def test_beam_almost_equal_same():
    "Ray almost equal are tested with beam_almost_equal"
    wavelength = 0.1234
    label = "dummy_label"

    beam1 = [
        ray.Ray(
            origin=np.array([i, 0.0, 0.0000001]),
            direction=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-1, 1)
    ]

    beam2 = [
        ray.Ray(
            origin=np.array([i, 0.0, 0.0]),
            direction=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-1, 1)
    ]

    print(beam1)
    print(beam2)
    assert beam_almost_equal(beam1, beam2)


def test_beam_less_almost_equal_not_same():
    "Ray almost equal are tested with beam_almost_equal"
    wavelength = 0.1234
    label = "dummy_label"

    beam1 = [
        ray.Ray(
            origin=np.array([i, 0.0, 0.000001]),
            direction=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-1, 1)
    ]

    beam2 = [
        ray.Ray(
            origin=np.array([i, 0.0, 0.0]),
            direction=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-1, 1)
    ]

    print(beam1)
    print(beam2)
    assert not beam_almost_equal(beam1, beam2)


def test_parallel_beam_c():
    # parallel_beam_c(origin=(0.,0.,0.),direction=(0.,0.,0.),size=(1.,1.),num_rays=(10,10),wavelength=0.58929, label="")
    origin = (0.0, 0.0, 0.0)
    direction = (0.0, 0.0, 0.0)
    size = (1.0, 1.0)
    num_rays = (2, 2)
    wavelength = 0.61234
    label = "test"

    expected = [
        ray.Ray(
            origin=np.array([i, j, 0.0]),
            direction=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-0.5, 0.5)
        for j in (-0.5, 0.5)
    ]

    calculated = ray_source.parallel_beam_c(
        origin=origin,
        direction=direction,
        size=size,
        num_rays=num_rays,
        wavelength=wavelength,
        label=label,
    )

    assert beam_equal(calculated, expected)


def test_parallel_beam_p():
    # def parallel_beam_p(origin=(0.,0.,0.),direction=(0.,0.,0),radius=0.5, num_rays=(5,10),wavelength=0.58929, label=""):
    pos_origin_ray = np.array([1.0, 1.5, 3.1416])
    direction = (5.2, 6.3, 10.3)
    radius = 1.43
    num_rays = (2, 3)
    wavelength = 0.61234
    label = "test13"

    expected_dir = np.array([0.67316743, -0.5721886, 0.46845044])
    expected_pos = [
        pos_origin_ray,
        np.array([1.52794548, 1.08696679, 1.87843843]),
        np.array([-0.05747045, 0.75593636, 3.75235817]),
        np.array([1.52952497, 2.65709685, 3.79400339]),
    ]
    expected = [
        ray.Ray(
            origin=i,
            direction=expected_dir,
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in expected_pos
    ]

    calculated = ray_source.parallel_beam_p(
        origin=pos_origin_ray,
        direction=direction,
        radius=radius,
        num_rays=num_rays,
        wavelength=wavelength,
        label=label,
    )

    assert beam_almost_equal(calculated, expected)


def test_point_source_c():
    # point_source_c(origin=(0.,0.,0.),direction=(0.,0.,0),span=(pi/8,pi/8)\
    #                      ,num_rays=(10,10),wavelength=0.58929, label="")
    pos_origin_ray = np.array([1.0, 1.5, 3.1416])
    direction = np.array([5.2, 6.3, 10.3])
    span = (np.pi / 8, np.pi / 8)
    num_rays = (2, 3)
    wavelength = 0.6328
    label = "another_test"

    expected_dir = [
        np.array([0.84217193, -0.45989021, 0.28150922]),
        np.array([0.73225862, -0.61754295, 0.28712021]),
        np.array([0.59697294, -0.75362932, 0.27507482]),
        np.array([0.69812011, -0.34719261, 0.6261674]),
        np.array([0.5882068, -0.50484536, 0.63177839]),
        np.array([0.45292112, -0.64093172, 0.619733]),
    ]
    expected = [
        ray.Ray(
            origin=pos_origin_ray,
            direction=i,
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in expected_dir
    ]

    calculated = ray_source.point_source_c(
        origin=pos_origin_ray,
        direction=direction,
        span=span,
        num_rays=num_rays,
        wavelength=wavelength,
        label=label,
    )

    assert beam_almost_equal(calculated, expected)


def test_point_source_p():
    # point_source_p(origin=(0., 0., 0.), direction=(0., 0., 0), span=pi / 8, num_rays=(10, 10), wavelength=0.58929,
    # label="")
    pos_origin_ray = np.array([1.0, 1.5, 3.1416])
    direction = np.array([0.2, 0.3, 0.5])
    span = np.pi / 16
    num_rays = (2, 3)
    wavelength = 0.432
    label = "one_more_test"

    expected_dir = [
        np.array([0.34942093, -0.03549297, 0.93629336]),
        np.array([0.38874339, -0.12238449, 0.91318159]),
        np.array([0.39840259, 0.04708771, 0.91600116]),
        np.array([0.25606913, -0.03066941, 0.96617182]),
    ]
    expected = [
        ray.Ray(
            origin=pos_origin_ray,
            direction=i,
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in expected_dir
    ]

    calculated = ray_source.point_source_p(
        origin=pos_origin_ray,
        direction=direction,
        span=span,
        num_rays=num_rays,
        wavelength=wavelength,
        label=label,
    )

    assert beam_almost_equal(calculated, expected)


def test_point_source_r():
    # point_source_r(origin=(0., 0., 0.), direction=(0., 0., 0), span=pi / 8, num_rays=100, wavelength=0.58929,
    #                            label=""):
    np.random.seed(0)

    pos_origin_ray = np.array([1.0, 1.234, 3.1416])
    direction = np.array([0.2, 0.3, 0.5])
    span = np.pi / 16
    num_rays = 5
    wavelength = 0.5460
    label = "test_with_random"

    expected_dir = [
        np.array([0.53157453, -0.30101208, 0.7917198]),
        np.array([0.74074805, -0.00967543, 0.67171326]),
        np.array([0.3209326, -0.43252757, 0.84256879]),
        np.array([0.39631705, -0.21297594, 0.89307001]),
        np.array([0.4072408, 0.01950802, 0.91311246]),
    ]
    expected = [
        ray.Ray(
            origin=pos_origin_ray,
            direction=i,
            intensity=1.0,
            wavelength=wavelength,
            n=nan,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in expected_dir
    ]

    calculated = ray_source.point_source_r(
        origin=pos_origin_ray,
        direction=direction,
        span=span,
        num_rays=num_rays,
        wavelength=wavelength,
        label=label,
    )

    assert beam_almost_equal(calculated, expected)
