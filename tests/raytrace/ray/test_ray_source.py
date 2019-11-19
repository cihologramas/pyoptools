# standard imports

# third-party imports
import numpy as np

# local imports
import pyoptools.raytrace.ray.ray as ray
import pyoptools.raytrace.ray.ray_source as ray_source


def beam_equal(beam1, beam2):
    "Beams are equals i.e. each Ray in them are equal. Rays are compared in the same order in the two beams"
    are_equal = True
    for ray1, ray2 in zip(beam1, beam2):
        are_equal = are_equal and (ray1 == ray2)
    print(are_equal)
    return are_equal


def beam_almost_equal(beam1, beam2):
    "Beams are almost equals i.e. each Ray in them are equal. Rays are compared in the same order in the two beams"
    are_almost_equal = True
    for ray1, ray2 in zip(beam1, beam2):
        are_almost_equal = are_almost_equal and (ray1.almost_equal(ray2))
    print(are_almost_equal)
    return are_almost_equal


def test_beam_equal_same():
    "Ray must have the same .pos to be equal"
    wavelength = 0.1234
    label = "dummy_label"

    beam1 = [
        ray.Ray(
            pos=np.array([i, 0.0, 0.0]),
            dir=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=None,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-0.5, 0.5)
    ]

    beam2 = [
        ray.Ray(
            pos=np.array([i, 0.0, 0.0]),
            dir=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=None,
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
    "Ray must have the same .pos to be equal"
    wavelength = 0.1234
    label = "dummy_label"

    beam1 = [
        ray.Ray(
            pos=np.array([i, 0.0, 0.0]),
            dir=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=None,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-0.5, 0.5)
    ]

    beam2 = [
        ray.Ray(
            pos=np.array([i, 0.0, 0.0]),
            dir=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=None,
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
            pos=np.array([i, 0.0, 0.0000001]),
            dir=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=None,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-1, 1)
    ]

    beam2 = [
        ray.Ray(
            pos=np.array([i, 0.0, 0.0]),
            dir=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=None,
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
            pos=np.array([i, 0.0, 0.000001]),
            dir=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=None,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-1, 1)
    ]

    beam2 = [
        ray.Ray(
            pos=np.array([i, 0.0, 0.0]),
            dir=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=None,
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
    wavelength = 0.61234
    label = "test"

    expected = [
        ray.Ray(
            pos=np.array([i, j, 0.0]),
            dir=np.array([0.0, 0.0, 1.0]),
            intensity=1.0,
            wavelength=wavelength,
            n=None,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in (-0.5, 0.5)
        for j in (-0.5, 0.5)
    ]

    calculated = ray_source.parallel_beam_c(
        origin=(0.0, 0.0, 0.0),
        direction=(0.0, 0.0, 0.0),
        size=(1.0, 1.0),
        num_rays=(2, 2),
        wavelength=wavelength,
        label=label,
    )

    assert beam_equal(calculated, expected)


def test_parallel_beam_p():
    # def parallel_beam_p(origin=(0.,0.,0.),direction=(0.,0.,0),radius=0.5, num_rays=(5,10),wavelength=0.58929, label=""):
    wavelength = 0.61234
    label = "test13"
    pos_origin_ray = np.array([1.0, 1.5, 3.1416])

    expected_dir = np.array([0.67316743, -0.5721886, 0.46845044])
    expected_pos = [
        pos_origin_ray,
        np.array([1.52794548, 1.08696679, 1.87843843]),
        np.array([-0.05747045, 0.75593636, 3.75235817]),
        np.array([1.52952497, 2.65709685, 3.79400339]),
    ]
    expected = [
        ray.Ray(
            pos=i,
            dir=expected_dir,
            intensity=1.0,
            wavelength=wavelength,
            n=None,
            label=label,
            orig_surf=None,
            order=0,
        )
        for i in expected_pos
    ]

    calculated = ray_source.parallel_beam_p(
        origin=pos_origin_ray,
        direction=(5.2, 6.3, 10.3),
        radius=1.43,
        num_rays=(2, 3),
        wavelength=wavelength,
        label="test13",
    )

    assert beam_almost_equal(calculated, expected)