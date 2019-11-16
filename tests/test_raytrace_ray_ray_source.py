# standard imports

# third-party imports
import numpy as np

# local imports
import pyoptools.raytrace.ray.ray as ray
import pyoptools.raytrace.ray.ray_source as ray_source

#TODO Write a function to compare beams


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

    for calculated_ray, expected_ray in zip(calculated, expected):
        assert calculated_ray == expected_ray


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

    for calculated_ray, expected_ray in zip(calculated, expected):
        assert calculated_ray.almost_equal(expected_ray)