# standard imports

# third-party imports
import numpy as np

# local imports
import pyoptools.raytrace.ray.ray as ray
import pyoptools.raytrace.ray.ray_source as ray_source


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
