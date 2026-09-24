# pyOpTools

**pyOpTools** is an open-source Python and Cython library for simulating optical systems using 3D non-sequential ray tracing, component modeling, and wavefront calculations.

Originally developed by the technological development team at [Combustión Ingenieros S.A.S.](http://www.cihologramas.com) and [Colombian Imaging Technologies S.A.S.](http://www.citech.com.co).

---

## Features

- **3D Ray Tracing**: Fast 3D non-sequential ray tracing with Eigen-accelerated coordinate transformations and intersection solvers.
- **Optical Component Library**: Spherical lenses, aspheric lenses, prisms, mirrors, apertures, stops, detectors, and compound optical systems.
- **Vendor Catalogs**: Built-in access to commercial optics catalogs (Thorlabs, Edmund Optics).
- **Wavefront Propagation & Diffraction**: Rayleigh-Sommerfeld, Angular Spectrum, and Fractional Fourier Transform methods.
- **Interactive Visualization**: 3D interactive system visualization in JupyterLab and Marimo using Plotly and WebGL.

---

## Installation

pyOpTools requires **Python 3.10+**.

### Using uv (Recommended)

```bash
# Add to an existing uv project:
uv add pyoptools

# Or install in your virtual environment:
uv pip install pyoptools
```

### Using pip

```bash
pip install pyoptools
```

---

## Quickstart

Here is a simple example showing how to build an optical system with a spherical lens and trace a ray:

```python
from pyoptools.raytrace.comp_lib import SphericalLens
from pyoptools.raytrace.ray import Ray
from pyoptools.raytrace.system import System

# 1. Define optical components
lens = SphericalLens(
    radius=25.0,
    thickness=5.0,
    curvature_s1=1.0 / 100.0,
    curvature_s2=-1.0 / 100.0,
    material=1.5168,
)

# 2. Assemble the optical system
system = System(complist=[(lens, (0, 0, 50), (0, 0, 0))], n=1.0)

# 3. Define an incident ray and propagate
ray = Ray(origin=(0, 10, 0), direction=(0, 0, 1), wavelength=0.589)
system.ray_add(ray)
system.propagate()

# 4. Inspect traced ray and component hits
final_rays = ray.get_final_rays()
print("Final ray direction:", final_rays[0].direction)
print("Impacts on lens:", len(lens.hit_list))
```

---

## Documentation

Comprehensive documentation, tutorials, and API reference are available at:
**[https://pyoptools.readthedocs.io/](https://pyoptools.readthedocs.io/)**

---

## Contributing & Development

We welcome contributions! Please see:
- **[CONTRIBUTING.md](CONTRIBUTING.md)** for developer setup, testing, and pull request guidelines.
- **[AGENTS.md](AGENTS.md)** for automated agent and internal architectural guidelines.

To set up a local development environment:

```bash
git clone https://github.com/cihologramas/pyoptools.git
cd pyoptools
uv sync
uv run python setup.py build_ext --inplace
uv run pytest
```
