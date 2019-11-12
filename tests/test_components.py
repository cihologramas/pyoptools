# standard imports

# third-party imports

# local imports
from pyoptools.raytrace import shape
from pyoptools.raytrace import surface
from pyoptools.raytrace import component
from pyoptools.raytrace.mat_lib import material


def test_component_surflist():
    S0 = surface.Spherical(shape=shape.Circular(radius=50),
                           curvature=1. / 200.)
    S1 = surface.Spherical(shape=shape.Circular(radius=50),
                           curvature=1. / 200.)
    S2 = surface.Cylinder(radius=50, length=10)
    L1 = component.Component(surflist=[(S0, (0, 0, -5), (0, 0, 0)),
                                       (S1, (0, 0, 5), (0, 0, 0)),
                                       (S2, (0, 0, 6.5), (0, 0, 0))
                                       ],
                             material=material.schott["BK7"])
    print(type(L1.surflist))
    for surf in L1.surflist:
        assert len(surf) == 3, 'Wrong surface length'