# standard imports

# third-party imports

# local imports
from pyoptools.raytrace.shape.circular import Circular


def test_circular():
    c = Circular()
    assert c.radius == 1.0, 'Unexpected radius'