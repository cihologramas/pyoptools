from pyoptools.raytrace import comp_lib as cl
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.system import System
from pyoptools.raytrace.mat_lib import material

from inspect import signature

def optic_factory(**kwargs):
    """Returns a new optical component, given dictionary descriptor

    Descriptor should countain all the keyword arguments for a
    optical component constructor.

    Descriptor should contain a 'type' element with string for the
    optical component class name.

    Only elements with keyword arguments explicitely defined by the
    optical component class will be passed through to the constructor.
    """
    try:
        klass = getattr(cl, kwargs['type'])
        if not issubclass(klass, (Component, System)):
            raise ValueError(f"{kwargs['type']} not a valid optic type.")
    except AttributeError:
        raise ValueError(f"{kwargs['type']} not a valid optic type.")

    # Include only elements with corresponding key in constructor
    sig = signature(klass)
    new_kwargs = {k:v for (k,v) in kwargs.items() if k in sig.parameters}

    # Lookup and convert any material types:
    for k, v in new_kwargs.items():
        if 'material' in k:
            new_kwargs[k] = material[v]

    return klass(**new_kwargs)

