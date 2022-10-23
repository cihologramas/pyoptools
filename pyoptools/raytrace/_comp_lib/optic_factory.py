
from pyoptools.raytrace import comp_lib as cl
from pyoptools.raytrace.component import Component
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
        if not issubclass(klass, Component):
            raise ValueError(f"{kwargs['type']} not a valid optic type.")
    except AttributeError:
        raise ValueError(f"{kwargs['type']} not a valid optic type.")

    sig = signature(klass)
    new_kwargs = {k:v for (k,v) in kwargs.items() if k in sig.parameters}

    return klass(**new_kwargs)

