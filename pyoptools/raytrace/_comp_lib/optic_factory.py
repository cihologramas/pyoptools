from pyoptools.raytrace import comp_lib as cl
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.system import System
from pyoptools.raytrace.mat_lib import material

from inspect import signature


def optic_factory(**kwargs):
    """Returns a new optical component, given dictionary descriptor

    Descriptor should contain all the keyword arguments for a
    optical component constructor.

    Descriptor should contain a 'type' element with string for the
    optical component class name.

    Only elements with keyword arguments explicitly defined by the
    optical component class will be passed through to the constructor.
    """
    try:
        klass = getattr(cl, kwargs["type"])
        if not issubclass(klass, (Component, System)):
            raise ValueError(f"{kwargs['type']} not a valid optic type.")
    except AttributeError:
        raise ValueError(f"{kwargs['type']} not a valid optic type.")

    # Lookup and convert any material types:
    for k, v in kwargs.items():
        if "material" in k:
            if "glass_catalogs" in kwargs:
                m = material.get_from(v, kwargs["glass_catalogs"])
            else:
                m = material[v]
            kwargs[k] = m

    # Include only elements with corresponding key in constructor,
    # or 'material' (since it is defined in the base class)
    sig = signature(klass)
    valid_keys = list(sig.parameters.keys()) + ["material"]
    new_kwargs = {k: v for (k, v) in kwargs.items() if k in valid_keys}

    return klass(**new_kwargs)
