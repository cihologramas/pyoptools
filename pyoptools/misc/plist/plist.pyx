
from pyoptools.misc.picklable.picklable cimport Picklable
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.component.component cimport Component
from pyoptools.raytrace.system.system cimport System



cdef class plist(Picklable):
    """
    Class representing a part or component list.

    This class is designed to facilitate access to the internals of a System
    or a component using a dictionary-like interface. It functions as a
    hybrid between a list and a dictionary, primarily operating as a
    dictionary but with an `append` method similar to a list. Each item
    must be a tuple `(O, P, D)`, where:

    - `O` is an object (limited to `Surfaces`, `Components`, or `Systems`
    in pyoptools).
    - `P` is a position vector (array, tuple, or list).
    - `D` is a direction vector (array, tuple, or list).

    When creating the list from a list or tuple, as no key is provided,
    the class automatically generates a key using `S#` for `Surfaces`
    or `C#` for `Components` and `Systems`, where `#` is an incrementing
    number. When created from a dictionary or another plist, the same
    keys from the original data are used.

    Parameters
    ----------
    items : list, tuple, dict, or plist
        A collection of `(O, P, D)` tuples. This can be provided as a list,
        tuple, dictionary, or plist.

    Notes
    -----
    This class uses the `dict` interface, allowing the use of string keys. It 
    does not inherit from `dict` due to Cython's limitations with multiple
    inheritance.
    """

    def __init__(self, items=None):
        """
        Create the list to save the objects, positions and rotations
        """

        # Register class attributes on Picklable superclass
        Picklable.__init__(self, "_buf")

        self._buf={}
        cdef tuple[double,double,double] ap, ar, a , r

        if isinstance(items, (list, tuple)):
            for i in items:
                if len(i)==3:
                    O, p, r=i
                    ap = p
                    ar = r
                    self.append((O, ap, ar))
                else: # if len = 4
                    O, p, r, k=i

                    ap = p
                    ar = r
                    self[k]=(O, ap, ar)

        elif isinstance(items, (dict, plist)):
            for k, v in items.iteritems():
                O, p, r = v
                ap=p
                ar=r
                self._buf[k]=(O, ap, ar)
        elif items is None:
            pass
        else:
            raise ValueError, \
                "plist can be initialized only with tuples, lists or " + \
                "dictionaries"

    # Expose DICT API

    def __len__(self):
        return self._buf.__len__()

    def __getitem__(self , x):
        return self._buf[x]

    def __setitem__(self, key, val):

        cdef tuple[double,double,double] ap, ar
        O, p, r =val
        ap= p
        ar= r
        self._buf[key]=(O, ap, ar)

    def __delitem__(self, key):
        self._buf.__delitem__(key)

    def __contains__(self, key):
        return self._buf.__contains__(key)

    def __repr__(self):
        return repr(self._buf)

    # Return an iterator so this can be used similar to a list
    def __iter__(self):
        return iter(self._buf.values())

    def iteritems(self):
        return self._buf.iteritems()

    def iter(self):
        return self._buf.iter()

    def clear(self):
        return self._buf.clear()

    def items(self):
        return self._buf.items()

    def iterkeys(self):
        return self._buf.iterkeys()

    def itervalues(self):
        return self._buf.itervalues()

    def keys(self):
        return self._buf.keys()

    def pop(self, *argv, **argkw):
        return self._buf.pop(*argv, **argkw)

    def popitem(self):
        return self._buf.popitem()

    def setdefault(self, *argv, **argkw):
        return self._buf.setdefault(*argv, **argkw)

    def update(self, *argv, **argkw):
        return self._buf.update(*argv, **argkw)

    def values(self):
        return self._buf.values()

    def viewitems(self):
        return self._buf.viewitems()

    def viewkeys(self):
        return self._buf.viewkeys()

    def viewvalues(self):
        return self._buf.viewvalues()

    def append(self, a):
        '''
        Append an instance in a similar way as they are append in a list.

        The key is auto generated
        '''

        cdef tuple[double, double, double] ap,ar

        O, p, r = a
        ap=p
        ar=r

        if isinstance(O, Surface):
            pf="S"
        elif isinstance(O, (Component, System)):
            pf="C"
        else:
            pf="A"  # Auto

        # Find an unused key
        k=0
        while pf+str(k) in self._buf:
            k=k+1

        self._buf[pf+str(k)]=(O, ap, ar)

