
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.component.component cimport Component
from pyoptools.raytrace.system.system cimport System


cdef class plist(dict):
    """
    Class representing a surface or component list.

    This class is designed to facilitate access to the internals of a system
    or component using a dictionary-like interface. It functions as a hybrid
    between a list and a dictionary, primarily operating as a dictionary but
    with an `append` method similar to a list. Each item must be a tuple `(O, P, D)`,
    where:

    - `O` is an object (limited to `Surfaces`, `Components`, or `Systems` in pyoptools).
    - `P` is a position vector (array, tuple, or list).
    - `D` is a direction vector (array, tuple, or list).

    When creating the list from a list or tuple, as no key is provided,
    the class automatically generates a key using `S#` for `Surfaces`
    or `C#` for `Components` and `Systems`, where `#` is an incrementing
    number. When created from a dictionary or another plist, the same
    keys from the original data are used.

    Parameters
    ----------
    items : list, tuple, dict, or plist, optional
        A collection of `(O, P, D)` tuples or `(O, P, D, Key)` tuples/lists. This can be
        provided as a list, tuple, dictionary, or plist.

    Notes
    -----
    This class uses the `dict` interface, allowing the use of string keys. It
    inherits from `dict` but adds additional functionality like automatic key
    generation and an `append` method.
    """

    def __init__(self, items=None):
        """
        Initialize the plist with optional items.

        Parameters
        ----------
        items : list, tuple, dict, or plist, optional
            A collection of items to initialize the plist. Each item should be in the
            form `(O, P, D)` or `(O, P, D, Key)`.

        Raises
        ------
        TypeError
            If any object in items is not an instance of `Surface`, `Component`,
            or `System`.
        ValueError
            If the items are not a list, tuple, dict, or plist, or if they do not
            have the correct format.
        """
        cdef tuple[double, double, double] position, rotation
        cdef str key

        super().__init__()  # Initialize the dictionary

        if isinstance(items, (list, tuple)):
            for item in items:
                if len(item) == 3:
                    obj, position, rotation = item
                    key = self._generate_key(obj)
                elif len(item) == 4:
                    obj, position, rotation, key = item
                else:
                    raise ValueError("Items must be a tuple or list of length 3 or 4")

                # The type check for obj is performed in self.__setitem__
                self[key] = (obj, position, rotation)

        elif isinstance(items, (dict, plist)):
            for key, value in items.items():
                obj, position, rotation = value

                # The type check for obj is performed in self.__setitem__
                self[key] = (obj, position, rotation)

        elif items is None:
            pass  # Do nothing if no items are provided

        else:
            raise ValueError("plist can be initialized only with tuples, "
                             "lists, dictionaries, or other plists")

    cdef str _generate_key(self, o):
        """
        Generate a unique key for the given object based on its type.

        Parameters
        ----------
        o : object
            The object for which the key is being generated. Must be an instance of
            `Surface`, `Component`, or `System`.

        Returns
        -------
        str
            A unique key string in the form 'S#' for `Surface` objects or 'C#' for
            `Component` or `System` objects.

        Raises
        ------
        ValueError
            If the object type is not recognized.
        """
        if isinstance(o, Surface):
            prefix = "S"
        elif isinstance(o, (Component, System)):
            prefix = "C"
        else:
            raise ValueError(f"Type {type(o).__name__} not recognized")

        # Find an unused key
        key_number = 0
        while prefix + str(key_number) in self:
            key_number += 1
        return prefix + str(key_number)

    def __setitem__(self, key, item):
        """
        Set an item in the dictionary.

        Parameters
        ----------
        key : str
            The key under which the item is stored. Must be a string.
        item : tuple
            A tuple containing the object, position, and rotation. The object
            must be an instance of `Surface`, `Component`, or `System`, and
            position and rotation must be tuples of doubles.

        Raises
        ------
        TypeError
            If the object is not an instance of `Surface`, `Component`, or
            `System`, or if the key is not a string.
        """

        cdef tuple[double, double, double] position, rotation

        # Unpack the item tuple
        obj, position, rotation = item

        # Validate the object type
        if not isinstance(obj, (Surface, Component, System)):
            raise TypeError(f"Object must be an instance of 'Surface', "
                            f"'Component', or 'System', got {type(obj).__name__}")

        # Validate the key type
        if not isinstance(key, str):
            raise TypeError("Key must be a string")

        # Assign the item to the dictionary
        super().__setitem__(key, (obj, position, rotation))

    def append(self, item):
        """
        Append an item to the plist, similar to how items are appended to a list.

        The key is auto-generated based on the type of the object.

        Parameters
        ----------
        item : tuple
            A tuple containing the object, position, and rotation. The object
            must be an instance of `Surface`, `Component`, or `System`.

        Raises
        ------
        TypeError
            If the object is not an instance of `Surface`, `Component`, or `System`.
        """

        cdef tuple[double, double, double] position, rotation
        cdef str key

        obj, position, rotation = item

        # Generate a key for the object
        key = self._generate_key(obj)

        # The type check for obj is performed in self.__setitem__
        self[key] = (obj, position, rotation)

    def __iter__(self):
        """
        Return an iterator over the values in the dictionary.

        Yields
        ------
        str
            The next key in the dictionary.
        """
        return iter(self.values())
