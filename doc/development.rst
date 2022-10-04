Developing Information
======================


Doc strings
-----------

The :doc:`old <development_old>` proposed way to document pyOpTools was confusing and should not be used any more. From now on all the documentation should be written in the `NumPy Style Python Docstrings <https://www.sphinx-doc.org/en/master/usage/extensions/example_numpy.html>`_. This style was chosen because it allows to easily document function and methods with multiple return values (return tuples).

See also: `<https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_

Creating :class:`Surface` Subclasses
------------------------------------

The following care must be taken when creating a :class:`Surface` subclass:

#. All surfaces super classes must be picklable. To do this follow:
    #. All subclasses inherit from pyoptools.picklable.picklable.Picklable
    #. For cython subclasses register all the attributes that affect the state using self.addkey("key")
        where key is changed by the attribute key.
    #. For python subclasses nothing needs to be done
        
