Developing Information
======================


Doc strings
-----------

The docstrings used to document the classes, methods and functions, should
contain the following information:

#. A paragraph containing a short description. (1 or 2 lines)
#. If needed a paragraph with a longer description
#. A description of the arguments received by the function, the method or
   the class constructor. If possible use a table as in the example.
#. For methods and functions description of the return value
#. A use example if possible
#. If further information is needed, a document should be attached to the
   end of the doc string using the include directive
The classes and the constructors should be documented in the class 
docstring. The __init__ docstring will not be used.

Example
>>>>>>>

The following is an example of how this should be written::

    def method( ..... ):
        '''
        Paragraph containing a short description of the method or function
        being documented
        
        A longer paragraph containing a more detailed description of the 
        method or function.
        
        **ARGUMENTS:**
            
            ==== ==============================
            arg1 Description of the argument 1
            arg2 Description of the argument 2.
                 It can be multiline.
            ==== ==============================
                        
        **RETURN VALUE:**
            
            Description of the return value. If the return value is a 
            tuple, it must show the tuple, and a table showing each of
            the values contained in the tuple.
            
        **EXAMPLE:**
            
            Code showing an example if possible. And a short description of it.
            
            
        .. include ../../../file.rst
        '''

Creating :class:`Surface` Subclasses
------------------------------------

The following care must be taken when creating a :class:`Surface` subclass:

#. All surfaces super classes must be picklable. To do this follow:
    #. All subclases inherit from pyoptools.picklable.picklable.Picklable
    #. For cython subclasses register all the attributes that affect the state using self.addkey("key")
        where key is changed by the attribute key.
    #. For python subclasses nothing needs to be done
        
#. D
#. D
#. D
