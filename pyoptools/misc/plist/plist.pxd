from pyoptools.misc.picklable.picklable cimport Picklable


cdef class plist(Picklable):
    cdef public dict _buf  # Buffer used to save the information

