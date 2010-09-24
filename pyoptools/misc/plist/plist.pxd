from pyoptools.misc.picklable.picklable cimport Picklable
cdef class plist(Picklable):
    cdef public dict _buf #Buffer used to save the information
    #~ cdef int _cnt  #counter used when plist is used as iterator
