
cdef class Picklable:
    cdef list __pkeys__
    cdef addkey(self, key)
