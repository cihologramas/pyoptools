from pyoptools.misc.picklable.picklable cimport Picklable
cdef class Shape(Picklable):
    cpdef hit(self, p)
    cpdef bint fhit(self,double px,double py,double pz)
    #cpdef polylist(self, topo)
    cpdef pointlist(self)
    cpdef limits(self)
       
       
