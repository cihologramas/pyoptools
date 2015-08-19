from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.misc.picklable.picklable cimport Picklable
cimport numpy as np
cdef class Surface(Picklable):
    cdef public object reflectivity
    cdef public Shape shape
    cdef public list _hit_list
    cdef public list id
    
    
    
    cpdef topo(self, x, y)
    cpdef int_nor(self,Ray iray)
    cpdef np.ndarray normal(self, int_p)
    cpdef np.ndarray intersection(self,Ray iray)
    cpdef _intersection(self,Ray iray)
    cpdef distance(self,Ray iray)
    cpdef distance_s(self,Ray iray)
    cpdef propagate(self,Ray ri,double ni,double nr)
    cpdef pw_propagate1(self, Ray ri,ni,nr, rsamples, isamples, knots)
    cpdef pw_propagate(self, Ray ri,ni,nr, rsamples, shape, order,z)
    cpdef pw_propagate_list(self, Ray ri,ni,nr, rsamples,z)
    cpdef wf_propagate(self, wf,ni,nr, samples,shape, knots)
    #def __repr__(self)
    cpdef reset(self)
    cpdef pw_cohef(self,ni,nr,ilimit, slimit, step, order, rsamples,zb)
