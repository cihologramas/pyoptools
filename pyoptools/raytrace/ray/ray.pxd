
cimport numpy as np
from cpython cimport bool


cdef class Ray:
    cdef np.ndarray cpos 
    cdef public object  n, parent, label
    cdef public int order # In surfaces that produce multiple rays, each ray
                          # should be assigned a different order. This is done
                          # in the ray.add_child method
                          
    cdef public list orig_surf # path of the originating surface
    cdef public double intensity, wavelength, pop
    cdef list __childs, 
    cdef np.ndarray _dir
    
    cpdef Ray ch_coord_sys(self, np.ndarray no,np. ndarray ae)
    cpdef Ray ch_coord_sys_inv_f(self,np.ndarray no ,np.ndarray ae,bool childs)
    
    #cpdef Ray ch_coord_sys(self, no, ae)
    

    #~ def ch_coord_sys_inv(self,no,ae,childs=False):
    #~ def ch_coord_sys(self,no,ae):
    #~ def get_final_rays(self, inc_zeros=True):
    #~ def copy(self):
    #~ def reverse(self):
    #~ def add_child(self, cr):
    #~ def optical_path_parent(self):
    #~ def optical_path(self):
cdef Ray Rayf(np.ndarray pos,np.ndarray dir,double intensity,double wavelength,
                n,label,parent,double pop,orig_surf,int order)
