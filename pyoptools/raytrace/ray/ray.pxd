from cpython cimport bool

from pyoptools.misc.cmisc.eigen cimport Vector3d, Matrix3d

from pyoptools.raytrace.surface.surface cimport Surface


cdef class Ray:
    # Se eliminó n de acá, hay que verificar por que era un objeto de python y
    # no un double
    cdef public object label, draw_color
    cdef public Ray parent

    # In surfaces that produce multiple rays, each ray should be assigned a
    # different order. This is done in the ray.add_child method

    cdef public int order

    cdef public list orig_surf  # path of the originating surface
    cdef public double intensity, wavelength, pop
    cdef list __childs

    cdef Vector3d _origin  # _origin is the backing attribute for the origin property
    cdef Vector3d _direction  # Baking of the dir property
    cdef public double n

    # amount of parents of this ray in the optical path. It counts the number of
    # times  a ray have been propagated through surfaces. It is used to stop
    # the propagation in resonant cavities.

    cdef public int _parent_cnt

    cdef Ray ch_coord_sys_inv_f(self, Vector3d& origin_coordinates ,
                                Vector3d& rotation_angles, bool childs)
    cdef Ray ch_coord_sys_f(self, Vector3d& origin_coordinates ,
                            Vector3d& rotation_angles)

    @staticmethod
    cdef Ray fast_init(Vector3d& origin, Vector3d& direction,
                       double intensity, double wavelength, double n,
                       object label, object draw_color, Ray parent,
                       double pop, list orig_surf,
                       int order, int parent_cnt)
