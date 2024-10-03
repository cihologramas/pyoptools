"""
Modulo con clases y funciones auxiliares.
"""

from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.cmisc.eigen cimport Vector3d, assign_nan_to_vector3d, \
    convert_vector3d_to_tuple

from libc.math cimport abs

cdef class IdealPPlanes(Surface):
    """Clase que representa un par de superficies principales ideales. Se
    utiliza para crear lentes ideales gruesas
    """
    cdef double f
    cdef double d

    def __init__(self, f=100, d=20, *args, **kwargs):
        """
        f: Focal length
        d: Distance between planes
        """
        Surface.__init__(self, *args, **kwargs)

        self.f = f
        self.d = d

        # Add attributes to the state list
        # self.addkey("f")

    cdef inline double topo_cy(self, double x, double y) noexcept nogil:
        return 0

    cdef void _calculate_intersection(self,
                                      Ray incident_ray,
                                      Vector3d& intersection_point) noexcept nogil:
        """Returns the intersection point between a ray and an the XY plane

        """
        # N_=array([0.,0.,1.])

        cdef Vector3d P1 = incident_ray._origin
        cdef Vector3d L1 = incident_ray._direction

        # if dot(N_,L1) ==0 : return inf_vect
        if L1(2) == 0:
            assign_nan_to_vector3d(intersection_point)

        # print N_,P1,L1
        # print dot(N_,-P1),dot(N_,L1)
        # u=dot(N_,-P1)/dot(N_,L1)

        cdef double u, u1, u2
        u1 = -(P1(2)-self.d/2.)/L1(2)
        u2 = -(P1(2)+self.d/2.)/L1(2)

        if abs(u1) < abs(u2):
            if u1 != 0:
                u = u1
            else:
                u = u2
        else:
            if u2 != 0:
                u = u2
            else:
                u = u1

        intersection_point = incident_ray._direction * u+ incident_ray._origin

    cdef void _calculate_normal(self, Vector3d& intersection_point,
                                Vector3d& normal) noexcept nogil:
        """Method that returns the normal to the surface
        """
        (<double*>(&normal(0)))[0] = 0.
        (<double*>(&normal(1)))[0] = 0.
        (<double*>(&normal(2)))[0] = 1.

    cpdef list propagate(self, Ray incident_ray, double ni, double nr):
        """
        This is overloaded, as the ray propagation here is different than the
        snell law
        """

        cdef Vector3d intersection_point

        self.intersection_cy(incident_ray, intersection_point)

        cdef double rz = incident_ray._direction(2)

        # Get the focussing point as the point where the principal ray hits
        # the focal plane

        cdef Vector3d FP = incident_ray._direction*self.f/abs(rz)
        cdef Vector3d d

        # Las ecuaciones de refraccion usadas funcionan para el caso donde el plano
        # esta en Z=0. Como aca el plano no esta en z=0, hay que moverlo para el
        # calculo, y luego moverlo nuevamente a las posiciones adecuadas
        if intersection_point(2) < 0:
            (<double*>(&intersection_point(2)))[0] = 0
            d = FP-intersection_point
            (<double*>(&intersection_point(2)))[0] = self.d/2.
        else:
            (<double*>(&intersection_point(2)))[0] = 0
            d = FP-intersection_point
            (<double*>(&intersection_point(2)))[0] = -self.d/2.
        ret = []
        cdef Vector3d temp_vect
        if self.reflectivity != 1:
            ret.append(Ray(origin=convert_vector3d_to_tuple(intersection_point),
                       direction=convert_vector3d_to_tuple(d),
                           intensity=incident_ray.intensity,
                           wavelength=incident_ray.wavelength, n=ni,
                           label=incident_ray.label,
                           orig_surf=self.id))
        if self.reflectivity != 0:
            # print "not 0"
            temp_vect = intersection_point-FP
            ret.append(Ray(origin=convert_vector3d_to_tuple(intersection_point),
                           direction=convert_vector3d_to_tuple(temp_vect),
                           intensity=incident_ray.intensity,
                           wavelength=incident_ray.wavelength,
                           n=ni, label=incident_ray.label,
                           orig_surf=self.id))

        # print self.reflectivity
        return ret
