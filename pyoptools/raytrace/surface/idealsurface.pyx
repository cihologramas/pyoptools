
"""
Modulo con clases y funciones auxiliares.
"""


from numpy import array, float64
from pyoptools.misc.definitions import inf_vect

from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray


cdef class IdealSurface(Surface):
    """Clase que representa una superficie ideal. Se utiliza para crear
    lentes ideales
    """

    cdef public double f

    def __init__(self, f=100, *args, **kwargs):
        """
        f represents the focal length of the ideal surface
        """
        Surface.__init__(self, *args, **kwargs)

        self.f = f

        # Add attributes to the state list
        self.addkey("f")

    cdef inline double topo_cy(self, double x, double y) noexcept nogil:
        return 0

    cpdef _intersection(self, Ray A):
        """Returns the intersection point between a ray and an the XY plane

        """
        # N_=array([0.,0.,1.])

        P1 = A.origin     # Punto que pertenece al rayo "Origen" del rayo
        L1 = A.direction  # Vector paralelo a la linea

        # if dot(N_,L1) ==0 : return inf_vect
        if L1[2] == 0:
            return inf_vect

        # print N_,P1,L1
        # print dot(N_,-P1),dot(N_,L1)
        # u=dot(N_,-P1)/dot(N_,L1)
        u = -P1[2]/L1[2]
        # Si u es muy grande, no hay intersección

        retval = P1+u*L1

        return retval

    cpdef normal(self, ri):
        """Method that returns the normal to the surface
        """
        N_ = array((0., 0., 1.)).astype(float64)
        return (N_)

    cpdef list propagate(self, Ray ri, double ni, double nr):

        PI, _P = self.int_nor(ri)
        _rx, _ry, rz = ri.direction
        # Get the focussing point as the point where the principal ray hits
        # the focal plane

        FP = ri.direction*self.f/abs(rz)

        ret = []

        if self.reflectivity != 1:
            # print "not 1"
            d = FP-PI
            if self.f < 0:
                d = -d
            ret.append(Ray(origin=PI, direction=d,
                           intensity=ri.intensity,
                           wavelength=ri.wavelength, n=ni, label=ri.label,
                           orig_surf=self.id))
        if self.reflectivity != 0:
            # print "not 0"
            ret.append(Ray(origin=PI, direction=PI-FP,
                           intensity=ri.intensity,
                           wavelength=ri.wavelength, n=ni, label=ri.label,
                           orig_surf=self.id))

        # print self.reflectivity
        return ret
