# distutils: language = c++
from pyoptools.misc.plist.plist cimport plist
from pyoptools.misc.picklable.picklable cimport Picklable
from pyoptools.raytrace.mat_lib import Material
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.misc.cmisc.eigen cimport (
    Vector3d,
    Matrix3d,
    compute_rotation_matrix,
    assign_to_vector3d,
    convert_vector3d_to_tuple,
)
from libc.math cimport INFINITY

__all__ = ["Component"]

cdef class Component(Picklable):
    '''
    Class used to define an optical component

    To create a component the the user must provide a list containing
    the surfaces bounding the component, the vertex coordinates of each
    surface, and the rotation angles describing its orientation. It must
    also provide information about the material used to build the component.

    The material information of the component is given in the "material"
    attribute. This information can be a **Material** class instance (some
    default materials are defined in the **mat_lib module**), or a
    refraction index for materials with constant refraction index.

    **ARGUMENTS:**

    ========  ==========================================================
    surflist  A list of tuples of the form *(surface, (PosX,PosY,PosZ),
              (RotX,RotY,RotZ))* where surface is an instance of the
              :class:surface class, *PosX,PosY,PosZ* are the surface's
              vertex coordinates, and *RotX,RotY,RotZ* are the rotation
              angles of the surface around the X , Y , and Z axes, in
              radians. The rotation about the Z axis if applied first,
              then the rotation about the Y axis, and finally the
              rotation about the X axis.
    material  Instance of the class Material with the material
              definition, or a floating point number to indicate a
              constant refraction index material
    ========  ==========================================================
    '''
    property surflist:
        def __get__(self):
            return self._surflist

        def __set__(self, list):
            self._surflist = plist(list)

    property material:
        def __get__(self):
            return self._material

        def __set__(self, material):
            assert isinstance(material, (float, Material)), \
                "material must be a floating point number or a Material instance"
            self._material = material

    property hit_list:
        """Return a tuple of (hit_point, incident_ray) for all hits on the component's
        surfaces, transformed to the Component's coordinate system.
        """
        def __get__(self):
            cdef list ret_list = []
            cdef Matrix3d tm
            cdef Vector3d sr_vec, sc_vec, pi_vec, pi_c
            cdef Surface S
            cdef tuple SC, SR

            for surf in self.surflist:
                S, SC, SR = surf
                HL = S.hit_list
                if not HL:
                    continue
                assign_to_vector3d(SR, sr_vec)
                compute_rotation_matrix(sr_vec, tm)
                assign_to_vector3d(SC, sc_vec)
                for item in HL:
                    PI, R = item
                    assign_to_vector3d(PI, pi_vec)
                    pi_c = tm * pi_vec + sc_vec
                    ret_list.append((convert_vector3d_to_tuple(pi_c), R))
            return tuple(ret_list)

    def __init__(self, surflist=None, material=1.):

        # surflist defaults to None to avoid mutable default argument issues.

        if surflist is None:
            self.surflist = []
        else:
            self.surflist = surflist

        self.material = material

        Picklable.__init__(self, "_surflist", "_material")

    def __repr__(self):
        '''Return an string with the representation of the optical component

        It must be overloaded in all subclasses
        '''

        retval = "OptComp(\n"
        for s in self.surflist:
            retval = retval+str(s)+"\n"

        retval = retval+")\n"

        return retval

    # Dict type and list type interface to expose surflist

    def __len__(self):
        return self._surflist.__len__()

    def __getitem__(self, x):
        return self._surflist[x]

    def __setitem__(self, key, val):
        self._surflist[key] = val

    def __delitem__(self, key):
        self._surflist.__delitem__(key)

    def __contains__(self, key):
        return self._surflist.__contains__(key)

    # Return an iterator so this can be used similar to a list
    def __iter__(self):
        return iter(self._surflist.values())

    def iteritems(self):
        return iter(self._surflist.items())

    def iter(self):
        return iter(self._surflist.values())

    def clear(self):
        return self._surflist.clear()

    def items(self):
        return self._surflist.items()

    def iterkeys(self):
        return iter(self._surflist.keys())

    def itervalues(self):
        return iter(self._surflist.values())

    def keys(self):
        return self._surflist.keys()

    def pop(self, *argv, **argkw):
        return self._surflist.pop(*argv, **argkw)

    def popitem(self):
        return self._surflist.popitem()

    def setdefault(self, *argv, **argkw):
        return self._surflist.setdefault(*argv, **argkw)

    def update(self, *argv, **argkw):
        return self._surflist.update(*argv, **argkw)

    def values(self):
        return self._surflist.values()

    def viewitems(self):
        return self._surflist.items()

    def viewkeys(self):
        return self._surflist.keys()

    def viewvalues(self):
        return self._surflist.values()

    def get_surf_paths(self):
        """
        Method that returns a list that contains the path for each surface.

        A path here is a list containing the keys needed to read each surface.

        This method is an auxiliary method so this works when called from a System
        """
        l = []
        keys = self.keys()
        for k in keys:
            l.append([k])
        return l

    def n(self, wavelength=0.58929):
        """
        Refraction index of the component at the specified wavelength

        The wavelength should be given in um. If the wavelength is not given, it
        is assumed wavelength=0.58929 um
        """

        if isinstance(self.material, Material):
            return self.material.n(wavelength)

        return self.material

    def surf_changed(self):
        """
        Increases changes, to indicate that any of the surfaces used to build
        the optical component was modified
        :param self:
        :return:
        """

        self.changes = self.changes+1

    cpdef distance(self, Ray ri_):
        """
        Distance length from a ray origin to a component, following the ray path.

        Method that calculates the distance traveled by a ray from its origin to
        the next surface of the component. It returns the physical distance, not
        the optical distance


        *Return value*
            A tuple with the distance, the point of intersection using the
            coordinate system of the surface, and a pointer to the surface
            that is closest to the ray (distance,point of intersection, surface)
        """

        cdef tuple[double, double, double] P, D
        cdef double min_d = INFINITY
        cdef double d
        cdef object min_pi = None
        cdef Surface min_surf = None
        cdef Surface S
        cdef Ray R

        for surf in self.surflist:
            S, P, D = surf
            # Change the coordinate system of the ray, From the Component
            # coordinate system to the surface component system

            R = ri_.ch_coord_sys(P, D)

            Dist = S.distance(R)
            d = Dist[0]

            if d < min_d or min_surf is None:
                min_d = d
                min_pi = Dist[1]
                min_surf = S

        return min_d, min_pi, min_surf

    def reset(self):
        """Reset the optical component

        Method that reset the optical surfaces that compose the Component.
        For example the detector surfaces should be reset before repeating the
        ray trace to erase the hit lists. Normally this method should not be
        used directly. It is called when the reset method of the system class is
        called
        """
        cdef Surface S
        for comp in self.surflist:
            S, _P, _D = comp
            S.reset()

    def propagate(self, Ray ri, n_m):
        """Returns the next ray in the propagation

        Taking into account the interaction (refraction, reflection, ...)
        at the Component surfaces, return a list containing the resulting rays.
        n_m is the refraction index of the media surrounding the component.
        """
        # Determine the refraction index of the incident and the refracted media
        # If ri.n is equal to the component refraction index, the ray is coming out
        # from the component, if not, it is goint in to the component.

        my_n = self.n(ri.wavelength)

        if ri.n == my_n:
            n_p = n_m
            n = my_n
        else:
            n = n_m
            n_p = my_n

        # Search for the next surface to be hit
        cdef double min_d = INFINITY
        cdef double d
        cdef object min_surf_item = None
        cdef Surface S
        cdef Ray R

        for i in self.surflist:
            S, P, D = i
            # Change the coordinate system of the ray, From the Component
            # coordinate system to the surface component system, and calculate
            # the distance to the next surface.

            R = ri.ch_coord_sys(P, D)
            d = S.distance(R)[0]
            if d < min_d or min_surf_item is None:
                min_d = d
                min_surf_item = i

        # Find the closest surface, and change the ray to its coordinate system
        # and calculate the refraction
        SR, PSR, DSR = min_surf_item
        R = ri.ch_coord_sys(PSR, DSR)
        ri_n = SR.propagate(R, n, n_p)

        ret_rays = []
        for i in ri_n:
            rr = i.ch_coord_sys_inv(PSR, DSR)
            ret_rays.append(rr)

        return ret_rays
