import warnings

from libc.math cimport isnan, NAN, INFINITY
from libc.string cimport memcpy

'''Ray class definition module
'''

cdef extern from "math.h":
    double sqrt(double)
cimport cython
from pyoptools.raytrace.surface.surface cimport Surface


#from pyoptools.misc.cmisc.cmisc cimport norm_vect, empty_vec, \
#     dot_product_3x3_matrix_vector, to_vector, rot_mat, rot_mat_i, \
#     norm_3d_vector, allclose_cython

from pyoptools.misc.cmisc.linalg cimport Vector3, Matrix3x3, \
     vector3_from_python_object, vector3_to_tuple, normalize_vector3, \
     matrix3x3_vector3_dot, compute_rotation_matrix, compute_rotation_matrix_i, \
     add_vector3, substract_vector3, vector3_magnitude, vector3_equals


cdef class Ray:
    """Class to define a ray

    **ARGUMENTS**

    =========== =======================================================
    pos         Tuple (x,y,z) containing the origin of the Ray
    dir         Tuple (x,y,z)  containing the direction vector of the Ray
    intensity   Floating point number representing the Intensity of the
                Ray.
                Warning: Check how can a physically correct definition
                can be made
    wavelength  Wavelength (in vacuum) of the ray in micrometers
                (.58929 by default)
    n           Refraction index of the point originating the ray.
                If the value is None, the ray was emitted from the media
                and its Refraction index is taken (not from inside a
                component)
    label       String used to follow the rays through the system.
    draw_color  Color used to render this ray. If None, the wavelength
                will be used to determine the color. Otherwise, can be
                any valid matplotlib color identifier.
    parent      Ray where this ray comes from (used to follow ray
                trajectory).
    childs      List of rays originated by the interaction of this ray
                with an optical surface.
    =========== =======================================================
    """

    def __init__(self, pos, dir, double intensity=1.,
                 double wavelength=.58929, double n=NAN, label="", draw_color=None,
                 parent=None, double pop=0., orig_surf=None, order=0, parent_cnt=0):

        # These are written through the properties, so any python object is valid
        self.origin = pos
        self.direction = dir

        self.intensity=intensity
        self.wavelength=wavelength

        #TODO: Check in all pyoptools where Ray.n is checked against None. It should be NaN
        self.n=n
        self.label=label
        self.draw_color=draw_color
        self.parent=parent
        self.pop=pop
        self.orig_surf=orig_surf
        self.__childs=[]
        self.order=order
        self._parent_cnt = parent_cnt

    @staticmethod

    cdef Ray fast_init(Vector3* origin_ptr, Vector3* direction_ptr, double intensity, double wavelength,
                double n, object label, object draw_color, Ray parent, double pop, list orig_surf,
                int order, int parent_cnt):
        """
        Function to create and initialize Ray instances fast from Cython.
        When changing the code, always check that the __init__ and the
        fast_init initialization do the same.

        Parameters
        ----------
        origin_ptr : Vector3*
            Pointer to a Vector3 struct representing the origin of the Ray.
            Must be a 3D point in space (x, y, z).
        direction_ptr : Vector3*
            Pointer to a Vector3 struct representing the direction vector of the Ray.
            Must be a 3D direction vector (x, y, z) with a magnitude of 1.
            Note: The magnitude is not checked for efficiency; ensure it is normalized before passing.
        intensity : double
            Floating point number representing the intensity of the Ray.
            Warning: Ensure a physically correct definition of intensity.
        wavelength : double
            Wavelength (in vacuum) of the ray in micrometers. Default is 0.58929 µm.
        n : double
            Refraction index of the point originating the ray.
            If None, the ray is emitted from the media and takes the refraction index
            of the surrounding environment (not from inside a component).
        label : object
            String or any identifier used to track the rays through the system.
        draw_color : object
            Color used to render this ray. If None, the wavelength is used to determine
            the color; otherwise, any valid matplotlib color identifier can be used.
        parent : object
            Ray instance or object representing the parent Ray from which this Ray originates.
        pop : double
            Represents a property related to the Ray's population.
        orig_surf : object
            Object representing the original optical surface related to the Ray.
        order : int
            Integer representing the order of the Ray in the system.
        parent_cnt : int
            Count of parent Rays for this Ray, used to manage complex ray trajectories.

        Returns
        -------
        Ray
            A new Ray instance with all properties initialized according to the provided arguments.

        Notes
        -----
        - This function directly initializes Ray attributes for fast construction, bypassing the
        usual Python __init__ method for performance reasons.
        - Be cautious when modifying the initialization parameters or the corresponding
        __init__ method to ensure consistency.
        - The direction vector's magnitude must be 1. However, for efficiency, this is not checked
        in the function. Ensure the direction is normalized before passing it to this function.
        """

        cdef Ray instance = Ray.__new__(Ray)

        # Efficiently copy the Vector3 structs using memcpy
        # _direction is not normalized here
        memcpy(&instance._origin, origin_ptr, sizeof(Vector3))
        memcpy(&instance._direction, direction_ptr, sizeof(Vector3))

        # Set the remaining properties
        instance.intensity = intensity
        instance.wavelength = wavelength
        instance.n = n
        instance.label = label
        instance.draw_color = draw_color
        instance.parent = parent
        instance.pop = pop
        instance.orig_surf = orig_surf
        instance.__childs = []
        instance.order = order
        instance._parent_cnt = parent_cnt

        return instance

    def __reduce__(self):
        args=(self.pos, self.dir, self.intensity, self.wavelength, self.n ,
              self.label, self.draw_color, self.parent, self.pop,
              self.orig_surf, self.order)
        return(type(self), args, self.__getstate__())

    def __getstate__(self):
        return self.__childs

    def __setstate__(self, state):
        self.__childs=state

    @property
    def childs(self):
        """
        Get a tuple of child Rays.

        This property returns a tuple containing all child Ray objects. The internal list
        is converted to a tuple to prevent accidental modifications (e.g., using `append`),
        ensuring that the child list cannot be modified directly through this property.

        Returns
        -------
        tuple of Ray
            A tuple containing all child Ray objects.
        """
        # The list is converted to a tuple to prevent accidental modification
        return tuple(self.__childs)

    @property
    def direction(self):
        """
        Get the direction vector of the ray.

        The direction vector is normalized to have a length of 1. However, due to
        floating-point rounding errors, the vector may occasionally deviate slightly
        from being perfectly normalized.

        Returns
        -------
        tuple
            A tuple (X, Y, Z) with the components of the direction vector of the ray.
        """
        return vector3_to_tuple( &self._direction)

    @direction.setter
    def direction(self, dir):
        """
        Set the direction vector of the ray, ensuring it is normalized.

        The input vector will be normalized to have a length of 1. Due to
        floating-point rounding errors, the magnitude of the input vector
        (`|dir|`) may occasionally be slightly different than 1.

        Parameters
        ----------
        dir : array-like
            The input direction vector and normalized to ensure its magnitude is 1.
        """
        # Convert input to memory view
        cdef Vector3 tmp_dir = vector3_from_python_object(dir)
        normalize_vector3(&tmp_dir)
        self._direction = tmp_dir


    @property
    def dir(self):
        warnings.warn(
        "The 'dir' attribute is deprecated and will be removed in a future"
        " version. Please use the 'direction' attribute instead.",
        DeprecationWarning,
        stacklevel=2)
        return self.direction

    @dir.setter
    def dir(self, dir):
        warnings.warn(
        "The 'dir' attribute is deprecated and will be removed in a future"
        " version. Please use the 'direction' attribute instead.",
        DeprecationWarning,
        stacklevel=2)
        self.direction = dir


    # Use property origin to access the field _origin from a python program.
    # It can receive lists, tuples or arrays
    # the cython only visible _origin, can only receive arrays
    @property
    def origin(self):
        return vector3_to_tuple(&self._origin)

    @origin.setter
    def origin(self, origin_coordinates):
        cdef int i
        self._origin = vector3_from_python_object(origin_coordinates)

    @property
    def pos(self):
        warnings.warn(
        "The 'pos' attribute is deprecated and will be removed in a future"
        " version. Please use the 'origin' attribute instead.",
        DeprecationWarning,
        stacklevel=2)
        return self.origin

    @pos.setter
    def pos(self, pos):
        warnings.warn(
        "The 'pos' attribute is deprecated and will be removed in a future"
        " version. Please use the 'origin' attribute instead.",
        DeprecationWarning,
        stacklevel=2)
        self.origin = pos


    def ch_coord_sys_inv(self, origin_coordinates, rotation_angles, bool childs=False):
        """
        Transform the coordinate system of the Ray.

        This method applies a transformation to the coordinate system of the Ray, 
        consisting of a rotation followed by a translation. The rotation is applied 
        in the order of Z, Y, and then X axes. Optionally, the transformation can 
        also be applied to any child objects.

        Parameters
        ----------
        origin_coordinates : array-like of float
            A list-like object (such as a tuple, list, or NumPy array) containing
            three elements (X, Y, Z) that represent the coordinates of the origin
            of the old coordinate system in the new coordinate system.
            
        rotation_angles : array-like of float
            A list-like object (such as a tuple, list, or NumPy array) containing
            three elements (RX, RY, RZ) that represent the rotation angles 
            in radians to be applied to the old coordinate system. The rotations 
            are applied in the order of RZ (around the Z-axis), then RY 
            (around the Y-axis), and finally RX (around the X-axis).

        childs : bool, optional
            If True, the coordinate system transformation is also applied to child
            objects. Default is False.

        Returns
        -------
        Ray
            A transformed Ray object with the new coordinate system applied.

        Notes
        -----
        The transformation consists of a rotation applied first, followed by a 
        translation to the new origin.

        Examples
        --------
        >>> ray = Ray()
        >>> new_ray = ray.ch_coord_sys_inv((0, 0, 0), (90, 0, 0), childs=True)
        >>> print(new_ray)
        """
        cdef Vector3 origin_coordinates_vector = \
            vector3_from_python_object(origin_coordinates)
        cdef Vector3 rotation_angles_vector = \
            vector3_from_python_object(rotation_angles)

        cdef Ray parent =  self.ch_coord_sys_inv_f(&origin_coordinates_vector,
                                                   &rotation_angles_vector,
                                                   childs)
        return parent
    
    cdef Ray ch_coord_sys_inv_f(self, Vector3 *origin_coordinates_ptr ,
                                Vector3 *rotation_angles_ptr, bool childs):
        """
        Transform the coordinate system of the Ray.

        Fast version to be used in Cython.

        Parameters
        ----------
        origin_coordinates : Vector3 representing the coordinates of the origin of the old
            coordinate system in the new coordinate system.

        rotation_angles : Vector3 representing the rotation angles to be applied to
            the old coordinate system.

        childs : bool
            If True, the coordinate system transformation is also applied to
            child objects.

        Returns
        -------
        Ray
            A transformed Ray object with the new coordinate system applied.
        """

        cdef Matrix3x3 tm
        compute_rotation_matrix(rotation_angles_ptr, &tm)

        # Calculate new ray origin by applying rotation and translation
        cdef Vector3 rotated_ray_origin
        matrix3x3_vector3_dot(&tm, &self._origin, &rotated_ray_origin)

        cdef Vector3 new_ray_origin
        add_vector3(&rotated_ray_origin, origin_coordinates_ptr, &new_ray_origin)


        #Calculate new ray direction
        cdef Vector3 new_ray_direction
        matrix3x3_vector3_dot(&tm, &self._direction, &new_ray_direction)

        cdef Ray parent = Ray.fast_init(&new_ray_origin,
                                        &new_ray_direction,
                                        self.intensity,
                                        self.wavelength,
                                        self.n,
                                        self.label,
                                        self.draw_color,
                                        None ,
                                        0,
                                        self.orig_surf,
                                        self.order,
                                        self._parent_cnt)

        # Calculate the transform of the childs and link them

        cdef Ray i

        if childs:
            for i in self.childs:
                it=i.ch_coord_sys_inv_f(origin_coordinates_ptr, rotation_angles_ptr,
                                        childs)
                parent.add_child(it)
        return parent

    def ch_coord_sys(self, origin_coordinates, rotation_angles):

        cdef Vector3 origin_coordinates_vector = \
            vector3_from_python_object(origin_coordinates)
        cdef Vector3 rotation_angles_vector = \
            vector3_from_python_object(rotation_angles)

        cdef Ray parent =  self.ch_coord_sys_f(&origin_coordinates_vector,
                                               &rotation_angles_vector)
        return parent

    cdef Ray ch_coord_sys_f(self, Vector3 *origin_coordinates_ptr, 
                            Vector3 *rotation_angles_ptr):
        
        cdef Matrix3x3 tm
        compute_rotation_matrix_i(rotation_angles_ptr, &tm)

        cdef Vector3 new_ray_origin, new_ray_direction, translated_ray_origin

        substract_vector3(&self._origin, origin_coordinates_ptr,&translated_ray_origin)

        matrix3x3_vector3_dot(&tm, &translated_ray_origin, &new_ray_origin)
        matrix3x3_vector3_dot(&tm, &self._direction, &new_ray_direction)

        return Ray.fast_init(&new_ray_origin,
                             &new_ray_direction,
                             self.intensity,
                             self.wavelength,
                             self.n,
                             self.label,
                             self.draw_color,
                             None ,
                             0,
                             self.orig_surf,
                             self.order, 
                             self._parent_cnt)

    def get_final_rays(self, inc_zeros=True):
        '''Find the final rays of the raytrace

        *inc_zeros*
            If inc_zeros == True, all the child rays are included.
            If set to false, the rays with intensity==0 are not
            included
        '''

        retval=[]

        if len(self.childs)==0:
            if(inc_zeros):
                retval=[self]
            elif self.intensity!=0:
                retval=[self]
            else:
                if len(self.parent.childs)==1:
                    retval=[self.parent]
        else:
            for i in self.childs:
                retval=retval+i.get_final_rays(inc_zeros)
        return retval

    def copy(self):
        '''Return a copy ray leaving the parent=None, and childs=[], and
        order=0
        '''
        return Ray(pos=self.pos, dir=self.dir, intensity=self.intensity,
                   wavelength=self.wavelength, n=self.n, label=self.label)

    def reverse(self):
        '''Return a copy ray leaving the parent=None, and childs=[], and
        order=0, and inverting the ray direction.
        '''
        return Ray(pos=self.pos, dir=-self.dir, intensity=self.intensity,
                   wavelength=self.wavelength, n=self.n, label=self.label,
                   orig_surf=self.orig_surf)

    def __repr__(self):
        return "Ray(" + \
                "," + repr(self.origin) + \
                "," + repr(self.dir) + \
                ",intensity=" + repr(self.intensity) + \
                ",wavelength=" + repr(self.wavelength) + \
                ",n = " + repr(self.n) + \
                ",label=" + repr(self.label) + \
                ",orig_surf=" + repr(self.orig_surf) + \
                ",order=" + repr(self.order) + ")"

    def add_child(self, cr):
        '''Add childs to the current ray, and create the appropriate links

        *cr*
            Ray to include in the child list
        '''

        # A child with intensity==0 is used to indicate a parent end point, for
        # example the intersection point with a stop, so it must be in the
        # list as a child.
        assert(isinstance(cr, Ray)), "A ray child must be a ray"
        cr.parent=self
        cr.order=len(self.__childs)
        # cr._parent_cnt = self._parent_cnt + 1
        self.__childs.append(cr)

    def optical_path_parent(self):
        ''' Return the optical path from the origin of the origin ray to the
       end of this ray parent (this ray origin)
        '''
        cdef Vector3 length_vector

        if self.parent is not None:
            if self.pop!=0:
                print("The pop attribute of the ray has a value of ", self.pop,
                      " instead the real parent optical path is being used")
            #path= norm_3d_vector(self.origin-self.parent.origin)*self.parent.n
            substract_vector3(&self._origin, &(self.parent._origin), &length_vector)
            path =  vector3_magnitude(&length_vector)*self.parent.n
            return path+self.parent.optical_path_parent()

        return self.pop

    def optical_path(self):
        ''' Return the optical path of the beam propagation from the origin of
        the origin ray, to the end of this ray (intersection with a surface)
        '''

        if self.intensity==0:
            return 0.
        elif len(self.childs)==0:
            return INFINITY
        else:
            return self.childs[0].optical_path_parent()

    def __eq__(self, other):

        cdef Ray other_ray = <Ray> other

        return (vector3_equals(&self._origin, &(other_ray._origin)) and
                vector3_equals(&self._direction, &(other_ray._direction)) and
                self.intensity == other.intensity and
                self.wavelength == other.wavelength and
                self.n == other.n and
                self.label == other.label and
                self.draw_color == other.draw_color and
                self.order == other.order and
                self.orig_surf == other.orig_surf)
        # TODO do we have to compare self.pop and other.pop ?
        # self.copy() indicate that we should actually only compare
        # fields pos, dir, intensity, wavelength, n and label

    @staticmethod
    def almost_equal(Ray ray1, Ray ray2, double tol=1e-7):
        return (vector3_equals(&(ray1._origin), &(ray2._origin), tol) and
                vector3_equals(&(ray1._direction), &(ray2._direction), tol) and
                ray1.intensity == ray2.intensity and
                ray1.wavelength == ray2.wavelength and
                ray1.n == ray2.n and
                ray1.label == ray2.label and
                ray1.draw_color == ray2.draw_color and
                ray1.order == ray2.order and
                ray1.orig_surf == ray2.orig_surf)
