# ------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amézquita Orozco
# Description:     Plane surface definition module
# Symbols Defined: Plane
# ------------------------------------------------------------------------------

import warnings

from libc.math cimport isnan, NAN, INFINITY

cdef extern from "math.h":
    double sqrt(double)

import warnings


from pyoptools.misc.cmisc.eigen cimport Vector3d, Matrix3d, convert_vector3d_to_tuple, \
    assign_to_vector3d, compute_rotation_matrix, compute_rotation_matrix_i, is_approx, \
    convert_vector3d_to_array

cdef class Ray:
    """
    Class to define a ray in an optical system.

    Parameters
    ----------
    origin : tuple of float
        A tuple (x, y, z) representing the origin of the ray.
    direction : tuple of float
        A tuple (x, y, z) representing the direction vector of the ray.
    intensity : float, optional
        The intensity of the ray. Defaults to 1.0.
        Note: The physical meaning of this intensity is currently not
        well-defined.
    wavelength : float, optional
        The wavelength of the ray in vacuum, in micrometers.
        Defaults to 0.58929 μm.
    n : float, optional
        The refractive index of the medium from which the ray originates.
        If `n` is `NaN`, the ray is assumed to have been emitted from
        the medium, and the medium's refractive index is used.
        Defaults to `NaN`.
    label : str or object, optional
        A string or identifier used to track the ray through the optical
        system. Defaults to an empty string.
    draw_color : object, optional
        The color used to render this ray. If `None`, the wavelength is
        used to determine the color. Otherwise, any valid Matplotlib
        color identifier can be used. Defaults to `None`.
    parent : Ray or None, optional
        The parent ray from which this ray originates, used to follow the
        ray's trajectory. Defaults to `None`.
    pop : float, optional
        The optical path length of the parent ray. This may be removed
        in future versions. Defaults to 0.0.
        TODO: Check if this may be removed

    orig_surf : list, optional
        A list representing the originating surface of the ray. Defaults to
        `None`.
        TODO: Explain what path means
    order : int, optional
        The index of the ray in the parent's child ray list. Defaults to 0.
    parent_cnt : int, optional
        The count of parent rays for this ray, used to manage complex
        ray trajectories. Defaults to 0.
    """
    def __init__(self, origin = (0, 0, 0), direction = (0, 0, 1),
                 double intensity=1., double wavelength=.58929, object n=NAN, label="",
                 draw_color=None, Ray parent=None, double pop=0., list orig_surf=None,
                 order=0, parent_cnt=0, **kwargs):

        # These are written through the properties, so any python object is valid

        if "pos" in kwargs:
            self.origin = kwargs["pos"]
            warnings.warn(
                "The 'pos' keyword argument is deprecated and will be removed in "
                "a future version. Please use 'origin' instead.",
                DeprecationWarning,
                stacklevel=2)
        else:
            self.origin = origin

        if "dir" in kwargs:
            self.direction = kwargs["dir"]
            warnings.warn(
                "The 'dir' keyword argument is deprecated and will be removed in "
                "a future version. Please use 'direction' instead.",
                DeprecationWarning,
                stacklevel=2)
        else:
            self.direction = direction

        self.intensity=intensity
        self.wavelength=wavelength

        if n is None:
            warnings.warn(
                "The refraction index given as None is deprecated and will be "
                "removed in a future version. Please use 'NaN' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            self.n = NAN
        else:
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
    cdef Ray fast_init(Vector3d& origin,
                       Vector3d& direction,
                       double intensity,
                       double wavelength,
                       double n,
                       object label,
                       object draw_color,
                       Ray parent,
                       double pop,
                       list orig_surf,
                       int order,
                       int parent_cnt):
        """
        Function to create and initialize Ray instances fast from Cython.
        When changing the code, always check that the __init__ and the
        fast_init initialization do the same.

        Parameters
        ----------
        origin : Vector3d&
            An eigen::Vector3d instance representing the origin of the Ray.
            Must be a 3D point in space (x, y, z).
        direction_ptr : Vector3&
            An eigen::Vector3d instance representing the direction vector of
            the Ray. Must be a 3D direction vector (x, y, z) with a magnitude
            of 1.
            Note: The magnitude is not checked for efficiency; ensure it is
            normalized before passing.
        intensity : double
            Floating point number representing the intensity of the Ray.
            Warning: Ensure a physically correct definition of intensity.
        wavelength : double
            Wavelength (in vacuum) of the ray in micrometers. Default is
            0.58929 µm.
        n : double
            Refraction index of the point originating the ray.
            If None, the ray is emitted from the media and takes the refraction
            index of the surrounding environment (not from inside a component).
        label : object
            String or any identifier used to track the rays through the system.
        draw_color : object
            Color used to render this ray. If None, the wavelength is used to
            determine the color; otherwise, any valid matplotlib color
            identifier can be used.
            TODO: This should be removed, is up to the plotting engine to
            decide the color.
        parent : Ray
            Ray instance representing the parent Ray from which this Ray
            originates.
        pop : double
            The optical path length of the parent ray. This may be removed
            in future versions. Defaults to 0.0.
            TODO: Check if this may be removed
        orig_surf : list
            A list representing the originating surface of the ray. Defaults
            to `None`.
            TODO: Explain what path means
        order : int
            The index of the ray in the parent's child ray list. Defaults to 0.
            parent_cnt : int

        Returns
        -------
        Ray
            A new Ray instance with all properties initialized according to the
            provided arguments.

        Notes
        -----
        - This function directly initializes Ray attributes for fast
        construction, bypassing the usual Python __init__ method for
        performance reasons.
        - Be cautious when modifying the initialization parameters or the
        corresponding __init__ method to ensure consistency.
        - The direction vector's magnitude must be 1. However, for efficiency,
        this is not checked in the function. Ensure the direction is
        normalized before passing it to this function.
        """

        cdef Ray instance = Ray.__new__(Ray)

        # _direction is not normalized here
        instance._origin = origin
        instance._direction = direction

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
        args=(self.origin, self.direction, self.intensity, self.wavelength, self.n ,
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

        This property returns a tuple containing all child Ray objects. The
        internal list is converted to a tuple to prevent accidental modifications
        (e.g., using `append`), ensuring that the child list cannot be modified
        directly through this property.

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
        array
            A numpy array (X, Y, Z) with the components of the direction vector
            of the ray.
        """
        return convert_vector3d_to_array(self._direction)

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
        cdef Vector3d tmp_dir
        assign_to_vector3d(dir, tmp_dir)
        self._direction = tmp_dir.normalized()

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
        """
        Get the origin of the ray as a tuple.

        This property returns the origin of the ray in 3D space, represented as
        a tuple (x, y, z). The origin is stored internally as an
        Eigen::Vector3d and is converted to a Python tuple when accessed
        through this property.

        Returns
        -------
        numpy array
            A numpy array containing the (x, y, z) coordinates of the ray's origin.
        """
        return convert_vector3d_to_array(self._origin)

    @origin.setter
    def origin(self, origin_coordinates):
        """
        Set the origin of the ray.

        This setter method assigns the given coordinates to the ray's origin.
        The input should be a Python object that supports indexing (such as a
        list, tuple, or NumPy array) with at least three elements representing
        the (x, y, z) coordinates. The internal origin is stored as an
        Eigen::Vector3d, and this method updates it with the provided values.

        Parameters
        ----------
        origin_coordinates : object
            A Python object that supports indexing (e.g., list, tuple, or
            NumPy array) and contains at least three elements representing
            the (x, y, z) coordinates of the new origin.

        Raises
        ------
        ValueError
            If `origin_coordinates` does not have at least three elements or
            if the elements cannot be converted to floating-point numbers.
        """
        assign_to_vector3d(origin_coordinates, self._origin)

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

        This method applies a transformation to the coordinate system of
        the Ray, consisting of a rotation followed by a translation. The
        rotation is applied in the order of the X, Y, and then Z axes. The
        rotation matrix is calculated as Rm = Rz * Ry * Rx. Optionally,
        the transformation can also be applied to any child objects.

        Parameters
        ----------
        origin_coordinates : array-like of float
            A list-like object (e.g., tuple, list, or NumPy array) containing
            three elements (X, Y, Z) that represent the coordinates of the
            origin of the old coordinate system in the new coordinate system.

        rotation_angles : array-like of float
            A list-like object (e.g., tuple, list, or NumPy array) containing
            three elements (Rx, Ry, Rz) that represent the rotation angles in
            radians to be applied to the old coordinate system. The rotations
            are applied in the order of Rz (around the Z-axis), then Ry (around
            the Y-axis), and finally Rx (around the X-axis).

        childs : bool, optional
            If True, the coordinate system transformation is also applied to
            child objects. Default is False.

        Returns
        -------
        Ray
            A transformed Ray object with the new coordinate system applied.

        Notes
        -----
        The transformation consists of a rotation applied first, followed
        by a translation to the new origin.
        """

        cdef Vector3d origin_coordinates_vector
        assign_to_vector3d(origin_coordinates, origin_coordinates_vector)
        cdef Vector3d rotation_angles_vector
        assign_to_vector3d(rotation_angles, rotation_angles_vector)

        cdef Ray parent = self.ch_coord_sys_inv_f(origin_coordinates_vector,
                                                  rotation_angles_vector,
                                                  childs)
        return parent

    cdef Ray ch_coord_sys_inv_f(self, Vector3d& origin_coordinates ,
                                Vector3d& rotation_angles, bool childs):
        """
        Transform the coordinate system of the Ray.

        This is a fast version intended for use in Cython. For a detailed
        explanation of the transformation process, see the `ch_coord_sys_inv`
        method.
        """

        cdef Matrix3d tm
        compute_rotation_matrix(rotation_angles, tm)

        # Calculate new ray origin by applying rotation and translation
        cdef Vector3d rotated_ray_origin = tm*self._origin
        # matrix3x3_vector3_dot(&tm, &self._origin, &rotated_ray_origin)

        cdef Vector3d new_ray_origin = rotated_ray_origin + origin_coordinates
        # add_vector3(&rotated_ray_origin, origin_coordinates_ptr, &new_ray_origin)

        # Calculate new ray direction
        cdef Vector3d new_ray_direction = tm*self._direction
        # matrix3x3_vector3_dot(&tm, &self._direction, &new_ray_direction)

        cdef Ray parent = Ray.fast_init(new_ray_origin,
                                        new_ray_direction,
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
                it=i.ch_coord_sys_inv_f(origin_coordinates, rotation_angles,
                                        childs)
                parent.add_child(it)
        return parent

    def ch_coord_sys(self, origin_coordinates, rotation_angles):
        """
        Transform the coordinate system of the Ray.

        This method applies a transformation to the coordinate system of
        the Ray, consisting of translation followed by a rotation. The
        rotation is applied in the order of the Z, Y, and then X axes. The
        rotation matrix is calculated as Rm = Rx * Ry * Rz.

        Parameters
        ----------
        origin_coordinates : array-like of float
            A list-like object (e.g., tuple, list, or NumPy array) containing
            three elements (X, Y, Z) that represent the coordinates of the
            origin of the new coordinate system in the old coordinate system.

        rotation_angles : array-like of float
            A list-like object (e.g., tuple, list, or NumPy array) containing
            three elements (Rx, Ry, Rz) that represent the rotation angles in
            radians to be applied to the old coordinate system. The rotations
            are applied in the order of Rz (around the Z-axis), then Ry (around
            the Y-axis), and finally Rx (around the X-axis).

        Returns
        -------
        Ray
            A transformed Ray object with the new coordinate system applied.

        Notes
        -----
        The transformation consists of a translation to the new origin first,
        followed by a rotation.
        """

        cdef Vector3d origin_coordinates_vector
        assign_to_vector3d(origin_coordinates, origin_coordinates_vector)
        cdef Vector3d rotation_angles_vector
        assign_to_vector3d(rotation_angles, rotation_angles_vector)

        cdef Ray parent = self.ch_coord_sys_f(origin_coordinates_vector,
                                              rotation_angles_vector)
        return parent

    cdef Ray ch_coord_sys_f(self, Vector3d& origin_coordinates,
                            Vector3d& rotation_angles):
        """
        Transform the coordinate system of the Ray.

        This is a fast version intended for use in Cython. For a detailed
        explanation of the transformation process, see the `ch_coord_sys`
        method.
        """

        cdef Matrix3d tm
        compute_rotation_matrix_i(rotation_angles, tm)

        cdef Vector3d new_ray_origin, new_ray_direction, translated_ray_origin

        # substract_vector3(&self._origin, origin_coordinates_ptr,
        # &translated_ray_origin)
        translated_ray_origin = self._origin - origin_coordinates

        # matrix3x3_vector3_dot(&tm, &translated_ray_origin, &new_ray_origin)
        new_ray_origin = tm*translated_ray_origin

        # matrix3x3_vector3_dot(&tm, &self._direction, &new_ray_direction)
        new_ray_direction = tm*self._direction

        return Ray.fast_init(new_ray_origin,
                             new_ray_direction,
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
        """
        Find the final rays of the raytrace.

        Parameters
        ----------
        inc_zeros : bool
            If True, all child rays are included in the result.
            If False, rays with `intensity == 0` are excluded.

        Returns
        -------
        list of Ray
            A list of final rays from the raytrace, optionally excluding rays
            with zero intensity.
        """
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
        """
        Return a copy of the Ray with specific attributes reset.

        This method returns a copy of the current Ray object with the `parent`
        attribute set to `None`, the `childs` attribute set to an empty list,
        and the `order` attribute set to 0.
        """

        return Ray(origin=self.origin, direction=self.direction,
                   intensity=self.intensity, wavelength=self.wavelength,
                   n=self.n, label=self.label)

    def reverse(self):
        """
        Return a copy of the Ray with specific attributes reset and the
        direction inverted.

        This method returns a copy of the current Ray object with the `parent`
        attribute set to `None`, the `childs` attribute set to an empty list,
        the `order` attribute set to 0, and the ray direction inverted.
        """
        cdef Vector3d new_dir = - self._direction
        return Ray(origin=self.origin,
                   direction=convert_vector3d_to_tuple(new_dir),
                   intensity=self.intensity,
                   wavelength=self.wavelength,
                   n=self.n,
                   label=self.label,
                   orig_surf=self.orig_surf)

    def __repr__(self):
        """
        Return a string representation of the Ray object.

        The string representation includes the origin, direction, intensity,
        wavelength, refractive index (n), label, originating surface, and order
        of the Ray.
        """
        return ("Ray(" + repr(self.origin) +
                ", " + repr(self.direction) +
                ", intensity=" + repr(self.intensity) +
                ", wavelength=" + repr(self.wavelength) +
                ", n=" + repr(self.n) +
                ", label=" + repr(self.label) +
                ", orig_surf=" + repr(self.orig_surf) +
                ", order=" + repr(self.order) + ")")

    def add_child(self, cr):
        """
        Add a child Ray to the current Ray and create the appropriate links.

        This method adds a Ray object to the child list of the current Ray.
        The child Ray's `parent` attribute is set to the current Ray, and
        its `order` attribute is set based on its position in the child list.

        Parameters
        ----------
        cr : Ray
            The Ray object to include in the child list.

        Notes
        -----
        A child Ray with `intensity == 0` is used to indicate a parent
        endpoint, such as an intersection point with a stop, so it must be
        included in the child list even if its intensity is zero.
        """

        # A child with intensity == 0 is used to indicate a parent endpoint,
        # for example the intersection point with a stop, so it must be in the
        # list as a child.
        assert isinstance(cr, Ray), "A ray child must be a Ray instance"
        cr.parent = self
        cr.order = len(self.__childs)
        self.__childs.append(cr)

    def optical_path_parent(self):
        """
        Return the optical path length from the origin of the original ray to
        the end of this ray's parent (i.e., the origin of this ray).

        This method calculates the optical path length recursively, starting
        from the origin of the first ray in the sequence and ending at the
        origin of the current ray. If the `pop` attribute of the current ray
        is non-zero, a warning is printed, and the actual parent optical path
        is used instead.

        Returns
        -------
        float
            The total optical path length from the original ray's origin to
            this ray's origin.

        Notes
        -----
        The optical path length is calculated as the product of the distance
        between two points and the refractive index (`n`) of the medium.
        If the ray has a parent, this calculation is performed recursively.
        """
        cdef Vector3d length_vector

        if self.parent is not None:
            if self.pop!=0:
                print("The pop attribute of the ray has a value of ", self.pop,
                      " instead the real parent optical path is being used")

            length_vector = self._origin - self.parent._origin
            path = length_vector.norm()*self.parent.n
            return path+self.parent.optical_path_parent()

        return self.pop

    def optical_path(self):
        """
        Return the optical path of the beam propagation from the origin of the
        original ray to the end of this ray (intersection with a surface).

        This method calculates the total optical path length from the origin
        of the original ray to the intersection point of this ray with a
        surface.

        Returns
        -------
        float
            The optical path length. Returns 0 if the ray's intensity is zero.
            If the ray has no child rays, returns infinity.
        """

        if self.intensity==0:
            return 0.
        elif len(self.childs)==0:
            return INFINITY
        else:
            return self.childs[0].optical_path_parent()

    def __eq__(self, other):

        cdef Ray other_ray = <Ray> other

        return Ray.almost_equal(self, other_ray, 1e-10)
        # TODO do we have to compare self.pop and other.pop ?
        # self.copy() indicate that we should actually only compare
        # fields pos, dir, intensity, wavelength, n and label

    @staticmethod
    def almost_equal(Ray ray1, Ray ray2, double tol=1e-7):

        cdef bint is_n_eq
        if isnan(ray1.n) and isnan(ray2.n):
            is_n_eq = True
        else:
            is_n_eq = (ray1.n == ray2.n)

        return (is_approx(ray1._origin, ray2._origin, tol) and
                is_approx(ray1._direction, ray2._direction, tol) and
                ray1.intensity == ray2.intensity and
                ray1.wavelength == ray2.wavelength and
                is_n_eq and
                ray1.label == ray2.label and
                ray1.draw_color == ray2.draw_color and
                ray1.order == ray2.order and
                ray1.orig_surf == ray2.orig_surf)
