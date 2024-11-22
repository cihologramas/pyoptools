# ------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amezquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amezquita Orozco
# Description:     Shape definition module
# Symbols Defined: Shape
# ------------------------------------------------------------------------------
#

# from pyoptools.misc.picklable.picklable cimport Picklable
from pyoptools.misc.cmisc.eigen cimport Vector3d, assign_tuple_to_vector3d

cdef class Shape:  # (Picklable):
    """
    Abstract superclass for all optical surface shapes.

    `Shape` is an abstract superclass that defines the interface for different
    surface shapes (e.g., circular, rectangular). This class provides an API that
    all subclasses must implement to define specific shapes and behaviors.

    Subclasses are required to implement methods for checking whether a point is
    inside the shape, defining the shape's limits, and generating lists of points
    that sample the shape adequately.

    Attributes
    ----------
    topo : callable
        A function Z(x, y) that describes the topography of the surface.
        This attribute should be initialized as needed.
    """

    def __init__(self):
        # Picklable.__init__(self)
        # TODO: Check if this needs to be picklable
        pass

    def hit(self, point):
        """
        Determine if a point is within the surface aperture.

        This method checks whether a given point `(x, y, z)` lies inside the
        surface aperture. It returns `True` if the point is within the aperture
        and `False` otherwise.

        The method relies on the `hit_cy` method, which must be implemented
        in all subclasses of `Shape`.

        Parameters
        ----------
        point : tuple of float
            A tuple representing the coordinates `(x, y, z)` of the point to
            check.

        Returns
        -------
        bool
            `True` if the point is within the surface aperture, `False` otherwise.
        """
        cdef Vector3d point_vect
        assign_tuple_to_vector3d(point, point_vect)
        return self.hit_cy(point_vect)

    cdef bint hit_cy(self, Vector3d& point)  noexcept nogil:
        """
        Determine if a point is inside the aperture defined by the shape.

        This method must be implemented in subclasses of `Shape`. It contains
        the algorithm necessary to determine whether a given point `(x, y, z)`
        is inside or outside of the aperture defined by the specific shape.

        Parameters
        ----------
        point : Vector3d&
            A reference to a `Vector3d` object representing the coordinates
            `(x, y, z)` of the point to be checked.

        Returns
        -------
        bool
            `True` if the point is within the aperture, `False` otherwise.

        Notes
        -----
        This method is abstract in the `Shape` class and must be implemented
        in all subclasses. The implementation should include the specific
        algorithm for determining whether the point lies within the shape's
        aperture.
        """
        with gil:
            raise NotImplementedError("The hit_cy method must be "
                                      "overloaded in the subclass " +
                                      self.__class__.__name__)

    cpdef pointlist(self):
        """
        Generate a list of points that adequately sample the shape.

        This method should return a tuple `(X, Y)` where `X` contains the X
        coordinates of the points and `Y` contains the Y coordinates. The
        points must be sampled adequately to be used as vertices for generating
        the triangles required to plot the surface defined by `(X, Y, self.topo(X, Y))`.

        The method must be implemented in each subclass to ensure that the
        sampling is appropriate for the specific shape.

        Returns
        -------
        tuple of lists
            A tuple `(X, Y)` where `X` is a list of X coordinates and `Y` is a
            list of Y coordinates, sampled adequately for surface plotting.

        Raises
        ------
        NotImplementedError
            If this method is not implemented in a subclass, an error is raised
            indicating that the method must be overloaded.
        """
        raise NotImplementedError("The `pointlist` method must be overloaded in " +
                                  "class " + self.__class__.__name__)

    # def mesh(self, size=None, grid_size=(100, 100), topo=None):
    #     """
    #     Generate a mesh grid for the surface.

    #     Creates a mesh grid for the surface, optionally limited by the `size`
    #     parameter. If `size` is not provided, the size will be determined by
    #     the aperture limits.

    #     Parameters
    #     ----------
    #     size : tuple of float, optional
    #         The limits `(xi, xf, yi, yf)` of the mesh grid. If not provided,
    #         the limits are derived from the aperture size.
    #     grid_size : tuple of int, optional
    #         The number of points in the mesh along each axis `(nx, ny)`.
    #         Defaults to `(100, 100)`.
    #     topo : callable, optional
    #         A function `Z(x, y)` that describes the topography of the surface.
    #         If not provided, the function will return 0 for points outside the
    #         aperture and 1 for points inside.

    #     Returns
    #     -------
    #     tuple of arrays
    #         A tuple `(XM, YM, Z)` where:
    #         - `XM` is the X coordinates of the mesh grid.
    #         - `YM` is the Y coordinates of the mesh grid.
    #         - `Z` is the Z values computed using the `topo` function.
    #     """
    #     if size is None:
    #         xi, xf, yi, yf = self.limits()
    #     else:
    #         xi, xf, yi, yf = size

    #     nx, ny = grid_size
    #     X_M = linspace(xi, xf, int(nx))
    #     Y_M = linspace(yi, yf, int(ny))
    #     XM, YM = meshgrid(X_M, Y_M)

    #     if topo is None:
    #         Z = self.hit((XM, YM, 0))
    #     else:
    #         Z = topo(XM, YM)
    #     return XM, YM, Z

    cpdef limits(self):
        """
        Return the minimum and maximum limits for the aperture.

        This method should be overridden in each subclass to return the
        specific limits of the shape.

        Returns
        -------
        tuple of float
            The minimum and maximum limits `(xi, xf, yi, yf)` of the aperture.
        """
        raise NotImplementedError("The `limits` method must be overloaded " +
                                  "in the subclass " + self.__class__.__name__)
