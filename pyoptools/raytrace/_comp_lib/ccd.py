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
# Description:     CCD definition module
# Symbols Defined: CCD
# ------------------------------------------------------------------------------
#
"""
Definition of a CCD like object and helper functions
"""
from PIL.Image import fromarray
from scipy.interpolate import interp2d, bisplrep, bisplev
from numpy import arange, ma, meshgrid, linspace

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import ArrayDetector, Plane
from pyoptools.misc.pmisc import wavelength2RGB
from pyoptools.misc.lsq import polyfit2d
from pyoptools.raytrace.shape import Shape
from pyoptools.raytrace.shape import Rectangular

# from gui.plotutils import plot, figure, cm,  legend


class CCD(Component):
    """Class to define a CCD like detector

    :param size: Tuple with the physical size (sx,sy)  of the CCD chip
    :type size: tuple(float, float)
    :param transparent: Boolean to set the detector transparent characteristic.
        Not implemented
    :type transparent: bool

    Using the same CCD, images of different resolutions can be simulated. See
    the im_show and spot_diagram methods
    """

    # Geometrical size of the CCD chip
    # size = Tuple(Float(5),Float(5))

    # Setting this attribute to *False*, make the CCD detector opaque
    # transparent=Bool(True)

    # Private attributes

    # detector surface

    # __d_surf = Instance(ArrayDetector)
    # __d_surf = Instance(Plane)

    def _get_hitlist(self):
        return tuple(self.__d_surf.hit_list)

    hit_list = property(_get_hitlist)
    """List containing a tuple for each ray hitting the CCD. The first component
    of the tuple is the coordinates of intersection of the ray with the CCD
    (in its coordinate system). The second component of each tuple points to the
    :class:`~pyoptools.raytrace.ray.Ray` that intersected the CCD.
    """

    def __init__(self, size=(10, 10), transparent=True, *args, **kwargs):
        Component.__init__(self, *args, **kwargs)
        self.__d_surf = Plane(
            shape=Rectangular(size=size)
        )  # ArrayDetector (size=self.size, transparent=self.transparent)
        self.size = size
        self.surflist["S1"] = (self.__d_surf, (0, 0, 0), (0, 0, 0))
        self.material = 1.0

    # ~ def __reduce__(self):
    # ~ args=() #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
    # ~ return(type(self),args,self.__getstate__())
    # ~
    # ~
    # ~ #TODO: Check if there is a better way to do this, because we are
    # ~ #rewriting the constructor values here
    # ~
    # ~ def __getstate__(self):
    # ~ return self.__d_surf,self.size,self.surflist,self.material
    # ~
    # ~ def __setstate__(self,state):
    # ~ self.__d_surf,self.size,self.surflist,self.material=state

    def get_image(self, size=(256, 256)):
        """
        Returns the ccd hit_list as a grayscale PIL image

        :param size: Tuple (dx,dy) containing the image size in pixels. Use this
            a ttribute to set the simulated resolution.
        """
        data = self.__d_surf.get_histogram(size)
        return fromarray(data)

    def get_color_image(self, size=(256, 256)):
        """
        Returns the CCD hit_list as a color image, using the rays wavelength.

        :param size: Tuple (dx,dy) containing the image size in pixels. Use this
            attribute to set the simulated resolution.
        """

        data = self.__d_surf.get_color_histogram(size)
        return fromarray(data, high=255, low=0)

    # ~ def im_show(self,fig=None, size=(256,256),cmap=cm.gray,title='Image',color=False):
    # ~ """Shows a simulated image
    # ~
    # ~ *Attributes:*
    # ~
    # ~ *size*
    # ~ Tuple (dx,dy) containing the image size in pixels. Use this
    # ~ attribute to set the simulated resolution.
    # ~ *cmap*
    # ~ Color map to use in the image simulation. See the matplotlib.cm
    # ~ module for information about colormaps.
    # ~ *fig*
    # ~ Pylab figure where the plot will be made. If set to None
    # ~ a new figure will be created.
    # ~ """
    # ~ if fig == None:
    # ~ fig=figure()
    # ~
    # ~ self.__d_surf.im_show(size,cmap,title,color)
    # ~
    # ~
    # ~ def spot_diagram(self,fig=None, style="o",  label=None):
    # ~ '''Plot a spot diagram in a pylab figure
    # ~
    # ~ Method that plots a spot diagram of the rays hitting the CCD.
    # ~
    # ~ *Attributes:*
    # ~
    # ~ *fig*
    # ~ Pylab figure where the plot will be made. If set to None
    # ~ a new figure will be created.
    # ~
    # ~ *style*
    # ~ Symbol to be used to represent the spot. See the pylab plot
    # ~ documentation for more information.
    # ~
    # ~ *label*
    # ~ String containing the label to show in the figure for this spot diagram.
    # ~ Can be used to identify different spot diagrams on the same figure.
    # ~ '''
    # ~
    # ~ if fig == None:
    # ~ fig=figure()
    # ~ X=[]
    # ~ Y=[]
    # ~ COL=[]
    # ~ if len(self.__d_surf._hit_list) >0:
    # ~ for i in self.__d_surf._hit_list:
    # ~ p=i[0]
    # ~ # Hitlist[1] points to the incident ray
    # ~ col=wavelength2RGB(i[1].wavelength)
    # ~ X.append(p[0])
    # ~ Y.append(p[1])
    # ~ COL.append(col)
    # ~ if label== None:
    # ~ plot(X, Y, style,  figure=fig)
    # ~ else:
    # ~ plot(X, Y, style,label=label,figure=fig)
    # ~ legend()
    # ~ return fig

    def get_optical_path_map(self, size=(20, 20), mask=None):
        """Return the optical path of the rays hitting the detector.

        This method uses the optical path of the rays hitting the surface to
        create a optical path map. The returned value is an interpolation of the
        values obtained by the rays.

        .. warning::

            If the rays hitting the surface are produced by more than one
            optical source, the returned map might not be valid.

        :param size: Tuple (nx,ny) containing the number of samples of the
            returned map. The map size will be the same as the CCD
        :param mask: :class:`~pyoptools.raytrace.shape.Shape`
            instance containing the mask of the aperture. If not given,
            the mask will be automatically calculated.

        :return: A masked array as defined in the numpy.ma module, containing
            the optical paths
        """

        X, Y, Z = self.get_optical_path_data()

        rv = bisplrep(X, Y, Z)
        nx, ny = size
        xs, ys = self.size
        xi = -xs / 2.0
        xf = -xi
        yi = -ys / 2.0
        yf = -yi

        xd = linspace(xi, xf, nx)
        yd = linspace(yi, yf, ny)
        data = bisplev(xd, yd, rv)

        if mask is not None:
            assert isinstance(mask, Shape)
            X, Y = meshgrid(xd, yd)
            m = ~mask.hit((X, Y, 0))
            retval = ma.array(data, mask=m)
        else:
            retval = data
        return retval

    def get_optical_path_map_lsq(self, order=10):
        """Return a 2D polinomial describing the the optical path of the
           rays hitting the detector.


        :param order: Order of the polynomial used to fit the data

        :return: tuple (e, p) where e es the rms error of the data when compared
           with the returned polynomial, and p is a
           :class:`~pyoptools.misc.Poly2D.poly2d` instance.
        """

        X, Y, Z = self.get_optical_path_data()
        e, p = polyfit2d(X, Y, Z, order=order)
        return e, p

    def get_optical_path_data(self):
        """Return the optical path of the rays hitting the detector.

        This method returns a tuple X,Y,D, containing the X,Y hit points, and
        D containing the optical path data

        .. warning::

            If the rays hitting the surface are produced by more than one
            optical source, the information may not be valid.

        """

        X = []
        Y = []
        Z = []
        for ip, r in self.hit_list:
            x, y, z = ip
            d = r.optical_path()
            X.append(x)
            Y.append(y)
            Z.append(d)

        return X, Y, Z
