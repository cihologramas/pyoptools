# cython: profile=True

'''Ray class definition module
'''

cdef extern from "math.h":
    double sqrt(double)
cimport cython

from sys import exit
# Se van a usar los traits para el desarrollo del programa
# Se cambio el traits por el HasStrictTraits, para verificar la asignación
# de los atributos en la creacion de los scripts

# from enthought.traits.api import Float, List, Trait, TraitHandler,\
#     TraitList,Array,Property, HasPrivateTraits, String, This, \
#     Any, Instance

#from enthought.tvtk.api import tvtk

import numpy as np
from numpy import dot, array, float_, inf, float64, empty, zeros, sqrt as npsqrt
cimport numpy as np

#from misc import rot_mat, rot_mat_i,mvdot, dot_test
from pyoptools.misc.cmisc.cmisc cimport *  # rot_mat, rot_mat_i

#import tokyo
# cimport tokyo
from cpython.ref cimport PyObject, Py_INCREF, Py_DECREF

cimport numpy as np

np.import_array()

# class Ray(HasPrivateTraits):

# Trick to create instances of ray fast from a cython code. From python
# use the standard
# Note: when changing the code, always check that the __init__, and the
# Rayf initialization do the same.


cdef Ray Rayf(np.ndarray pos, np.ndarray dir, double intensity, double wavelength,
              n, label, draw_color, parent, double pop, orig_surf, int order):
    """Function to create and initialize Ray instances Fast from cython
    when changing the code, always check that the __init__, and the
    Rayf initialization do the same.
    """
    # Code taken from:
    # http://wiki.cython.org/FAQ#CanCythoncreateobjectsorapplyoperatorstolocallycreatedobjectsaspureCcode.3F
    cdef Ray instance = Ray.__new__(Ray)

    instance.cpos=pos

    # Be careful no normalization is made here
    instance._dir=dir

    instance.intensity=intensity
    instance.wavelength=wavelength
    instance.n=n
    instance.label=label
    instance.draw_color=draw_color
    instance.parent=parent
    instance.pop=pop
    instance.orig_surf=orig_surf
    instance.__childs=[]
    instance.order=order
    return instance


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
    # ~ # (x,y,z) tuple containing the ray origin
    # ~ pos=Array('d',shape=(3,),value=(0.,0.,0.))
    # ~
    # ~ # (x,y,z) tuple containing the unit vector representing the direction of
    # ~ # propagation of the ray (this vector gets normalized automatically)
    # ~ dir=UnitVector()
    # ~
    # ~ # Intensity of the ray.
    # ~ intensity=Float(1)
    # ~
    # ~ # Wavelength of the ray in micrometers (.58929 by default)
    # ~ wavelength=Float(.58929)
    # ~
    # ~ # Refraction index of the point originating the ray. If the value is None,
    # ~ # the ray was emitted from the media (not from inside a component)
    # ~ n=Trait(None,None,Float)
    # ~
    # ~ # Label to follow the rays through the system. This attribute propagates with
    # ~ # the ray trace.
    # ~ label=String("")
    # ~
    # ~ # Ray where this ray comes from. It is the ray before the interaction with the
    # ~ # previous optical surface in the propagation
    # ~ parent= This
    # ~
    # ~ # If this ray has no parent, use this value as the parent optical path. If
    # ~ # the ray has a parent, this value must be ignored and instead the real
    # ~ # must be used.
    # ~ pop= Float(0)
    # ~
    # ~ # Surface where the ray originates. All the rays not originated by interaction
    # ~ # with a Surface have this attribute set to None
    # ~ orig_surf=Any()
    # ~

    # in order for the autodocumentation to work, the method definition mut be in a single line

    def __init__(self, pos=(0, 0, 0), dir=(0, 0, 1), double intensity=1., double wavelength=.58929, n=None, label="", draw_color=None, parent=None, double pop=0., orig_surf=None, order=0):

            self.pos=pos
            self.dir=dir
            self.intensity=intensity
            self.wavelength=wavelength
            self.n=n
            self.label=label
            self.draw_color=draw_color
            self.parent=parent
            self.pop=pop
            self.orig_surf=orig_surf
            self.__childs=[]
            self.order=order

        # HasPrivateTraits.__init__(self,**traits)
    def __reduce__(self):
        args=(self.pos, self.dir, self.intensity, self.wavelength, self.n , self.label, self.draw_color, self.parent, self.pop, self.orig_surf, self.order)
        return(type(self), args, self.__getstate__())

    def __getstate__(self):
        return self.__childs

    def __setstate__(self, state):
        self.__childs=state

    property childs:
        def __get__(self):
            # The list is converted to a tuple, so the internal list does not
            # get exposed, and it can not be modified by accident using an
            # append

            return tuple(self.__childs)

            # All internal lists that uses Properties are converted to
            # tuples to protect internal data.

    # ~ __childs= Trait([],TraitList(Trait(This)))
    # ~

    # Correct way to define properties in cython
    property dir:
        def __get__(self):
            return self._dir

        @cython.boundscheck(False)
        @cython.wraparound(False)
        def __set__(self, dir):
            """
            This should create a normalized direction. Some times due to
            rounding errors |dir|>1. Caution must be taken so this do not create
            a problem
            """
            cdef np.ndarray[np.float64_t, ndim=1] adir=array(dir, dtype=float64)
            cdef double n
            n =sqrt(1./(adir[0]*adir[0]+adir[1]*adir[1]+adir[2]*adir[2]))
            adir[0]=n*adir[0]
            adir[1]=n*adir[1]
            adir[2]=n*adir[2]
            self._dir=adir
            # if (adir[0]*adir[0]+adir[1]*adir[1]+adir[2]*adir[2])>1:
            #    print "warning: The direction of a ray was not normalized correctly"

    # Use property pos to access the field cpos from a python program.
    # It can receive lists, tuples or arrays
    # the cython only visible cpos, can only receive arrays
    property pos:
        def __get__(self):
            return self.cpos

        def __set__(self, pos):
            self.cpos=array(pos, dtype=float64)

    def ch_coord_sys_inv(self, no, ae, childs=False):
        '''Transform the coordinate system of the Ray

        Parameters

        *no*
            Tuple (X,Y,Z) containing the coordinates of the origin of the
            old coordinate system in the new coordinate system

        *ae*
            Tuple (RX,RY,RZ) containing the rotation angles to be applied to
            the  old coordinate system. The rotations are applied RZ first,
            then RY and last RX

        *childs*
            Transform also the coordinate system of the childs. By default (False)
            don't do it.

        The rotation is made first, and then the translation is made.
        Note, this has to be checked
        '''

        cdef np.ndarray tm=rot_mat(array(ae, dtype=float64))
        npos=dot(tm, self.pos)+no
        ndir=dot(tm, self.dir)
        # parent=Ray(pos=npos,dir=ndir,intensity=self.intensity,
        #        wavelength=self.wavelength,n=self.n, label=self.label,
        #        orig_surf=self.orig_surf)

        cdef Ray parent=Rayf(npos, ndir, self.intensity, self.wavelength, self.n, self.label, self.draw_color, None , 0, self.orig_surf, self.order)

        # Calculate the transform of the childs and link them

        if childs==True:
            for i in self.childs:
                it=i.ch_coord_sys_inv(no, ae, childs)
                parent.add_child(it)
        return parent

    cpdef Ray ch_coord_sys_inv_f(self, np.ndarray no , np.ndarray ae, bool childs):
        '''Transform the coordinate system of the Ray

        Fast version

        Parameters

        *no*
            Tuple (X,Y,Z) containing the coordinates of the origin of the
            old coordinate system in the new coordinate system

        *ae*
            Tuple (RX,RY,RZ) containing the rotation angles to be applied to
            the  old coordinate system. The rotations are applied RZ first,
            then RY and last RX

        *childs*
            Transform also the coordinate system of the childs. By default (False)
            don't do it.

        The rotation is made first, and then the translation is made.
        Note, this has to be checked
        '''

        cdef np.ndarray tm=rot_mat(ae.astype(float64))
        cdef np.ndarray npos  # =dot(tm,self.pos)+no
        cdef np.ndarray t

        t=dot(tm, self.pos)
        npos = empty_vec(3)

        cdef np.float64_t* nposd= <np.float64_t*>(np.PyArray_DATA(npos))
        cdef np.float64_t* nod= <np.float64_t*>(np.PyArray_DATA(no))
        cdef np.float64_t* td= <np.float64_t*>(np.PyArray_DATA(t))

        nposd[0]=td[0]+nod[0]
        nposd[1]=td[1]+nod[1]
        nposd[2]=td[2]+nod[2]

        ndir=dot(tm, self.dir)

        cdef Ray parent=Rayf(npos, ndir, self.intensity, self.wavelength, self.n, self.label, self.draw_color, None , 0, self.orig_surf, self.order)

        # Calculate the transform of the childs and link them

        if childs:
            for i in self.childs:
                it=i.ch_coord_sys_inv_f(no, ae, childs, True)
                parent.add_child(it)
        return parent

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef Ray ch_coord_sys(self, np.ndarray no, np.ndarray ae):
    # cpdef Ray ch_coord_sys(self, no, ae):
        '''Transform the coordinate system of the Ray

        Parameters
        no
        Tuple (X,Y,Z) containing the coordinates of the origin of the  new coordinate
        system in the old coordinate system
        ae
        Tuple (RX,RY,RZ) containing the rotation angles to be applied to the  old
        coordinate system. The rotations are applied RZ first, then RY and last RX

        Note this has to be checked
        '''
        # cdef np.ndarray tm=rot_mat_i(array(ae,dtype=float64))
        cdef np.ndarray tm=rot_mat_i(ae)
        cdef np.ndarray npos  # =empty_vec( 3 )
        cdef np.ndarray ndir  # =empty_vec( 3 )

        # mvdotf(<np.float64_t*>np.PyArray_DATA(npos),
        #        <np.float64_t*>np.PyArray_DATA(tm),<np.float64_t*>np.PyArray_DATA(array(self.pos-no,dtype=float64)))

        # mvdotf(<np.float64_t*>np.PyArray_DATA(ndir),
        #        <np.float64_t*>np.PyArray_DATA(tm),<np.float64_t*>np.PyArray_DATA(array(self.dir,dtype=float64)))

        cdef np.ndarray t = empty_vec(3)

        # t=self.cpos-no

        # (<np.float64_t*>np.PyArray_DATA(t))[0]=(<np.float64_t*>np.PyArray_DATA(self.cpos))[0]-(<np.float64_t*>np.PyArray_DATA(no))[0]
        # (<np.float64_t*>np.PyArray_DATA(t))[1]=(<np.float64_t*>np.PyArray_DATA(self.cpos))[1]-(<np.float64_t*>np.PyArray_DATA(no))[1]
        # (<np.float64_t*>np.PyArray_DATA(t))[2]=(<np.float64_t*>np.PyArray_DATA(self.cpos))[2]-(<np.float64_t*>np.PyArray_DATA(no))[2]

        cdef np.float64_t* td= <np.float64_t*>(np.PyArray_DATA(t))
        cdef np.float64_t* cd= <np.float64_t*>(np.PyArray_DATA(self.cpos))
        cdef np.float64_t* nd= <np.float64_t*>(np.PyArray_DATA(no))

        td[0]=cd[0]-nd[0]
        td[1]=cd[1]-nd[1]
        td[2]=cd[2]-nd[2]

        npos=dot(tm, t)        # if blas can be called directly, this can be improved
        ndir=dot(tm, self._dir)  # Using tokio this tokio works a little faster, but have to see
                              # How to install it right, or better to do something similar using inline snd blas

        #tokyo.dgemv3( tm, t, npos )
        #tokyo.dgemv3( tm, self.dir, ndir )
        # This resulted to be slower
        # npos=mvdot1(tm,t)
        # ndir=mvdot1(tm,self.dir)

        return Rayf(npos, ndir, self.intensity, self.wavelength, self.n, self.label, self.draw_color, None , 0, self.orig_surf, self.order)

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
        # TODO: why 'direc' and not 'dir' ? This is confusing
        return "Ray(pos="+repr(self.pos)+",direc="+repr(self.dir)+\
            ",intensity="+repr(self.intensity)+",wavelength="+\
            repr(self.wavelength)+",n="+repr(self.n)+",label="+\
            repr(self.label)+",orig_surf="+repr(self.orig_surf)+", order="+repr(self.order)+")"

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
        self.__childs.append(cr)

    def optical_path_parent(self):
        ''' Return the optical path from the origin of the origin ray to the
       end of this ray parent (this ray origin)
        '''

        if self.parent is not None:
            if self.pop!=0:
                print "The pop attribute of the ray has a value of ", self.pop, \
                " instead the real parent optical path is being used"
            path= sqrt(dot(self.pos-self.parent.pos, self.pos-self.parent.pos))*\
            self.parent.n
            return path+self.parent.optical_path_parent()

        return self.pop

    def optical_path(self):
        ''' Return the optical path of the beam propagation from the origin of
        the origin ray, to the end of this ray (intersection with a surface)
        '''

        if self.intensity==0:
            return 0.
        elif len(self.childs)==0:
            return inf
        else:
            return self.childs[0].optical_path_parent()

    def __eq__(self, other):
        return np.array_equal(self.pos, other.pos) and\
        np.array_equal(self.dir, other.dir) and\
        self.intensity == other.intensity and\
        self.wavelength == other.wavelength and\
        self.n == other.n and\
        self.label == other.label and\
        self.draw_color == other.draw_color and\
        self.order == other.order and\
        self.orig_surf == other.orig_surf
        # TODO do we have to compare self.pop and other.pop ?
        # self.copy() indicate that we should actually only compare
        # fields pos, dir, intensity, wavelength, n and label

    @staticmethod
    def almost_equal(ray1, ray2, decimal=7):
        """
        Test if two rays are equal up to desired precision.

        **ARGUMENTS**

            ======= ===============================================
            other   Ray other ray to be compared against.
            decimal int, optional, Desired precision, default is 7.
            ======= ===============================================


        **RETURN VALUE**
            bool   True if two rays are equal to the desired precision. i.e.
                   if abs(self - other) < 1.5 * 10**(-decimal) for attributes
                   pos, dir, wavelength and n.
        """
        atol = 1.5 * 10 ** (-decimal)
        rtol = 0

        return np.allclose(ray1.pos, ray2.pos, rtol=rtol, atol=atol) and \
        np.allclose(ray1.dir, ray2.dir, rtol=rtol, atol=atol) and \
        ray1.intensity == ray2.intensity and \
        np.allclose(ray1.wavelength, ray2.wavelength, rtol=rtol, atol=atol) and \
        ray1.n == ray2.n and \
        ray1.label == ray2.label and \
        ray1.draw_color == ray2.draw_color and \
        ray1.order == ray2.order and \
        ray1.orig_surf == ray2.orig_surf
        # TODO do we have to compare self.pop and other.pop ?
        # self.copy() indicate that we should actually only compare
        # fields pos, dir, intensity, wavelength, n and label
