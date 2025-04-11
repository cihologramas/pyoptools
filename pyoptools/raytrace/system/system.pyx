# ------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author: Ricardo Amézquita Orozco
# Description: Definition of the System class.
# ------------------------------------------------------------------------------


"""Module that defines the optical system class System()
"""

from numpy import asarray, array, all, isinf as npisinf

from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.plist.plist cimport plist

from pyoptools.misc.picklable.picklable cimport Picklable
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.component.component cimport Component

from libc.math cimport isinf, INFINITY, isnan

cdef class System(Picklable):
    """
    Class to define an optical system.

    The System class defines an optical system as a list of optical
    components, the coordinates of the component origin, and the rotation
    angles of such components. To define a system, the refraction index, or
    the material surrounding the components must also be given in the `n`
    attribute.

    To avoid recursion errors, two other parameters are given:
    `max_ray_parent` and `intensity_threshold`. See the parameters
    description to understand their use.

    Parameters
    ----------
    complist : list of tuples
        Contains a tuple list that defines the optical system. The first
        component of the tuple is an instance of the component to include in
        the system. The second is a tuple with the component location (X, Y,
        Z). The third is a tuple with the rotation around each axis (rx, ry,
        rz). The rotations are applied in the order: rz, then ry, and then rx.
    n : float
        Contains the refraction index of the media where the system is
        immersed.
    max_ray_parent : int, optional
        Limits the number of parents and grandparents a ray can have,
        effectively limiting the number of times a ray can be propagated.
        This is used to avoid recursion errors in case a resonant cavity is
        simulated. Default is None, meaning there is no limit for the
        propagation.
    intensity_threshold : float, optional
        Another way to limit the number of times a ray is propagated. This is
        used to avoid recursion errors in case a resonant cavity is simulated.
        If a ray has an intensity less than the threshold, the ray is not
        propagated. The intensity value of a ray is reduced each time it is
        reflected in a beam splitter. Default is 0, meaning there is no limit
        for the propagation.

    Examples
    --------
    Example of a system containing a doublet and a CCD detector::

        # Definition of a doublet type lens
        DB1 = Doublet(radius=12.5,
                    curvature_as=1./61.47,
                    curvature_ms=-1./44.64,
                    curvature_ps=-1./129.94,
                    thickness_al=6.,
                    thickness_pl=2.5,
                    material_al=N_BK7,
                    material_pl=SF5)

        # Definition of a CCD type detector
        ccd = CCD()

        # Definition of a system
        os = System(complist=[(DB1, (20, 0, 200), (0, 0, 0)),
                            (ccd, (20, 0, 400), (0, 0, 0))],
                    n=1)
    """

    @property
    def max_ray_parent(self):
        if self._max_ray_parent > 0:
            return self._max_ray_parent
        else:
            return None

    @max_ray_parent.setter
    def max_ray_parent(self, val):

        if val is None:
            self._max_ray_parent = 0
        else:
            self._max_ray_parent = int(val)

    @property
    def prop_ray(self):
        return tuple(self._p_rays)

    ###############################################################

    @property
    def complist(self):
        return self._complist

    @complist.setter
    def complist(self, list):
        self._complist=plist(list)

    def __init__(self, complist=None, n=1., max_ray_parent=None,
                 intensity_threshold=0):

        # Look in the init os component to see why this in done this way
        if complist is None:
            self.complist=[]
        else:
            self.complist=complist
        self.n=n
        self._np_rays=[]  # not propagated rays
        self._p_rays=[]  # propagated rays

        # self.propagation_limit is integer, so if it is 0, then there is no limit

        self.max_ray_parent = max_ray_parent

        self.intensity_threshold = intensity_threshold

        # Flag that indicates if a ray propagation was truncated or not by the
        # intensity_threshold or the max_ray_parent condition
        # if 0 no truncation was done

        self._exit_status_flag = 0

        for i in self.complist:
            comp=i[0]
            # If the component is a subsystem, the refraction index must be the
            # same of the system, also the propagation limits must be the same

            if isinstance(comp, System):
                comp.n=self.n
                comp.max_ray_parent = self.max_ray_parent
                comp.intensity_threshold = self.intensity_threshold

        # Add to the keys to the state key list
        Picklable.__init__(self, "complist", "n", "_np_rays", "_p_rays",
                           "_max_ray_parent", "intensity_threshold")

    # Dict type and list type interface to expose complist

    def __len__(self):
        return self._complist.__len__()

    def __getitem__(self , x):
        return self._complist[x]

    def __setitem__(self, key, val):
        self._complist[key]=val

    def __delitem__(self, key):
        self._complist.__delitem__(key)

    def __contains__(self, key):
        return self._complist.__contains__(key)

    # Return an iterator so this can be used similar to a list
    def __iter__(self):
        return self._complist.itervalues()

    def iteritems(self):
        return self._complist.iteritems()

    def iter(self):
        return self._complist.iter()

    def clear(self):
        return self._complist.clear()

    def items(self):
        return self._complist.items()

    def iterkeys(self):
        return self._complist.iterkeys()

    def itervalues(self):
        return self._complist.itervalues()

    def keys(self):
        return self._complist.keys()

    def pop(self, *argv, **argkw):
        return self._complist.pop(*argv, **argkw)

    def popitem(self):
        return self._complist.popitem()

    def setdefault(self, *argv, **argkw):
        return self._complist.setdefault(*argv, **argkw)

    def update(self, *argv, **argkw):
        return self._complist.update(*argv, **argkw)

    def values(self):
        return self._complist.values()

    def viewitems(self):
        return self._complist.viewitems()

    def viewkeys(self):
        return self._complist.viewkeys()

    def viewvalues(self):
        return self._complist.viewvalues()

    def clear_ray_list(self):
        """ Clear the ray lists of the system
        """
        self._np_rays=[]
        self._p_rays =[]

    def ray_add(self, ray):
        """
        Rutina que adiciona un rayo a la lista de rayos primarios del sistema
        optico.
        Recibe como parametro una rayo o una lista de rayos. Genera un error
        si se le pasa algo diferente a un rayo (instancia de la clase Ray,
        genera una excepcion        '''
        """

        if isinstance(ray, (list, tuple)):
            for i in ray:
                if isinstance(i, Ray):
                    self._np_rays.append(i)
                else:
                    raise Exception, "Not a valid Ray"
        elif isinstance(ray, Ray):
            self._np_rays.append(ray)
        else:
            raise Exception, "Not a valid Ray"

    def propagate(self, update_ids=True):
        """ Propagates all the rays in the non propagated list.
        """

        # This is not necessary for all propagations, but is safer to do it
        # When propagating a sub system this must not be done
        if update_ids:
            self.update_ids()

        while len(self._np_rays)>0:
            ri=self._np_rays.pop(0)
            self.propagate_ray(ri)
            self._p_rays.append(ri)

    def get_surf_paths(self):
        '''Method that returns a list that contains the path for each surface.

        A path here is a list containing the keys needed to read each surface.

        This method is an auxiliary method so this works when called from a System

        '''
        l=[]
        keys=self.keys()
        for k in keys:
            a=self[k][0].get_surf_paths()
            for k1 in a:
                l.append([k]+k1)
            if len(a)==0:
                l.append([k])
        return l

    def get_surface(self, path):
        '''
        Return a surface, given a path.

        A path is given as a list of keys
        '''

        O=self
        for k in path:
            try:
                O=O[k][0]
            except KeyError:
                raise KeyError, "Invalid path.  Key %s does not exist" %k
            except TypeError:
                raise TypeError, "Invalid path. Path too long, key %s does not exist" %k
        assert isinstance(O, Surface), "Error in path: Path too short"

        return O

    def get_component(self, path):
        '''
        Return the component thatis defined using the surface described by path
        '''
        O=self
        for k in path:
            try:
                C=O
                O=O[k][0]
            except KeyError:
                raise KeyError, "Invalid path.  Key %s does not exist" %k
            except TypeError:
                raise TypeError, "Invalid path. Path too long, key %s does not exist" %k
        assert isinstance(C, Component), "Error in path: Path too short"
        return C

    def update_ids(self):
        """
        Update the ids for all the surfaces in the system.

        This should be done before running a propagation.
        """
        paths=self.get_surf_paths()
        for path in paths:
            S=self.get_surface(path)
            S.id=path

    def __repr__(self):
        '''Return an string with the representation of the optical system

        It must be overloaded in all subclasses
        '''

        ret="OpSys(\n"
        for i in self.complist:
            ret=ret+repr(i)+",\n"

        ret=ret+")"

        return ret

    def reset(self):
        """
        Run the reset method on all the components used to create the system

        """
        self._np_rays=[]
        self._p_rays=[]
        for comp in self.complist:
            S, _P, _D=comp
            S.reset()

    cpdef propagate_ray(self, Ray ri):
        """
        Method to propagate the ray in the system. It creates the nexts rays
        and links them using the Ray.parent, and Ray.childs attributes. It calls
        itself recurrently.

        Arguments:


        *ri*
            Ray to propagate

        Return Value

        It returns *ri*
        """
        # TODO: All surfaces, elements and subsystems should know their
        # coordinates, so the transformations can be made a lot faster
        # the only problem of this approach is that the surfaces can not be
        # reused in a design. They must be copied to be reused. The same will
        # Happen to the components and subsystems.

        # Check if the ray comes from the media

        # These are defined as tuples, because they need to be python objects
        # to be included in lists

        cdef tuple[double, double, double] P, D, PSR, DSR, PSR0, DSR0, PSR1, \
            DSR1

        cdef int j, j1
        cdef double d0, d1

        if isnan(ri.n):
            ri.n=self.n

        cdef list dist_list=[]
        cdef list surf_list=[]
        cdef list comp_list=[]
        cdef list pi_list=[]
        # Calculate the path length followed by the ray until it intersects all
        # the components and subsystems

        # Note: C can be component or subsystem, so for the moment we will
        # leave it as a python object

        cdef object C
        for i in self.complist:
            C, P, D = i
            comp_list.append((C, P, D))
            # Reorientar el rayo, al sistema de coordenadas del elemento
            # y calcular el recorrido del rayo hasta chocar con la
            # el elemento
            R=ri.ch_coord_sys(P, D)

            Dist=C.distance(R)

            dist_list.append(Dist[0])

            pi_list.append(Dist[1])

            surf_list.append(Dist[2])

        # Check if there are more components in front of the ray
        # if not, return the original ray

        if all(npisinf(array(dist_list))):
            return ri

        # Sort the components by distance
        sort_list=asarray(dist_list).argsort()

        # Take the 2 nearest components. If there is only one component assume the 2nd
        # component at infinitum
        j=sort_list[0]
        d0=dist_list[j]

        # Add ray to the hit list
        surf_list[j]._hit_list.append((pi_list[j], ri))

        # TODO: The hitlists of the surfaces inside a subsystem are not accurate
        # because the rays are in the subsystem coordinate system, and not in
        # world coordinate system.
        if len(sort_list)>1:
            j1=sort_list[1]
            d1=dist_list[j1]
        else:
            d1 = INFINITY
        # Si las compomentes mas cercanas no estan en contacto, calcular la
        # propagacion a travez de la componente mas cercana
        # Nota_: La comparacion de punto flotante no esta funcionando. Para
        # que funcione toca definir un epsilon, en el que se consideran
        # nulas las diferencias

        # If the closest components are not in contact calculate the propagation
        # using the closest surface.

        cdef double N_EPS=1.e-12  # Used to check the zero

        # Check if you are propagating in a subsystem

        # Note: SR can be a System or a Component, so for the moment we will
        # leave it as a standard python object
        cdef object SR

        if isinstance(comp_list[j][0], System):
            # Leer el elemento que primero intersecta el rayo, asi como
            # su posicion y orientacion
            SR, PSR, DSR=comp_list[j]
            # SR.reset()
            SR.clear_ray_list()

            R=ri.ch_coord_sys(PSR, DSR)
            SR.ray_add(R)
            # Ids must not be updated when propagating in subsystems
            SR.propagate(update_ids=False)

            # Change the coordinate system of the propagated ray and its childs
            RT=R.ch_coord_sys_inv(PSR, DSR, childs=True)

            # Link the subsystem rays to the original ray
            for i in RT.childs:
                ri.add_child(i)

        # Verificar si no hay componentes en contacto
        elif abs(d0-d1)>N_EPS:

            # Get the nearest element to the ray origin, as well as its
            # position and orientation
            SR, PSR, DSR=comp_list[j]

            # Change the ray to the coordinate system of the element

            R=ri.ch_coord_sys(PSR, DSR)

            # as there are no components in contact the refraction index outside the
            # is the media's

            ri_n=SR.propagate(R, self.n)

            # Change the coordinate system of the propagated rays to the
            # system coordinate system
            for i in ri_n:
                ri_=i.ch_coord_sys_inv(PSR, DSR)
                # put the rays in the childs list
                ri.add_child(ri_)
        else:

            # There are 2 objects in contactt
            # Object 1

            SR0, PSR0, DSR0=comp_list[j]
            # Object 2
            SR1, PSR1, DSR1=comp_list[j1]
            # Add ray to the hit list
            surf_list[j1]._hit_list.append((pi_list[j1], ri))
            n0=SR0.n(ri.wavelength)
            n1=SR1.n(ri.wavelength)
            # print 1
            # Calculate the refraction for both components

            R0=ri.ch_coord_sys(PSR0, DSR0)

            ri_n0=SR0.propagate(R0, n1)

            R1=ri.ch_coord_sys(PSR1, DSR1)
            ri_n1=SR1.propagate(R1, n0)

            # TODO: Need to find a solution when the two surfaces return more
            # than one ray.
            if (len(ri_n0)>1)and(len(ri_n1)>1):
                raise Exception, "The two surfaces in contact, can not produce "\
                                 "both more than one propagated ray"
            elif len(ri_n0)>1:
                for i in ri_n0:
                    ri_=i.ch_coord_sys_inv(PSR0, DSR0)
                    # put the rays in the childs list
                    ri.add_child(ri_)
            elif len(ri_n1)>1:
                for i in ri_n1:
                    ri_=i.ch_coord_sys_inv(PSR1, DSR1)
                    # put the rays in the childs list
                    ri.add_child(ri_)
            else:
                ri_0=ri_n0[0].ch_coord_sys_inv(PSR0, DSR0)
                # ri_1=ri_n1[0].ch_coord_sys_inv(PSR1, DSR1)
                # TODO: ri_0 and ri_1 must be equal. Needs to be checked
                ri.add_child(ri_0)

        # Propagate childs

        for i in ri.get_final_rays():
            if (i!=ri):

                # stop propagating if the ray has an intensity below the
                # intensity propagation threshold or if it has been propagated
                # more than propagation_limit times

                if i.intensity>self.intensity_threshold:
                    if (self._max_ray_parent == 0 or
                       i._parent_cnt<self.max_ray_parent):
                        self.propagate_ray(i)
                    else:
                        self._exit_status_flag = 1
                else:
                    self._exit_status_flag = 1

            else:
                raise Exception, \
                    "Error, a a ray can not be parent and child at the same time"

        return ri

    cpdef propagate_ray_ns(self, Ray gr, dpath):
        '''Method to propagate the ray in the system.

        **Arguments**

            ===== ========================================================
            gr    Guide ray previously propagated in the system using the
                  non sequential algorithm. This ray contains the surface
                  sequence that the rays must follow.
            dpath Path (key) of the destination surface.
            ===== ========================================================

        This method uses the same n as the calculated in the non sequential
        propagation. If the wavelength change, assertion error is raised

        '''
        from pyoptools.raytrace.calc import ray_paths
        cdef list spath, Olist, paths
        cdef int i
        cdef Ray ri
        # Get all the paths traveled by the ray
        paths=ray_paths(gr)

        # Find the path of interest
        rp=[]

        for p in paths:
            if p[-1].orig_surf==dpath:
                rp=p
                break
        assert rp!=[], "Guide ray does not intersect the given surface"

        while len(self._np_rays)>0:
            ri=self._np_rays.pop(0)
            assert ri.wavelength==gr.wavelength, \
                "Propagated rays, and guide ray wavelength must match"
            # self.propagate_ray(ri)
            # ~ # Check if the ray comes from the media
            if isnan(ri.n):
                ri.n=self.n

            self._p_rays.append(ri)

            for i in range(len(rp)-1):

                # This uses the same n as the calculated in the non sequential
                # propagation. If the wavelength change, there is problem
                ni=rp[i].n
                nr=rp[i+1].n

                # Transform the ray to the surfaces coordinate system
                R=ri
                spath=rp[i+1].orig_surf
                # Nota, esto se puede hacer una sola vez antes, y meterlo
                # En una lista
                S=self.get_surface(spath)

                # Transform the ray to the surfaces coordinate system
                O=self
                Olist=[]
                for si in spath:
                    O=O[si]
                    Olist.append(O)
                    C, P, D=O
                    R=R.ch_coord_sys(P, D)
                    O=C

                # Get the distance to the next surface
                Dist, p_i=S.distance(R)

                if isinf(Dist):
                    break  # No intersection continue with next ray

                # Add ray to the hit list
                S._hit_list.append((p_i, ri))

                # TODO: The hitlists of the surfaces inside a subsystem are not accurate
                # because the rays are in the subsystem coordinate system, and not in
                # world coordinate system.

                # Need to check which is the real one to take
                ri_n=S.propagate(R, ni, nr)[rp[i+1].order]

                # Change the coordinate system of the propagated rays to the
                # system coordinate system

                for O in reversed(Olist):
                    C, P, D=O
                    ri_n=ri_n.ch_coord_sys_inv_f(P, D, False)

                ri.add_child(ri_n)
                ri=ri_n
                if ri.intensity==0:
                    break
        return rp

    cpdef distance(self, Ray ri):
        """Distance length from a ray origin to a subsystem, following the ray path.

        Method that calculates the distance traveled by a ray from its origin to
        the next surface of the component. It returns the physical distance, not
        the optical distance


        *Return value*
            A tuple with the distance, the point of intersection using the
            coordinate system of the surface, and a pointer to the surface
            that is closest to the ray (distance,point of intersection, surface)
        """
        # cdef np.ndarray P,D
        cdef list dist_list=[]
        cdef list pi_list=[]
        cdef list surf_list=[]
        # print self.complist
        for comp in self.complist:
            C, P, D =comp
            # C,P,D = i
            # Reorientar el rayo, al sistema de coordenadas del elemento
            # y calcular el recorrido del rayo hasta chocar con la
            # el elemento
            R=ri.ch_coord_sys(P, D)

            Dist=C.distance(R)

            dist_list.append(Dist[0])

            pi_list.append(Dist[1])
            surf_list.append(Dist[2])

        mini=asarray(dist_list).argmin()

        return dist_list[mini], pi_list[mini], surf_list[mini]

    def merge(self, os):
        """
        Method to merge simulation systems. Useful for joining raytraces that
        where split for parallel processing.
        """

        # Check that the 2 optical systems are equal. Right now the check is that
        # the surfaces paths of both systems are the same. There is no other check
        # TODO: Fix this

        cdef list path1=self.get_surf_paths()
        cdef list path2=os.get_surf_paths()

        assert path1==path2, "Different optical systems can not be merged"

        # Transfer the non propagated rays
        self._np_rays=self._np_rays+os._np_rays
        os._np_rays=[]

        # Transfer the propagated rays
        self._p_rays= self._p_rays+os._p_rays
        os._p_rays=[]

        # Transfer the hitlist information
        for p in path1:
            S1=self.get_surface(p)
            S2=os.get_surface(p)
            S1._hit_list=S1._hit_list+S2._hit_list
            S2._hit_list=[]
        # del(os)
