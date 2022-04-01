#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Method collection to obtain optical system information

This module contains a method collection to obtain information, and analyze 
optical systems
"""

__all__ = [
    "intersection",
    "nearest_points",
    "chief_ray_search",
    "pupil_location",
    "paraxial_location",
    "find_aperture",
    "find_ppp",
    "get_optical_path_ep",
    "find_reference_sphere_radius",
    "parallel_propagate",
    "parallel_propagate_ns",
    "ray_paths",
]

import numpy as np

from pyoptools.raytrace.ray import Ray
from pyoptools.misc.pmisc import cross
from pyoptools.raytrace.system import System
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.comp_lib import CCD
from pyoptools.raytrace.surface import Spherical

# from gui.plot_frame import PlotFrame
from pyoptools.raytrace.shape import Circular

from numpy import (
    inf,
    sqrt,
    square,
    pi,
    dot,
    array,
    arctan2,
    alltrue,
    isnan,
    nan,
    mgrid,
    where,
)
from scipy.optimize.minpack import fsolve
from numpy.random import normal
import multiprocessing as mp

# ******Logger definition *******#
# import logging

# log= logging.getLogger("ray_trace.calc")


def intersection(ray1, ray2, atol=1e-8):
    """
    Return the point of intersection between the rays ray1 and ray2.


    Parameters
    ----------

    ray1, ray2: :class:`~pyoptools.raytrace.ray.Ray`
        Rays to test for intersection.
    atol: absolute tolerance (maximum distance between ray lines)

    Returns
    -------

    ip : tuple(float, float, float)
        Intersection point coordinates. If the rays do not intersect
        ip=(nan,nan,nan)
    rv : bool
        Indicates if the intersection point represent a real image (rv=true),
        or a virtual image (rv=false). In this case virtual has the same meaning
        as in virtual image i.e. intersection point is not in the actual path,
        or is behind the ray's origin.
    """
    c1, c2, dist, rv = nearest_points(ray1, ray2)

    if dist < atol:
        return (c1 + c2) / 2, rv
    else:
        return array((nan, nan, nan)), False


def nearest_points(ray1, ray2):
    """
    Return the nearest points between 2 rays.

    The image point locations in optical systems are usually found by
    calculating the intersection between rays coming from a single object
    point, but in aberrated systems, the 2 rays will not really intersect.
    This function is used to find the point in space where the rays
    are closest to each other. If the rays intersect the values returned
    will be the intersection point.

    The solution was taken from:

        http://homepage.univie.ac.at/Franz.Vesely/notes/hard_sticks/hst/hst.html

    Parameters
    ----------

    r1,r2 : :class:`~pyoptools.raytrace.ray.Ray`
        Rays to test for intersection.

    Returns
    -------

    c1 : tuple(float, float, float)
        Coordinates of the point living on ray 1 closest to ray 2
    c2 : tuple(float, float, float)
        Coordinates of the point living on ray 2 closest to ray 1
    d : float
        The distance between p1 and p2
    rv : bool
        Indicates if the intersection is real or virtual. rv=True for
        real, rv=False for virtual. In this case virtual has the same meaning
        as in virtual image i.e. c1 and c2 are not in the actual path, or are
        behind the ray's origin.
    """
    e1 = ray1.dir
    e2 = ray2.dir
    p1 = ray1.pos
    p2 = ray2.pos

    # the directions are unit vectors
    denominator = 1 - np.dot(e1, e2) ** 2

    # check if the rays are parallel
    if np.isclose(denominator, 0, rtol=1e-8, atol=1e-15, equal_nan=False):
        return np.array((nan, nan, nan)), np.array((nan, nan, nan)), np.NaN, False

    p2p1 = p2 - p1

    e1e2 = np.dot(e1, e2)
    p2p1e1 = np.dot(p2p1, e1)
    p2p1e2 = np.dot(p2p1, e2)
    t1 = (p2p1e1 - p2p1e2 * e1e2) / denominator
    t2 = (-p2p1e2 + p2p1e1 * e1e2) / denominator
    # closest points
    c1 = p1 + t1 * e1  # on ray 1
    c2 = p2 + t2 * e2  # on ray 2

    # log.info("nearest points"+str(p1)+" "+str(p2))
    # log.info("tvalues "+str(t1)+" "+str(t2))

    rv = (t1 >= 0) and (t2 >= 0)

    return c1, c2, np.linalg.norm(c2 - c1), rv


def chief_ray_search(
    opsys,
    ccds,
    o=(0.0, 0.0, 0.0),
    rt=(0.0, 0.0, 0.0),
    er=0.1,
    w=pi / 2.0,
    maxiter=1000,
    wavelength=0.58929,
):
    """
    This function uses a random search algorithm to find the chief_ray for a
    given optical system and object point.

    Parameters
    ----------

    opsys : :class:`pyoptools.raytrace.system.System`
        Optical system that will be used to find the chief ray
    ccds : :class:`pyoptools.raytrace.comp_lib.CCD`
       Detector placed in the aperture plane. Must be centred in the optical
       axis
    o : tuple(float, flo0at, float)
        coordinates of the object point used to find the chief ray
    rt : tuple(float, float, float)
         rotations made to a ray propagating in the z direction to obtain the
         first test ray
    er : float
        Maximum acceptable distance between the ray and the center of the
        aperture
    w : float
        Gaussian width in radians
    wavelength : float
        Wavelength of the ray used to find the principal ray given in
        micrometers (.58929 by default).


    Returns
    -------
    :class:`~pyoptools.raytrace.ray.Ray`
        Chief ray found. (Ray instance)


    Notes
    -----

    The algorithm starts using a given ray, propagating it in the optical
    system, and finding the intersection point of this test ray and the
    aperture plane. The distance from this point and the optical axis is
    recorded.

    Using a gaussian random generator, two rotation angles are calculated,
    to generate a new test ray that is propagated in the optical system,
    and its distance to the optical axis is found at the aperture plane.
    If this distance is less than the distance found for the previous ray,
    this ray is taken as the new *chief ray* candidate, and  the algorithm
    is repeated until the number of iterations reaches *maxiter*, or until
    the distance is less than *er*.

    the *rt* parameter gives the rotations made to a ray originating in
    *o*, and propagating in the *Z* direction, to find the first test ray.

    A detector object *ccds* should be placed at the aperture plane. It is used
    to find the point where the ray intersects the aperture. To increase the
    convergense speed of the algorithm, it is better to make sure that the first
    test ray intersects the detector.



    .. todo::
        Implement a function similar to this one, using a minimization
        algorithm
    """

    # log.info("Entering chief_ray_search function")
    test_ray = Ray(wavelength=wavelength)
    opsys.clear_ray_list()

    btx, bty, btz = rt  # btz is not used
    ntry = 0
    nt = 0
    # Check the initial test ray
    retray = test_ray.ch_coord_sys_inv(o, (btx, bty, 0))
    # log.info("Calculating test_ray")
    opsys.clear_ray_list()
    opsys.reset()
    opsys.ray_add(retray)
    opsys.propagate()
    try:
        x, y, z = ccds.hit_list[0][0]
        dist = sqrt(square(x) + square(y))
    except:
        dist = inf

    p_dist = dist

    while (p_dist > er) and (ntry < maxiter):
        ntry = ntry + 1
        nt = nt + 1
        rx = normal(btx, w)
        ry = normal(bty, w)
        tray = test_ray.ch_coord_sys_inv(o, (rx, ry, 0))

        opsys.clear_ray_list()
        opsys.reset()
        opsys.ray_add(tray)
        opsys.propagate()
        try:
            x, y, z = ccds.hit_list[0][0]
            dist = sqrt(square(x) + square(y))
        except:
            # log.info("CCD not hitted by ray")
            dist = inf

        if p_dist > dist:
            # Select this ray as new generator ray

            btx = rx
            bty = ry
            p_dist = dist
            nt = 0
            retray = tray
            # log.info("distance to aperture center="+str(dist))
        if (nt > 10) and p_dist < inf:
            nt = 0
            w = w / 2
        # limit the minimum value of w
        if w < 0.0000001:
            w = 0.0000001
    # print p_dist,ntry
    return retray


def pupil_location(opsys, ccds, opaxis):
    """
    Function to find the optical system pupils position

    Note
    ----
        For this function to operate, the system should have a rotational
        symmetry around the optical axis.

    Parameters
    ----------
    opsys : :class:`pyoptools.raytrace.system.System`
        Optical system to use.
    opaxis :  :class:`~pyoptools.raytrace.ray.Ray`
        Ray representing the optical axis
    ccds : :class:`pyoptools.raytrace.comp_lib.CCD`
        Surface that represents a detector in the aperture plane

    Returns
    -------
        enpl : tuple(float, float, float)
           (xen,yen,zen) containing the entrance pupil coordinates
        expl : tuple(float, float, float)
           (xex,yex,zex) containing the exit pupil coordinates
    """

    # log.info("Propagate Optical axis ray")
    opsys.clear_ray_list()
    opsys.reset()
    # opsys.ray_add(cray)
    opsys.ray_add(opaxis)

    opsys.propagate()

    if len(ccds.hit_list) == 0:
        raise Exception("The optical axis did not intersect the aperture")
    if len(ccds.hit_list) > 1:
        raise Exception("The optical axis intersected the aperture more than once")

    aip = ccds.hit_list[0][0]
    air = ccds.hit_list[0][1]

    # log.info("Optical Axis Intersection point= "+str(aip))
    # log.info("Intersection Ray= "+str(air))

    # Getting Intersection point in global coordinates

    if len(air.childs) != 1:
        raise Exception("The intersected ray can only have one child")

    ip = air.childs[0].pos
    d = air.childs[0].dir
    # log.info("Intersection point in world coordinates= "+str(ip))
    # log.info("Direction of the optical axis at the intersection point"+str(d))

    # Todo: Check if the optical axis and the aperture are perpendicular

    # Calculate vectors perpendicular to the optical axis and to the XYZ axes
    pv1 = cross(d, (0, 0, 1))
    pv2 = cross(d, (0, 1, 0))
    pv3 = cross(d, (1, 0, 0))

    pv = [pv1, pv2, pv3]

    # Search for the longest pv
    pvn = array((dot(pv1, pv1), dot(pv2, pv2), dot(pv3, pv3)))

    pvm = pv[pvn.argmax()]

    # log.info("Displacement vector found: "+str(pvm))

    # Create ray to calculate the exit pupil
    expuray = air.childs[0].copy()
    expuray.dir = expuray.dir + pvm * 0.0001

    # Create the ray to calculate the entrance pupil
    enpuray = expuray.reverse()

    opsys.clear_ray_list()
    opsys.reset()
    opsys.ray_add(enpuray)
    opsys.ray_add(expuray)

    opsys.propagate()

    enp = enpuray.get_final_rays(inc_zeros=False)
    exp = expuray.get_final_rays(inc_zeros=False)
    oax = opaxis.get_final_rays(inc_zeros=False)
    # log.info("enp="+str(enp))
    # log.info("exp="+str(exp))
    # log.info("oax="+str(oax))
    if len(enp) != 1 or len(exp) != 1 or len(oax) != 1:
        raise Exception(
            "The principal ray or the optical axis ray have more" " than one final ray"
        )
    # log.info("Calculating entrance pupil location")

    # Find the nearest points between the rays.
    # Some times because of numerical errors, or some aberrations in the optical
    # system, the rays do not truly intersect.
    # Use instead the nearest points and issue a warning when the rays do not truly
    # intersect.

    enpl = intersection(opaxis, enp[0])[0]
    if (isnan(enpl)).all():
        p1, p2, d, rv = nearest_points(opaxis, enp[0])
        print(
            "Warning: The optical axis does not intersect the principal ray at the entrance"
        )
        print("pupil. The minimum distance is:", d)
        enpl = (p1 + p2) / 2

    # log.info("Calculating exit pupil location")
    expl = intersection(oax[0], exp[0])[0]
    if (isnan(expl)).all():
        p1, p2, d, rv = nearest_points(oax[0], exp[0])
        print(
            "Warning: The optical axis does not intersect the principal ray at the exit"
        )
        print("pupil. The minimum distance is:", d)
        expl = (p1 + p2) / 2

    return enpl, expl


def paraxial_location(opsys, opaxis):
    """Function to find the paraxial image location

    This function finds the paraxial image location of a point
    located in the optical axis, and a boolean indicating if the image
    is real or virtual (image_location, real_virtual).
    The origin of the opaxis location is taken as the object location

    Parameters
    ----------

    opsys : :class:`~pyoptools.raytrace.system.System`
        Optical system to use.
    opaxis: :class:`~pyoptools.raytrace.ray.Ray`
        Ray representing the optical axis

    Returns
    -------

    image_location : tuple(float, float, float)
        Image location coordinates
    real : bool
        Indicates if the intersection point represent a real image (real=True),
        or a virtual image (real=False).


    Note
    ----

    For this function to operate, the system should have a rotational symmetry
    around the optical axis.
    """

    # log.info("Propagate Optical axis ray")
    opsys.clear_ray_list()
    opsys.reset()
    # opsys.ray_add(cray)
    opsys.ray_add(opaxis)

    opsys.propagate()

    # Calculate vectors perpendicular to the optical axis and to the XYZ axes
    d = opaxis.dir
    pv1 = cross(d, (0, 0, 1))
    pv2 = cross(d, (0, 1, 0))
    pv3 = cross(d, (1, 0, 0))

    pv = [pv1, pv2, pv3]

    # Search for the longest pv
    pvn = array((dot(pv1, pv1), dot(pv2, pv2), dot(pv3, pv3)))

    pvm = pv[pvn.argmax()]

    # log.info("Displacement vector found: "+str(pvm))

    # Create paraxial ray

    par_ray = opaxis.copy()
    par_ray.dir = par_ray.dir + pvm * 0.001

    opsys.clear_ray_list()
    opsys.reset()
    opsys.ray_add(par_ray)
    opsys.propagate()

    par = par_ray.get_final_rays(inc_zeros=False)
    oax = opaxis.get_final_rays(inc_zeros=False)
    # log.info("par="+str(par))
    # log.info("oax="+str(oax))

    if len(par) != 1 or len(oax) != 1:
        raise Exception(
            "The paraxial ray or the optical axis ray have more" " than one final ray"
        )

    # log.info("Calculating object location")
    expl = intersection(oax[0], par[0])
    return expl


def find_aperture(ccd, size=(50, 50)):
    """Function to find a mask representing the aperture

    This function returns a array containing 1's and 0's representing
    the aperture shape. The aperture shape will be approximated from
    the CCD hit_list

    Parameters
    ----------

    ccd : :class:`~pyoptools.raytrace.comp_lib.CCD`
        CCD object that will be used to get the shape information from
    size : tuple(int, int)
        Array shape


    Returns
    -------
    array
        Array with the image of the aperture


    .. todo::
        please describe better

    Notes
    -----

    Right now only works for round apertures.

    .. todo::
        please be more specific
    """

    hl = ccd.hit_list
    sx, sy = ccd.size
    tx, ty = size
    dx, dy = sx / (tx - 1), sy / (ty - 1)
    CG = mgrid[
        float(-sx / 2.0) : float(sx / 2.0 + dx) : float(dx),
        float(-sy / 2.0) : float(sy / 2.0 + dy) : float(dy),
    ]

    rm = sqrt(CG[0] ** 2 + CG[1] ** 2)

    maxr = 0.0
    for i in hl:
        X, Y, Z = i[0]
        r = sqrt(X * X + Y * Y)
        if maxr < r:
            maxr = r

    return where(rm < maxr, 1.0, 0.0)


def find_ppp(opsys, opaxis):
    """Function to find the primary principal plane location of a lens or an
    optical component

    Parameters
    ----------

    opsys : :class:`~pyoptools.raytrace.system.System`
        Optical system or optical component whose principal planes are to be
        found
    opaxis : :class:`~pyoptools.raytrace.ray.Ray`
        Ray defining the optical axis of the system


    For this function to operate, the system should have a rotational symmetry
    around the optical axis.

    Notes
    -----
        This function is returns the intersection point of the optical axis and
        the principal plane.
    """

    # Create a system with the component
    if isinstance(opsys, (Component)):
        c = opsys
        opsys = System(
            complist=[
                (c, (0, 0, 0), (0, 0, 0)),
            ],
            n=1,
        )

    # To create a ray parallel to the optical axis, find a displacement vector
    # perpendicular to the optical axis, and to the XYZ axes

    d = opaxis.dir
    pv1 = cross(d, (0, 0, 1))
    pv2 = cross(d, (0, 1, 0))
    pv3 = cross(d, (1, 0, 0))

    pv = [pv1, pv2, pv3]

    # Search for the longest pv
    pvn = array((dot(pv1, pv1), dot(pv2, pv2), dot(pv3, pv3)))

    pvm = pv[pvn.argmax()]

    # Create parallel ray

    par_ray = opaxis.copy()
    par_ray.pos = par_ray.pos + pvm * 0.0001

    opsys.clear_ray_list()
    opsys.ray_add([opaxis, par_ray])
    opsys.propagate()
    par_ray_end = par_ray.get_final_rays(inc_zeros=False)

    if len(par_ray_end) != 1:
        raise Exception("The paraxial ray has more than one final ray")

    pppl = intersection(par_ray, par_ray_end[0])

    # Move the intersection point toward the optical axis

    ppp = pppl[0] - pvm * 0.0001
    return ppp  # , pppl[1])


def get_optical_path_ep(opsys, opaxis, raylist, stop=None, r=None):
    """Returns the optical path traveled by a ray up to the exit pupil

    The optical path is measured from the ray origin until it crosses the
    exit pupil of the system.
    If a stop (aperture) is not given, the measurement is made up to the primary
    principal plane.


    Parameters
    ----------
    opsys : :class:`~pyoptools.raytrace.system.System`
        Optical system under analysis
    opaxis : :class:`pyoptools.raytrace.ray.Ray`
        Ray indicating the optical axis. The origin of the optical axis must be
        the position of the object used in the image formation. This is needed
        to be able to calculate the radius of the reference sphere.

    raylist: list(:class:`pyoptools.raytrace.ray.Ray`)
        List of rays that will be used to sample the optical path

    stop : :class:`~pyoptools.raytrace.comp_lib.Stop`
        Stop aperture of the system. It must belong to opsys. In not given it
        will be assumed that the exit pupil is at the primary principal plane.
    r :
        If None, measure up to the exit pupil plane. If given, use a reference
        sphere with a vertex coinciding with the optical vertex.
        .. todo::
            Need to check the function and fix this documentation

    Returns
    -------

    hcl : list
        List containing the coordinates of the hits in the pupil coordinate
        system.
    opl : list
        list containing the optical paths measured
    pc : tuple(float, float, float)
        intersection point between the optical axis, and the pupil plane.


    hcl[i] corresponds to opl[i]

    Note: This method only works if the optical axis coincides with the Z axis.
    This must be corrected.
    """
    if stop != None:
        enp, exp = pupil_location(opsys, stop, opaxis)
    else:
        exp = find_ppp(opsys, opaxis)

    # Reset the system
    opsys.clear_ray_list()
    opsys.reset()

    # Propagate the rays
    # print "***", raylist
    opsys.ray_add(raylist)
    opsys.propagate()
    # pf=PlotFrame(opsys=opsys)
    rl = []
    l = []

    # Get the optical path up to the final element in the system
    for i in raylist:
        a = i.get_final_rays()
        if a[0].intensity != 0:
            # Reverse the rays to calculate the optical path from the final element
            # to the exit pupil

            nray = a[0].reverse()
            rl.append(nray)
            # TODO: This should not be done using the label
            nray.label = str(a[0].optical_path_parent())

    # Create a dummy system to calculate the wavefront at the exit pupil
    if r == None:
        # TODO: This ccd should be infinitely big. Have to see how this can be done
        ccd = CCD(size=(1000, 1000))
    else:
        ccds = Spherical(shape=Circular(radius=0.9 * r), curvature=1.0 / r)
        ccd = Component(
            surflist=[
                (ccds, (0, 0, 0), (0, 0, 0)),
            ]
        )
    # print rl

    dummy = System(
        complist=[
            (ccd, exp, (0, 0, 0)),
        ],
        n=1.0,
    )

    # Calculate the optical path from the final element to the exit pupil plane
    dummy.ray_add(rl)
    dummy.propagate()
    # PlotFrame(opsys=dummy)
    hcl = []
    opl = []
    for ip, r in ccd.hit_list:
        # print ip
        x, y, z = ip
        # TODO: This should not be done using the label
        d = float(r.label) - r.optical_path()
        hcl.append((x, y, z))
        opl.append(d)
    return (hcl, opl, exp)

    # rv=bisplrep(X,Y,Z)
    # data=bisplev(array(range(-20,20)),array(range(-20,20)),rv)

    # data=(data-data.mean())

    # print "Gaussian reference sphere radius =",sqrt(dot(impos-exp,impos-exp))


def find_reference_sphere_radius(ip, pl):
    """Find the radius os the reference sphere that best fits the input data.

    This method assumes that the optical axis coincides with the z axis. This
    means that the center of the sphere, has coordinates (0,0,r).

    Parameters
    ----------

    ip : list
        list of the points where the optical path is measured, that are being
        fitted. Each point is (XYZ) tuple. It can be also an array with a shape
        n,3 where n is the number of points.
    pl : list
        List of path lengths. pl[i] corresponds to the point ip[i].

    Returns
    -------
    float
        Reference sphere radius
    """

    ipa = array(ip)
    pla = array(pl)
    n, t = ipa.shape

    # Find the point closest to the center of the aperture.
    rm = sqrt(dot(ipa[0], ipa[0]))
    im = 0
    for i in range(n):
        if rm > sqrt(dot(ipa[i], ipa[i])):
            rm = sqrt(dot(ipa[i], ipa[i]))
            im = i

    # Make the OPL 0 at the center of the aperture
    pla = pla - pla[im]

    # Encontrar el radio de la esfera de mejor ajuste
    def F(z):
        dist = pla - (sqrt(ipa[:, 0] ** 2 + ipa[:, 1] ** 2 + (ipa[:, 2] - z) ** 2) - z)
        u = sqrt((dist ** 2).sum())
        # print "*", u
        # u=dist[-1]
        # print u
        return u

    r = fsolve(F, -10.0)

    return r


def aux_paral_f(x):
    """
    Auxiliary function needed in parallel propagate
    """

    os, rb = x
    os.ray_add(rb)
    os.propagate()
    return os


def parallel_propagate(os, r, np=None):
    """Perform a propagation of the rays in the system using all cores
    present on a computer

    os gets reset before beginning the propagation, so the only rays
    used in the simulation are the rays given in r

    Parameters
    ----------
    os : :class:`~pyoptools.raytrace.system.System`
        Optical system used in the simulation
    r : list(:class:`pyoptools.raytrace.ray.Ray`)
        List containing the rays to propagate
    np : int or None
        Number of processes used in the simulation. If not given use one
        process per cpu
    """

    if np == None:
        cpus = mp.cpu_count()
    else:
        cpus = np
    pool = mp.Pool(cpus)
    os.reset()
    # Split the ray list in the number of CPUS
    nr = len(r)
    r_list = []
    r_list.append((os, r[: nr / cpus]))
    for i in range(2, cpus):
        r_list.append((os, r[(nr / cpus) * (i - 1) : (nr / cpus) * (i)]))
    r_list.append((os, r[(nr / cpus) * (cpus - 1) :]))
    osi = pool.map(aux_paral_f, r_list)

    pool.close()
    pool.join()

    for osp in osi:
        os.merge(osp)
    return os


def aux_paral_f_ns(x):
    """
    Auxiliary function needed in parallel propagate
    """

    # os optical system
    # rg guide ray
    # dp Path (key) of the destination surface.
    # rb rays to propagate
    os, rg, dp, rb = x
    os.ray_add(rb)
    os.propagate_ray_ns(rg, dp)
    return os


def parallel_propagate_ns(os, rg, dp, r, np=None):
    """Perform a propagation of the rays in the system using all cores
    present on a computer

    os gets reset before beginning the propagation, so the only rays
    used in the simulation are the rays given in r

    Parameters
    ----------
    os :
        Optical system used in the simulation
    rg :
        Guide ray
    dp :
        Destination path
    r :
        List containing the rays to propagate
    np : int or None
        Number if processes used in the simulation. If not given use one
        process per cpu
    """

    if np == None:
        cpus = mp.cpu_count()
    else:
        cpus = np
    pool = mp.Pool(cpus)
    os.reset()
    # Split the ray list in the number of CPUS
    nr = len(r)
    r_list = []
    r_list.append((os, rg, dp, r[: nr / cpus]))
    for i in range(2, cpus):
        # os,rg,dp,rb=x
        r_list.append((os, rg, dp, r[(nr / cpus) * (i - 1) : (nr / cpus) * (i)]))
    r_list.append((os, rg, dp, r[(nr / cpus) * (cpus - 1) :]))
    osi = pool.map(aux_paral_f_ns, r_list)

    pool.close()
    pool.join()

    for osp in osi:
        os.merge(osp)
    return os


def ray_paths(r):
    """
    Return lists with all the possible paths traveled by the ray r.

    r must be previously propagated in an optical system

    When there are beam splitters, there is more than one path

    Parameters
    ----------
        r : :class:`pyoptools.raytrace.ray.Ray`


    .. todo::
        Finish documentation
    """

    def rt(r):
        l = []
        rays = r.childs
        for ray in rays:
            a = rt(ray)
            for ray1 in a:
                l.append([ray] + ray1)
            if len(a) == 0:
                l.append([ray])
        return l

    A = rt(r)
    B = []
    for rp in A:
        t = [r] + rp
        B.append(t)
    return B
