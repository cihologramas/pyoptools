"""
Definition of a radially-symmetric aspheric lens.
"""

from math import pi

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Aspherical, Cylinder, Plane, Aperture
from pyoptools.raytrace.shape import Circular
from pyoptools.misc.function_2d.poly_r.poly_r import PolyR

from types import SimpleNamespace


class AsphericLens(Component):
    """Helper class to define radially-symmetric aspherical lenses.

    Both single sided and double sided aspherics are supported.
    Lenses with a outer brim larger than the aspheric surface are supported.

    Parameters
    ----------

    thickness : float
                The total thickness in mm of the lens at the maximum
    outer_diameter : float, optional
                     The outer diameter of the lens in mm.
                     If unspecified or None, will be found from the largest
                     of either aspheric surface.
    material : Material or float
               Material defining the refraction index of the lens.
               Can be a material instance or floating point number if the
               refraction index is constant.
    origin : str
             Where along the optical axis to place the origin of the
             coordinate system. By default will the the geometric center. Can
             also be:

             * 's1_max' : Point of maximum thickness on s1
             * 's1_min' : Point of minimum thickness on s1
             * 's2_max' : Point of maximum thickness on s2
             * 's2_min' : Point of minimum thickness on s2
             * 'center' : Geometric center

             These options can be convenient for placing the origin
             at mounting face.
    s1 : dict
         Radially symmetric aspheric surface definition dict
         for the anterior surface, as described below.
    s2 : dict, optional
         Either None for a plano-aspheric lens or a symmetric
         aspheric surface definition dict for the posterior
         surface, as described below.

    Surface Definitions
    -------------------

    diameter : float
               Diameter of the aspheric surface in mm
    roc  : float
           Radius of curvature parameter for the surface in mm
    k : float
        Conic constant
    polycoefficents : tuple of float
                    Tuple listing the higher order aspheric coefficients
                    The first element corresponds to index zero.
                    Typically, the elements used in a surface description
                    are even numbered, starting at index 4, so the first
                    four elements are typically zero.
    max_thickness : float, optional
                    Maximum thickness of the surface in mm or None.
                    If None, the maximum thickness will be found from the
                    thickness at zero radius, however not all possible
                    surfaces have the maximum here.
    """

    def __init__(
        self,
        outer_diameter=8.0,
        thickness=3.0,
        material=1.5,
        origin="center",
        s1={
            "diameter": 6.0,
            "roc": 3.0,
            "k": -1.5,
            "polycoefficents": (0, 0, 0, 0, 3e-3, 0, -10e-6),
            "max_thickness": None,
        },
        s2=None,
        *args,
        **kwargs
    ):
        Component.__init__(self, *args, **kwargs)

        self.material = material
        self.thickness = thickness
        self.outer_diameter = outer_diameter

        # Fill defaults and put surface definitions into namespaces
        if "max_thickness" not in s1:
            s1["max_thickness"] = None
        if s2 is not None and "max_thickness" not in s2:
            s2["max_thickness"] = None

        s1_defn = SimpleNamespace(**s1)
        if s2 is not None:
            s2_defn = SimpleNamespace(**s2)
        else:
            s2_defn = None
        self.s1_defn, self.s2_defn = (s1_defn, s2_defn)

        # Auto select outer diameter if None
        if self.outer_diameter is None:
            candidates = [s1_defn.diameter]
            if s2 is not None:
                candidates.append(s2_defn.diameter)
            self.outer_diameter = max(candidates)

        # Start side thickness calculation, surfaces will be subtracted from this
        side_thickness = thickness

        # First surface
        s1_surf = Aspherical(
            shape=Circular(radius=0.5 * s1_defn.diameter),
            Ax=1.0 / s1_defn.roc,
            Ay=1.0 / s1_defn.roc,
            Kx=s1_defn.k,
            Ky=s1_defn.k,
            poly=PolyR(s1_defn.polycoefficents),
        )
        if s1_defn.max_thickness is None:
            s1_defn.max_thickness = s1_surf.get_z_at_point(s1_defn.diameter / 2.0, 0)
        side_thickness -= s1_defn.max_thickness
        self.surflist.append((s1_surf, (0, 0, 0), (0, 0, pi / 2)))

        # Second surface
        if s2_defn is None:
            s2_surf = Plane(shape=Circular(radius=0.5 * self.outer_diameter))
        else:
            s2_surf = Aspherical(
                shape=Circular(radius=0.5 * s2_defn.diameter),
                Ax=1.0 / s2_defn.roc,
                Ay=1.0 / s2_defn.roc,
                Kx=s2_defn.k,
                Ky=s2_defn.k,
                poly=PolyR(s2_defn.polycoefficents),
            )
            if s2_defn.max_thickness is None:
                s2_defn.max_thickness = s2_surf.get_z_at_point(s2_defn.diameter / 2.0, 0)
            side_thickness -= s2_defn.max_thickness
        self.surflist.append((s2_surf, (0, 0, thickness), (0, pi, pi / 2)))

        # Outer edge
        if side_thickness > 0:
            outer_edge = Cylinder(radius=self.outer_diameter / 2, length=side_thickness)
            self.surflist.append(
                (
                    outer_edge,
                    (0, 0, s1_defn.max_thickness + 0.5 * side_thickness),
                    (0, 0, pi / 2),
                )
            )
        elif side_thickness < 0:
            raise InvalidGeometryException(
                "Lens is not thick enough to support surfaces."
            )

        # Brim
        self._add_brim(s1_defn, s1_defn.max_thickness)
        if s2_defn is not None:
            self._add_brim(s2_defn, thickness - s2_defn.max_thickness)

        # Apply offset
        self._translate_origin(origin)

    def _add_brim(self, defn, z_position):
        """Adds an outer brim surface defined by surface definition defn.
        Only applies is the specified outer diameter of this lens is larger
        than the maximum diameter of the defined surface.
        Brim will be located at position z_position.
        """

        if self.outer_diameter > defn.diameter:
            brim = Aperture(
                shape=Circular(radius=0.5 * self.outer_diameter),
                ap_shape=Circular(radius=0.5 * defn.diameter),
            )
            self.surflist.append((brim, (0, 0, z_position), (0, 0, pi / 2)))
        elif self.outer_diameter < defn.diameter:
            raise InvalidGeometryException(
                "Lens outer diameter can not be smaller than the diameter of an aspheric surface."
            )

    def _translate_origin(self, origin):
        """Move all elements so that the origin coincides with the origin
        descriptor string, as documented in constructor.
        """
        # offset along the optical axis
        selector = {
            "s1_min": -1 * self.s1_defn.max_thickness,
            "s1_max": 0,
            "center": -0.5 * self.thickness,
        }
        if self.s2_defn is not None:
            selector["s2_min"] = (-1 * self.thickness + self.s2_defn.max_thickness,)
            selector["s2_max"] = (-1 * self.thickness,)

        try:
            delta = selector[origin]
        except KeyError:
            raise ValueError("Invalid origin offset string.")

        for k, v in self.surflist.items():
            surface, translation, rotation = v
            tx, ty, tz = translation
            self.surflist[k] = (surface, (tx, ty, tz + delta), rotation)


# TODO: move to a dedicated module
class InvalidGeometryException(Exception):
    pass
