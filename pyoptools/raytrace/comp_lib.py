"""Helper module with definition of standard optical components

The classes imported in this module are defined in the _comp_lib module,
and are used to create the surfaces and the components needed to define
some standard optical components, and return an instance to such components.
This module hides from the end user un needed stuff.
"""

# sorted is used so the classes are ordered alphabetically in the sphinx docs
__all__ = sorted(
    [
        "SphericalLens",
        "AsphericLens",
        "CylindricalLens",
        "CCD",
        "RightAnglePrism",
        "PentaPrism",
        "DovePrism",
        "Block",
        "BeamSplitingCube",
        "BeamSplittingCube",
        "Doublet",
        "AirSpacedDoublet",
        "MultiLens",
        "Stop",
        "IdealLens",
        "IdealTLens",
        "RoundMirror",
        "RectMirror",
        "RectGratting",
        "PowellLens",
        "SimpleDMDDevice",
    ]
)

from ._comp_lib.aspheric_lens import AsphericLens
from ._comp_lib.ccd import CCD
from ._comp_lib.compound_lens import AirSpacedDoublet, Doublet, MultiLens
from ._comp_lib.cube import BeamSplitingCube, BeamSplittingCube, Block
from ._comp_lib.cylindrical_lens import CylindricalLens
from ._comp_lib.diffraction import RectGratting
from ._comp_lib.ideallens import IdealLens, IdealTLens
from ._comp_lib.mirror import RectMirror, RoundMirror
from ._comp_lib.powell_lens import PowellLens
from ._comp_lib.prism import DovePrism, PentaPrism, RightAnglePrism
from ._comp_lib.simple_dmd_device import SimpleDMDDevice
from ._comp_lib.spherical_lens import SphericalLens
from ._comp_lib.stop import Stop
