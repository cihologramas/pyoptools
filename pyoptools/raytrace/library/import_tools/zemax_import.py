import numpy as np
from struct import Struct
from pathlib import Path
from enum import Enum, auto
import json
import sys
import csv
import textwrap
import traceback
import argparse
import re

from pyoptools.raytrace.mat_lib import material


def convert(d):
    try:
        return int(d)
    except ValueError:
        try:
            return float(d)
        except ValueError:
            return d


def find_key(key, lines):
    try:
        rv = ""
        for s in list(filter(lambda x: x.startswith(key), lines)):
            rv = rv + s[len(key) + 1:]
    except IndexError:
        # key not found
        rv = ""
    return rv


def checktype(surflist, t):
    "Returns true if any surfaces in list have type id string t"
    for s in surflist:
        if "TYPE" in s and s["TYPE"][0] == t:
            return True
    return False


def checkglas(surflist, t):
    for s in surflist:
        if (s.get("GLAS", [None])[0]) == t:
            return True
    return False


def zmf_decode(data, a, b):
    iv = np.cos(6 * a + 3 * b)
    iv = np.cos(655 * (np.pi / 180) * iv) + iv
    p = np.arange(len(data))
    k = 13.2 * (iv + np.sin(17 * (p + 3))) * (p + 1)
    k = (int(("{:.8e}".format(_))[4:7]) for _ in k)
    # data = np.fromstring(data, np.uint8)
    data = np.frombuffer(data, np.uint8).copy()
    data ^= np.fromiter(k, np.uint8, len(data))
    return data.tobytes()


def flip_aspheric_defn(defn):
    "Flip the definition dict for an aspheric surface"
    return {
        "diameter": defn["diameter"],
        "roc": -1 * defn["roc"],
        "k": 1 * defn["k"],
        "polycoefficents": tuple([-1 * x for x in defn["polycoefficents"]]),
    }


class FailedImport(Enum):
    manual_exclusion = auto()
    mode_not_seq = auto()
    units_not_mm = auto()
    mirrors_unsupported = auto()
    aspheres_unsupported = auto()
    doublet_pairs_unsupported = auto()
    diffractive_optics_unsupported = auto()
    conical_unsupported = auto()
    fresnet_unsupported = auto()
    acylindrical_unsupported = auto()
    no_primary_surface = auto()
    diameter_undefined = auto()
    unequal_surface_radius = auto()
    unknown_material_type = auto()
    glass_undefined = auto()
    missing_aspheric_surface = auto()
    unknown = auto()


class OpticImporter:
    def __init__(self, description="", surflist=[], gcat="", key=""):
        self.description = description
        self.key = key
        self.surflist = surflist
        self.lens_data = {"description": self.description, "glass_catalogs": gcat}


class SingletImporter(OpticImporter):
    def valid(self):
        return len(self.surflist) == 2

    def definition(self):

        c0 = self.surflist[0]["CURV"][0]
        c1 = self.surflist[1]["CURV"][0]
        d0 = self.surflist[0]["DISZ"][0]
        # Not all lenses have DIAM in both surfaces
        if "DIAM" in self.surflist[0] and "DIAM" in self.surflist[1]:
            r0 = self.surflist[0]["DIAM"][0]
            r1 = self.surflist[1]["DIAM"][0]
        elif "DIAM" in self.surflist[0]:
            r0 = self.surflist[0]["DIAM"][0]
            r1 = r0
        elif "DIAM" in self.surflist[1]:
            r0 = self.surflist[1]["DIAM"][0]
            r1 = r0
        else:
            return FailedImport.diameter_undefined

        g0 = self.surflist[0]["GLAS"][0]

        # Check that the surfaces are equal, if not issue an error
        if not r0 == r1:
            return FailedImport.unequal_surface_radius

        self.lens_data["material"] = g0
        self.lens_data["thickness"] = self.surflist[0]["DISZ"][0]
        self.lens_data["radius"] = r0
        self.lens_data["curvature_s2"] = self.surflist[1]["CURV"][0]
        self.lens_data["curvature_s1"] = self.surflist[0]["CURV"][0]
        self.lens_data["type"] = "SphericalLens"

        # Ball and half-ball lenses often have a numerical issue where
        # the aperture radius is slightly larger than allowed by the
        # surface curvature. This handles that special case
        if "ball" in self.lens_data["description"].lower():
            max_curve = max(
                abs(self.lens_data["curvature_s1"]), abs(self.lens_data["curvature_s2"])
            )
            if r0 * max_curve >= 1:
                self.lens_data["radius"] = 1.0 / max_curve
                print(
                    f"Fixed radius of item {self.key} : "
                    f"r0={r0} to {self.lens_data['radius']}"
                )

        return self.lens_data


class CylindricalImporter(OpticImporter):
    """Cylindrical lenses. It seems these are modeled as 'toroidal'
    Currently working for plano-convex and plano-concave, rectangular
    """

    def valid(self):
        return checktype(self.surflist, "TOROIDAL")

    def definition(self):

        self.lens_data["type"] = "CylindricalLens"

        # Get the size from the first surface
        all_sqaps = [tuple(s["SQAP"]) for s in self.surflist if "SQAP" in s]
        # Verify all SQAP elements same for all surfaces
        if len(set(all_sqaps)) != 1:
            return FailedImport.acylindrical_unsupported

        self.lens_data["size"] = tuple(reversed(all_sqaps[0][0:2]))

        # Get thickness. Seems the DISZ element in the surface with
        # COMM equal to the key is where this is found
        k = self.key
        surf = [s for s in self.surflist if "COMM" in s and s["COMM"][0] == k]
        if len(surf) != 1:
            return FailedImport.no_primary_surface

        self.lens_data["thickness"] = float(surf[0]["DISZ"][0])
        self.lens_data["curvature_s1"] = float(surf[0]["CURV"][0])
        self.lens_data["curvature_s2"] = 0.0
        self.lens_data["material"] = surf[0]["GLAS"][0]

        return self.lens_data


class AsphericImporter(OpticImporter):
    def valid(self):
        return checktype(self.surflist, "EVENASPH")

    def definition(self):
        # print('Importing aspheric')

        self.lens_data["type"] = "AsphericLens"
        self.lens_data["thickness"] = None
        self.lens_data["material"] = None
        self.lens_data["origin"] = "center"

        # find the maximum clear aperture
        max_clear_aperture = 0
        for s in self.surflist:
            if "MEMA" in s:
                max_clear_aperture = max(max_clear_aperture, s["MEMA"][0] * 2)
        if max_clear_aperture == 0:
            max_clear_aperture = None
        self.lens_data["max_clear_aperture"] = max_clear_aperture

        # find the diameter if given in NAME field
        match = re.search("Ø=(.*?)mm", self.description)
        if match is not None:
            diameter = float(match.group(0)[2:-2])
        else:
            diameter = None

        # find the aspheric surfaces, may be one or possibly two
        aspheric_surface_defs = []
        first_aspheric_index = None
        for i, s in enumerate(self.surflist):
            if "TYPE" in s and s["TYPE"][0] == "EVENASPH" and s["CURV"][0] != 0.0:

                if first_aspheric_index is None:
                    first_aspheric_index = i

                surface_def = {}

                # It is not trivial to find the surface diameter

                # if diameter directly specified for the surface, use that
                if "DIAM" in s:
                    surface_def["diameter"] = 2 * s["DIAM"][0]

                # otherwise, if clear aperture specified use that
                elif "MEMA" in s:
                    surface_def["diameter"] = 2 * s["MEMA"][0]

                # otherwise, if diameter specified globally, use that
                elif diameter is not None:
                    surface_def["diameter"] = diameter

                # otherwise, use the max clear aperture from all surfaces
                elif max_clear_aperture is not None:
                    surface_def["diameter"] = max_clear_aperture

                else:
                    return FailedImport.diameter_undefined

                surface_def["roc"] = 1.0 / s["CURV"][0]

                if "CONI" in s:
                    surface_def["k"] = s["CONI"][0]
                else:
                    surface_def["k"] = 0

                # Put the coefficents in the right place in polycoefficents
                coefficents = [0] * 17
                param_idx_to_coefficent = {2: 4, 3: 6, 4: 8, 5: 10, 6: 12, 7: 14, 8: 16}
                for p, value in s["PARM"].items():
                    if p in param_idx_to_coefficent:
                        coef_index = param_idx_to_coefficent[p]
                        coefficents[coef_index] = value

                if all(c == 0 for c in coefficents):
                    coefficents = [0, 0, 0, 0]
                else:
                    # Remove any trailing zeros
                    while coefficents[-1] == 0:
                        coefficents.pop(-1)

                surface_def["polycoefficents"] = coefficents

                aspheric_surface_defs.append(surface_def)

                # so far as I can tell, DISZ of first surface always
                # the total thickness
                if self.lens_data["thickness"] is None:
                    self.lens_data["thickness"] = s["DISZ"][0]

                if self.lens_data["material"] is None:
                    try:
                        self.lens_data["material"] = s["GLAS"][0]
                    except KeyError:
                        return FailedImport.glass_undefined

        try:
            self.lens_data["s1"] = aspheric_surface_defs[0]
        except IndexError:
            return FailedImport.missing_aspheric_surface

        if len(aspheric_surface_defs) > 1:
            # Double-aspheric lens
            self.lens_data["s2"] = flip_aspheric_defn(aspheric_surface_defs[1])
        else:
            # If only one aspheric surface found, decide what to do based
            # on the next defined surface
            posterior_surface = self.surflist[first_aspheric_index + 1]
            s1_data = self.lens_data["s1"]

            if posterior_surface["CURV"][0] == 0:
                # plano surface
                self.lens_data["s2"] = None
            else:
                # spherical surface
                self.lens_data["s2"] = {
                    "diameter": s1_data["diameter"],
                    "roc": -1.0 / posterior_surface["CURV"][0],
                    "k": 0,
                    "polycoefficents": (0, 0, 0, 0, 0, 0, 0),
                }

        # put in the diameter
        if diameter is None:
            if max_clear_aperture is None:
                # Default to s1 diameter if nothing known
                self.lens_data["outer_diameter"] = self.lens_data["s1"]["diameter"]
            else:
                self.lens_data["outer_diameter"] = max_clear_aperture
        else:
            # special case : both known but clear aperture largest, use that
            if max_clear_aperture is not None and max_clear_aperture > diameter:
                self.lens_data["outer_diameter"] = max_clear_aperture
            # in general, use the specified diameter
            else:
                self.lens_data["outer_diameter"] = diameter

        return self.lens_data


class DoubletImporter(OpticImporter):
    def valid(self):
        return len(self.surflist) == 3

    def definition(self):
        # Doublets
        self.lens_data["curvature_s1"] = self.surflist[0]["CURV"][0]
        self.lens_data["curvature_s2"] = self.surflist[1]["CURV"][0]
        self.lens_data["curvature_s3"] = self.surflist[2]["CURV"][0]
        self.lens_data["thickness_l1"] = self.surflist[0]["DISZ"][0]
        self.lens_data["thickness_l2"] = self.surflist[1]["DISZ"][0]
        r0 = self.surflist[0]["DIAM"][0]
        r1 = self.surflist[1]["DIAM"][0]
        r2 = self.surflist[2]["DIAM"][0]

        # Check that the surfaces are equal, if not issue an error
        if not (r0 == r1 and r1 == r2):
            return FailedImport.unequal_surface_radius

        self.lens_data["radius"] = r0

        self.lens_data["material_l1"] = self.surflist[0]["GLAS"][0]
        self.lens_data["material_l2"] = self.surflist[1]["GLAS"][0]
        self.lens_data["type"] = "Doublet"
        return self.lens_data


class AirSpacedDoubletImporter(OpticImporter):

    def valid(self):
        return len(self.surflist) == 4

    def definition(self):
        self.lens_data["curvature_s1"] = self.surflist[0]["CURV"][0]
        self.lens_data["curvature_s2"] = self.surflist[1]["CURV"][0]
        self.lens_data["curvature_s3"] = self.surflist[2]["CURV"][0]
        self.lens_data["curvature_s4"] = self.surflist[3]["CURV"][0]

        self.lens_data["thickness_l1"] = self.surflist[0]["DISZ"][0]
        self.lens_data["air_gap"] = self.surflist[1]["DISZ"][0]
        self.lens_data["thickness_l2"] = self.surflist[2]["DISZ"][0]

        # Verificar que las superficies son iguales, si no emitir un error
        # assert r0==r1 and r1== r2 and r2 ==r3 Esto no siempre se cumple

        r0 = self.surflist[0]["DIAM"][0]
        r1 = self.surflist[1]["DIAM"][0]
        r2 = self.surflist[2]["DIAM"][0]
        r3 = self.surflist[3]["DIAM"][0]

        # Verificar que las superficies son iguales, si no emitir un error
        # assert r0==r1 and r1== r2 and r2 ==r3 #Esto no siempre se cumple

        self.lens_data["radius"] = r0

        self.lens_data["material_l1"] = self.surflist[0]["GLAS"][0]
        # g1=surflist[1]["GLAS"][0] Este parametro no existe. Es aire
        self.lens_data["material_l2"] = self.surflist[2]["GLAS"][0]

        self.lens_data["type"] = "AirSpacedDoublet"
        return self.lens_data


class DoubletPairImporter(OpticImporter):
    """Doublet pair
    Thorlabs have the term 'Matched Achromatic Pair'
    in the description. Check for this to avoid
    collision with triplets
    """

    def valid(self):
        return len(self.surflist) == 6 and "Pair" in self.description

    def definition(self):

        self.lens_data["curvature_s1"] = self.surflist[0]["CURV"][0]
        self.lens_data["curvature_s2"] = self.surflist[1]["CURV"][0]
        self.lens_data["curvature_s3"] = self.surflist[2]["CURV"][0]
        self.lens_data["curvature_s4"] = self.surflist[3]["CURV"][0]
        self.lens_data["curvature_s5"] = self.surflist[4]["CURV"][0]
        self.lens_data["curvature_s6"] = self.surflist[5]["CURV"][0]

        self.lens_data["thickness_l1"] = self.surflist[0]["DISZ"][0]
        self.lens_data["thickness_l2"] = self.surflist[1]["DISZ"][0]
        self.lens_data["air_gap"] = self.surflist[2]["DISZ"][0]
        self.lens_data["thickness_l3"] = self.surflist[3]["DISZ"][0]
        self.lens_data["thickness_l4"] = self.surflist[4]["DISZ"][0]

        self.lens_data["material_l1"] = self.surflist[0]["GLAS"][0]
        self.lens_data["material_l2"] = self.surflist[1]["GLAS"][0]
        self.lens_data["material_l3"] = self.surflist[3]["GLAS"][0]
        self.lens_data["material_l4"] = self.surflist[4]["GLAS"][0]

        r0 = self.surflist[0]["DIAM"][0]
        r1 = self.surflist[1]["DIAM"][0]
        r2 = self.surflist[2]["DIAM"][0]
        r3 = self.surflist[3]["DIAM"][0]
        r4 = self.surflist[4]["DIAM"][0]
        r5 = self.surflist[5]["DIAM"][0]

        if not (r0 == r1 and r1 == r2 and r2 == r3 and r3 == r4 and r4 == r5):
            return FailedImport.unequal_surface_radius

        self.lens_data["radius"] = r0
        self.lens_data["type"] = "DoubletPair"

        # return lens_data
        # There isn't yet a component class which can load these
        return FailedImport.doublet_pairs_unsupported


class DoubletPairCentStopImporter(OpticImporter):
    """Doublet pair with stop in the middle, stop ignored
    Thorlabs have the term 'Matched Achromatic Pair'
    in the description. Check for this to avoid
    collision with triplets.
    """

    def valid(self):
        return len(self.surflist) == 7 and "Pair" in self.description

    def definition(self):
        self.lens_data["curvature_s1"] = self.surflist[0]["CURV"][0]
        self.lens_data["curvature_s2"] = self.surflist[1]["CURV"][0]
        self.lens_data["curvature_s3"] = self.surflist[2]["CURV"][0]
        self.lens_data["curvature_s4"] = self.surflist[4]["CURV"][0]
        self.lens_data["curvature_s5"] = self.surflist[5]["CURV"][0]
        self.lens_data["curvature_s6"] = self.surflist[6]["CURV"][0]

        self.lens_data["thickness_l1"] = self.surflist[0]["DISZ"][0]
        self.lens_data["thickness_l2"] = self.surflist[1]["DISZ"][0]
        self.lens_data["air_gap"] = (
            self.surflist[2]["DISZ"][0] + self.surflist[3]["DISZ"][0]
        )
        self.lens_data["thickness_l3"] = self.surflist[4]["DISZ"][0]
        self.lens_data["thickness_l4"] = self.surflist[5]["DISZ"][0]

        self.lens_data["material_l1"] = self.surflist[0]["GLAS"][0]
        self.lens_data["material_l2"] = self.surflist[1]["GLAS"][0]
        self.lens_data["material_l3"] = self.surflist[4]["GLAS"][0]
        self.lens_data["material_l4"] = self.surflist[5]["GLAS"][0]

        r0 = self.surflist[0]["DIAM"][0]
        r1 = self.surflist[1]["DIAM"][0]
        r2 = self.surflist[2]["DIAM"][0]
        r3 = self.surflist[4]["DIAM"][0]
        r4 = self.surflist[5]["DIAM"][0]
        r5 = self.surflist[1]["DIAM"][0]

        if not (r0 == r1 and r1 == r2 and r2 == r3 and r3 == r4 and r4 == r5):
            return FailedImport.unequal_surface_radius

        self.lens_data["radius"] = r0
        self.lens_data["type"] = "DoubletPair"

        # return lens_data
        # There isn't yet a component class which can load these
        return FailedImport.doublet_pairs_unsupported


class ZmfImporter:
    """Imports a .ZMF file. List of definition dicts available in
    .descriptors property after .import_parts has been called.
    """

    def __init__(self, zmf_filename):
        self.zmf_path = Path(zmf_filename)
        self.failed_imports = []

        self.manual_exclusions = [
            "Triplet",
            "GRIN Lens",
            "46227",
            "Fiber Collimation Pkg",
            "Collimator",
        ]

    def import_parts(self, match_name=None):
        """Import parts into local instance dict zmx_data.
        If match_name is None, will import all. Otherwise, only make available
        exact match to name.
        """
        self.read_zmf(match_name=match_name)
        self.make_pyot_descriptors()

    def read_zmf(self, match_name=None):
        """
        Reads a Zemax library (file ending zmf), and populates the
        self.zmx_data dict with the descriptions of each component.
        The zmx data, which is an ASCII description
        The key is the reference of each component.

        Code (GPL) taken and modified from:
        https://github.com/jordens/rayopt/blob/master/rayopt/zemax.py
        """
        with self.zmf_path.open("rb") as fp:
            self.zmx_data = {}
            head = Struct("<I")
            lens = Struct("<100sIIIIIIIdd")
            shapes = "?EBPM"
            (version,) = head.unpack(fp.read(head.size))
            assert version in (1001,)
            while True:
                li = fp.read(lens.size)
                if len(li) != lens.size:
                    if len(li) > 0:
                        print(fp, "additional data", repr(li))
                    break
                li = list(lens.unpack(li))
                li[0] = li[0].decode("latin1").strip("\0")
                li[3] = shapes[li[3]]
                description = fp.read(li[7])
                assert len(description) == li[7]
                description = zmf_decode(description, li[8], li[9])
                description = description.decode("latin1")
                assert description.startswith("VERS {:06d}\n".format(li[1]))

                if match_name is None:
                    self.zmx_data[li[0]] = description
                elif match_name == li[0]:
                    self.zmx_data[li[0]] = description
                    return

    def make_pyot_descriptors(self):
        self.descriptors = {}
        for k, v in self.zmx_data.items():
            try:
                descriptor = self.pyot_descriptor(k)
            except Exception as ex:
                print(f"Exception converting element {k}\n")
                print("Raw data : \n")
                print(textwrap.indent(v, " " * 4))
                print("Exception : \n")
                traceback.print_exc()
                sys.exit(0)
            if not isinstance(descriptor, FailedImport):
                self.descriptors[k] = descriptor
            else:
                self.failed_imports.append((k, str(descriptor)))

    def save_json(self):
        json_file = Path(self.zmf_path.parent / (self.zmf_path.stem + ".json"))
        with json_file.open("w") as jfp:
            json.dump(self.descriptors, jfp, indent=4)

    def save_failed_csv(self):
        failed_file = Path(self.zmf_path.parent / (self.zmf_path.stem + "_failed.csv"))
        with failed_file.open("w") as csv_fp:
            failed_writer = csv.writer(csv_fp)
            for row in self.failed_imports:
                failed_writer.writerow(row)

    def pyot_descriptor(self, key):
        """
        Returns a dictionary descriptor of a supported optical element, given
        key. If unsupported or unable to import, returns a FailedImport
        """

        data = self.zmx_data[key]
        surfaces = data.splitlines()

        # Separate header from data
        header = []

        # Interpret the header
        while True:
            if surfaces[0].startswith("SURF"):
                break
            line = surfaces.pop(0)
            header.append(line)

        lens_data = {}

        description = find_key("NOTE 0", header)
        name = find_key("NAME", header)

        # Thor uses the NAME key for description, while edmund puts there the
        # lens reference, and the description in NOTE 0
        if name != key:
            description = name

        mode = find_key("MODE", header)
        if mode != "SEQ":
            return FailedImport.mode_not_seq

        unit = find_key("UNIT", header).split()[0]

        if unit != "MM":
            return FailedImport.units_not_mm

        # Separate the surfaces into a list of dictionaries
        surflist = []
        for line in surfaces:
            # Ignore empty lines
            if line.split() == []:
                continue
            if line.startswith("SURF"):
                surflist.append(dict())
                continue
            line = line.lstrip()
            code = line[:4]
            data = line[5:].split()
            data = [convert(d) for d in data]

            # PARM items need to be in a sub-dict
            if code == "PARM":
                if "PARM" not in surflist[-1]:
                    surflist[-1]["PARM"] = {}
                surflist[-1]["PARM"][data[0]] = data[1]
            # otherwise, put data in array
            else:
                surflist[-1][code] = data

        # Flag
        if any(ex in description for ex in self.manual_exclusions):
            # print(f"\n *** \n Manually excluded : {key}")
            # print(libdata[key])
            # print('\n')
            # for s in surflist:
            #    print(s)
            # print('\n *** \n')
            return FailedImport.manual_exclusion

        # Delete the object plane and the image plane
        surflist.pop(0)
        surflist.pop()

        # Many lenses have a first surface with no glass, this mean air,
        # so they will be removed.
        # It seems this is user to mark the mounting tube
        if "GLAS" not in surflist[0]:
            surflist.pop(0)

        # Many lenses have a last surface with no glass, this mean air, so
        # they will be removed.
        # It seems this is user to mark the mounting tube location.
        if "GLAS" not in surflist[-1] and "GLAS" not in surflist[-2]:
            surflist.pop()

        # Remove mirrors
        if checkglas(surflist, "MIRROR"):
            return FailedImport.mirrors_unsupported

        # Don't include diffractive optics
        if checktype(surflist, "BINARY_2"):
            return FailedImport.diffractive_optics_unsupported

        # Don't include Fresnel lenses.
        if checktype(surflist, "FRESNELS"):
            return FailedImport.fresnet_unsupported

        if checktype(surflist, "CFRESNEL"):
            return FailedImport.fresnet_unsupported

        # Check that all the glass types are supported
        gcat = find_key("GCAT", header)
        for s in surflist:
            if "GLAS" in s:
                name = s["GLAS"][0]
                try:
                    material.get_from(name, gcat)
                except KeyError:
                    print("Unknown material", name)
                    return FailedImport.unknown_material_type

        # Select the importer to use and run the import
        args = {
            "description": description,
            "gcat": gcat,
            "key": key,
            "surflist": surflist,
        }
        optic_importers = (
            CylindricalImporter,
            AsphericImporter,
            SingletImporter,
            DoubletImporter,
            AirSpacedDoubletImporter,
            DoubletPairImporter,
            DoubletPairCentStopImporter,
        )
        for importer in optic_importers:
            i = importer(**args)
            if i.valid():
                return i.definition()
        return FailedImport.unknown


def main():

    parser = argparse.ArgumentParser(
        prog="Pyoptools ZMF importer",
        description=(
            "Finds optical component descriptions "
            "in ZMF files and produces a .json "
            "component descriptor file. Failed imports "
            "are logged to a .csv file."
        ),
        add_help=True,
    )

    parser.add_argument("ZMF_filename")
    parser.add_argument(
        "-p",
        "--part",
        type=str,
        help=("optional : for decoding just one part. " "Omit to decode all"),
    )
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        action="store_true",
        help=(
            "Output .json and .csv files. "
            "If not specified, output will be to terminal."
        ),
    )
    parser.add_argument(
        "-d",
        "--decode",
        required=False,
        action="store_true",
        help="Just print the raw decoded ZMX output",
    )

    args = parser.parse_args()

    imp = ZmfImporter(args.ZMF_filename)

    if not args.decode:
        imp.import_parts(match_name=args.part)
        print(
            f"Imported {len(imp.descriptors)} items. "
            f"Failed to import {len(imp.failed_imports)} items."
        )

        if args.output:
            imp.save_json()
            imp.save_failed_csv()
        else:
            print(json.dumps(imp.descriptors, indent=4))
            if len(imp.failed_imports) > 0:
                print("\n\nFailed imports:\n")
                for failed in imp.failed_imports:
                    print(f"{failed[0]} : {failed[1]}")

    else:
        imp.read_zmf(match_name=args.part)
        for k, v in imp.zmx_data.items():
            header = f"Part : {k}"
            print(header)
            print("─" * len(header) + "\n")
            print(v)
            print("\n")


if __name__ == "__main__":
    main()
