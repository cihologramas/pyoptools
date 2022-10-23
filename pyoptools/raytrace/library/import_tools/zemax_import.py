import numpy as np
from struct import Struct
from pathlib import Path
from enum import Enum, auto
import json
import sys
import csv
import textwrap
import traceback

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
            rv = rv+s[len(key)+1:]
    except IndexError:
        # key not found
        rv = ""
    return rv

def checktype(surflist, t):
    for s in surflist:
        if (t in s) or (s.get('TYPE', [None])[0]) == t:
            return True
    return False


def checkglas(surflist, t):
    for s in surflist:
        if (s.get('GLAS', [None])[0]) == t:
            return True
    return False

def zmf_decode(data, a, b):
    iv = np.cos(6*a + 3*b)
    iv = np.cos(655*(np.pi/180)*iv) + iv
    p = np.arange(len(data))
    k = 13.2*(iv + np.sin(17*(p + 3)))*(p + 1)
    k = (int(("{:.8e}".format(_))[4:7]) for _ in k)
    #data = np.fromstring(data, np.uint8)
    data = np.frombuffer(data, np.uint8).copy()
    data ^= np.fromiter(k, np.uint8, len(data))
    return data.tobytes()

#def zmx_read(fn):
#    f = open(fn, "rU")
#    data = f.read()
#    return zmx_parse(data)

class FailedImport(Enum):
    manual_exclusion = auto()
    mode_not_seq = auto()
    units_not_mm = auto()
    mirrors_unsupported = auto()
    aspheres_unsupported = auto()
    diffractive_optics_unsupported = auto()
    conical_unsupported = auto()
    fresnet_unsupported = auto()
    acylindrical_unsupported = auto()
    no_primary_surface = auto()
    diameter_undefined = auto()
    unequal_surface_radius = auto()
    unknown = auto()

class ZmfImporter:

    def __init__(self, zmf_filename):
        self.zmf_path = Path(zmf_filename)
        self.failed_imports = []

        self.manual_exclusions = ['Triplet']

    def import_all(self):
        self.read_zmf()
        self.make_pyot_descriptors()

    def read_zmf(self):
        """
        Reads a Zemax library (file ending zmf), and populates the
        self.zmx_data dict with the descriptions of each component.
        The zmx data, which is an ASCII description
        The key is the reference of each component.

        Code (GPL) taken and modified from:
        https://github.com/jordens/rayopt/blob/master/rayopt/zemax.py
        """
        with self.zmf_path.open('rb') as fp:
            self.zmx_data = {}
            head = Struct("<I")
            lens = Struct("<100sIIIIIIIdd")
            shapes = "?EBPM"
            version, = head.unpack(fp.read(head.size))
            assert version in (1001, )
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
                self.zmx_data[li[0]] = description

    def make_pyot_descriptors(self):
        self.descriptors = {}
        for k, v in self.zmx_data.items():
            try:
                descriptor = self.pyot_descriptor(k)
            except Exception as ex:
                print(f"Exception converting element {k}\n")
                print("Raw data : \n")
                print(textwrap.indent(v, ' '*4))
                print("Exception : \n")
                traceback.print_exc()
                sys.exit(0)
            if not isinstance(descriptor, FailedImport):
                self.descriptors[k] = descriptor
            else:
                self.failed_imports.append((k, str(descriptor)))

    def save_json(self):
        json_file = Path(self.zmf_path.parent/(self.zmf_path.stem + '.json'))
        with json_file.open('w') as jfp:
            json.dump(self.descriptors, jfp, indent = 4)

    def save_failed_csv(self):
        failed_file = Path(self.zmf_path.parent/(self.zmf_path.stem + '_failed.csv'))
        with failed_file.open('w') as csv_fp:
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
        if name == key:
            lens_data["description"] = description
        else:
            lens_data["description"] = name

        # lens_data["description"]=name+"\n"+description

        # I still don't know what to do with the mode
        mode = find_key("MODE", header)

        if mode != "SEQ":
            return FailedImport.mode_not_seq

        unit = find_key("UNIT", header).split()[0]

        if unit != "MM":
            return FailedImport.units_not_mm

        #gcat = find_key("GCAT", header)

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
            surflist[-1][code] = data

        # Flag
        if any(ex in lens_data['description'] for ex in self.manual_exclusions):
            #print(f"\n *** \n Manually excluded : {key}")
            #print(libdata[key])
            #print('\n')
            #for s in surflist:
            #    print(s)
            #print('\n *** \n')
            return FailedImport.manual_exclusion

        # Delete the object plane and the image plane
        surflist.pop(0)
        surflist.pop()

        # Many lenses have a first surface with no glass, this mean air, so they will be
        # removed. It seems this is user to mark the mounting tube
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

        # Don't include aspheres
        if checktype(surflist, 'EVENASPH'):
            return FailedImport.aspheres_unsupported

        # Don't include diffractive optics
        if checktype(surflist, 'BINARY_2'):
            return FailedImport.diffractive_optics_unsupported

        if checktype(surflist, 'CONI'):
            return FailedImport.conical_unsupported

        # Don't include Fresnel lenses.
        if checktype(surflist, 'FRESNELS'):
            return FailedImport.fresnet_unsupported

        if checktype(surflist, 'CFRESNEL'):
            return FailedImport.fresnet_unsupported

        # Cylindrical lenses. It seems these are modeled as 'toroidal'
        # Currently working for plano-convex and plano-concave, rectangular
        if checktype(surflist, 'TOROIDAL'):
            lens_data["type"] = "CylindricalLens"

            #print(f"Processing CylindricalLens {key}")

            # Get the size from the first surface
            all_sqaps = [tuple(s['SQAP']) for s in surflist if 'SQAP' in s]
            # Verify all SQAP elements same for all surfaces
            if len(set(all_sqaps)) != 1:
                return FailedImport.acylindrical_unsupported

            lens_data["size"] = tuple(reversed(all_sqaps[0][0:2]))

            # Get thickness. Seems the DISZ element in the surface with
            # COMM equal to the key is where this is found
            surf = [s for s in surflist if 'COMM' in s and s['COMM'][0] == key]
            if len(surf) != 1:
                return FailedImport.no_primary_surface

            lens_data["thickness"] = float(surf[0]['DISZ'][0])
            lens_data["curvature_s1"] = float(surf[0]['CURV'][0])
            lens_data["curvature_s2"] = 0.0
            lens_data["material"] = surf[0]['GLAS'][0]

            return lens_data

        # Identify the type of lenses from the number of surfaces
        ns = len(surflist)

        # Normal Singlets
        if ns == 2:
            c0 = surflist[0]["CURV"][0]
            c1 = surflist[1]["CURV"][0]
            d0 = surflist[0]["DISZ"][0]
            # Not all lenses have DIAM in both surfaces
            if "DIAM" in surflist[0] and "DIAM" in surflist[1]:
                r0 = surflist[0]["DIAM"][0]
                r1 = surflist[1]["DIAM"][0]
            elif "DIAM" in surflist[0]:
                r0 = surflist[0]["DIAM"][0]
                r1 = r0
            elif "DIAM" in surflist[1]:
                r0 = surflist[1]["DIAM"][0]
                r1 = r0
            else:
                return FailedImport.diameter_undefined

            g0 = surflist[0]["GLAS"][0]

            # Check that the surfaces are equal, if not issue an error
            if not r0 == r1:
                return FailedImport.unequal_surface_radius

            lens_data["material"] = g0
            #lens_data["glass_catalogs"] = gcat  # This is not used for the moment
            lens_data["thickness"] = surflist[0]["DISZ"][0]
            lens_data["radius"] = r0
            lens_data["curvature_s2"] = surflist[1]["CURV"][0]
            lens_data["curvature_s1"] = surflist[0]["CURV"][0]
            lens_data["type"] = "SphericalLens"

            return lens_data

            # return CL.SphericalLens(r0,d0,c0,c1,material=m0)

        elif ns == 3:
            lens_data["curvature_s1"] = surflist[0]["CURV"][0]
            lens_data["curvature_s2"] = surflist[1]["CURV"][0]
            lens_data["curvature_s3"] = surflist[2]["CURV"][0]
            lens_data["thickness_l1"] = surflist[0]["DISZ"][0]
            lens_data["thickness_l2"] = surflist[1]["DISZ"][0]
            r0 = surflist[0]["DIAM"][0]
            r1 = surflist[1]["DIAM"][0]
            r2 = surflist[2]["DIAM"][0]

            # Check that the surfaces are equal, if not issue an error
            if not (r0 == r1 and r1 == r2):
                return FailedImport.unequal_surface_radius

            lens_data["radius"] = r0

            lens_data["material_l1"] = surflist[0]["GLAS"][0]
            lens_data["material_l2"] = surflist[1]["GLAS"][0]
            #lens_data["glass_catalogs"] = gcat
            # return CL.Doublet(r0,c0,c1,c2, d0,d1,m0,m1)
            lens_data["type"] = "Doublet"
            return lens_data

        elif ns == 4:
            lens_data["curvature_s1"] = surflist[0]["CURV"][0]
            lens_data["curvature_s2"] = surflist[1]["CURV"][0]
            lens_data["curvature_s3"] = surflist[2]["CURV"][0]
            lens_data["curvature_s4"] = surflist[3]["CURV"][0]

            lens_data["thickness_l1"] = surflist[0]["DISZ"][0]
            lens_data["air_gap"] = surflist[1]["DISZ"][0]
            lens_data["thickness_l2"] = surflist[2]["DISZ"][0]

            # Verificar que las superficies son iguales, si no emitir un error
            # assert r0==r1 and r1== r2 and r2 ==r3 Esto no siempre se cumple

            r0 = surflist[0]["DIAM"][0]
            r1 = surflist[1]["DIAM"][0]
            r2 = surflist[2]["DIAM"][0]
            r3 = surflist[3]["DIAM"][0]

            # Verificar que las superficies son iguales, si no emitir un error
            # assert r0==r1 and r1== r2 and r2 ==r3 #Esto no siempre se cumple

            lens_data["radius"] = r0

            lens_data["material_l1"] = surflist[0]["GLAS"][0]
            # g1=surflist[1]["GLAS"][0] Este parametro no existe. Es aire
            lens_data["material_l2"] = surflist[2]["GLAS"][0]

            #lens_data["glass_catalogs"] = gcat

            lens_data["type"] = "AirSpacedDoublet"
            return lens_data

        # Doublet pair
        # Thorlabs have the term 'Matched Achromatic Pair'
        # in the description. Check for this to avoid
        # collision with triplets
        elif ns == 6 and 'Pair' in lens_data['description']:
            lens_data["curvature_s1"] = surflist[0]["CURV"][0]
            lens_data["curvature_s2"] = surflist[1]["CURV"][0]
            lens_data["curvature_s3"] = surflist[2]["CURV"][0]
            lens_data["curvature_s4"] = surflist[3]["CURV"][0]
            lens_data["curvature_s5"] = surflist[4]["CURV"][0]
            lens_data["curvature_s6"] = surflist[5]["CURV"][0]

            lens_data["thickness_l1"] = surflist[0]["DISZ"][0]
            lens_data["thickness_l2"] = surflist[1]["DISZ"][0]
            lens_data["air_gap"] = surflist[2]["DISZ"][0]
            lens_data["thickness_l3"] = surflist[3]["DISZ"][0]
            lens_data["thickness_l4"] = surflist[4]["DISZ"][0]

            lens_data["material_l1"] = surflist[0]["GLAS"][0]
            lens_data["material_l2"] = surflist[1]["GLAS"][0]
            lens_data["material_l3"] = surflist[3]["GLAS"][0]
            lens_data["material_l4"] = surflist[4]["GLAS"][0]

            r0 = surflist[0]["DIAM"][0]
            r1 = surflist[1]["DIAM"][0]
            r2 = surflist[2]["DIAM"][0]
            r3 = surflist[3]["DIAM"][0]
            r4 = surflist[4]["DIAM"][0]
            r5 = surflist[5]["DIAM"][0]

            if not (r0 == r1 and r1 == r2 and r2 == r3 and r3 == r4 and r4 == r5):
                return FailedImport.unequal_surface_radius

            lens_data["radius"] = r0
            #lens_data["glass_catalogs"] = gcat  # This is not used for the moment

            lens_data["type"] = "DoubletPair"

            #print('Imported a doublet pair with ns = 6\n')
            #print(self.zmx_data[key])

            return lens_data

        # Doublet pair with stop in the middle, stop ignored
        # Thorlabs have the term 'Matched Achromatic Pair'
        # in the description. Check for this to avoid
        # collision with triplets
        elif ns == 7 and 'Pair' in lens_data['description']:
            lens_data["curvature_s1"] = surflist[0]["CURV"][0]
            lens_data["curvature_s2"] = surflist[1]["CURV"][0]
            lens_data["curvature_s3"] = surflist[2]["CURV"][0]
            lens_data["curvature_s4"] = surflist[4]["CURV"][0]
            lens_data["curvature_s5"] = surflist[5]["CURV"][0]
            lens_data["curvature_s6"] = surflist[6]["CURV"][0]

            lens_data["thickness_l1"] = surflist[0]["DISZ"][0]
            lens_data["thickness_l2"] = surflist[1]["DISZ"][0]
            lens_data["air_gap"] = surflist[2]["DISZ"][0]+surflist[3]["DISZ"][0]
            lens_data["thickness_l3"] = surflist[4]["DISZ"][0]
            lens_data["thickness_l4"] = surflist[5]["DISZ"][0]

            lens_data["material_l1"] = surflist[0]["GLAS"][0]
            lens_data["material_l2"] = surflist[1]["GLAS"][0]
            lens_data["material_l3"] = surflist[4]["GLAS"][0]
            lens_data["material_l4"] = surflist[5]["GLAS"][0]

            r0 = surflist[0]["DIAM"][0]
            r1 = surflist[1]["DIAM"][0]
            r2 = surflist[2]["DIAM"][0]
            r3 = surflist[4]["DIAM"][0]
            r4 = surflist[5]["DIAM"][0]
            r5 = surflist[1]["DIAM"][0]

            if not(r0 == r1 and r1 == r2 and r2 == r3 and r3 == r4 and r4 == r5):
                return FailedImport.unequal_surface_radius

            lens_data["radius"] = r0
            #lens_data["glass_catalogs"] = gcat  # This is not used for the moment

            lens_data["type"] = "DoubletPair"

            #print('Imported a doublet pair with ns = 7\n')
            #print(self.zmx_data[key])

            return lens_data

        else:
            return FailedImport.unknown
            # Print elements present in the library that were not processed or explicitly excluded
            #print("***", key)
            #print(description)
            #print("Surfaces=", len(surflist))
            #for i in surflist:
            #    print("*", i)

def main(filename):
    imp = ZmfImporter(filename)
    imp.import_all()
    print(f"Imported {len(imp.descriptors)} items. "
          f"Failed to import {len(imp.failed_imports)} items.")


    imp.save_json()
    imp.save_failed_csv()


if __name__ == '__main__':
    main(sys.argv[1])
