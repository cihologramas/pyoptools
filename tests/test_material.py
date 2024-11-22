import unittest

from pyoptools.raytrace.mat_lib import material
from pyoptools.raytrace.mat_lib.mat_eq import Material as OpticalMaterial


class TestMaterial(unittest.TestCase):

    def test_access(self):

        n = 1.5214

        a = material["N-BK7"]
        b = material.schott["N-BK7"]
        c = material["BK7"]  # from alias

        self.assertAlmostEqual(a.n(0.5), n, places=3)
        self.assertAlmostEqual(b.n(0.5), n, places=3)
        self.assertAlmostEqual(c.n(0.5), n, places=3)

    def test_compounds(self):

        a = material.inorganic["Ag:"]
        b = material.inorganic["Ag"]
        c = material.inorganic["Ag:Yang"]

        self.assertIsInstance(a, OpticalMaterial)
        self.assertIsInstance(b, OpticalMaterial)
        self.assertIsInstance(c, OpticalMaterial)

        octane = material.organic["C8H18"]

        self.assertIsInstance(octane, OpticalMaterial)

    def test_multi_matches(self):

        with self.assertRaises(KeyError):
            material["SF11"]

    def test_find_material(self):

        matches = len(material.find_material("BK7", False))
        self.assertGreater(matches, 0)

    def test_find_material_compound(self):

        matches = len(material.find_material("styrene", False))
        self.assertGreater(matches, 0)

    def test_find_material_none(self):
        matches = len(material.find_material("not a material", False))
        self.assertEqual(matches, 0)

    def test_ailised_glass_wrong_catalog(self):
        with self.assertRaises(KeyError):
            material.ami["BK7"]

    def test_ailised_compounds_wrong_catalog(self):
        with self.assertRaises(KeyError):
            material.inorganic["PMMA"]

        with self.assertRaises(KeyError):
            material.organic["CAF2"]

    def test_ailised_organic_compounds(self):
        a = material["PMMA"]
        b = material.organic["PMMA"]
        c = material.organic["(C5H8O2)n:Zhang-Tomson"]
        self.assertEqual(a, b)
        self.assertEqual(a, c)

    def test_ailised_inorganic_compounds(self):
        a = material["GAAS"]
        b = material.inorganic["GAAS"]
        c = material.inorganic["GaAs:Skauli"]
        self.assertEqual(a, b)
        self.assertEqual(a, c)

    def test_ailas_glass(self):
        a = material.schott["N-BK7"]
        b = material["BK7"]
        self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main()
