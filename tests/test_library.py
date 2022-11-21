from pyoptools.raytrace.library import library
from pyoptools.raytrace._comp_lib.optic_factory import optic_factory
from pyoptools.raytrace._comp_lib.spherical_lens import SphericalLens

import unittest
import tempfile

class TestLibrary(unittest.TestCase):

    #@unittest.skip
    def test_all_optics(self):
        i = -1
        for part, descriptor in library.items():
            i += 1
            with self.subTest(part, i=i):
                optic = optic_factory(**descriptor)

    def test_direct_access(self):

        lens = library['LB1862']
        self.assertIsInstance(lens, SphericalLens)

    def test_attribute_access(self):

        lens = library.Thorlabs['LB1862']
        self.assertIsInstance(lens, SphericalLens)

        lens = library.Edmund['08068']
        self.assertIsInstance(lens, SphericalLens)

    def test_descriptor(self):

        self.assertEqual(library.Edmund.descriptor('08068')['radius'], 12.2)

        self.assertEqual(library.descriptor('08068')['radius'], 12.2)

        with self.assertRaises(KeyError):
            library.Edmund.descriptor('None')

        with self.assertRaises(KeyError):
            library.descriptor('None')

    def test_parts(self):

        self.assertTrue('LB1862' in library.Thorlabs.parts())
        self.assertTrue('LB1862' in library.parts())

    def test_user_lib(self):

        ntf = tempfile.NamedTemporaryFile('w', suffix='.json')
        with ntf as fp:
            fp.write("""
               {
                    "a_test_lens": {
                        "description": "test",
                        "material": "BK7",
                        "glass_catalogs": "SCHOTT",
                        "thickness": 3.5,
                        "radius": 12.7,
                        "curvature_s2": 0.0,
                        "curvature_s1": 0.02,
                        "type": "SphericalLens"
                    }
                }
            """)
            fp.flush()
            fp.seek(0)

            library.add(ntf.name)

            self.assertEqual(library.descriptor('a_test_lens')['thickness'],
                             3.5)

            self.assertEqual(library.descriptor('a_test_lens')['description'],
                             'test')

if __name__ == '__main__':
    unittest.main()
