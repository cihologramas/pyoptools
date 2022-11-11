from pyoptools.raytrace.library import library
from pyoptools.raytrace._comp_lib.optic_factory import optic_factory

import unittest

class TestAllOptics(unittest.TestCase):

    #@unittest.skip
    def test_all_optics(self):
        i = -1
        for part, descriptor in library.items():
            i += 1
            #print(part)
            with self.subTest(part, i=i):
                optic = optic_factory(**descriptor)

    def test_attribute_access(self):

        lens = library.Thorlabs['LB1862']
        lens = library.Edmund['08068']

    def test_descriptor(self):

        self.assertEqual(library.Edmund.descriptor('08068')['radius'], 12.2)
        with self.assertRaises(KeyError):
            library.Edmund.descriptor('None')

    #def test_user_lib(self):



if __name__ == '__main__':
    unittest.main()

# failes on:
"""
67395
67384
63228
63227
49566
49563
49561
48431
46230
46229
46228
46223
46122
43832
43822
43645
"""

