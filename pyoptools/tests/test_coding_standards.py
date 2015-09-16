import pep8
import os
import pyoptools
from nose.tools import assert_equal

PEP8_ADDITIONAL_IGNORE = []

# THIS IS THE LIST OF SHAME!!!!
# If you modify a single one of these files, take the opportunity,
# and clean it!.
# Once it's clean, remove it from this list
EXCLUDE_FILES = ['zernike.py',
                 'gs.py',
                 'calc.py',
                 'ray_source.py',
                 'material.py',
                 'genlenslib.py',
                 'library.py',
                 'poly_expansion.py',
                 'stop.py',
                 'spherical_lens.py',
                 'prism.py',
                 'ideallens.py',
                 'cube.py',
                 'compound_lens.py',
                 'ccd.py',
                 'frft.py',
                 'misc.py',
                 'shell_frame.py',
                 'wxrayos.py',
                 'comp_lib.py',
                 'plotutils.py',
                 'plot_frame.py',
                 'oglframe.py',
                 'logutils.py',
                 'glwindow.py',
                 'glplotframe2.py',
                 'glplotframe.py',
                 'all.py',
                 'field.pyx',
                 'fields.pyx',
                 'psurfrep.pyx',
                 'cpsurfrep.pyx',
                 'triangular.pyx',
                 'shape.pyx',
                 'rectangular.pyx',
                 'polygon.pyx',
                 'circular.pyx',
                 'ray.pyx',
                 'component.pyx',
                 'system.pyx',
                 'taylor_poly.pyx',
                 'surface.pyx',
                 'spherical.pyx',
                 'plane_mask.pyx',
                 'detector.pyx',
                 'cylindrical.pyx',
                 'plane.pyx',
                 'idealpplanes.pyx',
                 'cylinder.pyx',
                 'aspherical.pyx',
                 'idealsurface.pyx',
                 'aperture.pyx',
                 'resources.pyx',
                 'picklable.pyx',
                 'lsq.pyx',
                 'cmisc.pyx',
                 'definitions.pyx',
                 'plist.pyx',
                 'Poly2D.pyx',
                 ]


def test_pep8_conformance():

    dirname = os.path.dirname(pyoptools.__file__)
    pep8style = pep8.StyleGuide()

    # Extend the number of PEP8 guidelines which are not checked.
    pep8style.options.ignore = (pep8style.options.ignore +
                                tuple(PEP8_ADDITIONAL_IGNORE))
    pep8style.options.exclude.extend(EXCLUDE_FILES)
    pep8style.options.filename = ['*.py', '*.pyx']

    result = pep8style.check_files([dirname])
    msg = "Found code syntax errors (and warnings)."
    assert_equal(result.total_errors, 0, msg)

if __name__ == '__main__':
    import nose
    nose.runmodule()
