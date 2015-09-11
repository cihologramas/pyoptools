import pep8
import os
import pyoptools
from nose.tools import assert_equal

PEP8_ADDITIONAL_IGNORE = []
EXCLUDE_FILES = ['zernike.py',
                 'gs.py',
                 'calc.py',
                 'ray_source.py',
                 '__init__.py',
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
                 'ipynbplotutils.py',
                 'glwindow.py',
                 'glplotframe2.py',
                 'glplotframe.py',
                 'all.py'
                 ]


def test_pep8_conformance():

    dirname = os.path.dirname(pyoptools.__file__)
    pep8style = pep8.StyleGuide()

    # Extend the number of PEP8 guidelines which are not checked.
    pep8style.options.ignore = (pep8style.options.ignore +
                                tuple(PEP8_ADDITIONAL_IGNORE))
    pep8style.options.exclude.extend(EXCLUDE_FILES)

    result = pep8style.check_files([dirname])
    msg = "Found code syntax errors (and warnings)."
    assert_equal(result.total_errors, 0, msg)

if __name__ == '__main__':
    import nose
    nose.runmodule()
