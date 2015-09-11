from nose.tools import assert_equal


def test_circular():
    from pyoptools.raytrace.shape import Circular
    c = Circular()
    assert_equal(c.radius, 1.0, 'Unexpected radius')

if __name__ == '__main__':
    import nose
    nose.runmodule()
