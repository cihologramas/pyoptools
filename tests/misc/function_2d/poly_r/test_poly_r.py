import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal
from pyoptools.misc.function_2d.poly_r.poly_r import PolyR

def test_poly_r_init():
    """Test PolyR initialization"""
    coeff = [1.0, 2.0, 3.0]  # z = 1 + 2r + 3r^2
    poly = PolyR(coeff)
    assert poly.order == 3
    assert poly.dx is None
    assert poly.dy is None

def test_poly_r_eval_simple():
    """Test PolyR evaluation at origin and simple points"""
    # z = 1 + 2r + 3r^2
    poly = PolyR([1.0, 2.0, 3.0])
    
    # Test points as 2D arrays
    x = np.array([[0.0]], dtype=np.float64)
    y = np.array([[0.0]], dtype=np.float64)
    
    # At origin (0,0), r=0 so z=1
    result = poly.eval(x, y)
    assert_array_almost_equal(result, np.array([[1.0]]))
    
    # At point (1,0), r=1 so z=1+2+3=6
    x = np.array([[1.0]], dtype=np.float64)
    y = np.array([[0.0]], dtype=np.float64)
    result = poly.eval(x, y)
    assert_array_almost_equal(result, np.array([[6.0]]))

def test_poly_r_eval_grid():
    """Test PolyR evaluation on a grid of points"""
    # z = 1 + r  (simpler polynomial for easier testing)
    poly = PolyR([1.0, 1.0])
    
    # Create a 2x2 grid
    x = np.array([[0.0, 1.0],
                  [0.0, 1.0]], dtype=np.float64)
    y = np.array([[0.0, 0.0],
                  [1.0, 1.0]], dtype=np.float64)
    
    
    
    expected = np.array([[1.0, 2.0],           # r = 0, 1
                        [2.0, 1+np.sqrt(2)]])  # r = 1, sqrt(2)
    
    result = poly.eval(x, y)
    assert_array_almost_equal(result, expected)

def test_poly_r_dxdy():
    """Test partial derivatives of PolyR"""
    # z = 1 + 2r + 3r^2
    poly = PolyR([1.0, 2.0, 3.0])
    dx, dy = poly.dxdy()
    
    # Test derivatives at a point
    x = np.array([[1.0]], dtype=np.float64)
    y = np.array([[0.0]], dtype=np.float64)
    
    # For r = 1 (at point (1,0)):
    # dx = (x/r)(2 + 6r) = (1)(2 + 6(1)) = 8
    # dy = (y/r)(2 + 6r) = (0)(2 + 6(1)) = 0
    dx_result = dx.eval(x, y)
    dy_result = dy.eval(x, y)
    assert_array_almost_equal(dx_result, np.array([[8.0]]))
    assert_array_almost_equal(dy_result, np.array([[0.0]]))


    # Test derivatives at a point
    x = np.array([[1.0]], dtype=np.float64)
    y = np.array([[1.0]], dtype=np.float64)
    
    # For r = 1 (at point (1,0)):
    # dx = (x/r)(2 + 6r) = (1)(2 + 6(1)) = 8
    # dy = (y/r)(2 + 6r) = (0)(2 + 6(1)) = 0
    dx_result = dx.eval(x, y)
    dy_result = dy.eval(x, y)
    assert_array_almost_equal(dx_result, np.array([[np.sqrt(2)+6]]))
    assert_array_almost_equal(dy_result, np.array([[np.sqrt(2)+6]]))




def test_poly_r_eval_radial_symmetry():
    """Test that the polynomial is radially symmetric"""
    poly = PolyR([1.0, 2.0, 3.0])
    
    # Points with same radius but different angles
    r = 2.0
    angles = [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]
    
    results = []
    for angle in angles:
        x = np.array([[r * np.cos(angle)]], dtype=np.float64)
        y = np.array([[r * np.sin(angle)]], dtype=np.float64)
        results.append(poly.eval(x, y)[0,0])
    
    # All points should give same result since they're at same radius
    assert_array_almost_equal(results, [results[0]] * len(results))

# def test_poly_r_invalid_input():
#    """Test handling of invalid input"""
#    with pytest.raises(ValueError):
#        # Test with empty coefficient list
#        PolyR([])
    
#    poly = PolyR([1.0, 2.0])
    
#    with pytest.raises(ValueError):
#        # Test with mismatched x,y shapes
#        x = np.array([[1.0]], dtype=np.float64)
#        y = np.array([[1.0, 2.0]], dtype=np.float64)
#        poly.eval(x, y)

def test_poly_r_high_order():
    """Test high-order polynomial evaluation"""
    # z = 1 + r + r^2 + r^3 + r^4 + r^5
    coeff = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    poly = PolyR(coeff)
    
    x = np.array([[1.0]], dtype=np.float64)
    y = np.array([[0.0]], dtype=np.float64)
    
    # At (1,0), r=1 so sum should be 6
    result = poly.eval(x, y)
    assert_array_almost_equal(result, np.array([[6.0]]))

def test_poly_r_zero_polynomial():
    """Test zero polynomial (all coefficients zero)"""
    poly = PolyR([0.0, 0.0, 0.0])
    
    x = np.array([[1.0, 2.0]], dtype=np.float64)
    y = np.array([[1.0, 2.0]], dtype=np.float64)
    
    result = poly.eval(x, y)
    assert_array_almost_equal(result, np.zeros((1, 2)))