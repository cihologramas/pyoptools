import pytest
import numpy as np
from pyoptools.misc.poly_2d.poly_2d import Poly2D, index_to_powers, indices_to_powers, pxpy2i, ord2i

def test_index_to_powers():
    assert index_to_powers(0) == (0, 0)
    assert index_to_powers(1) == (1, 0)
    assert index_to_powers(2) == (0, 1)
    assert index_to_powers(3) == (2, 0)
    assert index_to_powers(4) == (1, 1)
    assert index_to_powers(5) == (0, 2)

def test_indices_to_powers():
    indices = np.array([0, 1, 2, 3, 4, 5], dtype=np.int32)
    x_powers, y_powers = indices_to_powers(indices)
    assert np.all(x_powers == [0, 1, 0, 2, 1, 0])
    assert np.all(y_powers == [0, 0, 1, 0, 1, 2])

def test_pxpy2i():
    assert pxpy2i(0, 0) == 0
    assert pxpy2i(1, 0) == 1
    assert pxpy2i(0, 1) == 2
    assert pxpy2i(2, 0) == 3
    assert pxpy2i(1, 1) == 4
    assert pxpy2i(0, 2) == 5

def test_ord2i():
    assert ord2i(0) == 1
    assert ord2i(1) == 3
    assert ord2i(2) == 6
    assert ord2i(3) == 10
    assert ord2i(4) == 15

def test_poly2d_initialization():
    coeff = [1.0, 2.0, 3.0]
    poly = Poly2D(coeff)

    assert poly.num_coefficients == 3
    assert np.all(poly.coefficients == np.array(coeff))

def test_poly2d_addition():
    poly1 = Poly2D([1, 2, 3])
    poly2 = Poly2D([4, 5, 6])
    poly_sum = poly1 + poly2

    assert poly_sum.coefficients.tolist() == [5, 7, 9]

def test_poly2d_subtraction():
    poly1 = Poly2D([1, 2, 3])
    poly2 = Poly2D([4, 5, 6])
    poly_diff = poly1 - poly2

    assert poly_diff.coefficients.tolist() == [-3, -3, -3]

def test_poly2d_multiplication():
    poly = Poly2D([1, 2, 3])
    poly_scaled = poly * 2

    assert poly_scaled.coefficients.tolist() == [2, 4, 6]

def test_poly2d_negation():
    poly = Poly2D([1, 2, 3])
    poly_neg = -poly

    assert poly_neg.coefficients.tolist() == [-1, -2, -3]

def test_poly2d_eval():
    poly = Poly2D([1, 1, 1])  # Simplest polynomial: 1 + x + y
    x = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    y = np.array([[1.0, 1.0], [1.0, 1.0]], dtype=np.float64)
    z = poly.eval(x, y)

    expected = np.array([[3.0, 4.0], [5.0, 6.0]])
    np.testing.assert_allclose(z, expected, rtol=1e-5)

def test_poly2d_str():
    poly = Poly2D([1.0, -1.0, 1.0, 2.0, 3.0, 4.0])
    assert str(poly) == "1.0 - x + y + 2.0x^2 + 3.0xy + 4.0y^2"

def test_poly2d_repr():
    poly = Poly2D([1.0, -1.0, 1.0])
    assert repr(poly) == "Poly2D(coeff=[1.0, -1.0, 1.0], order=1)"

def test_poly2d_dxdy():
    poly = Poly2D([1.0, -1.0, 1.0, 2.0, 3.0, 4.0])
    poly_dx, poly_dy = poly.dxdy()

    assert poly_dx.coefficients.tolist() == [-1, 4, 3, 0, 0, 0]
    assert poly_dy.coefficients.tolist() == [1, 3, 8, 0, 0, 0]
