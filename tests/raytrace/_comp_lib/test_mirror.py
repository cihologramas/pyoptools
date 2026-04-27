from pyoptools.raytrace.comp_lib import RoundMirror


def test_roundmirror_default_no_wedge():
    rm = RoundMirror()
    _, _, rotation = rm.surflist["S1"]
    assert rotation == (0, 0, 0)
    _, _, rotation = rm.surflist["S2"]
    assert rotation == (0, 0, 0)


def test_roundmirror_wedge_angle_s2_rotated():
    angle = 0.02
    rm = RoundMirror(wedge_angle=angle)
    _, _, rotation = rm.surflist["S1"]
    assert rotation == (0, 0, 0)
    _, pos, rotation = rm.surflist["S2"]
    assert pos == (0, 0, 10)
    assert rotation == (0, angle, 0)


def test_roundmirror_wedge_zero_backward_compat():
    rm1 = RoundMirror(wedge_angle=0.0)
    rm2 = RoundMirror()
    _, _, rot1 = rm1.surflist["S2"]
    _, _, rot2 = rm2.surflist["S2"]
    assert rot1 == rot2
