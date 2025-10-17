"""Plotting utilities that depend only on Matplotlib. To be used in
`jupyter notebooks <http://jupyter.org>`_

"""
from pylab import plot, axis, array
from pyoptools.misc.pmisc import wavelength2RGB


def spot_diagram(s):
    """Plot the spot diagram for the given surface, or element.

    Args:
        s: Object (usually :class:`~pyoptools.raytrace.comp_lib.CCD`)
            whose spot diagram will be plotted.
    """
    hl = s.hit_list
    X = []
    Y = []
    COL = []
    if len(hl) > 0:
        for i in hl:
            p = i[0]
            # Hitlist[1] points to the incident ray
            col = wavelength2RGB(i[1].wavelength)
            X.append(p[0])
            Y.append(p[1])
            COL.append(col)
    max = array(X + Y).max
    min = array(X + Y).min
    plot(
        X,
        Y,
        "o",
    )
    axis("equal")


def spot_diagram_c(s):
    """Plot the spot diagram for the given surface, or element using
    the rays colors.

    Args:
        s: Object (usually :class:`~pyoptools.raytrace.comp_lib.CCD`)
            whose spot diagram will be plotted.

    """
    hl = s.hit_list
    X = []
    Y = []
    COL = []
    if len(hl) > 0:
        for i in hl:
            p = i[0]
            # Hitlist[1] points to the incident ray
            col = wavelength2RGB(i[1].wavelength)
            plot(p[0], p[1], "o", color=col)
            # X.append(p[0])
            # Y.append(p[1])
            # COL.append(col)
    # max=array(X+Y).max
    # min=array(X+Y).min
    axis("equal")
