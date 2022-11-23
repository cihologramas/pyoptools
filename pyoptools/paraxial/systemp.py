from collections import namedtuple

#  from collections.abc import MutableSequence
from abc import abstractmethod
from math import isinf
from numpy import array, dot

PRay = namedtuple("PRay", ["y", "u"])


class PSurface:
    """
    Superclass to define paraxial surfaces
    """

    def __init__(self, d, n=1):
        self.d = d
        self.n = n

    def p_mat(self):
        """Method that returns the paraxial matrix representation used to
        calculate the (y,u) ray propagation a distance of d
        """
        return array(((1.0, self.d), (0, 1)))

    def matrix(self, na):
        """Return the full diffraction + propagation matrix related to the
        surface when the previous refractive index is na
        """
        return dot(self.p_mat(), self.r_mat(na))

    @abstractmethod
    def r_mat(self):
        """Method that returns the paraxial matrix representation used to
        calculate the refraction of the (y,u) ray on the surface.
        It must be overloaded by all the PSurface childrens.
        """
        pass

    def propagate(self, ray, na):
        """Calculates the full ray (u,v) diffraction+propagation.

        Returns:

        PRay instance
        """

        nray = dot(self.matrix(na), array(ray))
        return PRay(nray[0], nray[1])


class PPlane(PSurface):
    """
    Class to define a paraxial plane surface
    """

    def r_mat(self, na):
        """Matriz de la superficie correspondiente a la refracción"""
        return array(((1, 0), (0, na / self.n)))


class PSpherical(PSurface):
    def __init__(self, curv, d, n):
        super().__init__(d, n)
        self.curv = curv

    def r_mat(self, na):
        if self.curv != 0:
            r = 1 / self.curv
            m = array(((1, 0), ((na - self.n) / (r * self.n), na / self.n)))
            return m
        else:
            return array(((1, 0), (0, na / self.n)))


class PObj(PPlane):
    def __init__(self, h, d, n):
        super().__init__(d, n)
        self.h = h

    def r_mat(self, n):
        """Los "objetos" no refractan"""
        return array(((1, 0), (0, 1)))


class PImg(PObj):
    def __init__(self, h=None):
        super().__init__(h, 0, 1)


class PApe(PPlane):
    def __init__(self, r_app, d, n):
        super().__init__(d, n)
        self.r_app = r_app

    def r_mat(self, na):
        """Las "pupilas" no refractan"""
        return array(((1, 0), (0, 1)))

    def propagate(self, ray, na):
        """Calculates the full ray (u,v) diffraction+propagation.

        Returns:

        PRay instance or None
        """
        if ray.y <= self.r_app:
            nray = dot(self.matrix(na), array(ray))
            return PRay(nray[0], nray[1])
        else:
            return None


surfdict = {"SPH": PSpherical, "OBJ": PObj, "PLA": PPlane, "APE": PApe, "IMG": PImg}


class PSystem:  # (MutableSequence):
    def __init__(self, lst):

        self.surfaces = []

        for s in lst:
            tp, data = s
            if isinstance(tp, str):
                tp = surfdict[tp](**data)
                self.surfaces.append(tp)
            else:
                raise ValueError

    def propagate(self, ray):
        """Calculates the full ray (u,v) diffraction+propagation.

        Returns:

        PRay instance or None
        """
        na = self.surfaces[0].n
        Z = 0
        ZL = [Z]
        YL = [ray.y]
        UL = [ray.u]

        for s in self.surfaces:

            ray = s.propagate(ray, na)
            na = s.n
            Z = Z + s.d
            if not isinstance(s, PImg):
                ZL.append(Z)
                YL.append(ray.y)
                UL.append(ray.u)
        return ZL, YL, UL

    def get_matrix(self, ua=False):
        """Get total system matrix
        if ua == True, stop at the first system aperture found
        """
        na = self.surfaces[0].n
        mat = array(((1, 0), (0, 1)))
        for s in self.surfaces:
            if isinstance(s, PApe):
                break
            mat = dot(s.matrix(na), mat)
            na = s.n
        if ua:
            return mat, s
        else:
            return mat

    def get_principal_ray(self):

        h = self.surfaces[0].h
        MS, APP = self.get_matrix(ua=True)
        u = -h * MS[0, 0] / MS[0, 1]
        return PRay(h, u)

    def get_marginal_ray(self):
        """Obtener el rayo marginal del sistema S."""

        MS, APP = self.get_matrix(ua=True)
        h = APP.r_app
        u = h / MS[0, 1]
        return PRay(0, u)

    def reverse(self):

        rl = []

        for i in range(len(self.surfaces) - 1, 0, -1):
            s0 = self.surfaces[i]
            s1 = self.surfaces[i - 1]
            d = s1.d
            n = s1.n

            if isinstance(s0, PImg):
                if s0.h is None:
                    hobj = self.surfaces[0].h
                    if not isinf(hobj) and not (hobj is None):
                        pr = self.get_principal_ray()
                        Z, Y, U = self.propagate(pr)
                        h = Y[-1]
                    else:
                        h = hobj
                else:
                    h = s0.h
                s = ["OBJ", {"d": d, "n": n, "h": h}]
            elif isinstance(s0, PApe):
                s = ["APE", {"d": d, "n": n, "r_app": s0.r_app}]
            elif isinstance(s0, PSpherical):
                s = ["SPH", {"d": d, "n": n, "curv": -s0.curv}]

            rl.append(s)

        rl.append(["IMG", {}])

        return PSystem(rl)


if __name__ == "__main__":
    import pylab as pl
    from math import inf

    SP = [
        ["OBJ", {"d": inf, "n": 1.0, "h": None}],
        ["SPH", {"curv": 1 / 37.4, "d": 5.9, "n": 1.613586}],
        ["SPH", {"curv": -1 / 341.48, "d": 12.93, "n": 1.0}],
        ["SPH", {"curv": -1 / 42.65, "d": 2.50, "n": 1.648338}],
        ["SPH", {"curv": 1 / 36.40, "d": 2.0, "n": 1}],
        ["APE", {"r_app": 10.3, "d": 9.85, "n": 1}],
        ["SPH", {"curv": 1 / 204.52, "d": 5.90, "n": 1.613586}],
        ["SPH", {"curv": -1 / 37.05, "d": 77.405, "n": 1}],
        ["IMG", {"h": 21.248022}],
    ]

    PS1 = PSystem(SP)
    PS = PS1.reverse()
    PR = PS.get_principal_ray()
    MR = PS.get_marginal_ray()

    Z, Y, U = PS.propagate(MR)
    print(Z)
    print(Y)
    print(U)
    pl.plot(Z, Y, "-o")
    pl.show()
