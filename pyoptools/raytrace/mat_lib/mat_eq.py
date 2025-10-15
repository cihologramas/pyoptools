
# ------------------------------------------------------------------------------
# Copyright (c) 2007-2021, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amézquita Orozco
# Description:     Material definition helper class
# Symbols Defined: Material
#
#
# ------------------------------------------------------------------------------

"""
File that define the Material class and subclasses that follow
https://refractiveindex.info/database/doc/Dispersion%20formulas.pdf for the
dispersion formulas

"""
import yaml
import numpy
import io
from math import sqrt


class ModelNotImplemented(Exception):
    """Indicates that the model to calculate the index of refraction is not
    implemented yet"""

    pass


class Material:
    """Base class to define an optical material. It receives the YML text as
    defined in https://refractiveindex.info/about
    """

    def __init__(self, coef, cl, nd=None, vd=None):
        self.__nd__ = nd
        self.__vd__ = vd

        # Check if the number of coefficients agrees with the dispersion
        # formulas definition

        if cl is None:
            self.__coef__ = coef
        else:

            self.__coef__ = numpy.zeros(cl)
            self.__coef__[: coef.size] = coef

    @property
    def nd(self):
        if self.__nd__ is not None:
            return self.__nd__
        else:
            return self.n(0.58929)

    @property
    def vd(self):
        if self.__vd__ is not None:
            return self.__vd__
        else:
            return (self.n(0.5893) - 1) / (self.n(0.4861) - self.n(0.6563))

    def __eq__(self, other):
        return (
            isinstance(other, self.__class__)
            and numpy.array_equal(self.__coef__, other.__coef__)
            and self.__nd__ == other.__nd__
            and self.__vd__ == other.__vd__
        )


class Sellmeier(Material):
    """Class that define a material that complies with the Sellmeier dispersion
    model as defined at
    https://refractiveindex.info/database/doc/Dispersion%20formulas.pdf
    """

    def __init__(self, coef, nd=None, vd=None):
        Material.__init__(self, coef, 17, nd, vd)

    def n(self, wavelength=0.58929):
        n2 = 1 + self.__coef__[0]
        for i in range(1, 17, 2):
            n2 = n2 + (self.__coef__[i] * wavelength**2 / (wavelength**2 - self.__coef__[i + 1] ** 2))
        return sqrt(n2)


class Sellmeier_2(Material):
    """Class that define a material that complies with the Sellmeier_2
    dispersion model as defined at
    https://refractiveindex.info/database/doc/Dispersion%20formulas.pdf
    """

    def __init__(self, coef, nd=None, vd=None):
        Material.__init__(self, coef, 17, nd, vd)

    def n(self, wavelength=0.58929):
        n2 = 1 + self.__coef__[0]
        for i in range(1, 17, 2):
            n2 = n2 + (self.__coef__[i] * wavelength**2 / (wavelength**2 - self.__coef__[i + 1]))
        return sqrt(n2)


class Polynomial(Material):
    """Class that define a material that complies with the Polynomial
    dispersion model as defined at
    https://refractiveindex.info/database/doc/Dispersion%20formulas.pdf
    """

    def __init__(self, coef, nd=None, vd=None):
        Material.__init__(self, coef, 17, nd, vd)

    def n(self, wavelength=0.58929):
        n2 = self.__coef__[0]
        for i in range(1, 17, 2):
            n2 = n2 + self.__coef__[i] * wavelength ** self.__coef__[i + 1]
        return sqrt(n2)


class RefractiveIndex_Info(Material):
    def __init__(self, coef, nd=None, vd=None):
        raise ModelNotImplemented


class Cauchy(Material):
    def __init__(self, coef, nd=None, vd=None):
        Material.__init__(self, coef, 11, nd, vd)

    def n(self, wavelength=0.58929):
        n_ = self.__coef__[0]
        for i in range(1, 11, 2):
            n_ = n_ + self.__coef__[i] * wavelength ** self.__coef__[i + 1]
        return n_


class Gases(Material):
    def __init__(self, coef, nd=None, vd=None):
        raise ModelNotImplemented


class Herzberger(Material):
    def __init__(self, coef, nd=None, vd=None):
        raise ModelNotImplemented


class Retro(Material):
    def __init__(self, coef, nd=None, vd=None):
        raise ModelNotImplemented


class Exotic(Material):
    def __init__(self, coef, nd=None, vd=None):
        raise ModelNotImplemented


__models__ = [
    Sellmeier,
    Sellmeier_2,
    Polynomial,
    RefractiveIndex_Info,
    Cauchy,
    Gases,
    Herzberger,
    Retro,
    Exotic,
]


class Tabulated_N(Material):
    """Class that define a material that complies with the "tabulated n"
    format as provided by https://refractiveindex.info
    """

    def __init__(self, coef, nd=None, vd=None):
        Material.__init__(self, coef, None, nd, vd)

    def n(self, wavelength=0.58929):
        n_ = numpy.interp(wavelength, self.__coef__[:, 0], self.__coef__[:, 1])
        return n_


def from_yml(file_path):
    """Create a material instance from a YML file path as defined at
    https://refractiveindex.info/about
    """
    # with open(filename, encoding='utf-8') as f:

    # print('opening ', file_path)

    with file_path.open(encoding="utf-8") as f:
        mat = yaml.load(f, Loader=yaml.FullLoader)

    for c in mat["DATA"]:
        if "formula" in c["type"]:
            fn = int(c["type"].split()[1]) - 1
            coef = numpy.fromiter(c["coefficients"].split(), dtype=float)
            # Sometimes SPECS does not exists
            nd = mat.get("SPECS", {}).get("nd", None)
            vd = mat.get("SPECS", {}).get("Vd", None)
            return __models__[fn](coef, nd, vd)

        elif c["type"].startswith("tabulated n"):
            with io.StringIO(c["data"]) as data:
                coef = numpy.loadtxt(data, usecols=(0, 1))
            return Tabulated_N(coef)

        # Changed below to support tabulated nk etc formatted data
        # elif c["type"] == "tabulated n" :
        #    coef = numpy.fromstring(c["data"], sep=" ")
        #    l = len(coef)
        #    coef.shape = int(l / 2), 2
        #    return Tabulated_N(coef)

    else:
        # The else belongs to the for, and is used to check the break was
        # not used
        raise ModelNotImplemented
