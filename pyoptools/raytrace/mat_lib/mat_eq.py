#!/usr/bin/env python
# -*- coding: utf-8-*-

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
from math import sqrt


class ModelNotImplemented(Exception):
    """Indicates that the model to calculate the index of refraction is not
    implemented yet"""
    pass


class Material:
    ''' Base class to define an optical material. It receives the YML text as
    defined in https://refractiveindex.info/about
    '''
    def __init__(self, coef, cl, nd=None, vd=None):
        self.__nd__ = nd
        self.__vd__ = vd

        # Check if the number of coeficients agrees with the dispersion
        # formulas definition

        if cl is None:
            self.__coef__ = coef
        else:
            self.__coef__ = coef.copy()
            self.__coef__.resize(cl)

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
            return (self.n(0.5893)-1) / (self.n(0.4861)-self.n(0.6563))


class Sellmeier(Material):
    """Class that define a material that complies with the Sellmeier dispersion
    model as defined at
    https://refractiveindex.info/database/doc/Dispersion%20formulas.pdf
    """
    def __init__(self, coef, nd=None, vd=None):
        Material.__init__(self, coef, 17, nd, vd)

    def n(self, l=.58929):
        n2 = 1+self.__coef__[0]
        for i in range(1, 17, 2):
            n2 = n2+(self.__coef__[i]*l**2/(l**2-self.__coef__[i+1]**2))
        return sqrt(n2)


class Sellmeier_2(Material):
    """Class that define a material that complies with the Sellmeier_2
    dispersion model as defined at
    https://refractiveindex.info/database/doc/Dispersion%20formulas.pdf
    """
    def __init__(self, coef, nd=None, vd=None):
        Material.__init__(self, coef, 17, nd, vd)

    def n(self, l=.58929):
        n2 = 1+self.__coef__[0]
        for i in range(1, 17, 2):
            n2 = n2+(self.__coef__[i]*l**2/(l**2-self.__coef__[i+1]))
        return sqrt(n2)


class Polynomial(Material):
    """Class that define a material that complies with the Polynomial
    dispersion model as defined at
    https://refractiveindex.info/database/doc/Dispersion%20formulas.pdf
    """
    def __init__(self, coef, nd=None, vd=None):
        Material.__init__(self, coef, 17, nd, vd)

    def n(self, l=.58929):
        n2 = self.__coef__[0]
        for i in range(1, 17, 2):
            n2 = n2+self.__coef__[i]*l**self.__coef__[i+1]
        return sqrt(n2)


class RefractiveIndex_Info(Material):
    def __init__(self, coef, nd=None, vd=None):
        raise ModelNotImplemented


class Cauchy(Material):
    def __init__(self, coef, nd=None, vd=None):
        Material.__init__(self, coef, 11, nd, vd)

    def n(self, l=.58929):
        n_ = self.__coef__[0]
        for i in range(1, 11, 2):
            n_ = n_+self.__coef__[i]*l**self.__coef__[i+1]
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


__models__ = [Sellmeier, Sellmeier_2, Polynomial, RefractiveIndex_Info, Cauchy,
              Gases,  Herzberger, Retro,  Exotic]


def from_yml(filename):
    """Create a material instance from a YML file as defined at
       https://refractiveindex.info/about
    """
    with open(filename) as f:
        mat = yaml.load(f, Loader=yaml.FullLoader)

    for c in mat["DATA"]:
        if "formula" in c["type"]:
            break
    else:
        # The else belongs to the for, and is used to check the break was
        # not used
        raise ModelNotImplemented

    fn = int(c["type"].split()[1]) - 1

    coef = numpy.fromiter(c["coefficients"].split(), dtype=numpy.float)

    # Sometimes SPECS does not exists
    nd = mat.get("SPECS", {}).get("nd", None)
    vd = mat.get("SPECS", {}).get("Vd", None)
    return __models__[fn](coef, nd, vd)
