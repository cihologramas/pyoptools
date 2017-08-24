#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
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
#------------------------------------------------------------------------------

'''
Material class definition, and helper functions used to load the
constants of the dispersion formula to be used in the calculation of
the refraction index.

It also create libraries containing the materials defined in the directory data
so the materials can be used as a predefined.
'''

import six
from os import listdir
from os.path import join, split,  splitext, isfile, expanduser,isdir
#from enthought.traits.api import Float, HasTraits, String
from numpy import sqrt
from six.moves import configparser, reduce
from pkg_resources import resource_stream, resource_filename, \
                        resource_string,  resource_listdir

#mat_config = resource_filename(__name__, 'data/')
mat_config = resource_filename("pyoptools.raytrace.mat_lib", 'data')


# Get Datafiles from the system
libnames=[]

for fname in listdir(mat_config):
    fpath = join(mat_config,fname)
    if not isfile(fpath):
        continue
    if splitext(fpath)[1] != ".mat":
        continue
    libnames.append((splitext(fname)[0], fpath))

# get datafiles from the user repeated materials will be overwritten
local_mat_config = reduce(join,[expanduser("~"),".pyoptools","material"])
if isdir(local_mat_config):
    for fname in listdir(local_mat_config):
        fpath = join(local_mat_config,fname)
        if not isfile(fpath):
            continue
        if splitext(fpath)[1] != ".mat":
            continue
        libnames.append((splitext(fname)[0], fpath))




class Material:
    ''' Class to define an optical material.
    '''
    def __init__(self,ref="",nd=0.,vd=0):
        self.nd=nd
        self.vd=vd



class Sellmeier(Material):
    """
    Clase que define los materiales de acuerdo a la formula de Sellmeier
    """
    def __init__(self,nd=0.,vd=0.,B1=0.,B2=0.,B3=0.,C1=0.,C2=0.,C3=0.):
        Material.__init__(self,nd,vd)

        self.B1=B1
        self.B2=B2
        self.B3=B3

        self.C1=C1
        self.C2=C2
        self.C3=C3

    def n(self, l=.58929):

        ''' Returns the refraction index for a specific wavelength.

        if no wavelength is given, it returns the refraction index at
        l=.58929 um

        Note: the wavelength should be given in micrometers
        '''

        return  sqrt(1.+
                     (self.B1*l**2)/(l**2-self.C1)+
                     (self.B2*l**2)/(l**2-self.C2)+
                     (self.B3*l**2)/(l**2-self.C3))


class Schott(Material):
    """
    Clase que define los materiales de acuerdo a la formula de Schott
    Tomado de http://glassbank.ifmo.ru/eng/help.php#disp_formula
    """
    def __init__(self,nd=0.,vd=0.,C1=0.,C2=0.,C3=0.,C4=0.,C5=0.,C6=0.):
        Material.__init__(self,nd,vd)

        self.C1=C1
        self.C2=C2
        self.C3=C3

        self.C4=C4
        self.C5=C5
        self.C6=C6

    def n(self, l=.58929):

        ''' Returns the refraction index for a specific wavelength.

        if no wavelength is given, it returns the refraction index at
        l=.58929 um

        Note: the wavelength should be given in micrometers
        '''
        return  sqrt(self.C1+ self.C2*l**2+self.C3/(l**2)+
                     self.C4/(l**4)+ self.C5/(l**6)+
                     self.C6/(l**8))


class Cauchy(Material):
    def __init__(self,nd,vd,A,B,C,D):
        Material.__init__(self,nd,vd)
        #E n wikipedia no usan el termino A, pero en https://refractiveindex.info
        # Si lo usan y por eso se colocó

        self.A = A
        self.B = B
        self.C = C
        self.D = D

    def n(self,l=.58929):

        ''' Returns the refraction index for a specific wavelength.

        if no wavelength is given, it returns the refraction index at
        l=.58929 um

        Note: the wavelength should be given in micrometers
        '''

        return (self.A*l**2+
                self.B+
                self.C*l**(-2)+
                self.D*l**(-4))



# Creating the material dictionaries
__liblist__=set()

for libn, fname in libnames:
    #Crear el diccionario "libn" donde se van a guardar los materiales
    # Cada catalogo tiene su diccionario. Si los nombres de los catalogos se repiten,
    # entre el sistema y la carpeta del usuario, los materiales se adicionan a un mismo catalogo
    __liblist__.add(libn)
    if libn in globals():
        tdict =globals()[libn]
    else:
        tdict={}
        globals()[libn]=tdict
    matconfig = configparser.ConfigParser()
    matconfig.read(fname)

    for matref in matconfig.sections():
        matmodel = matconfig.get(matref,"model")

        #### Definir como se leen los diferentes modelos
        if matmodel == "Sellmeier":
            nd = matconfig.getfloat(matref,"nd")
            vd = matconfig.getfloat(matref,"vd")
            b1 = matconfig.getfloat(matref,"b1")
            b2 = matconfig.getfloat(matref,"b2")
            b3 = matconfig.getfloat(matref,"b3")
            c1 = matconfig.getfloat(matref,"c1")
            c2 = matconfig.getfloat(matref,"c2")
            c3 = matconfig.getfloat(matref,"c3")
            tmat=Sellmeier(nd,vd,b1,b2,b3,c1,c2,c3)
        elif matmodel == "Schott":
            nd = matconfig.getfloat(matref,"nd")
            vd = matconfig.getfloat(matref,"vd")
            c1 = matconfig.getfloat(matref,"c1")
            c2 = matconfig.getfloat(matref,"c2")
            c3 = matconfig.getfloat(matref,"c3")
            c4 = matconfig.getfloat(matref,"c4")
            c5 = matconfig.getfloat(matref,"c5")
            c6 = matconfig.getfloat(matref,"c6")
            tmat = Schott(nd,vd,c1,c2,c3,c4,c5,c6)

        elif matmodel == "Cauchy":
            nd = matconfig.getfloat(matref,"nd")
            vd = matconfig.getfloat(matref,"vd")
            a = matconfig.getfloat(matref,"a")
            b = matconfig.getfloat(matref,"b")
            c = matconfig.getfloat(matref,"c")
            d = matconfig.getfloat(matref,"d")
            tmat = Cauchy(nd,vd,a,b,c,d)

        else:
            print("Model not available")
            continue
        tdict[matref]=tmat
liblist=[]
for i in __liblist__:
    liblist.append((i,globals()[i]))
#del(libnames)

def find_material(material):
    """Search for a material in all the libraries

    This function prints all the libraries that contain the material

    Arguments:

    material
        String with the material name
    """
    retv=[]
    for libn,tdict in liblist:
       if material in tdict:
          retv.append(libn)
    return retv


def get_material(material):
    """Search for a material in all the libraries

    This function search in all the material libraries, and return the
    first instance found, that matches the name of the material requested.
    If no material found, returns None
    Arguments:

    material
        String with the material name
    """
    for libn,tdict in liblist:

        if material in tdict:
            return tdict[material]
    print (material, " not found")
    raise KeyError



def mat_list():
    for libn in liblist:
        tdict=globals()[libn]
        print (libn,  tdict.keys())
