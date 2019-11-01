#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import six
from six.moves import configparser as cp
from six.moves import reduce

from os import listdir, walk
from os.path import join, split, isdir, splitext,basename, expanduser
from numpy import sqrt
import inspect

from pkg_resources import resource_stream, resource_filename, \
                        resource_string,  resource_listdir

#import all the predefined components to the library and materials to the
#library
import pyoptools.raytrace.comp_lib as CL
from pyoptools.raytrace.mat_lib import get_material
from pyoptools.raytrace.system.system import System
from pyoptools.raytrace.component.component import Component

lenslib = resource_filename("pyoptools.raytrace.library", '')


###############################################################################
#Save the available component classes names defined in CL. This are the classes
#defined in pyoptools.raytrace.comp_lib
_av_comp=[]

#for cl in dir(CL):
#    if cl[0].isupper():
#        _av_comp.append(cl)

# Replace the check for something more robust than checking if the firs character
# is a capital letter.
for cl in dir(CL):
    c=getattr(CL,cl)
    if not inspect.isclass(c):
        continue
    if issubclass(c,(Component,System)):
        _av_comp.append(cl)
################################################################################

#Class to manage component libraries libraries
class Library:
    def __init__(self, filename, libname=None):

        self.parser=cp.ConfigParser()
        self.parser.read(filename)

        if libname != None:
            globals()[libname]=self
    def get(self, cmp):
        global _av_comp
        assert self.parser.has_section(cmp), "Library error in component %s, type not defined"%(cmp,)
        type=self.parser.get(cmp,"type")

        # Check that the class name is defined in the comp_lib module
        if not type in _av_comp:
            raise TypeError("Component type "+type+" not defined in comp_lib")
        options=self.parser.options(cmp)

        #Create the constructor string
        strobj="CL."+type+"("
        for opt in options:
            #Remove options not needed for the class constructors
            if opt=="type":
                pass
            elif opt=="description":
                pass
            elif opt=="glass_catalogs":
                pass

            #Check if the material exists in the material library and
            elif opt[0:3]=="mat":
                mat=self.parser.get(cmp,opt).upper()
                m=get_material(mat)
                if m==None:
                    raise TypeError("material "+mat+" not defined in mat_lib")
                strobj=strobj+opt+"=get_material(\""+mat+"\") ,"
            else:
                strobj=strobj+opt+"="+self.parser.get(cmp,opt)+" ,"
        strobj=strobj[:-1]+")"
        # strobj
        return eval(strobj)

    def parts(self):
        return self.parser.sections()

# Read and load all the component libraries

listdir(lenslib)

dirs=[]
#Get the directories. Each directory is a library.
for  i in listdir(lenslib):
    di=join(lenslib,i)
    if isdir(di):
        dirs.append(di)

#Get the local directories
# Check if the local library path exists
locallenslib = reduce(join,[expanduser("~"),".pyoptools","library"])
if isdir(locallenslib):
    for i in listdir(locallenslib):
        di=join(locallenslib,i)
        if isdir(di):
            dirs.append(di)


for di in dirs:
    fnames=listdir(di)
    # Retain only the .cmp files
    libfiles=[]
    for f in fnames:
        n,ext=splitext(f)
        if ext==".cmp":
            libfiles.append(f)

    if len(libfiles)>0:
        libname=basename(di)
        filename=[ join(di,fn) for fn in libfiles]
        print ("Loading component library",libname," from files ",filename)
        Library(filename,libname=libname)


#########################################################################################################################
## Mientras se organiza voy a colocar todos los archivos correspondientes a librerias acá

# Se crearon para poder leer facilmente archivos de thorlabs
def convert(d):
    try:
        return int(d)
    except ValueError:
        try:
            return float(d)
        except ValueError:
            return d


def zmx_parse(data):
    """Función que lee e interpreta un archivo zmx de zemax.

    - Solo funciona para lentes esfericas y dobletes.
    - No tiene ne cuenta recubrimientos antireflectivos
    - Asume que las medidas están en milimetros. Hay que arreglar esto
    """

    lines=data.splitlines()

    # Interpretar el encabezado
    while True:
        line=lines.pop(0)

        if lines[0].startswith("SURF"):
            break

    # Separar las superficies en una lista de diccionarios
    surflist=[]
    for line in lines:
        if line.startswith("SURF"):
            surflist.append(dict())
            continue
        line=line.lstrip()
        code=line[:4]
        data=line[5:].split()
        data=[convert(d) for d in data]
        surflist[-1][code]=data


    #Eliminar el plano objeto y el plano imagen

    surflist.pop(0)
    surflist.pop()


    # Identificar el tipo de lentes a partir de el numero de superficies
    # validas

    ns=len(surflist)

    if ns==2: #Lentes normales
        c0=surflist[0]["CURV"][0]
        c1=surflist[1]["CURV"][0]
        d0=surflist[0]["DISZ"][0]
        r0=surflist[0]["DIAM"][0]
        r1=surflist[1]["DIAM"][0]
        g0=surflist[0]["GLAS"][0]
        m0=get_material(g0)


        ##c|Verificar que las superficies son iguales, si no emitir un error
        assert r0==r1

        return CL.SphericalLens(r0,d0,c0,c1,material=m0)

    if ns==3: #Dobletes
        c0=surflist[0]["CURV"][0]
        c1=surflist[1]["CURV"][0]
        c2=surflist[2]["CURV"][0]
        d0=surflist[0]["DISZ"][0]
        d1=surflist[1]["DISZ"][0]
        r0=surflist[0]["DIAM"][0]
        r1=surflist[1]["DIAM"][0]
        r2=surflist[2]["DIAM"][0]
        g0=surflist[0]["GLAS"][0]
        g1=surflist[1]["GLAS"][0]
        m0=get_material(g0)
        m1=get_material(g1)

        #Verificar que las superficies son iguales, si no emitir un error
        assert r0==r1 and r1== r2

        return CL.Doublet(r0,c0,c1,c2, d0,d1,m0,m1)

    elif ns==4 and "GLAS" not in surflist[1]: #Dobletes con espaciado en Aire
        c0=surflist[0]["CURV"][0]
        c1=surflist[1]["CURV"][0]
        c2=surflist[2]["CURV"][0]
        c3=surflist[3]["CURV"][0]

        d0=surflist[0]["DISZ"][0]
        d1=surflist[1]["DISZ"][0]
        d2=surflist[2]["DISZ"][0]



        r0=surflist[0]["DIAM"][0]
        r1=surflist[1]["DIAM"][0]
        r2=surflist[2]["DIAM"][0]
        r3=surflist[3]["DIAM"][0]

        g0=surflist[0]["GLAS"][0]
        # g1=surflist[1]["GLAS"][0] Este parametro no existe. Es aire
        g2=surflist[2]["GLAS"][0]

        m0=get_material(g0)
        #m1=get_material(g1)
        m1=get_material(g2)

        #Verificar que las superficies son iguales, si no emitir un error

        #assert r0==r1 and r1== r2 and r2 ==r3 Esto no siempre se cumple


        return CL.AirSpacedDoublet(r0,c0,c1,c2,c3, d0,d1,d2,m0,m1)

    else:
        for i in surflist:
            print("*", i)
        raise ValueError # Esto toca arreglarlo y generar un error que realmente indique que está pasando
def zmx_read(fn):

    f=open(fn,"rU")
    data=f.read()
    return zmx_parse(data)

##Codigo tomado y modificado de https://github.com/jordens/rayopt/blob/master/rayopt/zemax.py
from struct import Struct
import numpy as np

def zmf2dict(fn):
    """Función que lee una librería de Zemax (archivo con terminación zmf), y genera un diccionario con las descripciones
    de cada componente. La llave es la referencia de cada componente
    """
    f=open(fn,"r")
    rd={}
    head = Struct("<I")
    lens = Struct("<100sIIIIIIIdd")
    shapes = "?EBPM"
    version, = head.unpack(f.read(head.size))
    assert version in (1001, )
    while True:
        li = f.read(lens.size)
        if len(li) != lens.size:
            if len(li) > 0:
                print(f, "additional data", repr(li))
            break
        li = list(lens.unpack(li))
        li[0] = li[0].decode("latin1").strip("\0")
        li[3] = shapes[li[3]]
        description = f.read(li[7])
        assert len(description) == li[7]
        description = zmf_obfuscate(description, li[8], li[9])
        description = description.decode("latin1")
        assert description.startswith("VERS {:06d}\n".format(li[1]))
        rd[li[0]]=description
    return rd


def zmf_obfuscate(data, a, b):
    iv = np.cos(6*a + 3*b)
    iv = np.cos(655*(np.pi/180)*iv) + iv
    p = np.arange(len(data))
    k = 13.2*(iv + np.sin(17*(p + 3)))*(p + 1)
    k = (int(("{:.8e}".format(_))[4:7]) for _ in k)
    data = np.fromstring(data, np.uint8)
    data ^= np.fromiter(k, np.uint8, len(data))
    return data.tostring()
