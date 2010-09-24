#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
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

from os import listdir, walk
from os.path import join, split,  splitext
#from enthought.traits.api import Float, HasTraits, String
from numpy import sqrt

from pkg_resources import resource_stream, resource_filename, \
                        resource_string,  resource_listdir

#mat_config = resource_filename(__name__, 'data/')
mat_config = resource_filename("pyoptools.raytrace.mat_lib", 'data')
# Get Datafiles
# TODO: Include user files
datafiles=[]
for root, dirs, files in walk(mat_config):
    for name in files:
        datafiles.append(join(root, name))

# Store material library names
libnames=[]
for fname in datafiles:
    f=open(fname, 'r')
    a=f.readline()
    if a[0:6]=="RTGL-0":
        libnames.append((splitext((split(fname))[1])[0], fname))
    f.close()
del(datafiles)


print "Hay que revisar la libreria de materiales material.py"
#######################

class Material:
    ''' Class to define an optical material.

    When the module is loaded, the spreadsheet is read, and the
    attribures created at runtime
    '''

    #~ # Glass Reference
    #~ ref=String()
#~ 
    #~ # Refraction index for the d line
    #~ nd=Float()
    #~ 
    #~ # Abbe's number for the d line
    #~ vd=Float()
#~ 
    #~ #Dispersion constants     
    #~ B1=Float()
    #~ B2=Float()
    #~ B3=Float()
    #~ C1=Float()
    #~ C2=Float()
    #~ C3=Float()
    #~ 
    def __init__(self,ref="",nd=0.,vd=0.,B1=0.,B2=0.,B3=0.,C1=0.,C2=0.,C3=0.):
        self.nd=nd
        self.vd=vd
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




# Creating the material dictionaries
liblist=[]
for libn, file in libnames:
    tdict={}
    liblist.append((libn, tdict))
    globals()[libn]=tdict
    fo=open(file, 'r')
    data=fo.readlines()
    fo.close()
    keys=data[1].split()
    for inst in data[2:]:
        tmat=Material()
        idata= inst.split()
        tdict[idata[0].strip("\"")]=tmat
        d={}
        for n, k in enumerate(keys):
            if k in dir(tmat):
                if idata[n][0]!="\"":
                    val=float(idata[n])
                else:
                    val=idata[n].strip("\"")
                setattr(tmat,k,val)
            #else:
            #    print "Error loading the library %s from the file %s %s"%(libn, file, k)
    

del(libnames)
def find_material(material):
    """Search for a material in all the libraries
    
    This function prints all the libraries that contain the material 
    
    Arguments:
    
    material
        String with the material name
    """
    for libn,tdict in liblist:
       if tdict.has_key(material):
          print libn 

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
       if tdict.has_key(material):
          return tdict[material]
    return None


def mat_list():
    for libn,tdict in liblist:
       print libn,  tdict.keys()
    
