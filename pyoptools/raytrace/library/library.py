import ConfigParser as cp

from os import listdir, walk
from os.path import join, split, isdir, splitext,basename
from numpy import sqrt

from pkg_resources import resource_stream, resource_filename, \
                        resource_string,  resource_listdir
 
#import all the predefined components to the library and materials to the
#library
import pyoptools.raytrace.comp_lib as CL
from pyoptools.raytrace.mat_lib import get_material

lenslib = resource_filename("pyoptools.raytrace.library", '')

# Get Datafiles
# TODO: Include user files

#Save the available component classes names defined in CL
_av_comp=[]

for cl in dir(CL):
    if cl[0].isupper():
        _av_comp.append(cl)
        

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
        print "Loading component library",libname," from files ",filename
        Library(filename,libname=libname)
             
