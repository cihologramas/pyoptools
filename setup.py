import sys, os, commands
from distutils.core import setup
from distutils.extension import Extension
try:
    from Cython.Distutils import build_ext
except:
    print "You don't seem to have Cython installed. Please get a"
    print "copy from www.cython.org and install it"
    sys.exit(1)

#Look for paths containig arrayobject.h 
#Necessary for non-unix systems 	
def contains_arrayobject_h(path):
    """
    Returns True iff the python path string contains the arrayobject.h
    include file where it is supposed to be.
    """
    f=False
    try:
        s=os.stat(os.path.join(path, 'numpy', 'core', 'include', \
                               'numpy', 'arrayobject.h'))
        f=True
    except OSError:
        pass
    return f

# scan the directory for extension files, converting
# them to extension names in dotted notation
def scandir(dir, files=[]):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if os.path.isfile(path) and path.endswith(".pyx"):
            files.append(path.replace(os.path.sep, ".")[2:-4])
        elif os.path.isdir(path):
            scandir(path, files)
    return files


def findpackages(dir,files=[], ):
    for file in os.listdir(dir):
        if file!="build":
            path = os.path.join(dir, file)
            if os.path.isdir(path):# and path.endswith(".py"):
                for file1 in os.listdir(path):
                    if file1=="__init__.py":
                        files.append(path.replace(os.path.sep, ".")[2:])
                findpackages(path, files)
    return files

# generate an Extension object from its dotted name
def makeExtension(extName):
    extPath = extName.replace(".", os.path.sep)+".pyx"
    pxdPath = extName.replace(".", os.path.sep)+".pxd"
    if os.path.isfile(pxdPath):
        flist=[extPath,pxdPath]
    else:
        flist=[extPath]
    return Extension(
        extName,
        [extPath],
        include_dirs = [".",include_numpy_array],   # adding the '.' to include_dirs is CRUCIAL!!
        extra_compile_args = ["-O3", "-Wall"],
        #extra_link_args = ['-g'],
        #libraries = ["dv",],
        )

#Check the availability of arrayobject.h		
valid_paths = filter(contains_arrayobject_h, sys.path)
if len(valid_paths) == 0:
    print "No paths in the python path contain numpy/arrayobject.h"
    sys.exit(0)

# The base path is by default the first python path with arrayobject.h in it.
include_numpy_array=valid_paths[0]

if len(valid_paths) > 1:
    print "There are several valid include directories containing numpy/arrayobject.h"
    l=[('%d: %s' % (i+1, valid_paths[i])) for i in xrange(0, len(valid_paths))]
    s = -1
    print string.join(l, '\n')
    # Prompt the user with a list of selections.
    while not (s >= 1 and s <= len(valid_paths)):
        s = input('Selection [default=1]:' % s)
        if s == '':
            s = 1
        else:
            s = int(s)
    include_numpy_array=valid_paths[s-1]

# Add the children directory path suffix to the base path.
include_numpy_array=os.path.join(include_numpy_array, 'numpy', 'core', \
                                 'include')

#Need to create a pyOpTools package
extNames = scandir("./")

# and build up the set of Extension objects
extensions = [makeExtension(name) for name in extNames]



setup(
        name =  "pyOpTools",
        #~ #~#version = "0.1.0",
        packages=findpackages("./"),
#        scripts=['wxRayTra.py'],
        install_requires = ['python-numpy',
                            'cython'],
        package_data= {
        'pyoptools.raytrace.mat_lib':['data/*.rtgl'],
        'pyoptools.raytrace.library':['Edmund/*.cmp'],
        },
        author= 'Ricardo Amezquita Orozco',
        author_email='ramezquitao@unal.edu.co',
        description='Optical ray tracing simulation system',
        license='BSD',
        url='http://code.google.com/p/pyoptools/',
        ext_modules=extensions,
        cmdclass = {'build_ext': build_ext},
    )
