import sys, os, stat, commands
from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
except:
    print "You don't seem to have Cython installed. Please get a"
    print "copy from www.cython.org and install it"
    sys.exit(1)

#from setuptools import find_packages#,setup


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
        include_dirs = ["."],   # adding the '.' to include_dirs is CRUCIAL!!
        extra_compile_args = ["-O3", "-Wall"],
        #extra_link_args = ['-g'],
        #libraries = ["dv",],
        )


#Need to create a pyOpTools package
extNames = scandir("./")

# and build up the set of Extension objects
extensions = [makeExtension(name) for name in extNames]



setup(
        name =  "pyOpTools",
        #~ #~#version = "0.1.0",
        packages=findpackages("./"),
        scripts=['wxRayTra.py'],
        install_requires = ['python-numpy',
                            'cython'],
        package_data= {
        'pyoptools.raytrace.mat_lib':['data/*.rtgl'],
        'pyoptools.raytrace.library':['Edmund/*.cmp'],
        },
        author= 'Ricardo Amezquita Orozco',
        author_email='ramezquitao@cihologramas.com',
        description='Optical ray tracing simulation system',
        license='BSD',
        url='https://trac.cihologramas.com/trac/wxRayTrace',
        ext_modules=extensions,
        cmdclass = {'build_ext': build_ext},
    )

