Installing pyoptools
====================


UBUNTU 10.04 (Lucid Lynx)
------

1. Install ipython, build-essential, python-numpy, python-scipy, python-matplotlib, python-wxversion, python-dev, python-opengl from the repository::
    
    ~$ sudo apt-get install ipython build-essential python-numpy python-scipy python-matplotlib python-wxversion python-dev python-opengl


2. Install Cython 0.13 (www.cython.org)

3. Download the pyoptools repository::

    ~$ svn checkout http://pyoptools.googlecode.com/svn/trunk/ pyoptools

4. Enter the pyoptools folder, build and install the pyoptools::

    ~$ cd pyoptools
    ~$ sudo python setup.py install
    
    


