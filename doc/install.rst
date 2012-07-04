Installing pyoptools (stable release)
=====================================


Note the stable release is currently being prepared.


UBUNTU (Latest versions)
------

1. Install ipython, build-essential, python-numpy, python-scipy, python-matplotlib, python-wxversion, python-dev, python-opengl from the repository::
    
    ~$ sudo apt-get install ipython build-essential python-numpy python-scipy python-matplotlib python-wxversion python-dev python-opengl


2. Download and install Cython 0.13 following the instructions found on the webpage (http://www.cython.org)

3. Download the pyoptools repository::

    ~$ hg clone https://code.google.com/p/pyoptools/ -r stable

4. Enter the pyoptools folder, build and install the pyoptools::

    ~$ cd pyoptools
    ~$ sudo python setup.py install
    
    


