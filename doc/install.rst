Installing pyOpTools
====================

Getting it
----------

pyOpTools can be downloaded from the project GitHub repository at:

    https://github.com/cihologramas/pyoptools


Installing it in Debian 11 (should work on any Debian derivative)
-----------------------------------------------------------------

pyOpTools is being developed almost exclusively using "Debian derivative" Linux
distributions, so this installation procedure is the most tested and reliable. To
generate the pyoptool .deb package you just need to run the following command in the project root::

    make deb

This will create a .deb package outside the pyOpTool's source root tree. This can be installed using the command::

    dpkg -i python3-pyoptools_<version>_<platform>.deb

where version/platform should be adjusted accordingly.  


Installing pyOpTools inside a python virtualenv
-----------------------------------------------

After a Python3 virtualenv is created and activated, pyOpTools can be installed by running the following command, inside the project root::
  
    pip install -r requirements.txt
    python3 setup.py install
    
.. _visualizing_pyoptools_in_jupyter:

Visualizing pyOpTools simulations in a JupyterLab notebook
----------------------------------------------------------

The easiest way to  use pyOpTools together with
`jupyterlab <https://jupyter.org>`_, is to install pyOpTools
inside a python virtualenv, and then install jupyterlab on it also::

    pip install jupyterlab

After that the jupyter plugin pythreejs must be installed and enabled::

    pip3 install pythreejs
    jupyter labextension install jupyter-threejs

after this is done, you will be able to visualize the simulations using the
:func:`~pyoptools.gui.ipywidgets.Plot3D` command.

The plot window is interactive. Using the mouse, it is possible to rotate
the image by click and drag with the left button, zoom by using the scroll
wheel and translate by click and drag with the right button.

.. note ::
    When installing pyoptools and jupyterlab in a debian 11 machine inside a
    virtualenv, it was necesary to install also importlib-resources::
      
      pip install importlib-resources

