#!/usr/bin/env python
# -*- coding: utf-8 -*-
#   Copyright (c) 2007, 2008, 2009,2010 Ricardo Amézquita Orozco
#   <ramezquitao@unal.edu.co>,
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following disclaimer
#     in the documentation and/or other materials provided with the
#     distribution.
#   * Neither the name of the  nor the names of its
#     contributors may be used to endorse or promote products derived from
#     this software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

'''
Package containing modules and submodules defining an *API* for optical
raytracing and wave propagation calculations.
'''

# Import all pyoptools packages

from pyoptools.misc.cmisc import *
from pyoptools.misc.definitions import *
#from misc.frft import *
from pyoptools.misc.lsq import *
from pyoptools.misc.pmisc import *
from pyoptools.misc.plist import *
from pyoptools.misc.Poly2D import *
from pyoptools.misc.resources import *

from pyoptools.raytrace.calc import *
from pyoptools.raytrace.comp_lib import *
from pyoptools.raytrace.component import *
from pyoptools.raytrace.library import *
from pyoptools.raytrace.mat_lib import *
from pyoptools.raytrace.ray import *
from pyoptools.raytrace.shape import *
from pyoptools.raytrace.surface import *
from pyoptools.raytrace.system import *

from pyoptools.wavefront.field import *
from pyoptools.wavefront.calc import *
from pyoptools.wavefront.psurfrep import *
from pyoptools.wavefront.zernike import *


#
#
# Import graphic packages This should be imported somewhere else
from pyoptools.gui.plotutils import *

# This module has problems with MESA in buster is disabled for the moment
#from pyoptools.gui.ipynbplotutils import *

# Module implemented using pythreejs
from pyoptools.gui.ipywidgets import *
