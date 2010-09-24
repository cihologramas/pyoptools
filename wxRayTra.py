#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from pkg_resources import require

from matplotlib import use
use('WXAgg') 

# Use wxversion to import wxpython 
MIN_WX_VERSION  = '2.8'
GET_WXPYTHON    = 'Get it from http://www.wxpython.org!'

#Define a variable to indicate this is a wxRayTra console
import os
os.environ['shell'] = 'wxRayTra'

from sys import exit, modules
try:
    import wxversion
    if modules.has_key('wx') or modules.has_key('wxPython'):
        pass    #probably not the first call to this module: wxPython already loaded
    else:
        #wxversion.select('2.6')
        wxversion.ensureMinimal(MIN_WX_VERSION)
    
except ImportError:
    #the old fashioned way as not everyone seems to have wxversion installed
    try:
        import wx
        if wx.VERSION_STRING < MIN_WX_VERSION:
            print 'You need to upgrade wxPython to v%s (or higher) to run pytbl.'%MIN_WX_VERSION
            print GET_WXPYTHON
            exit()
    except ImportError:
            print "Error: pytbl requires wxPython, which doesn't seem to be installed."
            print GET_WXPYTHON
            exit()
    print 'Warning: the package python-wxversion was not found, please install it!'
    
    
import wx
from wxRayTrace.gui import wxrayos
from wxRayTrace.gui.events import register_frame, SendEvent
from os.path import abspath

#from gui.glplotframe import *

#import logging.config
        

# Inicio Callbacks
def evtTextEnter(self,event):
    command=self.text_ctrl_1.GetValue()
    self.text_ctrl_2.AppendText(command+"\n")
    interpreter.runcode(command)
    self.text_ctrl_1.SetValue("")

def evtText(self,event):
    pass

def evt_run_file(self, event):
    from wx import FileDialog,ID_OK
    FD=FileDialog(self, "Choose a file", defaultDir = "", defaultFile = "",\
                     wildcard = "(*.py)|*.py", style = 0)
    if FD.ShowModal()==ID_OK:
        path=abspath(FD.GetPath())
        mpath=path.replace("\\","/")
        execute(self.rtshell,"execfile(\""+mpath+"\")")
    event.Skip()



def execute(shell,order):
    shell.text_ctrl.write(order)
    shell.setCurrentState('DO_EXECUTE_LINE')
    shell.stateDoExecuteLine()
    

if __name__ == "__main__":
    app = wx.PySimpleApp(redirect=False)
    
    
    wx.InitAllImageHandlers()
    wxrayos.wxRTFrame.evtTextEnter=evtTextEnter
    wxrayos.wxRTFrame.evt_run_file=evt_run_file
    wxrayos.wxRTFrame.evtText=evtText
    
    frame1 = wxrayos.wxRTFrame(None, -1, "")
    app.SetTopWindow(frame1)
    frame1.Show()

    # Register the application and main frame into the event manager gui.events
    register_frame(frame1)
    #logging.config.fileConfig('logger.ini')
    
    shell=frame1.rtshell
    
    #Different versions of ipython use diferent method
    try:
        shell.IP.updateNamespace(globals())
    except:
        shell.IP.update_namespace(globals())
    execute(shell,"from pyoptools import *")
    #execute(shell,"from wxRayTrace.gui.plotutils import *")
    #glFrame()
    app.MainLoop()
#
