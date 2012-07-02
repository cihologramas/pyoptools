#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''Module defining the comunication between the GUI and the processing thread
    still need to see if this is necesary
'''

import wx

# List of possible events

# Frame used to manage the events
__Frame__=None

def register_frame(frame):
    """Function to be called to register the frame into the comunication
    Module
    """
    global __Frame__
    __Frame__=frame
    EVT_NEW_COMMAND( frame, event_handler )
    
        
def event_handler(event):
    ##
    event.func(*event.args, **event.kwrds)
    return True



              
################################################
#Definicion de eventos nuevos para la interfaz

#Evento que crea una nueva ventana para colocar imagenes

NEW_COMMAND = wx.NewEventType() 
 
def EVT_NEW_COMMAND( window, function ): 
    """Utilizado para comunicar el programa de visualizacion y el
    interprete de python """ 
    window.Connect( -1, -1, NEW_COMMAND, function ) 
 
class NewCommandEvent(wx.PyCommandEvent): 
    eventType = NEW_COMMAND 
    def __init__(self, windowID, event_type,event_data): 
        wxPyCommandEvent.__init__(self, self.eventType, windowID)
        
    def Clone( self ): 
        self.__class__( self.GetId() )
        

def SendEvent(func, *args, **kwds):
    event = NewCommandEvent(__Frame__.GetId(),event_type,event_data)
    event.func=func
    event.args=args
    event.kwrds=kwds
    __Frame__.GetEventHandler().AddPendingEvent( event )

#def wxExec(func, *args, **kwds):
#    SendEvent(f)
