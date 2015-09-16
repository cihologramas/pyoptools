# cython: profile=True
cimport numpy as np
from numpy import array,float64
import six

from pyoptools.misc.picklable.picklable cimport Picklable
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.component.component cimport Component
from pyoptools.raytrace.system.system cimport System
cdef class plist(Picklable):
    """
    Class used to define a part list. o in
    
    It is based on the dict interface.
    
    Each item must have a form (O,P,D), where O is an object, P is a position
    vector (array, tuple or list), and D is a direction vector (array,tuple,list).
    
    This class is based on the dict interface, to be able to use string keys instead 
    of just numbers. It does not inherit from dict, because of cython limitation in
    multiple inheritance.
    """
    
    def __init__(self, a=None):
        """
        Create the list to save the objects, positions and rotations
        """
        
        #Register class attributes on Picklable superclass 
        Picklable.__init__(self,"_buf")
        
        self._buf={}
        cdef np.ndarray[np.float64_t, ndim=1] ap,ar
        
        # If a is a list, each key will be a number
        if isinstance(a,(list,tuple)):
            for i in a:
                if len(i)==3:
                    O,p,r=i
                
                    ap=array(p,dtype=float64)
                    ar=array(r,dtype=float64)
                    assert len(ap)==3, "The position vector must be an array, or a list, or a tuple of len 3"
                    assert len(ar)==3, "The rotation must be an array, or a list, or a tuple of len 3"
                    self.append((O,ap,ar))
                else:
                    O,p,r,k=i
                    ap=array(p,dtype=float64)
                    ar=array(r,dtype=float64)
                    assert len(ap)==3, "The position vector must be an array, or a list, or a tuple of len 3"
                    assert len(ar)==3, "The rotation must be an array, or a list, or a tuple of len 3"
                    self[k]=(O,ap,ar)
                    
                    
        elif isinstance(a,(dict,plist)):
            for k, v in a.iteritems():
                O,p,r = v
                ap=array(p,dtype=float64)
                ar=array(r,dtype=float64)
                assert len(ap)==3, "The position vector must be an array, or a list, or a tuple of len 3"
                assert len(ar)==3, "The rotation must be an array, or a list, or a tuple of len 3"
                self._buf[k]=(O,ap,ar)
        elif a==None:
            pass
        else:
            raise ValueError, "plist can be initialized only with tuples, lists or dictionaries"
    
    #Expose DICT API

    def __len__(self):
        return self._buf.__len__()

    def __getitem__(self , x):
        return self._buf[x]
        
    def __setitem__(self,key,val):
        O,p,r =val
        ap=array(p,dtype=float64)
        ar=array(r,dtype=float64)
        assert len(ap)==3, "The position vector must be an array, or a list, or a tuple of len 3"
        assert len(ar)==3, "The rotation must be an array, or a list, or a tuple of len 3"
        self._buf[key]=(O,ap,ar)
    
    def __delitem__(self, key):
        self._buf.__delitem__(key)

    
    def __contains__(self, key):
        return self._buf.__contains__(key)
    
    
    def __repr__(self):
        return repr(self._buf)
        
    
    # Return an iterator so this can be used similar to a list    
    def __iter__(self):
        return six.itervalues(self._buf)
    
    def iteritems(self):
        return self._buf.iteritems()
    
    def iter(self):
        return self._buf.iter()
        
    def clear(self):
        return self._buf.clear()
    
    def items(self):
        return self._buf.items()
    
    def iterkeys(self):
        return self._buf.iterkeys()
         
    def itervalues(self):
        return self._buf.itervalues()
        
    def keys(self):
        return self._buf.keys()
        
    def pop(self,*argv,**argkw):
        return self._buf.pop(*argv,**argkw)
    
    def popitem(self):
        return self._buf.popitem()
    
    def setdefault(self,*argv,**argkw):
        return self._buf.setdefault(*argv,**argkw)
    
    def update(self, *argv, **argkw):
        return self._buf.update(*argv,**argkw)
 
    def values(self):
        return self._buf.values()
     
    def viewitems(self):
        return self._buf.viewitems()
    
    def viewkeys(self):
        return self._buf.viewkeys()
        
    def viewvalues(self):
        return self._buf.viewvalues()
        
    
    def append(self, a):
        '''
        Append an instance in a similar way as they are append in a list.
        
        The key is auto generated
        '''
        
        O,p,r = a
        ap=array(p,dtype=float64)
        ar=array(r,dtype=float64)
        assert len(ap)==3, "The position vector must be an array, or a list, or a tuple of len 3"
        assert len(ar)==3, "The rotation must be an array, or a list, or a tuple of len 3"
        
        if isinstance(O,Surface): pf="S"
        elif isinstance(O,(Component,System)): pf="C"
        else: pf="A" #Auto
    
        #Find an unused key
        k=0
        while pf+str(k) in self._buf: k=k+1
        
        self._buf[pf+str(k)]=(O,ap,ar)
