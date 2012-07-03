cdef class Picklable:
    '''
    Defines support for pickle.
    
    All cython extensions, must inherit from this class in order to support
    Pickling.
    
    '''

    def __init__(self,*argv):
        self.__pkeys__=[]
        for key in argv:
            self.addkey(key)

    def __reduce__(self):
        '''
        Method needed to be able to pickle a surface.
        '''
        #print "in SYSTEM reduce creating a ", type(self)," class" 
        return(type(self),(),self.__getstate__())
    
    def __getstate__(self):
        "Return the current state"
        state={}

        #Fill the cython classes and subclasses
        for key in self.__pkeys__:
            t=getattr(self,key)
            state[key]=t
        
        if hasattr(self,"__dict__"):
            for k, v in self.__dict__.iteritems():
                state[k]=v
        
        return state
        
    def __setstate__(self,state):
        "Fill the current state"
        for k, v in state.iteritems():
            setattr(self,k,v) 
        

    
    cdef addkey(self, key):
        '''Check if the key is valid, and add a key to an attribute to the state list.
        
        Because of python limitations, the attributes must be public
        '''
        
        assert hasattr(self,key), "Class %s has no attribute %s"%(str(type(self)),key)
        self.__pkeys__.append(key)
        

