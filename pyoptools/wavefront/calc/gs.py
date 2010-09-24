
from pyoptools.all import *
from time import time
from numpy import exp,pi
from numpy.random import random


def GScgh(z,target,reference=1.,iterations=20,error=None):
    '''
    Gerber Saxton Algorithm for Fourier CGH
    
    Function that calculates Fourier holographic optical elements using
    the Gerber Saxton algorithm. The algorithm uses the FFT implementation
    of the Fraunhoffer Transform.
    
    **ARGUMENTS:**
    
        ========== ===================================================
        z          Propagation distance. This is used to calculate the 
                   resolution needed in the CGH
        target     Field used to represent the amplitude distribution 
                   to be obtained in the target plane.
        reference  Field used as reference to reconstruct the hologram. If
                   not given, a unitary amplitude plane wave is used. 
        iterations Maximum number of iterations
        error      Expected error
        ========== ===================================================
    
    
    **RETURN VALUE:**
        (holo,err)
        
        ====  ===========================================================
        holo  Field instance, containing the phase mask obtained from the
              iterative algorithm. The holo.res attribute contains the
              resolution of the calculated hologram for the given 
              propagation distance. The holo.l attribute contains the 
              wavelenght used to calculate the hologram.
        
        err   Final error obtained
        ====  ===========================================================


    '''
    
    
    
    #assert reference.shape==target.shape,\
    #    "The reference field, and the target, must have the same shape"
    
    #assert target.l====reference.l,\
    #    "The wave lenghts for the reference beam, and the target must be equal"
    if type(reference) is float:
        ramp=reference
        rang=0.
    else: #Reference is a field
        ramp=reference.abs()
        rang=reference.angle

    ctarget=target*exp(-2.j*pi*random(target.shape))
    
    holo=target.propagate_fraunhofer(-z)
    ntarget=target.abs()/target.abs().max()
    for n in range(iterations):
            if n!=0: holo=imp.propagate_fraunhofer(-z)
            
            #Keep only the phase in the hologram plane
            holo.data=exp(1.j*holo.angle)*ramp
            #Calculate the new image plane
            imp=holo.propagate_fraunhofer(z)
            err=(ntarget-imp.abs()/imp.abs().max()).std()
            if err!=None and err<error: break
            d=exp(1.j*(imp.angle-rang))
            imp=Field(data=d, psize=imp.psize, l=imp.l)
            imp=imp*target.abs()
    
    #Take into account the reference beam phase        
    holo=holo*exp(-1.j*rang)
    return holo,err


def GScghFr(z,target,reference=1.,iterations=20,error=None):
    '''
    Gerber Saxton Algorithm for Fresnel CGS
    
    Function that calculates a Fresnel holographic optical elements using
    the Gerber Saxton algorithm. The algorithm uses the FFT implementation
    of the Fresnel transform.
    
    .. note::
        A Fresnel CGH can be interpreted as a Fourier CGH, where the
        Fourier transforming lens has been multiplied into the CGH. 
    
    **ARGUMENTS:**
    
        ========== ===================================================
        z          propagation distance
        target     Field used to represent the amplitude distribution 
                   to be obtained in the target plane.
        reference  Field used as reference to reconstruct the hologram. If
                   not given, a unitary amplitude plane wave is used. 
        iterations Maximum number of iterations
        error      Expected error
        ========== ===================================================
    
    
    **RETURN VALUE:**
        (holo,err)
        
        ====  ===========================================================
        holo  Field instance, containing the phase mask obtained from the
              iterative algorithm. The holo.res attribute contains the
              resolution of the calculated hologram for the given 
              propagation distance. The holo.l attribute contains the 
              wavelenght used to calculate the hologram.
        
        err   Final error obtained
        ====  ===========================================================


    '''
    
    
    
    #assert reference.shape==target.shape,\
    #    "The reference field, and the target, must have the same shape"
    
    #assert target.l====reference.l,\
    #    "The wave lenghts for the reference beam, and the target must be equal"
    if type(reference) is float:
        ramp=reference
        rang=0.
    else: #Reference is a field
        ramp=reference.abs()
        rang=reference.angle

    ctarget=target*exp(-2.j*pi*random(target.shape))
    
    holo=target.propagate_fresnel(-z)
    ntarget=target.abs()/target.abs().max()
    for n in range(iterations):
            if n!=0: holo=imp.propagate_fresnel(-z)
            
            #Keep only the phase in the hologram plane
            holo.data=exp(1.j*holo.angle)*ramp
            #Calculate the new image plane
            imp=holo.propagate_fresnel(z)
            err=(ntarget-imp.abs()/imp.abs().max()).std()
            if err!=None and err<error: break
            d=exp(1.j*(imp.angle-rang))
            imp=Field(data=d, psize=imp.psize, l=imp.l)
            imp=imp*target.abs()
            
    #Take into account the reference beam phase        
    holo=holo*exp(-1.j*rang)

    return holo,err


def GScghAE(z,target,reference=1.,iterations=20,error=None):
    '''
    Gerber Saxton Algorithm for to calculate CGS, using the Angular 
    Spectrum propagation method.
    
    Function that calculates a holographic optical elements using
    the Gerber Saxton algorithm.  has been multiplied into the CGH. 
    
    **ARGUMENTS:**
    
        ========== ===================================================
        z          propagation distance
        target     Field used to represent the amplitude distribution 
                   to be obtained in the target plane.
        reference  Field used as reference to reconstruct the hologram. If
                   not given, a unitary amplitude plane wave is used. 
        iterations Maximum number of iterations
        error      Expected error
        ========== ===================================================
    
    
    **RETURN VALUE:**
        (holo,err)
        
        ====  ===========================================================
        holo  Field instance, containing the phase mask obtained from the
              iterative algorithm. The holo.res attribute contains the
              resolution of the calculated hologram for the given 
              propagation distance. The holo.l attribute contains the 
              wavelenght used to calculate the hologram.
        
        err   Final error obtained
        ====  ===========================================================


    '''
    
    
    
    #assert reference.shape==target.shape,\
    #    "The reference field, and the target, must have the same shape"
    
    #assert target.l====reference.l,\
    #    "The wave lenghts for the reference beam, and the target must be equal"
    if type(reference) is float:
        ramp=reference
        rang=0.
    else: #Reference is a field
        ramp=reference.abs()
        rang=reference.angle

    ctarget=target*exp(-2.j*pi*random(target.shape))
    
    holo=target.propagate_ae(-z)
    ntarget=target.abs()/target.abs().max()
    for n in range(iterations):
            if n!=0: holo=imp.propagate_ae(-z)
            
            #Keep only the phase in the hologram plane
            holo.data=exp(1.j*holo.angle)*ramp
            #Calculate the new image plane
            imp=holo.propagate_ae(z)
            err=(ntarget-imp.abs()/imp.abs().max()).std()
            if err!=None and err<error: break
            d=exp(1.j*(imp.angle-rang))
            imp=Field(data=d, psize=imp.psize, l=imp.l)
            imp=imp*target.abs()
            
    #Take into account the reference beam phase        
    holo=holo*exp(-1.j*rang)

    return holo,err
