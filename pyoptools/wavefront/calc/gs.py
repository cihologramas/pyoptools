
from pyoptools.all import *
from time import time
from numpy import exp,pi,angle,ones
from numpy.random import random
from numpy.fft import fft2, ifft2, fftshift, ifftshift

def ffGS(z,target,estimate=None, iterations=20,error=None):
	'''
	Far field Gerchberg - Saxton Algorithm
	
	Calculates the phase distribution in a object plane (for a given 
	amplitude constrain) to obtain an specific amplitude distribution in
	the target plane.
	It uses the Gerchberg - Saxton algorithm for Fraunhoffer propagation. 
	A FFT implementation of the Fraunhoffer Transform is used.
	
	**ARGUMENTS:**
	
		========== ======================================================
		z          Propagation distance. This is used to calculate the
		           resolution needed in the object plane, for a given
		           target resolution.
		target     :class:`Field` instance whose amplitude distribution
		           is used to represent the amplitude constrain to be
		           applied in the target plane. The phase of this field
		           is not used.
		estimate   :class:`Field` instance used as initial estimate for
		           the problem. The amplitude of this field is taken as
		           the reference amplitude and the phase is obtained. The
		           resolution used to define this field must match the
		           value needed to obtain the required target resolution
		           when the FFT-Fraunhoffer transform is used. If the
		           wrong value is given an exception is raised.
		           If not given, a unitary amplitude wave, with random
		           phase and the correct resolution, is used.
		iterations Maximum number of iterations
		error      Expected error
		========== ======================================================
	
		.. note:: target and object must have the same wavelength 
	
	**RETURN VALUE:**
		(holo,err)
		
		====  ==========================================================
		holo  Field instance, containing the reference amplitude 
			  information and the phase obtained from the iterative 
			  algorithm. The holo.res attribute contains the
			  resolution of the calculated hologram for the given 
			  propagation distance. The holo.l attribute contains the 
			  wavelength used to calculate the hologram.
		
		err   Final error obtained
		====  ==========================================================


	'''
	
	
	if estimate==None:
		edata=exp(2.j*pi*random(target.shape))
		sx,sy=target.size 
		dxe=target.l*z/sx
		dye=target.l*z/sy
		estimate=Field(data=edata,psize=(dxe,dye),l=target.l)

   
	assert estimate.shape==target.shape,\
		"The estimate field, and the target field, must have the same shape"
	
	assert target.l==estimate.l,\
		"The wave lengths for the reference beam, and the target must be equal"
		
	sx,sy=target.size 
	dxe=target.l*z/sx
	dye=target.l*z/sy
	
	dx,dy=estimate.res
	
	assert (dxe==dx) and (dye==dy),\
		"The resolution for the reference beam, and the target must be equal" 
	 
	
	
	
	holo=estimate
	eabs=estimate.abs()
	
	#Normalized Target amplitude
	ntarget=target.abs()/target.abs().max()

	for n in range(iterations):
		
			if n!=0: holo=imp.propagate_fraunhofer(-z)

			#Keep only the phase in the hologram plane
			holo.data=exp(1.j*holo.angle)
			holo=holo*eabs
			
			#Calculate the new image plane
			imp=holo.propagate_fraunhofer(z)
			
			err=(ntarget-imp.abs()/imp.abs().max()).std()
			
			if error!=None and err<error: break
			
			d=exp(1.j*imp.angle)
			imp=Field(data=d, psize=imp.psize, l=imp.l)
			imp=imp*target.abs()
	return holo,err
	
def fftGS(z,target,estimate=None, iterations=20,error=None,flagRand=True):
	'''
	Far field Gerchberg - Saxton Algorithm
	
	Calculates the phase distribution in a object plane (for a given 
	amplitude constrain) to obtain an specific amplitude distribution in
	the target plane.
	It uses the Gerchberg - Saxton algorithm for far-field propagation,
	using a standard FFT. 
	
	
	**ARGUMENTS:**
	
		========== ======================================================
		z          Propagation distance. This is used to calculate the
		           resolution needed in the object plane, for a given
		           target resolution.
		target     :class:`Field` instance whose amplitude distribution
		           is used to represent the amplitude constrain to be
		           applied in the target plane. The phase of this field
		           is not used.
		estimate   :class:`Field` instance used as initial estimate for
		           the problem. The amplitude of this field is taken as
		           the reference amplitude and the phase is obtained. The
		           resolution used to define this field must match the
		           value needed to obtain the required target resolution
		           when the FFT-Fraunhoffer transform is used. If the
		           wrong value is given an exception is raised.
		           If not given, a unitary amplitude wave, with random
		           phase and the correct resolution, is used.
		iterations Maximum number of iterations
		error      Expected error
		========== ======================================================
	
		.. note:: target and object must have the same wavelength 
	
	**RETURN VALUE:**
		(holo,err)
		
		====  ==========================================================
		holo  Field instance, containing the reference amplitude 
			  information and the phase obtained from the iterative 
			  algorithm. The holo.res attribute contains the
			  resolution of the calculated hologram for the given 
			  propagation distance. The holo.l attribute contains the 
			  wavelength used to calculate the hologram.
		
		err   Final error obtained
		====  ==========================================================


	'''
	
	
	if estimate==None:
		if flagRand:
			edata=exp(2.j*pi*random(target.shape))
		else:
			edata=exp(2.j*pi*ones(target.shape))
		sx,sy=target.size 
		dxe=target.l*z/sx
		dye=target.l*z/sy
		estimate=Field(data=edata,psize=(dxe,dye),l=target.l)

   
	assert estimate.shape==target.shape,\
		"The estimate field, and the target field, must have the same shape"
	
	assert target.l==estimate.l,\
		"The wave lengths for the reference beam, and the target must be equal"
		
	sx,sy=target.size 
	dxe=target.l*z/sx
	dye=target.l*z/sy
	
	dx,dy=estimate.res
	
	assert (dxe==dx) and (dye==dy),\
		"The resolution for the reference beam, and the target must be equal" 
	 
	
	
	
	holo=estimate.data
	eabs=estimate.abs()
	
	#Normalized Target amplitude
	ntarget=target.abs()/target.abs().max()

	for n in range(iterations):
		
			if n!=0: holo=fftshift(fft2(ifftshift(imp)))

			#Keep only the phase in the hologram plane
			holo=exp(1.j*angle(holo))
			holo=holo*eabs
			
			#Calculate the new image plane
			imp=ifftshift(ifft2(fftshift(holo)))
			
			err=(ntarget-abs(imp)/abs(imp).max()).std()
			
			if error!=None and err<error: break
			
			d=exp(1.j*angle(imp))
			imp=d*target.abs()
			
	holo=Field(data=holo, psize=(dxe,dye), l=target.l)
	return holo,err


def frGS(z,target,estimate=None, iterations=20,error=None):
	'''
	Fresnel transform Gerchberg - Saxton Algorithm
	
	Calculates the phase distribution in a object plane (for a given 
	amplitude constrain) to obtain an specific amplitude distribution in
	the target plane.
	
	A FFT implementation of the Fresnel Transform is used.
	
	**ARGUMENTS:**
	
		========== ======================================================
		z          Propagation distance. This is used to calculate the
		           resolution needed in the object plane, for a given
		           target resolution.
		target     :class:`Field` instance whose amplitude distribution
		           is used to represent the amplitude constrain to be
		           applied in the target plane. The phase of this field
		           is not used.
		estimate   :class:`Field` instance used as initial estimate for
		           the problem. The amplitude of this field is taken as
		           the reference amplitude and the phase is obtained. The
		           resolution used to define this field must match the
		           value needed to obtain the required target resolution
		           when the FFT-Fresnel transform is used. If the
		           wrong value is given an exception is raised.
		           If not given, a unitary amplitude wave, with random
		           phase and the correct resolution, is used.
		iterations Maximum number of iterations
		error      Expected error
		========== ======================================================
	
		.. note:: target and object must have the same wavelength 
	
	**RETURN VALUE:**
		(holo,err)
		
		====  ===========================================================
		holo  Field instance, containing the reference amplitude 
			  information and the phase obtained from the iterative 
			  algorithm. The holo.res attribute contains the
			  resolution of the calculated hologram for the given 
			  propagation distance. The holo.l attribute contains the 
			  wavelength used to calculate the hologram.
		
		err   Final error obtained
		====  ===========================================================


	'''
	
	
	if estimate==None:
		edata=exp(2.j*pi*random(target.shape))
		sx,sy=target.size 
		dxe=target.l*z/sx
		dye=target.l*z/sy
		estimate=Field(data=edata,psize=(dxe,dye),l=target.l)

   
	assert estimate.shape==target.shape,\
		"The estimate field, and the target field, must have the same shape"
	
	assert target.l==estimate.l,\
		"The wave lengths for the reference beam, and the target must be equal"
		
	sx,sy=target.size 
	dxe=target.l*z/sx
	dye=target.l*z/sy
	
	dx,dy=estimate.res
	
	assert (dxe==dx) and (dye==dy),\
		"The resolution for the reference beam, and the target must be equal" 
	 
	
	
	
	holo=estimate
	eabs=estimate.abs()
	
	#Normalized Target amplitude
	ntarget=target.abs()/target.abs().max()

	for n in range(iterations):
		
			if n!=0: holo=imp.propagate_fresnel(-z)

			#Keep only the phase in the hologram plane
			holo.data=exp(1.j*holo.angle)
			#Apply the amplitude constain
			holo=holo*eabs
			
			#Calculate the new image plane
			imp=holo.propagate_fresnel(z)
			
			err=(ntarget-imp.abs()/imp.abs().max()).std()
			
			if error!=None and err<error: break
			
			d=exp(1.j*imp.angle)
			imp=Field(data=d, psize=imp.psize, l=imp.l)
			#Apply the amplitude constain
			imp=imp*target.abs()
	return holo,err



def asGS(z,target,estimate=None, iterations=20,error=None):
	'''
	Angular spectrum Gerchberg - Saxton Algorithm
	
	Calculates the phase distribution in a object plane (for a given 
	amplitude constrain) to obtain an specific amplitude distribution in
	the target plane.
	It uses the Gerchberg - Saxton algorithm for the angular spectrum 
	propagation. 
	
	
	**ARGUMENTS:**
	
		========== =====================================================
		z          Propagation distance. This is used to calculate the
		           resolution needed in the object plane, for a given
		           target resolution.
		target     :class:`Field` instance whose amplitude distribution
		           is used to represent the amplitude constrain to be
		           applied in the target plane. The phase of this field
		           is not used.
		estimate   :class:`Field` instance used as initial estimate for
		           the problem. The amplitude of this field is taken as
		           the reference amplitude and the phase is obtained. It
		           must have the same resolution as the target field.
		
		           If not given, a unitary amplitude wave, with random
		           phase and the correct resolution, is used.
		iterations Maximum number of iterations
		error      Expected error
		========== =====================================================
			
		.. note:: target and object must have the same wavelength 
	
	**RETURN VALUE:**
		(holo,err)
		
		====  ===========================================================
		holo  Field instance, containing the reference amplitude 
			  information and the phase obtained from the iterative 
			  algorithm. The holo.res attribute contains the
			  resolution of the calculated hologram for the given 
			  propagation distance. The holo.l attribute contains the 
			  wavelength used to calculate the hologram.
		
		err   Final error obtained
		====  ===========================================================


	'''
	
	
	if estimate==None:
		edata=exp(2.j*pi*random(target.shape))
		sx,sy=target.size 
		dxe=target.l*z/sx
		dye=target.l*z/sy
		estimate=Field(data=edata,psize=(dxe,dye),l=target.l)

   
	assert estimate.shape==target.shape,\
		"The estimate field, and the target field, must have the same shape"
	
	assert target.l==estimate.l,\
		"The wave lengths for the reference beam, and the target must be equal"
		

	dxe,dye=target.res
	
	dx,dy=estimate.res
	
	assert (dxe==dx) and (dye==dy),\
		"The resolution for the estimate beam, and the target must be equal" 
	 
	
	
	
	holo=estimate
	eabs=estimate.abs()
	
	#Normalized Target amplitude
	ntarget=target.abs()/target.abs().max()

	for n in range(iterations):
		
			if n!=0: holo=imp.propagate_ae(-z)

			#Keep only the phase in the hologram plane
			holo.data=exp(1.j*holo.angle)
			holo=holo*eabs
			
			#Calculate the new image plane
			imp=holo.propagate_ae(z)
			
			err=(ntarget-imp.abs()/imp.abs().max()).std()
			
			if error!=None and err<error: break
			
			d=exp(1.j*imp.angle)
			imp=Field(data=d, psize=imp.psize, l=imp.l)
			imp=imp*target.abs()
	return holo,err

	
