from numpy.fft import ifft2
from numpy import unwrap

#a=plane_wave(n=(0, 0, 1), l=.633, size=(64*10., 64*10, ), samples=(64, 64))
#figure(); imshow(abs(fftshift(fft2(a.data))));colorbar()


#a=plane_wave(n=(.1, 0, 1), l=.633, size=(64*10., 64*10, ), samples=(64, 64))
#figure(); imshow(abs(fftshift(fft2(a.data))));colorbar()

#a=plane_wave(n=(0, 0, 1), l=.633, size=(64*10., 64*10, ), samples=(64, 64))

si=64
data=zeros((si,si))#,complex)

#test con un slit
data[:,:]=0.
#data[0.5*si/10:9.5*si/10, 4*si/10:5*si/10]=1.
data[0.5*si/10:9.5*si/10, 4.9*si/10:5.1*si/10]=1.
a=Field(data=data,psize=10,l=.633)

#a=plane_wave(n=(0, 0, 1), l=.633, size=(64*10., 64*10, ), samples=(64, 64))
r=array((0, 0.04, 0))
b=a.tilt(r=r)   
fb=fftshift(fft2(b.data))
figure(); imshow(abs(a.data));colorbar()
figure(); imshow(abs(fb));colorbar()
figure(); imshow(angle(b.data));colorbar()
figure(); imshow(abs(b.data));colorbar()


#ro=[0, 100, 200, 300]    
#figure(); imshow(abs(fftshift(fft2(a.data))));colorbar()
#figure(); imshow(abs(a.data));colorbar()
#figure(); imshow(angle(a.data));colorbar()
#for rot in ro:
#    rot=-rot*1e-3
#    r=array((rot,rot, 0))
#    b=a.tilt(r=r)
#    figure(); imshow(abs(fftshift(fft2(b.data))));colorbar()
#    figure(); imshow(abs(b.data));colorbar()    
#    #figure(); imshow(angle(b.data));colorbar()


