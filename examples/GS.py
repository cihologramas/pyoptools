"""
Gerber Saxton example
"""


from pyoptools.all import *
from numpy import ones


#Working wavelength
l=.442e-3

# image pixel size
ips=.450e-3         

# Read the Image plane
imp=Field(psize=(ips,ips), l=l,amp_im="oso_anteojos.png")
imp=imp.resize((1024,1024))

nx,ny=imp.shape

#hologram pixel size
hps=.05


# Distance where the condition for the pixel size of the image, and the
# Hologram are met. The input image is assumed square

z=ips*hps*nx/l

print "Distance of propagation :", z

ref=Field(data=ones((nx,ny)),psize=(hps,hps),l=l)
ref=ref.resize((nx,ny))



holo,err=GScgh(z,imp,ref,100,error=0.001)
#holo,err=GScghFr(z,imp,ref,100,error=0.05)
#holo,err=GScghAE(z,imp,ref,100,error=0.05)

#imp=(holo*ref).propagate_fresnel(z)
imp=(holo*ref).propagate_fraunhofer(z)
#imp=(holo*ref).propagate_ae(z)

figure()
imshow(holo.angle,interpolation='nearest');colorbar()
figure()
imshow(imp.intensity(),interpolation='nearest');colorbar()
im1=Field(psize=(ips,ips), l=l,amp_im="oso_anteojos.png")
im1=im1.resize((1024,1024))

holo1=im1.propagate_fraunhofer(-z)

figure()
imshow(holo1.angle,interpolation='nearest');colorbar()
figure()
imshow(holo1.intensity(),interpolation='nearest');colorbar()

figure()
imshow(holo.angle-holo1.angle,interpolation='nearest');colorbar()
