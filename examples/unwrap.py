#import numpy.ma as ma
from numpy import *
from numpy.ma import MaskedArray
#from gui.RayTraMods import *
#from wavefront.field import Field
from misc.cmisc import unwrap as unwrapc
from pylab import imread
from gui.plotutils import  *
im=imread ("/media/Portable1/Proyectos/wxRayTrace_examples/REAL.PNG")
im=imread ("/media/Portable1/Proyectos/wxRayTrace_examples/Iteraciona.png")
im=imread ("/media/Portable1/Proyectos/wxRayTrace_examples/fase_holo.png")
#X,Y=indices(im[:,:,0].shape)
X,Y=indices(im.shape)
mask=((X-113)**2+(Y-120)**2)>100**2
#mim=MaskedArray(im[:,:,0], dtype=double)#,mask=mask)
#wim=unwrapc(mim, in_p=(100, 100),  uv=1., nn=4)
wim=unwrapc(im.astype(double), in_p=(100, 100),  uv=1., nn=1)
figure()
imshow(wim)
#figure()
#imshow(im)



#X,Y=meshgrid(linspace(-1,1,1000),linspace(-1,1,1000))
#mask=(X**2+Y**2)>1**2
#
#e=array((.3, 0, 1))
#e=e/(sqrt(dot(e, e)))
#l=.632e-3
#k=(2*pi/l)*e
#
##a=MaskedArray(exp(1.j*(k[0]*X+k[1]*Y)), mask=mask)
#a1=MaskedArray(exp(1.j*2*pi/l*sqrt(-X**2-Y**2+10**2)), mask=mask)
#
#
##a2=exp(1.j*(-k[0]*X+k[1]*Y))
#
###a=exp(1.j*pi*(X**2)/10.)
###mask=where ((X-5)**2+(Y-5)**2<9, False, True)
#b=Field(data=a1, psize=X[0, 1]-X[0, 0], l=l)
#
#
#imshow(b.angle)
#figure()
#c=b.phase
#imshow(c)
#colorbar()
#
#
#rl=b.rayrep(20, 20)
#
#
#os=System(complist=[ ],n=1)
#os.ray_add(rl)
#os.propagate()
#pf=PlotFrame(opsys=os)
