from gui.RayTraMods import *
from time import time
from numpy import *
#import pylab as pl

## Programa para calcular la integral de rayleigh sommerfeld, para plano objeto 
## y plano imagen con diferentes muestreos


import pyopencl as cl
pl=cl.get_platforms()[0]
dev=pl.get_devices(cl.device_type.GPU)
print dev


## Inicio del programa

prg_src = """
__kernel void rsk(
                  __global float *x,
                  __global float *y,
                  __global float *ur,
                  __global float *ui,
                  unsigned h,
                  unsigned w,
                  float Z,
                  float K,
                  float dx,
                  float dy,
                  __global float *resr,
                  __global float *resi)
{
    int nWidth = get_global_size(0);
    int nHeigh = get_global_size(1);
    
    int ox=get_global_id(0); // Toma los indices en X
    int oy=get_global_id(1); // Toma los indices en Y
    int oid= oy*nWidth+ox;
    
    float X=(float)(ox-nWidth/2)*dx;
    float Y=(float)(oy-nWidth/2)*dy;
    float pi=3.14159265358979323846264338327950288;
    resr[oid]=0;
    resi[oid]=0;
    float Z2=Z/2.;
    for(unsigned i=0;i<h*w;i++)
    {
        float wi=2.;
        float wj=2.;
        
        int ci=i%w;
        int cj=i/w;
        
        //Hay que verificar si los i y los j estan correctos y si los 
        //coheficientes estan al derecho
        if ((ci==0)||(ci==(w-1)))wi=1;
        else if(ci%2!=0) wi=4;
        
        if ((cj==0)||(cj==(h-1)))wj=1;
        else if(cj%2!=0) wj=4;
        
        float w=wi*wj/9.;
         
        
        
        float Xx=X+x[i];
        float Yy=Y+y[i];
        float R2=pow(Xx,2.f)+pow(Yy,2.f)+pow(Z,2.f);
        float R= pow(R2,0.5f);
        float R3=R2*R;
        float sinkr=sin(R*K);
        float coskr=cos(R*K);
        float piR2=pi*R2;
        float piR3=pi*R3;
        
        //real 1/2*k*z*sin(R*k)/(piR2) + 1/2*z*cos(R*k)/(piR3)
        //imag -1/2*k*z*cos(R*k)/(piR2) + 1/2*z*sin(R*k)/(piR3)
        float re= K*Z2*sinkr/(piR2) + Z2*coskr/(piR3);
        float im= -K*Z2*coskr/(piR2) + Z2*sinkr/(piR3);
        resr[oid] +=  w*(re*ur[i]-im*ui[i]);
        resi[oid] +=  w*(re*ui[i]+im*ur[i]);
    }
}
"""


ctx0 = cl.Context(devices=[dev[0]])#dev_type=cl.device_type.GPU)
#ctx1 = cl.Context(devices=[dev[2]])#dev_type=cl.device_type.GPU)

#self.ctx = cl.Context()

queue0 = cl.CommandQueue(ctx0)
#queue1 = cl.CommandQueue(ctx1)


prg0 = cl.Program(ctx0,prg_src ).build()
#prg1 = cl.Program(ctx1,prg_src ).build()

mf = cl.mem_flags

#???7.5um pinhole
phd=zeros((20, 20))
ph_o=Field(data=phd, psize=2.5e-3, l=.442e-3)
X,Y= ph_o.field_sample_coord
ph_o.data=exp((X**2+Y**2)/.005) 

#Resulting field


msize=1024
#rfd=zeros((msize, msize))
#ph_d=Field(data=rfd, psize=5e-3, l=.442e-3)


Xo,Yo=ph_o.field_sample_coord
Xo=Xo.astype(float32)
Yo=Yo.astype(float32)
Uor=ph_o.data.real.astype(float32)
Uoi=ph_o.data.imag.astype(float32)

phor0=empty((msize,msize),dtype=float32)
phoi0=empty((msize,msize),dtype=float32)

#phor1=empty((msize,msize),dtype=float32)
#phoi1=empty((msize,msize),dtype=float32)


#Crear los buferes de entrada a partir de los arreglos de numpy

#kz=array((k,x,y,z,dx,dy,nx,ny),dtype=float32)
Xob0 = cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Xo)
Yob0 = cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Yo)

Uobr0= cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Uor)

Uobi0= cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Uoi)

#Uobi0= cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Uoi)

#Xob1=cl.Buffer(ctx1, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Xo)
#Yob1=cl.Buffer(ctx1, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Yo)


ti=time()
phbr0=cl.Buffer(ctx0, mf.WRITE_ONLY,phor0.nbytes)
phbi0=cl.Buffer(ctx0, mf.WRITE_ONLY,phoi0.nbytes)

#phbr1=cl.Buffer(ctx1, mf.WRITE_ONLY,phor1.nbytes)
#phbi1=cl.Buffer(ctx1, mf.WRITE_ONLY,phoi1.nbytes)


x=0;y=0;Z=10;
K=2*pi/(.633e-3)
dx=2.5e-3
dy=2.5e-3
prg0.rsk(queue0, (msize,msize),
            Xob0,
            Yob0,
            Uobr0,
            Uobi0,
            uint32(Xo.shape[0]),
            uint32(Xo.shape[1]),
            float32(Z),
            float32(K),
            float32(dx),
            float32(dy), 
            phbr0,
            phbi0,
            #local_size=(16,16)
            )

##prg1.rsk(queue1, (msize,msize),
##            Xob1,
##            Yob1,
##            uint32(ps*ps),
##            float32(Z),
##            float32(K),
##            float32(dx),
##            float32(dy), 
##            phbr1,
##            phbi1,
##            #local_size=(1,1)
##            )

#****************^ Global Size

cl.enqueue_read_buffer(queue0, phbr0, phor0)
#cl.enqueue_read_buffer(queue1, phbr1, phor1)
cl.enqueue_read_buffer(queue0, phbi0, phoi0).wait()
#cl.enqueue_read_buffer(queue1, phbi1, phoi1).wait()

ph0=phor0+1.j*phoi0
print time()-ti
figure()
imshow(abs(ph0))
colorbar()
figure()
imshow(angle(ph0))



