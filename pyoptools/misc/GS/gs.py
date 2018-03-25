from mako.template import Template
from pyoptools.misc.resources import has_double_support, has_amd_double_support

### ojo, toca solucionar esta importacion en caso de que no exista pypencl

try:
	from pyfft.cl import Plan
	import pyopencl as cl
	import pyopencl.array as cl_array
except:
	pass



from numpy.fft import fft2,ifft2,fftshift,ifftshift
from numpy import angle,exp,pi, complex128, zeros, sqrt,int32, zeros_like,ones
from numpy.random import random


from pylab import imshow,colorbar

KERNEL= \
"""   
    //There are some operations that are not defined in the RV770 GPUs
    // for doubles, so a cast to float is needed
    
    % if double_support:
        #pragma OPENCL EXTENSION cl_khr_fp64: enable
        #define CAST (double)
    % elif amd_double_support:
        #pragma OPENCL EXTENSION cl_amd_fp64: enable
        #define CAST (float)
    % endif
    
    __kernel void norm(__global double2 *data)
    {
      int nWidth = get_global_size(0);
      int nHeight = get_global_size(1);
      int ox=get_global_id(0); // Toma los indices en X
      int oy=get_global_id(1); // Toma los indices en Y
      int i= oy*nWidth+ox;
      
      double norm=sqrt(CAST(data[i].x*data[i].x+data[i].y*data[i].y));
      if (norm>0)
      {
        data[i].x=data[i].x/norm;
        data[i].y=data[i].y/norm;
      }gs
      else
      {
        data[i].x=1;
        data[i].y=0;
      }
    }
  
    __kernel void norm1(__global double2 *data, __global double2 *idata, __global double *error, int cut)
    {
      int nWidth = get_global_size(0);
      int nHeight = get_global_size(1);
      int ox=get_global_id(0); // Toma los indices en X
      int oy=get_global_id(1); // Toma los indices en Y
      int i;
    
      double norm,intdata;
      i= oy*nWidth+ox;
      error[i]=0;
      ///OJO, aca las matrices vienen con fftshift
      if( ((ox<cut)          && (oy<cut)          ) ||
          ((ox>(nWidth-cut)) && (oy<cut)          ) ||
          ((ox<cut)          && (oy>(nHeight-cut))) ||
          ((ox>(nWidth-cut)) && (oy>(nHeight-cut)))  )
      {
          intdata=data[i].x*data[i].x+data[i].y*data[i].y;
          intdata=sqrt((float)intdata);
          error[i]=(intdata-idata[i].x)*(intdata-idata[i].x);
        
          norm=sqrt(CAST(data[i].x*data[i].x+data[i].y*data[i].y));
          if (norm>0)
          {
            data[i].x=(data[i].x/norm)*idata[i].x;
            data[i].y=(data[i].y/norm)*idata[i].x;
          }
          else
          {
            data[i].x=idata[i].x;
            data[i].y=0;
          }
      }
      
      
    }


    __kernel void norm2(__global double2 *data, __global double2 *idata)
    {
      int nWidth = get_global_size(0);
      int nHeight = get_global_size(1);
      int ox=get_global_id(0); // Toma los indices en X
      int oy=get_global_id(1); // Toma los indices en Y
      int i;

      double norm;
      i= oy*nWidth+ox;

      norm=sqrt(CAST(data[i].x*data[i].x+data[i].y*data[i].y));
      if (norm>0)
      {
        data[i].x=(data[i].x/norm)*idata[i].x;
        data[i].y=(data[i].y/norm)*idata[i].x;
      }
      else
      {
        data[i].x=idata[i].x;
        data[i].y=0;
      }
      
      
    }

     """


#TODO: The GS algorithm should also use an maximum error condition to stop
#      Not only the iteration condition
def gs(idata,itera=10, ia=None):
    """Gerchberg-Saxton algorithm to calculate DOEs
    
    Calculates the phase distribution in a object plane to obtain an 
    specific amplitude distribution in the target plane. It uses a 
    FFT to calculate the field propagation.
    The wavefront at the DOE plane is assumed as a plane wave.
    
    **ARGUMENTS:**
	
		========== ======================================================
		idata      numpy array containing the target amplitude distribution 
        itera      Maximum number of iterations
        ia         Illumination amplitude at the hologram plane if not given
                   it is assumed to be a constant amplitude with a value
                   of 1. If given it should be an array with the same shape
                   of idata
		========== ======================================================
    """
    
    if ia==None:
        inpa=ones(idata.shape)
    else:
        inpa=ia
    
    assert idata.shape==inpa.shape, "ia and idata must have the same dimensions"
    
    fdata=fftshift(fft2(ifftshift(idata)))
    e=1000
    ea=1000
    
    for i in range (itera):
        fdata=exp(1.j*angle(fdata))*inpa
        
        rdata=ifftshift(ifft2(fftshift(fdata)))
        e= (abs(rdata)-idata).std()
        if e>ea: 
            break
        ea=e
        rdata=exp(1.j*angle(rdata))*(idata)
        fdata=fftshift(fft2(ifftshift(rdata)))        
    
    fdata=exp(1.j*angle(fdata))
    return fdata*inpa
    

def gs_mod(idata,itera=10,osize=256):
    """Modiffied Gerchberg-Saxton algorithm to calculate DOEs
    
    Calculates the phase distribution in a object plane to obtain an 
    specific amplitude distribution in the target plane. It uses a 
    FFT to calculate the field propagation.
    The wavefront at the DOE plane is assumed as a plane wave.
    This algorithm leaves a window around the image plane to allow the 
    noise to move there. It only optimises the center of the image.
    
    **ARGUMENTS:**
	
		========== ======================================================
		idata      numpy array containing the target amplitude distribution 
        itera      Maximum number of iterations
        osize      Size of the center of the image to be optimized
                   It should be smaller than the image itself.
		========== ======================================================
    """
    M,N=idata.shape
    cut=osize//2
    
    
    zone=zeros_like(idata)
    zone[M/2-cut:M/2+cut,N/2-cut:N/2+cut]=1
    zone=zone.astype(bool)

    mask=exp(2.j*pi*random(idata.shape))
    mask[zone]=0
    
    #~ imshow(abs(mask)),colorbar()
    
    fdata=fftshift(fft2(ifftshift(idata+mask))) #Nota, colocar esta mascara es muy importante, por que si no  no converge tan rapido
    
    e=1000
    ea=1000
    for i in range (itera):
        fdata=exp(1.j*angle(fdata))

        rdata=ifftshift(ifft2(fftshift(fdata)))
        #~ e= (abs(rdata[zone])-idata[zone]).std()
        #~ if e>ea: 
           #~ 
            #~ break
        ea=e
        rdata[zone]=exp(1.j*angle(rdata[zone]))*(idata[zone])        
        fdata=fftshift(fft2(ifftshift(rdata)))   
    fdata=exp(1.j*angle(fdata))
    return fdata


def gs_gpu(idata,itera=100):
    """Gerchberg-Saxton algorithm to calculate DOEs using the GPU
    
    Calculates the phase distribution in a object plane to obtain an 
    specific amplitude distribution in the target plane. It uses a 
    FFT to calculate the field propagation.
    The wavefront at the DOE plane is assumed as a plane wave.
    
    **ARGUMENTS:**
	
		========== ======================================================
		idata      numpy array containing the target amplitude distribution 
        itera      Maximum number of iterations
		========== ======================================================
    """ 
    
    pl=cl.get_platforms()[0]
    devices=pl.get_devices(device_type=cl.device_type.GPU)
    ctx = cl.Context(devices=[devices[0]])
    queue = cl.CommandQueue(ctx)

    plan = Plan(idata.shape, queue=queue,dtype=complex128) #no funciona con "complex128"
    
    src = str(Template(KERNEL).render(
        double_support=all(
            has_double_support(dev) for dev in devices),
        amd_double_support=all(
            has_amd_double_support(dev) for dev in devices)
        ))
    prg = cl.Program(ctx,src).build() 
    
    idata_gpu=cl_array.to_device(queue, ifftshift(idata).astype("complex128"))
    fdata_gpu=cl_array.empty_like(idata_gpu)
    rdata_gpu=cl_array.empty_like(idata_gpu)
    plan.execute(idata_gpu.data,fdata_gpu.data)
    
    e=1000
    ea=1000
    for i in range (itera):
        prg.norm(queue, fdata_gpu.shape, None,fdata_gpu.data)
        plan.execute(fdata_gpu.data,rdata_gpu.data,inverse=True)
        tr=rdata_gpu.get()
        rdata=ifftshift(tr)
        
        
        #TODO: This calculation should be done in the GPU
        e= (abs(rdata)-idata).std()
        if e>ea: 
            break
        ea=e
        
        prg.norm2(queue, rdata_gpu.shape,None,rdata_gpu.data,idata_gpu.data)
        
        plan.execute(rdata_gpu.data,fdata_gpu.data)
    
    fdata=fdata_gpu.get()
    
    #~ prg.norm(queue, fdata_gpu.shape, None,fdata_gpu.data)
    fdata=ifftshift(fdata)
    fdata=exp(1.j*angle(fdata))
    
    #~ fdata=fdata_gpu.get()
    return fdata

def gs_mod_gpu(idata,itera=10,osize=256):
    
    
    cut=osize//2
    
    pl=cl.get_platforms()[0]
    devices=pl.get_devices(device_type=cl.device_type.GPU)
    ctx = cl.Context(devices=[devices[0]])
    queue = cl.CommandQueue(ctx)

    plan = Plan(idata.shape, queue=queue,dtype=complex128) #no funciona con "complex128"
    
    src = str(Template(KERNEL).render(
        double_support=all(
            has_double_support(dev) for dev in devices),
        amd_double_support=all(
            has_amd_double_support(dev) for dev in devices)
        ))
    prg = cl.Program(ctx,src).build() 
    

    idata_gpu=cl_array.to_device(queue, ifftshift(idata).astype("complex128"))
    fdata_gpu=cl_array.empty_like(idata_gpu)
    rdata_gpu=cl_array.empty_like(idata_gpu)
    plan.execute(idata_gpu.data,fdata_gpu.data)
    
    mask=exp(2.j*pi*random(idata.shape))
    mask[512-cut:512+cut,512-cut:512+cut]=0
    
    
    idata_gpu=cl_array.to_device(queue, ifftshift(idata+mask).astype("complex128"))
    fdata_gpu=cl_array.empty_like(idata_gpu)
    rdata_gpu=cl_array.empty_like(idata_gpu)
    error_gpu=cl_array.to_device(ctx, queue, zeros(idata_gpu.shape).astype("double"))
    plan.execute(idata_gpu.data,fdata_gpu.data)
    
    e=1000
    ea=1000
    for i in range (itera):
        prg.norm(queue, fdata_gpu.shape, None,fdata_gpu.data)
        plan.execute(fdata_gpu.data,rdata_gpu.data,inverse=True)
        #~ prg.norm1(queue, rdata_gpu.shape,None,rdata_gpu.data,idata_gpu.data,error_gpu.data, int32(cut))
        norm1=prg.norm1
        norm1.set_scalar_arg_dtypes([None, None, None, int32])
        norm1(queue, rdata_gpu.shape,None,rdata_gpu.data,idata_gpu.data,error_gpu.data, int32(cut))
        
        e= sqrt(cl_array.sum(error_gpu).get())/(2*cut)

        #~ if e>ea: 
           #~ 
            #~ break
        #~ ea=e
        plan.execute(rdata_gpu.data,fdata_gpu.data)
    
    fdata=fdata_gpu.get()
    fdata=ifftshift(fdata)
    fdata=exp(1.j*angle(fdata))
    return fdata
