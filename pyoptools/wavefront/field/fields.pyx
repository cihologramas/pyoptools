from pyoptools.wavefront.field import Field
from numpy import indices,array,exp,sqrt,dot,pi


def plane_wave(n=(0, 0, 1.),l=.633, size=(10000, 10000), samples=(256, 256), a=1.,  ph=0.):
    """
    Function that returns a plane wave (Field object) that describes a plane wave
    
    ** ARGUMENTS **

    ======= ============================================================
    n       Tuple (nx,ny,nz) indicating a vector normal to the plane wave
    l       Wavelength of the plane wave (usually given in microns)
    size    Tuple (sx,sy) indicating the physical size of the window used
            where the field is defined. It must have the same units as l
    samples Tuple (nx,ny) indicating the number of samples to be used when
            creating the plane wave
    a       Amplitude of the plane wave
    ph      Phase of the plane wave at the origin (center of the sampling)
    ======= ============================================================

    """
    nx, ny=samples
    X,Y=indices((nx,ny))
    sx, sy=size
    dx,dy=float(sx)/nx, float(sy)/ny
    ux=(X-nx/2)*dx; uy=(Y-ny/2)*dy
    n=array(n)
    kx, ky, kz=n/(sqrt(dot(n, n)))*2.*pi/l
    
    f=a*exp(1.j*(kx*ux+ky*uy+ph))
    
    return Field(data=f, psize=(dx, dy), l=l)
    
def spherical_wave(o=(0,0,-100),l=.633,size=(10000,10000),samples=(256,256),a=1.,ph=0.):
    '''
    Function that returns a spherical wave (Field instance)
    
    The spherical wave returned, is evaluated in the plane Z=0, assuming
    the source point at the coordinates X,Y,Z, given by ''o''. If o is
    a number, and not a vector, the source point location is the point
    (0,0,o)
    
    ** Arguments: **
    
    ======= ============================================================
    o       Location of the point source, assuming that the observation
            plane is given by Z=0
    l       Wave length of the spherical wave, usually given in microns
    size    Tuple (sx,sy) indicating the physical size of the window used
            where the field is defined. It must have the same units as l
    samples Tuple (nx,ny) indicating the number of samples to be used when
            creating the plane wave
    a       Amplitude of the plane wave
    ph      Phase of the spherical wave at the origin (center of the sampling)
    ======= ============================================================
    
    .. note::
        If z is positive, the spherical wave generated will be convergent
        when the wave is propagating in the positive direction of the Z
        axis. If z is negative, the wave will be divergent.
    '''
    
    nx, ny=samples
    X,Y=indices((nx,ny))
    sx, sy=size
    dx,dy=float(sx)/nx, float(sy)/ny
    ux=(X-nx/2)*dx; uy=(Y-ny/2)*dy
    try:
        x,y,z=o
    except:
      x,y,z=0.,0.,o
        
    sp= sqrt((ux-x)**2+(uy-y)**2+z**2)
    sp=-(sp-sp.min())*z/abs(z)
    
    f=a*exp(2.j*pi*sp/l+ ph)
    
    return Field(data=f, psize=(dx, dy), l=l)
    

    
    
    
    
    
