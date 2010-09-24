import numpy as N
from pyoptools import *
#Calculo de las coordenadas del ovoide para campo medio plano (Coddington)
#yqp: coordenada sobre el plano imagen virtual (Altura de los puntos
#en el patron de Hartmann)
#yq: coordenada sobre la esfera de referencia
#ZM y YM: Coordenadas en el ovoide medio
#Y. Mejia, 2000
#Version Python O. Olarte 2007

def Polar_Array(Rmax=1,DS=0.1):
    """
    Esta funcion crea rho y theta para evaluar funciones en coordenadas cilindricas
    """
    CG= N.mgrid[float(-Rmax):float(Rmax+DS):float(DS),float(-Rmax):float(Rmax+DS):float(DS)]
    r = N.sqrt(CG[0]**2+CG[1]**2)
    th=N.arccos(N.transpose(CG[0]*1./r))
    th=N.where(th<2*N.pi,th,0)
    th=N.where(CG[0]<0,2*N.pi-th,th)
    return r,th

def mask_psrc(yqp):
    #*********************************************************
    #Parametros iniciales
    #*********************************************************
    R=7.8
    D=180.
    d=1./(2./R+1./D)
    zi=R-d

    #*********************************************************
    #Calculo de yq
    #*********************************************************

    a=(D+d)**2+yqp**2
    b=-(D+R)*(D+d)*yqp #el signo negativo para la parte superior
    c=((D+R)**2-R**2)*yqp**2
    yq=(-b-N.sqrt(b**2-a*c))/a
    zq=N.sqrt(R**2-yq**2)   #coordenada correspondiente sobre la esfera de referencia
    h=D+R-zq
    tgthe=yq/h
    the=N.arctan(tgthe)
    tgalf=yq/zq
    alf=N.arctan(tgalf)
    tgphi=(tgalf+tgthe)/(1-tgalf*tgthe)
    phi=N.arctan(tgphi)

    #**********************************************
    #ovoide para el campo medio plano
    #**********************************************
    Lm=(zq-zi)/N.cos(the)
    A=8.*Lm*N.cos(phi)-2.*R*((N.cos(phi))**2+1.)
    B=2.*R**2*(N.cos(phi))-4.*R*Lm*((N.cos(phi))**2+1.)
    C=2.*R**2*Lm*N.cos(phi)
    LM=(-B+N.sqrt(B**2-4.*A*C))/(2.*A)
    YM=yq-LM*N.sin(2.*phi-the)
    ZM=zq-LM*N.cos(2.*phi-the)
    return ZM,YM


BK7=schott["BK7"]
SF5=schott["SF5"]

r,t=Polar_Array(3.2,0.8)
#r,t=Polar_Array(3.2,1.6)

mask=N.where(r<=3.2,1.,0.)    
Z,R=mask_psrc(r)
X=R*N.cos(t)
Y=R*N.sin(t)

# Colocar el ovoide de tal forma que el vertice de la superficie de
# referencie este en el origen
Z=Z-7.8

XX=N.ravel(X)
YY=N.ravel(Y)
ZZ=N.ravel(Z)

MM=N.ravel(mask)
l=len(XX)


fp=[]



# Espejo que simula el ojo

#S1=Spherical(radius=7,curvature=-1./7.8,reflectivity=1)
#oc=Component(surflist=[(S1,(0,0,0),(0,0,0))])#n=1.5)

# Paraboloide para simular ojo aberrado
RD=7.8
S1=TaylorPoly(shape=Circular(radius=5), cohef=[[0,0,1./(2.*RD)],[0,0,0],[1./(2.*RD),0,0]],
                reflectivity=1)
oc=Component(surflist=[(S1,(0,0,0),(pi,0,0))])#n=1.5)


#Dobletes 32-327 Edmund Scientific
DB1=Doublet(radius=12.5,
    curvature_s1 =1./61.47,
    curvature_s2 =-1./44.64,
    curvature_s3 =-1./129.94,
    thickness_l1 = 6.,
    thickness_l2 = 2.5,
    material_l1  = BK7,
    material_l2  = SF5)

#Dobletes 32-917 Edmund Scientific

DB2=Doublet(radius=12.5,
    curvature_s1 =1./124.12,
    curvature_s2 =-1./87.26,
    curvature_s3 =-1./253.1,
    thickness_l1 = 8.5,
    thickness_l2 = 4.,
    material_l1  = BK7,
    material_l2  = SF5)

BS=BeamSplitingCube(size=25,material=BK7,reflectivity=0.5)

ccd = CCD(size=(200,200))

ccd1= CCD(size=(5,5))

ccd2= CCD(size=(10,10))

ccd3= CCD(size=(10,10))

ccd4= CCD(size=(10,10))


stop=Stop(shape=Rectangular(size=(60,60)),ap_shape=Circular(radius=2.5))

os=System(complist=[(DB2,(0,0,190.37+6.25),(pi,0,0)),
                    (DB1,(0,0,248.92),(0,0,0)),
                    (oc,(0,0,0),(0,0,0)),
                    (ccd1,(0,0,190.37+12.5+8.4+25.02+8.38+8.5+95.95),(0,0,0)),
                    (BS, (0,0,190.37+12.5+8.4+25.02/2.),(0,0,0)),
                    (ccd2,(0,0,180),(0,0,0)),
                    #(stop,(0,0,181),(0,0,0))
                    ],n=1)


# Simplified Optical system used to calculate the principal rays
os1=System(complist=[(oc,(0,0,0),(0,0,0)),
                     (ccd,(0,0,180),(0,0,0)),
                     (stop,(0,0,181),(0,0,0))
                     ]
                     ,n=1)



rl=[]



# Encontrar la posicion de la pupila de entrada
# nota:
# El programa mete los rayos reflejados en la lista de rayos secundarios, por lo
# que se producen 2 trazos de rayo de salida.
# Para este caso el primer rayo del segundo trazo de rayo, es el rayo reflejado

r1=Ray(pos=(0,0,180),dir=(0,0,-1))
r2=Ray(pos=(0,0,180),dir=(0,0.01,-1))
os1.ray_add([r1,r2])
os1.propagate()

## Get the final rays to calculate the virtual image point_source_c
## Do not include the 0 intensity rays produced in the reflections
r3=r1.get_final_rays(inc_zeros=False)[0]
r4=r2.get_final_rays(inc_zeros=False)[0]

#Invert ray direction to calculate the virtual image point_source_c

### Posicion aproximada del centro de la pupila de entrada

#print r1,r5

(xp,yp,zp),vr=intersection(r4,r1)
if vr==True:
    print "Pupila de salida real"
else:
    print "pupila de salida virtual"

print xp,yp,zp

for i in range(l):
    #print i
    # Limitar la pupila a un circulo
    if MM[i]==1:
        cord=(XX[i],YY[i],ZZ[i])
        ry=-pi+N.arctan(XX[i]/(ZZ[i]-zp))
        rx=N.arctan(YY[i]/(ZZ[i]-zp))#
    
        arc=N.arcsin(5./N.sqrt(XX[i]**2+YY[i]**2+ZZ[i]**2))/20.
        ###OJO cambiar el error para mejorar la precision de la imagen
        #print "Ojo, disminuir el error"
        cray=chief_ray_search(os1,ccd,cord,(rx,ry,0),w=arc,er=0.0001,maxiter=1000)
        #cray=Ray(pos=cord, dir=array((xp,yp,zp))-cord)
        #fp.append(cray.copy())
        dx,dy,dz=cray.dir
        rx=N.arcsin(-dy)
        ry=N.arctan2(dx,dz)
        #fp.append(point_source_r(origin=cord,direction=(rx,ry,0),span=arc/15,
        #               num_rays=100,wavelength=0.540))
        fp.append(point_source_r(origin=cord,direction=(rx,ry,0),span=arc/15,
                       num_rays=100,wavelength=0.540))
        
        #fp.append(point_source_c(cord,(rx,ry,0),(arc/100,arc/100),(10,10)))


os.clear_ray_list()
os.reset()

for i in fp:
    os.ray_add(i)

os.propagate()
#pf=PlotFrame(opsys=os)
#pf=glPlotFrame(os)
#ccd2.spot_diagram()
#ccd1.im_show(size=(512,512))
#ccd1.spot_diagram()

#im=ccd1.get_image((5000,5000))

#im.save("/home/richi/Paraboloide.jpg")

#ccds2.im_show("ccds2")
#ccds2.spot_diagram("ccds2")
#ccds3.im_show("ccds3")
#ccds3.spot_diagram("ccds3")
#ccds4.im_show("ccds4")
#ccds4.spot_diagram("ccds4")


