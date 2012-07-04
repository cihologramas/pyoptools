'''Example of a spherical lens simulation

'''
# TODO: This should be changed. It should import only the calculus libraries
from pyoptools.all import *

# Definition of the ray sources at the origin

# Blue
#r_b= parallel_beam_c(size=(2,2),num_rays=(10,10), wavelength=.470)
r_b=parallel_beam_p(radius = 2., num_rays = (5, 10), wavelength = 0.470)
#r_b=point_source_r(origin=(0.,0.,0.),direction=(0.,0.,0),span=pi/512
#                   ,num_rays=10,wavelength=0.470)



# Green
#r_g= parallel_beam_c(size=(10,10),num_rays=(10,10), wavelength=.540)
#r_g=point_source_r(origin=(5.,0.,0.),direction=(0.,0.,0),span=pi/512
#                   ,num_rays=10,wavelength=0.540)

# Red
#r_r= parallel_beam_c(size=(10,10),num_rays=(10,10), wavelength=.670)
#r_r=point_source_r(origin=(-5.,0.,0.),direction=(0.,0.,0),span=pi/512
#                   ,num_rays=10,wavelength=0.670)


# Lente asferica 47725 edmund scientific
#asf=AsphericalHO(n=10,Kx=-1.076527, Ky=-1.076527, Ax=1/6.63,Ay=1/6.63,   shape=Circular(radius=15./2.), 
#                #ho_cohef=array(((0, 0, 0, 0, 2.396046E-04, 0, 6.414674E-07, 0, 7.685840E-09, 0, -6.476209E-11),
#                               (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#                               (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#                               (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#        "execfile(\""+abspath(FD.GetPath())+"\")"                       (2.396046E-04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#                               (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#                               (6.414674E-07, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#                               (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#                               (7.685840E-09, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#                               (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#                               (-6.476209E-11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))
#                )  


#pl=Plane(size=(16, 16))
#oc=Component(surflist=[(asf, (0, 0, 0), (0, 0, 0)), 
#                       (pl, (0, 0, 7), (0, 0, 0))
#                       ],
#                        material=1.58913) # Ojo el BAL 35 no esta definido

#Thorlabs C330TME-A

R=2.75
k=-0.6139160
A2=0
A4=5.8891900E-04
A6=-1.7660200E-05
A8=1.0102500E-05
A10=-3.9148700E-06


r2=poly2d((0,0,0,1.,0,1.))
r4=r2*r2
r6=r4*r2
r8=r4*r4
r10=r8*r2

poly=A2*r2+A4*r4+ A6*r6 +A8*r8 +A10*r10

asf2=Aspherical(Kx=k, Ky=k, Ax=1./R,Ay=1./R, shape=Circular(radius=2.5), 
                poly=poly)  






R=-3.1885400
k=-12.6638600

A2=0
A4=1.2458340e-02
A6=-3.7119450e-03
A8=5.1223910e-04
A10=-3.1085780e-05
poly=A2*r2+A4*r4+ A6*r6 +A8*r8 +A10*r10

asf1=Aspherical(Kx=k, Ky=k, Ax=1./R,Ay=1./R, shape=Circular(radius=2.5), 
                poly=poly, reflectivity=.5)  



oc=Component(surflist=[(asf2, (0, 0, 0), (0, 0, 0)), 
                       (asf1, (0, 0, 2.8+.35), (0,0, 0))
                       ],
                        material=1.58913) # Ojo el BAL 35 no esta definido

##########33
ccd=CCD(size=(3,3))



os=System(complist=[(oc,(0,0,20),(0,0,0)),
                    (ccd,(0,0,20+2.8+2.14),(0,0,0)),
                    ],n=1)





#Add the ray sources
#os.ray_add(Ray(pos=(0, 3, 0)))
os.ray_add(r_b)
#os.ray_add(r_g)
#os.ray_add(r_r)
#os.ray_add(Ray())
os.propagate()


#pf=PlotFrame(opsys=os)
#pf=glPlotFrame(os)
#AsyncCall(glFrame, os)
#pf.plot()
#spot_diagram(ccd)
#a=ccd.spot_diagram()
#a.show()
