'''Example of a spherical lens simulation

'''



# Definition of the ray sources at the origin

# Blue
#r_b= parallel_beam_c(size=(2,2),num_rays=(10,10), wavelength=.470)
r_b=parallel_beam_p(radius = 2.5, num_rays = (12, 16), wavelength = 0.470)
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
r2=array(((0, 0, 1), (0, 0, 0), (1, 0, 0)))
r4=polymul2d(r2, r2)
r6=polymul2d(r4, r2)
r8=polymul2d(r6, r2)
r10=polymul2d(r8, r2)
#resize to be able to sum 
rr2=zeros(r10.shape)
rr4=zeros(r10.shape)
rr6=zeros(r10.shape)
rr8=zeros(r10.shape)

rr2[0:3, 0:3]=r2
rr4[0:5, 0:5]=r4
rr6[0:7, 0:7]=r6
rr8[0:9, 0:9]=r8



A4=2.396046E-04; A6=6.414674E-07; A8=7.685840E-09; A10=-6.476209E-11;

cohef=A4*rr4+ A6*rr6 +A8*rr8 +A10*r10

asf=Aspherical(n=10,Kx=-1.076527, Ky=-1.076527, Ax=1/6.63,Ay=1/6.63,shape=Circular(radius=15./2.), 
                 cohef=cohef
                )  


pl=Plane(shape=Rectangular(size=(7.5, 7.5)))
oc=Component(surflist=[(asf, (0, 0, 0), (0, 0, 0)), 
                       (pl, (0, 0, 7), (0, 0, 0))
                       ],
                        material=1.58913) # Ojo el BAL 35 no esta definido


##########33
ccd=CCD(size=(3,3))



os=System(complist=[(oc,(0,0,20),(0,0,0)),
                    (ccd,(0,0,35),(0,0,0)),
                    ],n=1)





#Add the ray sources
#os.ray_add(Ray(pos=(0, 3, 0)))
os.ray_add(r_b)
#os.ray_add(r_g)
#os.ray_add(r_r)
#os.ray_add(Ray())
os.propagate()


#pf=PlotFrame(opsys=os)

glPlotFrame(os)
spot_diagram(ccd)

#ccd.spot_diagram()
