'''Example of a spherical lens simulation
'''



# Definition of the ray sources at the origin
r_b= parallel_beam_c(size=(.1,.1),num_rays=(10,10), wavelength=.470)


L1=Edmund.get("31859")

#Get the paraxial constants for the lens
print L1.paraxial_constants()

ccd=CCD()

os=System(complist=[(L1,(0,0,50),(0,0,0)),
                    (ccd,(0,0,float(54.5)),(0,0,0)),
                    ],n=1)





#Add the ray sources
os.ray_add(r_b)

#Propagate
os.propagate()

#Show 3D ray trace
glPlotFrame(os)

