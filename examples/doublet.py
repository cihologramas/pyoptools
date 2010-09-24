'''Example of a spherical lens simulation

'''


from pyoptools.all import *
from numpy import pi
# Definition of the ray sources at the origin

# Blue
#r_b= parallel_beam_c(size=(10,10),num_rays=(10,10), wavelength=.470)

r_b=point_source_r(origin=(0.,0.,0.),direction=(0.,0.,0),span=pi/256
                   ,num_rays=10,wavelength=0.470, label="blue")

r_b1=point_source_r(origin=(20.,0.,0.),direction=(0.,0.,0),span=pi/256
                   ,num_rays=10,wavelength=0.470, label="blue1")


# Green
#r_g= parallel_beam_c(size=(10,10),num_rays=(10,10), wavelength=.540)
#r_g= point_source_r(origin=(0.,0.,0.),direction=(0.,0.,0),span=pi/512
#                   ,num_rays=100,wavelength=0.540, label="green")

# Red
#r_r= parallel_beam_c(size=(10,10),num_rays=(10,10), wavelength=.670)
#r_r=point_source_r(origin=(0.,0.,0.),direction=(0.,0.,0),span=pi/512
#                   ,num_rays=10,wavelength=0.670, label="red")


N_BK7=schott['BK7']
SF5=schott['SF5']
#Dobletes 32-327 Edmund Scientific
DB1=Doublet(radius=12.5,
    curvature_s1 =1./61.47,
    curvature_s2 =-1./44.64,
    curvature_s3 =-1./129.94,
    thickness_l1 = 6.,
    thickness_l2 = 2.5,
    material_l1  = N_BK7,
    material_l2  = SF5)


DB2=Doublet(radius=12.5,
    curvature_s1 =1./61.47,
    curvature_s2 =-1./44.64,
    curvature_s3 =-1./129.94,
    thickness_l1 = 6.,
    thickness_l2 = 2.5,
    material_l1  = N_BK7,
    material_l2  = SF5)


#Definition of a detector plane

ccd=CCD()
ccd1=CCD()

# Place de tetectors at the focal planes of the lenses

os=System(complist=[(DB1,(20,0,200),(0,0,0)),
                    (DB2,(0,0,200),(pi,0,0)),
                    (ccd,(20,0,float(400)),(0,0,0)),
                    (ccd1,(0,0,400),(0,0,0))
                    ],n=1)





#Add the ray sources

os.ray_add(r_b1)
os.ray_add(r_b)
#os.ray_add(r_g)
#os.ray_add(r_r)
os.propagate()


spot_diagram(ccd)
#title("Spot Diagram")
#ccd.im_show(size=(512,512))
#title("image")
glPlotFrame(os)

