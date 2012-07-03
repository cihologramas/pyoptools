
from numpy import pi
from pyoptools.all import *

'''Example of a spherical lens simulation

'''
N_BK7=schott["BK7"]
N_BAK4=schott["BAK4"]
N_SF10=schott["SF10"]



# Definition of the ray sources at the origin

# Blue
#r_b= parallel_beam_c(size=(10,10),num_rays=(10,10), wavelength=.470)
r_b=point_source_r(origin=(0.,0.,0.),direction=(0.,0.,0),span=pi/512
                   ,num_rays=50,wavelength=0.470)

r_b=Ray(pos=(0, 0, 0),  dir = (0, 20, 500))

# Green
#r_g= parallel_beam_c(size=(10,10),num_rays=(10,10), wavelength=.540)
r_g=point_source_r(origin=(0.,0.,0.),direction=(0.,0.,0),span=pi/512
                   ,num_rays=50,wavelength=0.540)

# Red
#r_r= parallel_beam_c(size=(10,10),num_rays=(10,10), wavelength=.670)
r_r=point_source_r(origin=(0.,0.,0.),direction=(0.,0.,0),span=pi/512
                   ,num_rays=500,wavelength=0.670)


# Definition
oc=Doublet(radius=25,
    curvature_s1 =1./162.59,
    curvature_s2 =-1./123.82,
    curvature_s3 =-1./402.58,
    thickness_l1 = 9.75,
    thickness_l2 = 3.50,
    material_l1  = N_BAK4,
    material_l2  = N_SF10)


bs=BeamSplitingCube(size=50,material=N_BK7,reflectivity=0.5)

#f,af,pf = oc.paraxial_constants(wavelength=.58929)

#print "Paraxial constants @589.29 nm"
#print "Effective Focal Length ",f
#print "Anterior Focal Length  ",af
#print "Posterior Focal Length ",pf


#Definition of a detector plane

ccd=CCD()


# Place de tetectors at the focal planes of the lenses

os=System(complist=[(oc,(0,0,500),(.50,0,0)),
                    (ccd,(0,0,float(990)),(0,0,0)),
                    (bs,(0,0,float(750)),(0,0,0)),
                    ],n=1)

#Add the ray sources
#os.ray_add(r_b)
#os.ray_add(r_g)
os.ray_add(r_r)
os.propagate()


#ccd.spot_diagram(title="Spot Diagram",color=True)
#spot_diagram(ccd)
#ccd.im_show(size=(512,512),title="Image",color=True)

#glPlotFrame(os)

