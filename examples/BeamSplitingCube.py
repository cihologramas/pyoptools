'''Example of a Beam Splitting Cube Simulation

'''

from pyoptools.all import *


# Definition of the ray sources at the origin

# Blue
r_b= parallel_beam_c(size=(10,10),num_rays=(10,10), wavelength=.470)

N_BK7=schott['BK7']

bs=BeamSplitingCube(size=50,material=N_BK7,reflectivity=0.5)


#Definition of a detector plane

ccd=CCD(size=(20, 20))


# Place de detectors at the focal planes of the lenses

os=System(complist=[(bs,(0,0,50),(0,0,0)),
                    (ccd,(0,0,100),(0,0,0)),
                    ],n=1)


#Add the ray sources
os.ray_add(r_b)

#Propagate the rays
os.propagate()


#spot_diagram(ccd)
#title("spot_diagram")

#pf=PlotFrame(opsys=os)
#glPlotFrame(os)
