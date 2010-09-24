from gui.RayTraMods import *

"""Program to test the difrection grating, using diferent configurations"""




l=.45
g=RPPMask(shape=Rectangular(size=(10,10)), phm=poly2d([0,0,0,(2*pi/(l*1e-3))*(1./(2*100)),0,(2*pi/(l*1e-3))*(1./(2*100))]),M=[-1])
#2*pi*0.04
#g=Plane(shape=Rectangular(size=(25,15)))

oc=Component(surflist=[(g, (0, 0, 0), (0, 0, 0)), 
                       ])
ccd=CCD(size=(50,50))



os=System(complist=[(oc,(0,0,20),(0,0,0)),
                    (ccd,(0,0,120),(0,0,0)),
                    ],n=1)



r=parallel_beam_c(origin=(0.,0.,0.),direction=(0.,0.,0.),size=(10.,10.)\
                      ,num_rays=(5,5),wavelength=l, label="")

os.ray_add(r)
os.propagate()


#pf=glPlotFrame(os)
spot_diagram(ccd)
