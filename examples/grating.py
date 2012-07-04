from gui.RayTraMods import *

"""Program to test the difrection grating, using diferent configurations"""





g=RPPMask(shape=Rectangular(size=(10,10)), phm=poly2d([0,2*pi*0.04,0,]),M=[-1,-2,-3])
#2*pi*0.04
#g=Plane(shape=Rectangular(size=(25,15)))

oc=Component(surflist=[(g, (0, 0, 0), (0, 0, 0)), 
                       ])
ccd=CCD(size=(1000,50))



os=System(complist=[(oc,(0,0,20),(0,0,0)),
                    (ccd,(0,0,120),(0,0,0)),
                    ],n=1)



#r=parallel_beam_c(origin=(0.,0.,0.),direction=(0.,0.,0.),size=(10.,10.)\
#                      ,num_rays=(10,10), label="")
r=[]

#Test using many wavelenghts
for w in (.4,.45,.5,.53,.6,.65,.7):
    r.append(Ray(wavelength=w))

#Add the ray sources
#os.ray_add(Ray(pos=(0, 3, 0)))
os.ray_add(r)
#os.ray_add(r_g)
#os.ray_add(r_r)
#os.ray_add(Ray())
os.propagate()


#pf=PlotFrame(opsys=os)
pf=glPlotFrame(os)
#AsyncCall(glFrame, os)
#pf.plot()
#spot_diagram(ccd)
#a=ccd.spot_diagram()
#a.show()

for pi,rd in ccd.hit_list:
    #for normal incidence
    #d sin(theta)=m lambda
    #d=m lambda/sin(theta)
    print pi, rd.wavelength
    sintheta=pi[0]/sqrt(100.**2+pi[0]**2)
    print 1.* rd.wavelength /sintheta, sintheta
