from gui.RayTraMods import *
#ps=Plane(shape=Rectangular(size=(25,15)))
#cs=Plane(shape=Circular(radius=25))
cs=Spherical(shape=Circular(radius=15),curvature=-1/50.)

I1, d1,  er=cs.pw_propagate(Ray(dir=(0.7, 0, 1)),1,1.5, rsamples=(200, 200) , isamples=(100, 100), knots=9)
print er
#figure();imshow(I);colorbar()
#figure();imshow(d);colorbar()

#I2, d2=cs.pw_propagate(Ray(dir=(.3, 0, 1)),1.5,1,  rsamples=(200, 200) ,  isamples=(400, 400))
#print d
figure();imshow(I1);colorbar()
figure();imshow(d1);colorbar()

#figure();imshow(I2);colorbar()
#figure();imshow(d2);colorbar()

#figure();imshow(d+d1);colorbar()
#figure();imshow(I);colorbar()
#figure();imshow(I1);colorbar()
#figure();imshow(da);colorbar()
