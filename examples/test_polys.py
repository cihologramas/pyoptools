

R=2.75
k=-0.6139160
A2=0
A4=5.8891900E-04
A6=-1.7660200E-05
A8=1.0102500E-05
A10=-3.9148700E-06


r2=poly2d((0,0,0,1,0,1))
r4=r2*r2
r6=r4*r2
r8=r4*r4
r10=r8*r2

poly=A2*r2+A4*r4+ A6*r6 +A8*r8 +A10*r10



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



A4=5.8891900E-04
A6=-1.7660200E-05
A8=1.0102500E-05
A10=-3.9148700E-06

cohef=A4*rr4+ A6*rr6 +A8*rr8 +A10*r10

x,y=meshgrid(linspace(-10,10,256),linspace(-10,10,256))

p1=poly.eval(x,y)
p2=eval_poly(cohef,x,y)

p3,p4=poly.dxdy()
c3,c4=Poly_DyDx(cohef)
