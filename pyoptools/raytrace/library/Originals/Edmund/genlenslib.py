# Generate lens library from


#from ray_trace.comp_lib.spherical_lens import SphericalLens
from six.moves import configparser as cp

def tofloat(x):
    try:
        retval=float(x.replace(",","."))
    except:
        retval=0
    return retval

f=open("Edmund Optics spherical data Jan 2007.csv","r")


desc=f.readline()

lib=cp.ConfigParser()


for line in f:
        c=line.split("\t")
        ref=c[0]
        desc=c[1]
        lib.add_section(ref) # Add the references as Sections
        diam=float(c[3].replace(",","."))
        nsurf=int(c[5])

        th1=tofloat(c[10])
        th2=tofloat(c[13])
        th3=tofloat(c[16])
        th4=tofloat(c[19])
        th5=tofloat(c[21])

        cur1=tofloat(c[9])
        cur2=tofloat(c[12])
        cur3=tofloat(c[15])
        cur4=tofloat(c[18])
        cur5=tofloat(c[21])
        cur6=tofloat(c[24])

        if cur1!=0:cur1=1./cur1
        if cur2!=0:cur2=1./cur2
        if cur3!=0:cur3=1./cur3
        if cur4!=0:cur4=1./cur4
        if cur5!=0:cur5=1./cur5
        if cur6!=0:cur6=1./cur6

        mat1=c[11]
        mat2=c[14]
        mat3=c[17]
        mat4=c[20]
        mat5=c[23]


        if nsurf==1:
            lib.set(ref,"type","SphericalMirror")
            lib.set(ref,"radius", diam/2)
            lib.set(ref,"description", desc)

        if nsurf==2:
            lib.set(ref,"type","SphericalLens")
            lib.set(ref,"radius", diam/2)
            lib.set(ref,"description", desc)

            lib.set(ref,"thickness",th1)

            lib.set(ref,"curvature_s1",cur1)
            lib.set(ref,"curvature_s2",cur2)

            lib.set(ref,"material",mat1)

        elif nsurf==3:
            lib.set(ref,"type","Doublet")
            lib.set(ref,"radius", diam/2)
            lib.set(ref,"description", desc)

            lib.set(ref,"curvature_s1",cur1)
            lib.set(ref,"curvature_s2",cur2)
            lib.set(ref,"curvature_s3",cur3)

            lib.set(ref,"thickness_l1",th1)
            lib.set(ref,"thickness_l2",th2)

            lib.set(ref,"material_l1",mat1)
            lib.set(ref,"material_l2",mat2)

        elif nsurf==4:
            lib.set(ref,"type","Triplet")
            lib.set(ref,"radius", diam/2)
            lib.set(ref,"description", desc)

            lib.set(ref,"curvature_s1",cur1)
            lib.set(ref,"curvature_s2",cur2)
            lib.set(ref,"curvature_s2",cur3)

            lib.set(ref,"thickness_l1",th1)
            lib.set(ref,"thickness_l2",th2)
            lib.set(ref,"thickness_l3",th3)


            lib.set(ref,"material_l1",mat1)
            lib.set(ref,"material_l2",mat2)
            lib.set(ref,"material_l3",mat3)

        elif nsurf==6:
            lib.set(ref,"type","AchromatPair")
            lib.set(ref,"radius", diam/2)
            lib.set(ref,"description", desc)

            lib.set(ref,"curvature_s1",cur1)
            lib.set(ref,"curvature_s2",cur2)
            lib.set(ref,"curvature_s3",cur3)
            lib.set(ref,"curvature_s4",cur4)
            lib.set(ref,"curvature_s5",cur5)
            lib.set(ref,"curvature_s6",cur6)

            lib.set(ref,"thickness_l1",th1)
            lib.set(ref,"thickness_l2",th2)
            lib.set(ref,"thickness_l3",th3)
            lib.set(ref,"thickness_l4",th4)
            lib.set(ref,"thickness_l5",th5)

            lib.set(ref,"material_l1",mat1)
            lib.set(ref,"material_l2",mat2)
            lib.set(ref,"material_l3",mat3)
            lib.set(ref,"material_l4",mat4)
            lib.set(ref,"material_l5",mat5)




lib.write(file("SphOptics.cmp","w"))
