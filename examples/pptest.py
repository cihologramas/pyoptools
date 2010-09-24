import pp
from time import *
import numpy 

d=numpy.zeros((2*4096,2*4096),dtype=numpy.complex)
d[1024,1024]=1
def test1(d):
    
    return numpy.fft.fft2(d)


jobserver=pp.Server(secret="secret")
print "Starting pp with", jobserver.get_ncpus(), "workers"



#5um pinhole
ti=time()
#print test1(d)

print "enviando proceso 1",time()-ti
f1=jobserver.submit(test1,(d,),(),("numpy",))
print "enviando proceso 2",time()-ti
f2=jobserver.submit(test1,(d,),(),("numpy",))
print "enviando proceso 3",time()-ti
#ph.propagate_rsi_gpu(z, n=1.,dfield=rfd)
print f1()
print f2()

print "rsi",time()-ti
del jobserver
