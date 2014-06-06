import vigra
import numpy
import h5py

def mkCircle(shape, center, radius):
    a = numpy.zeros(shape, dtype=numpy.float32)
    a[center[0], center[1]] = 1
    dist = vigra.filters.distanceTransform2D(a)
    dist = dist < radius
    return dist

cylinder = numpy.zeros((100,101,202), dtype=numpy.uint32)
for i in range(50,150):
    print i
    img = mkCircle((100,101), (50,50), 25) 
    cylinder[:,:,i] = img

f = h5py.File("cylinder.h5", 'w')
f.create_dataset("seg", data=cylinder)
f.close()
