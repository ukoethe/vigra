import vigra
from vigra import graphs

filepath = '12003.jpg'
img = vigra.impex.readImage(filepath).astype('float32')[:,:,0]

res = vigra.filters.shockFilter(img,sigma=1.5, rho=10.0, updwindFactorH=1.0, iterations=5)
res = res.squeeze()

import numpy as np
import pylab
import matplotlib.cm as cm

f = pylab.figure()
for n, arr in enumerate([img,res]):
    arr= arr.squeeze().T
    #f.add_subplot(2, 1, n)  # this line outputs images on top of each other
    f.add_subplot(1, 2, n+1)  # this line outputs images side-by-side
    pylab.imshow(arr,cmap=cm.Greys_r)
pylab.title('( III x) image')
pylab.show()
