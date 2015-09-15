import vigra
from vigra import graphs

filepath = '12003.jpg' 
filepath = "/home/tbeier/Desktop/boundProb217.png_files/view.png"
img = vigra.impex.readImage(filepath).astype('float32')#[1000:2000,0:1000]
#img[img>150] = 255


res = vigra.filters.shockFilter(img,sigma=1.5, rho=10.0, updwindFactorH=1.0, iterations=5)

#for c in range(res.shape[2]):
#    rc = res[:,:,c]
#    rc -=rc.min()
#    rc /=rc.max()
#    rc*255.0
res = res.squeeze()


import numpy as np
import pylab
import matplotlib.cm as cm
import Image

f = pylab.figure()
for n, arr in enumerate([img,res,img-res]):
    arr= arr.squeeze()
    #f.add_subplot(2, 1, n)  # this line outputs images on top of each other
    f.add_subplot(1, 3, n)  # this line outputs images side-by-side
    pylab.imshow(arr,cmap=cm.Greys_r)
pylab.title('( III x) image')
pylab.show()