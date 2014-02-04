import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys
import matplotlib
import pylab




print "get input"
f       = '100075.jpg'
f       = '69015.jpg'
#f       = '12003.jpg'
sigma   = 1.0
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab  = vigra.colors.transform_RGB2Lab(img)
gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
gradmag = numpy.squeeze(gradmag)

print gradmag.shape,gradmag.dtype


weights = gradmag.copy()
weights = numpy.exp(-1.0*weights)

path = vigra.graphs.shortestImagePathAStar(weights, (58,54),(237,105) )
#path = vigra.graphs.shortestImagePathDijkstra(weights, (58,54),(237,105) )
#print path
coords = (path[:,0],path[:,1])


gradmag[coords]=0

visuimg =img.copy()

visuimg[:,:,0]=gradmag[:,:]
visuimg[:,:,1]=gradmag[:,:]
visuimg[:,:,2]=gradmag[:,:]

iR=visuimg[:,:,0]
iG=visuimg[:,:,1]
iB=visuimg[:,:,2]

iR[coords]=visuimg.max()
iG[coords]=0
iB[coords]=0


vigra.imshow(visuimg)
pylab.show()