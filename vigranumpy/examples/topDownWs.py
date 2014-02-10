import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys
import matplotlib
import pylab as plt
import math
from collections import deque as Queue


print "get input"
f       = '100075.jpg'
#f       = '69015.jpg'
#f       = '12003.jpg'
sigma   = 3.0
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
img     = vigra.filters.gaussianSmoothing(img,0.8)
imgLab  = vigra.colors.transform_RGB2Lab(img)
gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
gradmag = numpy.squeeze(gradmag)

shape=gradmag.shape

seeds       = numpy.zeros(shape,dtype=numpy.uint32)
labels       = numpy.zeros(shape,dtype=numpy.uint32)
numIter  = 10000


def randPoint(shape):


    x=numpy.random.random_integers(low=0,high=shape[0]-1)
    y=numpy.random.random_integers(low=0,high=shape[1]-1)


    xx = x
    yy = y
    while (xx==x or yy==y ):
            xx=numpy.random.random_integers(low=0,high=shape[0]-1)
            yy=numpy.random.random_integers(low=0,high=shape[1]-1)

    return x,y,xx,yy

def intersectLabels(a,b):
    #a=numpy.array(a)
    #b=numpy.array(b)
    r=a.copy()
    ml  = numpy.max([a.max()+1,b.max()+1])
    b*=ml
    r+=a
    r+=b
    #r = a + b*ml
    r = vigra.analysis.labelImage(r)    
    return r


newColor = 3
print "uniqueAFTER FIRST",numpy.unique(labels)


crackedEdges = numpy.array(vigra.analysis.regionImageToCrackEdgeImage(labels,edgeLabel=0))
whereNoEdge    =  numpy.where(crackedEdges!=0)
crackedEdges[whereNoEdge]=0
crackedEdgesSum = crackedEdges.astype(numpy.float32)

for x in range(numIter):

    x1,y1,x2,y2=randPoint(shape)
    seeds[:]=0
    seeds[x1,y1]=1
    seeds[x2,y2]=2

    x1,y1,x2,y2=randPoint(shape)
    seeds[x1,y1]=3
    seeds[x2,y2]=4

    x1,y1,x2,y2=randPoint(shape)
    seeds[x1,y1]=5
    seeds[x2,y2]=6

    sigma = float(numpy.random.rand())*3.0+1.0
    print "simga",sigma
    gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
    gradmag = numpy.squeeze(gradmag)
    labels ,nSeg = vigra.analysis.watersheds(gradmag,seeds=seeds)

    crackedEdges = numpy.array(vigra.analysis.regionImageToCrackEdgeImage(labels,edgeLabel=0))
    whereEdge    =  numpy.where(crackedEdges==0)
    crackedEdgesSum[whereEdge]+=1.0
    print x
    if (x+1) % 500 == 0:
        
        vigra.imshow(numpy.swapaxes(numpy.log(crackedEdgesSum+1),0,1))
        vigra.show()