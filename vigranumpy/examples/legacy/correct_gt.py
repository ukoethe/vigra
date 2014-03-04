import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys
import matplotlib
import pylab
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import math
import os
import scipy.io as sio

import numpy as np, h5py 


def bImg(labels):
    crackedEdges = vigra.analysis.regionImageToCrackEdgeImage(labels)
    whereNoEdge    =  numpy.where(crackedEdges!=0)
    whereEdge    =  numpy.where(crackedEdges==0)
    crackedEdges[whereNoEdge]=0
    crackedEdges[whereEdge]=1
    return numpy.squeeze(crackedEdges)



gtFolder = '/home/tbeier/datasets/BSR/BSDS500/data/groundTruth/train/'
imageFolder = '/home/tbeier/datasets/BSR/BSDS500/data/images/train/'


imageFilesNumbers  = [ f.split('.')[0] for f in os.listdir(imageFolder) if f.endswith('.jpg')]


sigma = 1.0

for imgNr in imageFilesNumbers[40:50] :

    imageFile = imageFolder+imgNr+'.jpg'
    gtFile    = gtFolder+imgNr+'.mat'

    image=vigra.impex.readImage(imageFile)
    matContents = sio.loadmat(gtFile)


    gradmag = vigra.filters.gaussianGradientMagnitude(image,sigma)

    #print matContents['groundTruth']

    ngt = len(matContents['groundTruth'][0])
    print "ngts",ngt

    for gti in range(ngt):
        gt =  matContents['groundTruth'][0][gti][0]['Segmentation'][0]
        gt = numpy.swapaxes(gt,0,1)
        gt = gt.astype(numpy.uint32)
        print gt.min(),gt.max()

        print "shrink"
        seg,nseg =  vigra.analysis.watershedsReoptimization(gt,gradmag,10)
        seg=numpy.squeeze(seg)
        edgeGt  = bImg(gt)
        edgeSeg = bImg(seg)
        edgeGt*=1
        edgeSeg*=2

        mixed=edgeGt.copy()
        mixed+=edgeSeg


        whereAB = numpy.where(mixed==3) # both on (GREEN)
        whereA  = numpy.where(mixed==1) # gt on   (BLUE)
        whereB  = numpy.where(mixed==2) # seg on  (RED)

        colorAB = [0,255,0]
        colorA  = [0,0,255]
        colorB  = [255,0,0]

        imgBig = vigra.resize(image,mixed.shape)
        for c in range(3):
            imgBigC = imgBig[:,:,c]
            imgBigC[whereAB]=colorAB[c]
            imgBigC[whereA]=colorA[c]
            imgBigC[whereB]=colorB[c]

        print seg.shape
        vigra.segShow(image,gt)
        vigra.segShow(image,vigra.taggedView(seg,'xy'),[0,1,0])

        vigra.imshow(imgBig)
        vigra.show()
        break

        