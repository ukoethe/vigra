import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys

f       = '100075.jpg'
f       = '69015.jpg'
#f       = '12003.jpg'
sigma   = 3.0

print "prepare inpueeet"
img                 = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab              = vigra.colors.transform_RGB2Lab(img)
imgLabInterpolated  = vigra.resize(imgLab,[imgLab.shape[0]*2-1,imgLab.shape[1]*2-1 ])
gradmagInterpolated = vigra.filters.gaussianGradientMagnitude(imgLabInterpolated,sigma)
labels ,nseg        = vigra.analysis.slicSuperpixels(imgLab,10.0,5)
labels              = vigra.analysis.labelImage(labels)


graph0,graph1 = vigraph.gridRegionAdjacencyGraph(labels=labels,ignoreLabel=None)




# get grid graph and edge weights
graph0EdgeWeights =  vigraph.edgeFeaturesFromInterpolatedImage(graph0,gradmagInterpolated)
graph1EdgeWeights = graph1.accumulateEdgeFeatures(graph0EdgeWeights,acc='mean')
graph1NodeFeatures = graph1.accumulateNodeFeatures(img,acc='mean')
graph1Labels = vigraph.felzenszwalbSegmentation(graph1,graph1EdgeWeights,k=1)


graph2       = vigraph.regionAdjacencyGraph(graph=graph1,labels=graph1Labels,ignoreLabel=None)
graph2EdgeWeights = graph2.accumulateEdgeFeatures(graph1EdgeWeights,acc='mean')
graph2Labels = vigraph.felzenszwalbSegmentation(graph2,graph2EdgeWeights,k=5)

graph3       = vigraph.regionAdjacencyGraph(graph=graph2,labels=graph2Labels,ignoreLabel=None)









graph2.show(img=img)
vigra.show()

