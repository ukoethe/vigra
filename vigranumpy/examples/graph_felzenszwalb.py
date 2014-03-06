import vigra
import vigra.graphs as vigraph
import numpy

print "get input"
f       = '100075.jpg'
#f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 4.0

img          = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab       = vigra.colors.transform_RGB2Lab(img)
imgLabInterpolated  = vigra.resize(imgLab,[imgLab.shape[0]*2-1,imgLab.shape[1]*2-1 ])
gradmag      = numpy.squeeze(vigra.filters.gaussianGradientMagnitude(imgLabInterpolated,sigma))
labels ,nseg = vigra.analysis.slicSuperpixels(imgLab,10.0,2)
labels       = numpy.squeeze(vigra.analysis.labelImage(labels))


gridGraph = vigraph.gridGraph(img.shape[0:2])
gridGraphEdgeIndicator =  vigraph.edgeFeaturesFromInterpolatedImage(gridGraph,gradmag)



rag = vigraph.regionAdjacencyGraph(gridGraph, labels)
edgeWeights = rag.accumulateEdgeFeatures(gridGraphEdgeIndicator)
labels   = vigraph.felzenszwalbSegmentation(rag, edgeWeights, k=50)


rag.show(img, labels)
vigra.show()
