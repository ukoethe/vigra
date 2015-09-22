import vigra
import vigra.graphs as graphs
import pylab
import numpy
import matplotlib

# parameter:
filepath = '/media/tbeier/data/datasets/hhess/2x2x2nm_normalized/2x2x2nm.0005.tif.png'   # input image path
sigmaGradMag = 3.0       # sigma Gaussian gradient
superpixelDiameter = 20  # super-pixel size
slicWeight = 100.0        # SLIC color - spatial weight

# load image and convert to LAB
img = vigra.impex.readImage(filepath)[0:800,0:800].squeeze()



labels, nseg = vigra.analysis.watershedsNew(vigra.filters.gaussianSmoothing(-1.0*img,2.0))


labels = vigra.analysis.labelImage(labels)
labels = labels.squeeze()





# get 2D grid graph and  edgeMap for grid graph
# from gradMag of interpolated image
gridGraph = graphs.gridGraph(img.shape[0:2])

# get region adjacency graph from super-pixel labels
rag = graphs.regionAdjacencyGraph(gridGraph, labels)



featureExtractor = graphs.gridRagFeatureExtractor(rag, labels)
featureExtractor.labels = labels
featureExtractor.graph  = rag



accFeat = featureExtractor.accumulatedFeatures(img,float(img.min()),float(img.max()))
geoFeat = featureExtractor.geometricFeatures()

a = numpy.ones_like(img)
a[0,0] = 0
for i in range(5):
    print i
    rag.showEdgeFeature(a, geoFeat[:,i])
    vigra.show()