import vigra
import vigra.graphs as graphs
import pylab
import numpy
import matplotlib

import sklearn.decomposition
# parameter:
filepath = '100075.jpg'   # input image path
sigmaGradMag = 3.0       # sigma Gaussian gradient
superpixelDiameter = 20  # super-pixel size
slicWeight = 100.0        # SLIC color - spatial weight

# load image and convert to LAB
img = vigra.impex.readImage(filepath)[0:800,0:800].squeeze()


gradMag = vigra.filters.gaussianGradientMagnitude(img,7.0).squeeze()
labels, nseg = vigra.analysis.watershedsNew(gradMag)


labels = vigra.analysis.labelImage(labels)
labels = labels.squeeze()




print vigra.__file__

# get 2D grid graph and  edgeMap for grid graph
# from gradMag of interpolated image
gridGraph = graphs.gridGraph(img.shape[0:2])

# get region adjacency graph from super-pixel labels
rag = graphs.regionAdjacencyGraph(gridGraph, labels)



featureExtractor = graphs.gridRagFeatureExtractor(rag, labels)
featureExtractor.labels = labels
featureExtractor.graph  = rag



gradMag = vigra.filters.gaussianGradientMagnitude(img,3.0).squeeze()
miMa = float(gradMag.min()),float(gradMag.max()) 
accFeat = featureExtractor.accumulatedFeatures(gradMag,miMa[0],miMa[1])
geoFeat = featureExtractor.geometricFeatures()
topoFeat = featureExtractor.topologicalFeatures()
meanGrad = accFeat[:,0]


ucm = featureExtractor.ucmTransformFeatures(meanGrad[:,None])[:,0]

rag.showEdgeFeature(img, meanGrad)
vigra.show()


rag.showEdgeFeature(img, ucm)
vigra.show()


sys.exit(0)


cycleFeats = None
for x in range(1000):
    if cycleFeats is None:
        cycleFeats = featureExtractor.cyclePropergationFeatures(meanGrad)[:,0]
    else:
        cycleFeats = featureExtractor.cyclePropergationFeatures(cycleFeats)[:,0]

    if x%100 ==0 or x==0:
        rag.showEdgeFeature(img, cycleFeats)
        vigra.show()



sys.exit(0)

feat = numpy.concatenate([accFeat],axis=1)



dimRed = sklearn.decomposition.PCA(n_components=3)

dimRedFeat = dimRed.fit_transform(feat)

for i in range(dimRedFeat.shape[1]):
    print i
    rag.showEdgeFeature(img, dimRedFeat[:,i])
    vigra.show()
