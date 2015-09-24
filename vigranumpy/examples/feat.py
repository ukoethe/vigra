import vigra
import vigra.graphs as vigraph
import vigra.graphs as graphs
import pylab
import numpy
import numpy as np

import sys
import matplotlib
import pylab as plt
import math
from matplotlib.widgets import Slider, Button, RadioButtons






# parameter:
imPath = ('holyRegion.h5', 'im')   # input image path
labPath = ('segMaskOnly.h5', 'data')   # labeled image path

# load volume
labels = vigra.impex.readHDF5(*labPath).astype(np.uint32)#[:,:,0:20]
volume = vigra.impex.readHDF5(*imPath)#[:,:,0:20]


gridGraph = graphs.gridGraph(labels.shape)
rag = graphs.regionAdjacencyGraph(gridGraph, labels)

featureExtractor = graphs.gridRagFeatureExtractor(rag, labels)
featureExtractor.labels = labels
featureExtractor.graph  = rag

miMa = float(volume.min()),float(volume.max()) 
accFeat = featureExtractor.accumulatedFeatures(volume.astype('float32'),miMa[0],miMa[1])
geoFeat = featureExtractor.geometricFeatures()
topoFeat = featureExtractor.topologicalFeatures()

feat = numpy.concatenate([accFeat, geoFeat, topoFeat],axis=1).astype('float32')



gui = vigra.graphs.TinyEdgeLabelGui(rag=rag, img=volume, edgeLabels=None, labelMode=True)
gui.startGui()
labels =  gui.edgeLabels
labels[labels==0] = 19
labels+=1
labels/=2

whereLabels = numpy.where(labels<=1)[0]

X = feat[whereLabels,:]
Y = labels[whereLabels]
Y = Y[:,None].astype('uint32')
print X.shape,Y.shape


rf = vigra.learning.RandomForest(treeCount=1000)
rf.learnRF(X,Y)

p = rf.predictProbabilities(feat)[:,0]


gui = vigra.graphs.TinyEdgeLabelGui(rag=rag, img=volume, edgeLabels=p, labelMode=False)
gui.startGui()

sys.exit(0)
