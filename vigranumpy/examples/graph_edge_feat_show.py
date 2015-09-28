import vigra
import vigra.graphs as graphs
import pylab
import numpy
import matplotlib

import sklearn.decomposition
# parameter:
# 
# 
if False:
    filepath = "/media/tbeier/data/datasets/hhess/2x2x2nm/2x2x2nm.0002.tif"
    sigmaGradMag = 3.0       # sigma Gaussian gradient"
    superpixelDiameter = 20  # super-pixel size
    slicWeight = 100.0        # SLIC color - spatial weight

    # load image and convert to LAB
    img = vigra.impex.readImage(filepath)[0:1000,0:1000].squeeze().astype('float32')
    # make overseg
    emap = vigra.filters.hessianOfGaussianEigenvalues(img,5.6)[:,:,0].squeeze()
    emap = vigra.filters.gaussianSmoothing(emap, 4.0)
    labels, nseg = vigra.analysis.watershedsNew(emap)
else:
    filepath = "12003.jpg"
    sigmaGradMag = 3.0       # sigma Gaussian gradient"
    superpixelDiameter = 20  # super-pixel size
    slicWeight = 100.0        # SLIC color - spatial weight

    # load image and convert to LAB
    img = vigra.impex.readImage(filepath)[0:1000,0:1000].squeeze().astype('float32')
    # make overseg
    emap = vigra.filters.gaussianGradientMagnitude(img,2.0)[:,:,0].squeeze()
    emap = vigra.filters.gaussianSmoothing(emap, 2.0)
    labels, nseg = vigra.analysis.watershedsNew(emap)

    img = img[:,:,0]

labels = vigra.analysis.labelImage(labels)
labels = labels.squeeze()

# make rag
gridGraph = graphs.gridGraph(img.shape[0:2])
rag = graphs.regionAdjacencyGraph(gridGraph, labels)





def smoothedGradMag(img, sigma, simgaOuter=3.0):
    gradMag = vigra.filters.gaussianGradientMagnitude(img, sigma=sigma).squeeze()
    gradMag = vigra.gaussianSmoothing(gradMag,simgaOuter)
    return gradMag


def stEv(img, sigma, sigmaOuter=None):
    if sigmaOuter is None:
        sigmaOuter = sigma * 2.0
    res = vigra.filters.structureTensorEigenvalues(img, sigma, sigmaOuter)[:,:,0]
    return res

def hessianEv(img, sigma, sigmaOuter=None):
    if sigmaOuter is None:
        sigmaOuter = sigma * 2.0
    res = vigra.filters.hessianOfGaussianEigenvalues(img, sigma)[:,:,0]
    return res


def hessianEv2(img, sigma, sigmaOuter=None):
    if sigmaOuter is None:
        sigmaOuter = sigma * 2.0
    res = vigra.filters.hessianOfGaussian(img, sigma)
    res = vigra.filters.gaussianSmoothing(res, sigmaOuter)
    res = vigra.filters.tensorEigenvalues(res)[:,:,0]
    return res
    

def makeFeat(rag, raw,labels, show=True):

    featureExtractor = graphs.gridRagFeatureExtractor(rag, labels)
    featureExtractor.labels = labels
    featureExtractor.graph  = rag



    geoFeat = featureExtractor.geometricFeatures()
    topoFeat = featureExtractor.topologicalFeatures()

    features = [geoFeat, topoFeat]

    # ward facs
    wardness = numpy.array([0.0, 0.1, 0.15, 0.25, 0.5, 1.0], dtype='float32')

    accFeat = featureExtractor.accumulatedFeatures(raw)
    features.append(accFeat)


    for s in [2.0, 3.0, 4.0]:
        res = vigra.filters.gaussianGradientMagnitude(raw, s)
        accFeat = featureExtractor.accumulatedFeatures(res)
        features.append(accFeat)

    for s in [2.0, 3.0, 4.0]:
        res = vigra.filters.laplacianOfGaussian(raw, s)
        accFeat = featureExtractor.accumulatedFeatures(res)
        features.append(accFeat)



    for s in [2.0, 3.0, 4.0]:
        res = vigra.gaussianSmoothing(raw, s)
        accFeat = featureExtractor.accumulatedFeatures(res)
        features.append(accFeat)
        

    #  st ev
    for s in [2.0, 3.0, 4.0]:
        res = stEv(raw, s)
        accFeat = featureExtractor.accumulatedFeatures(res)
        mean = accFeat[:,0]
        ucm = featureExtractor.ucmTransformFeatures(mean[:,None],wardness)
        features.extend([accFeat,ucm])


        if False:
            for x in range(ucm.shape[1]):
                rag.showEdgeFeature(raw, ucm[:,x])
                vigra.show()


    #  hessian ev
    for s in [1.0,  3.0,  4.0, 6.0, 8.0]:


        img = hessianEv2(raw, s, 2.0)
        #vigra.imshow(img)
        #vigra.show()
        accFeat = featureExtractor.accumulatedFeatures(img)
        mean = accFeat[:,0]
        ucm = featureExtractor.ucmTransformFeatures(mean[:,None],wardness)

        features.extend([accFeat,ucm])
        if False:
            for x in range(ucm.shape[1]):
                print "x",x
                rag.showEdgeFeature(raw, ucm[:,x])
                vigra.show()
                #break
                # 


    features = numpy.concatenate(features,axis=1)
    return features

features = makeFeat(rag, raw=img, labels=labels)

print features.shape




gui = vigra.graphs.TinyEdgeLabelGui(rag=rag, img=img.squeeze(), edgeLabels=None, labelMode=True)
gui.startGui()
Y =  gui.edgeLabels
print Y.min(), Y.max()
whereLabels = numpy.where(Y!=0)[0]

Y = (Y[whereLabels].astype('int32')+1)/2
Y = Y.astype('uint32')
X = features[whereLabels]


rf = vigra.learning.RandomForest(treeCount=10000)
oob = rf.learnRF(X,Y[:,None])
print "OOB",oob
p  = rf.predictProbabilities(features)[:,1]

rag.showEdgeFeature(img, p)
vigra.show()
