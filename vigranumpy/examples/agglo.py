import vigra
import pylab
import numpy



def visualizeEdgeWeights(rag,coordMap,edgeWeightMap):

	edgeWeightImg = numpy.zeros(rag.topolocgialShape,dtype=numpy.float32)
	edgeWeightImg[:]=-1;
	#edgeWeightImg  = vigra.taggedView(edgeWeightImg,"xy")

	edgeWeightImg[:]=-1
	vigra.rag.edgeValueImage(
		rag,coordMap,edgeWeightMap,edgeWeightImg
	)
	edgeVal = numpy.where(edgeWeightImg>=-0.5)
	edgeWeightImg[edgeVal]-=edgeWeightImg[edgeVal].min()
	edgeWeightImg[edgeVal]/=edgeWeightImg[edgeVal].max()
	nonEdgeVal = numpy.where(edgeWeightImg<-0.5)
	edgeWeightImg[nonEdgeVal]=0.0
	pylab.imshow(numpy.swapaxes(edgeWeightImg,0,1), interpolation='nearest',cmap='hot')
	pylab.show()


f       = '100075.jpg'
sigma   = 1.5

print "get input"
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
gradmag = vigra.filters.gaussianGradientMagnitude(img,sigma)
gradmag = numpy.squeeze(gradmag).astype(numpy.float32)

if False:
    vigra.imshow(gradmag)
    pylab.show()

print "slic"
labels ,nseg = vigra.analysis.slicSuperpixels(vigra.colors.transform_RGB2Lab(img),10.0,5)
labels 		 = vigra.analysis.labelImage(labels)

print "get RAG"
rag = vigra.rag.Rag(labels)



print "edgeNum,nodeNum " ,rag.edgeNum,rag.nodeNum
print "maxEdgeId,maxNodeId" ,rag.maxEdgeId ,rag.maxNodeId


print "alloc maps"
edgeCoordMap     = vigra.rag.RagEdgeCoordinatesMap(rag)
edgeIndicatorMap = vigra.rag.RagEdgeFloatMap(rag)
edgeSizeMap      = vigra.rag.RagEdgeFloatMap(rag)

print "extract edge values"
vigra.rag.extractEdgeCoordinates(rag,edgeCoordMap)
vigra.rag.extractEdgeSizeFromCoords(rag,edgeCoordMap,edgeSizeMap)
vigra.rag.extractEdgeFeaturesFromImage(rag,edgeCoordMap,gradmag,edgeIndicatorMap)


print "visulize weights bevore ucm"
visualizeEdgeWeights(rag,edgeCoordMap,edgeIndicatorMap)
visualizeEdgeWeights(rag,edgeCoordMap,edgeSizeMap)

print "do ucm transform (inplace)"
vigra.rag.ucmTransform(rag,edgeIndicatorMap,edgeSizeMap)


print "visulize weights after ucm"
visualizeEdgeWeights(rag,edgeCoordMap,edgeSizeMap)
visualizeEdgeWeights(rag,edgeCoordMap,edgeIndicatorMap)

#visualizeEdgeWeights(rag,edgeSizeMap)
#visualizeEdgeWeights(rag,edgeSizeMap)
#visualizeEdgeWeights(rag,edgeSizeMap)



