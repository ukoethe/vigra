import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys


def showSeg(img,labels):
	crackedEdges = vigra.analysis.regionImageToCrackEdgeImage(labels)
	whereEdge    =  numpy.where(crackedEdges==0)
	imgToDisplay = vigra.resize(img,crackedEdges.shape)
	imgToDisplay-=imgToDisplay.min()
	imgToDisplay/=imgToDisplay.max()
	for c in range(img.ndim):
		ic = imgToDisplay[:,:,c]
		ic[whereEdge]=0

	vigra.imshow(imgToDisplay)



class MyOperator(object):
	def __init__(self,mergeGraph,graph,edgeIndicator,edgeSize,nodeSize):
		self.mergeGraph		= mergeGraph
		self.graph     		= graph
		self.edgeIndicator  = edgeIndicator
		self.edgeSize 		= edgeSize
		self.nodeSize       = nodeSize
		self.pq 			= vigra.utilities.ChangeablePriorityQueueFloat32Min(graph.maxEdgeId+1)

		edgeIds = mergeGraph.edgeIds()
		self.pq.push(edgeIds,self.edgeIndicator[edgeIds])

		assert len(self.pq)==mergeGraph.edgeNum

	def contractionEdge(self):
		eid = self.pq.top()
		while self.mergeGraph.edgeFromId(eid)==vigraph.INVALID:
			self.pq.deleteItem(eid)
			eid=self.pq.top()
		return self.mergeGraph.edgeFromId(eid)

	def contractionWeight(self):
		return float(self.pq.topPriority())

	def mergeNodes(self,a,b):
		# do the merging
		ia  = self.mergeGraph.id(a)
		ib  = self.mergeGraph.id(b)
		self.nodeSize[ia]+=self.nodeSize[ib]

	def mergeEdges(self,a,b):
		# do the merging
		ia  = self.mergeGraph.id(a)
		ib  = self.mergeGraph.id(b)
		sizeA = self.edgeSize[ia]
		sizeB = self.edgeSize[ib]
		self.edgeIndicator[ia] =  (self.edgeIndicator[ia]*sizeA + self.edgeIndicator[ib]*sizeB)/(sizeA+sizeB)
		self.edgeSize[ia]=sizeA+sizeB
		# delete b from pq
		self.pq.deleteItem(ib)
		print "ib",ib
		assert ia <= self.graph.maxEdgeId
		assert self.pq.contains(ia)
		self.pq.push(ia,float(self.edgeIndicator[ia]))


	def eraseEdge(self,a):
		eid  = self.mergeGraph.id(a)
		self.pq.deleteItem(eid)

print "get input"
f       = '100075.jpg'
f       = '69015.jpg'
f 		= '12003.jpg'
sigma   = 3.0
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab  = vigra.colors.transform_RGB2Lab(img)
newShape = [img.shape[0]*2-1,img.shape[1]*2-1]
#imgLab 	= vigra.resize(imgLab,newShape)
#img 	= vigra.resize(img,newShape)
gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
gradmag = numpy.squeeze(gradmag).astype(numpy.float32)
labels ,nseg = vigra.analysis.slicSuperpixels(imgLab,10.0,10)
labels 		 = numpy.squeeze(vigra.analysis.labelImage(labels))

numLabels = labels.max()





gridGraph 	   = vigraph.gridGraph(img.shape[0:2])


n1 = gridGraph.nodeFromId(500)

print "iteration test"

for nn in gridGraph.neighbourNodeIter(n1):
	print gridGraph.id(nn)


sys.exit(1)

rag,hyperEdges = vigraph.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)
hyperEdgeSizes = vigraph.hyperEdgeSizes(rag,hyperEdges)
hyperNodeSizes = vigraph.hyperNodeSizes(rag,gridGraph,labels)

edgeIndicator  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)
edgeMinWeight  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)
edgeIndicator  = vigraph.hyperEdgeImageFeatures(rag,gridGraph,hyperEdges,gradmag,edgeIndicator)

nodeFeatures   = vigraph.graphMap(graph=rag,item='node',dtype=numpy.float32,channels=3)
nodeFeatures   = vigraph.hyperNodeImageFeatures(rag,gridGraph,labels,imgLab,nodeFeatures)
mergeGraph 	   = vigraph.mergeGraph(rag)


print "get the pure python op"

myop 			 = MyOperator(mergeGraph,rag,edgeIndicator=edgeIndicator,edgeSize=hyperEdgeSizes,nodeSize=hyperNodeSizes)			
clusterOperator  = vigraph.pythonClusterOperator(mergeGraph=mergeGraph, opertator=myop)
print "get the cluster class"
hc = vigraph.hierarchicalClustering(clusterOperator,100)

print "do clustering"
hc.cluster()
newLabels = labels.copy()
newLabels = newLabels.reshape(-1)
# this is inplace
hc.reprNodeIds(newLabels)
newLabels = newLabels.reshape(labels.shape)
labels    = vigra.analysis.labelImage(newLabels)

showSeg(img,labels)
pylab.show()