import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys
import matplotlib
import pylab

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






def chiSquared(a,b):
    pass

def eukledian(a,b):
    return numpy.sqrt(numpy.sum((a-b)**2))

def cityblock(a,b):
    pass




























class MyOperator(object):
    def __init__(self,mergeGraph,graph,edgeIndicator,edgeSize,nodeSize,nodeFeatures,nodeDist=eukledian):
        self.mergeGraph     = mergeGraph
        self.graph          = graph
        self.edgeIndicator  = edgeIndicator
        self.edgeSize       = edgeSize
        self.nodeSize       = nodeSize
        self.nodeFeatures   = nodeFeatures
        self.pq             = vigra.utilities.ChangeablePriorityQueueFloat32Min(graph.maxEdgeId+1)

        self.nodeDist       = nodeDist

        edgeIds = mergeGraph.edgeIds()

        for i in edgeIds :
            ui,vi = self.mergeGraph.uvId(long(i))
            self.pq.push(long(i),self._getEdgeWeight(i,ui,vi))

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

    def _mergeItems(self,ia,ib,features,size):
        # do the merging
        sizeA = size[ia]
        sizeB = size[ib]
        features[ia] =  (features[ia]*sizeA + features[ib]*sizeB)/(sizeA+sizeB)
        size[ia]=sizeA+sizeB

    def _getEdgeWeight(self,ei,ui,vi):
        fromEdge  = self.edgeIndicator[ei]
        fromNodes = numpy.sqrt(self.nodeDist(self.nodeFeatures[ui],self.nodeFeatures[vi]))

        w = float(fromEdge)+float(fromNodes)
        ward  = 1.0/numpy.log(self.nodeSize[ui]) + 1.0/numpy.log(self.nodeSize[vi])
        return w
        return float(w/ward)

    def mergeNodes(self,a,b):
        self._mergeItems(self.mergeGraph.id(a),self.mergeGraph.id(b),self.nodeFeatures,self.nodeSize)

    def mergeEdges(self,a,b):
        # do the merging
        ia  = self.mergeGraph.id(a)
        ib  = self.mergeGraph.id(b)
        self._mergeItems(ia,ib,self.edgeIndicator,self.edgeSize)
        # delete b from pq
        self.pq.deleteItem(ib)
        
    def eraseEdge(self,deletedEdge):
        mg = self.mergeGraph
        pq = self.pq
        deletedEdgeId  = mg.id(deletedEdge)
        self.pq.deleteItem(deletedEdgeId)

        # get the node the erased edge is within (u==v for inactive/dead edges)
        node   = self.mergeGraph.inactiveEdgesNode(deletedEdge)
        nodeId = mg.id(node)
        # iterate over all neighbour nodes (and edges connecting)
        for otherNode,edge in zip(mg.neighbourNodeIter(node),mg.incEdgeIter(node)):
            # compute new weight
            edgeId = mg.id(edge)
            newWeight = self._getEdgeWeight(mg.id(edge),nodeId,mg.id(otherNode))
            pq.push(edgeId,newWeight)












print "get input"
f       = '100075.jpg'
f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 3.0
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab  = vigra.colors.transform_RGB2Lab(img)
newShape = [img.shape[0]*2-1,img.shape[1]*2-1]
#imgLab     = vigra.resize(imgLab,newShape)
#img    = vigra.resize(img,newShape)
gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
gradmag = numpy.squeeze(gradmag).astype(numpy.float32)
labels ,nseg = vigra.analysis.slicSuperpixels(imgLab,10.0,10)
labels       = numpy.squeeze(vigra.analysis.labelImage(labels))

numLabels = labels.max()





gridGraph      = vigraph.gridGraph(img.shape[0:2])
rag,hyperEdges = vigraph.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)
hyperEdgeSizes = vigraph.hyperEdgeSizes(rag,hyperEdges)
hyperNodeSizes = vigraph.hyperNodeSizes(rag,gridGraph,labels)

edgeIndicator  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)
edgeMinWeight  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)
edgeIndicator  = vigraph.hyperEdgeImageFeatures(rag,gridGraph,hyperEdges,gradmag,edgeIndicator)

nodeFeatures   = vigraph.graphMap(graph=rag,item='node',dtype=numpy.float32,channels=3)
nodeFeatures   = vigraph.hyperNodeImageFeatures(rag,gridGraph,labels,imgLab,nodeFeatures)
mergeGraph     = vigraph.mergeGraph(rag)


print "get the pure python op"

myop             = MyOperator(  mergeGraph,rag,edgeIndicator=edgeIndicator,edgeSize=hyperEdgeSizes,nodeSize=hyperNodeSizes,nodeFeatures=nodeFeatures,
                                nodeDist = eukledian
                                )         
clusterOperator  = vigraph.pythonClusterOperator(mergeGraph=mergeGraph, opertator=myop)
print "get the cluster class"
hc = vigraph.hierarchicalClustering(clusterOperator,50)

print "do clustering"
hc.cluster()
newLabels = labels.copy()
newLabels = newLabels.reshape(-1)
# this is inplace
hc.reprNodeIds(newLabels)
newLabels = newLabels.reshape(labels.shape)
labels    = vigra.analysis.labelImage(newLabels)

assert labels.min()==1


showSeg(img,labels)
pylab.show()


labeling ,nseg = vigra.analysis.watershedsReoptimization(labels,gradmag,6,visu=True)

showSeg(img,labeling)
pylab.show()