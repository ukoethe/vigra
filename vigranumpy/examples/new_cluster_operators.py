


def chiSquared(a,b):
    pass

def eukledian(a,b):
    pass

def cityblock(a,b):
    pass





def linearScaling(val,scale):
    return scale*nodeDist

def logScaleing(val,scale):
    return scale*math.log(val)

def sqrtScaling(val,scale):
    return scale*math.sqrt(val):






class MyOperator(object):
	def __init__(self,mergeGraph,graph,edgeIndicator,edgeSize,nodeSize,nodeFeatures):
		self.mergeGraph		= mergeGraph
		self.graph     		= graph
		self.edgeIndicator  = edgeIndicator
		self.edgeSize 		= edgeSize
		self.nodeSize       = nodeSize
		self.nodeFeatures   = nodeFeatures
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

    def _mergeItems(ia,ib,features,size):
        # do the merging
        ia  = self.mergeGraph.id(a)
        ib  = self.mergeGraph.id(b)
        sizeA = size[ia]
        sizeB = size[ib]
        self.edgeIndicator[ia] =  (features[ia]*sizeA + features[ib]*sizeB)/(sizeA+sizeB)
        size[ia]=sizeA+sizeB

    def _getEdgeWeight(self,ei,ui,vi):
        return self.edgeIndicator[ei]

	def mergeNodes(self,a,b):
		self._mergeItems(self.mergeGraph.id(a),self.mergeGraph.id(b),self.nodeFeatures,self.nodeSize)

	def mergeEdges(self,a,b):
		# do the merging
		ia  = self.mergeGraph.id(a)
		ib  = self.mergeGraph.id(b)
		self._mergeItems(ia,ib,self.nodeFeatures,self.nodeSize)
		# delete b from pq
		self.pq.deleteItem(ib)
		

    def eraseEdge(self,deletedEdge):
        mg = self.mergeGraph
        pq = self.pq
        deletedEdgeId  = mg.id(deletedEdge)
        self.pq.deleteItem(deletedEdgeId)

        # get the node the erased edge is within (u==v for inactive/dead edges)
        node   = self.inactiveEdgesNode(edge)
        nodeId = mg.id(node)
        # iterate over all neighbour nodes (and edges connecting)
        for otherNode,edge in zip(mg.neighbourNodeIter(node),mg.incEdgeIter(node)):
            # compute new weight
            edgeId = mg.id(edge)
            newWeight = self._getEdgeWeight(mg.id(edge),nodeId,mg.id(otherNode))
            pq.push(edgeId,newWeight)

