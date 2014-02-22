



def makeEdgeIndicator(
	graph,
	edgeIndicator = None,
	nodeFeatures  = None,
	**kwargs
):
	pass


def graphSegmentation(
	graph,
	edgeIndicator = None,
	nodeFeatures  = None
	method,
	**kwargs
) :
	
	if method in ("fw","felzenszwalb") :

		# get the mixed edge indicator
		edgeInd = makeEdgeIndicator(graph,edgeIndicator,nodeFeatures,**kwargs)
		

	if method in ("mc","multicut"):
		# get the mixed edge indicator
		edgeInd = makeEdgeIndicator(graph,edgeIndicator,nodeFeatures,**kwargs)

	if method in ("ws","watersheds"):
		# get the mixed edge indicator
		edgeInd = makeEdgeIndicator(graph,edgeIndicator,nodeFeatures,**kwargs)

	if method in ("hc","hierarchicalClustering"):
		pass
	