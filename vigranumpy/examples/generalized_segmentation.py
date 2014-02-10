


def makeEdgeWeights(graph,edgeIndicators,nodeFeatures,metric,beta):
    if nodeFeatures is None:
        return edgeIndicators
    else :
        nodeDiff = vigraph.nodeFeatureDistToEdgeWeight(graph,nodeFeatures,metric=metric)
        return (1.0-beta)*nodeDiff + beta*edgeIndicators



def segment(
    graph
    edgeIndicators  = None,
    nodeFeatures    = None,
    edgeSizes       = None,
    nodeSizes       = None,
    seeds           = None,
    metric          = None,
    beta            = 0.5,
    nodeNum         = None,
    method          = None,
    segParam        = None,
    refinement      = True,
    refinementParam = None 
):
    """
    edgeIndicators --  a floating point number for each edge where 
        high indicated to keep the nodes of the edges seperated

    nodeFeatures -- feature vector for each nodes

    edgeSize -- length of the edge 

    nodeSize -- size of the nodes

    seeds    -- seeds / segmetation constraints

    metric   -- metric to get an edge indicator from two adjacent nodes:
        The metric can be:
            - 'l1','manhatten'
            - 'eucledian','norm'
            - 'squaredNorm'
            - 'chiSquared'

    beta -- mixture between nodeDistance edge Indicator and edge indicator:
        Zero means no nodeDistance edge indicator,one means only nodeDistance.
        Any beta between zero and one is a linear combination between both. 

    nodeNum -- result segmetation node number


    method -- segmetation method:
        - 'wardAgglomeration'
        - 'felzenszwalb'
        - 'mst'
        - 'multicut'
        - 'watersheds'

    segParam -- parameter of the selected segmetation method:

        - 'wardAgglomeration':
            Parameter of Ward Agglomeration:
                - useLogWardnes -- scale node sizes linear of logarithmic
                - wardness -- 0 means no wardness 1 means full wardness

        - 'felzenszwalb':
            Parameter of Felzenszwalb
                - k -- a high k will lead to more equal sized region

        - 'mst'
            Parameter free

        - 'multicut'
            Parameter of multicut
                - threshold -- an edge weight below threshold indicates 
                    that nodes of edges should be in same segment.
                    Threshold can be uninformative since it will be
                    tuned to fit nodeNum parameter.
                    But a well choosen threshold will speed up the algorithm.
                - mcAlg -- algorithm for Multicut segmentation:
                    - 'cgc' -- fast approximative multicut solver
                    - 'multicut' -- gloval optimal multicut solver


        - 'watersheds'

        - 'cgc'

    refinement -- refine the resuting segmetation with seeded wathersheds ? :

    """
    
    out = None

    if method == 'felzenszwalb':
        mixedEdgeWeights = makeEdgeWeights(graph,edgeIndicators,nodeFeatures,metric,beta)
        out = vigraph.felzenszwalbSegmentation(edgeWeights=mixedEdgeWeights,nodeSizes=nodeSizes,
            k=param.k,nodeNum= nodeNum)

    elif :
        pass

    elif 'cgc':
        mixedEdgeWeights = makeEdgeWeights(graph,edgeIndicators,nodeFeatures,metric,beta)

    # refinement
    if refinement :

        # segmetation to seeds
