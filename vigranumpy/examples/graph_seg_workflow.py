import vigra
import vigra.graphs as vigraph
import vigra.graphs as graphs
import pylab
import numpy
import sys
import matplotlib
import pylab
import matplotlib.pyplot as plt
from   matplotlib.widgets import Button
import math
import pylab

f       = '100075.jpg'
f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 2.0

print "prepare input"
img                 = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab              = vigra.colors.transform_RGB2Lab(img)
imgLabInterpolated  = vigra.resize(imgLab,[imgLab.shape[0]*2-1,imgLab.shape[1]*2-1 ])
gradmagInterpolated = vigra.filters.gaussianGradientMagnitude(imgLabInterpolated,sigma)
labels ,nseg        = vigra.analysis.slicSuperpixels(imgLab,20.0,20)
labels              = vigra.analysis.labelImage(labels)



graph0,graph1 = vigraph.gridRegionAdjacencyGraph(labels=labels,ignoreLabel=None)





# get grid graph and edge weights
graph0EdgeWeights =  vigraph.edgeFeaturesFromInterpolatedImage(graph0,gradmagInterpolated)






def getNodeSizes(g):
    return g.nodeSize()
def getEdgeLength(g):
    return g.edgeLength()

def makeEdgeIndicator(  graph,edgeIndicator=None,nodeFeatures=None,
                        beta = 0.5,metric='l1',wardness=1.0
):

    assert edgeIndicator is not None or nodeFeatures is not None
    eI = None
    if nodeFeatures is None:
        eI = edgeIndicator
    else:
        nodeDist = graphs.nodeFeatureDistToEdgeWeight(graph,nodeFeatures,metric='l1')
        if edgeIndicator is not None :
            eI =edgeIndicator.copy()
            eI *=(1.0-beta)
            eI += beta*nodeDist
        else :
            eI =  nodeDist
    if wardness <=0.00000001:
        return eI
    else:
        nodeSize = getNodeSizes(graph)
        #print wardness
        #print "MAX NODE SIZE ",nodeSize.max()
        #print nodeSize
        #print nodeSize.shape,nodeSize.dtype
        #print eI.shape,eI.dtype
        eWard=graphs.wardCorrection(graph,eI,nodeSize,wardness)

        #print "diff" ,eWard-eI

        return eWard


def makeWeights(eI,gamma):
    e1 = numpy.exp(-gamma*eI)
    e0 = 1.0 - e1
    weights = e1-e0
    return weights

def multicutSegmentation(graph,edgeIndicator,nodeNumStop=None,gamma=0.1):



    # check for opengm version
    try:
        import opengm
    except:
        raise RuntimeError("'import opengm' failed : multicutSegmentation needs 'opengm' ")
    Cgc = None
    try :
        Cgc = opengm.inference.Cgc
    except:
        raise RuntimeError("'wrong opengm version/branch, Cgc solver missing ")


    # make graphical model
    vis,eInd = graphs.opengmMulticutDataStructure(graph,edgeIndicator)

    weights=makeWeights(eInd,0.1)


    numVar = graph.nodeNum
    gm = opengm.gm(numpy.ones(numVar,opengm.label_type)*numVar)

    pf = opengm.pottsFunctions([numVar,numVar],numpy.zeros(1),weights)

    fids = gm.addFunctions(pf)
    gm.addFactors(fids,vis)
    cgc = Cgc(gm)   

    #cgc.infer(visitor=cgc.verboseVisitor())
    def run(eInd,gamma,arg=None):
        weights=makeWeights(eInd,gamma)
        cgc.changeWeights(weights)

        if arg is not None:
            cgc.setStartingPoint(arg)

        #cgc.infer()
        cgc.infer(cgc.verboseVisitor())
        arg=cgc.arg()+1
        numSeg = len(numpy.bincount(arg.astype(numpy.int32)))
        return arg,numSeg

    bestUpper = [10.00000, None ] 
    bestLower = [0.0000000  , None ]
    bestMixed = [(bestUpper[0]+bestLower[0])/2.0  , None ]

    cgc.infer()
    arg=cgc.arg()+1
    

    numSeg = len(numpy.bincount(arg.astype(numpy.int32)))
    #print "numSeg" ,numSeg

    bestUpper[1]=run(bestUpper[0],eInd)[1]
    bestLower[1]=run(bestLower[0],eInd)[1]

    assert gamma < bestUpper[0]
    assert gamma > bestLower[0]

    toMuchSegs=False
    converged=False
    while bestMixed[1] is None or bestMixed[1]!=nodeNumStop:

        assert bestLower[1]<nodeNumStop
        assert bestUpper[1]>nodeNumStop

        if bestMixed[1] is None :
            bestMixed = [gamma , None ]
        else :
            oldGamma = bestMixed[0]
            bestMixed = [(bestUpper[0]+bestLower[0])/2.0  , None ]
            # check for convergence
            if numpy.abs(oldGamma-bestMixed[0])<0.00000001:
                converged=True
                ##print "converged"
                if numSeg < nodeNumStop:
                    #print "more segmented needed new gamma will be" 
                    arg,numSeg=run(eInd,bestUpper[0],arg)
                    assert numSeg>nodeNumStop
                elif  numSeg > nodeNumStop:
                    toMuchSegs=True
                labels=graphs.opengmArgToLabeling(graph,arg.astype(numpy.uint32))
                break

        arg,numSeg=run(eInd,bestMixed[0])
        bestMixed[1]=numSeg
       
        print numSeg ,"gamma",bestMixed[0],'(', bestLower[0],bestUpper[0],')'


        
        if numSeg < nodeNumStop:
            bestLower=[bestMixed[0],numSeg]
            #print "more segmented needed new gamma will be",(bestUpper[0]+bestLower[0])/2.0
           
        elif  numSeg > nodeNumStop:
            bestUpper=[bestMixed[0],numSeg]
            #print "less segmented needed new gamma will be",(bestUpper[0]+bestLower[0])/2.0
            
        else :
            labels=graphs.opengmArgToLabeling(graph,arg.astype(numpy.uint32))
            break

    if toMuchSegs :
        labels=graphs.opengmArgToLabeling(graph,arg.astype(numpy.uint32))
        g2  = graphs.regionAdjacencyGraph(graph,labels=labels)
        eF = g2.accumulateEdgeFeatures(edgeIndicator,acc='mean')
        labels2 = graphSegmentation(g2,eF,method='hc',nodeNumStop=nodeNumStop)
        lagels=g2.projectLabelsBack(steps=1)

    return labels

def graphSegmentation(  graph,edgeIndicator=None,nodeFeatures=None,
                        method='fw',beta= 0.1,metric = 'squaredNorm',
                        nodeNumStop=None,k=0.1,wardness=1.0,
                        makeHierarchical=True) :

    ##############################################
    # preprocessing
    ##############################################
    assert edgeIndicator is not None or nodeFeatures is not None

    if nodeNumStop is None :
        nodeNumStop = -1
    labels = None 

    if method in ('hc','hierarchicalClustering'):
        makeHierarchical=False
    ##############################################
    # segmentation
    ##############################################

    if makeHierarchical==False:
        if method in ("fw","felzenszwalb") :

            edgeInd = makeEdgeIndicator(graph=graph,edgeIndicator=edgeIndicator,
                                        nodeFeatures=nodeFeatures,beta=beta,
                                        metric=metric,wardness=wardness)
            nodeSize = getNodeSizes(graph)
            labels = graphs.felzenszwalbSegmentation(graph,edgeInd,nodeSizes=nodeSize,k=k,nodeNumStop=nodeNumStop)
        elif method in ("mc","multicut"):
            edgeInd = makeEdgeIndicator(graph=graph,edgeIndicator=edgeIndicator,
                                        nodeFeatures=nodeFeatures,beta=beta,
                                        metric=metric,wardness=wardness)
            labels=multicutSegmentation(graph=graph,edgeIndicator=edgeInd,gamma=0.2,nodeNumStop=nodeNumStop)


        elif method in ("hc","hierarchicalClustering"):
            nodeSize  = getNodeSizes(graph)
            edgeLength = getEdgeLength(graph)
            mg = graphs.mergeGraph(graph)
            clusterOp = graphs.minEdgeWeightNodeDist(mg,edgeIndicator,edgeLength,nodeFeatures,nodeSize,
                beta=float(beta),nodeDistType=metric,wardness=wardness)
            hc = graphs.hierarchicalClustering(clusterOp,nodeNumStopCond=nodeNumStop)
            hc.cluster()
            labels = hc.resultLabels()

        else :
            raise RuntimeError("unknown method :'%s'  try: 'hc','fw','mc' "%str(method))

    else :
        g1 = graph
        done=False

        eF = edgeIndicator
        nF = nodeFeatures
        counter=0
        while not done:
            sCond = max(int(g1.nodeNum*0.9),1)
            print "sCond",sCond
            if sCond <=  nodeNumStop:
                sCond = nodeNumStop
                done=True
            labels = graphSegmentation( graph=g1,edgeIndicator = eF,nodeFeatures  = nF,
                                        method=method,beta=beta,metric = metric,
                                        nodeNumStop   = sCond ,k=k,wardness=wardness,
                                        makeHierarchical=False)
            g2  = graphs.regionAdjacencyGraph(g1,labels=labels)
            eF = g2.accumulateEdgeFeatures(eF,acc='mean')
            nF = g2.accumulateNodeFeatures(nF,acc='mean')
            g1=g2
            counter+=1
        labels = g2.projectLabelsBack(steps=counter)

    ##############################################
    # postprocessing
    ##############################################
    return labels


graph1EdgeWeights = graph1.accumulateEdgeFeatures(graph0EdgeWeights,acc='mean')
graph1NodeFeatures = graph1.accumulateNodeFeatures(img,acc='mean')




l_hc = graphSegmentation(   graph1,graph1EdgeWeights,graph1NodeFeatures,
                                    method='hc',nodeNumStop=5,beta=0.01,k=100.0,
                                    wardness=1.0,makeHierarchical=True)   

l_mc = graphSegmentation(   graph1,graph1EdgeWeights,graph1NodeFeatures,
                                    method='mc',nodeNumStop=5,beta=0.01,k=100.0,
                                    wardness=1.0,makeHierarchical=True)   

f = pylab.figure()
f.add_subplot(2, 1, 0)
graph1.showNested(img,l_hc)
f.add_subplot(2, 1, 1)
graph1.showNested(img,l_mc)
pylab.show()
