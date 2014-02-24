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


f       = '100075.jpg'
f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 3.0

print "prepare input"
img                 = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab              = vigra.colors.transform_RGB2Lab(img)
imgLabInterpolated  = vigra.resize(imgLab,[imgLab.shape[0]*2-1,imgLab.shape[1]*2-1 ])
gradmagInterpolated = vigra.filters.gaussianGradientMagnitude(imgLabInterpolated,sigma)
labels ,nseg        = vigra.analysis.slicSuperpixels(imgLab,10.0,5)
labels              = vigra.analysis.labelImage(labels)



graph0,graph1 = vigraph.gridRegionAdjacencyGraph(labels=labels,ignoreLabel=None)




# get grid graph and edge weights
graph0EdgeWeights =  vigraph.edgeFeaturesFromInterpolatedImage(graph0,gradmagInterpolated)







def makeEdgeIndicator(
    graph,
    edgeIndicator = None,
    nodeFeatures  = None,
    beta          = 0.5,
    metric        = 'l1',
    normalize     = False,   
    **kwargs
):
    assert edgeIndicator is not None or nodeFeatures is not None

    if nodeFeatures is None:
        if normalize:
            pass
        else:
            return edgeIndicator
    else:
        nodeDist = graphs.nodeFeatureDistToEdgeWeight(graph,nodeFeatures,metric='l1')
        if edgeIndicator is not None :
            ei =edgeIndicator.copy()
            ei *=(1.0-beta)
            ei += beta*nodeDist
            return ei
        else :
            return nodeDist


def getNodeSizes(g):
    return g.nodeSize()
def getEdgeLength(g):
    return g.edgeLength()

def graphSegmentation(  graph,edgeIndicator = None,nodeFeatures  = None,
                        method='fw',beta= 0.1,metric = 'squaredNorm',
                        nodeNumStop   = None,**kwargs
) :
    global img
    ##############################################
    # preprocessing
    ##############################################
    assert edgeIndicator is not None or nodeFeatures is not None

    if nodeNumStop is None :
        nodeNumStop = -1

    labels = None 

    ##############################################
    # segmentation
    ##############################################
    if method in ("fw","felzenszwalb") :
        edgeInd = makeEdgeIndicator(graph,edgeIndicator,nodeFeatures,beta,metric,**kwargs)
        nodeSize = getNodeSizes(graph)
        labels = graphs.felzenszwalbSegmentation(graph,edgeInd,nodeSizes=nodeSize,k=0.1,nodeNumStop=nodeNumStop)


    if method in ("shfw",):


        g1 = graph
        done=False

        eF = edgeIndicator
        nF = nodeFeatures
        counter=0
        while not done:
            sCond = max(int(g1.nodeNum*0.75),1)
            #print "sCond",sCond
            if sCond <=  nodeNumStop:
                sCond = nodeNumStop
                done=True
            labels = graphSegmentation(g1,eF,nF,'fw',beta,metric,sCond,**kwargs)
            g2  = graphs.regionAdjacencyGraph(g1,labels=labels)
            eF = g2.accumulateEdgeFeatures(eF,acc='mean')
            nF = g2.accumulateNodeFeatures(nF,acc='mean')
            g1=g2
            counter+=1
        labels = g2.projectLabelsBack(steps=counter)

           





    if method in ("mc","multicut"):
        edgeInd = makeEdgeIndicator(graph,edgeIndicator,nodeFeatures,**kwargs)

    if method in ("ws","watersheds"):
        edgeInd = makeEdgeIndicator(graph,edgeIndicator,nodeFeatures,**kwargs)

    if method in ("hc","hierarchicalClustering"):
        nodeSize  = getNodeSizes(graph)
        edgeLength = getEdgeLength(graph)
        mg = graphs.mergeGraph(graph)
        clusterOp = graphs.minEdgeWeightNodeDist(mg,edgeIndicator,edgeLength,nodeFeatures,nodeSize,
            beta=float(beta),nodeDistType=metric,wardness=1.0)
        hc = graphs.hierarchicalClustering(clusterOp,nodeNumStopCond=nodeNumStop)
        hc.cluster()
        labels = hc.resultLabels()

    ##############################################
    # postprocessing
    ##############################################
    return labels


graph1EdgeWeights = graph1.accumulateEdgeFeatures(graph0EdgeWeights,acc='mean')
graph1NodeFeatures = graph1.accumulateNodeFeatures(img,acc='mean')




graph1Labels = graphSegmentation(graph1,graph1EdgeWeights,graph1NodeFeatures,method='fw',nodeNumStop=10)   
graph2       = vigraph.regionAdjacencyGraph(graph=graph1,labels=graph1Labels,ignoreLabel=None)

graph2EdgeWeights = graph2.accumulateEdgeFeatures(graph1EdgeWeights,acc='mean')
graph2NodeFeatures = graph2.accumulateNodeFeatures(graph1NodeFeatures,acc='mean')


graph2.show(img)
vigra.show()

