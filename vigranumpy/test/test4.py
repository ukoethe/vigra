#######################################################################
#                                                                      
#         Copyright 2009-2010 by Ullrich Koethe                        
#                                                                      
#    This file is part of the VIGRA computer vision library.           
#    The VIGRA Website is                                              
#        http://hci.iwr.uni-heidelberg.de/vigra/                       
#    Please direct questions, bug reports, and contributions to        
#        ullrich.koethe@iwr.uni-heidelberg.de    or                    
#        vigra@informatik.uni-hamburg.de                               
#                                                                      
#    Permission is hereby granted, free of charge, to any person       
#    obtaining a copy of this software and associated documentation    
#    files (the "Software"), to deal in the Software without           
#    restriction, including without limitation the rights to use,      
#    copy, modify, merge, publish, distribute, sublicense, and/or      
#    sell copies of the Software, and to permit persons to whom the    
#    Software is furnished to do so, subject to the following          
#    conditions:                                                       
#                                                                      
#    The above copyright notice and this permission notice shall be    
#    included in all copies or substantial portions of the             
#    Software.                                                         
#                                                                      
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   
#    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       
#    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      
#    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      
#    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     
#    OTHER DEALINGS IN THE SOFTWARE.                                   
#                                                                      
#######################################################################

import sys
print >> sys.stderr, "\nexecuting test file", __file__
execfile('set_paths.py')


from nose.tools import assert_equal, raises
from vigra import numpy as np
from vigra import numpy as numpy
from vigra import graphs as vigraph
from vigra import graphs,taggedView
import vigra

def testGridGraphSegmentationFelzenszwalbSegmentation():
    dataRGB  = numpy.random.random([3,3,3]).astype(numpy.float32)
    dataRGB  = taggedView(dataRGB,'xyc')
    data  = numpy.random.random([3,3]).astype(numpy.float32)
    edata = numpy.random.random([3*2-1,3*2-1]).astype(numpy.float32)
    g0 = graphs.gridGraph(data.shape)

    ew = graphs.edgeFeaturesFromInterpolatedImage(g0,edata)

    labels = graphs.felzenszwalbSegmentation(graph=g0,edgeWeights=ew,k=1.0,nodeNumStop=5)

    g1  = graphs.regionAdjacencyGraph(graph=g0,labels=labels)
    assert g1.nodeNum == 5


    
    data  = numpy.random.random([3,3,3]).astype(numpy.float32)
    edata = numpy.random.random([3*2-1,3*2-1,3*2-1]).astype(numpy.float32)
    g0 = graphs.gridGraph(data.shape)

    ew = graphs.edgeFeaturesFromInterpolatedImage(g0,edata)

    labels = graphs.felzenszwalbSegmentation(graph=g0,edgeWeights=ew,k=1.0,nodeNumStop=15)

    g1  = graphs.regionAdjacencyGraph(graph=g0,labels=labels)
    assert g1.nodeNum == 15

def testGridGraphWatersheds():

    data  = numpy.random.random([10,10,10]).astype(numpy.float32)
    edata = numpy.random.random([10*2-1,10*2-1,10*2-1]).astype(numpy.float32)
    g0 = graphs.gridGraph(data.shape)


    ew = graphs.edgeFeaturesFromInterpolatedImage(graph=g0,image=edata)

    # generate seeds
    seeds = graphs.nodeWeightedWatershedsSeeds(graph=g0,nodeWeights=data)
    # node weighted watershed seeds
    labelsNodeWeightedA  = graphs.nodeWeightedWatersheds(graph=g0,nodeWeights=data,seeds=seeds)
    # node weighted watershed seeds
    labelsNodeWeightedB  = graphs.nodeWeightedWatersheds(graph=g0,nodeWeights=data)
    # edge weighted watershed seeds
    seeds = graphs.nodeWeightedWatershedsSeeds(graph=g0,nodeWeights=data)
    labelsEdgeWeighted  = graphs.edgeWeightedWatersheds(graph=g0,edgeWeights=ew,seeds=seeds)

    assert numpy.array_equal(labelsNodeWeightedA,labelsNodeWeightedB)

    data  = numpy.random.random([10,10]).astype(numpy.float32)
    edata = numpy.random.random([10*2-1,10*2-1]).astype(numpy.float32)
    g0 = graphs.gridGraph(data.shape)


    ew = graphs.edgeFeaturesFromInterpolatedImage(graph=g0,image=edata)

    # generate seeds
    seeds = graphs.nodeWeightedWatershedsSeeds(graph=g0,nodeWeights=data)
    # node weighted watershed seeds
    labelsNodeWeightedA  = graphs.nodeWeightedWatersheds(graph=g0,nodeWeights=data,seeds=seeds)
    # node weighted watershed seeds
    labelsNodeWeightedB  = graphs.nodeWeightedWatersheds(graph=g0,nodeWeights=data)
    # edge weighted watershed seeds
    labelsEdgeWeighted  = graphs.edgeWeightedWatersheds(graph=g0,edgeWeights=ew,seeds=seeds)

    assert numpy.array_equal(labelsNodeWeightedA,labelsNodeWeightedB)


def testGridGraphAgglomerativeClustering():
    dataRGB  = numpy.random.random([10,10,3]).astype(numpy.float32)
    dataRGB  = vigra.taggedView(dataRGB,'xyc')
    data  = numpy.random.random([10,10]).astype(numpy.float32)
    edata = numpy.random.random([10*2-1,10*2-1]).astype(numpy.float32)
    g0 = graphs.gridGraph(data.shape)


    ew = graphs.edgeFeaturesFromInterpolatedImage(graph=g0,image=edata)
    #ew = taggedView(ew,'xyz')

    
    labels = graphs.agglomerativeClustering(graph=g0,edgeWeights=ew,nodeFeatures=dataRGB,nodeNumStop=5)
    g1  = graphs.regionAdjacencyGraph(graph=g0,labels=labels)
    assert g1.nodeNum == 5
    
    labels = graphs.agglomerativeClustering(graph=g0,edgeWeights=ew,nodeNumStop=5)
    g1  = graphs.regionAdjacencyGraph(graph=g0,labels=labels)
    assert g1.nodeNum == 5



    dataRGB  = numpy.random.random([10,10,10,3]).astype(numpy.float32)
    dataRGB  = vigra.taggedView(dataRGB,'xyzc')
    data  = numpy.random.random([10,10,10]).astype(numpy.float32)
    edata = numpy.random.random([10*2-1,10*2-1,10*2-1]).astype(numpy.float32)
    g0 = graphs.gridGraph(data.shape)


    ew = graphs.edgeFeaturesFromInterpolatedImage(graph=g0,image=edata)
    #ew = taggedView(ew,'xyz')

    
    labels = graphs.agglomerativeClustering(graph=g0,edgeWeights=ew,nodeFeatures=dataRGB,nodeNumStop=5)
    g1  = graphs.regionAdjacencyGraph(graph=g0,labels=labels)
    assert g1.nodeNum == 5
    
    labels = graphs.agglomerativeClustering(graph=g0,edgeWeights=ew,nodeNumStop=5)
    g1  = graphs.regionAdjacencyGraph(graph=g0,labels=labels)
    assert g1.nodeNum == 5

class TestGraph(object):

    def testAddNodesWithIds(self):

        IV = vigraph.INVALID

        g  = vigraph.listGraph()
        assert g.nodeNum == 0
        assert g.edgeNum == 0
        assert g.nodeFromId(0)==IV
        assert g.nodeFromId(1)==IV
        assert g.nodeFromId(2)==IV
        assert g.nodeFromId(3)==IV
        assert g.nodeFromId(4)==IV
        assert g.nodeFromId(5)==IV
        assert g.nodeFromId(6)==IV
        assert g.nodeFromId(7)==IV

        n5 =g.addNode(5)

        assert g.id(n5)==5

        assert n5!=IV
        assert g.nodeNum == 1
        assert g.edgeNum == 0
        assert g.nodeFromId(0)==IV
        assert g.nodeFromId(1)==IV
        assert g.nodeFromId(2)==IV
        assert g.nodeFromId(3)==IV
        assert g.nodeFromId(4)==IV
        assert g.nodeFromId(5)!=IV
        assert g.nodeFromId(6)==IV
        assert g.nodeFromId(7)==IV

        n2=g.addNode(2)
        assert n2!=IV
        assert g.nodeNum == 2
        assert g.nodeFromId(0)==IV
        assert g.nodeFromId(1)==IV
        assert g.nodeFromId(2)!=IV
        assert g.nodeFromId(3)==IV
        assert g.nodeFromId(4)==IV
        assert g.nodeFromId(5)!=IV
        assert g.nodeFromId(6)==IV
        assert g.nodeFromId(7)==IV

        n6=g.addNode(6)
        assert n6!=IV
        assert g.nodeNum == 3
        assert g.nodeFromId(0)==IV
        assert g.nodeFromId(1)==IV
        assert g.nodeFromId(2)!=IV
        assert g.nodeFromId(3)==IV
        assert g.nodeFromId(4)==IV
        assert g.nodeFromId(5)!=IV
        assert g.nodeFromId(6)!=IV
        assert g.nodeFromId(7)==IV

        n6=g.addNode(5)
        n2=g.addNode(2)
        n6=g.addNode(6)
        assert g.nodeNum == 3
        assert g.nodeFromId(0)==IV
        assert g.nodeFromId(1)==IV
        assert g.nodeFromId(2)!=IV
        assert g.nodeFromId(3)==IV
        assert g.nodeFromId(4)==IV
        assert g.nodeFromId(5)!=IV
        assert g.nodeFromId(6)!=IV
        assert g.nodeFromId(7)==IV

    def testAddEdges(self):
        IV = vigraph.INVALID
        elist = [
            [1,3],
            [1,5],
            [5,7],
            [3,4]
        ] 

        edges = np.array(elist,dtype=np.uint32)
        nodeIds = np.unique(edges.reshape(-1))


        g = vigraph.listGraph()

        edgeIds = g.addEdges(edges)

        assert g.edgeNum == len(elist)
        assert g.edgeNum == len(edgeIds)
        assert g.nodeNum == len(nodeIds)
        assert g.maxNodeId == nodeIds.max()

        for ui,vi in elist :
            assert g.findEdge(ui,vi)!=IV
            assert g.findEdge(g.nodeFromId(ui),g.nodeFromId(vi))!=IV

        for nId in nodeIds :
            nId = int(nId)
            assert g.nodeFromId(nId)!=IV

        for eId in edgeIds :
            eId = int(eId)
            assert g.edgeFromId(eId)!=IV

        findEdges = g.findEdges(edges)

        assert np.array_equal(findEdges,edgeIds)

    def testIters(self):
        g  = vigraph.listGraph()
        
        nodes = [n for n in g.nodeIter()]
        assert len(nodes)==0

        g.addNode(3)
        nodes = [n for n in g.nodeIter()]
        assert len(nodes)==1
        assert g.id(nodes[0]) == 3 

        g.addNode(6)
        nodes = [n for n in g.nodeIter()]
        assert len(nodes)==2
        assert g.id(nodes[0]) == 3 
        assert g.id(nodes[1]) == 6 

        g.addNode(2)
        nodes = [n for n in g.nodeIter()]
        assert len(nodes)==3
        assert g.id(nodes[0]) == 2 
        assert g.id(nodes[1]) == 3 
        assert g.id(nodes[2]) == 6 