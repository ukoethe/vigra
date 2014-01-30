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
from vigra import graphs as vigraph




def testGridGraph():
    pass



class TestGraph(object):

    def testAddNodesWithIds(self):

        IV = vigraph.INVALID

        g  = vigraph.graph()
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


        g = vigraph.graph()

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