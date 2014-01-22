from vigra import *
import numpy






class NumpyArrayEdgeMap(numpy.ndarray):
    def __new__(cls, graph,dtype):
        #print 'In __new__ with class %s' % cls
        return numpy.ndarray.__new__(cls, shape=[graph.maxEdgeId+1],dtype=dtype)

    def __init__(self, *args, **kwargs):
        # in practice you probably will not need or want an __init__
        # method for your subclass
        #print 'In __init__ with class %s' % self.__class__
        pass

    def __array_finalize__(self, obj):
        pass
        #print 'In array_finalize:'
        #print '   self type is %s' % type(self)
        #print '   obj type is %s' %  type(obj)

    def __getitem__(self,key):
        return super(NumpyArrayEdgeMap, self).__getitem__(key.id)

    def __setitem__(self,key,value):
        super(NumpyArrayEdgeMap, self).__setitem__(key.id,value)



def graphMap(graph,item,dtype):
    if item=='edge':
        return NumpyArrayEdgeMap(graph,dtype)





def testGridGraph():


    INVALID = graphs.Invalid()

    g = graphs.GridGraphUndirected2d([10,10])
    n1 = g.nodeFromId(1)

    edgeMap = graphMap(g,'edge',numpy.float32)
    print edgeMap[n1]
    edgeMap[n1]=1.1
    print edgeMap[n1]


    returnedMap = g.maptest(edgeMap)


    print "retMap",returnedMap[n1]

    
    print "returnedMapType",type(returnedMap)


    #print g.edgeIds()

    edge = graphs.EdgeGridGraphUndirected2d()



    assert edge==INVALID

    n1 = g.nodeFromId(1)
    assert n1.id==1
    n2 = g.nodeFromId(2)
    assert n2.id==2
    n99 = g.nodeFromId(99)
    assert n99.id==99
    n999 = g.nodeFromId(999)
    assert n1!=INVALID
    assert n2!=INVALID
    assert n99!=INVALID
    #assert n999==INVALID
    assert INVALID!=n1
    assert INVALID!=n2
    assert INVALID!=n99
    print "bug for ulli"#assert INVALID==n999

    e12 = g.findEdge(n1,n2)
    assert e12!=INVALID
    assert INVALID!=e12

    e199 = g.findEdge(n1,n99)
    assert e199==INVALID
    assert INVALID==e199



    



testGridGraph()