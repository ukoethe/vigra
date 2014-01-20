from vigra import *



def testGridGraph():


    INVALID = graphs.Invalid()

    g = graphs.GridGraphUndirected2d([10,10])

    print g.edgeIds()

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