import vigra
import vigra.graphs as vigraph
import vigra.graphs as graphs
import pylab
import numpy as np
import numpy
import sys
import math
import matplotlib
import pylab as plt
import math
import h5py
from matplotlib.widgets import Slider, Button, RadioButtons
from vigra import blockwise as bw
from functools import partial
import matplotlib.lines as lines
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
import pyqtgraph as pg
import pyqtgraph.ptime as ptime
import threading
import opengm

from bv_view_box import *
from bv_layer    import *
from bv_feature_selection import *


import zmq
import random
import sys
import time


    
if False:
    ragB = vigra.graphs.loadGridRagHDF5("ragb.h5",'data')
    labelsB=ragB.labels

    labels = labelsB[0:500,0:500,0:500]
    gridGraph = graphs.gridGraph(labels.shape)
    rag = graphs.regionAdjacencyGraph(gridGraph, labels)
    rag.writeHDF5("rags.h5",'data')
elif False:
    rag = vigra.graphs.loadGridRagHDF5("rag.h5",'data')
    labels = rag.labels
    raw = vigra.impex.readHDF5("/home/tbeier/Desktop/hhes/pmap_pipe/data_sub.h5","data")
    raw = raw[0:500,0:500,0:500].astype('float32')



    extractor = graphs.gridRagFeatureExtractor(rag, labels)


    print "geo feat"
    geoFeat, geoFeatNames = extractor.geometricFeatures()

    print "top feat"
    topoFeat, topoFeatNames = extractor.topologicalFeatures()

    print "acc "
    accFeat, accFeatNames = extractor.accumulatedFeatures(raw)

    features = [geoFeat, topoFeat, accFeat]
    features = numpy.concatenate(features,axis=1)
    vigra.impex.writeHDF5(features, "features.h5", "data")

else:

    if False:
        rag = vigra.graphs.loadGridRagHDF5("ragb.h5",'data')
        labels=rag.labels[0:500,0:500,0:500]
        gridGraph = graphs.gridGraph(labels.shape)
        rag2 = graphs.regionAdjacencyGraph(gridGraph, labels)
        rag2.writeHDF5("rag_saalfeld.h5",'data')
    else:
        
        rag = vigra.graphs.loadGridRagHDF5("rag_saalfeld.h5",'data')

        raw = vigra.impex.readHDF5("/home/tbeier/Desktop/hhes/pmap_pipe/data_sub.h5","data")
        raw = raw[0:500,0:500,0:500].astype('float32')

        labels=rag.labels[0:500,0:500,0:500]
        extractor = graphs.gridRagFeatureExtractor(rag, labels)
        
        if True:
            geoFeat, geoFeatNames = extractor.geometricFeatures()
            topoFeat, topoFeatNames = extractor.topologicalFeatures()
            accFeat, accFeatNames = extractor.accumulatedFeatures(raw)

            features = [geoFeat, topoFeat, accFeat]
            names = graphs.StringVector()
            names.extend(geoFeatNames)
            names.extend(topoFeatNames)
            names.extend(accFeatNames)
    
            features = numpy.concatenate(features,axis=1)
            vigra.impex.writeHDF5(features, "features2.h5", "data")
            #else:
            #features = vigra.impex.readHDF5("features2.h5", "data")
            print  features
            for fi in range(features.shape[1]):
                feat = features[:,fi]
                print  "mima",names[fi],feat.min(), feat.max(), feat

            features = numpy.nan_to_num(features)


        labels = rag.labels

        
        print features.shape

        userLabels = dict()

        port = "5556"
        context = zmq.Context()
        socket = context.socket(zmq.PAIR)
        socket.bind("tcp://*:%s" % port)



        print "RAG STUFF",rag.maxNodeId, rag.maxEdgeId

        arg = None
        def predict():
            global features,rag, sp,arg

            trainingInstances = numpy.array(userLabels.keys(),dtype='uint64')
            labels = numpy.array([userLabels[k] for k in trainingInstances],dtype='uint32')[:,None]
            
            trainingInstances = numpy.clip(trainingInstances, 0, rag.maxEdgeId-2)

            assert trainingInstances.max() <= rag.maxEdgeId

            if len(labels)>10 and labels.min()==0 and labels.max()==1: 
                feat = features[trainingInstances,:]
                rf = vigra.learning.RandomForest(treeCount=255)
                oob = rf.learnRF(feat, labels)
                print "oob", oob
                probs = rf.predictProbabilities(features)[:,1]
                p1 = probs.copy()
                p1 = numpy.clip(p1, 0.005, 1.0-0.005)
                p0 = 1.0 - p1

                weights = numpy.log(p0/p1)
                nVar = rag.maxNodeId + 1
                nos = numpy.ones(nVar)*nVar
                gm = opengm.gm(nos)


                uv = rag.uvIds()

                if weights.shape[0] < uv.shape[1]:
                    diff  = uv.shape[1] - weights.shape[0]
                    val = numpy.zeros(diff)
                    weights = numpy.concatenate([weights,val])

                
                uv = numpy.sort(uv,axis=1)
                pf = opengm.pottsFunctions([nVar,nVar], numpy.array([0]),weights)
                fid = gm.addFunctions(pf)
                gm.addFactors(fid,uv)

                param = opengm.InfParam(planar=False)
                inf = opengm.inference.Cgc(gm,parameter=param)
                if arg is not None:
                    inf.setStartingPoint(arg)
                visitor = inf.timingVisitor(timeLimit=60.0)
                inf.infer(visitor)
                arg = inf.arg()
                eArg = arg[uv[:,0]]!=arg[uv[:,1]]


                return arg
                
            else:
                return None



        print "MAX RAG ID",rag.maxNodeId

        while True:
            msg = socket.recv()

            print "THE MESSAGE FROM SAALFELD\n",msg
            if msg.startswith("merge"):
                isMerge = True
            else:
                isMerge = False



            if isMerge:
                msg = msg.replace("merge","")
                msg = msg.replace("(","")
                msg = msg.replace(")","")
                print "clean msg",msg
                mergeList = eval(msg)

                for a in range(len(mergeList)):
                    for b in range(a+1, len(mergeList)):

                        u = mergeList[a]
                        v = mergeList[b]

                        u = int(u)
                        v = int(v)
                        s = int(0)
                        
                        #print  "u",u,"v",v

                        edge = rag.findEdge(u,v)
                        if edge.id == -1:
                            pass
                            #print "air edges",edge.id
                            #un = rag.nodeFromId(u)
                            #vn = rag.nodeFromId(v)
                            #print "un vn",un.id, vn.id
                            ##edge = rag.addEdge(un,vn)
                            #print "new edge",edge.id
                        else:
                            userLabels[rag.id(edge)] = s
            else:
                msg = msg.replace("detach","")
                msg = msg.replace("(","")
                msg = msg.replace(")","")
                splitList = eval(msg)
                toSplitAway = splitList[0][0]
                rest = splitList[1]

                for r in rest:
                    if r != toSplitAway:

                        edge = rag.findEdge(r,toSplitAway)
                        if edge.id == -1:
                            pass
                        else:
                            userLabels[rag.id(edge)] = 1


            res = predict()
            if res is None:
                socket.send("NEED MORE")
            else:
                print "make the string"
                strRes = "["
                for v in res:
                    strRes +=str(v)+" "
                strRes+="]"
                print "send it"
                socket.send(strRes)
                print "done send"

            time.sleep(1)

