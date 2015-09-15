import vigra
import vigra.graphs as vigraph
import vigra.graphs as graphs
import pylab
import numpy
np=numpy
import sys
import matplotlib
import pylab as plt
import math
from matplotlib.widgets import Slider, Button, RadioButtons






# parameter:
filepath = '100075.jpg'   # input image path
sigmaGradMag = 3.0       # sigma Gaussian gradient
superpixelDiameter = 20  # super-pixel size
slicWeight = 10.0        # SLIC color - spatial weight

# load image and convert to LAB
img = vigra.impex.readImage(filepath)

# get super-pixels with slic on LAB image
imgLab = vigra.colors.transform_RGB2Lab(img)
labels, nseg = vigra.analysis.slicSuperpixels(imgLab, slicWeight,
                                              superpixelDiameter)
labels = vigra.analysis.labelImage(labels)
gridGraph = graphs.gridGraph(img.shape[0:2])
rag = graphs.regionAdjacencyGraph(gridGraph, labels)



gui = vigra.graphs.TinyEdgeLabelGui(rag=rag, img=img)
gui.startGui()





if False:


    ragEdgeLabels = numpy.zeros(rag.edgeNum, dtype='float32')



    print "bar"
    imgWithEdges = rag.showEdgeFeature(img, ragEdgeLabels,returnImg=True,  labelMode=True)
    print imgWithEdges.shape






    visuimg =imgWithEdges.copy()
    ax = plt.gca()
    fig = plt.gcf()
    visuimg-=visuimg.min()
    visuimg/=visuimg.max()
    implot = ax.imshow(numpy.swapaxes(visuimg,0,1))




    def onclick(event):
        global img
        global labels
        global ragEdgeLabels
        shape = img.shape
        if event.xdata != None and event.ydata != None:
            xRaw,yRaw = event.xdata,event.ydata
            if xRaw >=0.0 and yRaw>=0.0 and xRaw<img.shape[0] and yRaw<img.shape[1]:
                x,y = long(math.floor(event.xdata)),long(math.floor(event.ydata))
                  
                print "X,Y",x,y
                l = labels[x,y]


                other  = None

                for xo in [-2,-1,0,1,2]:
                    for yo in [-2,-1,0,1,2]:
                        xx = x+xo
                        yy = y+yo
                        if xo is not 0 or yo is not 0:
                            if  xx >=0 and xx<shape[0] and \
                                yy >=0 and yy<shape[0]:
                                otherLabel = labels[xx, yy]
                                if l != otherLabel:
                                    edge = rag.findEdge(long(l), long(otherLabel))
                                    print edge
                                    other = (xx,yy,edge)
                                    break
                    if other is not None:
                        break
                
                if other is not None:
                    eid = other[2].id
                    oldLabel  = ragEdgeLabels[eid]

                    newLabel = 0 
                    if oldLabel == 0: 
                        newLabel = 1
                    elif oldLabel ==1 :
                        newLabel = -1
                    elif oldLabel == -1:
                        newLabel = 0

                    print "old label",oldLabel
                    print "new label",newLabel

                    ragEdgeLabels[eid] = newLabel
                    imgWithEdges = rag.showEdgeFeature(img, ragEdgeLabels,returnImg=True, labelMode=True)
                    implot.set_data(numpy.swapaxes(imgWithEdges,0,1))
                    plt.draw()


    def onslide(event):
        global img,gradmag,weights,clickList,sgamma
        weights  = makeWeights(sgamma.val)
        print "onslide",clickList
        #if len(clickList)>=2:
        #    print "we have  path"
        #    source = gridGraph.coordinateToNode(clickList[0])
        #    target = gridGraph.coordinateToNode(clickList[1])
        #    path = pathFinder.run(weights, source,target).path(pathType='coordinates')
        #    visuimg = makeVisuImage(path,img)
        #    implot.set_data(numpy.swapaxes(visuimg,0,1))
        #    plt.draw()




    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
