import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys
import matplotlib
import pylab
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import math

def labelImageFromRagLabels(ragLabels,spLabels,shape,out):

    for x in range(shape[0]):
        for y in range(shape[1]):
            ragId    = spLabels[x,y]
            out[x,y] = ragLabels[ragId]
    out = vigra.analysis.labelImage(out)
    return out


def segImg(img,labels):
    crackedEdges = vigra.analysis.regionImageToCrackEdgeImage(labels)
    whereEdge    =  numpy.where(crackedEdges==0)
    imgToDisplay = vigra.resize(img,crackedEdges.shape)
    imgToDisplay-=imgToDisplay.min()
    imgToDisplay/=imgToDisplay.max()
    for c in range(img.ndim):
        ic = imgToDisplay[:,:,c]
        ic[whereEdge]=0

    return imgToDisplay




print "get input"
f       = '100075.jpg'
#f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 4.0
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab  = vigra.colors.transform_RGB2Lab(img)
newShape = [img.shape[0]*2-1,img.shape[1]*2-1]
#imgLab     = vigra.resize(imgLab,newShape)
#img    = vigra.resize(img,newShape)
gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
gradmag = numpy.squeeze(gradmag).astype(numpy.float32)
labels ,nseg = vigra.analysis.slicSuperpixels(imgLab,10.0,25)
labels       = numpy.squeeze(vigra.analysis.labelImage(labels))

numLabels = labels.max()
newLabel  = labels.copy()




gridGraph      = vigraph.gridGraph(img.shape[0:2])
rag,hyperEdges = vigraph.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)
edgeIndicator  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)
edgeIndicator  = vigraph.hyperEdgeImageFeatures(rag,gridGraph,hyperEdges,gradmag,edgeIndicator)


seeds = vigraph.graphMap(graph=rag,item='node',dtype=numpy.uint32,channels=1)


seeds[10]=1
seeds[20]=2
seeds[30]=3

visuImg = segImg(img,labels)
ax = plt.gca()
fig = plt.gcf()
implot = ax.imshow(numpy.swapaxes(visuImg,0,1))




#fig, ax = plt.subplots()
#plt.subplots_adjust(bottom=0.2)
#t = np.arange(0.0, 1.0, 0.001)
#s = np.sin(2*np.pi*freqs[0]*t)
#l, = plt.plot(t, s, lw=2)


class WsCallsback:
    ind = 0
    def __init__(self,seeds,labels,img):
        self.seeds = seeds
        #self.clearSeeds()
        self.labels=labels
        self.newLabel = self.labels.copy()
        self.img = img
        self.seedNr = 1
    def clearSeeds(self, event=None):
        print "clear segs"
        self.seeds[:]=0
        self.seedNr = 1

    def nextSeed(self, event):
        print "next seed "
        self.seedNr +=1
        print self.seedNr

    def segment(self, event):

        print "seeds", numpy.unique(self.seeds)

        ragLabels = vigra.graphs.watershedsSegmentation(rag,edgeIndicator,self.seeds)
        self.newLabel  = labelImageFromRagLabels(ragLabels,self.labels,img.shape,out=self.newLabel)
        implot.set_data(numpy.swapaxes(segImg(self.img,self.newLabel),0,1))
        ax.imshow(numpy.swapaxes(segImg(self.img,self.newLabel),0,1))
        plt.draw()

    def onclick(self,event):
        if event.xdata != None and event.ydata != None:
            
            if event.xdata>1 and event.ydata>1:
                x,y = long(math.floor(event.xdata)),long(math.floor(event.ydata))
                print "ONLICKRAW",event.xdata,event.ydata
                x/=2 
                y/=2
                print "ragId",x,y,self.labels[x,y]
                self.seeds[self.labels[x,y]]=self.seedNr
                #self.seedNr+=1

callback = WsCallsback(seeds,labels,img)
axclear = plt.axes([0.5, 0.05, 0.1, 0.075])
axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
axnext = plt.axes([0.81, 0.05, 0.1, 0.075])


clear = Button(axclear, 'Clear')
clear.on_clicked(callback.clearSeeds)

segment = Button(axprev, 'Segment')
segment.on_clicked(callback.segment)

nextSeed = Button(axnext, 'NextSeed')
nextSeed.on_clicked(callback.nextSeed)

cid = fig.canvas.mpl_connect('button_press_event', callback.onclick)

pylab.show()