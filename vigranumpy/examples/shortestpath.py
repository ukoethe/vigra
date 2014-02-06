import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys
import matplotlib
import pylab as plt
import math



print "get input"
f       = '100075.jpg'
#f       = '69015.jpg'
#f       = '12003.jpg'
sigma   = 2.0
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab  = vigra.colors.transform_RGB2Lab(img)
gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
gradmag = numpy.squeeze(gradmag)

print gradmag.shape,gradmag.dtype
    

if False:
    weights = gradmag.copy()
    weights = vigra.resize(weights,(100,200))
    weights = weights**2
    vigra.graphs.allPairsOfshortestImagePath(weights)

if True:
    weights = gradmag.copy()
    weights = numpy.exp(-5.0*weights)
    visuimg =img.copy()
    ax = plt.gca()
    fig = plt.gcf()
    visuimg-=visuimg.min()
    visuimg/=visuimg.max()
    implot = ax.imshow(numpy.swapaxes(visuimg,0,1))

    clickList=[]

    def onclick(event):
        global clickList
        global weights
        global img
        if event.xdata != None and event.ydata != None:
            x,y = long(math.floor(event.xdata)),long(math.floor(event.ydata))
            clickList.append((x,y))
            if len(clickList)==1:
                print "source ",clickList[0]
            elif len(clickList)==2:
                print "target ",clickList[1]

                path = vigra.graphs.shortestImagePathDijkstra(weights, clickList[0],clickList[1] )
                coords = (path[:,0],path[:,1])
                visuimg =img.copy()
                iR=visuimg[:,:,0]
                iG=visuimg[:,:,1]
                iB=visuimg[:,:,2]
                iR[coords]=255
                iG[coords]=0
                iB[coords]=0

                visuimg-=visuimg.min()
                visuimg/=visuimg.max()
                implot.set_data(numpy.swapaxes(visuimg,0,1))
                plt.draw()


                clickList=[]


    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    #vigra.imshow(visuimg)
    #pylab.show()

if True:

    weights = gradmag.copy()
    weights = weights**2

    img = vigra.graphs.shortestImagePathDistanceImageDijkstra(weights,(58,54))


    img-=img.min()
    img/=img.max()


    ax = plt.gca()
    fig = plt.gcf()
    implot = ax.imshow(numpy.swapaxes(img,0,1))

    fig2=plt.figure()

    gradmag-=gradmag.min()
    gradmag/=gradmag.max()

    ax2 = plt.gca()
    fig2 = plt.gcf()
    implot2 = plt.imshow(numpy.swapaxes(gradmag,0,1))

    img=None

    def onclick1(event):
        global weights
        global img
        if event.xdata != None and event.ydata != None:
            x,y = long(math.floor(event.xdata)),long(math.floor(event.ydata))

            # distance image
            img = vigra.graphs.shortestImagePathDistanceImageDijkstra(weights,(x,y))
            img-=img.min()
            img/=img.max()

            img =vigra.filters.nonlinearDiffusion(img,scale=20,edgeThreshold=0.004)

            implot.set_data(numpy.swapaxes(img,0,1))
            plt.draw()
            print x,y


    def onclick2(event):
        global weights
        global img
        if event.xdata != None and event.ydata != None:
            grad = vigra.filters.gaussianGradientMagnitude(img,2.0)
            grad-=grad.min()
            grad/=grad.max()
            implot2.set_data(numpy.swapaxes(grad,0,1))
            plt.draw()

    cid = fig.canvas.mpl_connect('button_press_event', onclick1)
    cid = fig2.canvas.mpl_connect('button_press_event', onclick2)
    plt.show()