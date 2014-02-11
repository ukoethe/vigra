import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import time


def showSeg(img,labels):
	crackedEdges = vigra.analysis.regionImageToCrackEdgeImage(labels)
	whereEdge    =  numpy.where(crackedEdges==0)
	imgToDisplay = vigra.resize(img,crackedEdges.shape)
	imgToDisplay-=imgToDisplay.min()
	imgToDisplay/=imgToDisplay.max()
	for c in range(img.ndim):
		ic = imgToDisplay[:,:,c]
		ic[whereEdge]=0

	vigra.imshow(imgToDisplay)



if True:
	print "get input"
	f       = '100075.jpg'
	f       = '69015.jpg'
	f 		= '12003.jpg'
	sigma   = 3.0
	img     = vigra.impex.readImage(f)#[0:100,0:100,:]
	imgLab  = vigra.colors.transform_RGB2Lab(img)
	newShape = [img.shape[0]*2-1,img.shape[1]*2-1]
	imgLab 	= vigra.resize(imgLab,newShape)
	img 	= vigra.resize(img,newShape)
	gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
	gradmag = numpy.squeeze(gradmag).astype(numpy.float32)
	labels ,nseg = vigra.analysis.slicSuperpixels(imgLab,10.0,10)
	labels 		 = numpy.squeeze(vigra.analysis.labelImage(labels))

	numLabels = labels.max()

	

	
	#print "get histogram"
	#hist = vigra.histogram.gaussianHistogram(img,
	#	minVals=[0.0,0.0,0.0],maxVals=[255.01,255.01,255.01],
	#	bins=30,sigma=7.0).reshape([img.shape[0],img.shape[1],-1])
	#hsum = numpy.sum(hist,axis=2)
	#hist/=hsum[:,:,numpy.newaxis]


	t0=time.time()
	print "start hierarchicalSuperpixels"
	labels = vigraph.hierarchicalSuperpixels(labels=labels,edgeIndicatorImage=gradmag,nodeFeaturesImage=img,
		nSuperpixels=5,verbose=True,beta=0.,nodeDistType='chiSquared',wardness=0.0)#0.001)
	t1=time.time()
	print "hcluster time",t1-t0
	#showSeg(img,labels)
	#pylab.show()


if False:
	print "get input"
	f   = 'data.h5'
	img = vigra.impex.readHDF5(f,'data')[0:200,0:100,0:100].astype(numpy.float32)
	gradmag = vigra.filters.gaussianGradientMagnitude(img,1.0)
	labels,nlabels = vigra.analysis.slicSuperpixels(img,10.0,2)
	labels = vigra.analysis.labelVolume(labels)
	print "start hierarchicalSuperpixels"
	labels = vigraph.hierarchicalSuperpixels(labels=labels,edgeIndicatorImage=gradmag,nodeFeaturesImage=img,
		nSuperpixels=10,verbose=True,beta=0.2,nodeDistType='squaredNorm',wardness=1)#0.001)
