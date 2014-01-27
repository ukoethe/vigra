import vigra
import vigra.graphs as vigraph
import pylab
import numpy



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



if False:
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
	labels ,nseg = vigra.analysis.slicSuperpixels(imgLab,10.0,2)
	labels 		 = numpy.squeeze(vigra.analysis.labelImage(labels))

	numLabels = labels.max()

	print "start hierarchicalSuperpixels"
	labels = vigraph.hierarchicalSuperpixels(labels=labels,edgeIndicatorImage=gradmag,nodeFeaturesImage=imgLab,
		nSuperpixels=40,verbose=True,beta=0.2,nodeDistType='squaredNorm',degree1Fac=1.0,wardness=1)#0.001)

	#showSeg(img,labels)
	#pylab.show()


if True:
	print "get input"
	f   = 'data.h5'
	img = vigra.impex.readHDF5(f,'data')[0:200,0:100,0:100].astype(numpy.float32)
	gradmag = vigra.filters.gaussianGradientMagnitude(img,1.0)
	labels,nlabels = vigra.analysis.slicSuperpixels(img,10.0,2)
	labels = vigra.analysis.labelVolume(labels)
	print "start hierarchicalSuperpixels"
	labels = vigraph.hierarchicalSuperpixels(labels=labels,edgeIndicatorImage=gradmag,nodeFeaturesImage=img,
		nSuperpixels=10,verbose=True,beta=0.2,nodeDistType='squaredNorm',degree1Fac=1.0,wardness=1)#0.001)
