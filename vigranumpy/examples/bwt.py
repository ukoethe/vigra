import vigra
from vigra import numpy


produceBug = False

# parameter
filepath = '12003.jpg'  # input image path
img = vigra.impex.readImage(filepath)



if produceBug:
	roiBegin = (img.shape[0]-100, img.shape[1]-100)
	roiEnd  = (img.shape[0] , img.shape[1])
else:
	roiBegin = (100, 100)
	roiEnd = (200, 200)



res = vigra.filters.gaussianSmoothing(img, 4.0)
resRoi = res.copy()
resRoiS = vigra.filters.gaussianSmoothing(img, 4.0, roi=(roiBegin,roiEnd))
resRoi[roiBegin[0]:roiEnd[0], roiBegin[1]:roiEnd[1]] = resRoiS

diff = numpy.abs(res[roiBegin[0]:roiEnd[0], roiBegin[1]:roiEnd[1]] - resRoi[roiBegin[0]:roiEnd[0], roiBegin[1]:roiEnd[1]])

if diff.sum()>0.00001:
    print diff.sum()/(diff.size*3)
    print "there is a difference between the roi end the non roi version"
    vigra.imshow(diff)
    vigra.show()
