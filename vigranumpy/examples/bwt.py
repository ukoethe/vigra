import vigra
from vigra import numpy


produceBug = True

# parameter
filepath = '12003.jpg'  # input image path
img = vigra.impex.readImage(filepath)


roiBegin = (100, 100)

if produceBug:
    roiEnd  = (img.shape[0] - 3, img.shape[1] - 3)
else:
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
