from volumina.api import Viewer
from PyQt4.QtGui import QApplication, QColor
import h5py
import numpy
import sys
import vigra

app = QApplication([])

#f = h5py.File("seg3d_small.h5")
#f = h5py.File("cylinder/cylinder.h5")

f = h5py.File("vectorialDist.h5")
seg = f["seg"].value[:,:,:].astype(numpy.uint32)
print "seg.shape=",seg.shape

ew = f["ew"].value
curl = f["curl"].value
magnitude = f["magnitude"].value
div = f["div"].value
f.close()
print "magnitude.shape=",magnitude.shape
print "ew.shape=",ew.shape

x = ew[1,:,:,:]
x = (x-x.min())/(x.max()-x.min())
x = x > 0.2

y = seg*x

'''
f = h5py.File("maxDist.h5")
vmd = f["vectorialMaxDistance"].value
vmd = vmd[::2, ::2, ::2, :]/2.0
vmd = numpy.sqrt(numpy.sum(numpy.square(vmd), axis=3))
print "vmd.shape=",vmd.shape
f.close()

vmd_mag = magnitude/vmd
'''

def normalize(a):
    m,M = a.min(), a.max()
    a = (255*(a-m)/float(M-m)).astype(numpy.uint8)
    return a

def normalizeIndividually(a, seg, dtype=None):
    return (a-a.min())/(a.max()-a.min())

    assert a.shape == seg.shape, "%r <-> %r" % (a.shape, seg.shape)

    feat= vigra.analysis.extractRegionFeatures(a.astype(numpy.float32), seg,
        features=["Minimum", "Maximum"])

    m = feat["Minimum"]
    M = feat["Maximum"]

    m2 = m[seg]
    M2 = M[seg]

    M2m2 = M2-m2
    M2m2[ M2m2 < 1E-5 ] = 1.0

    res = (a-m2)/M2m2

    if dtype is None:
        return res
    elif dtype==numpy.uint8:
        return (255*res).astype(numpy.uint8)
    else:
        raise RuntimeError("unhandled dtype")

normalizedM = (magnitude-magnitude.min())/float(magnitude.max() - magnitude.min())

def makeTintedColortable(tintColor):
    ctable = [None]*256
    for i in range(256):
        ctable[i] = QColor(tintColor.red(), tintColor.green(), tintColor.blue(),i).rgba()
    return ctable

res = []

v = Viewer()
v.addRandomColorsLayer(seg, name="seg")

v.addClickableSegmentationLayer(seg[numpy.newaxis, :,:,:, numpy.newaxis], name="seg/click")

tintColor = [QColor(255,0,0), QColor(0,255,0), QColor(0,0,255)]
for i in range(3):
    a = ew[i,:,:,:]
    a = numpy.sqrt(numpy.abs(a))

    a = normalizeIndividually(a, seg);

    l = v.addGrayscaleLayer(normalizeIndividually(a, seg, dtype=numpy.uint8), name="ch=%d" % i)
    l.visible = False

    col = tintColor[i]    

    x = normalizeIndividually(a*normalizedM, seg)

    res.append(x)

    l = v.addColorTableLayer((255*x).astype(numpy.uint8),
        colortable=makeTintedColortable(col), name="ch=%d * magnitude" % i)
    l.visible = False

    if i == 1:
        col = tintColor[i]    
        l = v.addColorTableLayer(normalizeIndividually(a*res[i-1], seg, dtype=numpy.uint8),
            colortable=makeTintedColortable(col), name="ch=%d * ch=%d" % (i, i-1))
        l.visible = False


l = v.addGrayscaleLayer(normalizeIndividually(magnitude, seg, dtype=numpy.uint8), name="magnitude")
l.visible = False

'''
l = v.addGrayscaleLayer(normalizeIndividually(vmd, seg, dtype=numpy.uint8), name="vmd")
l.visible = False

l = v.addGrayscaleLayer(normalizeIndividually(vmd_mag, seg, dtype=numpy.uint8), name="vmd/mag")
l.visible = False
'''

l = v.addGrayscaleLayer(normalizeIndividually(numpy.sqrt(numpy.sum(numpy.square(curl), axis=3)), seg, dtype=numpy.uint8), name="curl")
l.visible = False

l = v.addGrayscaleLayer(normalizeIndividually(div, seg, dtype=numpy.uint8), name="div")
l.visible = False

print ew.shape
l = v.addGrayscaleLayer(normalize(ew[1,:,:,:]-ew[2,:,:,:]), name="dd")
l.visible = False

l = v.addRandomColorsLayer(y, name="yyyy")
l.visible = False

'''
diff = ew[1,:,:,:]-ew[2,:,:,:]
diff = (diff-diff.min())/(diff.max()-diff.min())
x = diff > 0.2
l = v.addGrayscaleLayer(255*x, name="dd x")
l.visible = False

x = ew[1,:,:,:]
x = (x-x.min())/(x.max()-x.min())
x = x > 0.2
l = v.addGrayscaleLayer(255*x, name="dd y")
l.visible = False
'''

v.showMaximized()
app.exec_()
