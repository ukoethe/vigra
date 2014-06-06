import vigra
import numpy
import h5py
import sys
import random
from matplotlib import pyplot as plot

from PyQt4.QtGui import QApplication, QImage, QPainter, QPen, QPainterPath, QColor, QPolygonF
from PyQt4.QtCore import Qt, QRectF, QPointF, QPoint, QRect, QLineF
from qimage2ndarray import gray2qimage, array2qimage

def drawArrow(painter, start, offset, color=(255,0,0), sz=1):
    arrowH = 2*sz #height of arrow tip
    arrowW = 1*sz #width of arrow tip
    tailR  = sz/4.0 #radius of tail circle
    
    myPen = QPen(QColor(*color)) 
    myColor = QColor(*color)
    arrowSize = sz;
    painter.setPen(myPen);
    painter.setBrush(myColor);

    line = QLineF(start, start+offset)
    if line.length() == 0.0:
        return

    a1 = numpy.asarray( [line.p1().x(), line.p1().y()] ) #start 
    a2 = numpy.asarray( [line.p2().x(), line.p2().y()] ) #end
    b  = a2-a1
    B  = b/numpy.linalg.norm(b) #unit vector in arrow direction
    Bp = numpy.asarray( [-B[1], B[0]] ) #unit vector perpendicular to arrow direction
    #d,e: points defining arrow head
    c = a2 - arrowH*B
    d = c+arrowW/2.0*Bp
    e = c-arrowW/2.0*Bp
   
    #draw arrow line
    painter.drawLine(line)
    
    #draw arrow head
    poly = QPolygonF( [QPointF(d[0],d[1]),QPointF(e[0],e[1]),QPointF(a2[0], a2[1]), QPointF(d[0],d[1])] )
    painter.drawPolygon(poly)
   
    #draw arrow tail
    painter.drawEllipse(QRectF(line.p1().x()-tailR, line.p1().y()-tailR, 2*tailR, 2*tailR))


def drawLine(p, N, x,y, dx, dy):
    #r = N*random.random()/5
    r = 0
    drawArrow(p, QPointF(N*x+N/2.0+r, N*y+N/2.0+r), QPointF(N*dx, N*dy), sz=N/5.0)


def makePlotNew(a, img, fname=None):
    N =30 
    img = img.astype(numpy.uint8)
    if img.ndim == 2:
        qimg = gray2qimage(img.swapaxes(0,1))
    else:
        qimg = array2qimage(img.swapaxes(0,1))
    
    canvas = QImage(N*qimg.width(), N*qimg.height(), QImage.Format_ARGB32)
    p = QPainter(canvas)
    p.setRenderHints(QPainter.Antialiasing)
    p.drawImage(QRectF(0,0,N*qimg.width(), N*qimg.height()), qimg)

    p.setPen(QPen(Qt.blue))
    for i in range(a.shape[0]):
        p.drawLine(QPointF(N*i,0), QPointF(N*i, N*a.shape[1]))
    for j in range(a.shape[1]):
        p.drawLine(QPointF(0,N*j), QPointF(N*a.shape[0], N*j))

    p.setPen(QPen(Qt.red))
    #drawLine(0,1, 2,2)
    k = 0
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            I, J = float(a[i,j,0]), float(a[i,j,1])
            drawLine(p,N, i,j, I,J)
            k+=1
    del p
    canvas.save(fname)

######################################################################################################################

#seg = numpy.asarray([[1,2,1,1],[1,1,1,1],[1,1,1,2]], dtype=numpy.uint8)
seg = numpy.asarray([
    [7,8,7,7], 
    #[8,8,8,8],
    #[8,8,8,8]
    ],
    dtype=numpy.uint8).T
    
#seg = vigra.impex.readImage("seg.tiff")[400:415,400:420,0].astype(numpy.uint32).T
#vigra.impex.writeImage(seg, "seg_small.tiff")
#seg = vigra.impex.readImage("seg_small.tiff")[:,:,0].astype(numpy.uint32).T

#############

cei = vigra.analysis.regionImageToCrackEdgeImage(seg.astype(numpy.uint32))
img = numpy.zeros(cei.shape, numpy.uint8)
img[cei == 0] = 1
a = vigra.filters.vectorialDistanceTransform(img.astype(numpy.float64))
img = 255*img
makePlotNew(a, img, fname="/tmp/A.png")

#############

a = vigra.filters.vectorialBoundaryDistanceTransform(seg.astype(numpy.float64))

print seg

for i in range(a.shape[0]):
    for j in range(a.shape[1]):
        sys.stdout.write("(%1.2f, %1.2f) " % (a[i,j,0], a[i,j,1]))
    sys.stdout.write("\n")

relabel = (255*numpy.random.random((seg.max()+1, 3))).astype(numpy.uint8)
img = relabel[seg]
makePlotNew(a, img, fname="/tmp/B.png")


