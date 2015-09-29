import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np

class BvViewBox(pg.ViewBox):

    sigScrolled = QtCore.Signal(object)


    def __init__(self):
        super(BvViewBox,self).__init__()
        self.setAspectLocked(True)

    #def mouseClickEvent(self, ev):
    #    print "vb click"

    def wheelEvent(self, ev, axis=None):
        kmods = ev.modifiers()
        if kmods & pg.QtCore.Qt.ControlModifier:
            super(BvViewBox,self).wheelEvent(ev, axis)
        else:
            d = (ev.delta() * self.state['wheelScaleFactor'])
            self.sigScrolled.emit(d)

    def mouseDragEvent(self, ev, axis=None):
        kmods = ev.modifiers()
        if kmods & pg.QtCore.Qt.ShiftModifier:
            print "shift"
        if kmods & pg.QtCore.Qt.ControlModifier:
            super(BvViewBox,self).mouseDragEvent(ev, axis)

class BvImageItem(pg.ImageItem):
    def __init__(self,**kwargs):
        super(BvImageItem,self).__init__(**kwargs)

    def mouseClickEvent(self, ev):
        print "image click"



class BvPlotCurveItem(pg.PlotCurveItem):
    def __init__(self,**kwargs):
        super(BvPlotCurveItem,self).__init__(**kwargs)

        self.setAcceptHoverEvents(True)
        self.setAcceptDrops(True)
        self.setHandlesChildEvents(True)
        #self.acceptDrags(True)
 
    def mouseOverEvent(self, ev):
        print "mouse over"

    def mouseDragEvent(self, ev):
        print "drag click"

    def mouseMoveEvent(self, ev):
        print "move"

    def mouseHoverEvent(self, ev):
        print "mouse move"

    def mouseClickEvent(self, ev):
        print "clicked"
        self.setPen(pg.mkPen({'color': (150,150,50), 'width': 6}))

    def mouseReleaseEvent(self, ev):
        print "released"

    #def hoverMoveEvent(self, ev):
    #    print "fooo"

    def hoverEnterEvent(self, ev):
        #if ev.button() == QtCore.Qt.RightButton:
        #print "enter",ev
        if self.isUnderMouse():
            if self.contains(ev.pos()):
                self.setPen(pg.mkPen({'color': (255,0,1), 'width': 6}))
            else:
                self.setPen(pg.mkPen({'color': (0,0,255), 'width': 6}))
        else:
            self.setPen(pg.mkPen({'color': (0,255,0), 'width': 6}))
    #def hoverLeaveEvent(self, ev):
    #    self.setPen(pg.mkPen({'color': (0,0,1), 'width': 3}))

    #def shape(self):
    #    print "shape"
    #    return self.path

    #def boundingRect(self):
    #    print "JOO"

if __name__ == '__main__':

    import vigra
    import numpy

    app = QtGui.QApplication([])
    mw = QtGui.QMainWindow()    
    mw.show()
    mw.resize(800, 600)

    gv = pg.GraphicsView()
    mw.setCentralWidget(gv)
    l = QtGui.QGraphicsGridLayout()
    l.setHorizontalSpacing(0)
    l.setVerticalSpacing(0)
    vb = BvViewBox()

    l.addItem(vb, 0, 1)
    gv.centralWidget.setLayout(l)


    imgItem = BvImageItem(border='w')



    # parameter:
    imPath = ('../holyRegion.h5', 'im')   # input image path
    labPath = ('../segMaskOnly.h5', 'data')   # labeled image path



    # load volume
    labels = vigra.impex.readHDF5(*labPath).astype(numpy.uint32)[:,:,0:20]
    volume = vigra.impex.readHDF5(*imPath)[:,:,0:20]
    imgItem.setImage(volume[:,:,0])
    vb.addItem(imgItem)

    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()