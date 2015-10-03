import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import numpy#








class BvGridViewBox(pg.ViewBox):

    sigScrolled = QtCore.Signal(object)

    #                                   

    sigBlocksAppeared = QtCore.Signal(object)
    sigBlocksDisappeared = QtCore.Signal(object)




    def __init__(self, blocking2d):


        super(BvGridViewBox,self).__init__()
        self.setAspectLocked(True)
        self.setMenuEnabled(False)
        self.blocking2d = blocking2d
        self.blockVisibility = numpy.zeros(len(self.blocking2d),dtype='bool')
        self.visibleBlocks = None
        self.setAcceptHoverEvents(True)



        #proxy = pg.SignalProxy(self.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)

        self.sigXRangeChanged.connect(self.rangeChanged)
        self.sigYRangeChanged.connect(self.rangeChanged)




    def rangeChanged(self):
        #print '\n\n\n\n'
        rect =  self.viewRect()

        self.visibleBlocks =  self.blocking2d.iBlocks(rect)

        tmp = self.blockVisibility.copy()
        #print self.blockVisibility
        self.blockVisibility[:] = 0
        self.blockVisibility[self.visibleBlocks] = 1

        changingBlocks = numpy.where(tmp!=self.blockVisibility)[0]

        if(len(changingBlocks)>0):
            newState = self.blockVisibility[changingBlocks]
            self.sigBlocksAppeared.emit( changingBlocks[numpy.where(newState==1)] )
            self.sigBlocksDisappeared.emit( changingBlocks[numpy.where(newState==0)] )


        #print self.blockVisibility
    def wheelEvent(self, ev, axis=None):
        kmods = ev.modifiers()
        if kmods & pg.QtCore.Qt.ControlModifier:
            super(BvGridViewBox,self).wheelEvent(ev, axis)
        else:
            d = (ev.delta() * self.state['wheelScaleFactor'])
            self.sigScrolled.emit(d)


    def mouseDragEvent(self, ev, axis=None):

        kmods = ev.modifiers()
        if kmods & pg.QtCore.Qt.ControlModifier:
            super(BvGridViewBox,self).mouseDragEvent(ev, axis)
        else:
            ev.accept()
            #print self.mapSceneToView(ev.pos())

            view = self.scene().views()[0]
            tr = view.viewportTransform()
            point = ev.scenePos()
            items = self.scene().items(point, QtCore.Qt.IntersectsItemShape, QtCore.Qt.DescendingOrder, tr)
            

            for item in items:
                #print item
                if isinstance(item, BvPlotCurveItem):
                    if(item.mouseShape().contains(item.mapFromScene(point))):
                        item.viewer.edgeClicked(item, ev)
                        #item.mouseClickEvent(ev)










            #ev.accept()
class BvImageItem(pg.ImageItem):
    def __init__(self,**kwargs):
        super(BvImageItem,self).__init__(**kwargs)




class BvPlotCurveItem(pg.PlotCurveItem):
    def __init__(self,**kwargs):
        super(BvPlotCurveItem,self).__init__(**kwargs)
        self.bRect = QtCore.QRectF()

        self.setClickable(None, 8)

    def mouseClickEvent(self, ev):
        print "click event"
        if self.mouseShape().contains(ev.pos()):
            self.viewer.edgeClicked(self, ev)
        else:
            print "outside"
  
    def boundingRect(self):
        return self.bRect



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
    imPath = ('../holyRegion.h5', 'im')  #ut image path
    labPath = ('../segMaskOnly.h5', 'data')   # labeled image path



    # load volume
    labels = vigra.impex.readHDF5(*labPath).astype(numpy.uint32)[:,:,0:20]
    volume = vigra.impex.readHDF5(*imPath)[:,:,0:20]
    imgItem.setImage(volume[:,:,0])
    vb.addItem(imgItem)

    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
