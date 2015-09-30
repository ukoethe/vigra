import vigra
import vigra.graphs as vigraph
import vigra.graphs as graphs
import pylab
import numpy as np
import numpy
import sys
import math
import matplotlib
import pylab as plt
import math
import h5py
from matplotlib.widgets import Slider, Button, RadioButtons

from functools import partial


import matplotlib.lines as lines


from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
import pyqtgraph as pg
import pyqtgraph.ptime as ptime

from bv_view_box import *

# parameter:
imPath = ('../holyRegion.h5', 'im')   # input image path
labPath = ('../segMaskOnly.h5', 'data')   # labeled image path


labPath=  ("/home/tbeier/Desktop/hhes/pmap_pipe/superpixels_10000.h5",'data')
imPath = ("/home/tbeier/Desktop/hhes/pmap_pipe/data_sub.h5",'data')


# load volume
labels = vigra.impex.readHDF5(*labPath).astype(np.uint32)
volume = vigra.impex.readHDF5(*imPath).T

print labels.shape, volume.shape

labels = labels[0:500,0:500,0:20]
volume = volume[0:500,0:500,0:20]

gridGraph = graphs.gridGraph(labels.shape)
rag = graphs.regionAdjacencyGraph(gridGraph, labels)








app = QtGui.QApplication([])





class DownCtrl(QtGui.QWidget):
    def __init__(self,*args,**kwargs):
        super(DownCtrl,self).__init__(*args,**kwargs)

        self.mainLayout = QtGui.QHBoxLayout()
        self.setLayout(self.mainLayout)


        # MODE SELECTOR
        self.sharedCtrlLayout = QtGui.QVBoxLayout()
        self.mainLayout.addLayout(self.sharedCtrlLayout)
        self.modeSelectorComboBox = QtGui.QComboBox()
        self.modeSelectorComboBox.addItem('LabelMode')
        self.modeSelectorComboBox.addItem('None')
        self.modeSelectorComboBox.addItem('Features')
        self.modeSelectorComboBox.addItem('Black')
        self.modeSelectorComboBox.addItem('Probabilities')
        self.sharedCtrlLayout.addWidget(self.modeSelectorComboBox)
        
        # color change button
        self.colorChangeButton = pg.ColorButton('Edge Color')
        self.sharedCtrlLayout.addWidget(self.colorChangeButton)

        # Edge Size Slide
        self.brushSizeSlider = QtGui.QSlider(orientation=QtCore.Qt.Horizontal)
        self.sharedCtrlLayout.addWidget(self.brushSizeSlider)


        #Save load LABELS
        self.saveLoadLabelsLayout = QtGui.QVBoxLayout()
        self.mainLayout.addLayout(self.saveLoadLabelsLayout)
        self.saveLabelsButton = QtGui.QPushButton('Save Labels')
        self.loadLabelsButton = QtGui.QPushButton('Load Labels')
        self.saveLoadLabelsLayout.addWidget(self.saveLabelsButton)
        self.saveLoadLabelsLayout.addWidget(self.loadLabelsButton)


        # Compute / load / save features?
        self.featureLayout = QtGui.QVBoxLayout()
        self.mainLayout.addLayout(self.featureLayout)
        self.computeFeaturesButton = QtGui.QPushButton('Save Features')
        self.saveFeaturesButton = QtGui.QPushButton('Save Features')
        self.loadFeaturesButton = QtGui.QPushButton('Load Features')
        self.featureLayout.addWidget(self.computeFeaturesButton)
        self.featureLayout.addWidget(self.saveFeaturesButton)
        self.featureLayout.addWidget(self.loadFeaturesButton)

    def mode(self):
        return self.modeSelectorComboBox.currentText()
    def edgeSize(self):
        return self.brushSizeSlider.value()

class AllCurves(pg.GraphItem):
    def __init__(self):
        pg.GraphItem.__init__(self)
        self.curves = []

    def setCurves(self, curves):
        self.curves = curves
        for c in self.curves:
            c.setParentItem(self)


class EdgeGui(object):
    def __init__(self, rag, ndim=3, axis=2):
        self.rag = rag
        self.labels = rag.labels
        self.shape = self.labels.shape
        self.ndim = ndim
        self.axis = axis
        assert len(self.shape) == ndim

        self.nSlices = 1
        if ndim == 3:
            self.nSlices = self.shape[axis]

        self.sliceDict = None

        # qt gui
        self.win = QtGui.QMainWindow()
        self.win.setWindowTitle('pyqtgraph example: ImageItem')
        self.win.show()
        self.win.resize(800, 600)


        self.cw = QtGui.QWidget()
        self.cw.setMouseTracking(True)
        self.win.setCentralWidget(self.cw)
        self.layout = QtGui.QGridLayout()
        self.cw.setLayout(self.layout)
        self.gv = pg.GraphicsLayoutWidget()
        self.gv.setMouseTracking(True)
        self.layout.addWidget(self.gv)
        self.viewBox = BvViewBox()
        #self.viewBox.setMouseTracking(True)
        self.gv.addItem(self.viewBox)
        self.viewBox.setAspectLocked(True)


        self.edgeClickLabels = dict()

        def scrolled(d):
            if d>0: 
                d=1
            else :
                d=-1
            if self.ndim == 3:
                newSlice = min(self.nSlices-1, self.currentSlice - d)
                newSlice = max(0, newSlice)
                self.currentSlice = newSlice
                self.setZ(self.currentSlice)
        self.viewBox.sigScrolled.connect(scrolled)






        self.ctrlWidget = DownCtrl()
        self.layout.addWidget(self.ctrlWidget)



        def sliderMoved(val):
            self.updatePens()
        self.ctrlWidget.brushSizeSlider.sliderMoved.connect(sliderMoved)


        self.ctrlWidget.saveLabelsButton.clicked.connect(self.onClickedSaveLabels)
        self.ctrlWidget.loadLabelsButton.clicked.connect(self.onClickedLoadLabels)


        self.imgItem = pg.ImageItem(border='w')
        self.viewBox.addItem(self.imgItem)
        self.dataDict = dict()

        self.curves = []
        self.allCurves = None
        self.currentSlice = 0

        self.pathHint = None



    def onClickedComputeFeatures(self):
        print "compute features"

    def onClickedComputeFeatures(self):
        print "compute features"
    def onClickedComputeFeatures(self):
        print "compute features"


    def onClickedSaveLabels(self):
        
        keys = self.edgeClickLabels.keys()
        if(len(keys) == 0):
            raise Exception("has no labels to save")
        vals = [self.edgeClickLabels[k] for k in keys]

        
        keys = numpy.array(keys,dtype='int64')
        vals = numpy.array(keys,dtype='int64')

        path = str(pg.QtGui.QFileDialog.getSaveFileName(caption='Save file',directory='/home'))
        f = h5py.File(path,'w')
        f['edgeIds'] = keys
        f['labels'] = vals
        f.close()
        
    def onClickedLoadLabels(self):
        path = str(pg.QtGui.QFileDialog.getOpenFileName(caption='Open file',directory='/home'))
        f = h5py.File(path,'r')
        keys = f['edgeIds'][:]
        vals = f['labels'][:]

        for k,v in zip(keys,vals):
            self.edgeClickLabels[k] = v

        self.updatePens()

    def mode(self):
        return self.ctrlWidget.mode()

    def setData(self, data, key):
        self.dataDict[key] = data


    def changeLineSize(self, size):
        for curve in self.curves:
            curve.setPen(pg.mkPen({'color': (0,0,1), 'width': size+1}))
        self.viewBox.update()


    def updatePens(self):
        if self.allCurves is not None:
            for curve in self.allCurves.curves:
                curve.setPen(self.getPen(curve.edge))

    def getPen(self, edge):
        w = self.edgeWidth()
        if self.mode() == 'LabelMode':
            if edge in self.edgeClickLabels:
                label = self.edgeClickLabels[edge]
                if label == 0 :
                    return pg.mkPen({'color': (255,0,0,50), 'width':w})
                else:
                    return pg.mkPen({'color': (0,255,0), 'width':w})
            else:
                return pg.mkPen({'color': (0,0,255), 'width':w})
        else:
            return pg.mkPen({'color': (0,0,1), 'width':w})

    def edgeWidth(self):
        return self.ctrlWidget.edgeSize()

    def edgeClicked(self,curve, ev):
        mode = self.mode()
        print curve.edge, "clicked",self.mode()

        if(mode == 'LabelMode'):
            if ev.button() == 1:
                self.edgeClickLabels[curve.edge] = 1
            elif ev.button() == 2:
                self.edgeClickLabels[curve.edge] = 0
            elif ev.button() == 4:
                if  curve.edge in self.edgeClickLabels:
                    self.edgeClickLabels.pop(curve.edge)
            curve.setPen(self.getPen(curve.edge))    


    def setZ(self, z):


            

        vb = self.viewBox

      

        dataSlice = self.dataDict['raw'][:,:,z]
        labelSlice = self.labels[:,:,z]  

        self.imgItem.setImage(dataSlice)
        self.imgItem.update()
        #self.viewBox.update()


        #for curve in self.curves:
        #    self.viewBox.removeItem(curve)
        self.curves = []

        slicesEdges = vigra.graphs.SliceEdges(self.rag)


        with vigra.Timer("find slices"):
            slicesEdges.findSlicesEdges(labelSlice)


        with vigra.Timer("build curves"):
            visibleEdges = slicesEdges.visibleEdges()
            for edge in visibleEdges:

                edge = long(edge)
                lineIds = slicesEdges.lineIds(edge)
                totalLine = []
                for lid in lineIds:
                    line = slicesEdges.line(long(lid)).astype('float32')/2.0
                    totalLine.append(line)
                    totalLine.append([[float('nan'),float('nan')]])
                totalLine = numpy.concatenate(totalLine,axis=0)
                lx = totalLine[:,0]
                ly = totalLine[:,1]
                #with vigra.Timer("get curve"):
                curve = BvPlotCurveItem(clickable=False,parent=vb.childGroup)
                curve.edge = edge
                curve.viewer = self
                curve.setPen(self.getPen(edge))
                curve.setData(lx,ly, connect="finite")
                #curve.sigClicked.connect(self.edgeClicked)
                self.curves.append(curve)

        with vigra.Timer("add"):
            if self.allCurves is not None:
                self.viewBox.removeItem(self.allCurves)
                self.allCurves = None
            self.allCurves = AllCurves()
      
            self.allCurves.setCurves(self.curves)
            self.viewBox.addItem(self.allCurves)
            vb.updateAutoRange()

    def show(self):
        self.win.show()


gui = EdgeGui(rag)
gui.setData(volume,'raw')


gui.show()
gui.setZ(0)




## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
