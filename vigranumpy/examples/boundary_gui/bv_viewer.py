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
if True:
    imPath = ('../holyRegion.h5', 'im')   # input image path
    labPath = ('../segMaskOnly.h5', 'data')   # labeled image path
    labels = vigra.impex.readHDF5(*labPath).astype(np.uint32)
    volume = vigra.impex.readHDF5(*imPath).astype('float32')
else:
    labPath=  ("/home/tbeier/Desktop/hhes/pmap_pipe/superpixels_10000.h5",'data')
    imPath = ("/home/tbeier/Desktop/hhes/pmap_pipe/data_sub.h5",'data')


    # load volume
    labels = vigra.impex.readHDF5(*labPath).astype(np.uint32)
    volume = vigra.impex.readHDF5(*imPath).astype('float32').T
    labels = labels[0:500,0:500,0:20]
    volume = volume[0:500,0:500,0:20]

print labels.shape, volume.shape



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
        self.modeSelectorComboBox.addItem('LabelProbabilities')
        self.sharedCtrlLayout.addWidget(self.modeSelectorComboBox)
        
      

        # Edge Size Slide
        self.brushSizeSlider = QtGui.QSlider(orientation=QtCore.Qt.Horizontal)
        self.brushSizeSlider.setValue(2)
        self.sharedCtrlLayout.addWidget(self.brushSizeSlider)


        #Save load LABELS
        self.saveLoadLabelsLayout = QtGui.QVBoxLayout()
        self.mainLayout.addLayout(self.saveLoadLabelsLayout)
        self.saveLabelsButton = QtGui.QPushButton('Save Labels')
        self.loadLabelsButton = QtGui.QPushButton('Load Labels')
        self.saveLoadLabelsLayout.addWidget(self.saveLabelsButton)
        self.saveLoadLabelsLayout.addWidget(self.loadLabelsButton)


        # Compute / load / save features?
        self.featureGetLayout = QtGui.QVBoxLayout()
        self.mainLayout.addLayout(self.featureGetLayout)
        self.computeFeaturesButton = QtGui.QPushButton('Compute Features')
        self.saveFeaturesButton = QtGui.QPushButton('Save Features')
        self.loadFeaturesButton = QtGui.QPushButton('Load Features')
        self.featureGetLayout.addWidget(self.computeFeaturesButton)
        self.featureGetLayout.addWidget(self.saveFeaturesButton)
        self.featureGetLayout.addWidget(self.loadFeaturesButton)


        # show the features
        self.featureShowLayout = QtGui.QVBoxLayout()
        self.mainLayout.addLayout(self.featureShowLayout)
        self.featureSpinBox = QtGui.QSpinBox()
        self.featureSpinBox.setEnabled(False)

        self.featureGradientWidget = pg.GradientWidget(orientation='top')
        self.featureGradientWidget.setEnabled(False)
        self.featureShowLayout.addWidget(self.featureSpinBox)
        self.featureShowLayout.addWidget(self.featureGradientWidget)


        # train rf 
        self.rfLayout = QtGui.QHBoxLayout()
        self.mainLayout.addLayout(self.rfLayout)
        self.rfLayout1 = QtGui.QVBoxLayout()
        self.rfLayout2 = QtGui.QVBoxLayout()
        self.rfLayout.addLayout(self.rfLayout1)
        self.rfLayout.addLayout(self.rfLayout2)

        self.trainRfButton = QtGui.QPushButton('TrainRf')
        self.predictProbsButton = QtGui.QPushButton('Predict Probs')

        self.saveRfButton = QtGui.QPushButton('Save Rf')
        self.loadRfButton = QtGui.QPushButton('Load Rf')
        
        self.rfLayout1.addWidget(self.trainRfButton)
        self.rfLayout1.addWidget(self.predictProbsButton)

        self.rfLayout2.addWidget(self.saveRfButton)
        self.rfLayout2.addWidget(self.loadRfButton)
        


    def getColor(self, x, toQColor=True):
        return self.featureGradientWidget.item.getColor(x, toQColor)

    def setFeatures(self, nFeatures):
        self.featureSpinBox.setMinimum(0)
        self.featureSpinBox.setMaximum(nFeatures-1)
        self.featureSpinBox.setEnabled(True)
        self.featureGradientWidget.setEnabled(True)

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

        self.featureExtractor = graphs.gridRagFeatureExtractor(self.rag, self.labels)
        self.featureExtractor.labels = self.labels
        self.featureExtractor.graph  = self.rag

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
        self.layout = QtGui.QVBoxLayout()
        self.cw.setLayout(self.layout)
        self.gv = pg.GraphicsLayoutWidget()
        self.gv.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Expanding)
        self.gv.setMouseTracking(True)
        
        self.layout.addWidget(self.gv,3)
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
        self.ctrlWidget.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum)

        self.layout.insertStretch(1,0)
        self.layout.addWidget(self.ctrlWidget)




        self.ctrlWidget.modeSelectorComboBox.currentIndexChanged.connect(self.updatePens)
        self.ctrlWidget.brushSizeSlider.sliderMoved.connect(self.updatePens)


        self.ctrlWidget.saveLabelsButton.clicked.connect(self.onClickedSaveLabels)
        self.ctrlWidget.loadLabelsButton.clicked.connect(self.onClickedLoadLabels)

        self.ctrlWidget.computeFeaturesButton.clicked.connect(self.onClickedComputeFeatures)
        self.ctrlWidget.saveFeaturesButton.clicked.connect(self.onClickedSaveFeatures)
        self.ctrlWidget.loadFeaturesButton.clicked.connect(self.onClickedLoadFeatures)



        self.ctrlWidget.featureSpinBox.valueChanged.connect(self.onFeatureToShowChanged)
        self.ctrlWidget.featureGradientWidget.sigGradientChanged.connect(self.onGradientChanged)
        

        self.ctrlWidget.trainRfButton.clicked.connect(self.onClickedTrainRf)
        self.ctrlWidget.predictProbsButton.clicked.connect(self.onClickedPredictProbs)
        self.ctrlWidget.saveRfButton.clicked.connect(self.onClickedSaveRf)
        self.ctrlWidget.loadRfButton.clicked.connect(self.onClickedLoadRf)


        self.imgItem = pg.ImageItem(border='w')
        self.viewBox.addItem(self.imgItem)
        self.dataDict = dict()

        self.curves = []
        self.allCurves = None
        self.currentSlice = 0

        self.pathHint = None


        self.currentFeatures = None
        self.featureMinMax = None
        self.currentFi = self.ctrlWidget.featureSpinBox.value()

        self.rf = None
        self.probs = None

    def onClickedTrainRf(self):
        assert self.currentFeatures is not None
        trainingInstances = numpy.array(self.edgeClickLabels.keys(),dtype='uint64')
        labels = numpy.array([self.edgeClickLabels[k] for k in trainingInstances],dtype='uint32')[:,None]
        features = self.currentFeatures[trainingInstances,:]

        self.rf = vigra.learning.RandomForest(treeCount=255)
        oob = self.rf.learnRF(features, labels)
        print "oob error ",oob
        self.onClickedPredictProbs()

    def onClickedPredictProbs(self):
        assert self.currentFeatures is not None
        assert self.rf is not None

        self.probs = self.rf.predictProbabilities(self.currentFeatures)[:,1]
        print "predict probs done"
        self.ctrlWidget.modeSelectorComboBox.setCurrentIndex(4)

    def onClickedSaveRf(self):
        print "save rf"

    def onClickedLoadRf(self):
        print "load rf"



    def onModeChanged(self):
        self.updatePens()

    def onGradientChanged(self):
        if self.mode() == 'Features' or self.mode() == 'Probabilities' or self.mode() == 'LabelProbabilities':
            self.updatePens()

    def onFeatureToShowChanged(self, value):
        self.currentFi = value
        if(self.mode()=='Features'):
            self.updatePens()

    def onClickedComputeFeatures(self):
        extractor =  self.featureExtractor
        rawData = self.dataDict['raw']

        print "compute features"
        geoFeat = extractor.geometricFeatures()
        topoFeat = extractor.topologicalFeatures()
        features = [geoFeat, topoFeat]

        for s in [2.0, 3.0, 4.0]:
            res = vigra.gaussianSmoothing(rawData, s)
            accFeat = extractor.accumulatedFeatures(res)
            features.append(accFeat)
        

        features = numpy.concatenate(features,axis=1)
        self.currentFeatures = features
        print "compute done"

        fMin = numpy.min(features, axis=0)
        fMax = numpy.max(features, axis=0)
        self.featureMinMax = (fMin, fMax)

        self.ctrlWidget.setFeatures(self.currentFeatures.shape[1])
        self.ctrlWidget.modeSelectorComboBox.setCurrentIndex(2)
    def onClickedSaveFeatures(self):
        if self.currentFeatures is None:
            raise RuntimeError("Has no features to save")
        path = str(pg.QtGui.QFileDialog.getSaveFileName(caption='Save file',directory='/home'))
        f = h5py.File(path,'w')
        f['features'] = self.currentFeatures
        f.close()
        
    def onClickedLoadFeatures(self):
        print "load features"
        path = str(pg.QtGui.QFileDialog.getOpenFileName(caption='Open file',directory='/home'))
        f = h5py.File(path,'r')
        self.currentFeatures = f['features'][:]
        fMin = numpy.min(self.currentFeatures, axis=0)
        fMax = numpy.max(self.currentFeatures, axis=0)
        self.featureMinMax = (fMin, fMax)

        self.ctrlWidget.setFeatures(self.currentFeatures.shape[1])
        self.ctrlWidget.modeSelectorComboBox.setCurrentIndex(2)
        #self.updatePens()

    def onClickedSaveLabels(self):
        
        keys = self.edgeClickLabels.keys()
        if(len(keys) == 0):
            #errorMsg = QErrorMessage()
            #errorMsg.showMsg("has no labels to save")
            raise RuntimeError("has no labels to save")
        else:
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
        m = self.mode() 
        if m == 'LabelMode':
            if edge in self.edgeClickLabels:
                label = self.edgeClickLabels[edge]
                if label == 0 :
                    return pg.mkPen({'color': (255,0,0,50), 'width':w})
                else:
                    return pg.mkPen({'color': (0,255,0), 'width':w})
            else:
                return pg.mkPen({'color': (0,0,255), 'width':w})

        elif m == "Features":
            assert self.currentFeatures is not None
            fi = self.currentFi
            fVal = self.currentFeatures[edge, fi]
            fVal -= self.featureMinMax[0][fi]
            fVal /= (self.featureMinMax[1][fi]-self.featureMinMax[0][fi])
            #fVal *= 255.0

            color = self.ctrlWidget.getColor(fVal)
            return pg.mkPen({'color': color, 'width':w})
        elif m == "Black":
            return pg.mkPen({'color': (0,0,1), 'width':w})
        elif m == "Probabilities":
            #assert self.probs is not None
            color = self.ctrlWidget.getColor(self.probs[edge])
            return pg.mkPen({'color': color, 'width':w})
        elif m == "LabelProbabilities":

            if edge in self.edgeClickLabels:
                label = self.edgeClickLabels[edge]
                if label == 0 :
                    return pg.mkPen({'color': (255,0,0,50), 'width':w})
                else:
                    return pg.mkPen({'color': (0,255,0), 'width':w})
            else:
                #assert self.probs is not None
                color = self.ctrlWidget.getColor(self.probs[edge])
                return pg.mkPen({'color': color, 'width':w})

    def edgeWidth(self):
        return self.ctrlWidget.edgeSize()

    def edgeClicked(self,curve, ev):
        mode = self.mode()
        print curve.edge, "clicked",self.mode()

        if(mode == 'LabelMode' or mode == 'LabelProbabilities') :
            if ev.button() == 1:
                self.edgeClickLabels[curve.edge] = 1
            elif ev.button() == 2:
                self.edgeClickLabels[curve.edge] = 0
            elif ev.button() == 4:
                if  curve.edge in self.edgeClickLabels:
                    self.edgeClickLabels.pop(curve.edge)
            curve.setPen(self.getPen(curve.edge))    

        ev.accept()

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
                curve = BvPlotCurveItem(clickable=True,parent=vb.childGroup)
                curve.edge = edge
                curve.viewer = self
                curve.setPen(self.getPen(edge))
                curve.setData(lx,ly, connect="finite")
                curve.sigClicked.connect(self.edgeClicked)
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
