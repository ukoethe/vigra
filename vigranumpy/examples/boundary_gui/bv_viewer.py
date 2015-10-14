import vigra
import vigra.graphs as graphs
import vigra.blockwise as bw
import numpy
import sys
import math
import math
import h5py
from vigra import blockwise as bw
from functools import partial
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import threading
import opengm

from bv_view_box import *
from bv_layer    import *
from bv_feature_selection import *





imPath = ("/media/tbeier/data/datasets/hhess/data_sub.h5",'data')
imPath = ("/home/tbeier/Desktop/hhes/pmap_pipe/data_sub.h5",'data')
volume = vigra.impex.readHDF5(*imPath).astype('float32')
volume = volume[:,:,0:600]
volume = vigra.taggedView(volume,'xyz')

if True:

    if False:
        options = bw.BlockwiseConvolutionOptions3D()
        sigma = 2.7
        options.stdDev = (sigma, )*3
        options.blockShape = (64, )*3

        print "hessianEv"
        ev = bw.hessianOfGaussianFirstEigenvalue(volume, options)
        print ev.shape, ev.dtype

        options.stdDev = (4.5, )*3
        print "smooth"
        ev = bw.gaussianSmooth(ev, options)

        vigra.impex.writeHDF5(ev,"growmap.h5",'data')
    
    else:
        ev = vigra.impex.readHDF5("growmap.h5",'data')

    with vigra.Timer("watershedsNew"):
        labels, nseg = vigra.analysis.unionFindWatershed3D(ev,(100,100,100))
   

    print "gridGraph"
    gridGraph = graphs.gridGraph(labels.shape)
    rag = graphs.regionAdjacencyGraph(gridGraph, labels)
    rag.writeHDF5("rag.h5",'data')
else:
    rag = vigra.graphs.loadGridRagHDF5("rag.h5",'data')
    labels=rag.labels










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
        self.modeSelectorComboBox.addItem('McRes')
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


        self.featureGradientWidget = pg.GradientWidget(orientation='top')
    
        d = {'ticks': [(0.0, (165, 0, 60, 255)), (1.0, (0, 170, 60, 255))], 'mode': 'rgb'}

        self.featureGradientWidget.restoreState(d)

        self.featureGradientWidget.setEnabled(False)
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
        self.predictMcButton = QtGui.QPushButton('MC')


        self.saveRfButton = QtGui.QPushButton('Save Rf')
        self.loadRfButton = QtGui.QPushButton('Load Rf')
        
        self.rfLayout1.addWidget(self.trainRfButton)
        self.rfLayout1.addWidget(self.predictProbsButton)
        self.rfLayout1.addWidget(self.predictMcButton)

        self.rfLayout2.addWidget(self.saveRfButton)
        self.rfLayout2.addWidget(self.loadRfButton)
        


    def getColor(self, x, toQColor=True):
        return self.featureGradientWidget.item.getColor(x, toQColor)

    def setFeatures(self, nFeatures):
        self.featureGradientWidget.setEnabled(True)

    def mode(self):
        return self.modeSelectorComboBox.currentText()
    def edgeSize(self):
        return self.brushSizeSlider.value()

class AllCurves(pg.GraphItem):
    def __init__(self):
        pg.GraphItem.__init__(self)
        self.curves = []

    #def setCurves(self, curves):
    #    self.curves = curves
    #    self.curves[0].setParentItem(self)
    #    #for c in self.curves:
    #    #    c.setParentItem(self)




class Blocking2d(object):
    def __init__(self, shape, blockShape):
        self.blocking = vigra.blockwise.Blocking2D(shape, blockShape)

    def iBlocks(self, rect):
        tl = rect.topLeft()
        br = rect.bottomRight()

        begin = (  int(tl.x()-0.5), int(tl.y()-0.5))
        end = (  int(br.x()+0.5), int(br.y()+0.5))
        return self.blocking.intersectingBlocks(begin,end)

    def __len__(self):
        return len(self.blocking)
class EdgeGui(object):
    def __init__(self, rag, ndim=3, axis=2):
        self.rag = rag
        self.labels = rag.labels

        self.featureExtractor = graphs.gridRagFeatureExtractor(self.rag, self.labels)
        self.featureExtractor.labels = self.labels
        self.featureExtractor.graph  = self.rag



        self.shape = self.labels.shape
        self.blocking2d = Blocking2d(self.shape[0:2], (40,40))
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
        self.layout = QtGui.QHBoxLayout()
        self.layout2 = QtGui.QVBoxLayout()
        self.cw.setLayout(self.layout)
        self.gv = pg.GraphicsLayoutWidget()


        self.gv.scene()
    
        self.gv.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Expanding)
        self.gv.setMouseTracking(True)
        
        self.layerCtrl = BvLayerCtrl()
        self.layout.addLayout(self.layout2)
        self.layout.addWidget(self.layerCtrl)


        self.layout2.addWidget(self.gv,3)
        self.viewBox = BvGridViewBox(blocking2d=self.blocking2d)
        #self.viewBox.setMouseTracking(True)
        self.gv.addItem(self.viewBox)
        self.viewBox.setAspectLocked(True)


        self.edgeClickLabels = dict()

        def scrolled(d):
            if d>0: 
                d=5
            else :
                d=-5
            if self.ndim == 3:
                newSlice = min(self.nSlices-1, self.currentSlice - d)
                newSlice = max(0, newSlice)
                self.currentSlice = newSlice
                with vigra.Timer("scroll"):
                    self.setZ(self.currentSlice)
        self.viewBox.sigScrolled.connect(scrolled)






        self.ctrlWidget = DownCtrl()
        self.ctrlWidget.setSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum)

        self.layout2.insertStretch(1,0)
        self.layout2.addWidget(self.ctrlWidget)




        self.ctrlWidget.modeSelectorComboBox.currentIndexChanged.connect(self.updatePens)
        self.ctrlWidget.brushSizeSlider.sliderMoved.connect(self.changeLineSize)


        self.ctrlWidget.saveLabelsButton.clicked.connect(self.onClickedSaveLabels)
        self.ctrlWidget.loadLabelsButton.clicked.connect(self.onClickedLoadLabels)

        self.ctrlWidget.computeFeaturesButton.clicked.connect(self.onClickedComputeFeatures)
        self.ctrlWidget.saveFeaturesButton.clicked.connect(self.onClickedSaveFeatures)
        self.ctrlWidget.loadFeaturesButton.clicked.connect(self.onClickedLoadFeatures)




        self.ctrlWidget.featureGradientWidget.sigGradientChanged.connect(self.onGradientChanged)
        

        self.ctrlWidget.trainRfButton.clicked.connect(self.onClickedTrainRf)
        self.ctrlWidget.predictProbsButton.clicked.connect(self.onClickedPredictProbs)
        self.ctrlWidget.predictMcButton.clicked.connect(self.onClickedMulticut)
        self.ctrlWidget.saveRfButton.clicked.connect(self.onClickedSaveRf)
        self.ctrlWidget.loadRfButton.clicked.connect(self.onClickedLoadRf)


        self.layerCtrl.sigFeatureSelected.connect(self.onFeatureSelected)

        self.imgItem = pg.ImageItem(border='w')
        self.viewBox.addItem(self.imgItem)
        self.dataDict = dict()

        self.curves = []
        self.allCurves = None
        self.currentSlice = 0

        self.pathHint = None


        self.currentFeatures = None
        self.featureMinMax = None
        self.currentFi = 0

        self.rf = None
        self.probs = None
        self.nCurves = 0


        self.ScrollLock = threading.Lock()
        self.scrollingIsLocked = False


    def onFeatureSelected(self, fIndex):
        self.currentFi = fIndex
        self.updatePens()

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


        #extractor =  self.featureExtractor
        #wardness = numpy.array([0.0, 0.1, 0.15, 0.25, 0.5, 1.0], dtype='float32')
        #ucmProbs = extractor.ucmTransformFeatures(self.probs[:,None],wardness)
        #newFeat = numpy.concatenate([self.currentFeatures, ucmProbs ],axis=1)
        #trainingInstances = numpy.array(self.edgeClickLabels.keys(),dtype='uint64')
        #labels = numpy.array([self.edgeClickLabels[k] for k in trainingInstances],dtype='uint32')[:,None]
        #features = newFeat[trainingInstances,:]
        #rf2 = vigra.learning.RandomForest(treeCount=255)
        #oob = rf2.learnRF(features, labels)
        #self.probs = rf2.predictProbabilities(newFeat)[:,1]
        #print "predict probs done"
        self.ctrlWidget.modeSelectorComboBox.setCurrentIndex(4)

    def onClickedMulticut(self):

        p1 = self.probs.copy()
        p1 = numpy.clip(p1, 0.005, 1.0-0.005)
        p0 = 1.0 - p1

        weights = numpy.log(p0/p1)
        nVar = self.rag.maxNodeId + 1
        nos = numpy.ones(nVar)*nVar
        gm = opengm.gm(nos)

        uv = self.rag.uvIds()
        uv = numpy.sort(uv,axis=1)
        pf = opengm.pottsFunctions([nVar,nVar], numpy.array([0]),weights)
        fid = gm.addFunctions(pf)
        gm.addFactors(fid,uv)

        pparam = opengm.InfParam(seedFraction=0.02)
        parameter = opengm.InfParam(generator='randomizedWatershed',proposalParam=pparam,numStopIt=20,numIt=3000)
        inf = opengm.inference.IntersectionBased(gm, parameter=parameter)
        inf.infer(inf.verboseVisitor())
        arg = inf.arg()

        self.eArg = arg[uv[:,0]]!=arg[uv[:,1]]

        self.ctrlWidget.modeSelectorComboBox.setCurrentIndex(6)

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




        d = FeatureSelectionDialog(viewer=self,parent=self.ctrlWidget)
        d.show()

    def onClickedComputeFeaturesImpl(self, param):





        lCtrl = self.layerCtrl

        extractor =  self.featureExtractor
        rawData = self.dataDict['raw']

        print "compute features"

        print param


        options = bw.BlockwiseConvolutionOptions3D()
        options.blockShape = (128, )*3


        fSize = [0] 
        def fRange(feat):
            old = long(fSize[0])
            fSize[0] += feat.shape[1]
            return (old, fSize[0])


        geoFeat,geoFeatNames = extractor.geometricFeatures()
        lCtrl.addFeature("GeometricFeature",fRange(geoFeat),subNames=geoFeatNames)

        topoFeat, topoFeatNames = extractor.topologicalFeatures()
        lCtrl.addFeature("TopologicalFeatures",fRange(topoFeat),subNames=topoFeatNames)

        features = [geoFeat, topoFeat]



        for order in range(3):
            orderName = '%d-Order Filter'%order
            doFilter = param[('RawData',orderName,'computeFilter')]
            if(doFilter):
                sigmas = param[('RawData',orderName,'sigma')]
                sigmas = list(eval(sigmas))

                doUcm = param[('RawData',orderName,'UCM','ucmFilters')]
                ucmMeanSign = float(param[('RawData',orderName,'UCM','meanSign')])
                print "ucmMeanSign",ucmMeanSign
                print "doUcm",doUcm
                wardness = param[('RawData',orderName,'UCM','wardness')]
                wardness = list(eval(wardness))
                wardness = numpy.array(wardness, dtype='float32')
                wnames =[]
                for w in wardness :
                    wnames.append(" w"+str(w))
                    wnames.append(" r"+str(w))

                print "sigma",sigmas, "w",wardness


                for sigma in sigmas:
                    print sigma
                    orderSigmaName = orderName + str(sigma)
                    options.stdDev = (float(sigma), )*3
                    if order == 0:
                        if sigma < 0.01:
                            res = rawData
                        else:
                            res = bw.gaussianSmooth(rawData, options)
                            accFeat, accFeatNames = extractor.accumulatedFeatures(res)
                            lCtrl.addFeature(orderSigmaName,fRange(accFeat),subNames=accFeatNames)
                            features.append(accFeat)

                            if doUcm :
                                mean = accFeat[:,0]
                                mean *= ucmMeanSign
                                ucmName = orderSigmaName+"MeanUcm"
                                ucm = extractor.ucmTransformFeatures(mean[:,None], wardness) 
                                lCtrl.addFeature(ucmName,fRange(ucm),subNames=wnames)
                                features.append(ucm)

                    if order == 1 :
                        if sigma < 0.5:
                            continue
                        else:
                            res = bw.gaussianGradientMagnitude(rawData, options)
                            accFeat, accFeatNames = extractor.accumulatedFeatures(res)
                            lCtrl.addFeature(orderSigmaName,fRange(accFeat),subNames=accFeatNames)
                            features.append(accFeat)

                            if doUcm :
                                mean = accFeat[:,0]
                                print "ucmMeanSign",ucmMeanSign
                                mean *= ucmMeanSign
                                ucmName = orderSigmaName+"MeanUcm"
                                ucm = extractor.ucmTransformFeatures(mean[:,None], wardness) 
                                lCtrl.addFeature(ucmName,fRange(ucm),subNames=wnames)
                                features.append(ucm)

                    if order == 2 :
                        if sigma < 0.5:
                            continue
                        else:
                            res = bw.hessianOfGaussianFirstEigenvalue(rawData, options)
                            accFeat, accFeatNames = extractor.accumulatedFeatures(res)
                            lCtrl.addFeature(orderSigmaName,fRange(accFeat),subNames=accFeatNames)
                            features.append(accFeat)

                            if doUcm :
                                mean = accFeat[:,0]
                                mean *= ucmMeanSign
                                ucmName = orderSigmaName+"MeanUcm"
                                ucm = extractor.ucmTransformFeatures(mean[:,None], wardness) 
                                lCtrl.addFeature(ucmName,fRange(ucm),subNames=wnames)
                                features.append(ucm)



        if False:
            for s in [2.0]:
                options.stdDev = (s, )*3
                res = bw.gaussianSmooth(rawData, options)
                accFeat, accFeatNames = extractor.accumulatedFeatures(res)
                lCtrl.addFeature("RawGaussianSmooth",fRange(accFeat),subNames=accFeatNames)
                features.append(accFeat)
            
            # 6
            wardness = numpy.array([0.0, 0.1, 0.15, 0.25, 0.5, 1.0], dtype='float32')
            wnames =[]
            for w in wardness :
                wnames.append(" w"+str(w))
                wnames.append(" r"+str(w))
            for s in [2.0]:
                options.stdDev = (s, )*3
                res = bw.hessianOfGaussianFirstEigenvalue(rawData, options)
                accFeat, accFeatNames = extractor.accumulatedFeatures(res)
                lCtrl.addFeature("hessianOfGaussianFirstEigenvalue",fRange(accFeat),subNames=accFeatNames)
                features.append(accFeat)
            
                print "ucm"
                mean = accFeat[:,0]
                ucm = extractor.ucmTransformFeatures(mean[:,None], wardness)
                print  "ucm done"
                lCtrl.addFeature("hessianMeanUcm",fRange(ucm),subNames=wnames)
                features.append(ucm)
            
        #for s in [1.0,  3.0,  4.0]:
        #    img = hessianEv2(rawData, s, 2.0)
        #    #vigra.imshow(img)
        #    #vigra.show()
        #    accFeat = extractor.accumulatedFeatures(img)
        #    features.append(accFeat)
        #    mean = accFeat[:,0]
        #    ucm = extractor.ucmTransformFeatures(mean[:,None],wardness)
        #    features.extend([accFeat,ucm])


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
        path = str(pg.QtGui.QFileDialog.getSaveFileName(caption='Save file',directory='/home/tbeier'))
        f = h5py.File(path,'w')
        f['features'] = self.currentFeatures
        f.close()
        
    def onClickedLoadFeatures(self):
        print "load features"
        path = str(pg.QtGui.QFileDialog.getOpenFileName(caption='Open file',directory='/home/tbeier'))
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
            vals = numpy.array(vals,dtype='int64')

            path = str(pg.QtGui.QFileDialog.getSaveFileName(caption='Save file',directory='/home/tbeier'))
            f = h5py.File(path,'w')
            f['edgeIds'] = keys
            f['labels'] = vals
            f.close()
        
    def onClickedLoadLabels(self):
        path = str(pg.QtGui.QFileDialog.getOpenFileName(caption='Open file',directory='/home/tbeier'))
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
        
        self.updatePens()

    def updatePens(self):
        if self.allCurves is not None:
            for curve in self.allCurves.curves[0:self.nCurves]:
                curve.setPen(self.getPen(curve.edge))
        self.viewBox.update()

    def getPen(self, edge):
        w = self.edgeWidth()
        m = self.mode() 
        if m == 'LabelMode':
            if edge in self.edgeClickLabels:
                label = self.edgeClickLabels[edge]
                if label == 0 :
                    return pg.mkPen({'color': (255,0,0,255), 'width':w})
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
                    return pg.mkPen({'color': (255,0,0,255), 'width':w})
                else:
                    return pg.mkPen({'color': (0,255,0), 'width':w})
            else:
                #assert self.probs is not None
                color = self.ctrlWidget.getColor(self.probs[edge])
                return pg.mkPen({'color': color, 'width':w})
        elif m == "McRes":
            state = self.eArg[edge]
            if(state == 0):
                return pg.mkPen({'color': (255,0,0,50), 'width':w})
            else:
                return pg.mkPen({'color': (0,255,0), 'width':w})

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

        
        
        slicesEdges = vigra.graphs.SliceEdges(self.rag)


        with vigra.Timer("find slices"):
            slicesEdges.findSlicesEdges(labelSlice)



        if self.allCurves is  None:
            self.allCurves = AllCurves()
            #self.viewBox.addItem(self.allCurves)
        
        self.curves = self.allCurves.curves
            

        ci = 0

        #with vigra.Timer("build curves"):
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
            

            
            try:
                curve = self.curves[ci]
            except:
                curve = BvPlotCurveItem(clickable=False,parent=self.imgItem)
                self.curves.append(curve)
                curve.viewer = self
                #self.viewBox.addItem(curve,ignoreBounds=True)
            ci += 1

            leftTop = numpy.nanmin(lx),numpy.nanmin(ly)
            rightBottom =  numpy.nanmax(lx),numpy.nanmax(ly)
            curve.setToolTip("id%d"%edge)
            curve.bRect.setCoords(leftTop[0],leftTop[1],rightBottom[0],rightBottom[1])
            curve.edge = edge
            curve.setPen(self.getPen(edge))
            curve.setData(lx,ly, connect="finite")
            curve.setVisible(True)

        self.nCurves = ci
        print self.nCurves
        for cii in range(ci, len(self.curves)):
            self.curves[cii].setVisible(False)



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
