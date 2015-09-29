import vigra
import vigra.graphs as vigraph
import vigra.graphs as graphs
import pylab
import numpy as np
import numpy
import sys
import matplotlib
import pylab as plt
import math
from matplotlib.widgets import Slider, Button, RadioButtons

from functools import partial


import matplotlib.lines as lines


from pyqtgraph.Qt import QtCore, QtGui
import numpy as np
import pyqtgraph as pg
import pyqtgraph.ptime as ptime



# parameter:
imPath = ('holyRegion.h5', 'im')   # input image path
labPath = ('segMaskOnly.h5', 'data')   # labeled image path

# load volume
labels = vigra.impex.readHDF5(*labPath).astype(np.uint32)[:,:,0:20]
volume = vigra.impex.readHDF5(*imPath)[:,:,0:20]

gridGraph = graphs.gridGraph(labels.shape)
rag = graphs.regionAdjacencyGraph(gridGraph, labels)








app = QtGui.QApplication([])





class DownCtrl(QtGui.QWidget):
    def __init__(self,*args,**kwargs):
        super(DownCtrl,self).__init__(*args,**kwargs)

        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)


        # MODE SELECTOR
        modeSelector = QtGui.QComboBox()
        modeSelector.addItem('Black')
        modeSelector.addItem('None')
        modeSelector.addItem('Features')
        modeSelector.addItem('Labels')
        modeSelector.addItem('Probabilities')
        self.layout.addWidget(modeSelector)
        
        # color change button
        self.colorChangeButton = pg.ColorButton('Edge Color')
        self.layout.addWidget(self.colorChangeButton)

        # Edge Size Slide
        self.brushSizeSlider = QtGui.QSlider(orientation=QtCore.Qt.Horizontal)
        self.layout.addWidget(self.brushSizeSlider)


        # Save load LABELS
        self.saveLabelsButton = QtGui.QPushButton('save Labels')
        self.layout.addWidget(self.saveLabelsButton)
        self.loadLabelsButton = QtGui.QPushButton('load Labels')
        self.layout.addWidget(self.loadLabelsButton)


class AllCurves(pg.GraphItem):
    def __init__(self):
        pg.GraphItem.__init__(self)
        self.curves = []

    def setCurves(self, curves):
        self.curves = curves
        for c in self.curves:
            c.setParentItem(self)


class EdgeGui(object):
    def __init__(self, rag):
        self.rag = rag
        self.labels = rag.labels
        self.nZ = self.labels.shape[2]
        self.shapeXY = self.labels.shape[0:2]
        self.sliceDict = None

        # qt gui
        self.win = QtGui.QMainWindow()
        self.win.setWindowTitle('pyqtgraph example: ImageItem')
        self.win.show()
        self.win.resize(800, 600)


        self.cw = QtGui.QWidget()
        self.win.setCentralWidget(self.cw)
        self.layout = QtGui.QGridLayout()
        self.cw.setLayout(self.layout)
        self.gv = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.gv)
        self.viewBox
        self.viewBox = self.gv.addViewBox()
        self.viewBox.setAspectLocked(True)


        self.downCtrlWidget = DownCtrl()
        self.layout.addWidget(self.downCtrlWidget)



        def sliderMoved(val):
            self.changeLineSize(val)
        self.downCtrlWidget.brushSizeSlider.sliderMoved.connect(sliderMoved)




        button = QtGui.QPushButton()
        self.layout.addWidget(button)
        def pressed():
            self.cz +=1
            #self.setZ(self.cz%self.nZ)
            if(self.cz % 2 ==0):
                for curve in self.curves:
                    curve.setPen(pg.mkPen({'color': (0,0,0,0), 'width': 3}))
            else:
                for curve in self.curves:
                    curve.setPen(pg.mkPen({'color': (0,0,0,255), 'width': 3}))
        button.clicked.connect(pressed)


        self.imgItem = pg.ImageItem(border='w')
        self.viewBox.addItem(self.imgItem)
        self.dataDict = dict()

        self.curves = []
        self.allCurves = AllCurves()
        self.cz = 0

    def setData(self, data, key):
        self.dataDict[key] = data


    def changeLineSize(self, size):
        for curve in self.curves:
            curve.setPen(pg.mkPen({'color': (0,0,1), 'width': size+1}))
        self.viewBox.update()


    def setZ(self, z):

        def handleClic(curve, edge):
            print curve,"edge",edge
            curve.setPen(pg.mkPen({'color': (0,0,1), 'width': 6}))

        pen = pg.mkPen({'color': "FF0", 'width': 6})
        vb = self.viewBox

        def vbAdd(item):
            if item.zValue() < vb.zValue():
                item.setZValue(vb.zValue()+1)

            #item.setParentItem(vb.childGroup)
            vb.addedItems.append(item)

        dataSlice = self.dataDict['raw'][:,:,z]
        labelSlice = self.labels[:,:,z]  

        self.imgItem.setImage(dataSlice)
        self.imgItem.update()
        #self.viewBox.update()


        #for curve in self.curves:
        #    self.viewBox.removeItem(curve)
        self.curves = []

        slicesEdges = vigra.graphs.SliceEdges(self.rag)


        with vigra.Timer("findSlicesEdges"):
            slicesEdges.findSlicesEdges(labelSlice)

        with vigra.Timer("render them"):
            
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
                curve = pg.PlotCurveItem(clickable=True,parent=vb.childGroup)
                curve.setPen(pen)
                curve.setData(lx,ly, connect="finite")
                curve.sigClicked.connect(partial(handleClic,edge=edge))
                    
                self.curves.append(curve)
                #with vigra.Timer("add curve"):
                #vbAdd(curve)
                #self.viewBox.update()
        self.allCurves.setCurves(self.curves)
        self.viewBox.addItem(self.allCurves)
        with vigra.Timer("update auto range"):
            vb.updateAutoRange()

    def show(self):
        self.win.show()


gui = EdgeGui(rag)
gui.setData(volume,'raw')


gui.show()
gui.setZ(10)


i=0
updateTime = ptime.time()
fps = 0

def updateData():
    global gui,i, updateTime, fps

    ## Display the data
    i = (i+1) % gui.nZ
    print i
    gui.setZ(i)
    print i

    QtCore.QTimer.singleShot(1, updateData)
    now = ptime.time()
    fps2 = 1.0 / (now-updateTime)
    updateTime = now
    fps = fps * 0.9 + fps2 * 0.1
    
    #print "%0.1f fps" % fps
    

#updateData()




## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
