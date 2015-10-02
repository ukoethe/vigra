import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np



class BvLayerCtrl(QtGui.QWidget):


    sigFeatureSelected = QtCore.Signal(object)

    def __init__(self):
        super(BvLayerCtrl, self).__init__()

        self.layout = QtGui.QHBoxLayout()
        self.setLayout(self.layout)




        
        self.edgeLayerTree  = pg.TreeWidget(parent=self)
        self.edgeLayerTree.setColumnCount(1)
        self.edgeLayerTree.setAcceptDrops(False)
        self.edgeLayerTree.setDragEnabled(False)

        self.edgeLayerTree.setColumnCount(1)

        self.labeslItem  = QtGui.QTreeWidgetItem(["Labels"])
        self.featuresItem  = QtGui.QTreeWidgetItem(["Features"])
        self.predictionsItem  = QtGui.QTreeWidgetItem(["Predictions"])

        self.edgeLayerTree.addTopLevelItem(self.labeslItem)
        self.edgeLayerTree.addTopLevelItem(self.featuresItem)
        self.edgeLayerTree.addTopLevelItem(self.predictionsItem)


        self.layout.addWidget(self.edgeLayerTree)


        self.edgeLayerTree.itemClicked.connect(self.onItemActivated)
        self.edgeLayerTree.itemActivated.connect(self.onItemActivated)
        self.edgeLayerTree.itemEntered.connect(self.onItemActivated)

    def addFeature(self, iname, indexRange, subNames ):
        
        item = QtGui.QTreeWidgetItem([iname])
        self.featuresItem.addChild(item)
        for i,sn in  enumerate(subNames):
            sitem = QtGui.QTreeWidgetItem([sn])
            sitem.featureIndex = i+indexRange[0]
            item.addChild(sitem)


    def onItemActivated(self, item):
        if item.childCount() ==0:
            print "item",item,"activated"
            fIndex = None
            try:
                fIndex = item.featureIndex
            except:
                fIndex = None

            if fIndex is not None:
                self.sigFeatureSelected.emit(fIndex)
