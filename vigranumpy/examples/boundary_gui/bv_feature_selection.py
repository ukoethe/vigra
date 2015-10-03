import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import numpy#







class FeatureSelectionDialog(QtGui.QDialog):

    def __init__(self,viewer, parent):
        super(FeatureSelectionDialog, self).__init__(parent)
        #self.resize(300,200)
        self.viewer = viewer
        self.layout = QtGui.QHBoxLayout()
        self.setLayout(self.layout)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.layout.addWidget(self.buttonBox)

        self.buttonBox.accepted.connect(self.onPressAccepted)
    #def showEvent(self, event):
    #    geom = self.frameGeometry()
    #    geom.moveCenter(QtGui.QCursor.pos())
    #    self.setGeometry(geom)
    #    super(Dialog, self).showEvent(event)

    def onPressAccepted(self):
        self.hide()
        self.viewer.onClickedComputeFeaturesImpl()
    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.hide()
            event.accept()
        else:
            super(Dialog, self).keyPressEvent(event)
