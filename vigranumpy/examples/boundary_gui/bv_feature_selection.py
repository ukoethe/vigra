import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import numpy#
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType










class FeatureSelectionDialog(QtGui.QDialog):

    def __init__(self,viewer, parent):
        super(FeatureSelectionDialog, self).__init__(parent)
        
        self.resize(800,600)
        self.viewer = viewer
        self.layout = QtGui.QVBoxLayout()
        self.setLayout(self.layout)

        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        

        self.buttonBox.accepted.connect(self.onPressAccepted)


        def makeCheckBox(name, val=True):
            return {
                'name': name, 
                'type': 'bool', 
                'value': val, 
                #'tip': "This is a checkbox",
            }


        sigmaOpt = {'name': 'sigma', 'type': 'str', 'value': '[0.0, 1.0, 2.0, 4.0]' }
        wardOpts = {'name': 'wardness', 'type': 'str', 'value': '[0.0, 0.1, 0.2]' }

        filterChild = [
            makeCheckBox("computeFilter"),
            sigmaOpt,
            {
                'name':'UCM',
                'children': [
                    makeCheckBox("ucmFilters"),
                    wardOpts,
                    {'name': 'meanSign', 'type': 'float', 'value': '1.0' }
                ]
            }
        ] 

        params = [
            {
                'name' : "RawData",
                'type' : 'group',
                'children' : [
                    {
                        'name': 'Compute Features On Raw Data', 
                        'type': 'bool', 
                        'value': True, 
                        'tip': "This is a checkbox",
                    },
                    {
                        'name' : "0-Order Filter",
                        'type' : 'group',
                        'children' : filterChild
                    },
                    {
                        'name' : "1-Order Filter",
                        'type' : 'group',
                        'children' : filterChild
                    },
                    {
                        'name' : "2-Order Filter",
                        'type' : 'group',
                        'children' : filterChild
                    }
                ]
            },
            #ComplexParameter(name='Custom parameter group (reciprocal values)'),
            #ScalableGroup(name="Expandable Parameter Group", children=[
            #    {'name': 'ScalableParam 1', 'type': 'str', 'value': "default param 1"},
            #    {'name': 'ScalableParam 2', 'type': 'str', 'value': "default param 2"},
            #]),
        ]

        ## Create tree of Parameter objects
        self.p = Parameter.create(name='params', type='group', children=params)
        self.t = ParameterTree()
        self.t.setParameters(self.p, showTop=False)

        self.layout.addWidget(self.t)
        self.layout.addWidget(self.buttonBox)

        ## If anything changes in the tree, print a message
        def change(param, changes):
            print("tree changes:")
            for param, change, data in changes:
                path = self.p.childPath(param)
                if path is not None:
                    childName = '.'.join(path)
                else:
                    childName = param.name()
                print('  parameter: %s'% childName)
                print('  change:    %s'% change)
                print('  data:      %s'% str(data))
                print('  ----------')
            
        self.p.sigTreeStateChanged.connect(change)










    def onPressAccepted(self):
        self.hide()
        self.viewer.onClickedComputeFeaturesImpl(self.p)
    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Escape:
            self.hide()
            event.accept()
        else:
            super(QtGui.QDialog, self).keyPressEvent(event)
