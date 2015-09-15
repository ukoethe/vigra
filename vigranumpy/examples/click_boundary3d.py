import vigra
import vigra.graphs as vigraph
import vigra.graphs as graphs
import pylab
import numpy
np=numpy
import sys
import matplotlib
import pylab as plt
import math
from matplotlib.widgets import Slider, Button, RadioButtons


from libs.viewerModule import view3d




# parameter:
imPath = ('holyRegion.h5', 'im')   # input image path
labPath = ('segMaskOnly.h5', 'data')   # labeled image path

# load volume
labels = vigra.impex.readHDF5(*labPath).astype(np.uint32)
volume = vigra.impex.readHDF5(*imPath)

gridGraph = graphs.gridGraph(labels.shape)
rag = graphs.regionAdjacencyGraph(gridGraph, labels)

# view3d(grayData=((volume, 'vol'),), segData=((labels, 'lab'),))
# from IPython import embed; embed()

gui = vigra.graphs.TinyEdgeLabelGui(rag=rag, img=volume)
gui.startGui()
