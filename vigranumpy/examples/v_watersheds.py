import vigra
from vigra import numpy
import vigra.graphs as graphs
import pylab


# parameter:
filepath = '100075.jpg'   # input image path
filepath = '/home/tbeier/Desktop/HOLES-x13250-y5800/iso.007.png'
sigmaGradMag = 2.0       # sigma Gaussian gradient
superpixelDiameter = 10  # super-pixel size
slicWeight = 10.0        # SLIC color - spatial weight

# load image and convert to LAB
img = vigra.impex.readImage(filepath)

# get super-pixels with slic on LAB image
imgLab = vigra.colors.transform_RGB2Lab(img)

# compute gradient
imgLabBig = vigra.resize(imgLab, [imgLab.shape[0]*2-1, imgLab.shape[1]*2-1])
gradMag    = vigra.filters.gaussianGradientMagnitude(imgLab, sigmaGradMag)
gradMagBig = vigra.filters.gaussianGradientMagnitude(imgLabBig, sigmaGradMag*2.0)


print "get graphs"

# get 2D grid graph and  edgeMap for grid graph
# from gradMag of interpolated image
gridGraph = graphs.gridGraph(img.shape[0:2])
gridGraphEdgeIndicator = graphs.edgeFeaturesFromInterpolatedImage(gridGraph,
                                                                  gradMagBig)

print "get seg seeds"

# generate seeds
seeds = graphs.nodeWeightedWatershedsSeeds(gridGraph, gradMag)


# make viscosity term
C = 10.0
sigma = 2.0

viscTerm = 1.0-numpy.exp(-1.0*C*gradMagBig) 
viscTerm = vigra.filters.gaussianSmoothing(viscTerm, sigma)
viscTerm = numpy.log(1.0 + viscTerm)*(1.0/C)
print viscTerm.shape
vigra.imshow(viscTerm)
vigra.show()



labelsNodeWeighted  = graphs.nodeWeightedWatersheds(gridGraph, gradMag, seeds)
labelsEdgeNodeWeighted  = graphs.shortestPathSegmentation(gridGraph, gridGraphEdgeIndicator, viscTerm, seeds)
labelsEdgeWeighted  = graphs.edgeWeightedWatersheds(gridGraph, gridGraphEdgeIndicator, seeds)


# get region adjacency graph from super-pixel labels
ragEN = graphs.regionAdjacencyGraph(gridGraph, labelsEdgeNodeWeighted)
ragN = graphs.regionAdjacencyGraph(gridGraph, labelsNodeWeighted)
ragE = graphs.regionAdjacencyGraph(gridGraph, labelsEdgeWeighted)

f = pylab.figure()

ax0 = f.add_subplot(1, 3, 0)
ragEN.show(img)
ax0.set_title("edge node weighted")

ax1 = f.add_subplot(1, 3, 1)
ragN.show(img)
ax1.set_title("node weighted")

ax2 = f.add_subplot(1, 3, 2)
ragE.show(img)
ax2.set_title("edge weighted")

pylab.show()