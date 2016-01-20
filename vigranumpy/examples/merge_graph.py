import vigra
from vigra import graphs
from vigra import numpy
import pylab
# parameter
filepath = '12003.jpg'  # input image path
sigmaGradMag = 2.0      # sigma Gaussian gradient
superpixelDiameter = 10 # super-pixel size
slicWeight = 10.0       # SLIC color - spatial weight
beta = 0.5              # node vs edge weight
nodeNumStop = 50        # desired num. nodes in result


# load image and convert to LAB
img = vigra.impex.readImage(filepath)

# get super-pixels with slic on LAB image
imgLab = vigra.colors.transform_RGB2Lab(img)
labels, nseg = vigra.analysis.slicSuperpixels(imgLab, slicWeight,
                                              superpixelDiameter)
labels = vigra.analysis.labelImage(labels)


gridGraph = graphs.gridGraph(img.shape[0:2])
rag = graphs.regionAdjacencyGraph(gridGraph, labels)

# get the merge graph
mg = graphs.mergeGraph(rag)






# do n runs where we erase k edges in each run

n = 3
k = 200
for r in range(n):

    erased = 0 

    while(erased<k):

        # get a random edge
        randEdgeId = numpy.random.randint(rag.edgeNum)

        print "random edge:",randEdgeId
        # edge could be gone 
        # -since we could have merged it already
        # - or due to transitivity of other merges
        if mg.hasEdgeId(randEdgeId):
            mg.contractEdge(mg.edgeFromId(randEdgeId))
            erased+=1

    # get the segmentation
    # (brand new c++ function for max speed)
    labels = mg.graphLabels()

    # view results
    rag.show(img,labels)
    vigra.show()

    # get the result as pixels wise labeling
    asImage = rag.projectLabelsToGridGraph(labels)
    asImage = vigra.taggedView(asImage, "xy")
    
    
