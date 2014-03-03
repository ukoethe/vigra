import vigra
import numpy
import vigra.graphs as vigraph
import matplotlib.pyplot as plt
import scipy.misc
import sys


gamma = 0.0001
percentage = 100
f = "figure_2.png"





## img: image segment with 0: inside, 1: outside
## distFunc: function applied after distance transform, must be one of "exponential", "linear", "inverse"
## showPathImage: if True, the image with distance transform and paths will be shown
## percentageOfPaths: percentage of computed paths
def eccentricity( img, distFunc = "exponential", showPathImage = False, percentageOfPaths = 100, imgSaveName = "" ):

    img = img.astype(numpy.uint8)
    dim = len(img.shape)

    ## Enlarge image by one pixel on each side
    bigShape = []
    for s in img.shape:
        bigShape.append(s + 2)
    bigImg = numpy.ones(bigShape)
    slices = []
    for i in range(dim):
        slices.append(slice(1, bigImg.shape[i]-1))
    bigImg[ slices ] = img
    inside = numpy.where(bigImg==0)
    outside = numpy.where(bigImg==1)


    ## Apply distanceTransform and modify (outside: high values, inside: low values)
    distImage = vigra.filters.distanceTransform2D(bigImg.astype(numpy.float32))
    imgp = distImage.copy()
    # if showPathImage:
    #     imgp = distImage.copy()
    if distFunc == "exponential":
        distImage = numpy.exp(distImage*-gamma)
    elif distFunc == "linear":
        maxDist = distImage.max()
        distImage = maxDist - distImage
    elif distFunc == "inverse":
        w = numpy.where(distImage!=0)
        distImage[w] = 1/distImage[w]
    else:
        print "wrong parameters for distFunc in eccentricity"

    ## Distance in the inside between two pixels is 1.0
    #distImage = bigImg.copy().astype(numpy.float32)
    #distImage[inside]=1.0

    ## Set the outside to a very high value
    distImage[outside]=1000.0

    ## Image copy to draw the paths in
    #imgp = distImage.copy()
    imgp[outside] = 100

    ## Get image graph and its path finder
    gridGraph = vigraph.gridGraph(bigImg.shape[0:dim],False)
    graphShape = []
    for s in distImage.shape:
        graphShape.append(s*2-1)
    edgeWeights = vigra.resize(distImage, graphShape, order=0)
    #edgeWeights = vigra.graphs.edgeFeaturesFromInterpolatedImage(gridGraph, edgeWeights)
    edgeWeights = vigra.graphs.edgeFeaturesFromInterpolatedImageCorrected(gridGraph,edgeWeights)
    pathFinder = vigraph.ShortestPathPathDijkstra(gridGraph)

    ## Find borders in img
    if dim == 2:
        bigLblImg = vigra.analysis.labelImage(bigImg.astype(numpy.uint8))
    else:
        bigLblImg = vigra.analysis.labelVolume(bigImg.astype(numpy.uint8))
    rag = vigraph.GridRegionAdjacencyGraph(gridGraph, bigLblImg)
    node = vigraph.GridRegionAdjacencyGraph.nodeFromId(rag, long( bigLblImg[inside][0] ))
    edges = vigraph._ragFindEdges(rag, gridGraph, rag.affiliatedEdges, bigLblImg, node)

    borderImg = numpy.zeros(bigImg.shape)
    for edge in edges:
        slices = []
        for d in range(dim):
            slices.append( slice(edge[d], edge[d]+1) )
        borderImg[slices] = 1

    ## End points for paths (all points on the border)
    targets = numpy.where(borderImg==1)
    nTargets = len(targets[0])

    ## Find the diameter (longest path)
    eccLength = numpy.empty(nTargets)
    eccLength.fill(-1)
    eccTargetPath = {}
    vpIndex = 0
    vpGraphIndex = []
    for d in range(dim):
        vpGraphIndex.append(int(targets[d][vpIndex]))
    vp = gridGraph.coordinateToNode(vpGraphIndex)
    visited = numpy.zeros(nTargets)
    while True:
        visited[vpIndex] = 1
        pathFinder.run(edgeWeights, vp)
        eccChanged = False
        for j in range(nTargets):
            targetIndex = []
            for d in range(dim):
                targetIndex.append(int(targets[d][j]))
            target = gridGraph.coordinateToNode(targetIndex)
            pathLength = pathFinder.distance(target)
            m = max(eccLength[j], pathLength)
            if m > eccLength[j]:
                eccChanged = True
                eccLength[j] = m
                eccTargetPath[j] = pathFinder.path(pathType='coordinates', target=target)
        vpIndex = numpy.argmax(eccLength)
        vpGraphIndex = []
        for d in range(dim):
            vpGraphIndex.append(int(targets[d][vpIndex]))
        vp = gridGraph.coordinateToNode(vpGraphIndex)
        if visited[vpIndex] or not eccChanged:
            break

    ## Find the length of the diameter (non-weighted)
    path = eccTargetPath[vpIndex]
    dMax = 0
    for k in range(1, len(path)):
        diffCount = 0
        for d in range(dim):
            if path[k][d] != path[k-1][d]:
                diffCount += 1
        dMax += numpy.sqrt(diffCount)

    ## Find the midpoint of the diameter
    dMax = dMax/2
    if len(path) == 0:
        path = numpy.empty( (1, 2), numpy.uint8 )
        path[0][0] = targets[0][0]
        path[0][1] = targets[1][0]
    p1 = path[0]
    d1 = 0
    for k in range(1, len(path)):
        p2 = path[k]
        d2 = d1 + numpy.linalg.norm(p2-p1)
        if d2 > dMax:
            if (abs(d2-dMax) < abs(d1-dMax)):
                p1 = p2
            break
        p1 = p2
        d1 = d2

    ## Compute eccentricity from center (p1) to all points on border
    sourceIndex = []
    for d in range(dim):
        sourceIndex.append(int(p1[d]))
    source = gridGraph.coordinateToNode(sourceIndex)
    pathFinder.run(edgeWeights, source)
    maxPathLength = 0
    for j in range(nTargets):
        targetIndex = []
        for d in range(dim):
            targetIndex.append(int(targets[d][j]))
        target = gridGraph.coordinateToNode(targetIndex)
        pathLength = pathFinder.distance(target)
        maxPathLength = max(maxPathLength, pathLength)
        imgp[targetIndex[0], targetIndex[1]] = 40*pathLength

    imgp[ path[:,0], path[:,1] ] = 12*maxPathLength
    imgp[sourceIndex[0], sourceIndex[1]] = 40*maxPathLength
    plt.figure(distFunc)
    plt.imshow(numpy.swapaxes(imgp, 1, 0), interpolation='none')
    if len(imgSaveName)>1:
         scipy.misc.imsave(imgSaveName, numpy.swapaxes(imgp, 1, 0))
    #plt.show()




loadImage = True
if loadImage:
    ## Read image from file
    img = vigra.impex.readImage(f)
    labels = numpy.squeeze(vigra.analysis.labelImage(img))
else:
    ## Random 3d image
    img = numpy.ones( (100, 100, 20) )
    img = numpy.random.rand(100, 100, 20)

    ## Compute slic superpixels
    labels ,nseg = vigra.analysis.slicSuperpixels(img.astype(numpy.float32),100.0,50)
    labels       = numpy.squeeze(vigra.analysis.labelVolume(labels))

## Dimension
dim = len(labels.shape)

## Compute bounding boxes
regionFeatures = vigra.analysis.extractRegionFeatures(img.astype(numpy.float32), labels)
upperLeftBBs = regionFeatures["Coord<Minimum>"]
lowerRightBBs = regionFeatures["Coord<Maximum>"]
nBoxes = len(upperLeftBBs)-1

## Get segment inside its bounding box
segments = []
nonEmptyBoxIndices = []
for i in range(nBoxes):
    slices = []
    for d in range(dim):
        slices.append(slice(upperLeftBBs[i+1][d], lowerRightBBs[i+1][d]))
    subImg = labels[slices].copy()
    where = numpy.where(subImg==i+1)
    if len(where[0]) > 0:
        subImg[where] = 0
        subImg[numpy.where(subImg!=0)] = 1
        segments.append(subImg)
        nonEmptyBoxIndices.append(i+1)

## Apply eccentricity transform
pathLengths = []
counter = 0
for seg in segments:
    saveName = "ecc_transform_from_center_with_distance_transform/ecc_seg_"+`counter`+".png"
    #saveName = ""
    pathLength = eccentricity(seg, distFunc="linear", showPathImage=False, percentageOfPaths=percentage, imgSaveName=saveName)
    pathLengths.append(pathLength)
    counter = counter+1

# ## Testimage: map longest path to color
# maxPath = 0
# for i in range(len(pathLengths)):
#     m = max(pathLengths[i])
#     if m > maxPath:
#         maxPath = m
# labelCopy = labels.copy()
# for i in range(len(pathLengths)):
#     val = max(pathLengths[i]) * 255.0/maxPath
#     j = nonEmptyBoxIndices[i]
#     labelCopy[numpy.where(labels == j)] = val
#
# vigra.imshow(labelCopy)
# vigra.show()
