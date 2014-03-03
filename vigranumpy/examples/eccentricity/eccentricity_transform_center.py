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

    ## Distance in the inside between two pixels is 1.0
    distImage = bigImg.copy().astype(numpy.float32)
    distImage[numpy.where(bigImg==0)]=1.0

    ## Set the outside to a very high value
    distImage[numpy.where(bigImg==1)]=10000.0

    ## Image copy to draw the paths in
    imgp = distImage.copy()
    imgp[numpy.where(bigImg==1)] = 100

    ## Get image graph and its path finder
    gridGraph = vigraph.gridGraph(bigImg.shape[0:dim],False)
    graphShape = []
    for s in distImage.shape:
        graphShape.append(s*2-1)
    edgeWeights = vigra.resize(distImage, graphShape, order=0)
    edgeWeights = vigra.graphs.edgeFeaturesFromInterpolatedImageCorrected(gridGraph,edgeWeights)
    pathFinder = vigraph.ShortestPathPathDijkstra(gridGraph)

    ## Find borders in img (TODO: replace with graph functions, maybe regionImageToCrackEdgeImage ( labelImage ) )
    borderImg = numpy.zeros(bigImg.shape)
    if dim == 2:
        for y in range(bigImg.shape[1]-1):
            for x in range(bigImg.shape[0]-1):
                if bigImg[x,y] == 0:
                    if bigImg[x+1,y] == 1 or bigImg[x,y+1] == 1:
                        borderImg[x, y] = 1
                else:
                    if bigImg[x+1,y] == 0:
                        borderImg[x+1, y] = 1
                    if bigImg[x,y+1] == 0:
                        borderImg[x, y+1] = 1
    else:
        for z in range(bigImg.shape[2]-1):
            for y in range(bigImg.shape[1]-1):
                for x in range(bigImg.shape[0]-1):
                    if bigImg[x,y,z] == 0:
                        if bigImg[x+1,y,z] == 1 or bigImg[x,y+1,z] == 1 or bigImg[x,y,z+1] == 1:
                            borderImg[x, y, z] = 1
                    else:
                        if bigImg[x+1,y,z] == 0:
                            borderImg[x+1,y,z] = 1
                        if bigImg[x,y+1,z] == 0:
                            borderImg[x,y+1,z] = 1
                        if bigImg[x,y,z+1] == 0:
                            borderImg[x,y,z+1] = 1

    #vigra.imshow(numpy.swapaxes(borderImg, 1, 0))
    #vigra.show()
    #exit()

    ## End points for paths (all points on the border)
    targets = numpy.where(borderImg==1)
    nTargets = len(targets[0])








    ## Find the diameter (longest path)
    eccLength = numpy.zeros(nTargets)
    #eccTargetIndex = numpy.zeros(nTargets)
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
                #eccTargetIndex[j] = vpIndex
                eccTargetPath[j] = pathFinder.path(pathType='coordinates', target=target)
        vpIndex = numpy.argmax(eccLength)
        vpGraphIndex = []
        for d in range(dim):
            vpGraphIndex.append(int(targets[d][vpIndex]))
        vp = gridGraph.coordinateToNode(vpGraphIndex)
        if visited[vpIndex] or not eccChanged:
            break

    ## Find the midpoint of the diameter
    dMax = eccLength[vpIndex]/2
    path = eccTargetPath[vpIndex]
    p1 = path[0]
    d = 0
    for k in range(1, len(path)):
        p2 = path[k]
        d = d + numpy.linalg.norm(p2-p1)
        if d > dMax:
            break
        p1 = p2

    print eccTargetPath[vpIndex]
    print "------------"
    print p1
    exit()
    ## TODO: Den Mittelpunkt des laengsten Pfades finden.


        #sourceIndex = []
        #for d in range(dim):
        #    sourceIndex.append(int(targets[d][i]))
        #source = gridGraph.coordinateToNode(sourceIndex)
        #pathFinder.run(edgeWeights, source)
        #for k in range(nTargets):











    ## Indices of start points for paths (random)
    nPoints = int(numpy.ceil(percentageOfPaths * nTargets / 100.0))
    numpy.random.seed(42)
    starts = numpy.random.permutation(range(nTargets))[:nPoints]

    ## Compute paths
    maxPaths = []
    maxPathLengths = []
    for i in range(nPoints):
        sourceIndex = []
        for d in range(dim):
            sourceIndex.append(int(targets[d][starts[i]]))
        source = gridGraph.coordinateToNode(sourceIndex)
        pathFinder.run(edgeWeights, source)
        maxPathLength = 0
        for j in range(nTargets):
            targetIndex = []
            for d in range(dim):
                targetIndex.append(int(targets[d][j]))
            target = gridGraph.coordinateToNode(targetIndex)
            path = pathFinder.path(pathType='coordinates', target=target)
            pathLength = pathFinder.distance(target)
            if pathLength > maxPathLength or maxPathLength == 0:
                maxPathLength = pathLength
                maxPath = path
        maxPaths.append(maxPath)
        maxPathLengths.append(maxPathLength)
        imgp[sourceIndex[0], sourceIndex[1]] = 40*maxPathLength

    if showPathImage:
        val = (imgp.max()+imgp.min())/2
        for p in maxPaths:
            imgp[p[:,0], p[:,1]] = val
    if showPathImage:
        plt.figure(distFunc)
        plt.imshow(numpy.swapaxes(imgp, 1, 0), interpolation='none')
    if len(imgSaveName)>1:

        scipy.misc.imsave(imgSaveName, numpy.swapaxes(imgp, 1, 0))

    return maxPathLengths



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
    saveName = "ecc_seg_"+`counter`+".png"
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
