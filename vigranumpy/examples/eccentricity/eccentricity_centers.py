import vigra
import numpy

f = "figure_1.png"

# Shuffle the colors of a labeled image and return the result.
def shuffleColors(labels):
    from random import shuffle
    maxLabel = labels.max()
    mp = numpy.arange(0, maxLabel+1)
    keys = range(maxLabel)
    vals = range(maxLabel)
    shuffle(vals)
    mp[keys] = vals
    return mp[labels]

def drawX(img, pos, color):
    img[pos] = color
    if pos[0] > 0 and pos[1] > 0:
        img[pos[0]-1, pos[1]-1] = color
    if pos[0] > 0 and pos[1] < img.shape[1]-1:
        img[pos[0]-1, pos[1]+1] = color
    if pos[0] < img.shape[0]-1 and pos[1] > 0:
        img[pos[0]+1, pos[1]-1] = color
    if pos[0] < img.shape[0]-1 and pos[1] < img.shape[1]-1:
        img[pos[0]+1, pos[1]+1] = color

# Find the eccentricity centers of each connected component.
img = numpy.squeeze(vigra.readImage(f))
lbl = numpy.array(vigra.analysis.labelImage(img))
centers = vigra.graphs.eccentricityCenters(lbl)

# Draw the centers into the image and show it.
centerImg = shuffleColors(lbl)
maxLabel = centers.shape[0]
for i in range(centers.shape[0]):
    coords = tuple(centers[i, :])
    if centerImg[coords] < maxLabel/2:
        drawX(centerImg, coords, maxLabel)
        # centerImg[coords] = maxLabel
    else:
        drawX(centerImg, coords, 0)
        # centerImg[coords] = 0
vigra.imshow(numpy.swapaxes(centerImg, 1, 0))
vigra.show()
