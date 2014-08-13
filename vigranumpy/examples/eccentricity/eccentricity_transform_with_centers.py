import vigra
import numpy

f = "figure_1.png"

# Draw 'X' in image img at position pos with color color.
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

# Find eccentricity and centers of each connected component.
img = numpy.squeeze(vigra.readImage(f))
lbl = numpy.array(vigra.analysis.labelImage(img))
ecc, centers = vigra.graphs.eccentricityTransformWithCenters(lbl)

# Draw the centers into the image and show it.
maxLabel = centers.shape[0]
for i in range(centers.shape[0]):
    coords = tuple(centers[i, :])
    drawX(ecc, coords, 150)
vigra.imshow(numpy.swapaxes(ecc, 1, 0))
vigra.show()
