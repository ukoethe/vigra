import vigra
import numpy
import matplotlib.pyplot as plt

f = "islands.png"
cmap = plt.get_cmap("afmhot")

# Compute labeled iamge.
img = numpy.squeeze(vigra.readImage(f))
lbl = numpy.array(vigra.analysis.labelImage(img))

# Compute the eccentricity transform.
ecc = vigra.filters.eccentricityTransform(lbl)
plt.imshow(numpy.swapaxes(ecc, 1, 0), cmap=cmap)
plt.show()

# Compute the eccentricity centers and draw them into the image.
centers = vigra.filters.eccentricityCenters(lbl)
m = ecc.max()
for c in centers[1:]:
    ecc[c] = m
plt.imshow(numpy.swapaxes(ecc, 1, 0), cmap=cmap)
plt.show()

# # Compute the transformation and the centers in one step:
# ecc, centers = vigra.filters.eccentricityTransformWithCenters(lbl)
