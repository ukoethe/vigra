import vigra
import numpy

f = "figure_1.png"

img = numpy.squeeze(vigra.readImage(f))
lbl = numpy.array(vigra.analysis.labelImage(img))
ecc = vigra.graphs.eccentricityTransform(lbl)
vigra.imshow(numpy.swapaxes(ecc, 1, 0))
vigra.show()
