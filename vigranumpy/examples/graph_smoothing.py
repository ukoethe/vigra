import vigra
from vigra import graphs

# parameter:
filepath = '12003.jpg'   # input image path
sigmaGradMag = 2.0       # sigma Gaussian gradient
superpixelDiameter = 10  # super-pixel size
slicWeight = 10.0        # SLIC color - spatial weight
gamma = 0.15             # exp(-gamma * edgeIndicator)
edgeThreshold = 2.5      # values higher are considered as edges
scale = 1.0              # how much smoothing
iterations = 10          # how man smoothing iterations

# load image and convert to LAB
img = vigra.impex.readImage(filepath)

res = vigra.filters.nonlinearDiffusion(img, 1.9, 20.0)

vigra.imshow(res)
vigra.show()