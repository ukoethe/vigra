import vigra
from vigra import graphs

filepath = '12003.jpg' 
img = vigra.impex.readImage(filepath).astype('float32')
imgM = img.copy()

sigmaS = 8.0
sigmaB = 1.0


for c in range(3):
    imgC = img[:, :, c].squeeze()
    imgM[:,:,c] = vigra.histogram.gaussianRankOrder(imgC,sigmas=(sigmaS, sigmaS, sigmaB), ranks=(0.5,), bins=255).squeeze()
vigra.imshow(vigra.taggedView(imgM,'xyc'))
vigra.show()
