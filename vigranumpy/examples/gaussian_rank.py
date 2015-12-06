import vigra
from vigra import graphs

filepath = '12003.jpg' 
img = vigra.impex.readImage(filepath).astype('float32')
imgM = img.copy()

sigmaS = 1.0
sigmaB = 5.0

with vigra.Timer("compute rank"):
    for c in range(3):
        print "channel",c
        imgC = img[:, :, c].squeeze()
        imgM[:,:,c] = vigra.histogram.gaussianRankOrder(imgC,sigmas=(sigmaS, sigmaS, sigmaB), ranks=(0.5,), bins=100).squeeze()
vigra.imshow(vigra.taggedView(imgM,'xyc'))
vigra.show()
