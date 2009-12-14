import vigranumpycmodule
from vigranumpycmodule import *

import arraytypes
from arraytypes import *

try:
    import vigranumpyfourier as fourier
    from vigranumpyfourier import *
except:
    print "cannot import vigra.fourier"
    
try:
    import pyqt
    from pyqt import showImage
except:
    print "cannot import vigra.pyqt"
    

def imshow(image):
    import matplotlib
    
    if image.ndim == 3:
        if image.shape[2] != 3:
            raise RuntimeError("vigra.imshow(): Multi channel image must have 3 channels.")
        if image.dtype != uint8:
            mi, ma = image.min(), image.max()
            if mi >= ma:
                image = image.__class__(image.shape, dtype=uint8)
            else:
                image = (image - mi) / (ma - mi)
        return matplotlib.pyplot.imshow(image.swapaxes(0,1))
    elif image.ndim == 2:
        return matplotlib.pyplot.imshow(image.swapaxes(0,1), cmap=matplotlib.cm.gray, \
                                     norm=matplotlib.cm.colors.Normalize())
    else:
        raise RuntimeError("vigra.imshow(): ndim must be 2 or 3.")

