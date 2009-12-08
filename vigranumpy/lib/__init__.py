from vigranumpycmodule import *
from arraytypes import *

try:
    from vigranumpyfourier import *
except:
    print "vigranumpyfourier not found"
    

def imshow(image):
    import matplotlib
    
    if image.ndim == 3:
        if image.shape[2] != 3:
            raise RuntimeError("showImage(): Multi channel image must have 3 channels.")
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
        raise RuntimeError("showImage(): ndim must be 2 or 3.")

