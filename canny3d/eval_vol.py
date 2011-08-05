import numpy as npy, canny, pdb, math, scipy
from numpy import asfortranarray as fort
from numpy import *
from canny import *
try:
    from scipy import spatial
except:
    print "WARNING: cannot import scipy.spatial"

def filterPoints(arr):
    o = (SHAPE - 1)/2.
    d = sqrt(((arr - o)**2).sum(1))
    return array(arr)[d<5]

def expectedErrorStd(SNR,sigma):
    return 1./SNR * \
           sqrt(3)/(4.*pi**(1/4.)*sqrt(sigma))


def cannySurfelsWithNoises(vol,thresh,sigma,n_scales,optimizer,
                           error=None,shift=0.0):
    if type(vol) is not npy.ndarray:
        vol = sample(vol)
    vol = fort(vol)
    sts = []
    for n_scale in n_scales:
        ##n_scale = npy.exp(-n_scale/20.)*sigma
        if n_scale > 0.0:
            n_scale = exp(-n_scale/20.)*vol.mean()
            print "noise: %f" % n_scale
            noi = npy.random.normal(scale=n_scale,size=vol.shape)
            vol = fort(vol+noi)

        if type(optimizer) is not str:
            optimizer = optimizer.__name__
            
        ret = canny(vol,sigma,thresh,optimizer,
                    True,False,NUMBER_OFCGITS,
                    DIST if FILTER else SHAPE[0]+1,
                    shift)
            
        if error is None:
            sts.append(ret)
        else:
            sts.append(error(ret))

        if DEBUG:
            pdb.vol = vol
            pdb.ret = ret
            pdb.set_trace()

    return sts

