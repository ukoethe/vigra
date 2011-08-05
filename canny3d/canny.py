import vigra, numpy as npy, pdb, scipy, math

try:
    from scipy.optimize import brent
except:
    print "WARNING: cannot import scipy.optimize.brent"
try:
    from pygsl import _numobj as numx
    from pygsl import multiminimize, errno
    from pygsl import minimize
    from pygsl import multiroots
    from pygsl import multifit
except:
    print "WARNING: cannot import pygsl"

from hourglasscmodule import *

def cgMore1Iter(vol,pos):
    return cgMore(vol,pos,1,False)

def mcgMore1Iter(vol,mgn,pos):
    return mcgMore(vol,mgn,pos,1,False)

def sample(vol,s=None,oversample=False):
    os = npy.array(vol.shape)
    if oversample:
        assert s is not None
    if oversample:
        d = (os-1.)/(npy.array(s)-1.)
    else:
        d = [1,1,1]
    m = npy.asfortranarray(\
            npy.mgrid[0:os[0]:d[0],0:os[1]:d[1],0:os[2]:d[2]].astype('uint32'))
    
    return vol[m]

def indices(vol):
    return npy.asfortranarray(scipy.indices(vol.shape,dtype='uint32'))

def pythinning(surfels,shape):
    shape = [int(i) for i in shape]
    surfels_thin = Surfels()
    thinning(surfels,surfels_thin,shape)
    return npy.array(surfels_thin)

def getSurfels(arr,thresh,mgn=None,grd=None,sigma=2.0,noise=None,thin=True):
    ##arr = sample(vol)
    if noise != None:
        arr += npy.random.normal(scale=noise,size=arr.shape)
    if mgn is None:
        mgn = gradientMagnitude(arr,sigma)
    if grd is None:
        grd = gaussianGradient(arr,sigma)
    grd = npy.asfortranarray(grd.transpose(3,0,1,2))
    mgn = npy.asfortranarray(mgn)
    surfels = Surfels()
    findSurfels(grd,mgn,thresh,surfels)
    if thin:
        return pythinning(surfels,mgn.shape)
    return npy.array(surfels)

def parabolafit(p,vol,d=None):
    p = vigra.Vector3(p[0],p[1],p[2])
    if d is None:
        d = vol.d(p)
    d = d / d.norm()
##    x = npy.r_[-1:1.1:1]
    x = npy.array([-1,0,1])
    y = npy.array([vol[tuple(p-i*d)] for i in x])
    n = len(x)
    X = numx.ones((n,3),)*1.
    X[:,0] = 1.0
    X[:,1] = x
    X[:,2] = x ** 2
    work = multifit.linear_workspace(n,3)
    c, cov, chisq = multifit.linear(X, y, work)
    return p - d*(-c[1]/(2.*c[2]))##FIXME

def parabolafits(vol,points):
    ret = []
    if points.shape[1] != 3:
        points = points.T
    for p in points:
        rp = parabolafit(p,vol)
        ##assert vol[rp] >= vol[p]
        ret.append(rp)
    return ret
