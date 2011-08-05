import numpy as npy, canny, pdb, math, scipy
from numpy import asfortranarray as fort
from numpy import array
from canny import *
try:
    from scipy import spatial
except:
    print "WARNING: cannot import scipy.spatial"

### std experiments parameters
USE_ANA_SPHERE = True

FSIGMAS     = [0.85,1.0,2.0,3.0]
RADI        = [2.,4.,10.,20.]
STD_SIGMA_F = 2.
STD_SIGMA_P = 0.9
STD_SIGMA   = 2*npy.sqrt(2.)
STD_Q = 1
STD_P = 0

##SHAPE     = npy.array((75,75,75))
SHAPE     = npy.array((50,50,50))
O = (SHAPE-1)/2.
TRANS = 0.0
THRESH = 0.02
ONLY_NEAREST=False
NUMBER_OFCGITS=10
N_NOISE_SCALES=2
DEBUG = False
FILTER=True
DIST = 20000
PRINT=False
THIN=False
##def SNR()

def getTresh(s):
    d = {0.85:0.2,1.0:0.21,2.0:0.01,3.0:0.01}
    return d[s]

### helpers

def predictedDisloaction(R,s):
##    func = lambda r: -(r*R + s**2 + math.exp((2*r*R)/s**2)*(r*R - s**2))/ \
##           (math.exp((r + R)**2/(2*s**2))*math.sqrt(2*math.pi)*r**2*s)

    def func(r):
        t1 = (R*(npy.exp(-(r + R)**2/(2*s**2))))/(npy.sqrt(2*npy.pi)*r*s)
        t2 = (s*(npy.exp(-(r + R)**2/(2*s**2))))/(npy.sqrt(2*npy.pi)*r**2)
        t3 = (R*(npy.exp(-(r - R)**2/(2*s**2))))/(npy.sqrt(2*npy.pi)*r*s)
        t4 = -((s*(npy.exp(-(r - R)**2/(2*s**2)))))/(npy.sqrt(2*npy.pi)*r**2)
        ##print R,s,r,t1,t2,t3,t4,t1+t2+t3+t4
        return -(t1+t2+t3+t4)
    
    return npy.asscalar(scipy.optimize.fmin(func,R,disp=False))

    ##return (1.064*npy.sqrt(0.108**2 + (s/R - 0.444)**2) - 0.4835 )*s + R

def norm(a,b=None):
    if b is None:
        ret = (npy.array(a)**2).sum(1)
    else:
        ret = a**2 + b**2
    return npy.sqrt(ret)

### experiments
def cannySurfelsWithNoises(vol,thresh,sigma,n_scales,optimizer,error=None,shift=0):
    if type(vol) is not npy.ndarray:
        vol = sample(vol)
    vol = fort(vol)
    print "mean: %s" % vol.mean()
    sts = []
    for n_scale in n_scales:
        ##n_scale = npy.exp(-n_scale/20.)*sigma
        if n_scale > 0.0:
            if PRINT:
                print "noise: %f" % n_scale
            noi = npy.random.normal(scale=n_scale,size=vol.shape)
            vol = fort(vol+noi)

        if type(optimizer) is not str:
            optimizer = optimizer.__name__

        ret = canny(vol,sigma,thresh,optimizer,
                    THIN,PRINT,NUMBER_OFCGITS,
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


def estimateRadiusFromShiftedSphere(r,thresh,s_p,s_f,n_scales,optimizer,shift=(0,0,0), error=True):
    shape = SHAPE
    o = (shape-1.)/2. + shift
    errorf = lambda arr:\
             npy.sqrt(((arr - o)**2).sum(1))

    if USE_ANA_SPHERE:
        gpv  = GaussianSphere(r,s_p,STD_Q,STD_P,shape.astype('int32'),o)
    else:
        gpv  = STD_P + (-STD_P + STD_Q)*createSphere(r,s_p,shape,o)

    assert shift[0]==shift[1]==shift[2]
    sts  = cannySurfelsWithNoises(gpv,thresh,s_f,n_scales,optimizer,
                                  errorf if error else None,shift[0])
    return npy.array(sts)

def estimateRadiusFromShiftedSpheres(r,thresh,s_p,s_f,n_scales,optimizer,n_shifts=0):
    shift = npy.linspace(0.0,npy.sqrt(3),n_shifts)
    shift = zip(shift,shift,shift)
    ret = []
    for sh in list(shift)+[(0,0,0)]:
        sts = estimateRadiusFromShiftedSphere(r,thresh,s_p,s_f,n_scales,optimizer,sh)
        sts = npy.array(sts)
        ret.append(sts.mean())
    return npy.array(ret)

def exp2(optimizer,SNR=0):
    """ first example on sphere radius estimation
    with various radi and sigmas"""
    
    s_p = STD_SIGMA_P
    rs  = RADI
    sfs = FSIGMAS
    
    ret = []
    pre = []
    for R in RADI:
        for s_f in FSIGMAS:
            s = norm(s_f,s_p)
            r = predictedDisloaction(R,s)
            print s,R,r
            radi = estimateRadiusFromShiftedSpheres(R,getTresh(s_f),s_p,s_f,[SNR],optimizer)
            ret.append(((radi.mean()-R)/s,radi.std()))
            pre.append((r-R)/s)

            print radi.mean()

##    sfs = npy.array(zip(sfs,len(sfs)*[s_p]))
##    sfs = npy.sqrt((sfs**2).sum(1))
##    rs  = 1./npy.array(rs)
##    params = (rs[:,npy.newaxis]*sfs).ravel()
    pre = array(pre)
    ret = array(ret)
    return pre,ret
            

def testSNRSphereRadius(thresh,s_p,s_f,n_scales,optimizer):
    rs = npy.linspace(2,10,20)
    ret = []
    for r in rs:
        R = realRadius(r,math.sqrt(s_p**2+s_f**2))
        ###print r,thresh,sigma,n_scales,optimizer.__name__,sh
        sts = testSNRSphere(r,thresh,s_p,s_f,n_scales,optimizer)
        sts = npy.array(sts)
        ret.append((sts.mean()-R,sts.std()))
    return npy.array(ret)
    
def meanDistNearestNeigh(points,r=1e-06):
    tree = spatial.kdtree.KDTree(points)
    ret  = tree.query_ball_point(points,r)
    c = npy.array([len(l) for l in ret])
    return c.mean()

def plotExp2(optimizer):
    import pgfplots
    
    the,est = exp2(optimizer)
    est = array(est)
    pgfplots.plot([the,est.T[0],est.T[1]],"latex/plotSphere_%s_dis"%optimizer)

def plotPredislo():
    import pgfplots
    
    vals = []
    for s in npy.linspace(0.001,1):
        vals.append((s,predictedDisloaction(1,s)-1/s))
    vals = npy.array(vals)
    pgfplots.plot(vals.T,"latex/predis_sphere")
             

def plotExp4(optimizer,SNR):
    import pgfplots
    
    the,est = exp2(optimizer,SNR)
    est = array(est)
    pgfplots.plot([the,est.T[0],est.T[1]],"latex/plotSphere_%s_dis_%s"%(optimizer,str(SNR)))

    return the,est
