import numpy as npy, canny, pdb, math, scipy
from numpy import asfortranarray as fort
from numpy import *
from canny import *
try:
    from scipy import spatial
except:
    print "WARNING: cannot import scipy.spatial"

### std experiments parameters
USE_ANA_PLANE = True

FSIGMAS     = [0.7,1.0,2.0,3.0]

STD_SIGMA_P = 0.9
STD_SIGMA_F = 2.0

STD_SIGMA   = 2*npy.sqrt(2.)
STD_Q = 1
STD_P = 0

SHAPE     = npy.array((30,30,30))
#SHAPE     = npy.array((50,50,50))
#SHAPE     = npy.array((100,100,100))
O = (SHAPE-1)/2.
TRANS = 0.1
TRESH = 0.1  
ONLY_NEAREST=False
NUMBER_OFCGITS=100
N_SPIRAL_POINTS=50
N_NOISE_SCALES=10
DEBUG = False
FILTER=True
DIST=5
##def SNR()

def filterPoints(arr):
    o = (SHAPE - 1)/2.
    d = sqrt(((arr - o)**2).sum(1))
    return array(arr)[d<5]

def expectedErrorStd(N,S,p,s):
    assert STD_P == 0
    return N/S*sqrt(3.)/(4.*pi**(1/4.)*sqrt(s))
##    
##    return N/S*(sqrt(3)*(p**2+s**2)**(3./2))/(4.*pi**(1/4.)*s**(7./2))


def cannySurfelsWithNoises(vol,thresh,sigma,n_scale,optimizer,error=None,shift=0.0):
    if type(vol) is not npy.ndarray:
        vol = sample(vol)
    vol = fort(vol)
    sts = []

    if n_scale > 0.0:
        noi = npy.random.normal(scale=n_scale,size=vol.shape)
        vol = fort(vol+noi)

    if type(optimizer) is not str:
        optimizer = optimizer.__name__
    print shift
    ret = canny(vol,sigma,thresh,optimizer,
                True,False,NUMBER_OFCGITS,
                DIST if FILTER else SHAPE[0]+1,
                shift)
            
    if DEBUG:
        pdb.vol = vol
        pdb.ret = ret
        pdb.set_trace()

    return ret

def estimateErrorFromPlane(angles,trans,thresh,s_p,s_f,
                           n_scales,optimizer,error=True):
    shape = SHAPE

    if USE_ANA_PLANE:
        ##print angles,trans,s_p,STD_Q,STD_P,shape
        gpv  = GaussianPlane(angles,trans,s_p,STD_Q,STD_P,shape.astype('int32'))
        if DEBUG:
            pdb.gpv = gpv
    else:
        raise NotImplementedError

    def errorf(arr):
        arr = array(arr)
        if ONLY_NEAREST:
            dis = sqrt(((arr-O)**2).sum(1))
            ret = array([arr[dis.argmin()]])
        else:
            ret = arr
        ret = npy.array(planeDistance(gpv,ret.T))
        ##pdb.set_trace()
        print "std: %f"%ret.std()
        print "max: %f"%ret.max()
        print "min: %f"%ret.min()
        ## ##if len(ret) > 1 else ret[0]
        if ONLY_NEAREST:
            return ret[0]
        else:
            return (ret**2).mean(), (ret).std()
    
    ret  = cannySurfelsWithNoises(gpv,thresh,s_f,n_scales,optimizer,
                                  errorf if error else None,
                                  trans)
    return errorf(ret) if error else ret


def exp1(optimizer):
    
    s_p = STD_SIGMA_P
    s_f = STD_SIGMA_F
    thresh = 0.01
    trans = TRANS

    rads1 = npy.linspace(0,npy.pi,100)
    rads2 = npy.linspace(0,npy.pi,100)
    
    ret = []
    for angles in zip([0]*len(rads1),rads1,rads2):
        ret.append(float(estimateErrorFromPlane(angles,trans,thresh,s_p,s_f,
                                                [0.0],optimizer,error=True)))

    return rads1,array(ret)

def exp2(optimizer):
    
    s_p = STD_SIGMA_P
    s_f = STD_SIGMA_F
    thresh = 0.001
    trans = TRANS

    rads1 = npy.linspace(0,npy.pi,10)
    rads2 = npy.linspace(0,npy.pi,10)
    
    ret = []
    for angles in zip([0]*len(rads1),rads1,rads2):
        ret.append(estimateErrorFromPlane(angles,trans,thresh,s_p,s_f,
                                          [20]*20,optimizer,error=True))

    return rads1,array(ret)

def exp3(trans,optimizer,SNR=0.0):
    s_p = STD_SIGMA_P
    s_f = STD_SIGMA_F
    thresh = TRESH

    ts,dirs = spiral(N_SPIRAL_POINTS)
    ret = []
    for d in dirs:
        est = estimateErrorFromPlane(d,trans,thresh,s_p,s_f,SNR,\
                                     optimizer,error=True)
        ret.append(est)
        
    return ts,ret

def testAngles():
    ts,ps = spiral(100)
    ret = []
    for n in ps:
        n = array(n)
        n = n/sqrt((n**2).sum())
        gpv = GaussianPlane(n,0.0,STD_SIGMA_P,\
                            STD_Q,STD_P,SHAPE.astype('int32'))
        vol = sample(gpv)
        svv = SplineVolumeView5(vol)
        grd = gaussianGradient(vol,1)
        d = array(svv.d(O))
        d = d / sqrt((d**2).sum())
        ret.append(dot(d,n))
    return ret

##def spiral(n,k=None):
####    if k is None:
####        k = 0.5*(3.-sqrt(5))*n
####    ts = linspace(-1,1,n)
####    ts = arange(1-(2.0/n)/2.,0,-2./n)
####    return ts,[(cos(pi*t*k)*sqrt(1.-t**2),\
####                sin(pi*t*k)*sqrt(1.-t**2),\
####                t)  for t in ts]


##    dLong = pi*(3.-sqrt(5.))
##    Long = 0
##    dz = 2.0/(n*2)
##    z = 1 - dz/2.
##    pt = [None]*n
##    ts = [None]*n
##    kr = range(0,n)
##    kr.reverse()
##    for k in kr:
##        r = sqrt(1-z*z)
##        pt[k]=(cos(Long)*r, sin(Long)*r, z)
##        ts[k]=z
##        z = z - dz
##        Long = Long + dLong

##    return array(ts),array(pt)

def spiral(N):
    s = 3.6/sqrt(N*2)
    Long = 0
    dz = 2.0/(N*2)
    z = 1 - dz/2
    pt = []
    ts = []
    for k in range(0,N):
        r = sqrt(1-z*z)
        pt.append((cos(Long)*r, sin(Long)*r, z))
        ts.append(z)
        z = z - dz
        Long = Long + s/r

    ts.reverse()
    pt.reverse()
    return array(ts),array(pt)

def plotExp3All(trans=0.0,SNR=0):
    import pgfplots

    end = str(trans).replace('.','_') + "_" + str(SNR)
    
    ts,vals = exp3(trans,"parabolafits",SNR)
    vals = array(vals)
    pgfplots.plot([ts,vals.T[0],vals.T[1]],"latex/plotPlane_para_%s"%end)

    ts,vals = exp3(trans,"brend",SNR)
    vals = array(vals)
    pgfplots.plot([ts,vals.T[0],vals.T[1]],"latex/plotPlane_brend_%s"%end)

    ts,vals = exp3(trans,"mcgBrend",SNR)
    vals = array(vals)
    pgfplots.plot([ts,vals.T[0],vals.T[1]],"latex/plotPlane_cg_%s"%end)

def plotExp4All(SNR=0):
    import pgfplots

    end = str(SNR)
    trans = [0.1,0.25,0.5,0.75]

    def func(opti):
        ret = []
        for t in trans:
            ts,vals = exp3(t,opti,SNR)
            ret.append(array(vals).T[0])
        return ts,array(ret)

    ts,ret = func("parabolafits")
    pgfplots.plot([ts,ret.mean(0),ret.std(0)],"latex/plotPlane_para_%s"%end)

    ts,ret = func("brend")
    pgfplots.plot([ts,ret.mean(0),ret.std(0)],"latex/plotPlane_brend_%s"%end)

    ts,ret = func("mcgBrend")
    pgfplots.plot([ts,ret.mean(0),ret.std(0)],"latex/plotPlane_cg_%s"%end)


def plotExp5All(SNR=0):
    import pgfplots

    trans = [0.1,0.25,0.5,0.75]

    ret = []
    def func(opti):
        for t in trans:
            end = str(t).replace('.','_') + "_" + str(SNR)
            ts,vals = exp3(t,opti,SNR)
            vals = array(vals)
            ret.append(vals)
            pgfplots.plot([ts,vals.T[0],vals.T[1]],"latex/plotPlane_%s_%s"%(opti,end))

    func("parabolafits")
    func("brend")
    func("mcgBrend")

    return ret

def plotExp6All():
    import pgfplots

    trans = [0.1,0.25,0.5,0.75]

    t = trans[0]

    def func(opti):
        ret = []
        for SNR in linspace(0.01,0.08,10):
            end = str(t).replace('.','_') + "_" + str(SNR)
            ts,vals = exp3(t,opti,SNR)
            vals = array(vals)
            ret.append((expectedErrorStd(SNR,STD_Q,STD_SIGMA_P,\
                                         STD_SIGMA_F),\
                        vals.T[1].mean(),vals.T[1].std()))
        return array(ret)


    ret = func("parabolafits")
    pgfplots.plot(ret.T,"latex/exp_std_parabolafits")
    
    ret = func("brend")
    pgfplots.plot(ret.T,"latex/exp_std_brend")
    
    ret = func("mcgBrend")
    pgfplots.plot(ret.T,"latex/exp_std_mcgBrend")

    return ret
