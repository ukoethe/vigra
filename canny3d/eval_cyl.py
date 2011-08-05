from eval_vol import *
### std experiments parameters
USE_ANA = True

FSIGMAS     = [4.0,3.5,3.0,2.5,2.0,1.0,0.7,0.5,0.2]
RADI        = [10.]
STD_SIGMA_P = 0.9
STD_SIGMA_F = 1.
STD_SIGMA   = 2*npy.sqrt(2.)
STD_Q = 1
STD_P = 0

SHAPE     = npy.array((50,50,50),dtype='int32')
#SHAPE     = npy.array((50,50,50))
#SHAPE     = npy.array((100,100,100))
O = (SHAPE-1)/2.
TRANS = 0.1
TRESH = 0.01
ONLY_NEAREST=False
NUMBER_OFCGITS=10
N_SPIRAL_POINTS=50
N_NOISE_SCALES=10
DEBUG = False
FILTER=True
PRINT=True
DIST=15

##def SNR()

def init():
    import eval_vol
    reload(eval_vol)
    for name,val in globals().iteritems():
        if name.isupper():
            setattr(eval_vol,name,val)

init()

def predictedDisloaction(R,s):
    assert s/R < 0.5
    return 0.52*sqrt(0.12**2+(s/R-0.476)**2)-0.255

def norm(a,b=None):
    if b is None:
        ret = (npy.array(a)**2).sum(1)
    else:
        ret = a**2 + b**2
    return npy.sqrt(ret)

def estimateRadii(r,angles,shift,thresh,s_p,s_f,
                  n_scales,optimizer,error=True):
    shape = SHAPE
    o = (shape-1.)/2. + shift

    if USE_ANA:
        ##print angles,trans,s_p,STD_Q,STD_P,shape
        gpv  = GaussianCylinder(r,s_p,\
                                STD_Q,STD_P,shape.astype('int32'),
                                angles[0],angles[1],shift)
        if DEBUG:
            pdb.gpv = gpv
    else:
        gpv  = STD_P + (-STD_P + STD_Q)*createCylinder(r,s_p,shape,o,
                                                       angles[0],angles[1])

    def errorf(arr):
        arr = array(arr)
        if USE_ANA:
            dis = cylinderDistance(gpv,arr.T)
        else:
            raise NotImplementedError, "impl shift"
            a,b=angles
            dis = []
            for p in arr:
                p = p-o;
                p[0] = p[0]*cos(b)+p[2]*sin(b);
                p[1] = p[1]*cos(a)-p[2]*cos(b)*sin(a)+p[0]*sin(b)*sin(a);
                p[2] = p[2]*cos(b)*cos(a)-p[0]*cos(a)*sin(b)+p[1]*sin(a);
                dis.append(sqrt(p[0]*p[0]+p[1]*p[1]))

        if DEBUG:
            pdb.set_trace()
        if ONLY_NEAREST:
            return array([dis[dis.argmin()]])
        else:
            return dis 
    
    sts  = cannySurfelsWithNoises(gpv,thresh,s_f,n_scales,optimizer,
                                  errorf if error else None,shift)
    return npy.array(sts)

def spiral(N):
    s = 3.6/sqrt(N*2)
    Long = 0
    dz = 2.0/(N*2)
    z = 1 - dz/2
    pt = []
    ts = []
    for k in range(0,N):
        r = sqrt(1-z*z)
        pt.append((Long, arcsin(z)))
        ts.append(z)
        z = z - dz
        Long = Long + s/r

    ts.reverse()
    pt.reverse()
    return array(ts),array(pt)

def testSpiral(N):

    import view

    ts,ps = spiral(N)
    x,y,z=(0,0,1)
    ret = []
    for (a,b) in ps:
        ret.append(\
            (x*cos(b)+z*sin(b),
             y*cos(a)-z*cos(b)*sin(a)+x*sin(b)*sin(a),
             z*cos(b)*cos(a)-x*cos(a)*sin(b)+y*sin(a)))

    view.clf()
    view.points(ret)

def dislocation(thresh,s_p,s_f,n_scales,optimizer):
    rs = npy.linspace(2,10,20)
    ret = []
    for r in rs:
        R = realRadius(r,math.sqrt(s_p**2+s_f**2))
        ###print r,thresh,sigma,n_scales,optimizer.__name__,sh
        sts = estimateRadii(r,(0,0),0,thresh,s_p,s_f,[0.0],optimizer)
        sts = npy.array(sts)
        ret.append((sts.mean()-R,sts.std()))
    return npy.array(ret)

def exp2(optimizer):
    """ first example on cylinder radius estimation
    with various radi and sigmas"""
    
    s_p = STD_SIGMA_P
    rs  = RADI
    sfs = FSIGMAS
    thresh = TRESH
    
    est = []
    the = []
    for R in RADI:
        for s_f in FSIGMAS:
            s = norm(s_f,s_p)
            dis = predictedDisloaction(R,s)
            radi = estimateRadii(R,(0,0),0,thresh,s_p,s_f,[0.0],optimizer)
            est.append(((radi.mean()-R)/s,radi.std()))
            the.append(dis)

##    sfs = npy.array(zip(sfs,len(sfs)*[s_p]))
##    sfs = npy.sqrt((sfs**2).sum(1))
##    rs  = 1./npy.array(rs)
##    params = (rs[:,npy.newaxis]*sfs).ravel()
    the = array(the)
    est = array(est)
    return the,est

def exp3(R,trans,optimizer,SNR=0.0):
    s_p = STD_SIGMA_P
    s_f = STD_SIGMA_F
    thresh = TRESH

    ts,angles = spiral(N_SPIRAL_POINTS)
    ret = []
    for (a,b) in angles:
        est = estimateRadii(R,(a,b),trans,thresh,s_p,s_f,[0.0],optimizer)
        ret.append(est.mean())
        
    return ts,ret

def exp5(R,trans,optimizer,SNR=0.0):
    s_p = STD_SIGMA_P
    s_f = STD_SIGMA_F
    thresh = TRESH

    ts,angles = spiral(N_SPIRAL_POINTS)
    ret = []
    for (a,b) in angles:
        est = estimateRadii(R,(a,b),trans,thresh,s_p,s_f,[0.0],
                            optimizer,error=False)
        ret.append(est)
        
    return ts,ret


def plotExp2(optimizer="brend"):
    import pgfplots
    
    the,est = exp2(optimizer)
    pgfplots.plot([the,est.T[0],est.T[1]],"latex/plotCyl_%s_dis"%optimizer)

    return the,est

def plotExp3(optimizer="brend",SNR=0):
    import pgfplots
    R = 5
    trans=TRANS
    end = str(trans).replace('.','_') + "_" + str(SNR)
    
    s = norm(STD_SIGMA_F,STD_SIGMA_F)
    ts,est = exp3(R,trans,optimizer,SNR=SNR)
    ts = list(ts)
    ts.reverse()
    est = array(est)-R-predictedDisloaction(R,s)
    pgfplots.plot(array([ts,est]),"latex/plotCyl_%s_%s"%(optimizer,end))

    return ts,est

def plotExp3All(R=5,SNR=0):
    import pgfplots

    trans = [0.0]

    s = norm(STD_SIGMA_F,STD_SIGMA_F)

    ret = []
    def func(optimizer):
        for t in trans:
            end = str(t).replace('.','_') + "_" + str(SNR)
            ts,est = exp3(R,t,optimizer,SNR=SNR)
            ts = list(ts)
            ts.reverse()
            est = (array(est)-R)+predictedDisloaction(R,s)
            ret.append(est)
            pgfplots.plot(array([ts,est]),"latex/plotCyl_%s_%s"%(optimizer,end))

    func("parabolafits")
    func("brend")
    func("mcgBrend")

    return ret
