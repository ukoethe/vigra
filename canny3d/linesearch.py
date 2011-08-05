##try:
##    from minfx.line_search.more_thuente import more_thuente
##    ##from minfx.line_search.nocedal_wright_wolfe import nocedal_wright_wolfe
##except:
##    print "WARNING: import 'minfx.line_search.more_thuente' failed"

from nocedal_wright_wolfe import nocedal_wright_wolfe
from more_thuente import more_thuente
    
from pygsl import minimize
import numpy as npy, pdb, sys
from hourglasscmodule import brentLineSearch as cbrentLineSearchImpl
from hourglasscmodule import moreLineSearch as cmoreLineSearchImpl

this_module = sys.modules[__name__]

def func1(pos, u, vol):#
    def func(alpha,param=None):
        p = pos + alpha*u
        return -vol[p]
    return func

def dfunc1(pos, u, vol):
    def func(alpha):
        p = pos + alpha*u
        d = npy.array(-vol.d(p))
        return npy.dot(d,u)
    return func

def brentLineSearch(p, vol,itmax=10):
    p = npy.array(p)
    u = vol.d(p)
    u/= u.norm()
    u = npy.array(u)
    f = minimize.gsl_function(func1(p,u,vol),None)
    minimizer = minimize.brent(f)
    minimizer.set(0,-10,10)
    for i in range(4):
        status = minimizer.iterate()
        a      = minimizer.x_lower()
        b      = minimizer.x_upper()
        m      = minimizer.minimum()
        status = minimize.test_interval(a, b, 0.001, 0)
        if status == 0:
            break
    p = p + u*m
    return p, None

def cbrentLineSearch(p, vol,itmax=10):
    p = npy.array(p)
    u = vol.d(p)
    u/= u.norm()
    a = cbrentLineSearchImpl(p,u,vol)
    p = p + u*a
    return p, None


def moreLineSearch(x, vol, mu=0.001, eta=0.1,a0=1e-6):
    u = vol.d(x)
    u = npy.array(u / u.norm())
    f = lambda v: -vol[v]
    g = lambda v: -vol.d(v)
    a = more_thuente(f,g,x,f(x),g(x),u,
                     a_init=a0, mu=mu, eta=eta)
    return x + a*u, None

def cmoreLineSearch(p, vol,itmax=10):
    p = npy.array(p)
    u = vol.d(p)
    u/= u.norm()
    a = cmoreLineSearchImpl(p,u,vol)
    p = p + u*a
    return p, None

def nocedalWolfeLineSearch(x, vol, c1=0.001, c2=0.1, a0=1e-10):
    u = vol.d(x)
    u = npy.array(u / u.norm())
    f = lambda v: -vol[v]
    g = lambda v: -vol.d(v)
    
    a = nocedal_wright_wolfe(f,g,x,f(x),g(x),u,
                             max_a=10,
                             a_init=a0)##,mu=c1,eta=c2,print_flag=0)

    if vol[x] > vol[x + a*u]: print 'false'
    return x + a*u, None

def newtonLineSearchStep(p,u,vol):
    H = npy.diag([vol.dxx(p),vol.dyy(p),vol.dzz(p)])
    H[0,1]=H[1,0]=vol.dxy(p)
    H[0,2]=H[2,0]=vol.dxz(p)
    H[1,2]=H[2,1]=vol.dyz(p)
    ##print H
    du = npy.dot(vol.d(p),u)
    a = du/npy.dot(u,npy.dot(H,u))
    return tuple(p - a*u), du

def newtonLineSearch(p,vol,tol=1e-8,maxit=20):
    p = tuple(p)
    u = vol.d(p)
    u = u / u.norm()
    tra = [p]
    for i in range(maxit):
        p,du = newtonLineSearchStep(p,u,vol)
        tra.append(p)
        if du < tol:
            break
    return p, tra

def lineSearch(points,vol,method="newton"):
    method = getattr(this_module, method + 'LineSearch')
    ret = []

    for p in points:
        ##try:
        rp,tra = method(p,vol)
        ##assert vol[rp] >= vol[p]
        ret.append(rp)
        ##except: pass
    
    return ret
