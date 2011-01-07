/************************************************************************/
/*                                                                      */
/*               Copyright 2010-2011 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */                
/*                                                                      */
/************************************************************************/

// code has been adapted from C.R. Bond's implementation
// see http://www.crbond.com/math.htm

#ifndef VIGRA_BESSEL_HXX
#define VIGRA_BESSEL_HXX

#include "mathutil.hxx"

namespace vigra {

/** \addtogroup MathFunctions 
*/
//@{

namespace detail {

template <class REAL>
int msta1(REAL x, int mp)
{
    double a0,f0,f1,f;
    int i,n0,n1,nn;

    a0 = abs(x);
    n0 = (int)(1.1*a0)+1;
    f0 = 0.5*std::log10(6.28*n0) - n0*std::log10(1.36*a0/n0)-mp;
    n1 = n0+5;
    f1 = 0.5*std::log10(6.28*n1) - n1*std::log10(1.36*a0/n1)-mp;
    for(i=0;i<20;i++) 
    {
        nn = int(n1-(n1-n0)/(1.0-f0/f1));
        f = 0.5*std::log10(6.28*nn) - nn*std::log10(1.36*a0/nn)-mp;
        if(abs(nn-n1) < 1) 
            break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn;
}

template <class REAL>
int msta2(REAL x, int n, int mp)
{
    double a0,ejn,hmp,f0,f1,f,obj;
    int i,n0,n1,nn;

    a0 = abs(x);
    hmp = 0.5*mp;
    ejn = 0.5*std::log10(6.28*n) - n*std::log10(1.36*a0/n);
    if (ejn <= hmp) 
    {
        obj = mp;
        n0 = (int)(1.1*a0);
        if (n0 < 1) 
            n0 = 1;
    }
    else 
    {
        obj = hmp+ejn;
        n0 = n;
    }
    f0 = 0.5*std::log10(6.28*n0) - n0*std::log10(1.36*a0/n0)-obj;
    n1 = n0+5;
    f1 = 0.5*std::log10(6.28*n1) - n1*std::log10(1.36*a0/n1)-obj;
    for (i=0;i<20;i++) 
    {
        nn = int(n1-(n1-n0)/(1.0-f0/f1));
        f = 0.5*std::log10(6.28*nn) - nn*std::log10(1.36*a0/nn)-obj;
        if (abs(nn-n1) < 1) 
            break;
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn+10;
}

template <class REAL>
int bessjyv(double v,REAL x,double &vm,double *jv,double *yv,
            double *djv,double *dyv)
{
    double v0,vl,vg,vv,a,a0,r,x2,bjv0 = 0.0,bjv1 = 0.0,bjvl = 0.0,f = 0.0,f0,f1,f2;
    double r0,r1,ck,cs,cs0,cs1,sk,qx,px,byv0 = 0.0,byv1 = 0.0,rp,xk,rq;
    double b,ec,w0,w1,bju0 = 0.0,bju1 = 0.0,byvk = 0.0;
    double el = 0.5772156649015329; // Euler's gamma
    int j,k,l,m,n,kz;

    x2 = x*x;
    n = (int)v;
    v0 = v-n;
    if ((x < 0.0) || (v < 0.0)) 
        return 1;
    if (x < 1e-15) 
    {
        for (k=0;k<=n;k++) 
        {
            jv[k] = 0.0;
            yv[k] = -1e308;
            djv[k] = 0.0;
            dyv[k] = 1e308;
            if (v0 == 0.0) 
            {
                jv[0] = 1.0;
                djv[1] = 0.5;
            }
            else 
                djv[0] = 1e308;
        }
        vm = v;
        return 0;
    }
    if (x <= 12.0) 
    {
        for (l=0;l<2;l++) 
        {
            vl = v0 + l;
            bjvl = 1.0;
            r = 1.0;
            for (k=1;k<=40;k++) 
            {
                r *= -0.25*x2/(k*(k+vl));
                bjvl += r;
                if (abs(r) < abs(bjvl)*1e-15) 
                    break;
            }
            vg = 1.0 + vl;
            a = pow(0.5*x,vl)/gamma(vg);
            if (l == 0) 
                bjv0 = bjvl*a;
            else 
                bjv1 = bjvl*a;
        }
    }
    else 
    {
        if (x >= 50.0) kz = 8;
        else if (x >= 35.0) kz = 10;
        else kz = 11;
        for (j=0;j<2;j++) 
        {
            vv = 4.0*(j+v0)*(j+v0);
            px = 1.0;
            rp = 1.0;
            for (k=1;k<=kz;k++) 
            {
                rp *= (-0.78125e-2)*(vv-pow(4.0*k-3.0,2.0))*
                    (vv-pow(4.0*k-1.0,2.0))/(k*(2.0*k-1.0)*x2);
                px += rp;
            }
            qx = 1.0;
            rq = 1.0;
            for (k=1;k<=kz;k++) 
            {
                rq *= (-0.78125e-2)*(vv-pow(4.0*k-1.0,2.0))*
                    (vv-pow(4.0*k+1.0,2.0))/(k*(2.0*k+1.0)*x2);
                qx += rq;
            }
            qx *= 0.125*(vv-1.0)/x;
            xk = x-(0.5*(j+v0)+0.25)*M_PI;
            a0 = std::sqrt(M_2_PI/x);
            ck = std::cos(xk);
            sk = std::sin(xk);

            if (j == 0) 
            {
                bjv0 = a0*(px*ck-qx*sk);
                byv0 = a0*(px*sk+qx*ck);
            }
            else if (j == 1) 
            {
                bjv1 = a0*(px*ck-qx*sk);
                byv1 = a0*(px*sk+qx*ck);
            }
        }
    }
    jv[0] = bjv0;
    jv[1] = bjv1;
    djv[0] = v0*jv[0]/x-jv[1];
    djv[1] = -(1.0+v0)*jv[1]/x+jv[0];
    if ((n >= 2) && (n <= (int)(0.9*x))) 
    {
        f0 = bjv0;
        f1 = bjv1;
        for (k=2;k<=n;k++) 
        {
            f = 2.0*(k+v0-1.0)*f1/x-f0;
            jv[k] = f;
            f0 = f1;
            f1 = f;
        }
    }
    else if (n >= 2) 
    {
        m = msta1(x,200);
        if (m < n) n = m;
        else m = msta2(x,n,15);
        f2 = 0.0;
        f1 = 1.0e-100;
        for (k=m;k>=0;k--) 
        {
            f = 2.0*(v0+k+1.0)*f1/x-f2;
            if (k <= n) jv[k] = f;
            f2 = f1;
            f1 = f;
        }
        if (abs(bjv0) > abs(bjv1)) 
            cs = bjv0/f;
        else 
            cs = bjv1/f2;
        for (k=0;k<=n;k++) 
        {
            jv[k] *= cs;
        }
    }
    for (k=2;k<=n;k++) 
    {
        djv[k] = -(k+v0)*jv[k]/x+jv[k-1];
    }
    if (x <= 12.0) 
    {
        if (v0 != 0.0) 
        {
            for (l=0;l<2;l++) 
            {
                vl = v0 +l;
                bjvl = 1.0;
                r = 1.0;
                for (k=1;k<=40;k++) 
                {
                    r *= -0.25*x2/(k*(k-vl));
                    bjvl += r;
                    if (abs(r) < abs(bjvl)*1e-15) 
                        break;
                }
                vg = 1.0-vl;
                b = pow(2.0/x,vl)/gamma(vg);
                if (l == 0) 
                    bju0 = bjvl*b;
                else 
                    bju1 = bjvl*b;
            }
            byv0 = (bjv0*cos_pi(v0)-bju0)/sin_pi(v0);
            byv1 = (bjv1*cos_pi(v0+1.0)-bju1)/sin_pi(v0+1.0);
        }
        else 
        {
            ec = std::log(0.5*x)+el;
            cs0 = 0.0;
            w0 = 0.0;
            r0 = 1.0;
            for (k=1;k<=30;k++) 
            {
                w0 += 1.0/k;
                r0 *= -0.25*x2/(k*k);
                cs0 += r0*w0;
            }
            byv0 = M_2_PI*(ec*bjv0-cs0);
            cs1 = 1.0;
            w1 = 0.0;
            r1 = 1.0;
            for (k=1;k<=30;k++) 
            {
                w1 += 1.0/k;
                r1 *= -0.25*x2/(k*(k+1));
                cs1 += r1*(2.0*w1+1.0/(k+1.0));
            }
            byv1 = M_2_PI*(ec*bjv1-1.0/x-0.25*x*cs1);
        }
    }
    yv[0] = byv0;
    yv[1] = byv1;
    for (k=2;k<=n;k++) 
    {
        byvk = 2.0*(v0+k-1.0)*byv1/x-byv0;
        yv[k] = byvk;
        byv0 = byv1;
        byv1 = byvk;
    }
    dyv[0] = v0*yv[0]/x-yv[1];
    for (k=1;k<=n;k++) 
    {
        dyv[k] = -(k+v0)*yv[k]/x+yv[k-1];
    }
    vm = n + v0;
    return 0;
}

template <class REAL>
int bessikv(double v,REAL x,double &vm,
            double *iv,double *kv, double *ivp,double *kvp)
{
    double x2,v0,piv,vt,a1,v0p,gap = 0.0,r,bi0,ca,sum;
    double f = 0.0,f1,f2,ct,cs,wa,gan,ww,w0,v0n;
    double r1,r2,bk0,bk1,bk2,a2,cb;
    double el = 0.5772156649015329; // Euler's gamma
    double eps = 1e-15;
    int n,k,kz,m;

    if ((v < 0.0) || (x < 0.0)) 
        return 1;
    x2 = x*x;
    n = (int)v;
    v0 = v-n;
    if (n == 0) n = 1;
    if (x == 0.0) 
    {
        for (k=0;k<=n;k++) 
        {
            iv[k] = 0.0;
            kv[k] = -1e308;
            ivp[k] = 0.0;
            kvp[k] = 1e308;
        }
        if (v0 == 0.0) 
        {
            iv[0] = 1.0;
            ivp[1] = 0.5;
        }
        vm = v;
        return 0;
    }
    piv = M_PI*v0;
    vt = 4.0*v0*v0;
    if (v0 == 0.0) 
    {
        a1 = 1.0;
    }
    else 
    {
        v0p = 1.0+v0;
        gap = gamma(v0p);
        a1 = pow(0.5*x,v0)/gap;
    }
    if (x >= 50.0) 
        kz = 8;
    else if (x >= 35.0) 
        kz = 10;
    else 
        kz = 14;
    if (x <= 18.0) 
    {
        bi0 = 1.0;
        r = 1.0;
        for (k=1;k<=30;k++) 
        {
            r *= 0.25*x2/(k*(k+v0));
            bi0 += r;
            if (abs(r/bi0) < eps) 
                break;
        }
        bi0 *= a1;
    }
    else 
    {
        ca = std::exp(x)/std::sqrt(2.0*M_PI*x);
        sum = 1.0;
        r = 1.0;
        for (k=1;k<=kz;k++) 
        {
            r *= -0.125*(vt-pow(2.0*k-1.0,2.0))/(k*x);
            sum += r;
        }
        bi0 = ca*sum;
    }
    m = msta1(x,200);
    if (m < n) 
        n = m;
    else 
        m = msta2(x,n,15);
    f2 = 0.0;
    f1 = 1.0e-100;
    for (k=m;k>=0;k--) 
    {
        f = 2.0*(v0+k+1.0)*f1/x+f2;
        if (k <= n) iv[k] = f;
        f2 = f1;
        f1 = f;
    }
    cs = bi0/f;
    for (k=0;k<=n;k++) 
    {
        iv[k] *= cs;
    }
    ivp[0] = v0*iv[0]/x+iv[1];
    for (k=1;k<=n;k++) 
    {
        ivp[k] = -(k+v0)*iv[k]/x+iv[k-1];
    }
    ww = 0.0;
    if (x <= 9.0) 
    {
        if (v0 == 0.0) 
        {
            ct = -log(0.5*x)-el;
            cs = 0.0;
            w0 = 0.0;
            r = 1.0;
            for (k=1;k<=50;k++) 
            {
                w0 += 1.0/k;
                r *= 0.25*x2/(k*k);
                cs += r*(w0+ct);
                wa = abs(cs);
                if (abs((wa-ww)/wa) < eps) 
                    break;
                ww = wa;
            }
            bk0 = ct+cs;
        }
        else 
        {
            v0n = 1.0-v0;
            gan = gamma(v0n);
            a2 = 1.0/(gan*pow(0.5*x,v0));
            a1 = pow(0.5*x,v0)/gap;
            sum = a2-a1;
            r1 = 1.0;
            r2 = 1.0;
            for (k=1;k<=120;k++) 
            {
                r1 *= 0.25*x2/(k*(k-v0));
                r2 *= 0.25*x2/(k*(k+v0));
                sum += a2*r1-a1*r2;
                wa = abs(sum);
                if (abs((wa-ww)/wa) < eps) 
                    break;
                ww = wa;
            }
            bk0 = M_PI_2*sum/std::sin(piv);
        }
    }
    else 
    {
        cb = std::exp(-x)*std::sqrt(M_PI_2/x);
        sum = 1.0;
        r = 1.0;
        for (k=1;k<=kz;k++) 
        {
            r *= 0.125*(vt-pow(2.0*k-1.0,2.0))/(k*x);
            sum += r;
        }
        bk0 = cb*sum;
    }
    bk1 = (1.0/x-iv[1]*bk0)/iv[0];
    kv[0] = bk0;
    kv[1] = bk1;
    for (k=2;k<=n;k++) 
    {
        bk2 = 2.0*(v0+k-1.0)*bk1/x+bk0;
        kv[k] = bk2;
        bk0 = bk1;
        bk1 = bk2;
    }
    kvp[0] = v0*kv[0]/x-kv[1];
    for (k=1;k<=n;k++) 
    {
        kvp[k] = -(k+v0)*kv[k]/x-kv[k-1];
    }
    vm = n+v0;
    return 0;
}
 
} // namespace detail

    /*! Bessel function of the first kind. 

        Computes the value of BesselJ of order <tt>n</tt> and argument <tt>x</tt>.
        Negative <tt>x</tt> are currently unsupported and will result in an exception.

        This function implements the algorithm from<br>
        Zhang and Jin: "Computation of Special Functions", John Wiley and Sons, 1996.

        <b>\#include</b> \<vigra/bessel.hxx\><br>
        Namespace: vigra
    */
inline double besselj(double n, double x)
{
    vigra_precondition(x >= 0.0,
        "besselj(): negative x are currently unsupported, sorry.");
        
    if(x == 0.0)
        return n == 0.0 ? 1.0 : 0.0;

    double an = abs(n);
    int in = (int)an, s = in+1;
    
    ArrayVector<double> tmp(4*s);
    double * t = tmp.begin();
    
    vigra_postcondition(
		detail::bessjyv(an, x, an, t, t+s, t+2*s, t+3*s) == 0,
       "besselj(): internal error.");
    if(n < 0.0)
        return tmp[in] * cos_pi(an) - tmp[in+s] * sin_pi(an);
    else
        return tmp[in];
}

    /*! Bessel function of the second kind. 

        Computes the value of BesselY of order <tt>n</tt> and argument <tt>x</tt>.
        Negative <tt>x</tt> are currently unsupported and will result in an exception.

        This function implements the algorithm from<br>
        Zhang and Jin: "Computation of Special Functions", John Wiley and Sons, 1996.

        <b>\#include</b> \<vigra/bessel.hxx\><br>
        Namespace: vigra
    */
inline double bessely(double n, double x)
{
    vigra_precondition(x >= 0.0,
        "bessely(): negative x are currently unsupported, sorry.");
        
    if(x == 0.0)
        return -std::numeric_limits<double>::infinity();
        
    double an = abs(n);
    int in = (int)an, s = in+1;
    
    ArrayVector<double> tmp(4*s);
    double * t = tmp.begin();
    
    vigra_postcondition(
		detail::bessjyv(an, x, an, t, t+s, t+2*s, t+3*s) == 0,
       "bessely(): internal error.");
    if(n < 0.0)
        return tmp[in] * sin_pi(an) + tmp[in+s] * cos_pi(an);
    else
        return tmp[in+s];
}

    /*! Modified Bessel function of the first kind. 

        Computes the value of BesselI of order <tt>n</tt> and argument <tt>x</tt>.
        Negative <tt>x</tt> are currently unsupported and will result in an exception.

        This function implements the algorithm from<br>
        Zhang and Jin: "Computation of Special Functions", John Wiley and Sons, 1996.

        <b>\#include</b> \<vigra/bessel.hxx\><br>
        Namespace: vigra
    */
inline double besseli(double n, double x)
{
    vigra_precondition(x >= 0.0,
        "besseli(): negative x are currently unsupported, sorry.");
        
    if(x == 0.0)
        return n == 0.0 ? 1.0 : 0.0;

        double an = abs(n);
    int in = (int)an, s = in+1;
    
    ArrayVector<double> tmp(4*s);
    double * t = tmp.begin();
    
    vigra_postcondition(
		detail::bessikv(an, x, an, t, t+s, t+2*s, t+3*s) == 0,
       "besseli(): internal error.");
    if(n < 0.0)
        return tmp[in] + M_2_PI * tmp[in+s] * sin_pi(an);
    else
        return tmp[in];
}

    /*! Modified Bessel function of the second kind. 

        Computes the value of BesselK of order <tt>n</tt> and argument <tt>x</tt>.
        Negative <tt>x</tt> are currently unsupported and will result in an exception.

        This function implements the algorithm from<br>
        Zhang and Jin: "Computation of Special Functions", John Wiley and Sons, 1996.

        <b>\#include</b> \<vigra/bessel.hxx\><br>
        Namespace: vigra
    */
inline double besselk(double n, double x)
{
     vigra_precondition(x >= 0.0,
        "besselk(): negative x are currently unsupported, sorry.");

    if(x == 0.0)
        return std::numeric_limits<double>::infinity();
        
    double an = abs(n);
    int in = (int)an, s = in+1;
    
    ArrayVector<double> tmp(4*s);
    double * t = tmp.begin();
    
    vigra_postcondition(
		detail::bessikv(an, x, an, t, t+s, t+2*s, t+3*s) == 0,
       "besselk(): internal error.");
    return tmp[in+s];
}

//@}

} // namespace vigra

#endif // VIGRA_BESSEL_HXX