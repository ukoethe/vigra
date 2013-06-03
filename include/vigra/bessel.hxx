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

#ifndef VIGRA_BESSEL_HXX
#define VIGRA_BESSEL_HXX

#include "mathutil.hxx"

#ifdef HasBoostMath
#include <boost/math/special_functions/bessel.hpp>
#endif

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

//
//  INPUT:
//  double x    -- argument of Bessel function of 1st and 2nd kind.
//  int n       -- order
//
//  OUPUT:
//
//  int nm      -- highest order actually computed (nm <= n)
//  double jn[] -- Bessel function of 1st kind, orders from 0 to nm
//  double yn[] -- Bessel function of 2nd kind, orders from 0 to nm
//
//  Computes Bessel functions of all order up to 'n' using recurrence
//  relations. If 'nm' < 'n' only 'nm' orders are returned.
//
// code has been adapted from C.R. Bond's implementation
// see http://www.crbond.com/math.htm
//
template <class REAL>
void bessjyn(int n, REAL x,int &nm, double *jn, double *yn)
{
    double t1,t2,f,f1,f2,bj0,bj1,bjk,by0,by1,cu,s0,su,sv;
    double ec,bs,byk,p0,p1,q0,q1;
    double a[] = {
        -0.7031250000000000e-1,
         0.1121520996093750,
        -0.5725014209747314,
         6.074042001273483};
    double b[] = {
         0.7324218750000000e-1,
        -0.2271080017089844,
         1.727727502584457,
        -2.438052969955606e1};
    double a1[] = {
         0.1171875,
        -0.1441955566406250,
         0.6765925884246826,
        -6.883914268109947};
    double b1[] = {
       -0.1025390625,
        0.2775764465332031,
       -1.993531733751297,
        2.724882731126854e1};
        
    int i,k,m;
    nm = n;
    if (x < 1e-15) 
    {
        for (i=0;i<=n;i++) 
        {
            jn[i] = 0.0;
            yn[i] = -1e308;
        }
        jn[0] = 1.0;
        return;
    }
    if (x <= 300.0 || n > (int)(0.9*x)) 
    {
        if (n == 0) 
            nm = 1;
        m = msta1(x,200);
        if (m < nm) 
            nm = m;
        else 
            m = msta2(x,nm,15);
        bs = 0.0;
        su = 0.0;
        sv = 0.0;
        f2 = 0.0;
        f1 = 1.0e-100;
        for (k = m;k>=0;k--) 
        {
            f = 2.0*(k+1.0)/x*f1 - f2;
            if (k <= nm) 
                jn[k] = f;
            if ((k == 2*(int)(k/2)) && (k != 0)) 
            {
                bs += 2.0*f;
                su += (-1)*((k & 2)-1)*f/(double)k;
            }
            else if (k > 1) 
            {
                sv += (-1)*((k & 2)-1)*(double)k*f/(k*k-1.0);
            }
            f2 = f1;
            f1 = f;
        }
        s0 = bs+f;
        for (k=0;k<=nm;k++) 
        {
            jn[k] /= s0;
        }
        ec = std::log(0.5*x) + M_EULER_GAMMA;
        by0 = M_2_PI*(ec*jn[0]-4.0*su/s0);
        yn[0] = by0;
        by1 = M_2_PI*((ec-1.0)*jn[1]-jn[0]/x-4.0*sv/s0);
        yn[1] = by1;
    }
    else 
    {
        t1 = x-M_PI_4;
        p0 = 1.0;
        q0 = -0.125/x;
        for (k=0;k<4;k++) 
        {
            p0 += a[k]*std::pow(x,-2*k-2);
            q0 += b[k]*std::pow(x,-2*k-3);
        }
        cu = std::sqrt(M_2_PI/x);
        bj0 = cu*(p0*std::cos(t1)-q0*std::sin(t1));
        by0 = cu*(p0*std::sin(t1)+q0*std::cos(t1));
        jn[0] = bj0;
        yn[0] = by0;
        t2 = x-0.75*M_PI;
        p1 = 1.0;
        q1 = 0.375/x;
        for (k=0;k<4;k++) 
        {
            p1 += a1[k]*std::pow(x,-2*k-2);
            q1 += b1[k]*std::pow(x,-2*k-3);
        }
        bj1 = cu*(p1*std::cos(t2)-q1*std::sin(t2));
        by1 = cu*(p1*std::sin(t2)+q1*std::cos(t2));
        jn[1] = bj1;
        yn[1] = by1;
        for (k=2;k<=nm;k++) 
        {
            bjk = 2.0*(k-1.0)*bj1/x-bj0;
            jn[k] = bjk;
            bj0 = bj1;
            bj1 = bjk;
        }
    }
    for (k=2;k<=nm;k++) 
    {
        byk = 2.0*(k-1.0)*by1/x-by0;
        yn[k] = byk;
        by0 = by1;
        by1 = byk;
    }
}


 
} // namespace detail

    /** \brief Bessel function of the first kind. 

        Computes the value of BesselJ of integer order <tt>n</tt> and argument <tt>x</tt>.
        Negative <tt>x</tt> are unsupported and will result in a <tt>std::domain_error</tt>.

        This function wraps a number of existing implementations and falls back to 
        a rather slow algorithm if none of them is available. In particular,
        it uses boost::math when <tt>HasBoostMath</tt> is \#defined, or native 
        implementations on gcc and MSVC otherwise.

        <b>\#include</b> \<vigra/bessel.hxx\><br>
        Namespace: vigra
    */
inline double besselJ(int n, double x)
{
    if(x < 0.0)
        throw std::domain_error("besselJ(n, x): x cannot be negative");
    if(x < 1e-15)
        return n == 0 ? 1.0 : 0.0;
#if defined(HasBoostMath)
    return boost::math::cyl_bessel_j((double)n, x);
#elif defined(__GNUC__)
    return ::jn(n, x);
#elif defined(_MSC_VER)
    return _jn(n, x);
#else
    int an = abs(n), nr = n, s = an+2;
    ArrayVector<double> t(2*s);
    detail::bessjyn(an, x, nr, &t[0], &t[s]);
    if(n < 0 && odd(an))
        return -t[an];
    else
        return  t[an];
#endif
}

    /** \brief Bessel function of the second kind. 

        Computes the value of BesselY of integer order <tt>n</tt> and argument <tt>x</tt>.
        Negative <tt>x</tt> are unsupported and will result in a <tt>std::domain_error</tt>.

        This function wraps a number of existing implementations and falls back to 
        a rather slow algorithm if none of them is available. In particular,
        it uses boost::math when <tt>HasBoostMath</tt> is \#defined, or native 
        implementations on gcc and MSVC otherwise.

        <b>\#include</b> \<vigra/bessel.hxx\><br>
        Namespace: vigra
    */
inline double besselY(int n, double x)
{
    if(x < 0.0)
        throw std::domain_error("besselY(n, x): x cannot be negative");
    if(x == 0.0 )
        return -std::numeric_limits<double>::infinity();
#if defined(HasBoostMath)
    return boost::math::cyl_neumann((double)n, x);
#elif defined(__GNUC__)
    return ::yn(n, x);
#elif defined(_MSC_VER)
    return _yn(n, x);
#else
    int an = abs(n), nr = n, s = an+2;
    ArrayVector<double> t(2*s);
    detail::bessjyn(an, x, nr, &t[0], &t[s]);
    if(n < 0.0 && odd(n))
        return -t[an+s];
    else
        return  t[an+s];
#endif
}

//@}

} // namespace vigra

#endif // VIGRA_BESSEL_HXX
