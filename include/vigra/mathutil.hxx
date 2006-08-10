/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2005 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#ifndef VIGRA_MATHUTIL_HXX
#define VIGRA_MATHUTIL_HXX

#include <cmath>
#include <cstdlib>
#include "vigra/config.hxx"
#include "vigra/sized_int.hxx"
#include "vigra/numerictraits.hxx"

/*! \page MathConstants Mathematical Constants

    <TT>M_PI, M_SQRT2</TT>

    <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"

    Since <TT>M_PI</TT> and <TT>M_SQRT2</TT> are not officially standardized,
    we provide definitions here for those compilers that don't support them.

    \code
    #ifndef M_PI
    #    define M_PI     3.14159265358979323846
    #endif

    #ifndef M_SQRT2
    #    define M_SQRT2  1.41421356237309504880
    #endif
    \endcode
*/
#ifndef M_PI
#    define M_PI     3.14159265358979323846
#endif

#ifndef M_SQRT2
#    define M_SQRT2  1.41421356237309504880
#endif

namespace vigra {

#ifndef __sun 

/** \addtogroup MathFunctions Mathematical Functions

    Useful mathematical functions and functors.
*/
//@{
    /*! The error function.

        If <tt>erf()</tt> is not provided in the C standard math library (as it should according to the
        new C99 standard ?), VIGRA implements <tt>erf()</tt> as an approximation of the error 
        function
        
        \f[
            \mbox{erf}(x) = \int_0^x e^{-t^2} dt
        \f]
        
        according to the formula given in Press et al. "Numerical Recipes".

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
template <class T>
double erf(T x)
{
    double t = 1.0/(1.0+0.5*VIGRA_CSTD::fabs(x));
    double ans = t*VIGRA_CSTD::exp(-x*x-1.26551223+t*(1.00002368+t*(0.37409196+
                                    t*(0.09678418+t*(-0.18628806+t*(0.27886807+
                                    t*(-1.13520398+t*(1.48851587+t*(-0.82215223+
                                    t*0.17087277)))))))));
    if (x >= 0.0)
        return 1.0 - ans;
    else
        return ans - 1.0;
}

#else

using VIGRA_CSTD::erf;

#endif

// import functions into namespace vigra which VIGRA is going to overload

using VIGRA_CSTD::pow;  
using VIGRA_CSTD::floor;  
using VIGRA_CSTD::ceil;  

// import abs(float), abs(double), abs(long double) from <cmath>
//    and abs(int), abs(long), abs(long long) from <cstdlib>
using std::abs;  

// define the missing variants of abs() to avoid 'ambigous overload'
// errors in template functions
#define VIGRA_DEFINE_UNSIGNED_ABS(T) \
    inline T abs(T t) { return t; }

VIGRA_DEFINE_UNSIGNED_ABS(bool)
VIGRA_DEFINE_UNSIGNED_ABS(unsigned char)
VIGRA_DEFINE_UNSIGNED_ABS(unsigned short)
VIGRA_DEFINE_UNSIGNED_ABS(unsigned int)
VIGRA_DEFINE_UNSIGNED_ABS(unsigned long)
VIGRA_DEFINE_UNSIGNED_ABS(unsigned long long)

#undef VIGRA_DEFINE_UNSIGNED_ABS

#define VIGRA_DEFINE_MISSING_ABS(T) \
    inline T abs(T t) { return t < 0 ? -t : t; }

VIGRA_DEFINE_MISSING_ABS(signed char)
VIGRA_DEFINE_MISSING_ABS(signed short)

#undef VIGRA_DEFINE_MISSING_ABS

    /*! The rounding function.

        Defined for all floating point types. Rounds towards the nearest integer 
        such that <tt>abs(round(t)) == round(abs(t))</tt> for all <tt>t</tt>.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
inline float round(float t)
{
     return t >= 0.0
                ? floor(t + 0.5)
                : ceil(t - 0.5);
}

inline double round(double t)
{
     return t >= 0.0
                ? floor(t + 0.5)
                : ceil(t - 0.5);
}

inline long double round(long double t)
{
     return t >= 0.0
                ? floor(t + 0.5)
                : ceil(t - 0.5);
}

    /*! The square function.

        sq(x) is needed so often that it makes sense to define it as a function.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
template <class T>
inline 
typename NumericTraits<T>::Promote sq(T t)
{
    return t*t;
}

namespace detail {

template <class T>
class IntSquareRoot
{
  public:
    static int sqq_table[];
    static UInt32 exec(UInt32 v);
};

template <class T>
int IntSquareRoot<T>::sqq_table[] = {
	       0,  16,  22,  27,  32,  35,  39,  42,  45,  48,  50,  53,  55,  57,
	      59,  61,  64,  65,  67,  69,  71,  73,  75,  76,  78,  80,  81,  83,
	      84,  86,  87,  89,  90,  91,  93,  94,  96,  97,  98,  99, 101, 102,
	     103, 104, 106, 107, 108, 109, 110, 112, 113, 114, 115, 116, 117, 118,
	     119, 120, 121, 122, 123, 124, 125, 126, 128, 128, 129, 130, 131, 132,
	     133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 144, 145,
	     146, 147, 148, 149, 150, 150, 151, 152, 153, 154, 155, 155, 156, 157,
	     158, 159, 160, 160, 161, 162, 163, 163, 164, 165, 166, 167, 167, 168,
	     169, 170, 170, 171, 172, 173, 173, 174, 175, 176, 176, 177, 178, 178,
	     179, 180, 181, 181, 182, 183, 183, 184, 185, 185, 186, 187, 187, 188,
	     189, 189, 190, 191, 192, 192, 193, 193, 194, 195, 195, 196, 197, 197,
	     198, 199, 199, 200, 201, 201, 202, 203, 203, 204, 204, 205, 206, 206,
	     207, 208, 208, 209, 209, 210, 211, 211, 212, 212, 213, 214, 214, 215,
	     215, 216, 217, 217, 218, 218, 219, 219, 220, 221, 221, 222, 222, 223,
	     224, 224, 225, 225, 226, 226, 227, 227, 228, 229, 229, 230, 230, 231,
	     231, 232, 232, 233, 234, 234, 235, 235, 236, 236, 237, 237, 238, 238,
	     239, 240, 240, 241, 241, 242, 242, 243, 243, 244, 244, 245, 245, 246,
	     246, 247, 247, 248, 248, 249, 249, 250, 250, 251, 251, 252, 252, 253,
	     253, 254, 254, 255
};

template <class T>
UInt32 IntSquareRoot<T>::exec(UInt32 x) 
{
    unsigned long xn;
	if (x >= 0x10000)
	    if (x >= 0x1000000)
	        if (x >= 0x10000000)
	            if (x >= 0x40000000) {
	                if (x >= (UInt32)65535*(UInt32)65535)
	                    return 65535;
	                xn = sqq_table[x>>24] << 8;
	            } else
	                xn = sqq_table[x>>22] << 7;
	        else
	            if (x >= 0x4000000)
	                xn = sqq_table[x>>20] << 6;
	            else
	                xn = sqq_table[x>>18] << 5;
	    else {
	        if (x >= 0x100000)
	            if (x >= 0x400000)
	                xn = sqq_table[x>>16] << 4;
	            else
	                xn = sqq_table[x>>14] << 3;
	        else
	            if (x >= 0x40000)
	                xn = sqq_table[x>>12] << 2;
	            else
	                xn = sqq_table[x>>10] << 1;

	        goto nr1;
	    }
	else
	    if (x >= 0x100) {
	        if (x >= 0x1000)
	            if (x >= 0x4000)
	                xn = (sqq_table[x>>8] >> 0) + 1;
	            else
	                xn = (sqq_table[x>>6] >> 1) + 1;
	        else
	            if (x >= 0x400)
	                xn = (sqq_table[x>>4] >> 2) + 1;
	            else
	                xn = (sqq_table[x>>2] >> 3) + 1;

	        goto adj;
	    } else
	        return sqq_table[x] >> 4;

    /* Run two iterations of the standard convergence formula */

	xn = (xn + 1 + x / xn) / 2;
  nr1:
	xn = (xn + 1 + x / xn) / 2;
  adj:

	if (xn * xn > x) /* Correct rounding if necessary */
	    xn--;

	return xn;
}

} // namespace detail

using VIGRA_CSTD::sqrt;

    /*! Signed integer square root.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
inline Int32 sqrti(Int32 v)
{
	if(v < 0)
	    throw std::domain_error("sqrti(Int32): negative argument.");
    return (Int32)detail::IntSquareRoot<UInt32>::exec((UInt32)v);
}

    /*! Unsigned integer square root.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
inline UInt32 sqrti(UInt32 v)
{
    return detail::IntSquareRoot<UInt32>::exec(v);
}

#ifdef VIGRA_NO_HYPOT
    /*! Compute the Euclidean distance (length of the hypothenuse of a right-angled triangle).

        The  hypot()  function  returns  the  sqrt(a*a  +  b*b).
        It is implemented in a way that minimizes round-off error.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
inline double hypot(double a, double b) 
{ 
    double absa = VIGRA_CSTD::fabs(a), absb = VIGRA_CSTD::fabs(b);
    if (absa > absb) 
        return absa * VIGRA_CSTD::sqrt(1.0 + sq(absb/absa)); 
    else 
        return absb == 0.0
                   ? 0.0
                   : absb * VIGRA_CSTD::sqrt(1.0 + sq(absa/absb)); 
}

#else

using ::hypot;

#endif

namespace detail {

template <class T>
T ellipticRD(T x, T y, T z)
{
    double f = 1.0, s = 0.0, X, Y, Z, m;
    while(true)
    {
        m = (x + y + 3.0*z) / 5.0;
        X = 1.0 - x/m;
        Y = 1.0 - y/m;
        Z = 1.0 - z/m;
        if(std::max(std::max(std::fabs(X), std::fabs(Y)), std::fabs(Z)) < 0.01)
            break;
        double l = std::sqrt(x*y) + std::sqrt(x*z) + std::sqrt(y*z);
        s += f / (std::sqrt(z)*(z + l));
        f /= 4.0;
        x = (x + l)/4.0;
        y = (y + l)/4.0;
        z = (z + l)/4.0;
    }
    double a = X*Y;
    double b = sq(Z);
    double c = a - b;
    double d = a - 6.0*b;
    double e = d + 2.0*c;
    return 3.0*s + f*(1.0+d*(-3.0/14.0+d*9.0/88.0-Z*e*4.5/26.0)
                      +Z*(e/6.0+Z*(-c*9.0/22.0+a*Z*3.0/26.0))) / std::pow(m,1.5);
}

template <class T>
T ellipticRF(T x, T y, T z)
{
    double X, Y, Z, m;
    while(true)
    {
        m = (x + y + z) / 3.0;
        X = 1.0 - x/m;
        Y = 1.0 - y/m;
        Z = 1.0 - z/m;
        if(std::max(std::max(std::fabs(X), std::fabs(Y)), std::fabs(Z)) < 0.01)
            break;
        double l = std::sqrt(x*y) + std::sqrt(x*z) + std::sqrt(y*z);
        x = (x + l)/4.0;
        y = (y + l)/4.0;
        z = (z + l)/4.0;
    }
    double d = X*Y - sq(Z);
    double p = X*Y*Z;
    return (1.0 - d/10.0 + p/14.0 + sq(d)/24.0 - d*p*3.0/44.0) / std::sqrt(m);
}

} // namespace detail

    /*! The incomplete elliptic integral of the first kind.

        Computes
        
        \f[
            \mbox{F}(x, k) = \int_0^x \frac{1}{\sqrt{1 - k^2 \sin(t)^2}} dt
        \f]
        
        according to the algorithm given in Press et al. "Numerical Recipes". The
        complete elliptic integral of the first kind is simply <tt>ellipticIntegralF(M_PI/2, k)</TT>.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
inline double ellipticIntegralF(double x, double k)
{
    double c2 = sq(std::cos(x));
    double s = std::sin(x);
    return s*detail::ellipticRF(c2, 1.0 - sq(k*s), 1.0);
}

    /*! The incomplete elliptic integral of the second kind.

        Computes
        
        \f[
            \mbox{E}(x, k) = \int_0^x \sqrt{1 - k^2 \sin(t)^2} dt
        \f]
        
        according to the algorithm given in Press et al. "Numerical Recipes". The
        complete elliptic integral of the second kind is simply <tt>ellipticIntegralE(M_PI/2, k)</TT>.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
inline double ellipticIntegralE(double x, double k)
{
    double c2 = sq(std::cos(x));
    double s = std::sin(x);
    k = sq(k*s);
    return s*(detail::ellipticRF(c2, 1.0-k, 1.0) - k/3.0*detail::ellipticRD(c2, 1.0-k, 1.0));
}

    /*! The sign function.

        Returns 1, 0, or -1 depending on the sign of \a t.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
template <class T>
T sign(T t) 
{ 
    return t > NumericTraits<T>::zero()
               ? NumericTraits<T>::one()
               : t < NumericTraits<T>::zero()
                    ? -NumericTraits<T>::one()
                    : NumericTraits<T>::zero();
}

    /*! The binary sign function.

        Transfers the sign of \a t2 to \a t1.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
template <class T1, class T2>
T1 sign(T1 t1, T2 t2) 
{ 
    return t2 >= NumericTraits<T2>::zero()
               ? abs(t1)
               : -abs(t1);
}

#define VIGRA_DEFINE_NORM(T) \
    inline NormTraits<T>::SquaredNormType squaredNorm(T t) { return sq(t); } \
    inline NormTraits<T>::NormType norm(T t) { return abs(t); }

VIGRA_DEFINE_NORM(bool)
VIGRA_DEFINE_NORM(signed char)
VIGRA_DEFINE_NORM(unsigned char)
VIGRA_DEFINE_NORM(short)
VIGRA_DEFINE_NORM(unsigned short)
VIGRA_DEFINE_NORM(int)
VIGRA_DEFINE_NORM(unsigned int)
VIGRA_DEFINE_NORM(long)
VIGRA_DEFINE_NORM(unsigned long)
VIGRA_DEFINE_NORM(float)
VIGRA_DEFINE_NORM(double)
VIGRA_DEFINE_NORM(long double)

#undef VIGRA_DEFINE_NORM

template <class T>
inline typename NormTraits<std::complex<T> >::SquaredNormType
squaredNorm(std::complex<T> const & t)
{
    return sq(t.real()) + sq(t.imag());
}

#ifdef DOXYGEN // only for documentation
    /*! The squared norm of a numerical object.

        For scalar types: equals <tt>vigra::sq(t)</tt><br>.
        For vectorial types: equals <tt>vigra::dot(t, t)</tt><br>.
        For complex types: equals <tt>vigra::sq(t.real()) + vigra::sq(t.imag())</tt><br>.
        For matrix types: results in the squared Frobenius norm (sum of squares of the matrix elements).
    */
NormTraits<T>::SquaredNormType squaredNorm(T const & t);

#endif

    /*! The norm of a numerical object.

        For scalar types: implemented as <tt>abs(t)</tt><br>
        otherwise: implemented as <tt>sqrt(squaredNorm(t))</tt>.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
template <class T>
inline typename NormTraits<T>::NormType 
norm(T const & t)
{
    typedef typename NormTraits<T>::SquaredNormType SNT;
    return sqrt(static_cast<typename SquareRootTraits<SNT>::SquareRootArgument>(squaredNorm(t)));
}

namespace detail {

// both f1 and f2 are unsigned here
template<class FPT>
inline
FPT safeFloatDivision( FPT f1, FPT f2 )
{
    return  f2 < NumericTraits<FPT>::one() && f1 > f2 * NumericTraits<FPT>::max()
                ? NumericTraits<FPT>::max() 
                : (f2 > NumericTraits<FPT>::one() && f1 < f2 * NumericTraits<FPT>::smallestPositive()) || 
                   f1 == NumericTraits<FPT>::zero()
                     ? NumericTraits<FPT>::zero() 
                     : f1/f2;
}

} // namespace detail
    
    /*! Tolerance based floating-point comparison.

        Check whether two floating point numbers are equal within the given tolerance.
        This is useful because floating point numbers that should be equal in theory are
        rarely exactly equal in practice. If the tolerance \a epsilon is not given,
        twice the machine epsilon is used.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
template <class T1, class T2>
bool closeAtTolerance(T1 l, T2 r, typename PromoteTraits<T1, T2>::Promote epsilon)
{
    typedef typename PromoteTraits<T1, T2>::Promote T;
    if(l == 0.0)
        return VIGRA_CSTD::fabs(r) <= epsilon;
    if(r == 0.0)
        return VIGRA_CSTD::fabs(l) <= epsilon;
    T diff = VIGRA_CSTD::fabs( l - r );
    T d1   = detail::safeFloatDivision<T>( diff, VIGRA_CSTD::fabs( r ) );
    T d2   = detail::safeFloatDivision<T>( diff, VIGRA_CSTD::fabs( l ) );

    return (d1 <= epsilon && d2 <= epsilon);
}

template <class T1, class T2>
bool closeAtTolerance(T1 l, T2 r)
{
    typedef typename PromoteTraits<T1, T2>::Promote T;
    return closeAtTolerance(l, r, 2.0 * NumericTraits<T>::epsilon());
}

//@}

} // namespace vigra

#endif /* VIGRA_MATHUTIL_HXX */
