/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2005 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_MATHUTIL_HXX
#define VIGRA_MATHUTIL_HXX

#include <cmath>
#include <cstdlib>
#include "vigra/config.hxx"
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
            \mbox{erf}(x) = \int_0^x e^{-x^2} dx
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
using std::abs;  

#define VIGRA_DEFINE_UNSIGNED_ABS(T) \
    inline T abs(T t) { return t; }

VIGRA_DEFINE_UNSIGNED_ABS(bool)
VIGRA_DEFINE_UNSIGNED_ABS(unsigned char)
VIGRA_DEFINE_UNSIGNED_ABS(unsigned short)
VIGRA_DEFINE_UNSIGNED_ABS(unsigned int)
VIGRA_DEFINE_UNSIGNED_ABS(unsigned long)

#undef VIGRA_DEFINE_UNSIGNED_ABS

    /*! The rounding function.

        Defined for all floating point types. Rounds towards the nearest integer for both 
        positive and negative inputs.

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

#ifdef VIGRA_NO_HYPOT
    /*! Compute the Euclidean distance (length of the hypothenuse of a right-angled triangle).

        The  hypot()  function  returns  the  sqrt(a*a  +  b*b).
        It is implemented in a way that minimizes round-off error.

        <b>\#include</b> "<a href="mathutil_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
        Namespace: vigra
    */
template <class T>
T hypot(T a, T b) 
{ 
    T absa = abs(a), absb = abs(b);
    if (absa > absb) 
        return absa * VIGRA_CSTD::sqrt(1.0 + sq(absb/absa)); 
    else 
        return absb == NumericTraits<T>::zero()
                   ? NumericTraits<T>::zero()
                   : absb * VIGRA_CSTD::sqrt(1.0 + sq(absa/absb)); 
}

#else

using ::hypot;

#endif

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
    return VIGRA_CSTD::sqrt(static_cast<typename NormTraits<T>::NormType>(squaredNorm(t)));
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
    if(l == 0.0 && r != 0.0)
        return VIGRA_CSTD::fabs(r) <= epsilon;
    if(l != 0.0 && r == 0.0)
        return VIGRA_CSTD::fabs(r) <= epsilon;
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
