/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2011 by Ullrich Koethe                  */
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

#ifndef VIGRA_MATHUTIL_HXX
#define VIGRA_MATHUTIL_HXX

#ifdef _MSC_VER
# pragma warning (disable: 4996) // hypot/_hypot confusion
#endif

#include <cmath>
#include <cstdlib>
#include <complex>
#include "config.hxx"
#include "error.hxx"
#include "tuple.hxx"
#include "sized_int.hxx"
#include "numerictraits.hxx"
#include "algorithm.hxx"

/** \page MathConstants Mathematical Constants

    <TT>M_PI, M_SQRT2 etc.</TT>

    <b>\#include</b> \<vigra/mathutil.hxx\>

    Since mathematical constants such as <TT>M_PI</TT> and <TT>M_SQRT2</TT>
    are not officially standardized, we provide definitions here for those
    compilers that don't support them.

    \code
    #ifndef M_PI
    #    define M_PI     3.14159265358979323846
    #endif

    #ifndef M_SQRT2
    #    define M_2_PI   0.63661977236758134308
    #endif

    #ifndef M_PI_2
    #    define M_PI_2   1.57079632679489661923
    #endif

    #ifndef M_PI_4
    #    define M_PI_4   0.78539816339744830962
    #endif

    #ifndef M_SQRT2
    #    define M_SQRT2  1.41421356237309504880
    #endif

    #ifndef M_EULER_GAMMA
    #    define M_EULER_GAMMA  0.5772156649015329
    #endif
    \endcode
*/
#ifndef M_PI
#    define M_PI     3.14159265358979323846
#endif

#ifndef M_2_PI
#    define M_2_PI   0.63661977236758134308
#endif

#ifndef M_PI_2
#    define M_PI_2   1.57079632679489661923
#endif

#ifndef M_PI_4
#    define M_PI_4   0.78539816339744830962
#endif

#ifndef M_SQRT2
#    define M_SQRT2  1.41421356237309504880
#endif

#ifndef M_E
#    define M_E      2.71828182845904523536
#endif

#ifndef M_EULER_GAMMA
#    define M_EULER_GAMMA  0.5772156649015329
#endif

namespace vigra {

/** \addtogroup MathFunctions Mathematical Functions

    Useful mathematical functions and functors.
*/
//@{
// import functions into namespace vigra which VIGRA is going to overload

using VIGRA_CSTD::pow;
using VIGRA_CSTD::floor;
using VIGRA_CSTD::ceil;
using VIGRA_CSTD::exp;

// import abs(float), abs(double), abs(long double) from <cmath>
//        abs(int), abs(long), abs(long long) from <cstdlib>
//        abs(std::complex<T>) from <complex>
using std::abs;

// define the missing variants of abs() to avoid 'ambiguous overload'
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
    inline T abs(T t) { return t < 0 ? static_cast<T>(-t) : t; }

VIGRA_DEFINE_MISSING_ABS(signed char)
VIGRA_DEFINE_MISSING_ABS(signed short)

#if defined(_MSC_VER) && _MSC_VER < 1600
VIGRA_DEFINE_MISSING_ABS(signed long long)
#endif

#undef VIGRA_DEFINE_MISSING_ABS

#ifndef _MSC_VER

using std::isinf;
using std::isnan;

#else

template <class REAL>
inline bool isinf(REAL v)
{
    return _finite(v) == 0;
}

template <class REAL>
inline bool isnan(REAL v)
{
    return _isnan(v) != 0;
}

#endif

// scalar dot is needed for generic functions that should work with
// scalars and vectors alike

#define VIGRA_DEFINE_SCALAR_DOT(T) \
    inline NumericTraits<T>::Promote dot(T l, T r) { return l*r; }

VIGRA_DEFINE_SCALAR_DOT(unsigned char)
VIGRA_DEFINE_SCALAR_DOT(unsigned short)
VIGRA_DEFINE_SCALAR_DOT(unsigned int)
VIGRA_DEFINE_SCALAR_DOT(unsigned long)
VIGRA_DEFINE_SCALAR_DOT(unsigned long long)
VIGRA_DEFINE_SCALAR_DOT(signed char)
VIGRA_DEFINE_SCALAR_DOT(signed short)
VIGRA_DEFINE_SCALAR_DOT(signed int)
VIGRA_DEFINE_SCALAR_DOT(signed long)
VIGRA_DEFINE_SCALAR_DOT(signed long long)
VIGRA_DEFINE_SCALAR_DOT(float)
VIGRA_DEFINE_SCALAR_DOT(double)
VIGRA_DEFINE_SCALAR_DOT(long double)

#undef VIGRA_DEFINE_SCALAR_DOT

using std::pow;

// support 'double' exponents for all floating point versions of pow()

inline float pow(float v, double e)
{
    return std::pow(v, (float)e);
}

inline long double pow(long double v, double e)
{
    return std::pow(v, (long double)e);
}

    /** \brief The rounding function.

        Defined for all floating point types. Rounds towards the nearest integer
        such that <tt>abs(round(t)) == round(abs(t))</tt> for all <tt>t</tt>.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
#ifdef DOXYGEN // only for documentation
REAL round(REAL v);
#endif

inline float round(float t)
{
     return t >= 0.0
                ? floor(t + 0.5f)
                : ceil(t - 0.5f);
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


    /** \brief Round and cast to integer.

        Rounds to the nearest integer like round(), but casts the result to
        <tt>long long</tt> (this will be faster and is usually needed anyway).

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline long long roundi(double t)
{
     return t >= 0.0
                ? (long long)(t + 0.5)
                : (long long)(t - 0.5);
}

    /** \brief Round up to the nearest power of 2.

        Efficient algorithm for finding the smallest power of 2 which is not smaller than \a x
        (function clp2() from Henry Warren: "Hacker's Delight", Addison-Wesley, 2003,
         see http://www.hackersdelight.org/).
        If \a x > 2^31, the function will return 0 because integer arithmetic is defined modulo 2^32.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline UInt32 ceilPower2(UInt32 x)
{
    if(x == 0) return 0;

    x = x - 1;
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    x = x | (x >> 8);
    x = x | (x >>16);
    return x + 1;
}

    /** \brief Round down to the nearest power of 2.

        Efficient algorithm for finding the largest power of 2 which is not greater than \a x
        (function flp2() from Henry Warren: "Hacker's Delight", Addison-Wesley, 2003,
         see http://www.hackersdelight.org/).

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline UInt32 floorPower2(UInt32 x)
{
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    x = x | (x >> 8);
    x = x | (x >>16);
    return x - (x >> 1);
}

namespace detail {

template <class T>
struct IntLog2
{
    static Int32 table[64];
};

template <class T>
Int32 IntLog2<T>::table[64] = {
         -1,  0,  -1,  15,  -1,  1,  28,  -1,  16,  -1,  -1,  -1,  2,  21,
         29,  -1,  -1,  -1,  19,  17,  10,  -1,  12,  -1,  -1,  3,  -1,  6,
         -1,  22,  30,  -1,  14,  -1,  27,  -1,  -1,  -1,  20,  -1,  18,  9,
         11,  -1,  5,  -1,  -1,  13,  26,  -1,  -1,  8,  -1,  4,  -1,  25,
         -1,  7,  24,  -1,  23,  -1,  31,  -1};

} // namespace detail

    /** \brief Compute the base-2 logarithm of an integer.

        Returns the position of the left-most 1-bit in the given number \a x, or
        -1 if \a x == 0. That is,

        \code
        assert(k >= 0 && k < 32 && log2i(1 << k) == k);
        \endcode

        The function uses Robert Harley's algorithm to determine the number of leading zeros
        in \a x (algorithm nlz10() at http://www.hackersdelight.org/). But note that the functions
        \ref floorPower2() or \ref ceilPower2() are more efficient and should be preferred when possible.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline Int32 log2i(UInt32 x)
{
    // Propagate leftmost 1-bit to the right.
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    x = x | (x >> 8);
    x = x | (x >>16);
    x = x*0x06EB14F9; // Multiplier is 7*255**3.
    return detail::IntLog2<Int32>::table[x >> 26];
}

    /** \brief The square function.

        <tt>sq(x) = x*x</tt> is needed so often that it makes sense to define it as a function.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T>
inline
typename NumericTraits<T>::Promote sq(T t)
{
    return t*t;
}

namespace detail {

template <class V, unsigned>
struct cond_mult
{
    static V call(const V & x, const V & y) { return x * y; }
};
template <class V>
struct cond_mult<V, 0>
{
    static V call(const V &, const V & y) { return y; }
};

template <class V, unsigned n>
struct power_static
{
    static V call(const V & x)
    {
        return n / 2
            ? cond_mult<V, n & 1>::call(x, power_static<V, n / 2>::call(x * x))
            : n & 1 ? x : V();
    }
};
template <class V>
struct power_static<V, 0>
{
    static V call(const V & /* x */)
    {
        return V(1);
    }
};

} // namespace detail

    /** \brief Exponentiation to a positive integer power by squaring.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <unsigned n, class V>
inline V power(const V & x)
{
    return detail::power_static<V, n>::call(x);
}
//doxygen_overloaded_function(template <unsigned n, class V> power(const V & x))

namespace detail {

template <class T>
struct IntSquareRoot
{
    static UInt32 sqq_table[];
    static UInt32 exec(UInt32 v);
};

template <class T>
UInt32 IntSquareRoot<T>::sqq_table[] = {
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
    UInt32 xn;
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

    /** \brief Signed integer square root.

        Useful for fast fixed-point computations.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline Int32 sqrti(Int32 v)
{
    if(v < 0)
        throw std::domain_error("sqrti(Int32): negative argument.");
    return (Int32)detail::IntSquareRoot<UInt32>::exec((UInt32)v);
}

    /** \brief Unsigned integer square root.

        Useful for fast fixed-point computations.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline UInt32 sqrti(UInt32 v)
{
    return detail::IntSquareRoot<UInt32>::exec(v);
}

#ifdef VIGRA_NO_HYPOT
    /** \brief Compute the Euclidean distance (length of the hypotenuse of a right-angled triangle).

        The  hypot()  function  returns  the  sqrt(a*a  +  b*b).
        It is implemented in a way that minimizes round-off error.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
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

    /** \brief The sign function.

        Returns 1, 0, or -1 depending on the sign of \a t, but with the same type as \a t.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T>
inline T sign(T t)
{
    return t > NumericTraits<T>::zero()
               ? NumericTraits<T>::one()
               : t < NumericTraits<T>::zero()
                    ? -NumericTraits<T>::one()
                    : NumericTraits<T>::zero();
}

    /** \brief The integer sign function.

        Returns 1, 0, or -1 depending on the sign of \a t, converted to int.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T>
inline int signi(T t)
{
    return t > NumericTraits<T>::zero()
               ? 1
               : t < NumericTraits<T>::zero()
                    ? -1
                    : 0;
}

    /** \brief The binary sign function.

        Transfers the sign of \a t2 to \a t1.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T1, class T2>
inline T1 sign(T1 t1, T2 t2)
{
    return t2 >= NumericTraits<T2>::zero()
               ? abs(t1)
               : -abs(t1);
}


#ifdef DOXYGEN // only for documentation
    /** \brief Check if an integer is even.

        Defined for all integral types.
    */
bool even(int t);

    /** \brief Check if an integer is odd.

        Defined for all integral types.
    */
bool odd(int t);

#endif

#define VIGRA_DEFINE_ODD_EVEN(T) \
    inline bool even(T t) { return (t&1) == 0; } \
    inline bool odd(T t)  { return (t&1) == 1; }

VIGRA_DEFINE_ODD_EVEN(char)
VIGRA_DEFINE_ODD_EVEN(short)
VIGRA_DEFINE_ODD_EVEN(int)
VIGRA_DEFINE_ODD_EVEN(long)
VIGRA_DEFINE_ODD_EVEN(long long)
VIGRA_DEFINE_ODD_EVEN(unsigned char)
VIGRA_DEFINE_ODD_EVEN(unsigned short)
VIGRA_DEFINE_ODD_EVEN(unsigned int)
VIGRA_DEFINE_ODD_EVEN(unsigned long)
VIGRA_DEFINE_ODD_EVEN(unsigned long long)

#undef VIGRA_DEFINE_ODD_EVEN

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
VIGRA_DEFINE_NORM(long long)
VIGRA_DEFINE_NORM(unsigned long long)
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
    /** \brief The squared norm of a numerical object.

        <ul>
        <li>For scalar types: equals <tt>vigra::sq(t)</tt>.
        <li>For vectorial types (including TinyVector): equals <tt>vigra::dot(t, t)</tt>.
        <li>For complex number types: equals <tt>vigra::sq(t.real()) + vigra::sq(t.imag())</tt>.
        <li>For array and matrix types: results in the squared Frobenius norm (sum of squares of the matrix elements).
        </ul>
    */
NormTraits<T>::SquaredNormType squaredNorm(T const & t);

#endif

    /** \brief The norm of a numerical object.

        For scalar types: implemented as <tt>abs(t)</tt><br>
        otherwise: implemented as <tt>sqrt(squaredNorm(t))</tt>.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T>
inline typename NormTraits<T>::NormType
norm(T const & t)
{
    typedef typename NormTraits<T>::SquaredNormType SNT;
    return sqrt(static_cast<typename SquareRootTraits<SNT>::SquareRootArgument>(squaredNorm(t)));
}

    /** \brief Compute the eigenvalues of a 2x2 real symmetric matrix.

        This uses the analytical eigenvalue formula
        \f[
           \lambda_{1,2} = \frac{1}{2}\left(a_{00} + a_{11} \pm \sqrt{(a_{00} - a_{11})^2 + 4 a_{01}^2}\right)
        \f]

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T>
void symmetric2x2Eigenvalues(T a00, T a01, T a11, T * r0, T * r1)
{
    double d  = hypot(a00 - a11, 2.0*a01);
    *r0 = static_cast<T>(0.5*(a00 + a11 + d));
    *r1 = static_cast<T>(0.5*(a00 + a11 - d));
    if(*r0 < *r1)
        std::swap(*r0, *r1);
}

    /** \brief Compute the eigenvalues of a 3x3 real symmetric matrix.

        This uses a numerically stable version of the analytical eigenvalue formula according to
        <p>
        David Eberly: <a href="http://www.geometrictools.com/Documentation/EigenSymmetric3x3.pdf">
        <em>"Eigensystems for 3 × 3 Symmetric Matrices (Revisited)"</em></a>, Geometric Tools Documentation, 2006

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T>
void symmetric3x3Eigenvalues(T a00, T a01, T a02, T a11, T a12, T a22,
                             T * r0, T * r1, T * r2)
{
    double inv3 = 1.0 / 3.0, root3 = std::sqrt(3.0);

    double c0 = a00*a11*a22 + 2.0*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01;
    double c1 = a00*a11 - a01*a01 + a00*a22 - a02*a02 + a11*a22 - a12*a12;
    double c2 = a00 + a11 + a22;
    double c2Div3 = c2*inv3;
    double aDiv3 = (c1 - c2*c2Div3)*inv3;
    if (aDiv3 > 0.0)
        aDiv3 = 0.0;
    double mbDiv2 = 0.5*(c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1));
    double q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3;
    if (q > 0.0)
        q = 0.0;
    double magnitude = std::sqrt(-aDiv3);
    double angle = std::atan2(std::sqrt(-q),mbDiv2)*inv3;
    double cs = std::cos(angle);
    double sn = std::sin(angle);
    *r0 = static_cast<T>(c2Div3 + 2.0*magnitude*cs);
    *r1 = static_cast<T>(c2Div3 - magnitude*(cs + root3*sn));
    *r2 = static_cast<T>(c2Div3 - magnitude*(cs - root3*sn));
    if(*r0 < *r1)
        std::swap(*r0, *r1);
    if(*r0 < *r2)
        std::swap(*r0, *r2);
    if(*r1 < *r2)
        std::swap(*r1, *r2);
}

namespace detail {

template <class T>
T ellipticRD(T x, T y, T z)
{
    double f = 1.0, s = 0.0, X, Y, Z, m;
    for(;;)
    {
        m = (x + y + 3.0*z) / 5.0;
        X = 1.0 - x/m;
        Y = 1.0 - y/m;
        Z = 1.0 - z/m;
        if(std::max(std::max(VIGRA_CSTD::fabs(X), VIGRA_CSTD::fabs(Y)), VIGRA_CSTD::fabs(Z)) < 0.01)
            break;
        double l = VIGRA_CSTD::sqrt(x*y) + VIGRA_CSTD::sqrt(x*z) + VIGRA_CSTD::sqrt(y*z);
        s += f / (VIGRA_CSTD::sqrt(z)*(z + l));
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
                      +Z*(e/6.0+Z*(-c*9.0/22.0+a*Z*3.0/26.0))) / VIGRA_CSTD::pow(m,1.5);
}

template <class T>
T ellipticRF(T x, T y, T z)
{
    double X, Y, Z, m;
    for(;;)
    {
        m = (x + y + z) / 3.0;
        X = 1.0 - x/m;
        Y = 1.0 - y/m;
        Z = 1.0 - z/m;
        if(std::max(std::max(VIGRA_CSTD::fabs(X), VIGRA_CSTD::fabs(Y)), VIGRA_CSTD::fabs(Z)) < 0.01)
            break;
        double l = VIGRA_CSTD::sqrt(x*y) + VIGRA_CSTD::sqrt(x*z) + VIGRA_CSTD::sqrt(y*z);
        x = (x + l)/4.0;
        y = (y + l)/4.0;
        z = (z + l)/4.0;
    }
    double d = X*Y - sq(Z);
    double p = X*Y*Z;
    return (1.0 - d/10.0 + p/14.0 + sq(d)/24.0 - d*p*3.0/44.0) / VIGRA_CSTD::sqrt(m);
}

} // namespace detail

    /** \brief The incomplete elliptic integral of the first kind.

        This function computes

        \f[
             \mbox{F}(x, k) = \int_0^x \frac{1}{\sqrt{1 - k^2 \sin(t)^2}} dt
        \f]

        according to the algorithm given in Press et al. "Numerical Recipes".

        Note: In some libraries (e.g. Mathematica), the second parameter of the elliptic integral
        functions must be k^2 rather than k. Check the documentation when results disagree!

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline double ellipticIntegralF(double x, double k)
{
    double c2 = sq(VIGRA_CSTD::cos(x));
    double s = VIGRA_CSTD::sin(x);
    return s*detail::ellipticRF(c2, 1.0 - sq(k*s), 1.0);
}

    /** \brief The incomplete elliptic integral of the second kind.

        This function computes

        \f[
            \mbox{E}(x, k) = \int_0^x \sqrt{1 - k^2 \sin(t)^2} dt
        \f]

        according to the algorithm given in Press et al. "Numerical Recipes". The
        complete elliptic integral of the second kind is simply <tt>ellipticIntegralE(M_PI/2, k)</TT>.

        Note: In some libraries (e.g. Mathematica), the second parameter of the elliptic integral
        functions must be k^2 rather than k. Check the documentation when results disagree!

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline double ellipticIntegralE(double x, double k)
{
    double c2 = sq(VIGRA_CSTD::cos(x));
    double s = VIGRA_CSTD::sin(x);
    k = sq(k*s);
    return s*(detail::ellipticRF(c2, 1.0-k, 1.0) - k/3.0*detail::ellipticRD(c2, 1.0-k, 1.0));
}

#if defined(_MSC_VER) && _MSC_VER < 1800

namespace detail {

template <class T>
double erfImpl(T x)
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

} // namespace detail

    /** \brief The error function.

        If <tt>erf()</tt> is not provided in the C standard math library (as it should according to the
        new C99 standard ?), VIGRA implements <tt>erf()</tt> as an approximation of the error
        function

        \f[
            \mbox{erf}(x) = \int_0^x e^{-t^2} dt
        \f]

        according to the formula given in Press et al. "Numerical Recipes".

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline double erf(double x)
{
    return detail::erfImpl(x);
}

#else

using ::erf;

#endif

namespace detail {

template <class T>
double noncentralChi2CDFApprox(unsigned int degreesOfFreedom, T noncentrality, T arg)
{
    double a = degreesOfFreedom + noncentrality;
    double b = (a + noncentrality) / sq(a);
    double t = (VIGRA_CSTD::pow((double)arg / a, 1.0/3.0) - (1.0 - 2.0 / 9.0 * b)) / VIGRA_CSTD::sqrt(2.0 / 9.0 * b);
    return 0.5*(1.0 + erf(t/VIGRA_CSTD::sqrt(2.0)));
}

template <class T>
void noncentralChi2OneIteration(T arg, T & lans, T & dans, T & pans, unsigned int & j)
{
    double tol = -50.0;
    if(lans < tol)
    {
        lans = lans + VIGRA_CSTD::log(arg / j);
        dans = VIGRA_CSTD::exp(lans);
    }
    else
    {
        dans = dans * arg / j;
    }
    pans = pans - dans;
    j += 2;
}

template <class T>
std::pair<double, double> noncentralChi2CDF(unsigned int degreesOfFreedom, T noncentrality, T arg, T eps)
{
    vigra_precondition(noncentrality >= 0.0 && arg >= 0.0 && eps > 0.0,
        "noncentralChi2P(): parameters must be positive.");
    if (arg == 0.0 && degreesOfFreedom > 0)
        return std::make_pair(0.0, 0.0);

    // Determine initial values
    double b1 = 0.5 * noncentrality,
           ao = VIGRA_CSTD::exp(-b1),
           eps2 = eps / ao,
           lnrtpi2 = 0.22579135264473,
           probability, density, lans, dans, pans, sum, am, hold;
    unsigned int maxit = 500,
        i, m;
    if(degreesOfFreedom % 2)
    {
        i = 1;
        lans = -0.5 * (arg + VIGRA_CSTD::log(arg)) - lnrtpi2;
        dans = VIGRA_CSTD::exp(lans);
        pans = erf(VIGRA_CSTD::sqrt(arg/2.0));
    }
    else
    {
        i = 2;
        lans = -0.5 * arg;
        dans = VIGRA_CSTD::exp(lans);
        pans = 1.0 - dans;
    }

    // Evaluate first term
    if(degreesOfFreedom == 0)
    {
        m = 1;
        degreesOfFreedom = 2;
        am = b1;
        sum = 1.0 / ao - 1.0 - am;
        density = am * dans;
        probability = 1.0 + am * pans;
    }
    else
    {
        m = 0;
        degreesOfFreedom = degreesOfFreedom - 1;
        am = 1.0;
        sum = 1.0 / ao - 1.0;
        while(i < degreesOfFreedom)
            detail::noncentralChi2OneIteration(arg, lans, dans, pans, i);
        degreesOfFreedom = degreesOfFreedom + 1;
        density = dans;
        probability = pans;
    }
    // Evaluate successive terms of the expansion
    for(++m; m<maxit; ++m)
    {
        am = b1 * am / m;
        detail::noncentralChi2OneIteration(arg, lans, dans, pans, degreesOfFreedom);
        sum = sum - am;
        density = density + am * dans;
        hold = am * pans;
        probability = probability + hold;
        if((pans * sum < eps2) && (hold < eps2))
            break; // converged
    }
    if(m == maxit)
        vigra_fail("noncentralChi2P(): no convergence.");
    return std::make_pair(0.5 * ao * density, std::min(1.0, std::max(0.0, ao * probability)));
}

} // namespace detail

    /** \brief Chi square distribution.

        Computes the density of a chi square distribution with \a degreesOfFreedom
        and tolerance \a accuracy at the given argument \a arg
        by calling <tt>noncentralChi2(degreesOfFreedom, 0.0, arg, accuracy)</tt>.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline double chi2(unsigned int degreesOfFreedom, double arg, double accuracy = 1e-7)
{
    return detail::noncentralChi2CDF(degreesOfFreedom, 0.0, arg, accuracy).first;
}

    /** \brief Cumulative chi square distribution.

        Computes the cumulative density of a chi square distribution with \a degreesOfFreedom
        and tolerance \a accuracy at the given argument \a arg, i.e. the probability that
        a random number drawn from the distribution is below \a arg
        by calling <tt>noncentralChi2CDF(degreesOfFreedom, 0.0, arg, accuracy)</tt>.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline double chi2CDF(unsigned int degreesOfFreedom, double arg, double accuracy = 1e-7)
{
    return detail::noncentralChi2CDF(degreesOfFreedom, 0.0, arg, accuracy).second;
}

    /** \brief Non-central chi square distribution.

        Computes the density of a non-central chi square distribution with \a degreesOfFreedom,
        noncentrality parameter \a noncentrality and tolerance \a accuracy at the given argument
        \a arg. It uses Algorithm AS 231 from Appl. Statist. (1987) Vol.36, No.3 (code ported from
        http://lib.stat.cmu.edu/apstat/231). The algorithm has linear complexity in the number of
        degrees of freedom.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline double noncentralChi2(unsigned int degreesOfFreedom,
              double noncentrality, double arg, double accuracy = 1e-7)
{
    return detail::noncentralChi2CDF(degreesOfFreedom, noncentrality, arg, accuracy).first;
}

    /** \brief Cumulative non-central chi square distribution.

        Computes the cumulative density of a chi square distribution with \a degreesOfFreedom,
        noncentrality parameter \a noncentrality and tolerance \a accuracy at the given argument
        \a arg, i.e. the probability that a random number drawn from the distribution is below \a arg
        It uses Algorithm AS 231 from Appl. Statist. (1987) Vol.36, No.3 (code ported from
        http://lib.stat.cmu.edu/apstat/231). The algorithm has linear complexity in the number of
        degrees of freedom (see noncentralChi2CDFApprox() for a constant-time algorithm).

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline double noncentralChi2CDF(unsigned int degreesOfFreedom,
              double noncentrality, double arg, double accuracy = 1e-7)
{
    return detail::noncentralChi2CDF(degreesOfFreedom, noncentrality, arg, accuracy).second;
}

    /** \brief Cumulative non-central chi square distribution (approximate).

        Computes approximate values of the cumulative density of a chi square distribution with \a degreesOfFreedom,
        and noncentrality parameter \a noncentrality at the given argument
        \a arg, i.e. the probability that a random number drawn from the distribution is below \a arg
        It uses the approximate transform into a normal distribution due to Wilson and Hilferty
        (see Abramovitz, Stegun: "Handbook of Mathematical Functions", formula 26.3.32).
        The algorithm's running time is independent of the inputs, i.e. is should be used
        when noncentralChi2CDF() is too slow, and approximate values are sufficient. The accuracy is only
        about 0.1 for few degrees of freedom, but reaches about 0.001 above dof = 5.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline double noncentralChi2CDFApprox(unsigned int degreesOfFreedom, double noncentrality, double arg)
{
    return detail::noncentralChi2CDFApprox(degreesOfFreedom, noncentrality, arg);
}

namespace detail  {

// computes (l+m)! / (l-m)!
// l and m must be positive
template <class T>
T facLM(T l, T m)
{
    T tmp = NumericTraits<T>::one();
    for(T f = l-m+1; f <= l+m; ++f)
        tmp *= f;
    return tmp;
}

} // namespace detail

    /** \brief Associated Legendre polynomial.

        Computes the value of the associated Legendre polynomial of order <tt>l, m</tt>
        for argument <tt>x</tt>. <tt>x</tt> must be in the range <tt>[-1.0, 1.0]</tt>,
        otherwise an exception is thrown. The standard Legendre polynomials are the
        special case <tt>m == 0</tt>.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class REAL>
REAL legendre(unsigned int l, int m, REAL x)
{
    vigra_precondition(abs(x) <= 1.0, "legendre(): x must be in [-1.0, 1.0].");
    if (m < 0)
    {
        m = -m;
        REAL s = odd(m)
                   ? -1.0
                   :  1.0;
        return legendre(l,m,x) * s / detail::facLM<REAL>(l,m);
    }
    REAL result = 1.0;
    if (m > 0)
    {
        REAL r = std::sqrt( (1.0-x) * (1.0+x) );
        REAL f = 1.0;
        for (int i=1; i<=m; i++)
        {
            result *= (-f) * r;
            f += 2.0;
        }
    }
    if((int)l == m)
        return result;

    REAL result_1 = x * (2.0 * m + 1.0) * result;
    if((int)l == m+1)
        return result_1;
    REAL other = 0.0;
    for(unsigned int i = m+2; i <= l; ++i)
    {
        other = ( (2.0*i-1.0) * x * result_1 - (i+m-1.0)*result) / (i-m);
        result = result_1;
        result_1 = other;
    }
    return other;
}

    /** \brief \brief Legendre polynomial.

        Computes the value of the Legendre polynomial of order <tt>l</tt> for argument <tt>x</tt>.
        <tt>x</tt> must be in the range <tt>[-1.0, 1.0]</tt>, otherwise an exception is thrown.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class REAL>
REAL legendre(unsigned int l, REAL x)
{
    return legendre(l, 0, x);
}

    /** \brief sin(pi*x).

        Essentially calls <tt>std::sin(M_PI*x)</tt> but uses a more accurate implementation
        to make sure that <tt>sin_pi(1.0) == 0.0</tt> (which does not hold for
        <tt>std::sin(M_PI)</tt> due to round-off error), and <tt>sin_pi(0.5) == 1.0</tt>.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class REAL>
REAL sin_pi(REAL x)
{
    if(x < 0.0)
        return -sin_pi(-x);
    if(x < 0.5)
        return std::sin(M_PI * x);

    bool invert = false;
    if(x < 1.0)
    {
        invert = true;
        x = -x;
    }

    REAL rem = std::floor(x);
    if(odd((int)rem))
        invert = !invert;
    rem = x - rem;
    if(rem > 0.5)
        rem = 1.0 - rem;
    if(rem == 0.5)
        rem = NumericTraits<REAL>::one();
    else
        rem = std::sin(M_PI * rem);
    return invert
              ? -rem
              : rem;
}

    /** \brief cos(pi*x).

        Essentially calls <tt>std::cos(M_PI*x)</tt> but uses a more accurate implementation
        to make sure that <tt>cos_pi(1.0) == -1.0</tt> and <tt>cos_pi(0.5) == 0.0</tt>.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class REAL>
REAL cos_pi(REAL x)
{
    return sin_pi(x+0.5);
}

namespace detail {

template <class REAL>
struct GammaImpl
{
    static REAL gamma(REAL x);
    static REAL loggamma(REAL x);

    static double g[];
    static double a[];
    static double t[];
    static double u[];
    static double v[];
    static double s[];
    static double r[];
    static double w[];
};

template <class REAL>
double GammaImpl<REAL>::g[] = {
    1.0,
    0.5772156649015329,
   -0.6558780715202538,
   -0.420026350340952e-1,
    0.1665386113822915,
   -0.421977345555443e-1,
   -0.9621971527877e-2,
    0.7218943246663e-2,
   -0.11651675918591e-2,
   -0.2152416741149e-3,
    0.1280502823882e-3,
   -0.201348547807e-4,
   -0.12504934821e-5,
    0.1133027232e-5,
   -0.2056338417e-6,
    0.6116095e-8,
    0.50020075e-8,
   -0.11812746e-8,
    0.1043427e-9,
    0.77823e-11,
   -0.36968e-11,
    0.51e-12,
   -0.206e-13,
   -0.54e-14,
    0.14e-14
};

template <class REAL>
double GammaImpl<REAL>::a[] = {
    7.72156649015328655494e-02,
    3.22467033424113591611e-01,
    6.73523010531292681824e-02,
    2.05808084325167332806e-02,
    7.38555086081402883957e-03,
    2.89051383673415629091e-03,
    1.19270763183362067845e-03,
    5.10069792153511336608e-04,
    2.20862790713908385557e-04,
    1.08011567247583939954e-04,
    2.52144565451257326939e-05,
    4.48640949618915160150e-05
};

template <class REAL>
double GammaImpl<REAL>::t[] = {
    4.83836122723810047042e-01,
    -1.47587722994593911752e-01,
    6.46249402391333854778e-02,
    -3.27885410759859649565e-02,
    1.79706750811820387126e-02,
    -1.03142241298341437450e-02,
    6.10053870246291332635e-03,
    -3.68452016781138256760e-03,
    2.25964780900612472250e-03,
    -1.40346469989232843813e-03,
    8.81081882437654011382e-04,
    -5.38595305356740546715e-04,
    3.15632070903625950361e-04,
    -3.12754168375120860518e-04,
    3.35529192635519073543e-04
};

template <class REAL>
double GammaImpl<REAL>::u[] = {
    -7.72156649015328655494e-02,
    6.32827064025093366517e-01,
    1.45492250137234768737e+00,
    9.77717527963372745603e-01,
    2.28963728064692451092e-01,
    1.33810918536787660377e-02
};

template <class REAL>
double GammaImpl<REAL>::v[] = {
    0.0,
    2.45597793713041134822e+00,
    2.12848976379893395361e+00,
    7.69285150456672783825e-01,
    1.04222645593369134254e-01,
    3.21709242282423911810e-03
};

template <class REAL>
double GammaImpl<REAL>::s[] = {
    -7.72156649015328655494e-02,
    2.14982415960608852501e-01,
    3.25778796408930981787e-01,
    1.46350472652464452805e-01,
    2.66422703033638609560e-02,
    1.84028451407337715652e-03,
    3.19475326584100867617e-05
};

template <class REAL>
double GammaImpl<REAL>::r[] = {
    0.0,
    1.39200533467621045958e+00,
    7.21935547567138069525e-01,
    1.71933865632803078993e-01,
    1.86459191715652901344e-02,
    7.77942496381893596434e-04,
    7.32668430744625636189e-06
};

template <class REAL>
double GammaImpl<REAL>::w[] = {
    4.18938533204672725052e-01,
    8.33333333333329678849e-02,
    -2.77777777728775536470e-03,
    7.93650558643019558500e-04,
    -5.95187557450339963135e-04,
    8.36339918996282139126e-04,
    -1.63092934096575273989e-03
};

template <class REAL>
REAL GammaImpl<REAL>::gamma(REAL x)
{
    int i, k, m, ix = (int)x;
    double ga = 0.0, gr = 0.0, r = 0.0, z = 0.0;

    vigra_precondition(x <= 171.0,
        "gamma(): argument cannot exceed 171.0.");

    if (x == ix)
    {
        if (ix > 0)
        {
            ga = 1.0;               // use factorial
            for (i=2; i<ix; ++i)
            {
               ga *= i;
            }
        }
        else
        {
            vigra_precondition(false,
                 "gamma(): gamma function is undefined for 0 and negative integers.");
        }
     }
     else
     {
        if (abs(x) > 1.0)
        {
            z = abs(x);
            m = (int)z;
            r = 1.0;
            for (k=1; k<=m; ++k)
            {
                r *= (z-k);
            }
            z -= m;
        }
        else
        {
            z = x;
        }
        gr = g[24];
        for (k=23; k>=0; --k)
        {
            gr = gr*z+g[k];
        }
        ga = 1.0/(gr*z);
        if (abs(x) > 1.0)
        {
            ga *= r;
            if (x < 0.0)
            {
                ga = -M_PI/(x*ga*sin_pi(x));
            }
        }
    }
    return ga;
}

/*
 * the following code is derived from lgamma_r() by Sun
 *
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 *
 */
template <class REAL>
REAL GammaImpl<REAL>::loggamma(REAL x)
{
    vigra_precondition(x > 0.0,
        "loggamma(): argument must be positive.");

    vigra_precondition(x <= 1.0e307,
        "loggamma(): argument must not exceed 1e307.");

    double res;

    if (x < 4.2351647362715017e-22)
    {
        res = -std::log(x);
    }
    else if ((x == 2.0) || (x == 1.0))
    {
        res = 0.0;
    }
    else if (x < 2.0)
    {
        const double tc  =  1.46163214496836224576e+00;
        const double tf  = -1.21486290535849611461e-01;
        const double tt  = -3.63867699703950536541e-18;
        if (x <= 0.9)
        {
            res = -std::log(x);
            if (x >= 0.7316)
            {
                double y = 1.0-x;
                double z = y*y;
                double p1 = a[0]+z*(a[2]+z*(a[4]+z*(a[6]+z*(a[8]+z*a[10]))));
                double p2 = z*(a[1]+z*(a[3]+z*(a[5]+z*(a[7]+z*(a[9]+z*a[11])))));
                double p  = y*p1+p2;
                res  += (p-0.5*y);
            }
            else if (x >= 0.23164)
            {
                double y = x-(tc-1.0);
                double z = y*y;
                double w = z*y;
                double p1 = t[0]+w*(t[3]+w*(t[6]+w*(t[9] +w*t[12])));
                double p2 = t[1]+w*(t[4]+w*(t[7]+w*(t[10]+w*t[13])));
                double p3 = t[2]+w*(t[5]+w*(t[8]+w*(t[11]+w*t[14])));
                double p  = z*p1-(tt-w*(p2+y*p3));
                res += (tf + p);
            }
            else
            {
                double y = x;
                double p1 = y*(u[0]+y*(u[1]+y*(u[2]+y*(u[3]+y*(u[4]+y*u[5])))));
                double p2 = 1.0+y*(v[1]+y*(v[2]+y*(v[3]+y*(v[4]+y*v[5]))));
                res += (-0.5*y + p1/p2);
            }
        }
        else
        {
            res = 0.0;
            if (x >= 1.7316)
            {
                double y = 2.0-x;
                double z = y*y;
                double p1 = a[0]+z*(a[2]+z*(a[4]+z*(a[6]+z*(a[8]+z*a[10]))));
                double p2 = z*(a[1]+z*(a[3]+z*(a[5]+z*(a[7]+z*(a[9]+z*a[11])))));
                double p  = y*p1+p2;
                res  += (p-0.5*y);
            }
            else if(x >= 1.23164)
            {
                double y = x-tc;
                double z = y*y;
                double w = z*y;
                double p1 = t[0]+w*(t[3]+w*(t[6]+w*(t[9] +w*t[12])));
                double p2 = t[1]+w*(t[4]+w*(t[7]+w*(t[10]+w*t[13])));
                double p3 = t[2]+w*(t[5]+w*(t[8]+w*(t[11]+w*t[14])));
                double p  = z*p1-(tt-w*(p2+y*p3));
                res += (tf + p);
            }
            else
            {
                double y = x-1.0;
                double p1 = y*(u[0]+y*(u[1]+y*(u[2]+y*(u[3]+y*(u[4]+y*u[5])))));
                double p2 = 1.0+y*(v[1]+y*(v[2]+y*(v[3]+y*(v[4]+y*v[5]))));
                res += (-0.5*y + p1/p2);
            }
        }
    }
    else if(x < 8.0)
    {
        double i = std::floor(x);
        double y = x-i;
        double p = y*(s[0]+y*(s[1]+y*(s[2]+y*(s[3]+y*(s[4]+y*(s[5]+y*s[6]))))));
        double q = 1.0+y*(r[1]+y*(r[2]+y*(r[3]+y*(r[4]+y*(r[5]+y*r[6])))));
        res = 0.5*y+p/q;
        double z = 1.0;
        while (i > 2.0)
        {
            --i;
            z *= (y+i);
        }
        res += std::log(z);
    }
    else if (x < 2.8823037615171174e+17)
    {
        double t = std::log(x);
        double z = 1.0/x;
        double y = z*z;
        double yy = w[0]+z*(w[1]+y*(w[2]+y*(w[3]+y*(w[4]+y*(w[5]+y*w[6])))));
        res = (x-0.5)*(t-1.0)+yy;
    }
    else
    {
        res =  x*(std::log(x) - 1.0);
    }

    return res;
}


} // namespace detail

    /** \brief The gamma function.

        This function implements the algorithm from<br>
        Zhang and Jin: "Computation of Special Functions", John Wiley and Sons, 1996.

        The argument must be <= 171.0 and cannot be zero or a negative integer. An
        exception is thrown when these conditions are violated.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline double gamma(double x)
{
    return detail::GammaImpl<double>::gamma(x);
}

    /** \brief The natural logarithm of the gamma function.

        This function is based on a free implementation by Sun Microsystems, Inc., see
        <a href="http://www.sourceware.org/cgi-bin/cvsweb.cgi/~checkout~/src/newlib/libm/mathfp/er_lgamma.c?rev=1.6&content-type=text/plain&cvsroot=src">sourceware.org</a> archive. It can be removed once all compilers support the new C99
        math functions.

        The argument must be positive and < 1e30. An exception is thrown when these conditions are violated.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
inline double loggamma(double x)
{
    return detail::GammaImpl<double>::loggamma(x);
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

    /** \brief Tolerance based floating-point equality.

        Check whether two floating point numbers are equal within the given tolerance.
        This is useful because floating point numbers that should be equal in theory are
        rarely exactly equal in practice. If the tolerance \a epsilon is not given,
        twice the machine epsilon is used.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T1, class T2>
bool
closeAtTolerance(T1 l, T2 r, typename PromoteTraits<T1, T2>::Promote epsilon)
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
inline bool closeAtTolerance(T1 l, T2 r)
{
    typedef typename PromoteTraits<T1, T2>::Promote T;
    return closeAtTolerance(l, r, T(2.0) * NumericTraits<T>::epsilon());
}

    /** \brief Tolerance based floating-point less-or-equal.

        Check whether two floating point numbers are less or equal within the given tolerance.
        That is, \a l can actually be greater than \a r within the given \a epsilon.
        This is useful because floating point numbers that should be equal in theory are
        rarely exactly equal in practice. If the tolerance \a epsilon is not given,
        twice the machine epsilon is used.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T1, class T2>
inline bool
lessEqualAtTolerance(T1 l, T2 r, typename PromoteTraits<T1, T2>::Promote epsilon)
{
    return l < r || closeAtTolerance(l, r, epsilon);
}

template <class T1, class T2>
inline bool lessEqualAtTolerance(T1 l, T2 r)
{
    typedef typename PromoteTraits<T1, T2>::Promote T;
    return lessEqualAtTolerance(l, r, T(2.0) * NumericTraits<T>::epsilon());
}

    /** \brief Tolerance based floating-point greater-or-equal.

        Check whether two floating point numbers are greater or equal within the given tolerance.
        That is, \a l can actually be less than \a r within the given \a epsilon.
        This is useful because floating point numbers that should be equal in theory are
        rarely exactly equal in practice. If the tolerance \a epsilon is not given,
        twice the machine epsilon is used.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T1, class T2>
inline bool
greaterEqualAtTolerance(T1 l, T2 r, typename PromoteTraits<T1, T2>::Promote epsilon)
{
    return r < l || closeAtTolerance(l, r, epsilon);
}

template <class T1, class T2>
inline bool greaterEqualAtTolerance(T1 l, T2 r)
{
    typedef typename PromoteTraits<T1, T2>::Promote T;
    return greaterEqualAtTolerance(l, r, T(2.0) * NumericTraits<T>::epsilon());
}

//@}

#define VIGRA_MATH_FUNC_HELPER(TYPE) \
    inline TYPE clipLower(const TYPE t){ \
        return t < static_cast<TYPE>(0.0) ? static_cast<TYPE>(0.0) : t; \
    } \
    inline TYPE clipLower(const TYPE t,const TYPE valLow){ \
        return t < static_cast<TYPE>(valLow) ? static_cast<TYPE>(valLow) : t; \
    } \
    inline TYPE clipUpper(const TYPE t,const TYPE valHigh){ \
        return t > static_cast<TYPE>(valHigh) ? static_cast<TYPE>(valHigh) : t; \
    } \
    inline TYPE clip(const TYPE t,const TYPE valLow, const TYPE valHigh){ \
        if(t<valLow) \
            return valLow; \
        else if(t>valHigh) \
            return valHigh; \
        else  \
            return t; \
    } \
    inline TYPE sum(const TYPE t){ \
        return t; \
    }\
    inline NumericTraits<TYPE>::RealPromote mean(const TYPE t){ \
        return t; \
    }\
    inline TYPE isZero(const TYPE t){ \
        return t==static_cast<TYPE>(0); \
    } \
    inline NumericTraits<TYPE>::RealPromote sizeDividedSquaredNorm(const TYPE t){ \
        return  squaredNorm(t); \
    } \
    inline NumericTraits<TYPE>::RealPromote sizeDividedNorm(const TYPE t){ \
        return  norm(t); \
    }


VIGRA_MATH_FUNC_HELPER(unsigned char)
VIGRA_MATH_FUNC_HELPER(unsigned short)
VIGRA_MATH_FUNC_HELPER(unsigned int)
VIGRA_MATH_FUNC_HELPER(unsigned long)
VIGRA_MATH_FUNC_HELPER(unsigned long long)
VIGRA_MATH_FUNC_HELPER(signed char)
VIGRA_MATH_FUNC_HELPER(signed short)
VIGRA_MATH_FUNC_HELPER(signed int)
VIGRA_MATH_FUNC_HELPER(signed long)
VIGRA_MATH_FUNC_HELPER(signed long long)
VIGRA_MATH_FUNC_HELPER(float)
VIGRA_MATH_FUNC_HELPER(double)
VIGRA_MATH_FUNC_HELPER(long double)



#undef VIGRA_MATH_FUNC_HELPER


} // namespace vigra

#endif /* VIGRA_MATHUTIL_HXX */
