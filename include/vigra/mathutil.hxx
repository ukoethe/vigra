/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
#include "vigra/config.hxx"

/*! \page MathConstants Mathematical Constants

    <TT>M_PI, M_SQRT2</TT>

    <b>\#include</b> "<a href="utilities_8hxx-source.html">vigra/mathutil.hxx</a>"

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

    Useful math functions that didn't make it into the C++ standard.
*/
//@{
/*! The error function.

    With the exception of Solaris (where <tt>erf()</tt> is provided as an extension of the 
    C math library), VIGRA implements <tt>erf()</tt> as an approximation of the error 
    function
    
    \f[
        \text{erf}(x) = \int_0^x e^{-x^2} dx
    \f]
    
    according to the formula given in Press et al. "Numerical Recipes".

    <b>\#include</b> "<a href="utilities_8hxx-source.html">vigra/mathutil.hxx</a>"<br>
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

template <class T>
inline
typename NumericTraits<T>::Promote
sq(T const & t)
{
    return t*t;
}

//@}

} // namespace vigra


#endif /* VIGRA_MATHUTIL_HXX */
