/************************************************************************/
/*                                                                      */
/*        Copyright 2009-2010 by Ullrich Koethe and Janis Fehr          */
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

#ifndef VIGRA_INVARIANT_FEATURES3D_HXX
#define VIGRA_INVARIANT_FEATURES3D_HXX

#include <complex>
#include "config.hxx"
#include "error.hxx"
#include "utilities.hxx"
#include "mathutil.hxx"
#include "array_vector.hxx"
#include "matrix.hxx"
#include "tinyvector.hxx"
#include "quaternion.hxx"

namespace vigra {

namespace detail  {

// computes the normalization for SH base functions
inline double realSH(double l, double m)
{
    return std::sqrt((2.0*l + 1.0) / (4.0*M_PI*facLM(l,m)));

}

template<int N, class T, class C>
TinyVector<float, N> centerOfBB(MultiArrayView<N, T, C> const & A)
{
    return TinyVector<float, N>(A.shape()) /= 2.0;                        
}

} // namespace detail


} // namespace vigra 

#endif // VIGRA_INVARIANT_FEATURES3D_HXX
