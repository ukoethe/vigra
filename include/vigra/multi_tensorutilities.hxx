/************************************************************************/
/*                                                                      */
/*               Copyright 2009-2010 by Ullrich Koethe                  */
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

#ifndef VIGRA_MULTI_TENSORUTILITIES_HXX
#define VIGRA_MULTI_TENSORUTILITIES_HXX

#include <cmath>
#include "utilities.hxx"
#include "mathutil.hxx"
#include "metaprogramming.hxx"
#include "multi_shape.hxx"
#include "multi_pointoperators.hxx"

namespace vigra {

namespace detail {

template <int N, class ArgumentVector, class ResultVector>
class OuterProductFunctor
{
public:
    typedef ArgumentVector argument_type;
    typedef ResultVector result_type;
    typedef typename ArgumentVector::value_type ValueType;
    
    result_type operator()(argument_type const & in) const
    {
        result_type res;
        for(int b=0, i=0; i<N; ++i)
        {
            for(int j=i; j<N; ++j, ++b)
            {
                res[b] = detail::RequiresExplicitCast<ValueType>::cast(in[i]*in[j]);
            }
        }
        return res;
    }
};

template <int N, class ArgumentVector>
class TensorTraceFunctor
{
public:

    typedef ArgumentVector argument_type;
    typedef typename ArgumentVector::value_type result_type;
    
    result_type exec(argument_type const & v, MetaInt<1>) const
    {
        return v[0];
    }
    
    result_type exec(argument_type const & v, MetaInt<2>) const
    {
        return v[0] + v[2];
    }
    
    result_type exec(argument_type const & v, MetaInt<3>) const
    {
        return v[0] + v[3] + v[5];
    }
    
    template <int N2>
    void exec(argument_type const & v, result_type & r, MetaInt<N2>) const
    {
        vigra_fail("tensorTraceMultiArray(): Sorry, can only handle dimensions up to 3.");
    }

    result_type operator()( const argument_type & a ) const
    {
        return exec(a, MetaInt<N>());
    }
};

template <int N, class ArgumentVector, class ResultVector>
class EigenvaluesFunctor
{
public:

    typedef ArgumentVector argument_type;
    typedef ResultVector result_type;
    
    void exec(argument_type const & v, result_type & r, MetaInt<1>) const
    {
        symmetric2x2Eigenvalues(v[0], &r[0]);
    }
    
    void exec(argument_type const & v, result_type & r, MetaInt<2>) const
    {
        symmetric2x2Eigenvalues(v[0], v[1], v[2], &r[0], &r[1]);
    }
    
    void exec(argument_type const & v, result_type & r, MetaInt<3>) const
    {
        symmetric3x3Eigenvalues(v[0], v[1], v[2], v[3], v[4], v[5], &r[0], &r[1], &r[2]);
    }
    
    template <int N2>
    void exec(argument_type const & v, result_type & r, MetaInt<N2>) const
    {
        vigra_fail("tensorEigenvaluesMultiArray(): Sorry, can only handle dimensions up to 3.");
    }

    result_type operator()( const argument_type & a ) const
    {
        result_type res;
        exec(a, res, MetaInt<N>());
        return res;
    }
};


template <int N, class ArgumentVector>
class DeterminantFunctor
{
public:

    typedef ArgumentVector argument_type;
    typedef typename ArgumentVector::value_type result_type;
    
    result_type exec(argument_type const & v, MetaInt<1>) const
    {
        return v[0];
    }
    
    result_type exec(argument_type const & v, MetaInt<2>) const
    {
        return v[0]*v[2] - sq(v[1]);
    }
    
    result_type exec(argument_type const & v, MetaInt<3>) const
    {
        result_type r0, r1, r2;
        symmetric3x3Eigenvalues(v[0], v[1], v[2], v[3], v[4], v[5], &r0, &r1, &r2);
        return r0*r1*r2;
    }
    
    template <int N2>
    void exec(argument_type const & v, result_type & r, MetaInt<N2>) const
    {
        vigra_fail("tensorDeterminantMultiArray(): Sorry, can only handle dimensions up to 3.");
    }

    result_type operator()( const argument_type & a ) const
    {
        return exec(a, MetaInt<N>());
    }
};

} // namespace detail


/** \addtogroup MultiPointoperators
*/
//@{

/********************************************************/
/*                                                      */
/*                vectorToTensorMultiArray              */
/*                                                      */
/********************************************************/

/** \brief Calculate the tensor (outer) product of a N-D vector with itself.

    This function is useful to transform vector arrays into a tensor representation 
    that can be used as input to tensor based processing and analysis functions
    (e.g. tensor smoothing). When the input array has N dimensions, the input value_type 
    must be a vector of length N, whereas the output value_type mus be vectors of length 
    N*(N-1)/2 which will represent the upper triangular part of the resulting (symmetric) 
    tensor. That is, for 2D arrays the output contains the elements 
    <tt>[t11, t12 == t21, t22]</tt> in this order, whereas it contains the elements
    <tt>[t11, t12, t13, t22, t23, t33]</tt> for 3D arrays.
    
    Currently, <tt>N <= 3</tt> is required.
    
    <b> Declarations:</b>

    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T1, class S1,
                                  class T2, class S2>
        void 
        vectorToTensorMultiArray(MultiArrayView<N, T1, S1> const & source,
                                 MultiArrayView<N, T2, S2> dest);
    }
    \endcode

    \deprecatedAPI{vectorToTensorMultiArray}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        vectorToTensorMultiArray(SrcIterator  si, SrcShape const & shape, SrcAccessor src,
                                 DestIterator di, DestAccessor dest);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        vectorToTensorMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> s,
                                 pair<DestIterator, DestAccessor> d);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_tensorutilities.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<3, float>                  vol(shape);
    MultiArray<3, TinyVector<float, 3> >  gradient(shape);
    MultiArray<3, TinyVector<float, 6> >  tensor(shape);
    
    gaussianGradientMultiArray(vol, gradient, 2.0);
    vectorToTensorMultiArray(gradient, tensor);
    \endcode
*/
doxygen_overloaded_function(template <...> void vectorToTensorMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
vectorToTensorMultiArray(SrcIterator  si, SrcShape const & shape, SrcAccessor src,
                         DestIterator di, DestAccessor dest)
{
    static const int N = SrcShape::static_size;
    static const int M = N*(N+1)/2;
    
    typedef typename SrcAccessor::value_type  SrcType;
    typedef typename DestAccessor::value_type DestType;

    for(int k=0; k<N; ++k)
        if(shape[k] <=0)
            return;

    vigra_precondition(N == (int)src.size(si),
        "vectorToTensorMultiArray(): Wrong number of channels in input array.");
    vigra_precondition(M == (int)dest.size(di),
        "vectorToTensorMultiArray(): Wrong number of channels in output array.");

    transformMultiArray(si, shape, src, di, dest, 
                        detail::OuterProductFunctor<N, SrcType, DestType>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
vectorToTensorMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> s,
                         pair<DestIterator, DestAccessor> d)
{
    vectorToTensorMultiArray(s.first, s.second, s.third, d.first, d.second);
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2>
inline void 
vectorToTensorMultiArray(MultiArrayView<N, T1, S1> const & source,
                         MultiArrayView<N, T2, S2> dest)
{
    vigra_precondition(source.shape() == dest.shape(),
        "vectorToTensorMultiArray(): shape mismatch between input and output.");
    vectorToTensorMultiArray(srcMultiArrayRange(source), destMultiArray(dest));
}

/********************************************************/
/*                                                      */
/*                tensorTraceMultiArray                 */
/*                                                      */
/********************************************************/

/** \brief Calculate the tensor trace for every element of a N-D tensor array.

    This function turns a N-D tensor (whose value_type is a vector of length N*(N+1)/2, 
    see \ref vectorToTensorMultiArray()) representing the upper triangular part of a 
    symmetric tensor into a scalar array holding the tensor trace.
    
    Currently, <tt>N <= 3</tt> is required.
    
    <b> Declarations:</b>

    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T1, class S1,
                                  class T2, class S2>
        void 
        tensorTraceMultiArray(MultiArrayView<N, T1, S1> const & source,
                              MultiArrayView<N, T2, S2> dest);
    }
    \endcode

    \deprecatedAPI{tensorTraceMultiArray}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        tensorTraceMultiArray(SrcIterator si,  SrcShape const & shape, SrcAccessor src,
                              DestIterator di, DestAccessor dest);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        tensorTraceMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> s,
                              pair<DestIterator, DestAccessor> d);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_tensorutilities.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<3, float>                  vol(shape);
    MultiArray<3, TinyVector<float, 6> >  hessian(shape);
    MultiArray<3, float>                  trace(shape);
    
    hessianOfGaussianMultiArray(vol, hessian, 2.0);
    tensorTraceMultiArray(hessian, trace);
    \endcode

    <b> Preconditions:</b>

    <tt>N == 2</tt> or <tt>N == 3</tt>
*/
doxygen_overloaded_function(template <...> void tensorTraceMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
tensorTraceMultiArray(SrcIterator si,  SrcShape const & shape, SrcAccessor src,
                      DestIterator di, DestAccessor dest)
{
    static const int N = SrcShape::static_size;
    typedef typename SrcAccessor::value_type  SrcType;

    transformMultiArray(si, shape, src, di, dest, 
                        detail::TensorTraceFunctor<N, SrcType>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
tensorTraceMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> s,
                      pair<DestIterator, DestAccessor> d)
{
    tensorTraceMultiArray(s.first, s.second, s.third, d.first, d.second);
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2>
inline void 
tensorTraceMultiArray(MultiArrayView<N, T1, S1> const & source,
                      MultiArrayView<N, T2, S2> dest)
{
    vigra_precondition(source.shape() == dest.shape(),
        "tensorTraceMultiArray(): shape mismatch between input and output.");
    tensorTraceMultiArray(srcMultiArrayRange(source), destMultiArray(dest));
}


/********************************************************/
/*                                                      */
/*             tensorEigenvaluesMultiArray              */
/*                                                      */
/********************************************************/

/** \brief Calculate the tensor eigenvalues for every element of a N-D tensor array.

    This function turns a N-D tensor (whose value_type is a vector of length N*(N+1)/2, 
    see \ref vectorToTensorMultiArray()) representing the upper triangular part of a 
    symmetric tensor into a vector-valued array holding the tensor eigenvalues (thus,
    the destination value_type must be vectors of length N).
    
    Currently, <tt>N <= 3</tt> is required.
    
    <b> Declarations:</b>

    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T1, class S1,
                                  class T2, class S2>
        void 
        tensorEigenvaluesMultiArray(MultiArrayView<N, T1, S1> const & source,
                                    MultiArrayView<N, T2, S2> dest);
    }
    \endcode

    \deprecatedAPI{tensorEigenvaluesMultiArray}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        tensorEigenvaluesMultiArray(SrcIterator si,  SrcShape const & shape, SrcAccessor src,
                                    DestIterator di, DestAccessor dest);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        tensorEigenvaluesMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> s,
                                    pair<DestIterator, DestAccessor> d);
    }
    \endcode
    \deprecatedEnd

    <b> Usage (MultiArrayView API):</b>

    <b>\#include</b> \<vigra/multi_tensorutilities.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<3, float>                  vol(shape);
    MultiArray<3, TinyVector<float, 6> >  hessian(shape);
    MultiArray<3, TinyVector<float, 3> >  eigenvalues(shape);
    
    hessianOfGaussianMultiArray(vol, hessian, 2.0);
    tensorEigenvaluesMultiArray(hessian, eigenvalues);
    \endcode

    <b> Preconditions:</b>

    <tt>N == 2</tt> or <tt>N == 3</tt>
*/
doxygen_overloaded_function(template <...> void tensorEigenvaluesMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
tensorEigenvaluesMultiArray(SrcIterator si,  SrcShape const & shape, SrcAccessor src,
                            DestIterator di, DestAccessor dest)
{
    static const int N = SrcShape::static_size;
    static const int M = N*(N+1)/2;
    
    typedef typename SrcAccessor::value_type  SrcType;
    typedef typename DestAccessor::value_type DestType;

    for(int k=0; k<N; ++k)
        if(shape[k] <=0)
            return;

    vigra_precondition(M == (int)src.size(si),
        "tensorEigenvaluesMultiArray(): Wrong number of channels in input array.");
    vigra_precondition(N == (int)dest.size(di),
        "tensorEigenvaluesMultiArray(): Wrong number of channels in output array.");

    transformMultiArray(si, shape, src, di, dest, 
                        detail::EigenvaluesFunctor<N, SrcType, DestType>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
tensorEigenvaluesMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> s,
                            pair<DestIterator, DestAccessor> d)
{
    tensorEigenvaluesMultiArray(s.first, s.second, s.third, d.first, d.second);
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2>
inline void 
tensorEigenvaluesMultiArray(MultiArrayView<N, T1, S1> const & source,
                            MultiArrayView<N, T2, S2> dest)
{
    vigra_precondition(source.shape() == dest.shape(),
        "tensorEigenvaluesMultiArray(): shape mismatch between input and output.");
    tensorEigenvaluesMultiArray(srcMultiArrayRange(source), destMultiArray(dest));
}

/********************************************************/
/*                                                      */
/*             tensorDeterminantMultiArray              */
/*                                                      */
/********************************************************/

/** \brief Calculate the tensor determinant for every element of a ND tensor array.

    This function turns a N-D tensor (whose value_type is a vector of length N*(N+1)/2, 
    see \ref vectorToTensorMultiArray()) representing the upper triangular part of a 
    symmetric tensor into the a scalar array holding the tensor determinant.
    
    Currently, <tt>N <= 3</tt> is required.
    
    <b> Declarations:</b>

    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T1, class S1,
                                  class T2, class S2>
        void 
        tensorDeterminantMultiArray(MultiArrayView<N, T1, S1> const & source,
                                    MultiArrayView<N, T2, S2> dest);
    }
    \endcode

    \deprecatedAPI{tensorDeterminantMultiArray}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        tensorDeterminantMultiArray(SrcIterator si,  SrcShape const & shape, SrcAccessor src,
                                    DestIterator di, DestAccessor dest);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        tensorDeterminantMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> s,
                                    pair<DestIterator, DestAccessor> d);
    }
    \endcode
    \deprecatedEnd

    <b> Usage (MultiArrayView API):</b>

    <b>\#include</b> \<vigra/multi_tensorutilities.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<3, float>                  vol(shape);
    MultiArray<3, TinyVector<float, 6> >  hessian(shape);
    MultiArray<3, float>                  determinant(shape);
    
    hessianOfGaussianMultiArray(vol, hessian, 2.0);
    tensorDeterminantMultiArray(hessian, determinant);
    \endcode

    <b> Preconditions:</b>

    <tt>N == 2</tt> or <tt>N == 3</tt>
*/
doxygen_overloaded_function(template <...> void tensorDeterminantMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
tensorDeterminantMultiArray(SrcIterator si,  SrcShape const & shape, SrcAccessor src,
                            DestIterator di, DestAccessor dest)
{
    typedef typename SrcAccessor::value_type  SrcType;

    static const int N = SrcShape::static_size;
    static const int M = N*(N+1)/2;
    
    for(int k=0; k<N; ++k)
        if(shape[k] <=0)
            return;

    vigra_precondition(M == (int)src.size(si),
        "tensorDeterminantMultiArray(): Wrong number of channels in output array.");

    transformMultiArray(si, shape, src, di, dest, 
                        detail::DeterminantFunctor<N, SrcType>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
tensorDeterminantMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> s,
                            pair<DestIterator, DestAccessor> d)
{
    tensorDeterminantMultiArray(s.first, s.second, s.third, d.first, d.second);
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2>
inline void 
tensorDeterminantMultiArray(MultiArrayView<N, T1, S1> const & source,
                            MultiArrayView<N, T2, S2> dest)
{
    vigra_precondition(source.shape() == dest.shape(),
        "tensorDeterminantMultiArray(): shape mismatch between input and output.");
    tensorDeterminantMultiArray(srcMultiArrayRange(source), destMultiArray(dest));
}

//@}

} // namespace vigra

#endif /* VIGRA_MULTI_TENSORUTILITIES_HXX */
