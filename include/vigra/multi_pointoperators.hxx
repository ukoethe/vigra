//-- -*- c++ -*-
/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Christian-Dennis Rahn                */
/*                        and Ullrich Koethe                            */
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

#ifndef VIGRA_MULTI_POINTOPERATORS_H
#define VIGRA_MULTI_POINTOPERATORS_H

#include <vigra/initimage.hxx>
#include <vigra/copyimage.hxx>
#include <vigra/transformimage.hxx>
#include <vigra/combineimages.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/metaprogramming.hxx>



namespace vigra
{

/** \addtogroup MultiPointoperators Point operators for multi-dimensional arrays.

    Copy, transform, and inspect arbitrary dimensional arrays which are represented
    by iterators compatible to \ref MultiIteratorPage. Note that are range is here
    specified by a pair: an iterator referring to the first point of the array 
    and a shape object specifying the size of the (rectangular) ROI.

    <b>\#include</b> "<a href="multi_pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"
*/
//@{

/********************************************************/
/*                                                      */
/*                    initMultiArray                    */
/*                                                      */
/********************************************************/

template <class Iterator, class Shape, class Accessor, 
          class VALUETYPE>
inline void
initMultiArrayImpl(Iterator s, Shape const & shape, Accessor a,  VALUETYPE v, MetaInt<0>)
{
    initLine(s, s + shape[0], a, v);
}
    
template <class Iterator, class Shape, class Accessor, 
          class VALUETYPE, int N>
void
initMultiArrayImpl(Iterator s, Shape const & shape, Accessor a,  
                   VALUETYPE v, MetaInt<N>)
{
    Iterator send = s + shape[N];
    for(; s != send; ++s)
    {
        initMultiArrayImpl(s.begin(), shape, a, v, MetaInt<N-1>());
    }
}
    
/** \brief Write a value to every pixel in a multi-dimensional array.

    This function can be used to init the array which must be represented by
    a pair of iterators compatible to \ref vigra::MultiIterator.
    It uses an accessor to access the data alements. Note that the iterator range 
    must be specified by a shape object, because otherwise we could not control
    the range simultaneously in all dimensions (this is a necessary consequence
    of the \ref vigra::MultiIterator design).
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class Iterator, class Shape, class Accessor, class VALUETYPE>
        void
        initMultiArray(Iterator s, Shape const & shape, Accessor a,  VALUETYPE v);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class Iterator, class Shape, class Accessor, class VALUETYPE>
        void
        initMultiArray(triple<Iterator, Shape, Accessor> const & s, VALUETYPE v);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="multi_pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"br>
    Namespace: vigra
    
    \code
    typedef vigra::MultiArray<3, int> Array;
    Array array(Array::size_type(100, 200, 50));
    
    // zero the array
    vigra::initMultiArray(destMultiArrayRange(array), 0);
    \endcode

    <b> Required Interface:</b>
    
    \code
    MultiIterator begin;
    
    Accessor accessor;
    VALUETYPE v;
    
    accessor.set(v, begin); 
    \endcode
    
*/
template <class Iterator, class Shape, class Accessor, class VALUETYPE>
inline void
initMultiArray(Iterator s, Shape const & shape, Accessor a,  VALUETYPE v)
{
    initMultiArrayImpl(s, shape, a, v, MetaInt<Iterator::level>());
}
    
template <class Iterator, class Shape, class Accessor, class VALUETYPE>
inline 
void
initMultiArray(triple<Iterator, Shape, Accessor> const & s, VALUETYPE v)
{
    initMultiArray(s.first, s.second, s.third, v);
}
    
/********************************************************/
/*                                                      */
/*                    copyMultiArray                    */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
copyMultiArrayImpl(SrcIterator s, SrcShape const & shape, SrcAccessor src,
               DestIterator d, DestAccessor dest, MetaInt<0>)
{
    copyLine(s, s + shape[0], src, d, dest);
}
    
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, int N>
void
copyMultiArrayImpl(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                   DestIterator d, DestAccessor dest, MetaInt<N>)
{
    SrcIterator send = s + shape[N];
    for(; s != send; ++s, ++d)
    {
        copyMultiArrayImpl(s.begin(), shape, src, d.begin(), dest, MetaInt<N-1>());
    }
}
    
/** \brief Copy a multi-dimensional array.

    If necessary, type conversion takes place. The arrays must be represented by
    iterators compatible with \ref vigra::MultiIterator.
    The function uses accessors to access the data elements. Note that the iterator range 
    must be specified by a shape object, because otherwise we could not control
    the range simultaneously in all dimensions (this is a necessary consequence
    of the \ref vigra::MultiIterator design).
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        copyMultiArray(SrcIterator s, 
                       SrcShape const & shape, SrcAccessor src,
                       DestIterator d, DestAccessor dest);
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        copyMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
                       pair<DestIterator, DestAccessor> const & dest);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="multi_pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"br>
    Namespace: vigra
    
    \code
    typedef vigra::MultiArray<3, int> Array;
    Array src(Array::size_type(100, 200, 50)),
          dest(Array::size_type(100, 200, 50));
    ...
    
    vigra::copyMultiArray(srcMultiArrayRange(src), destMultiArray(dest));
    \endcode

    <b> Required Interface:</b>
    
    \code
    MultiIterator src_begin, dest_begin;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    dest_accessor.set(src_accessor(src_begin), dest_begin);

    \endcode
    
*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
copyMultiArray(SrcIterator s, 
               SrcShape const & shape, SrcAccessor src,
               DestIterator d, DestAccessor dest)
{    
    copyMultiArrayImpl(s, shape, src, d, dest, MetaInt<SrcIterator::level>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
copyMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               pair<DestIterator, DestAccessor> const & dest)
{
    
    copyMultiArray(src.first, src.second, src.third, dest.first, dest.second);
}

/********************************************************/
/*                                                      */
/*                 transformMultiArray                  */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Functor>
inline void
transformMultiArrayImpl(SrcIterator s, SrcShape const & shape, SrcAccessor src,
               DestIterator d, DestAccessor dest, Functor const & f, MetaInt<0>)
{
    transformLine(s, s + shape[0], src, d, dest, f);
}
    
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Functor, int N>
void
transformMultiArrayImpl(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                   DestIterator d, DestAccessor dest, 
                   Functor const & f, MetaInt<N>)
{
    SrcIterator send = s + shape[N];
    for(; s != send; ++s, ++d)
    {
        transformMultiArrayImpl(s.begin(), shape, src, d.begin(), dest, 
                                f, MetaInt<N-1>());
    }
}
    
/** \brief Transform a multi-dimensional array with a unary function or functor.

    The transformation given by the functor is applied to every source
    element and the result written into the corresponding destination element.
    The arrays must be represented by
    iterators compatible with \ref vigra::MultiIterator.
    The function uses accessors to access the pixel data.
    Note that the unary functors of the STL can be used in addition to
    the functors specifically defined in \ref TransformFunctor.
    Creation of new functors is easiest by using \ref FunctorExpressions. Note that the iterator range 
    must be specified by a shape object, because otherwise we could not control
    the range simultaneously in all dimensions (this is a necessary consequence
    of the \ref vigra::MultiIterator design).

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, 
                  class Functor>
        void
        transformMultiArray(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                            DestIterator d, DestAccessor dest, Functor const & f);
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, 
                  class Functor>
        void
        transformMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
                       pair<DestIterator, DestAccessor> const & dest, Functor const & f);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi_pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"br>
    Namespace: vigra
    
    \code
    #include <cmath>         // for sqrt()

    typedef vigra::MultiArray<3, int> Array;
    Array src(Array::size_type(100, 200, 50)),
          dest(Array::size_type(100, 200, 50));
    ...
    
    vigra::transformMultiArray(srcMultiArrayRange(src),
                               destMultiArray(dest),
                               &std::sqrt );

    \endcode

    <b> Required Interface:</b>

    \code
    MultiIterator src_begin, src_end, dest_begin;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    Functor functor;

    dest_accessor.set(functor(src_accessor(src_begin)), dest_begin);

    \endcode

*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          class Functor>
inline void
transformMultiArray(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                    DestIterator d, DestAccessor dest, Functor const & f)
{    
    transformMultiArrayImpl(s, shape, src, d, dest, f, MetaInt<SrcIterator::level>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          class Functor>
inline void
transformMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               pair<DestIterator, DestAccessor> const & dest, Functor const & f)
{
    
    transformMultiArray(src.first, src.second, src.third, 
                        dest.first, dest.second, f);
}

/********************************************************/
/*                                                      */
/*                combineTwoMultiArrays                 */
/*                                                      */
/********************************************************/

template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestAccessor, 
          class Functor>
inline void
combineTwoMultiArraysImpl(SrcIterator1 s1, SrcShape const & shape, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               DestIterator d, DestAccessor dest, Functor const & f, MetaInt<0>)
{
    combineTwoLines(s1, s1 + shape[0], src1, s2, src2, d, dest, f);
}
    
template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestAccessor, 
          class Functor, int N>
void
combineTwoMultiArraysImpl(SrcIterator1 s1, SrcShape const & shape, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               DestIterator d, DestAccessor dest, 
               Functor const & f, MetaInt<N>)
{
    SrcIterator1 s1end = s1 + shape[N];
    for(; s1 != s1end; ++s1, ++s2, ++d)
    {
        combineTwoMultiArraysImpl(s1.begin(), shape, src1, 
                                  s2.begin(), src2, d.begin(), dest, 
                                  f, MetaInt<N-1>());
    }
}
    
/** \brief Combine two multi-dimensional arrays into one using a binary function or functor.

    The transformation given by the functor is applied to the source 
    array elements and the result written into the corresponding destination element.
    This is typically used for operations like add and subtract.
    The arrays must be represented by
    iterators compatible with \ref vigra::MultiIterator.
    The function uses accessors to access the pixel data.
    Note that the binary functors of the STL can be used in addition to
    the functors specifically defined in \ref CombineFunctor.
    Creation of new functors is easiest by using \ref FunctorExpressions. Note that the iterator range 
    must be specified by a shape object, because otherwise we could not control
    the range simultaneously in all dimensions (this is a necessary consequence
    of the \ref vigra::MultiIterator design).
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator1, class SrcShape, class SrcAccessor1,
                  class SrcIterator2, class SrcAccessor2,
                  class DestIterator, class DestAccessor, 
                  class Functor>
        inline void
        combineTwoMultiArrays(SrcIterator1 s1, SrcShape const & shape, SrcAccessor1 src1,
                       SrcIterator2 s2, SrcAccessor2 src2,
                       DestIterator d, DestAccessor dest, Functor const & f);
            }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator1, class SrcShape, class SrcAccessor1,
                  class SrcIterator2, class SrcAccessor2,
                  class DestIterator, class DestAccessor, class Functor>
        void
        combineTwoMultiArrays(triple<SrcIterator1, SrcShape, SrcAccessor1> const & src1,
                       pair<SrcIterator2, SrcAccessor2> const & src2,
                       pair<DestIterator, DestAccessor> const & dest, Functor const & f);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="multi_pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"br>
    Namespace: vigra
    
    \code
    #include <functional>     // for plus

    typedef vigra::MultiArray<3, int> Array;
    Array src1(Array::size_type(100, 200, 50)),
          src2(Array::size_type(100, 200, 50)),
          dest(Array::size_type(100, 200, 50));
    ...
    
    vigra::combineTwoMultiArrays(
                srcMultiArrayRange(src1), 
                srcMultiArray(src2), 
                destMultiArray(dest),  
                std::plus<SrcValueType>());
    
    \endcode
    
    Note that <TT>SrcValueType</TT> must be replaced with the appropriate type (e.g. 
    the promote type of the input images' pixel type, see also 
    \ref NumericPromotionTraits)
    
    <b> Required Interface:</b>
    
    \code
    MultiIterator src1_begin, src2_begin, dest_begin;
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;

    dest_accessor.set(
          functor(src1_accessor(src1_begin), src2_accessor(src2_begin)), 
          dest_begin);

    \endcode
    
    
*/
template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestAccessor, 
          class Functor>
inline void
combineTwoMultiArrays(SrcIterator1 s1, SrcShape const & shape, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               DestIterator d, DestAccessor dest, Functor const & f)
{    
    combineTwoMultiArraysImpl(s1, shape, src1, s2, src2, d, dest, f, 
                              MetaInt<SrcIterator1::level>());
}

template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestAccessor, class Functor>
inline void
combineTwoMultiArrays(triple<SrcIterator1, SrcShape, SrcAccessor1> const & src1,
               pair<SrcIterator2, SrcAccessor2> const & src2,
               pair<DestIterator, DestAccessor> const & dest, Functor const & f)
{
    
    combineTwoMultiArrays(
           src1.first, src1.second, src1.third, 
           src2.first, src2.second, dest.first, dest.second, f);
}

/********************************************************/
/*                                                      */
/*               combineThreeMultiArrays                */
/*                                                      */
/********************************************************/

template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class SrcIterator3, class SrcAccessor3,
          class DestIterator, class DestAccessor, 
          class Functor>
inline void
combineThreeMultiArraysImpl(SrcIterator1 s1, SrcShape const & shape, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               SrcIterator3 s3, SrcAccessor3 src3,
               DestIterator d, DestAccessor dest, Functor const & f, MetaInt<0>)
{
    combineThreeLines(s1, s1 + shape[0], src1, s2, src2, s3, src3, d, dest, f);
}
    
template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class SrcIterator3, class SrcAccessor3,
          class DestIterator, class DestAccessor, 
          class Functor, int N>
void
combineThreeMultiArraysImpl(SrcIterator1 s1, SrcShape const & shape, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               SrcIterator3 s3, SrcAccessor3 src3,
               DestIterator d, DestAccessor dest, 
                   Functor const & f, MetaInt<N>)
{
    SrcIterator1 s1end = s1 + shape[N];
    for(; s1 != s1end; ++s1, ++s2, ++s3, ++d)
    {
        combineThreeMultiArraysImpl(s1.begin(), shape, src1, 
                                  s2.begin(), src2, s3.begin(), src3, d.begin(), dest, 
                                  f, MetaInt<N-1>());
    }
}
    
    
/** \brief Combine three multi-dimensional arrays into one using a binary function or functor.

    Except for the fact that it operates on three input arrays, this function is
    identical to \ref combineTwoMultiArrays().
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator1, class SrcShape, class SrcAccessor1,
                  class SrcIterator2, class SrcAccessor2,
                  class SrcIterator3, class SrcAccessor3,
                  class DestIterator, class DestAccessor, 
                  class Functor>
        void
        combineThreeMultiArrays(SrcIterator1 s1, SrcShape const & shape, SrcAccessor1 src1,
                       SrcIterator2 s2, SrcAccessor2 src2,
                       SrcIterator3 s3, SrcAccessor3 src3,
                       DestIterator d, DestAccessor dest, Functor const & f);
                    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator1, class SrcShape, class SrcAccessor1,
                  class SrcIterator2, class SrcAccessor2,
                  class SrcIterator3, class SrcAccessor3,
                  class DestIterator, class DestAccessor, 
                  class Functor>
        inline void
        combineThreeMultiArrays(triple<SrcIterator1, SrcShape, SrcAccessor1> const & src1,
                       pair<SrcIterator2, SrcAccessor2> const & src2,
                       pair<SrcIterator3, SrcAccessor3> const & src3,
                       pair<DestIterator, DestAccessor> const & dest, Functor const & f);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="multi_pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"br>
    Namespace: vigra
    
    \code
    #include <functional>     // for plus

    typedef vigra::MultiArray<3, int> Array;
    Array src1(Array::size_type(100, 200, 50)),
          src2(Array::size_type(100, 200, 50)),
          src3(Array::size_type(100, 200, 50)),
          dest(Array::size_type(100, 200, 50));
    ...
    
    vigra::combineThreeMultiArrays(
                srcMultiArrayRange(src1), 
                srcMultiArray(src2), 
                srcMultiArray(src3), 
                destMultiArray(dest),  
                SomeThreeArgumentFunctor());
    
    \endcode
*/
template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class SrcIterator3, class SrcAccessor3,
          class DestIterator, class DestAccessor, 
          class Functor>
inline void
combineThreeMultiArrays(SrcIterator1 s1, SrcShape const & shape, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               SrcIterator3 s3, SrcAccessor3 src3,
               DestIterator d, DestAccessor dest, Functor const & f)
{    
    combineThreeMultiArraysImpl(s1, shape, src1, s2, src2, s3, src3, d, dest, f, 
                              MetaInt<SrcIterator1::level>());
}

template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class SrcIterator3, class SrcAccessor3,
          class DestIterator, class DestAccessor, 
          class Functor>
inline void
combineThreeMultiArrays(triple<SrcIterator1, SrcShape, SrcAccessor1> const & src1,
               pair<SrcIterator2, SrcAccessor2> const & src2,
               pair<SrcIterator3, SrcAccessor3> const & src3,
               pair<DestIterator, DestAccessor> const & dest, Functor const & f)
{
    
    combineThreeMultiArrays(
           src1.first, src1.second, src1.third, 
           src2.first, src2.second, src3.first, src3.second, dest.first, dest.second, f);
}

/********************************************************/
/*                                                      */
/*                  inspectMultiArray                   */
/*                                                      */
/********************************************************/

template <class Iterator, class Shape, class Accessor, class Functor>
inline void
inspectMultiArrayImpl(Iterator s, Shape const & shape, Accessor a,  Functor & f, MetaInt<0>)
{
    inspectLine(s, s + shape[0], a, f);
}
    
template <class Iterator, class Shape, class Accessor, class Functor, int N>
void
inspectMultiArrayImpl(Iterator s, Shape const & shape, Accessor a,  Functor & f, MetaInt<N>)
{
    Iterator send = s + shape[N];
    for(; s != send; ++s)
    {
        inspectMultiArrayImpl(s.begin(), shape, a, f, MetaInt<N-1>());
    }
}
    
/** \brief Apply read-only functor to every element of a multi-dimensional array.

    This function can be used to collect statistics of the array etc.
    The results must be stored in the functor, which serves as a return
    value. The arrays must be represented by
    iterators compatible with \ref vigra::MultiIterator.
    The function uses an accessor to access the pixel data. Note that the iterator range 
    must be specified by a shape object, because otherwise we could not control
    the range simultaneously in all dimensions (this is a necessary consequence
    of the \ref vigra::MultiIterator design).

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class Iterator, class Shape, class Accessor, class Functor>
        void
        inspectMultiArray(Iterator s, Shape const & shape, Accessor a,  Functor & f);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class Iterator, class Shape, class Accessor, class Functor>
        void
        inspectMultiArray(triple<Iterator, Shape, Accessor> const & s, Functor & f);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi_pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"br>
    Namespace: vigra

    \code
    typedef vigra::MultiArray<3, int> Array;
    Array array(Array::size_type(100, 200, 50));

    // init functor
    vigra::FindMinMax<int> minmax;

    vigra::inspectMultiArray(srcMultiArrayRange(array), minmax);

    cout << "Min: " << minmax.min << " Max: " << minmax.max;

    \endcode

    <b> Required Interface:</b>

    \code
    MultiIterator src_begin;

    Accessor accessor;
    Functor functor;

    functor(accessor(src_begin)); 
    \endcode

*/
template <class Iterator, class Shape, class Accessor, class Functor>
inline void
inspectMultiArray(Iterator s, Shape const & shape, Accessor a,  Functor & f)
{
    inspectMultiArrayImpl(s, shape, a, f, MetaInt<Iterator::level>());
}
    
template <class Iterator, class Shape, class Accessor, class Functor>
inline void
inspectMultiArray(triple<Iterator, Shape, Accessor> const & s, Functor & f)
{
    inspectMultiArray(s.first, s.second, s.third, f);
}
    
/********************************************************/
/*                                                      */
/*                  inspectTwoMultiArrays               */
/*                                                      */
/********************************************************/

template <class Iterator1, class Shape, class Accessor1, 
          class Iterator2, class Accessor2, 
          class Functor>
inline void
inspectTwoMultiArraysImpl(Iterator1 s1, Shape const & shape, Accessor1 a1,
                          Iterator2 s2, Accessor2 a2,
                          Functor & f, MetaInt<0>)
{
    inspectTwoLines(s1, s1 + shape[0], a1, s2, a2, f);
}
    
template <class Iterator1, class Shape, class Accessor1, 
          class Iterator2, class Accessor2, 
          class Functor, int N>
void
inspectTwoMultiArraysImpl(Iterator1 s1, Shape const & shape, Accessor1 a1,
                          Iterator2 s2, Accessor2 a2,
                          Functor & f, MetaInt<N>)
{
    Iterator1 s1end = s1 + shape[N];
    for(; s1 != s1end; ++s1, ++s2)
    {
        inspectTwoMultiArraysImpl(s1.begin(), shape, a1, 
                                  s2.begin(), a2, f, MetaInt<N-1>());
    }
}
    
/** \brief Apply read-only functor to every element of a multi-dimensional array.

    This function can be used to collect statistics of the array etc.
    The results must be stored in the functor, which serves as a return
    value. The arrays must be represented by
    iterators compatible with \ref vigra::MultiIterator.
    The function uses an accessor to access the pixel data. Note that the iterator range 
    must be specified by a shape object, because otherwise we could not control
    the range simultaneously in all dimensions (this is a necessary consequence
    of the \ref vigra::MultiIterator design).

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class Iterator1, class Shape, class Accessor1, 
                  class Iterator2, class Accessor2, 
                  class Functor>
        void
        inspectTwoMultiArrays(Iterator1 s1, Shape const & shape, Accessor1 a1,
                              Iterator2 s2, Accessor2 a2, Functor & f);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class Iterator1, class Shape1, class Accessor1, 
                  class Iterator2, class Accessor2, 
                  class Functor>
        void
        inspectTwoMultiArrays(triple<Iterator1, Shape1, Accessor1> const & s1, 
                              pair<Iterator2, Accessor2> const & s2, Functor & f);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi_pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"br>
    Namespace: vigra

    \code
    typedef vigra::MultiArray<3, int> Array;
    Array array1(Array::size_type(100, 200, 50)),
          array2(Array::size_type(100, 200, 50));

    // init functor
    SomeStatisticsFunctor stats(..);

    vigra::inspectTwoMultiArrays(srcMultiArrayRange(array1), srcMultiArray(array2), stats);

    \endcode

    <b> Required Interface:</b>

    \code
    MultiIterator src1_begin, src2_begin;

    Accessor a1, a2;
    Functor functor;

    functor(a1(src1_begin), a2(src2_begin)); 
    \endcode

*/
template <class Iterator1, class Shape, class Accessor1, 
          class Iterator2, class Accessor2, 
          class Functor>
inline void
inspectTwoMultiArrays(Iterator1 s1, Shape const & shape, Accessor1 a1,
                      Iterator2 s2, Accessor2 a2, Functor & f)
{
    inspectTwoMultiArraysImpl(s1, shape, a1, s2, a2, f, MetaInt<Iterator1::level>());
}
    
template <class Iterator1, class Shape, class Accessor1, 
          class Iterator2, class Accessor2, 
          class Functor>
inline 
void
inspectTwoMultiArrays(triple<Iterator1, Shape, Accessor1> const & s1, 
                      pair<Iterator2, Accessor2> const & s2, Functor & f)
{
    inspectTwoMultiArrays(s1.first, s1.second, s1.third, 
                          s2.first, s2.second, f);
}
    
//@}

};	//-- namespace vigra


#endif	//-- VIGRA_MULTI_POINTOPERATORS_H
