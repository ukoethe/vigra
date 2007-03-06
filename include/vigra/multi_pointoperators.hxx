//-- -*- c++ -*-
/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Ullrich Koethe                       */
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

#ifndef VIGRA_MULTI_POINTOPERATORS_H
#define VIGRA_MULTI_POINTOPERATORS_H

#include "initimage.hxx"
#include "copyimage.hxx"
#include "transformimage.hxx"
#include "combineimages.hxx"
#include "inspectimage.hxx"
#include "multi_array.hxx"
#include "metaprogramming.hxx"



namespace vigra
{

/** \addtogroup MultiPointoperators Point operators for multi-dimensional arrays.

    Copy, transform, and inspect arbitrary dimensional arrays which are represented
    by iterators compatible to \ref MultiIteratorPage. Note that are range is here
    specified by a pair: an iterator referring to the first point of the array 
    and a shape object specifying the size of the (rectangular) ROI.

    <b>\#include</b> "<a href="multi__pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"
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


        template <class Iterator, class Shape, class Accessor, class FUNCTOR>
        void
        initMultiArray(Iterator s, Shape const & shape, Accessor a,  FUNCTOR const & f);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class Iterator, class Shape, class Accessor, class VALUETYPE>
        void
        initMultiArray(triple<Iterator, Shape, Accessor> const & s, VALUETYPE v);


        template <class Iterator, class Shape, class Accessor, class FUNCTOR>
        void
        initMultiArray(triple<Iterator, Shape, Accessor> const & s, FUNCTOR const & f);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="multi__pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"<br>
    Namespace: vigra
    
    \code
    typedef vigra::MultiArray<3, int> Array;
    Array array(Array::size_type(100, 200, 50));
    
    // zero the array
    vigra::initMultiArray(destMultiArrayRange(array), 0);
    \endcode

    <b> Required Interface:</b>
    
    The function accepts either a value that is copied into every destination element: 
    
    \code
    MultiIterator begin;
    
    Accessor accessor;
    VALUETYPE v;
    
    accessor.set(v, begin); 
    \endcode
    
    or a functor that is called (without argument) at every location,
    and the result is written into the current element. Internally,
    functors are recognized by the meta function 
    <tt>FunctorTraits&lt;FUNCTOR&gt;::</tt><tt>isInitializer</tt> yielding <tt>VigraTrueType</tt>.
    Make sure that your functor correctly defines <tt>FunctorTraits</tt> because
    otherwise the code will not compile.
    
    \code
    MultiIterator begin;    
    Accessor accessor;
    
    FUNCTOR f;
    assert(typeid(FunctorTraits<FUNCTOR>::isInitializer) == typeid(VigraTrueType));
    
    accessor.set(f(), begin); 
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
/*                    initVolumeBorder                  */
/*                                                      */
/********************************************************/

/** \brief Write value to the specified border voxels in the array.

   This is a simple port of <TT>initImageBorder</TT> to 3d for the 
   Seeded Region Growing algorithm

*/
template <class Iterator, class Diff_type, class Accessor, class VALUETYPE>
inline 
void
initVolumeBorder(Iterator upperleft, Diff_type shape, 
                Accessor a,  int border_width, VALUETYPE v)
{
    int w = shape[0]; 
    int h = shape[1]; 
    int d = shape[2]; 
    
    int hb = (border_width > h) ? h : border_width;
    int wb = (border_width > w) ? w : border_width;
    int db = (border_width > d) ? d : border_width;
    
    initMultiArray(upperleft, Diff_type(w,hb,d), a, v);
    initMultiArray(upperleft, Diff_type(wb,h,d), a, v);
    
    initMultiArray(upperleft+Diff_type(0,h-hb,0), Diff_type(w,hb,d), a, v);
    initMultiArray(upperleft+Diff_type(w-wb,0,0), Diff_type(wb,h,d), a, v);
    
    initMultiArray(upperleft, Diff_type(w,h,db), a, v);
    initMultiArray(upperleft+Diff_type(0,0,d-db), Diff_type(w,h,db), a, v);
}
    
template <class Iterator, class Diff_type, class Accessor, class VALUETYPE>
inline 
void
initVolumeBorder(triple<Iterator, Diff_type, Accessor> img, 
                int border_width, VALUETYPE v)
{
    initVolumeBorder(img.first, img.second, img.third, border_width, v);
}

    
/********************************************************/
/*                                                      */
/*                    copyMultiArray                    */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor>
void
copyMultiArrayImpl(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
               DestIterator d, DestShape const & dshape, DestAccessor dest, MetaInt<0>)
{
    if(sshape[0] == 1)
    {
        initLine(d, d + dshape[0], dest, src(s));
    }
    else
    {
        copyLine(s, s + sshape[0], src, d, dest);
    }
}
    
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor, int N>
void
copyMultiArrayImpl(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
                   DestIterator d, DestShape const & dshape, DestAccessor dest, MetaInt<N>)
{
    DestIterator dend = d + dshape[N];
    if(sshape[N] == 1)
    {
        for(; d != dend; ++d)
        {
            copyMultiArrayImpl(s.begin(), sshape, src, d.begin(), dshape, dest, MetaInt<N-1>());
        }
    }
    else
    {
        for(; d != dend; ++s, ++d)
        {
            copyMultiArrayImpl(s.begin(), sshape, src, d.begin(), dshape, dest, MetaInt<N-1>());
        }
    }
}
    
/** \brief Copy a multi-dimensional array.

    This function can be applied in two modes:
    
    <DL>
    <DT><b>Standard Mode:</b>
        <DD>If the source and destination arrays have the same size, 
        the corresponding array elements are simply copied.
        If necessary, type conversion takes place.
    <DT><b>Expanding Mode:</b>
        <DD>If the source array has length 1 along some (or even all) dimensions,
        the source value at index 0 is used for all destination
        elements in those dimensions. For example, if we have single row of data
        (column length is 1), we can copy it into a 2D image of the same width:
        The given row is automatically repeated for every row of the destination image.
        Again, type conversion os performed if necessary.
    </DL>
        
    The arrays must be represented by
    iterators compatible with \ref vigra::MultiIterator, and the iteration range 
    is specified by means of shape objects. If only the source shape is given
    the destination array is assumed to have the same shape, and standard mode
    is applied. If two shapes are given, the size of corresponding dimensions
    must be either equal (standard copy), or the source length must be 1 
    (expanding copy). The function uses accessors to access the data elements. 
    
    <b> Declarations:</b>
    
    <b>\#include</b> "<a href="multi__pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"<br>
    Namespace: vigra
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        copyMultiArray(SrcIterator s, 
                       SrcShape const & shape, SrcAccessor src,
                       DestIterator d, DestAccessor dest);


        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestShape, class DestAccessor>
        void
        copyMultiArray(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
                       DestIterator d, DestShape const & dshape, DestAccessor dest);
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
                       
                       
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestShape, class DestAccessor>
        void
        copyMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
                       triple<DestIterator, DestShape, DestAccessor> const & dest);
    }
    \endcode
    
    <b> Usage - Standard Mode:</b>
    
    \code
    typedef vigra::MultiArray<3, int> Array;
    Array src(Array::size_type(100, 200, 50)),
          dest(Array::size_type(100, 200, 50));
    ...
    
    vigra::copyMultiArray(srcMultiArrayRange(src), destMultiArray(dest));
    \endcode

    <b> Usage - Expanding Mode:</b>
    
    The source array is only 2D (it has depth 1). Thus, the destination
    will contain 50 identical copies of this image. Note that the destination shape
    must be passed to the algorithm for the expansion to work, so we use 
    <tt>destMultiArrayRange()</tt> rather than <tt>destMultiArray()</tt>.
    
    \code
    typedef vigra::MultiArray<3, int> Array;
    Array src(Array::size_type(100, 200, 1)),
          dest(Array::size_type(100, 200, 50));
    ...
    
    vigra::copyMultiArray(srcMultiArrayRange(src), destMultiArrayRange(dest));
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
    copyMultiArrayImpl(s, shape, src, d, shape, dest, MetaInt<SrcIterator::level>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
copyMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               pair<DestIterator, DestAccessor> const & dest)
{
    
    copyMultiArray(src.first, src.second, src.third, dest.first, dest.second);
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor>
void
copyMultiArray(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
               DestIterator d, DestShape const & dshape, DestAccessor dest)
{    
    vigra_precondition(sshape.size() == dshape.size(),
        "copyMultiArray(): dimensionality of source and destination array differ");
    for(unsigned int i=0; i<sshape.size(); ++i)
        vigra_precondition(sshape[i] == 1 || sshape[i] == dshape[i],
            "copyMultiArray(): mismatch between source and destination shapes:\n"
            "length of each source dimension must either be 1 or equal to the corresponding "
            "destination length.");
    copyMultiArrayImpl(s, sshape, src, d, dshape, dest, MetaInt<SrcIterator::level>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor>
inline void
copyMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               triple<DestIterator, DestShape, DestAccessor> const & dest)
{
    
    copyMultiArray(src.first, src.second, src.third, dest.first, dest.second, dest.third);
}

/********************************************************/
/*                                                      */
/*                 transformMultiArray                  */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
void
transformMultiArrayReduceImpl(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               SrcShape const & reduceShape,
               Functor const & ff, MetaInt<0>)
{
    DestIterator dend = d + dshape[0];
    for(; d != dend; ++s.template dim<0>(), ++d)
    {
        Functor f = ff;
        inspectMultiArray(s, reduceShape, src, f);
        dest.set(f(), d);
    }
}
    
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor, int N>
void
transformMultiArrayReduceImpl(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
                   DestIterator d, DestShape const & dshape, DestAccessor dest, 
                   SrcShape const & reduceShape,
                   Functor const & f, MetaInt<N>)
{
    DestIterator dend = d + dshape[N];
    for(; d != dend; ++s.template dim<N>(), ++d)
    {
        transformMultiArrayReduceImpl(s, sshape, src, d.begin(), dshape, dest,
                                      reduceShape, f, MetaInt<N-1>());
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
void
transformMultiArrayImpl(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               Functor const & f, VigraTrueType)
{
    // reduce mode
    SrcShape reduceShape = sshape;
    for(unsigned int i=0; i<dshape.size(); ++i)
    {
        vigra_precondition(dshape[i] == 1 || sshape[i] == dshape[i],
            "transformMultiArray(): mismatch between source and destination shapes:\n"
            "In 'reduce'-mode, the length of each destination dimension must either be 1\n"
            "or equal to the corresponding source length.");
        if(dshape[i] != 1)
            reduceShape[i] = 1;
    }
    transformMultiArrayReduceImpl(s, sshape, src, d, dshape, dest, reduceShape,
                                  f, MetaInt<SrcIterator::level>());
}
    
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
void
transformMultiArrayExpandImpl(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               Functor const & f, MetaInt<0>)
{
    if(sshape[0] == 1)
    {
        initLine(d, d + dshape[0], dest, f(src(s)));
    }
    else
    {
        transformLine(s, s + sshape[0], src, d, dest, f);
    }
}
    
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor, int N>
void
transformMultiArrayExpandImpl(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
                   DestIterator d, DestShape const & dshape, DestAccessor dest, 
                   Functor const & f, MetaInt<N>)
{
    DestIterator dend = d + dshape[N];
    if(sshape[N] == 1)
    {
        for(; d != dend; ++d)
        {
            transformMultiArrayExpandImpl(s.begin(), sshape, src, d.begin(), dshape, dest,
                                          f, MetaInt<N-1>());
        }
    }
    else
    {
        for(; d != dend; ++s, ++d)
        {
            transformMultiArrayExpandImpl(s.begin(), sshape, src, d.begin(), dshape, dest,
                                          f, MetaInt<N-1>());
        }
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
void
transformMultiArrayImpl(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               Functor const & f, VigraFalseType)
{
    // expand mode
    for(unsigned int i=0; i<sshape.size(); ++i)
        vigra_precondition(sshape[i] == 1 || sshape[i] == dshape[i],
            "transformMultiArray(): mismatch between source and destination shapes:\n"
            "In 'expand'-mode, the length of each source dimension must either be 1\n"
            "or equal to the corresponding destination length.");
    transformMultiArrayExpandImpl(s, sshape, src, d, dshape, dest, 
                                  f, MetaInt<SrcIterator::level>());
}
    
/** \brief Transform a multi-dimensional array with a unary function or functor.

    This function can be applied in three modes:
    
    <DL>
    <DT><b>Standard Mode:</b>
        <DD>If the source and destination arrays have the same size, 
        the transformation given by the functor is applied to every source
        element and the result written into the corresponding destination element.
        Unary functions, unary functors from the STL and the functors specifically 
        defined in \ref TransformFunctor can be used in standard mode.
        Creation of new functors is easiest by using \ref FunctorExpressions. 
    <DT><b>Expanding Mode:</b>
        <DD>If the source array has length 1 along some (or even all) dimensions,
        the source value at index 0 is used for all destination
        elements in those dimensions. In other words, the source index is not
        incremented along these dimensions, but the transformation functor
        is applied as usual. So, we can expand a small array (e.g. a single row of data,
        column length is 1), into a larger one (e.g. a 2D image with the same width): 
        the given values are simply reused as necessary (e.g. for every row of the 
        destination image). The same functors as in standard mode can be applied.
    <DT><b>Reducing Mode:</b>
        <DD>If the destination array has length 1 along some (or even all) dimensions,
        the source values in these dimensions are reduced to single values by means
        of a suitable functor (e.g. \ref vigra::ReduceFunctor), which supports two 
        function call operators: one
        with a single argument to collect the values, and without argument to 
        obtain the final (reduced) result. This behavior is a multi-dimensional
        generalization of the C++ standard function <tt>std::accumulate()</tt>.
    </DL>
        
    The arrays must be represented by
    iterators compatible with \ref vigra::MultiIterator, and the iteration range 
    is specified by means of shape objects. If only the source shape is given
    the destination array is assumed to have the same shape, and standard mode
    is applied. If two shapes are given, the size of corresponding dimensions
    must be either equal (standard copy), or the source length must be 1 
    (expand mode), or the destination length must be 1 (reduce mode). However,
    reduction and expansion cannot be executed at the same time, so the latter
    conditions are mutual exclusive, even if they apply to different dimensions.
    
    The function uses accessors to access the data elements. 
    
    <b> Declarations:</b>

    <b>\#include</b> "<a href="multi__pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"<br>
    Namespace: vigra
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, 
                  class Functor>
        void
        transformMultiArray(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                            DestIterator d, DestAccessor dest, Functor const & f);


        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestShape, class DestAccessor, 
                  class Functor>
        void
        transformMultiArray(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
                            DestIterator d, DestShape const & dshape, DestAccessor dest, 
                            Functor const & f);
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


        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestShape, class DestAccessor, 
                  class Functor>
        void
        transformMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
                            triple<DestIterator, DestShape, DestAccessor> const & dest, 
                            Functor const & f)
    }
    \endcode

    <b> Usage - Standard Mode:</b>

    Source and destination array have the same size.
    
    \code
    #include <cmath>         // for sqrt()

    typedef vigra::MultiArray<3, float> Array;
    Array src(Array::size_type(100, 200, 50)),
          dest(Array::size_type(100, 200, 50));
    ...
    
    vigra::transformMultiArray(srcMultiArrayRange(src),
                               destMultiArray(dest),
                               (float(*)(float))&std::sqrt );

    \endcode

    <b> Usage - Expand Mode:</b>

    The source array is only 2D (it has depth 1). Thus, the destination
    will contain 50 identical copies of the transformed source array. 
    Note that the destination shape must be passed to the algorithm for 
    the expansion to work, so we use <tt>destMultiArrayRange()</tt> 
    rather than <tt>destMultiArray()</tt>.
    
    \code
    #include <cmath>         // for sqrt()

    typedef vigra::MultiArray<3, float> Array;
    Array src(Array::size_type(100, 200, 1)),
          dest(Array::size_type(100, 200, 50));
    ...
    
    vigra::transformMultiArray(srcMultiArrayRange(src),
                               destMultiArrayRange(dest),
                               (float(*)(float))&std::sqrt );

    \endcode

    <b> Usage - Reduce Mode:</b>

    The destination array is only 1D (it's width and height are 1). 
    Thus, it will contain accumulated data for every slice of the source volume
    (or for every frame, if the source is intepreted as an image sequence).
    In the example, we use the functor \ref vigra::FindAverage to calculate
    the average gray value of every slice. Note that the destination shape
    must also be passed for the reduction to work, so we use 
    <tt>destMultiArrayRange()</tt> rather than <tt>destMultiArray()</tt>.
    
    \code
    typedef vigra::MultiArray<3, float> Array;
    Array src(Array::size_type(100, 200, 50)),
          dest(Array::size_type(1, 1, 50));
    ...
    
    vigra::transformMultiArray(srcMultiArrayRange(src),
                               destMultiArrayRange(dest),
                               vigra::FindAverage<float>() );

    \endcode

    <b> Required Interface:</b>

    In standard and expand mode, the functor must be a model of UnaryFunction
    (i.e. support function call with one argument and a return value
    <tt>res = functor(arg)</tt>):
    
    \code
    MultiIterator src_begin, src_end, dest_begin;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    Functor functor;

    dest_accessor.set(functor(src_accessor(src_begin)), dest_begin);
    \endcode
    
    In reduce mode, it must be a model of UnaryAnalyser (i.e. support function call
    with one argument and no return vakue <tt>functor(arg)</tt>) and Initializer
    (i.e. support function call with no argument, but return value 
    <tt>res = functor()</tt>). Internally, such functors are recognized by the 
    meta functions <tt>FunctorTraits&lt;FUNCTOR&gt;::</tt><tt>isUnaryAnalyser</tt> and
    <tt>FunctorTraits&lt;FUNCTOR&gt;::</tt><tt>isInitializer</tt> which must both yield 
    <tt>VigraTrueType</tt>. Make sure that your functor correctly defines 
    <tt>FunctorTraits</tt> because otherwise reduce mode will not work. In addition,
    the functor must be copy constructible in order to start each reduction
    with a fresh functor.
    
    \code
    MultiIterator src_begin, src_end, dest_begin;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    FUNCTOR initial_functor, functor(initial_functor);
    assert(typeid(FunctorTraits<FUNCTOR>::isInitializer) == typeid(VigraTrueType));
    assert(typeid(FunctorTraits<FUNCTOR>::isUnaryAnalyser) == typeid(VigraTrueType));
    
    functor(src_accessor(src_begin));
    dest_accessor.set(functor(), dest_begin);
    \endcode

*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          class Functor>
inline void
transformMultiArray(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                    DestIterator d, DestAccessor dest, Functor const & f)
{    
    transformMultiArrayExpandImpl(s, shape, src, d, shape, dest, 
                                  f, MetaInt<SrcIterator::level>());
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

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
void
transformMultiArray(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               Functor const & f)
{    
    vigra_precondition(sshape.size() == dshape.size(),
        "transformMultiArray(): dimensionality of source and destination array differ");
    typedef FunctorTraits<Functor> FT;
    typedef typename 
        And<typename FT::isInitializer, typename FT::isUnaryAnalyser>::result
        isAnalyserInitializer;
    transformMultiArrayImpl(s, sshape, src, d, dshape, dest, 
                            f, isAnalyserInitializer());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
inline void
transformMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               triple<DestIterator, DestShape, DestAccessor> const & dest, 
               Functor const & f)
{
    transformMultiArray(src.first, src.second, src.third, 
                        dest.first, dest.second, dest.third, f);
}

/********************************************************/
/*                                                      */
/*                combineTwoMultiArrays                 */
/*                                                      */
/********************************************************/

template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
void
combineTwoMultiArraysReduceImpl(
               SrcIterator1 s1, SrcShape const & sshape, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               DestIterator d,  DestShape const & dshape, DestAccessor dest, 
               SrcShape const & reduceShape,
               Functor const & ff, MetaInt<0>)
{
    DestIterator dend = d + dshape[0];
    for(; d != dend; ++s1.template dim<0>(), ++s2.template dim<0>(), ++d)
    {
        Functor f = ff;
        inspectTwoMultiArrays(s1, reduceShape, src1, s2, src2, f);
        dest.set(f(), d);
    }
}
    
template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor, int N>
void
combineTwoMultiArraysReduceImpl(
               SrcIterator1 s1, SrcShape const & sshape, SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               DestIterator d,  DestShape const & dshape, DestAccessor dest, 
               SrcShape const & reduceShape,
               Functor const & f, MetaInt<N>)
{
    DestIterator dend = d + dshape[N];
    for(; d != dend; ++s1.template dim<N>(), ++s2.template dim<N>(), ++d)
    {
        combineTwoMultiArraysReduceImpl(s1, sshape, src1, s2, src2, 
                                        d.begin(), dshape, dest,
                                        reduceShape, f, MetaInt<N-1>());
    }
}

template <class SrcIterator1, class SrcShape1, class SrcAccessor1,
          class SrcIterator2, class SrcShape2, class SrcAccessor2,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
void
combineTwoMultiArraysImpl(
               SrcIterator1 s1, SrcShape1 const & sshape1, SrcAccessor1 src1,
               SrcIterator2 s2, SrcShape2 const & sshape2, SrcAccessor2 src2,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               Functor const & f, VigraTrueType)
{
    // reduce mode
    SrcShape1 reduceShape = sshape1;
    for(unsigned int i=0; i<dshape.size(); ++i)
    {
        vigra_precondition(sshape1[i] == sshape2[i] && 
                           (dshape[i] == 1 || sshape1[i] == dshape[i]),
            "combineTwoMultiArrays(): mismatch between source and destination shapes:\n"
            "In 'reduce'-mode, the two source shapes must be equal, and\n"
            "the length of each destination dimension must either be 1\n"
            "or equal to the corresponding source length.");
        if(dshape[i] != 1)
            reduceShape[i] = 1;
    }
    combineTwoMultiArraysReduceImpl(s1, sshape1, src1, s2, src2, 
                                    d, dshape, dest, reduceShape,
                                    f, MetaInt<SrcIterator1::level>());
}
    
template <class SrcIterator1, class SrcShape1, class SrcAccessor1,
          class SrcIterator2, class SrcShape2, class SrcAccessor2,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
void
combineTwoMultiArraysExpandImpl(
               SrcIterator1 s1, SrcShape1 const & sshape1, SrcAccessor1 src1,
               SrcIterator2 s2, SrcShape2 const & sshape2, SrcAccessor2 src2,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               Functor const & f, MetaInt<0>)
{
    DestIterator dend = d + dshape[0];
    if(sshape1[0] == 1 && sshape2[0] == 1)
    {
        initLine(d, dend, dest, f(src1(s1), src2(s2)));
    }
    else if(sshape1[0] == 1)
    {
        typename SrcAccessor1::value_type sv1 = src1(s1);
        for(; d != dend; ++d, ++s2)
            dest.set(f(sv1, src2(s2)), d);
    }
    else if(sshape2[0] == 1)
    {
        typename SrcAccessor2::value_type sv2 = src2(s2);
        for(; d != dend; ++d, ++s1)
            dest.set(f(src1(s1), sv2), d);
    }
    else
    {
        combineTwoLines(s1, s1 + sshape1[0], src1, s2, src2, d, dest, f);
    }
}
    
template <class SrcIterator1, class SrcShape1, class SrcAccessor1,
          class SrcIterator2, class SrcShape2, class SrcAccessor2,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor, int N>
void
combineTwoMultiArraysExpandImpl(
               SrcIterator1 s1, SrcShape1 const & sshape1, SrcAccessor1 src1,
               SrcIterator2 s2, SrcShape2 const & sshape2, SrcAccessor2 src2,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               Functor const & f, MetaInt<N>)
{
    DestIterator dend = d + dshape[N];
    int s1inc = sshape1[N] == 1
                    ? 0 
                    : 1;
    int s2inc = sshape2[N] == 1
                    ? 0 
                    : 1;
    for(; d != dend; ++d, s1 += s1inc, s2 += s2inc)
    {
        combineTwoMultiArraysExpandImpl(s1.begin(), sshape1, src1, 
                                        s2.begin(), sshape2, src2, 
                                        d.begin(), dshape, dest,
                                        f, MetaInt<N-1>());
    }
}

template <class SrcIterator1, class SrcShape1, class SrcAccessor1,
          class SrcIterator2, class SrcShape2, class SrcAccessor2,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
void
combineTwoMultiArraysImpl(
               SrcIterator1 s1, SrcShape1 const & sshape1, SrcAccessor1 src1,
               SrcIterator2 s2, SrcShape2 const & sshape2, SrcAccessor2 src2,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               Functor const & f, VigraFalseType)
{
    // expand mode
    for(unsigned int i=0; i<sshape1.size(); ++i)
        vigra_precondition((sshape1[i] == 1 || sshape1[i] == dshape[i]) &&
                           (sshape2[i] == 1 || sshape2[i] == dshape[i]),
            "combineTwoMultiArrays(): mismatch between source and destination shapes:\n"
            "In 'expand'-mode, the length of each source dimension must either be 1\n"
            "or equal to the corresponding destination length.");
    combineTwoMultiArraysExpandImpl(s1, sshape1, src1, s2, sshape2, src2, 
                                    d, dshape, dest, 
                                    f, MetaInt<SrcIterator1::level>());
}

/** \brief Combine two multi-dimensional arrays into one using a binary function or functor.

    This function can be applied in three modes:
    
    <DL>
    <DT><b>Standard Mode:</b>
        <DD>If the source and destination arrays have the same size, 
        the transformation given by the functor is applied to every pair of
        corresponding source elements and the result written into the corresponding 
        destination element.
        Binary functions, binary functors from the STL and the functors specifically 
        defined in \ref CombineFunctor can be used in standard mode.
        Creation of new functors is easiest by using \ref FunctorExpressions. 
    <DT><b>Expanding Mode:</b>
        <DD>If the source arrays have length 1 along some (or even all) dimensions,
        the source values at index 0 are used for all destination
        elements in those dimensions. In other words, the source index is not
        incremented along those dimensions, but the transformation functor
        is applied as usual. So, we can expand small arrays (e.g. a single row of data,
        column length is 1), into larger ones (e.g. a 2D image with the same width): 
        the given values are simply reused as necessary (e.g. for every row of the 
        destination image). It is not even necessary that the source array shapes
        are equal. For example, we can combine a small array with one that
        hase the same size as the destination array. 
        The same functors as in standard mode can be applied.
    <DT><b>Reducing Mode:</b>
        <DD>If the destination array has length 1 along some (or even all) dimensions,
        the source values in these dimensions are reduced to single values by means
        of a suitable functor which supports two function call operators: one
        with two arguments to collect the values, and one without argument to 
        obtain the final (reduced) result. This behavior is a multi-dimensional
        generalization of the C++ standard function <tt>std::accumulate()</tt>.
    </DL>
        
    The arrays must be represented by
    iterators compatible with \ref vigra::MultiIterator, and the iteration range 
    is specified by means of shape objects. If only a single source shape is given
    the destination array is assumed to have the same shape, and standard mode
    is applied. If three shapes are given, the size of corresponding dimensions
    must be either equal (standard copy), or the length of this dimension must
    be 1 in one or both source arrays
    (expand mode), or the destination length must be 1 (reduce mode). However,
    reduction and expansion cannot be executed at the same time, so the latter
    conditions are mutual exclusive, even if they apply to different dimensions.
    
    The function uses accessors to access the data elements. 
    
    <b> Declarations:</b>
    
    <b>\#include</b> "<a href="multi__pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"<br>
    Namespace: vigra
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator1, class SrcShape, class SrcAccessor1,
                  class SrcIterator2, class SrcAccessor2,
                  class DestIterator, class DestAccessor, 
                  class Functor>
        void combineTwoMultiArrays(
                       SrcIterator1 s1, SrcShape const & shape, SrcAccessor1 src1,
                       SrcIterator2 s2, SrcAccessor2 src2,
                       DestIterator d, DestAccessor dest, Functor const & f);


        template <class SrcIterator1, class SrcShape1, class SrcAccessor1,
                  class SrcIterator2, class SrcShape2, class SrcAccessor2,
                  class DestIterator, class DestShape, class DestAccessor, 
                  class Functor>
        void combineTwoMultiArrays(
                       SrcIterator1 s1, SrcShape1 const & sshape1, SrcAccessor1 src1,
                       SrcIterator2 s2, SrcShape2 const & sshape2, SrcAccessor2 src2,
                       DestIterator d, DestShape const & dshape, DestAccessor dest, 
                       Functor const & f);
            }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator1, class SrcShape, class SrcAccessor1,
                  class SrcIterator2, class SrcAccessor2,
                  class DestIterator, class DestAccessor, class Functor>
        void combineTwoMultiArrays(
                       triple<SrcIterator1, SrcShape, SrcAccessor1> const & src1,
                       pair<SrcIterator2, SrcAccessor2> const & src2,
                       pair<DestIterator, DestAccessor> const & dest, Functor const & f);


        template <class SrcIterator1, class SrcShape1, class SrcAccessor1,
                  class SrcIterator2, class SrcShape2, class SrcAccessor2,
                  class DestIterator, class DestShape, class DestAccessor, 
                  class Functor>
        void combineTwoMultiArrays(
                       triple<SrcIterator1, SrcShape1, SrcAccessor1> const & src1,
                       triple<SrcIterator2, SrcShape2, SrcAccessor2> const & src2,
                       triple<DestIterator, DestShape, DestAccessor> const & dest, 
                       Functor const & f);
    }
    \endcode
    
    <b> Usage - Standard Mode:</b>
    
    Source and destination arrays have the same size.
    
    \code
    #include <functional>     // for std::plus

    typedef vigra::MultiArray<3, int> Array;
    Array src1(Array::size_type(100, 200, 50)),
          src2(Array::size_type(100, 200, 50)),
          dest(Array::size_type(100, 200, 50));
    ...
    
    vigra::combineTwoMultiArrays(
                srcMultiArrayRange(src1), 
                srcMultiArray(src2), 
                destMultiArray(dest),  
                std::plus<int>());
    
    \endcode
    
    <b> Usage - Expand Mode:</b>

    One source array is only 2D (it has depth 1). This image will be added
    to every slice of the other source array, and the result
    if written into the corresponding destination slice. Note that the shapes
    of all arrays must be passed to the algorithm, so we use 
    <tt>srcMultiArrayRange()</tt> and <tt>destMultiArrayRange()</tt> 
    rather than <tt>srcMultiArray()</tt> and <tt>destMultiArray()</tt>.
    
    \code
    #include <functional>     // for std::plus

    typedef vigra::MultiArray<3, int> Array;
    Array src1(Array::size_type(100, 200, 1)),
          src2(Array::size_type(100, 200, 50)),
          dest(Array::size_type(100, 200, 50));
    ...
    
    vigra::combineTwoMultiArrays(
                srcMultiArrayRange(src1), 
                srcMultiArray(src2), 
                destMultiArray(dest),  
                std::plus<int>());

    \endcode

    <b> Usage - Reduce Mode:</b>

    The destination array is only 1D (it's width and height are 1). 
    Thus, it will contain accumulated data for every slice of the source volumes
    (or for every frame, if the sources are intepreted as image sequences).
    In the example, we use \ref vigra::ReduceFunctor together with a functor 
    expression (see \ref FunctorExpressions)
    to calculate the total absolute difference of the gray values in every pair of 
    source slices. Note that the shapes of all arrays must be passed 
    to the algorithm in order for the reduction to work, so we use 
    <tt>srcMultiArrayRange()</tt> and <tt>destMultiArrayRange()</tt> 
    rather than <tt>srcMultiArray()</tt> and <tt>destMultiArray()</tt>.
    
    \code
    #include <vigra/functorexpression.hxx>
    using namespace vigra::functor;
        
    typedef vigra::MultiArray<3, int> Array;
    Array src1(Array::size_type(100, 200, 50)),
          src2(Array::size_type(100, 200, 50)),
          dest(Array::size_type(1, 1, 50));
    ...
    
    vigra::combineTwoMultiArrays(
                srcMultiArrayRange(src1), 
                srcMultiArray(src2), 
                destMultiArray(dest),  
                reduceFunctor(Arg1() + abs(Arg2() - Arg3()), 0) );
                // Arg1() is the sum accumulated so far, initialzed with 0

    \endcode

    <b> Required Interface:</b>
    
    In standard and expand mode, the functor must be a model of BinaryFunction
    (i.e. support function call with two arguments and a return value
    <tt>res = functor(arg1, arg2)</tt>):
    
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
        
    In reduce mode, it must be a model of BinaryAnalyser (i.e. support function call
    with two arguments and no return vakue <tt>functor(arg1, arg2)</tt>) and Initializer
    (i.e. support function call with no argument, but return value 
    <tt>res = functor()</tt>). Internally, such functors are recognized by the 
    meta functions <tt>FunctorTraits&lt;FUNCTOR&gt;::</tt><tt>isBinaryAnalyser</tt> and
    <tt>FunctorTraits&lt;FUNCTOR&gt;::</tt><tt>isInitializer</tt> which must both yield 
    <tt>VigraTrueType</tt>. Make sure that your functor correctly defines 
    <tt>FunctorTraits</tt> because otherwise reduce mode will not work. In addition,
    the functor must be copy constructible in order to start each reduction
    with a fresh functor.
    
    \code
    MultiIterator src1_begin, src2_begin, dest_begin;
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    DestAccessor dest_accessor;
    
    FUNCTOR initial_functor, functor(initial_functor);
    assert(typeid(FunctorTraits<FUNCTOR>::isInitializer) == typeid(VigraTrueType));
    assert(typeid(FunctorTraits<FUNCTOR>::isBinaryAnalyser) == typeid(VigraTrueType));
    
    functor(src1_accessor(src1_begin), src2_accessor(src2_begin));
    dest_accessor.set(functor(), dest_begin);
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
    combineTwoMultiArraysExpandImpl(s1, shape, src1, s2, shape, src2, d, shape, dest, f, 
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

template <class SrcIterator1, class SrcShape1, class SrcAccessor1,
          class SrcIterator2, class SrcShape2, class SrcAccessor2,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
void
combineTwoMultiArrays(
               SrcIterator1 s1, SrcShape1 const & sshape1, SrcAccessor1 src1,
               SrcIterator2 s2, SrcShape2 const & sshape2, SrcAccessor2 src2,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               Functor const & f)
{    
    vigra_precondition(sshape1.size() == dshape.size() && sshape2.size() == dshape.size(),
        "combineTwoMultiArrays(): dimensionality of source and destination arrays differ");
    
    typedef FunctorTraits<Functor> FT;
    typedef typename 
        And<typename FT::isInitializer, typename FT::isBinaryAnalyser>::result
        isAnalyserInitializer;
    combineTwoMultiArraysImpl(s1, sshape1, src1, s2, sshape2, src2, d, dshape, dest, 
                              f, isAnalyserInitializer());
}

template <class SrcIterator1, class SrcShape1, class SrcAccessor1,
          class SrcIterator2, class SrcShape2, class SrcAccessor2,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
inline void
combineTwoMultiArrays(
               triple<SrcIterator1, SrcShape1, SrcAccessor1> const & src1,
               triple<SrcIterator2, SrcShape2, SrcAccessor2> const & src2,
               triple<DestIterator, DestShape, DestAccessor> const & dest, 
               Functor const & f)
{
    combineTwoMultiArrays(src1.first, src1.second, src1.third, 
                          src2.first, src2.second, src2.third, 
                          dest.first, dest.second, dest.third, f);
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
    
    
/** \brief Combine three multi-dimensional arrays into one using a 
           ternary function or functor.

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
    
    <b>\#include</b> "<a href="multi__pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"<br>
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
    
/** \brief Call an analyzing functor at every element of a multi-dimensional array.

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

    <b>\#include</b> "<a href="multi__pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"<br>
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
    
/** \brief Call an analyzing functor at all corresponding elements of 
           two multi-dimensional arrays.

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

    <b>\#include</b> "<a href="multi__pointoperators_8hxx-source.html">vigra/multi_pointoperators.hxx</a>"<br>
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

}	//-- namespace vigra


#endif	//-- VIGRA_MULTI_POINTOPERATORS_H
