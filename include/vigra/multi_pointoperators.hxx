//-- -*- c++ -*-
/************************************************************************/
/*                                                                      */
/*      Copyright 2003 by Ullrich Koethe, B. Seppke, F. Heinrich        */
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

#ifndef VIGRA_MULTI_POINTOPERATORS_H
#define VIGRA_MULTI_POINTOPERATORS_H

#include "initimage.hxx"
#include "copyimage.hxx"
#include "transformimage.hxx"
#include "combineimages.hxx"
#include "inspectimage.hxx"
#include "multi_array.hxx"
#include "metaprogramming.hxx"
#include "inspector_passes.hxx"



namespace vigra
{

/** \addtogroup MultiPointoperators Point operators for multi-dimensional arrays.

    Copy, transform, and inspect arbitrary dimensional arrays which are represented
    by iterators compatible to \ref MultiIteratorPage. Note that are range is here
    specified by a pair: an iterator referring to the first point of the array 
    and a shape object specifying the size of the (rectangular) ROI.

    <b>\#include</b> \<vigra/multi_pointoperators.hxx\><br/>
    Namespace: vigra
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
initMultiArrayImpl(Iterator s, Shape const & shape, Accessor a,  VALUETYPE const & v, MetaInt<0>)
{
    initLine(s, s + shape[0], a, v);
}
    
template <class Iterator, class Shape, class Accessor, 
          class VALUETYPE, int N>
void
initMultiArrayImpl(Iterator s, Shape const & shape, Accessor a,  
                   VALUETYPE const & v, MetaInt<N>)
{
    Iterator send = s + shape[N];
    for(; s < send; ++s)
    {
        initMultiArrayImpl(s.begin(), shape, a, v, MetaInt<N-1>());
    }
}
    
/** \brief Write a value to every element in a multi-dimensional array.

    The initial value can either be a constant of appropriate type (compatible with 
    the destination's value_type), or a functor with compatible result_type. These two 
    cases are automatically distinguished when <tt>FunctorTraits<FUNCTOR>::isInitializer</tt>
    yields <tt>VigraTrueType</tt>. Since the functor is passed by <tt>const</tt> reference, its 
    <tt>operator()</tt> must be const, and its internal state may need to be <tt>mutable</tt>.
    
    <b> Declarations:</b>
    
    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T, class S, class VALUETYPE>
        void
        initMultiArray(MultiArrayView<N, T, S> s, VALUETYPE const & v);
        
        template <unsigned int N, class T, class S, class FUNCTOR>
        void
        initMultiArray(MultiArrayView<N, T, S> s, FUNCTOR const & f);
    }
    \endcode
    
    \deprecatedAPI{initMultiArray}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
    \code
    namespace vigra {
        template <class Iterator, class Shape, class Accessor, class VALUETYPE>
        void
        initMultiArray(Iterator s, Shape const & shape, Accessor a,  VALUETYPE const & v);

        template <class Iterator, class Shape, class Accessor, class FUNCTOR>
        void
        initMultiArray(Iterator s, Shape const & shape, Accessor a,  FUNCTOR const & f);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class Iterator, class Shape, class Accessor, class VALUETYPE>
        void
        initMultiArray(triple<Iterator, Shape, Accessor> const & s, VALUETYPE const & v);


        template <class Iterator, class Shape, class Accessor, class FUNCTOR>
        void
        initMultiArray(triple<Iterator, Shape, Accessor> const & s, FUNCTOR const & f);
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/multi_pointoperators.hxx\><br>
    Namespace: vigra
    
    \code
    MultiArray<3, unsigned int> array(Shape3(100, 200, 50));
    
    // make an array of all ones
    initMultiArray(array, 1);
    
    // equivalent calls:
    array = 1;
    array.init(1);
    
    // fill the array with random numbers
    #include <vigra/random.hxx> 
    
    initMultiArray(array, MersenneTwister());
    \endcode

    \deprecatedUsage{initMultiArray}
    \code
    MultiArray<3, int> array(Shape3(100, 200, 50));
    
    // make an array of all twos
    vigra::initMultiArray(destMultiArrayRange(array), 2);
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
    <tt>FunctorTraits<FUNCTOR>::isInitializer</tt> yielding <tt>VigraTrueType</tt>.
    Make sure that your functor correctly defines <tt>FunctorTraits</tt> because
    otherwise the code will not compile.
    \code
    MultiIterator begin;    
    Accessor accessor;
    
    FUNCTOR f;
    assert(typeid(FunctorTraits<FUNCTOR>::isInitializer) == typeid(VigraTrueType));
    
    accessor.set(f(), begin); 
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void initMultiArray)

template <class Iterator, class Shape, class Accessor, class VALUETYPE>
inline void
initMultiArray(Iterator s, Shape const & shape, Accessor a,  VALUETYPE const & v)
{
    initMultiArrayImpl(s, shape, a, v, MetaInt<Iterator::level>());
}
    
template <class Iterator, class Shape, class Accessor, class VALUETYPE>
inline void
initMultiArray(triple<Iterator, Shape, Accessor> const & s, VALUETYPE const & v)
{
     initMultiArrayImpl(s.first, s.second, s.third, v, MetaInt<Iterator::level>());
}

template <unsigned int N, class T, class S, class VALUETYPE>
inline void
initMultiArray(MultiArrayView<N, T, S> s, VALUETYPE const & v)
{
    initMultiArray(destMultiArrayRange(s), v);
}

/********************************************************/
/*                                                      */
/*                  initMultiArrayBorder                */
/*                                                      */
/********************************************************/

/** \brief Write values to the specified border values in the array.

    This functions is similar to \ref initMultiArray(), but it initializes only 
    the array elements whose distance from any array border is at most \a border_width.
    
    <b> Declarations:</b>
    
    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
            // init equal borders on all array sides
        template <unsigned int N, class T, class S, 
                  class VALUETYPE>
        void 
        initMultiArrayBorder( MultiArrayView<N, T, S> array, 
                              MultiArrayIndex border_width, VALUETYPE const & v);
        
        template <unsigned int N, class T, class S, 
                  class FUNCTOR>
        void 
        initMultiArrayBorder( MultiArrayView<N, T, S> array, 
                              MultiArrayIndex border_width, FUNCTOR const & v);
        
            // specify border width individually for all array sides
        template <unsigned int N, class T, class S, 
                  class VALUETYPE>
        void 
        initMultiArrayBorder( MultiArrayView<N, T, S> array, 
                              typename MultiArrayShape<N>::type const & lower_border, 
                              typename MultiArrayShape<N>::type const & upper_border, 
                              VALUETYPE const & v);
    }
    \endcode
    
    \deprecatedAPI{initMultiArrayBorder}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
    \code
    namespace vigra {
        template <class Iterator, class Diff_type, class Accessor, 
                  class VALUETYPE>
        void
        initMultiArrayBorder(Iterator upperleft, Diff_type shape, Accessor a,
                             MultiArrayIndex border_width, VALUETYPE const & v);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class Iterator, class Diff_type, class Accessor, 
                  class VALUETYPE>
        inline void 
        initMultiArrayBorder( triple<Iterator, Diff_type, Accessor> multiArray, 
                              MultiArrayIndex border_width, VALUETYPE const & v);
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/multi_pointoperators.hxx\><br>
    Namespace: vigra
    
    \code
    MultiArray<3, unsigned int> array(Shape3(100, 200, 50));
    
    int border_width = 5;
    
    // init the array interior to 1, the border to 2
    initMultiArray(array.subarray(Shape3(border_width), Shape3(-border_width)), 1);
    initMultiArrayBorder(array, border_width, 2);
    \endcode
*/
doxygen_overloaded_function(template <...> void initMultiArrayBorder)

template <class Iterator, class Diff_type, class Accessor, 
          class VALUETYPE>
void
initMultiArrayBorder(Iterator upperleft, Diff_type shape, Accessor a,
                     Diff_type lower_border, Diff_type upper_border, 
                     VALUETYPE const & v)
{
    for(unsigned int dim=0; dim<shape.size(); dim++)
    {
        lower_border[dim] = (lower_border[dim] > shape[dim]) ? shape[dim] : lower_border[dim];
        upper_border[dim] = (upper_border[dim] > shape[dim]) ? shape[dim] : upper_border[dim];
    }

    for(unsigned int dim=0; dim<shape.size(); dim++)
    {
        Diff_type  start,
                   offset(shape);
        offset[dim] = lower_border[dim];

        initMultiArray(upperleft+start, offset, a, v);
 
        start[dim]  = shape[dim] - upper_border[dim];
        offset[dim] = upper_border[dim];
        initMultiArray(upperleft+start, offset, a, v);
    }
}
    
template <class Iterator, class Diff_type, class Accessor, 
          class VALUETYPE>
inline void
initMultiArrayBorder(Iterator upperleft, Diff_type shape, Accessor a,
                     MultiArrayIndex border_width, VALUETYPE const & v)
{
    initMultiArrayBorder(upperleft, shape, a,
                         Diff_type(border_width), Diff_type(border_width), v);
}
    
template <class Iterator, class Diff_type, class Accessor, 
          class VALUETYPE>
inline void 
initMultiArrayBorder( triple<Iterator, Diff_type, Accessor> multiArray, 
                      MultiArrayIndex border_width, VALUETYPE const & v)
{
    initMultiArrayBorder(multiArray.first, multiArray.second, multiArray.third, border_width, v);
}

template <class Iterator, class Diff_type, class Accessor, 
          class VALUETYPE>
inline void 
initMultiArrayBorder( triple<Iterator, Diff_type, Accessor> multiArray, 
                      Diff_type const & lower_border, Diff_type const & upper_border, 
                      VALUETYPE const & v)
{
    initMultiArrayBorder(multiArray.first, multiArray.second, multiArray.third, 
                         lower_border, upper_border, v);
}

template <unsigned int N, class T, class S, 
          class VALUETYPE>
inline void 
initMultiArrayBorder( MultiArrayView<N, T, S> array, 
                      MultiArrayIndex border_width, VALUETYPE const & v)
{
    initMultiArrayBorder(destMultiArrayRange(array), border_width, v);
}

template <unsigned int N, class T, class S, 
          class VALUETYPE>
inline void 
initMultiArrayBorder( MultiArrayView<N, T, S> array, 
                      typename MultiArrayShape<N>::type const & lower_border, 
                      typename MultiArrayShape<N>::type const & upper_border, 
                      VALUETYPE const & v)
{
    initMultiArrayBorder(destMultiArrayRange(array), lower_border, upper_border, v);
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
        for(; d < dend; ++d)
        {
            copyMultiArrayImpl(s.begin(), sshape, src, d.begin(), dshape, dest, MetaInt<N-1>());
        }
    }
    else
    {
        for(; d < dend; ++s, ++d)
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
    (expanding copy). 
    
    <b> Declarations:</b>
    
    <b>\#include</b> \<vigra/multi_pointoperators.hxx\><br>
    Namespace: vigra
    
    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T1, class S1,
                  class T2, class S2>
        void
        copyMultiArray(MultiArrayView<N, T1, S1> const & source,
                       MultiArrayView<N, T2, S2> dest);
    }
    \endcode
    
    \deprecatedAPI{copyMultiArray}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
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
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    \deprecatedEnd
    
    <b> Usage - Standard Mode:</b>
    
    \code
    MultiArray<3, int> src(Shape3(100, 200, 50)),
                       dest(Shape3(100, 200, 50));
    ...
    
    copyMultiArray(src, dest);
    
    // equivalent to
    dest = src;
    
    // copy only the red channel (i.e. channl 0) of an RGB array
    MultiArray<3, RGBValue<int> > rgb_src(Shape3(100, 200, 50));
    
    copyMultiArray(rgb_src.bindElementChannel(0), dest);
    
    // equivalent to 
    dest = rgb_src.bindElementChannel(0);
    \endcode

    <b> Usage - Expanding Mode:</b>
    
    The source array is effectively only a 2D image (it has a 3D shape, but 'depth' is a
    singleton dimension with length 1). Thus, the destination will contain 50 identical 
    copies of this image:
    
    \code
    MultiArray<3, int> src(Shape2(100, 200)),
                       dest(Shape3(100, 200, 50));
    ...
    
    copyMultiArray(src.insertSingletonDimension(2), dest);
    
    // create an RGB image with three identical color bands
    MultiArray<3, RGBValue<int> > rgb_dest(Shape2(100, 200));
    
    copyMultiArray(src.insertSingletonDimension(2), rgb_dest.expandElements(2));
    \endcode

    \deprecatedUsage{copyMultiArray}
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
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void copyMultiArray)

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
          class DestIterator, class DestAccessor>
inline void
copyMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               pair<DestIterator, DestAccessor> const & dest)
{
    
    copyMultiArray(src.first, src.second, src.third, dest.first, dest.second);
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor>
inline void
copyMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               triple<DestIterator, DestShape, DestAccessor> const & dest)
{
    
    copyMultiArray(src.first, src.second, src.third, dest.first, dest.second, dest.third);
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2>
inline void
copyMultiArray(MultiArrayView<N, T1, S1> const & source,
               MultiArrayView<N, T2, S2> dest)
{
    for(int k=0; k<N; ++k)
        vigra_precondition(source.shape(k) == dest.shape(k) || source.shape(k) == 1 || 1 == dest.shape(k),
            "copyMultiArray(): shape mismatch between input and output.");
    if(source.shape() == dest.shape())
        copyMultiArray(srcMultiArrayRange(source), destMultiArray(dest));
    else
        copyMultiArray(srcMultiArrayRange(source), destMultiArrayRange(dest));
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
transformMultiArrayReduceImpl(SrcIterator s, SrcShape const &, SrcAccessor src,
               DestIterator d, DestShape const & dshape, DestAccessor dest, 
               SrcShape const & reduceShape,
               Functor const & ff, MetaInt<0>)
{
    DestIterator dend = d + dshape[0];
    for(; d < dend; ++s.template dim<0>(), ++d)
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
    for(; d < dend; ++s.template dim<N>(), ++d)
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
        for(; d < dend; ++d)
        {
            transformMultiArrayExpandImpl(s.begin(), sshape, src, d.begin(), dshape, dest,
                                          f, MetaInt<N-1>());
        }
    }
    else
    {
        for(; d < dend; ++s, ++d)
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

    Note: The effect of this function can often be achieved in a simpler and
    more readable way by means of \ref MultiMathModule "array expressions".
    
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
        
    The arrays must be represented by MultiArrayViews. If source and destination shapes
    match, standard mode is applied. If the shapes differ, the size of corresponding 
    dimensions must either be equal, or the source length must be 1 
    (expand mode), or the destination length must be 1 (reduce mode). However,
    reduction and expansion cannot be executed at the same time, so the latter
    conditions are mutual exclusive, even if they apply to different dimensions.
    
    <b> Declarations:</b>

    <b>\#include</b> \<vigra/multi_pointoperators.hxx\><br>
    Namespace: vigra
    
    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T1, class S1,
                                  class T2, class S2, 
                  class Functor>
        void
        transformMultiArray(MultiArrayView<N, T1, S1> const & source,
                            MultiArrayView<N, T2, S2> dest, Functor const & f);
    }
    \endcode
    
    \deprecatedAPI{transformMultiArray}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
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
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    \deprecatedEnd

    <b> Usage - Standard Mode:</b>

    Source and destination array have the same size.
    
    \code
    #include <cmath>         // for sqrt()

    MultiArray<3, float>  src(Shape3(100, 200, 50)),
                          dest(Shape3(100, 200, 50));
    ...
    
    transformMultiArray(src, dest, &std::sqrt );
    \endcode

    <b> Usage - Expand Mode:</b>

    The source array is effectively only a 2D image(it has a 3D shape, but depth is a singleton dimension 
    with length 1). Thus, the destination will contain 50 identical copies of the transformed source image. 
    
    \code
    #include <cmath>         // for sqrt()

    MultiArray<3, float> src(Shape3(100, 200, 1)),
                         dest(Shape3(100, 200, 50));
    ...
    
    transformMultiArray(src, dest, &std::sqrt );
    \endcode

    <b> Usage - Reduce Mode:</b>

    The destination array is effectively only 1D (it's width and height are singleton dimensions). 
    Thus, it will contain accumulated data for every slice of the source volume
    (or for every frame, if the source is interpreted as an image sequence).
    In the example, we use the functor \ref vigra::FindAverage to calculate
    the average gray value of every slice. 
    
    \code
    MultiArray<3, float>  src(Shape3(100, 200, 50)),
                          dest(Shape3(1, 1, 50));
    ...
    
    transformMultiArray(src, dest,
                        FindAverage<float>() );
    \endcode
    
    Note that the functor must define the appropriate traits described below in order to be 
    recognized as a reduce functor. This is most easily achieved by deriving from 
    <tt>UnaryReduceFunctorTag</tt> (see \ref vigra::FunctorTraits).

    \deprecatedUsage{transformMultiArray}
    \code
    #include <cmath>         // for sqrt()

    typedef vigra::MultiArray<3, float> Array;
    Array src(Shape3(100, 200, 50)),
          dest(Shape3(100, 200, 50));
    ...
    
    vigra::transformMultiArray(srcMultiArrayRange(src),
                               destMultiArray(dest),
                               (float(*)(float))&std::sqrt );

    \endcode
    \deprecatedEnd

    <b> Required Interface:</b>

    In standard and expand mode, the functor must be a model of UnaryFunction
    (i.e. support one-argument function call which accepts values of type
    <tt>T1</tt> and a return value that is convertible into <tt>T2</tt>.
    
    In reduce mode, it must be a model of UnaryAnalyser (i.e. support function call
    with one argument and no return value <tt>functor(arg)</tt>) and Initializer
    (i.e. support function call with no argument, but return value 
    <tt>res = functor()</tt>). Internally, such functors are recognized by the 
    meta functions <tt>FunctorTraits<FUNCTOR>::isUnaryAnalyser</tt> and
    <tt>FunctorTraits<FUNCTOR>::isInitializer</tt> which must both yield 
    <tt>VigraTrueType</tt>. Make sure that your functor correctly defines 
    <tt>FunctorTraits</tt> because otherwise reduce mode will not work. 
    This is most easily achieved by deriving the functor from 
    <tt>UnaryReduceFunctorTag</tt> (see \ref vigra::FunctorTraits).
    In addition, the functor must be copy constructible in order to start each reduction
    with a fresh functor.
    
    \see TransformFunctor, MultiMathModule, \ref FunctorExpressions
*/
doxygen_overloaded_function(template <...> void transformMultiArray)

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
inline void
transformMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               triple<DestIterator, DestShape, DestAccessor> const & dest, 
               Functor const & f)
{
    transformMultiArray(src.first, src.second, src.third, 
                        dest.first, dest.second, dest.third, f);
}

template <unsigned int N, class T1, class S1,
          class T2, class S2, 
          class Functor>
inline void
transformMultiArrayImpl(MultiArrayView<N, T1, S1> const & source,
                        MultiArrayView<N, T2, S2> dest,
                        Functor const & f, VigraFalseType)
{
    if(source.shape() == dest.shape())
        transformMultiArray(srcMultiArrayRange(source), destMultiArray(dest), f);
    else
        transformMultiArray(srcMultiArrayRange(source), destMultiArrayRange(dest), f);
}

template <unsigned int N, class T1, class S1,
          class T2, class S2, 
          class Functor>
inline void
transformMultiArrayImpl(MultiArrayView<N, T1, S1> const & source,
                        MultiArrayView<N, T2, S2> dest,
                        Functor const & f, VigraTrueType)
{
    transformMultiArray(srcMultiArrayRange(source), destMultiArrayRange(dest), f);
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2, 
          class Functor>
inline void
transformMultiArray(MultiArrayView<N, T1, S1> const & source,
                    MultiArrayView<N, T2, S2> dest, Functor const & f)
{
    for(unsigned int k=0; k<N; ++k)
        vigra_precondition(source.shape(k) == dest.shape(k) || source.shape(k) == 1 || 1 == dest.shape(k),
            "transformMultiArray(): shape mismatch between input and output.");

    typedef FunctorTraits<Functor> FT;
    typedef typename 
        And<typename FT::isInitializer, typename FT::isUnaryAnalyser>::result
        isAnalyserInitializer;
    transformMultiArrayImpl(source, dest, f, isAnalyserInitializer());
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
               SrcIterator1 s1, SrcShape const & , SrcAccessor1 src1,
               SrcIterator2 s2, SrcAccessor2 src2,
               DestIterator d,  DestShape const & dshape, DestAccessor dest, 
               SrcShape const & reduceShape,
               Functor const & ff, MetaInt<0>)
{
    DestIterator dend = d + dshape[0];
    for(; d < dend; ++s1.template dim<0>(), ++s2.template dim<0>(), ++d)
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
    for(; d < dend; ++s1.template dim<N>(), ++s2.template dim<N>(), ++d)
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
        for(; d < dend; ++d, ++s2)
            dest.set(f(sv1, src2(s2)), d);
    }
    else if(sshape2[0] == 1)
    {
        typename SrcAccessor2::value_type sv2 = src2(s2);
        for(; d < dend; ++d, ++s1)
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
    for(; d < dend; ++d, s1 += s1inc, s2 += s2inc)
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

    Note: The effect of this function can often be achieved in a simpler and
    more readable way by means of \ref MultiMathModule "array expressions".
    
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
        
    The arrays must be represented by MultiArrayViews. If all shapes are identical, 
    standard mode is applied. If the shapes differ, the size of corresponding dimensions
    must either be equal, or the length of this dimension must be 1 in one or both source 
    arrays (expand mode), or the destination length must be 1 (reduce mode). However,
    reduction and expansion cannot be executed at the same time, so the latter
    conditions are mutual exclusive, even if they apply to different dimensions.
    
    <b> Declarations:</b>
    
    <b>\#include</b> \<vigra/multi_pointoperators.hxx\><br>
    Namespace: vigra
    
    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T11, class S11,
                                  class T12, class S12,
                                  class T2, class S2, 
                  class Functor>
        void
        combineTwoMultiArrays(MultiArrayView<N, T11, S11> const & source1,
                              MultiArrayView<N, T12, S12> const & source2,
                              MultiArrayView<N, T2, S2> dest, 
                              Functor const & f);
    }
    \endcode
    
    \deprecatedAPI{combineTwoMultiArrays}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
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
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    \deprecatedEnd
    
    <b> Usage - Standard Mode:</b>
    
    Source and destination arrays have the same size.
    
    \code
    #include <functional>     // for std::plus

    MultiArray<3, int>  src1(Shape3(100, 200, 50)),
                        src2(Shape3(100, 200, 50)),
                        dest(Shape3(100, 200, 50));
    ...
    
    combineTwoMultiArrays(src1, src2, dest,  
                          std::plus<int>());
    \endcode
    
    <b> Usage - Expand Mode:</b>

    One source array is effectively only a 2D image (it has depth 1). This image will be added
    to every slice of the other source array, and the result is written into the 
    corresponding destination slice. 
    
    \code
    #include <functional>     // for std::plus

    MultiArray<3, int> src1(Shape3(100, 200, 1)),
                       src2(Shape3(100, 200, 50)),
                       dest(Shape3(100, 200, 50));
    ...
    
    combineTwoMultiArrays(src1, src2, dest,  
                          std::plus<int>());
    \endcode

    <b> Usage - Reduce Mode:</b>

    The destination array is only 1D (it's width and height are singleton dimensions). 
    Thus, it will contain accumulated data for every slice of the source volumes
    (or for every frame, if the sources are interpreted as image sequences).
    In the example, we use \ref vigra::ReduceFunctor together with a functor 
    expression (see \ref FunctorExpressions) to calculate the total absolute difference 
    of the gray values in every pair of source slices.
    
    \code
    #include <vigra/functorexpression.hxx>
    using namespace vigra::functor;
        
    MultiArray<3, int> src1(Shape3(100, 200, 50)),
                       src2(Shape3(100, 200, 50)),
                       dest(Shape3(1, 1, 50));
    ...
    
    combineTwoMultiArrays(src1, src2, dest,  
                          reduceFunctor(Arg1() + abs(Arg2() - Arg3()), 0) );
                          // Arg1() is the sum accumulated so far, initialized with 0
    \endcode

    Note that the functor must define the appropriate traits described below in order to be 
    recognized as a reduce functor. This is most easily achieved by deriving from 
    <tt>BinaryReduceFunctorTag</tt> (see \ref vigra::FunctorTraits).

    \deprecatedUsage{combineTwoMultiArrays}
    \code
    #include <functional>     // for std::plus

    typedef vigra::MultiArray<3, int> Array;
    Array src1(Shape3(100, 200, 50)),
          src2(Shape3(100, 200, 50)),
          dest(Shape3(100, 200, 50));
    ...
    
    vigra::combineTwoMultiArrays(
                srcMultiArrayRange(src1), 
                srcMultiArray(src2), 
                destMultiArray(dest),  
                std::plus<int>());
    \endcode
    \deprecatedEnd
    
    <b> Required Interface:</b>
    
    In standard and expand mode, the functor must be a model of BinaryFunction
    (i.e. support function call with two arguments and a return value which is convertible 
    into <tt>T2</tt>:  <tt>T2 res = functor(arg1, arg2)</tt>):
    
    In reduce mode, it must be a model of BinaryAnalyser (i.e. support function call
    with two arguments and no return value <tt>functor(arg1, arg2)</tt>) and Initializer
    (i.e. support function call with no argument, but return value 
    <tt>res = functor()</tt>). Internally, such functors are recognized by the 
    meta functions <tt>FunctorTraits<FUNCTOR>::isBinaryAnalyser</tt> and
    <tt>FunctorTraits<FUNCTOR>::isInitializer</tt> which must both yield 
    <tt>VigraTrueType</tt>. Make sure that your functor correctly defines 
    <tt>FunctorTraits</tt> because otherwise reduce mode will not work. 
    This is most easily achieved by deriving the functor from 
    <tt>BinaryReduceFunctorTag</tt> (see \ref vigra::FunctorTraits).
    In addition, the functor must be copy constructible in order to start each reduction
    with a fresh functor.
    
    \see TransformFunctor, MultiMathModule, \ref FunctorExpressions
*/
doxygen_overloaded_function(template <...> void combineTwoMultiArrays)

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

template <class SrcIterator1, class SrcShape, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestAccessor, class Functor>
inline void
combineTwoMultiArrays(triple<SrcIterator1, SrcShape, SrcAccessor1> const & src1,
                      pair<SrcIterator2, SrcAccessor2> const & src2,
                      pair<DestIterator, DestAccessor> const & dest, 
                      Functor const & f)
{
    
    combineTwoMultiArrays(
           src1.first, src1.second, src1.third, 
           src2.first, src2.second, dest.first, dest.second, f);
}

template <class SrcIterator1, class SrcShape1, class SrcAccessor1,
          class SrcIterator2, class SrcShape2, class SrcAccessor2,
          class DestIterator, class DestShape, class DestAccessor, 
          class Functor>
inline void
combineTwoMultiArrays(triple<SrcIterator1, SrcShape1, SrcAccessor1> const & src1,
                      triple<SrcIterator2, SrcShape2, SrcAccessor2> const & src2,
                      triple<DestIterator, DestShape, DestAccessor> const & dest, 
                      Functor const & f)
{
    combineTwoMultiArrays(src1.first, src1.second, src1.third, 
                          src2.first, src2.second, src2.third, 
                          dest.first, dest.second, dest.third, f);
}

template <unsigned int N, class T11, class S11,
                          class T12, class S12,
                          class T2, class S2, 
          class Functor>
inline void
combineTwoMultiArraysImpl(MultiArrayView<N, T11, S11> const & source1,
                          MultiArrayView<N, T12, S12> const & source2,
                          MultiArrayView<N, T2, S2> dest, 
                          Functor const & f, VigraFalseType)
{
    
    if(source1.shape() == source2.shape() && source1.shape() == dest.shape())
        combineTwoMultiArrays(srcMultiArrayRange(source1), 
                              srcMultiArray(source2), destMultiArray(dest), f);
    else
        combineTwoMultiArrays(srcMultiArrayRange(source1), 
                              srcMultiArrayRange(source2), 
                              destMultiArrayRange(dest), f);
}

template <unsigned int N, class T11, class S11,
                          class T12, class S12,
                          class T2, class S2, 
          class Functor>
inline void
combineTwoMultiArraysImpl(MultiArrayView<N, T11, S11> const & source1,
                          MultiArrayView<N, T12, S12> const & source2,
                          MultiArrayView<N, T2, S2> dest, 
                          Functor const & f, VigraTrueType)
{
    
    combineTwoMultiArrays(srcMultiArrayRange(source1), 
                          srcMultiArrayRange(source2), 
                          destMultiArrayRange(dest), f);
}

template <unsigned int N, class T11, class S11,
                          class T12, class S12,
                          class T2, class S2, 
          class Functor>
inline void
combineTwoMultiArrays(MultiArrayView<N, T11, S11> const & source1,
                      MultiArrayView<N, T12, S12> const & source2,
                      MultiArrayView<N, T2, S2> dest, 
                      Functor const & f)
{
    for(unsigned int k=0; k<N; ++k)
        vigra_precondition((source1.shape(k) == source2.shape(k) || source1.shape(k) == 1 || 1 == source2.shape(k)) &&
                           (source1.shape(k) == dest.shape(k) || source1.shape(k) == 1 || 1 == dest.shape(k)),
            "combineTwoMultiArrays(): shape mismatch between inputs and/or output.");

    typedef FunctorTraits<Functor> FT;
    typedef typename 
        And<typename FT::isInitializer, typename FT::isBinaryAnalyser>::result
        isAnalyserInitializer;
    combineTwoMultiArraysImpl(source1, source2, dest, f, isAnalyserInitializer());
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
    for(; s1 < s1end; ++s1, ++s2, ++s3, ++d)
    {
        combineThreeMultiArraysImpl(s1.begin(), shape, src1, 
                                  s2.begin(), src2, s3.begin(), src3, d.begin(), dest, 
                                  f, MetaInt<N-1>());
    }
}
    
    
/** \brief Combine three multi-dimensional arrays into one using a 
           ternary function or functor.

    Note: The effect of this function can often be achieved in a simpler and
    more readable way by means of \ref MultiMathModule "array expressions".
    
    Except for the fact that it operates on three input arrays, this function is
    identical to the standard mode of \ref combineTwoMultiArrays() (reduce and expand 
    modes are not supported).
    
    <b> Declarations:</b>
    
    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T11, class S11,
                                  class T12, class S12, 
                                  class T13, class S13,
                                  class T2, class S2, 
                  class Functor>
        void
        combineThreeMultiArrays(MultiArrayView<N, T11, S11> const & source1,
                                MultiArrayView<N, T12, S12> const & source2,
                                MultiArrayView<N, T13, S13> const & source3,
                                MultiArrayView<N, T2, S2> dest,
                                Functor const & f);
    }
    \endcode
    
    \deprecatedAPI{combineThreeMultiArrays}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
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
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/multi_pointoperators.hxx\><br>
    Namespace: vigra
    
    \code
    #include <vigra/functorexpression.hxx>
    
    MultiArray<3, int> src1(Shape3(100, 200, 50)),
                       src2(Shape3(100, 200, 50)),
                       src3(Shape3(100, 200, 50)),
                       dest(Shape3(100, 200, 50));
    ...
    
    using namespace vigra::functor; // activate VIGRA's lambda library
    
    combineThreeMultiArrays(src1, src2, src3, dest,  
                            Arg1()*exp(-abs(Arg2()-Arg3())));
    \endcode
    
    \deprecatedUsage{combineThreeMultiArrays}
    \code
    #include <functional>     // for plus

    typedef vigra::MultiArray<3, int> Array;
    Array src1(Shape3(100, 200, 50)),
          src2(Shape3(100, 200, 50)),
          src3(Shape3(100, 200, 50)),
          dest(Shape3(100, 200, 50));
    ...
    
    vigra::combineThreeMultiArrays(
                srcMultiArrayRange(src1), 
                srcMultiArray(src2), 
                srcMultiArray(src3), 
                destMultiArray(dest),  
                SomeThreeArgumentFunctor());
    \endcode
    \deprecatedEnd
    
    \see TransformFunctor, MultiMathModule, \ref FunctorExpressions
*/
doxygen_overloaded_function(template <...> void combineThreeMultiArrays)

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

template <unsigned int N, class T11, class S11,
                          class T12, class S12, 
                          class T13, class S13,
                          class T2, class S2, 
          class Functor>
inline void
combineThreeMultiArrays(MultiArrayView<N, T11, S11> const & source1,
                        MultiArrayView<N, T12, S12> const & source2,
                        MultiArrayView<N, T13, S13> const & source3,
                        MultiArrayView<N, T2, S2> dest, Functor const & f)
{
    vigra_precondition(source1.shape() == source2.shape() && source1.shape() == source3.shape() && source1.shape() == dest.shape(),
        "combineThreeMultiArrays(): shape mismatch between inputs and/or output.");
    
    combineThreeMultiArrays(
           srcMultiArrayRange(source1), 
           srcMultiArray(source2), srcMultiArray(source3), destMultiArray(dest), f);
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
    for(; s < send; ++s)
    {
        inspectMultiArrayImpl(s.begin(), shape, a, f, MetaInt<N-1>());
    }
}
    
/** \brief Call an analyzing functor at every element of a multi-dimensional array.

    This function can be used to collect statistics of the array etc.
    The results must be stored in the functor, which serves as a return
    value (therefore, it is passed to the function by reference). The array must be 
    represented as a MultiArrayView.
    
    For many common statistics, the use of \ref  vigra::acc::extractFeatures() in combination with 
    \ref FeatureAccumulators is more convenient.

    <b> Declarations:</b>

    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T, class S, 
                  class Functor>
        void
        inspectMultiArray(MultiArrayView<N, T, S> const & s, 
                          Functor & f);
    }
    \endcode

    \deprecatedAPI{inspectMultiArray}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
    \code
    namespace vigra {
        template <class Iterator, class Shape, class Accessor, class Functor>
        void
        inspectMultiArray(Iterator s, Shape const & shape, Accessor a,  Functor & f);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class Iterator, class Shape, class Accessor, class Functor>
        void
        inspectMultiArray(triple<Iterator, Shape, Accessor> const & s, Functor & f);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_pointoperators.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<3, int>  array(Shape3(100, 200, 50));
    ... // fill array
    
    // init functor
    FindMinMax<int> minmax;

    inspectMultiArray(array, minmax);

    cout << "Min: " << minmax.min << " Max: " << minmax.max;
    \endcode
    The functor must support function call with one argument.

    \deprecatedUsage{inspectMultiArray}
    \code
    typedef vigra::MultiArray<3, int> Array;
    Array array(Shape3(100, 200, 50));

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
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void inspectMultiArray)

template <class Iterator, class Shape, class Accessor>
struct inspectMultiArray_binder
{
    Iterator      s;
    const Shape & shape;
    Accessor      a;
    inspectMultiArray_binder(Iterator s_, const Shape & shape_, Accessor a_)
        : s(s_), shape(shape_), a(a_) {}
    template <class Functor>
    void operator()(Functor & f)
    {
        inspectMultiArrayImpl(s, shape, a, f, MetaInt<Iterator::level>());
    }
};

template <class Iterator, class Shape, class Accessor, class Functor>
inline void
inspectMultiArray(Iterator s, Shape const & shape, Accessor a, Functor & f)
{
    inspectMultiArray_binder<Iterator, Shape, Accessor> g(s, shape, a);
    detail::extra_passes_select(g, f);
}
    
template <class Iterator, class Shape, class Accessor, class Functor>
inline void
inspectMultiArray(triple<Iterator, Shape, Accessor> const & s, Functor & f)
{
    inspectMultiArray(s.first, s.second, s.third, f);
}
    
template <unsigned int N, class T, class S, class Functor>
inline void
inspectMultiArray(MultiArrayView<N, T, S> const & s, Functor & f)
{
    inspectMultiArray(srcMultiArrayRange(s), f);
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
    for(; s1 < s1end; ++s1, ++s2)
    {
        inspectTwoMultiArraysImpl(s1.begin(), shape, a1, 
                                  s2.begin(), a2, f, MetaInt<N-1>());
    }
}
    
/** \brief Call an analyzing functor at all corresponding elements of 
           two multi-dimensional arrays.

    This function can be used to collect statistics over tow arrays.
    For example, one can holde data, and the other region labels or weights.
    The results must be stored in the functor, which serves as a return
    value (and is therefore passed by reference). The arrays must be represented by
    MultiArrayViews.
    
    For many common statistics, the use of \ref  vigra::acc::extractFeatures() in combination with 
    \ref FeatureAccumulators is more convenient.

    <b> Declarations:</b>

    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <unsigned int N, class T1, class S1, 
                                  class T2, class S2, 
                  class Functor>
        void
        inspectTwoMultiArrays(MultiArrayView<N, T1, S1> const & s1, 
                              MultiArrayView<N, T2, S2> const & s2,
                              Functor & f);
    }
    \endcode

    \deprecatedAPI{inspectTwoMultiArrays}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
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
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_pointoperators.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<3, int>  array1(Shape3(100, 200, 50)),
                        array2(Shape3(100, 200, 50));

    // init functor
    SomeStatisticsFunctor stats(..);

    inspectTwoMultiArrays(array1, array2, stats);
    \endcode
    The functor must support function call with two arguments.

    \deprecatedUsage{inspectTwoMultiArrays}
    \code
    MultiArray<3, int>  array1(Shape3(100, 200, 50)),
                        array2(Shape3(100, 200, 50));

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
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void inspectTwoMultiArrays)

template <class Iterator1, class Shape, class Accessor1,
          class Iterator2, class Accessor2>
struct inspectTwoMultiArrays_binder
{
    Iterator1     s1;
    const Shape & shape;
    Accessor1     a1;
    Iterator2     s2;
    Accessor2     a2;
    inspectTwoMultiArrays_binder(Iterator1 s1_, const Shape & shape_,
                                 Accessor1 a1_, Iterator2 s2_, Accessor2 a2_)
        : s1(s1_), shape(shape_), a1(a1_), s2(s2_), a2(a2_) {}
    template <class Functor>
    void operator()(Functor & f)
    {
        inspectTwoMultiArraysImpl(s1, shape, a1, s2, a2, f,
                                  MetaInt<Iterator1::level>());
    }
};
    
template <class Iterator1, class Shape, class Accessor1,
          class Iterator2, class Accessor2,
          class Functor>
inline void
inspectTwoMultiArrays(Iterator1 s1, Shape const & shape, Accessor1 a1,
                      Iterator2 s2, Accessor2 a2, Functor & f)
{
    inspectTwoMultiArrays_binder<Iterator1, Shape, Accessor1,
                                 Iterator2, Accessor2>
        g(s1, shape, a1, s2, a2);
    detail::extra_passes_select(g, f);
}
    
template <class Iterator1, class Shape, class Accessor1, 
          class Iterator2, class Accessor2, 
          class Functor>
inline void
inspectTwoMultiArrays(triple<Iterator1, Shape, Accessor1> const & s1, 
                      pair<Iterator2, Accessor2> const & s2, Functor & f)
{
    inspectTwoMultiArrays(s1.first, s1.second, s1.third, 
                          s2.first, s2.second, f);
}
    
template <unsigned int N, class T1, class S1, 
                          class T2, class S2, 
          class Functor>
inline void
inspectTwoMultiArrays(MultiArrayView<N, T1, S1> const & s1, 
                      MultiArrayView<N, T2, S2> const & s2, Functor & f)
{
    vigra_precondition(s1.shape() == s2.shape(),
        "inspectTwoMultiArrays(): shape mismatch between inputs.");
    
    inspectTwoMultiArrays(srcMultiArrayRange(s1), 
                          srcMultiArray(s2), f);
}
    
//@}

}  //-- namespace vigra


#endif  //-- VIGRA_MULTI_POINTOPERATORS_H
