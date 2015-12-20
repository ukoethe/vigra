/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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


#ifndef VIGRA_ITERATORTRAITS_HXX
#define VIGRA_ITERATORTRAITS_HXX

#include "accessor.hxx"
#include "imageiteratoradapter.hxx"

namespace vigra {

/** \addtogroup ImageIterators
*/
//@{
/** \brief Export associated information for each image iterator.

    The IteratorTraits class contains the following fields:

    \code
    template <class T>
    struct IteratorTraits
    {
        typedef T                                     Iterator;
        typedef Iterator                              iterator;
        typedef typename iterator::iterator_category  iterator_category;
        typedef typename iterator::value_type         value_type;
        typedef typename iterator::reference          reference;
        typedef typename iterator::index_reference    index_reference;
        typedef typename iterator::pointer            pointer;
        typedef typename iterator::difference_type    difference_type;
        typedef typename iterator::row_iterator       row_iterator;
        typedef typename iterator::column_iterator    column_iterator;
        typedef typename
         AccessorTraits<value_type>::default_accessor DefaultAccessor;
        typedef DefaultAccessor                       default_accessor;

        typedef VigraTrueType/VigraFalseType          hasConstantStrides;
    };
    \endcode

    By (partially) specializing this template for an iterator class
    the defaults given above can be changed as appropriate. For example, iterators
    for rgb images are associated with <TT>RGBAccessor<value_type></TT>
    instead of <TT>StandardAccessor<value_type></TT>. To get the accessor
    associated with a given iterator, use code like this:

    \code
    template <class Iterator>
    void foo(Iterator i)
    {
        typedef typename IteratorTraits<Iterator>::DefaultAccessor Accessor;
        Accessor a;
        ...
    }
    \endcode

    This technique is, for example, used by the
    \ref IteratorBasedArgumentObjectFactories. The possibility to retrieve the default accessor by means of a traits
    class is especially important since this information is not
    contained in the iterator directly.
    
    The member <tt>hasConstantStrides</tt> is useful for certain 
    optimizations: it helps to decide whether we can replace iterator
    operations such as <tt>iter++</tt> or <tt>iter += n</tt> with
    corresponding pointer operations (which may be faster), where
    the pointer is obtained as the address of iterator's pointee 
    (the object the iterator currently  refers to). 
    This flag would be <tt>VigraFalseType</tt> for a
    <tt>std::list<int>::iterator</tt>, but is <tt>VigraTrueType</tt> 
    for most VIGRA iterators.

    <b>\#include</b> \<vigra/iteratortraits.hxx\>
    Namespace: vigra
*/
template <class T>
struct IteratorTraits
{
    typedef T                                          Iterator;
    typedef Iterator                                   iterator;
    typedef typename iterator::iterator_category       iterator_category;
    typedef typename iterator::value_type              value_type;
    typedef typename iterator::reference               reference;
    typedef typename iterator::index_reference         index_reference;
    typedef typename iterator::pointer                 pointer;
    typedef typename iterator::difference_type         difference_type;
    typedef typename iterator::row_iterator            row_iterator;
    typedef typename iterator::column_iterator         column_iterator;
    typedef typename
        AccessorTraits<value_type>::default_accessor   DefaultAccessor;
    typedef DefaultAccessor                            default_accessor;

    // default: disable the constant strides optimization
    typedef VigraFalseType                             hasConstantStrides;
};

template <class T>
struct IteratorTraitsBase
{
    typedef T                                     Iterator;
    typedef Iterator                              iterator;
    typedef typename iterator::iterator_category  iterator_category;
    typedef typename iterator::value_type         value_type;
    typedef typename iterator::reference          reference;
    typedef typename iterator::index_reference    index_reference;
    typedef typename iterator::pointer            pointer;
    typedef typename iterator::difference_type    difference_type;
    typedef typename iterator::row_iterator       row_iterator;
    typedef typename iterator::column_iterator    column_iterator;
};


//@}


/***********************************************************/

/** \page ArgumentObjectFactories Argument Object Factories

    Factory functions to create argument objects which simplify long argument lists.

    <UL style="list-style-image:url(documents/bullet.gif)">
    <LI> \ref ImageBasedArgumentObjectFactories
    <LI> \ref MultiArrayBasedArgumentObjectFactories
    <LI> \ref IteratorBasedArgumentObjectFactories
    </UL>

    Long argument lists provide for greater flexibility of functions,
    but they are also tedious and error prone, when we don't need
    the flexibility. Thus, we define argument objects which
    automatically provide reasonable defaults for those arguments that we
    didn't specify explicitly.

    The argument objects are created via a number of factory functions.
    Since these functions have descriptive names, they also serve
    to improve readability: the name of each factory tells the purpose of its
    argument object.

    Consider the following example. Without argument objects we had to
    write something like this (cf. \ref copyImageIf()):

    \code
    vigra::BImage img1, img2, img3;

    // fill img1 and img2 ...

    vigra::copyImageIf(img1.upperLeft(), img1.lowerRight(), img1.accessor(),
                img2.upperLeft(), img2.accessor(),
                img3.upperLeft(), img3.accessor());
    \endcode

    Using the argument object factories, this becomes much shorter and
    more readable:

    \code
    vigra::copyImageIf(srcImageRange(img1),
                maskImage(img2),
                destImage(img3));
    \endcode

    The names of the factories clearly tell which image is source, mask,
    and destination. In addition, the suffix <TT>Range</TT> must be used
    for those argument objects that need to specify the lower right
    corner of the region of interest. Typically, this is only the first
    source argument, but sometimes the first destiniation argument must
    also contain a range.

    The factory functions come in two flavours: Iterator based and
    image based factories. Above we have seen the image based variant.
    The iterator based variant would look like this:

    \code
    vigra::copyImageIf(srcIterRange(img1.upperLeft(), img1.lowerRight()),
                maskIter(img2.upperLeft()),
                destIter(img3.upperLeft()));
    \endcode

    These factory functions contain the word <TT>Iter</TT> instead of the word
    <TT>Image</TT>,  They would normally be used if we couldn't access the
    images (for example, within a function which got passed iterators)
    or if we didn't want to operate on the entire image. The default
    accessor is obtained via \ref vigra::IteratorTraits.

    All factory functions also allow to specify accessors explicitly. This
    is useful if we can't use the default accessor. This variant looks
    like this:

    \code
    vigra::copyImageIf(srcImageRange(img1),
                maskImage(img2, MaskPredicateAccessor()),
                destImage(img3));
    \endcode

    or

    \code
    vigra::copyImageIf(srcIterRange(img1.upperLeft(), img1.lowerRight()),
                maskIter(img2.upperLeft(), MaskPredicateAccessor()),
                destIter(img3.upperLeft()));
    \endcode

    All versions can be mixed freely within one expression.
    Technically, the argument objects are simply defined as
    pairs and triples of iterators and accessor so that all algorithms
    should declare a call interface version based on pairs and triples
    (see for example \ref copyImageIf()).

  \section ImageBasedArgumentObjectFactories Image Based Argument Object Factories

    <b>Include:</b> automatically included with the image classes<br>
    Namespace: vigra

    These factories can be used to create argument objects when we
    are given instances or subclasses of \ref vigra::BasicImage (see
    \ref StandardImageTypes for instances defined per default).
    These factory functions access <TT>img.upperLeft()</TT>,
    <TT>img.lowerRight()</TT>, and <TT>img.accessor()</TT> to obtain the iterators
    and accessor for the given image (unless the accessor is
    given explicitly). The following factory functions are provided:

    <table>
    <tr><th bgcolor="#f0e0c0" colspan=2 align=left>
        <TT>\ref vigra::BasicImage "vigra::BasicImage<SomeType>" img;</TT> or <br>
         <TT>\ref vigra::BasicImageView "vigra::BasicImageView<SomeType>" img;</TT>
        </th>
    </tr>
    <tr><td>

    <TT>srcImageRange(img)</TT>
    </td><td>
        create argument object containing upper left, lower right, and
        default accessor of source image

    </td></tr>
    <tr><td>

    <TT>srcImageRange(img, Rect2D(...))</TT>
    </td><td>
        create argument object containing the ROI specified by <tt>\ref vigra::Rect2D</tt> and
        default accessor of source image

    </td></tr>
    <tr><td>

    <TT>srcImageRange(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing upper left, lower right
        of source image, and given accessor

    </td></tr>
    <tr><td>

    <TT>srcImageRange(img, Rect2D(...), SomeAccessor())</TT>
    </td><td>
        create argument object containing the ROI specified by <tt>\ref vigra::Rect2D</tt> and
        of source image, and given accessor

    </td></tr>
    <tr><td>

    <TT>srcImage(img)</TT>
    </td><td>
        create argument object containing upper left, and
        default accessor of source image

    </td></tr>
    <tr><td>

    <TT>srcImage(img, Point2D(...))</TT>
    </td><td>
        create argument object with upper left at point given by <tt>\ref vigra::Point2D</tt>, and
        default accessor of source image

    </td></tr>
    <tr><td>

    <TT>srcImage(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing upper left
        of source image, and given accessor

    </td></tr>
    <tr><td>

    <TT>srcImage(img, Point2D(...), SomeAccessor())</TT>
    </td><td>
        create argument object with upper left at point given by <tt>\ref vigra::Point2D</tt> of source image,
        and given accessor

    </td></tr>
    <tr><td>

    <TT>maskImage(img)</TT>
    </td><td>
        create argument object containing upper left, and
        default accessor of mask image

    </td></tr>
     <tr><td>

    <TT>maskImage(img, Point2D(...))</TT>
    </td><td>
        create argument object with upper left at point given by <tt>\ref vigra::Point2D</tt>, and
        default accessor of mask image

    </td></tr>
   <tr><td>

    <TT>maskImage(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing upper left
        of mask image, and given accessor

    </td></tr>
    <tr><td>

    <TT>maskImage(img, Point2D(...), SomeAccessor())</TT>
    </td><td>
        create argument object with upper left at point given by <tt>\ref vigra::Point2D</tt> of mask image,
        and given accessor

    </td></tr>
    <tr><td>

    <TT>destImageRange(img)</TT>
    </td><td>
        create argument object containing upper left, lower right, and
        default accessor of destination image

    </td></tr>
    <tr><td>

    <TT>destImageRange(img, Rect2D(...))</TT>
    </td><td>
        create argument object containing the ROI specified by <tt>\ref vigra::Rect2D</tt> and
        default accessor of destination image

    </td></tr>
    <tr><td>

    <TT>destImageRange(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing upper left, lower right
        of destination image, and given accessor

    </td></tr>
    <tr><td>

    <TT>destImageRange(img, Rect2D(...), SomeAccessor())</TT>
    </td><td>
        create argument object containing the ROI specified by <tt>\ref vigra::Rect2D</tt>
        of destination image, and given accessor

    </td></tr>
     <tr><td>

    <TT>destImage(img)</TT>
    </td><td>
        create argument object containing upper left, and
        default accessor of destination image

    </td></tr>
     <tr><td>

    <TT>destImage(img, Point2D(...))</TT>
    </td><td>
        create argument object with upper left at point given by <tt>\ref vigra::Point2D</tt>, and
        default accessor of destination image

    </td></tr>
    <tr><td>

    <TT>destImage(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing upper left
        of destination image, and given accessor

    </td></tr>
    <tr><td>

    <TT>destImage(img, Point2D(...), SomeAccessor())</TT>
    </td><td>
        create argument object with upper left at point given by <tt>\ref vigra::Point2D</tt> of destination image,
        and given accessor

    </td></tr>
    </table>


  \section MultiArrayBasedArgumentObjectFactories MultiArrayView Based Argument Object Factories

    <b>Include:</b> automatically included with 
       \<vigra/multi_array.hxx\><br>
    Namespace: vigra

    These factories can be used to create argument objects when we
    are given instances or subclasses of \ref vigra::MultiArrayView.
    These factory functions access <TT>array.traverser_begin()</TT>,
    <TT>array.traverser_end()</TT> to obtain the iterators. If no accessor is
    given, they use the <tt>AccessorTraits<T></tt> to determine the default 
    accessor associated with the array's value type <tt>T</tt>.
    The following factory functions are provided:

    <table>
    <tr><th bgcolor="#f0e0c0" colspan=2 align=left>
        <TT>\ref vigra::MultiArrayView "vigra::MultiArrayView<N, SomeType>" img;</TT>
        </th>
    </tr>
    <tr><td>

    <TT>srcMultiArrayRange(img)</TT>
    </td><td>
        create argument object containing a \ref vigra::MultiIterator 
        marking the begin of the array, a shape object giving the desired
        shape of the array (possibly a subarray) and the default const accessor for
        <tt>SomeType</tt>

    </td></tr>
    <tr><td>

    <TT>srcMultiArrayRange(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing a \ref vigra::MultiIterator 
        marking the begin of the array, a shape object giving the desired
        shape of the array (possibly a subarray) and the given accessor

    </td></tr>
    <tr><td>

    <TT>srcMultiArray(img)</TT>
    </td><td>
        create argument object containing a \ref vigra::MultiIterator
        marking the begin of the array, and the default const accessor for
        <tt>SomeType</tt>

    </td></tr>
    <tr><td>

    <TT>srcMultiArray(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing a \ref vigra::MultiIterator 
        marking the begin of the array and the given accessor

    </td></tr>
    <tr><td>

    <TT>destMultiArrayRange(img)</TT>
    </td><td>
        create argument object containing a \ref vigra::MultiIterator 
        marking the begin of the array, a shape object giving the desired
        shape of the array (possibly a subarray) and the default accessor for
        <tt>SomeType</tt>

    </td></tr>
    <tr><td>

    <TT>destMultiArrayRange(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing a \ref vigra::MultiIterator's 
        marking the begin of the array, a shape object giving the desired
        shape of the array (possibly a subarray) and the given accessor

    </td></tr>
    <tr><td>

    <TT>destMultiArray(img)</TT>
    </td><td>
        create argument object containing a \ref vigra::MultiIterator 
        marking the begin of the array and the default accessor for
        <tt>SomeType</tt>

    </td></tr>
    <tr><td>

    <TT>destMultiArray(img, SomeAccessor())</TT>
    </td><td>
        create argument object containing a \ref vigra::MultiIterator's 
        marking the begin of the array and the given accessor

    </td></tr>
    </table>


  \section IteratorBasedArgumentObjectFactories Iterator Based Argument Object Factories

    <b>\#include</b> \<vigra/iteratortraits.hxx\>
    Namespace: vigra

    These factories can be used to create argument objects when we
    are given \ref ImageIterators.
    These factory functions use \ref vigra::IteratorTraits to
    get the default accessor for the given iterator unless the
    accessor is given explicitly. The following factory functions
    are provided:

    <table>
    <tr><th bgcolor="#f0e0c0" colspan=2 align=left>
        <TT>\ref vigra::BasicImage::Iterator "vigra::BasicImage<SomeType>::Iterator" i1, i2;</TT>
        </th>
    </tr>
    <tr><td>

    <TT>srcIterRange(i1, i2)</TT>
    </td><td>
        create argument object containing the given iterators and
        corresponding default accessor (for source image)

    </td></tr>
    <tr><td>

    <TT>srcIterRange(i1, i2, SomeAccessor())</TT>
    </td><td>
        create argument object containing given iterators and
        accessor (for source image)

    </td></tr>
    <tr><td>

    <TT>srcIter(i1)</TT>
    </td><td>
        create argument object containing the given iterator and
        corresponding default accessor (for source image)

    </td></tr>
    <tr><td>

    <TT>srcIter(i1, SomeAccessor())</TT>
    </td><td>
        create argument object containing given iterator and
        accessor (for source image)

    </td></tr>
    <tr><td>

    <TT>maskIter(i1)</TT>
    </td><td>
        create argument object containing the given iterator and
        corresponding default accessor (for mask image)

    </td></tr>
    <tr><td>

    <TT>maskIter(i1, SomeAccessor())</TT>
    </td><td>
        create argument object containing given iterator and
        accessor (for mask image)

    </td></tr>
    <tr><td>

    <TT>destIterRange(i1, i2)</TT>
    </td><td>
        create argument object containing the given iterators and
        corresponding default accessor (for destination image)

    </td></tr>
    <tr><td>

    <TT>destIterRange(i1, i2, SomeAccessor())</TT>
    </td><td>
        create argument object containing given iterators and
        accessor (for destination image)

    </td></tr>
    <tr><td>

    <TT>destIter(i1)</TT>
    </td><td>
        create argument object containing the given iterator and
        corresponding default accessor (for destination image)

    </td></tr>
    <tr><td>

    <TT>destIter(i1, SomeAccessor())</TT>
    </td><td>
        create argument object containing given iterator and
        accessor (for destination image)

    </td></tr>
    </table>
*/

template <class Iterator, class Accessor>
inline triple<Iterator, Iterator, Accessor>
srcIterRange(Iterator const & upperleft, Iterator const & lowerright, Accessor a)
{
    return triple<Iterator, Iterator, Accessor>(upperleft, lowerright, a);
}

template <class Iterator, class Accessor>
inline pair<Iterator, Accessor>
srcIter(Iterator const & upperleft, Accessor a)
{
    return pair<Iterator, Accessor>(upperleft, a);
}

template <class Iterator, class Accessor>
inline pair<Iterator, Accessor>
maskIter(Iterator const & upperleft, Accessor a)
{
    return pair<Iterator, Accessor>(upperleft, a);
}

template <class Iterator, class Accessor>
inline pair<Iterator, Accessor>
destIter(Iterator const & upperleft, Accessor a)
{
    return pair<Iterator, Accessor>(upperleft, a);
}


template <class Iterator, class Accessor>
inline triple<Iterator, Iterator, Accessor>
destIterRange(Iterator const & upperleft, Iterator const & lowerright, Accessor a)
{
    return triple<Iterator, Iterator, Accessor>(upperleft, lowerright, a);
}

template <class Iterator>
inline pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>
srcIter(Iterator const & upperleft)
{
    return pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>(
                  upperleft,
                  typename IteratorTraits<Iterator>::DefaultAccessor());
}

template <class Iterator>
inline triple<Iterator, Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>
srcIterRange(Iterator const & upperleft, Iterator const & lowerright)
{
    return triple<Iterator, Iterator,
                  typename IteratorTraits<Iterator>::DefaultAccessor>(
                  upperleft, lowerright,
                  typename IteratorTraits<Iterator>::DefaultAccessor());
}

template <class Iterator>
inline pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>
maskIter(Iterator const & upperleft)
{
    return pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>(
                  upperleft,
                  typename IteratorTraits<Iterator>::DefaultAccessor());
}

template <class Iterator>
inline pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>
destIter(Iterator const & upperleft)
{
    return pair<Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>(
                  upperleft,
                  typename IteratorTraits<Iterator>::DefaultAccessor());
}

template <class Iterator>
inline triple<Iterator, Iterator, typename IteratorTraits<Iterator>::DefaultAccessor>
destIterRange(Iterator const & upperleft, Iterator const & lowerright)
{
    return triple<Iterator, Iterator,
                  typename IteratorTraits<Iterator>::DefaultAccessor>(
                  upperleft, lowerright,
                  typename IteratorTraits<Iterator>::DefaultAccessor());
}

} // namespace vigra

#endif // VIGRA_ITERATORTRAITS_HXX
