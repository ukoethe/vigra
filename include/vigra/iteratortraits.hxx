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
 
 
#ifndef VIGRA_ITERATORTRAITS_HXX
#define VIGRA_ITERATORTRAITS_HXX

#include <vigra/accessor.hxx>
#include <vigra/imageiteratoradapter.hxx>

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
        typedef StandardAccessor<value_type>          DefaultAccessor;
        typedef StandardAccessor<value_type>          default_accessor;
    };
    \endcode
    
    By (partially) specializing this template for an iterator class
    the defaults given above can be changed as approiate. For example, iterators
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
    
    <b>\#include</b> "<a href="iteratortraits_8hxx-source.html">vigra/iteratortraits.hxx</a>"
    Namespace: vigra
*/
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
    typedef StandardAccessor<value_type>          DefaultAccessor;
    typedef StandardAccessor<value_type>          default_accessor;
};

//@}

template <> 
struct IteratorTraits<Diff2D > 
{
    typedef Diff2D                               Iterator;
    typedef Iterator                             iterator;
    typedef iterator::iterator_category          iterator_category;
    typedef iterator::value_type                 value_type;
    typedef iterator::reference                  reference;
    typedef iterator::index_reference            index_reference;
    typedef iterator::pointer                    pointer;
    typedef iterator::difference_type            difference_type;
    typedef iterator::row_iterator               row_iterator;
    typedef iterator::column_iterator            column_iterator;
    typedef StandardConstValueAccessor<Diff2D>   DefaultAccessor;
    typedef StandardConstValueAccessor<Diff2D>   default_accessor;
    
};

} // namespace vigra

#endif // VIGRA_ITERATORTRAITS_HXX
