/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    ( Version 1.1.4a, Nov 20 2001 )                                   */
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
/** \brief Define the default accessor for each image iterator.

    With each image iterator, a default accessor is associated. The type of
    this accessor can be retrieved by means of IteratorTraits. This is, 
    for example, used by the \ref IteratorBasedArgumentObjectFactories.
    
    \code
    template <class Iterator>
    void foo(Iterator i)
    {
        typedef typename IteratorTraits<Iterator>::DefaultAccessor Accessor;
        Accessor a;
        ...
    }
    \endcode
    
    <b>\#include</b> "<a href="iteratortraits_8hxx-source.html">vigra/iteratortraits.hxx</a>"
    
    Namespace: vigra
*/
template <class T> 
struct IteratorTraits 
{
    struct IteratorTraitsNotDefinedForThisCase
    {};
    
    // must be defined due to Visual C++ bug
    typedef IteratorTraitsNotDefinedForThisCase DefaultAccessor;
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
