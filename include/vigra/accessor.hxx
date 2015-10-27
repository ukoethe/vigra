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

#ifndef VIGRA_ACCESSOR_HXX
#define VIGRA_ACCESSOR_HXX

#include "metaprogramming.hxx"
#include "numerictraits.hxx"
#include "tuple.hxx"

namespace vigra {

/** \addtogroup DataAccessors Data Accessors

Basic templates to encapsulate access to the data of an iterator.

Data accessors are used to allow for flexible access to the data
an iterator points to. When we access the data directly, we
are bound to what <TT>operator*()</TT> returns, if this method exists at
all. Encapsulating access in an accessor enables a better
decoupling of data structures and algorithms.
<a href="http://ukoethe.github.io/vigra/doc/vigra/documents/DataAccessors.ps">This paper</a> contains
a detailed description of the concept. Here is a brief list of the basic
accessor requirements:

<p>
<table border=2 cellspacing=0 cellpadding=2 width="100%">
<tr><th>
    Operation
    </th><th>
    Result
    </th><th>
    Semantics
    </th>
</tr>
<tr>
    <td><tt>accessor(iter)</tt></td><td>convertible to <br><tt>Accessor::value_type const &</tt></td>
    <td>read data at the current position of the iterator</td>
</tr>
<tr>
    <td><tt>accessor(iter, index)</tt></td><td>convertible to <br><tt>Accessor::value_type const &</tt></td>
    <td>read data at offset <tt>index</tt> relative to iterator's current position
    (random-access iterator only)</td>
</tr>
<tr>
    <td><tt>accessor.set(value, iter)</tt></td><td><tt>void</tt></td>
    <td>write data <tt>value</tt> at the current position of the iterator (mutable iterator only)</td>
</tr>
<tr>
    <td><tt>accessor.set(value, iter, index)</tt></td><td><tt>void</tt></td>
    <td>write data <tt>value</tt> at offset <tt>index</tt> relative to iterator's current position
    (mutable random-access iterator only)</td>
</tr>
<tr><td colspan=2>
    <tt>Accessor::value_type</tt></td>
    <td>type of the data field the accessor refers to</td>
</tr>
<tr><td colspan=3>
        <tt>iter</tt> is an iterator<br>
        <tt>index</tt> has the iterator's index type (<tt>Iterator::difference_type</tt>)<br>
        <tt>value</tt> is convertible to <tt>Accessor::value_type const &</tt>
    </td>
</tr>
</table>
</p>

The template <tt>AccessorTraits<T></tt> can be used to find the default accessor
associated with the type <tt>T</tt>, e.g.

\code
typedef typename AccessorTraits<typename Image::value_type>::default_accessor       Accessor;
typedef typename AccessorTraits<typename Image::value_type>::default_const_accessor ConstAccessor;
\endcode
*/
//@{

/********************************************************/
/*                                                      */
/*                     StandardAccessor                 */
/*                                                      */
/********************************************************/

/** \brief Encapsulate access to the values an iterator points to.

    StandardAccessor is a trivial accessor that simply encapsulates
    the iterator's operator*() and operator[]() in its
    read and write functions. It passes its arguments <em>by reference</em>.
    If you want to return items by value, you
    must use StandardValueAccessor instead of StandardAccessor.
    Both accessors have different optimization properties --
    StandardAccessor is usually faster for compound pixel types,
    while StandardValueAccessor is faster for the built-in types.

    When a floating point number is assigned by means of an accessor
    with integral value_type, the value is rounded and clipped as appropriate.

    <b>\#include</b> \<vigra/accessor.hxx\><br>
    Namespace: vigra
*/
template <class VALUETYPE>
class StandardAccessor
{
  public:
        /** the value_type
        */
    typedef VALUETYPE value_type;

        /** read the current data item
        */
    template <class ITERATOR>
    VALUETYPE const & operator()(ITERATOR const & i) const { return *i; }

    VALUETYPE const & operator()(VALUETYPE const * i) const { return *i; }

        /** read the data item at an offset (can be 1D or 2D or higher order difference).
        */
    template <class ITERATOR, class OFFSET>
    VALUETYPE const & operator()(ITERATOR const & i, OFFSET const & diff) const
    {
        return i[diff];
    }

        /** Write the current data item. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR>
    void set(V const & value, ITERATOR const & i) const
    { *i = detail::RequiresExplicitCast<VALUETYPE>::cast(value); }

        /* This overload is needed to make the accessor work with a std::back_inserter */
    template <class V, class ITERATOR>
    void set(V const & value, ITERATOR & i) const
    { *i = detail::RequiresExplicitCast<VALUETYPE>::cast(value); }

        /** Write the data item at an offset (can be 1D or 2D or higher order difference)..
            The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR, class OFFSET>
    void set(V const & value, ITERATOR const & i, OFFSET const & diff) const
    {
        i[diff]= detail::RequiresExplicitCast<VALUETYPE>::cast(value);
    }
};

/** \brief Encapsulate access to the values an iterator points to.

    StandardValueAccessor is a trivial accessor that simply encapsulates
    the iterator's operator*() and operator[]() in its
    read and write functions. It passes its arguments <em>by value</em>.
    If the iterator returns its items by reference (such as \ref vigra::ImageIterator),
    you can also use StandardAccessor.
    These accessors have different optimization properties --
    StandardAccessor is usually faster for compound pixel types,
    while StandardValueAccessor is faster for the built-in types.

    When a floating point number is assigned by means of an accessor
    with integral value_type, the value is rounded and clipped as appropriate.

    <b>\#include</b> \<vigra/accessor.hxx\><br>
    Namespace: vigra
*/
template <class VALUETYPE>
class StandardValueAccessor
{
  public:
        /** the value_type
        */
    typedef VALUETYPE value_type;

        /** Read the current data item. The type <TT>ITERATOR::reference</TT>
            is automatically converted to <TT>VALUETYPE</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class ITERATOR>
    VALUETYPE operator()(ITERATOR const & i) const
        { return detail::RequiresExplicitCast<VALUETYPE>::cast(*i); }

        /** Read the data item at an offset (can be 1D or 2D or higher order difference).
            The type <TT>ITERATOR::index_reference</TT>
            is automatically converted to <TT>VALUETYPE</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class ITERATOR, class OFFSET>
    VALUETYPE operator()(ITERATOR const & i, OFFSET const & diff) const
    {
        return detail::RequiresExplicitCast<VALUETYPE>::cast(i[diff]);
    }
        /** Write the current data item. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR>
    void set(V value, ITERATOR const & i) const
        { *i = detail::RequiresExplicitCast<VALUETYPE>::cast(value); }

        /* This overload is needed to make the accessor work with a std::back_inserter */
    template <class V, class ITERATOR>
    void set(V value, ITERATOR & i) const
        { *i = detail::RequiresExplicitCast<VALUETYPE>::cast(value); }

        /** Write the data item at an offset (can be 1D or 2D or higher order difference)..
            The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR, class OFFSET>
    void set(V value, ITERATOR const & i, OFFSET const & diff) const
    {
        i[diff]= detail::RequiresExplicitCast<VALUETYPE>::cast(value);
    }
};

/********************************************************/
/*                                                      */
/*                StandardConstAccessor                 */
/*                                                      */
/********************************************************/

/** \brief Encapsulate read access to the values an iterator points to.

    StandardConstAccessor is a trivial accessor that simply encapsulates
    the iterator's operator*() and operator[]() in its
    read functions. It passes its arguments <em>by reference</em>.
    If the iterator returns its items by value (such as \ref vigra::CoordinateIterator), you
    must use StandardConstValueAccessor instead of StandardConstAccessor.
    Both accessors also have different optimization properties --
    StandardConstAccessor is usually faster for compound pixel types,
    while StandardConstValueAccessor is faster for the built-in types.

    <b>\#include</b> \<vigra/accessor.hxx\><br>
    Namespace: vigra
*/
template <class VALUETYPE>
class StandardConstAccessor
{
  public:
    typedef VALUETYPE value_type;

        /** read the current data item
        */
    template <class ITERATOR>
    VALUETYPE const & operator()(ITERATOR const & i) const
        { return *i; }

        /** read the data item at an offset (can be 1D or 2D or higher order difference).
        */
    template <class ITERATOR, class OFFSET>
    VALUETYPE const & operator()(ITERATOR const & i, OFFSET const & diff) const
    {
        return i[diff];
    }
};

/** \brief Encapsulate access to the values an iterator points to.

    StandardConstValueAccessor is a trivial accessor that simply encapsulates
    the iterator's operator*() and operator[]() in its
    read functions. It passes its arguments <em>by value</em>.
    If the iterator returns its items by reference (such as \ref vigra::ConstImageIterator),
    you can also use StandardConstAccessor.
    These accessors have different optimization properties --
    StandardConstAccessor is usually faster for compound pixel types,
    while StandardConstValueAccessor is faster for the built-in types.

    When an iterator passes a floating point number to an accessor
    with integral value_type, the value is rounded and clipped as appropriate.

    <b>\#include</b> \<vigra/accessor.hxx\><br>
    Namespace: vigra
*/
template <class VALUETYPE>
class StandardConstValueAccessor
{
  public:
    typedef VALUETYPE value_type;

        /** Read the current data item. The type <TT>ITERATOR::reference</TT>
            is automatically converted to <TT>VALUETYPE</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class ITERATOR>
    VALUETYPE operator()(ITERATOR const & i) const
        { return detail::RequiresExplicitCast<VALUETYPE>::cast(*i); }

        /** Read the data item at an offset (can be 1D or 2D or higher order difference).
            The type <TT>ITERATOR::index_reference</TT>
            is automatically converted to <TT>VALUETYPE</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class ITERATOR, class OFFSET>
    VALUETYPE operator()(ITERATOR const & i, OFFSET const & diff) const
    {
        return detail::RequiresExplicitCast<VALUETYPE>::cast(i[diff]);
    }
};

/********************************************************/
/*                                                      */
/*                 VectorComponentAccessor              */
/*                                                      */
/********************************************************/

/** \brief Accessor for one component of a vector.

    This accessor allows to select a single component (a single 'band')
    of a vector valued pixel type. The pixel type must support
    <TT>operator[]</TT>. The index of the component to be selected
    is passed in the constructor. The accessor returns its items
    <em>by reference</em>. If you want to pass/return items by value,
    use VectorComponentValueAccessor. If a floating point number
    is assigned by means of an accessor with integral value_type, the
    value is rounded and clipped as appropriate.

    <b>Usage:</b>

    \code
    vigra::BRGBImage image(w,h);

    // init red channel with 255
    initImage(destImageRange(image,
                             VectorComponentAccessor<vigra::BRGBImage::value_type>(0)),
              255);
    \endcode

    <b>\#include</b> \<vigra/accessor.hxx\><br>
    Namespace: vigra

*/
template <class VECTORTYPE>
class VectorComponentAccessor
{
    int index_;
  public:
        /** the value_type
        */
    typedef typename VECTORTYPE::value_type value_type;

        /** determine the component to be accessed
        */
    VectorComponentAccessor(int index)
    : index_(index)
    {}

        /** read the current data item
        */
    template <class ITERATOR>
    value_type const & operator()(ITERATOR const & i) const
        { return (*i)[index_]; }

        /** read the data item at an offset (can be 1D or 2D or higher order difference).
        */
    template <class ITERATOR, class OFFSET>
    value_type const & operator()(ITERATOR const & i, OFFSET const & diff) const
    {
        return i[diff][index_];
    }

        /** Write the current data item. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR>
    void set(V const & value, ITERATOR const & i) const
    {
        (*i)[index_] = detail::RequiresExplicitCast<value_type>::cast(value);
    }

        /** Write the data item at an offset (can be 1D or 2D or higher order difference)..
            The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR, class OFFSET>
    void set(V const & value, ITERATOR const & i, OFFSET const & diff) const 
    { 
        i[diff][index_]= detail::RequiresExplicitCast<value_type>::cast(value); 
    }
    
        /** Reset the index to the given number.
        */
    void setIndex(int i)
    {
        index_ = i;
    }
};

/** \brief Accessor for one component of a vector.

    This accessor allows to select a single component (a single 'band')
    of a vector valued pixel type. The pixel type must support
    <TT>operator[]</TT>. The index of the component to be selected
    is passed in the constructor. The accessor returns its items
    <em>by value</em>. If you want to pass/return items by reference,
    use VectorComponentAccessor. If a floating point number
    is assigned by means of an accessor with integral value_type, the
    value is rounded and clipped as appropriate.

    <b>Usage:</b>

    \code
    vigra::BRGBImage image(w,h);

    // init red channel with 255
    initImage(destImageRange(image,
                             VectorComponentValueAccessor<vigra::BRGBImage::value_type>(0)),
              255);
    \endcode

    <b>\#include</b> \<vigra/accessor.hxx\><br>
    Namespace: vigra

*/
template <class VECTORTYPE>
class VectorComponentValueAccessor
{
    int index_;
  public:
        /** the value_type
        */
    typedef typename VECTORTYPE::value_type value_type;

        /** determine the component to be accessed
        */
    VectorComponentValueAccessor(int index)
    : index_(index)
    {}

        /** Read the current data item.
            The type <TT>ITERATOR::index_reference::value_type</TT>
            is automatically converted to <TT>value_type</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const
        { return detail::RequiresExplicitCast<value_type>::cast((*i)[index_]); }

        /** Read the data item at an offset (can be 1D or 2D or higher order difference).
            The type <TT>ITERATOR::index_reference::value_type</TT>
            is automatically converted to <TT>value_type</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class ITERATOR, class OFFSET>
    value_type operator()(ITERATOR const & i, OFFSET const & diff) const
    { 
        return detail::RequiresExplicitCast<value_type>::cast(i[diff][index_]); 
    }
    
        /** Write the current data item. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR>
    void set(V value, ITERATOR const & i) const 
    { 
        (*i)[index_] = detail::RequiresExplicitCast<value_type>::cast(value); 
    }
    
        /** Write the data item at an offset (can be 1D or 2D or higher order difference)..
            The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR, class OFFSET>
    void set(V value, ITERATOR const & i, OFFSET const & diff) const 
    { 
        i[diff][index_]= detail::RequiresExplicitCast<value_type>::cast(value); 
    }
    
        /** Reset the index to the given number.
        */
    void setIndex(int i)
    {
        index_ = i;
    }
};

/********************************************************/
/*                                                      */
/*                   VectorElementAccessor              */
/*                                                      */
/********************************************************/

/** \brief Accessor for one component of a vector.

    This works like VectorComponentAccessor, only the template parameters differ: 
    Here, we need a vector accessor type , whereas VectorComponentAccessor requires a vector type.

    <b>Usage:</b>
    
    \code
    vigra::BRGBImage image(w,h);
    
    // init red channel with 255
    initImage(destImageRange(image, 
                             VectorElementAccessor<vigra::BRGBImage::Accessor>(0)),
              255);
    \endcode
    
    <b>\#include</b> \<vigra/accessor.hxx\><br>
    Namespace: vigra
    
*/
template <class ACCESSOR>
class VectorElementAccessor
{
    int index_;
    ACCESSOR a_;
  public:
        /** the value_type
        */
    typedef typename ACCESSOR::component_type value_type;
    
        /** determine the component to be accessed
        */
    VectorElementAccessor(int index, ACCESSOR a = ACCESSOR())
    : index_(index),
      a_(a)
    {}
    
        /** read the current data item
        */
    template <class ITERATOR>
    value_type const & operator()(ITERATOR const & i) const 
        { return a_.getComponent(i, index_); }
    
        /** read the data item at an offset (can be 1D or 2D or higher order difference).
        */
    template <class ITERATOR, class OFFSET>
    value_type const & operator()(ITERATOR const & i, OFFSET const & diff) const
    { 
        return a_.getComponent(i, diff, index_); 
    }
    
        /** Write the current data item. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR>
    void set(V const & value, ITERATOR const & i) const 
    { 
        a_.setComponent(detail::RequiresExplicitCast<value_type>::cast(value), i, index_); 
    }

        /** Write the data item at an offset (can be 1D or 2D or higher order difference)..
            The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR, class OFFSET>
    void set(V const & value, ITERATOR const & i, OFFSET const & diff) const 
    { 
       a_.setComponent(detail::RequiresExplicitCast<value_type>::cast(value), i, diff, index_); 
    }
    
        /** Reset the index to the given number.
        */
    void setIndex(int i)
    {
        index_ = i;
    }
};

/********************************************************/
/*                                                      */
/*                    SequenceAccessor                  */
/*                                                      */
/********************************************************/

/** \brief Accessor for items that are STL compatible sequences.

    It encapsulates access to the sequences' begin() and end()
    functions.

    <b>Usage:</b>

    <b>\#include</b> \<vigra/accessor.hxx\><br>
    Namespace: vigra

    \code
    typedef std::list<std::list<int> > ListOfLists;

    ListOfLists ll;
    ...

    typedef vigra::SequenceAccessor<ListOfLists::value_type> ListOfListsAccessor;
    ListOfListsAccessor a;
    for(ListOfLists::iterator li = ll.begin(); li != ll.end(); ++li)
    {
        for(ListOfListsAccessor::iterator i = a.begin(li); i != a.end(li); ++i)
        {
            *i = 10;
        }
    }
    \endcode
*/
template <class SEQUENCE>
class SequenceAccessor
: public StandardAccessor<SEQUENCE>
{
  public:
        /** the sequence's value_type
        */
    typedef typename SEQUENCE::value_type component_type;

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION
    typedef typename
            If<typename TypeTraits<SEQUENCE>::isConst,
               typename SEQUENCE::const_iterator,
               typename SEQUENCE::iterator>::type 
            iterator;
#else
        /** the sequence's iterator type
        */
    typedef typename SEQUENCE::iterator iterator;
#endif

        /** get begin iterator for sequence at given iterator position
        */
    template <class ITERATOR>
    iterator begin(ITERATOR const & i) const
    {
        return (*i).begin();
    }

        /** get end iterator for sequence at given iterator position
        */
    template <class ITERATOR>
    iterator end(ITERATOR const & i)  const
    {
         return (*i).end();
    }

        /** get begin iterator for sequence at an offset
            of given iterator position
        */
    template <class ITERATOR, class OFFSET>
    iterator begin(ITERATOR const & i, OFFSET const & diff)  const
    {
        return i[diff].begin();
    }

        /** get end iterator for sequence at a 2D difference vector
            of given iterator position
        */
    template <class ITERATOR, class OFFSET>
    iterator end(ITERATOR const & i, OFFSET const & diff)  const
    {
        return i[diff].end();
    }

        /** get size of sequence at given iterator position
        */
    template <class ITERATOR>
    unsigned int size(ITERATOR const & i) const { return (*i).size(); }

        /** get size of sequence at 2D difference vector of given iterator position
        */
    template <class ITERATOR, class OFFSET>
    unsigned int size(ITERATOR const & i, OFFSET const & diff) const
    { return i[diff].size(); }
};

/********************************************************/
/*                                                      */
/*                     VectorAccessor                   */
/*                                                      */
/********************************************************/

/** \brief Accessor for items that are STL compatible vectors.

    It encapsulates access to a vector's access functionality.

    <b> Usage:</b>

    <b>\#include</b> \<vigra/accessor.hxx\><br>
    Namespace: vigra

    The accessor has two modes of operation:

    <ol>
    <li> Access the vector's iterator via the <TT>begin()</TT> and <TT>end()</TT>
    functions:

    \code
    typedef std::list<std::vector<int> > ListOfVectors;

    ListOfVectors ll;
    ...

    typedef vigra::SequenceAccessor<ListOfVectors::value_type> ListOfVectorsAccessor;
    ListOfVectorsAccessor a;
    for(ListOfVectors::iterator li = ll.begin(); li != ll.end(); ++li)
    {
        for(ListOfVectorsAccessor::iterator i = a.begin(li); i != a.end(li); ++i)
        {
            *i = 10;
        }
    }
    \endcode
    <li> Access the vector's components via an index (internally calls
    the vector's <TT>operator[]</TT> ):
    \code
    typedef std::list<std::vector<int> > ListOfVectors;

    ListOfVectors ll;
    ...

    typedef vigra::SequenceAccessor<ListOfVectors::value_type> ListOfVectorsAccessor;
    ListOfVectorsAccessor a;
    for(ListOfVectors::iterator li = ll.begin(); li != ll.end(); ++li)
    {
        for(int i = 0; i != a.size(li); ++i)
        {
            a.setComponent(10, li, i);
        }
    }
    \endcode
    </ol>

    <b> Required Interface:</b>

    \code
    VECTOR v;
    VECTOR::iterator i;
    value_type d;
    int index;

    d = v[index];
    v[index] = d;
    i = v.begin();
    i = v.end();
    v.size();
    \endcode
*/
template <class VECTOR>
class VectorAccessor
: public SequenceAccessor<VECTOR>
{
  public:
        /** the vector's value_type
        */
    typedef typename VECTOR::value_type component_type;

        /** the vector element accessor associated with this vector accessor
            (see \ref VectorElementAccessor)
        */
    typedef VectorElementAccessor<VectorAccessor<VECTOR> > ElementAccessor;

        /** Read the component data at given vector index
            at given iterator position
        */
    template <class ITERATOR>
    component_type const & getComponent(ITERATOR const & i, int idx) const
    {
        return (*i)[idx];
    }

        /** Set the component data at given vector index
            at given iterator position. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
        */
    template <class V, class ITERATOR>
    void setComponent(V const & value, ITERATOR const & i, int idx) const
    {
        (*i)[idx] = detail::RequiresExplicitCast<component_type>::cast(value);
    }

        /** Read the component data at given vector index
            at an offset of given iterator position
        */
    template <class ITERATOR, class OFFSET>
    component_type const & getComponent(ITERATOR const & i, OFFSET const & diff, int idx) const
    {
        return i[diff][idx];
    }

    /** Set the component data at given vector index
        at an offset of given iterator position. The type <TT>V</TT> of the passed
        in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
            In case of a conversion floating point -> integral this includes rounding and clipping.
    */
    template <class V, class ITERATOR, class OFFSET>
    void
    setComponent(V const & value, ITERATOR const & i, OFFSET const & diff, int idx) const
    {
        i[diff][idx] = detail::RequiresExplicitCast<component_type>::cast(value);
    }
};


/********************************************************/
/*                                                      */
/*                 MultiImageAccessor2                  */
/*                                                      */
/********************************************************/

/** \brief Access two images simultaneously.

    This accessor is used when two images need to be treated as one
    because an algorithm accepts only one image. For example,
    \ref seededRegionGrowing() uses only one image two calculate
    the cost for aggregating each pixel into a region. Sometimes, we
    need more information to calculate this cost, for example gray value
    and local gradient magnitude. These values can be stored in two images,
    which appear as only one when we pass a <TT>MultiImageAccessor2</TT> to
    the algorithms. Of course, the cost functor must accept a <TT>pair</TT>
    of values for this to work. Instead of an actual image iterator, we
    pass a <a href="CoordinateIterator.html">CoordinateIterator</a> which
    selects the right pixels form both images.

    <b> Usage:</b>

    <b>\#include</b> \<vigra/accessor.hxx\><br>
    Namespace: vigra

    \code
    using namespace vigra;

    FImage gray_values(w,h), gradient_magnitude(w,h);
    IImage seeds(w,h), labels(w,h);

    seededRegionGrowing(
        srcIterRange(CoordinateIterator(), CoordinateIterator(w,h),
           MultiImageAccessor2<FImage::iterator, FImage::Accessor,
                               FImage::iterator, FImage::Accessor>
                              (gray_values.upperLeft(), gray_values.accessor(),
                               gradient_magnitude.upperLeft(), gradient_magnitude.accessor())),
        srcImage(seeds),
        destImage(labels),
        SomeCostFunctor());
    \endcode
*/

template <class Iter1, class Acc1, class Iter2, class Acc2>
class MultiImageAccessor2
{
  public:
        /** The accessors value_type: construct a pair that contains
            the corresponding image values.
        */
    typedef pair<typename Acc1::value_type, typename Acc2::value_type>
            value_type;

        /** Construct from two image iterators and associated accessors.
        */
    MultiImageAccessor2(Iter1 i1, Acc1 a1, Iter2 i2, Acc2 a2)
    : i1_(i1), a1_(a1), i2_(i2), a2_(a2)
    {}

        /** read the current data item
        */
    template <class OFFSET>
    value_type operator()(OFFSET const & d) const
    {
        return std::make_pair(a1_(i1_, d), a2_(i2_, d));
    }

        /** read the data item at an offset
        */
    template <class OFFSET1, class OFFSET2>
    value_type operator()(OFFSET1 d1, OFFSET2 const & d2) const
    {
        d1 += d2;
        return std::make_pair(a1_(i1_, d1), a2_(i2_, d1));
    }

  private:
    Iter1 i1_;
    Acc1 a1_;
    Iter2 i2_;
    Acc2 a2_;
};

//@}

template <class T>
struct AccessorTraits
{
    typedef StandardAccessor<T>        default_accessor;
    typedef StandardConstAccessor<T>   default_const_accessor;
};

#define VIGRA_DEFINE_ACCESSOR_TRAITS(VALUE, ACCESSOR, CONST_ACCESSOR) \
    template <> \
    struct AccessorTraits<VALUE > \
    { \
        typedef ACCESSOR<VALUE >         default_accessor; \
        typedef CONST_ACCESSOR<VALUE >   default_const_accessor; \
    };

VIGRA_DEFINE_ACCESSOR_TRAITS(signed char, StandardValueAccessor, StandardConstValueAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(unsigned char, StandardValueAccessor, StandardConstValueAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(short, StandardValueAccessor, StandardConstValueAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(unsigned short, StandardValueAccessor, StandardConstValueAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(int, StandardValueAccessor, StandardConstValueAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(unsigned int, StandardValueAccessor, StandardConstValueAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(long, StandardValueAccessor, StandardConstValueAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(unsigned long, StandardValueAccessor, StandardConstValueAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(float, StandardValueAccessor, StandardConstValueAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(double, StandardValueAccessor, StandardConstValueAccessor)

template <class T, unsigned int RED_IDX, unsigned int GREEN_IDX, unsigned int BLUE_IDX> class RGBValue;
template <class T> class RGBAccessor;
template <class T, int SIZE> class TinyVector;

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class T, unsigned int RED_IDX, unsigned int GREEN_IDX, unsigned int BLUE_IDX>
struct AccessorTraits<RGBValue<T, RED_IDX, GREEN_IDX, BLUE_IDX> >
{
    typedef RGBAccessor<RGBValue<T, RED_IDX, GREEN_IDX, BLUE_IDX> >   default_accessor;
    typedef RGBAccessor<RGBValue<T, RED_IDX, GREEN_IDX, BLUE_IDX> >   default_const_accessor;
};

template <class T, int SIZE>
struct AccessorTraits<TinyVector<T, SIZE> >
{
    typedef VectorAccessor<TinyVector<T, SIZE> >   default_accessor;
    typedef VectorAccessor<TinyVector<T, SIZE> >   default_const_accessor;
};

#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION

VIGRA_DEFINE_ACCESSOR_TRAITS(RGBValue<unsigned char>, RGBAccessor, RGBAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(RGBValue<signed char>, RGBAccessor, RGBAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(RGBValue<short>, RGBAccessor, RGBAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(RGBValue<unsigned short>, RGBAccessor, RGBAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(RGBValue<int>, RGBAccessor, RGBAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(RGBValue<unsigned int>, RGBAccessor, RGBAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(RGBValue<long>, RGBAccessor, RGBAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(RGBValue<unsigned long>, RGBAccessor, RGBAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(RGBValue<float>, RGBAccessor, RGBAccessor)
VIGRA_DEFINE_ACCESSOR_TRAITS(RGBValue<double>, RGBAccessor, RGBAccessor)

#define VIGRA_PIXELTYPE TinyVector<unsigned char, 2>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<unsigned char, 3>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<unsigned char, 4>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<short, 2>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<short, 3>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<short, 4>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<int, 2>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<int, 3>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<int, 4>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<float, 2>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<float, 3>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<float, 4>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<double, 2>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<double, 3>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<double, 4>
VIGRA_DEFINE_ACCESSOR_TRAITS(VIGRA_PIXELTYPE, VectorAccessor, VectorAccessor)
#undef VIGRA_PIXELTYPE

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

#undef VIGRA_DEFINE_ACCESSOR_TRAITS

} // namespace vigra

#endif // VIGRA_ACCESSOR_HXX
