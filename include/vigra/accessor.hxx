/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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
 
#ifndef VIGRA_ACCESSOR_HXX
#define VIGRA_ACCESSOR_HXX

#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"

namespace vigra {

namespace detail {

template <class T>
struct RequiresExplicitCast {
    template <class U>
    static U const & cast(U const & v)
        { return v; }
};

template <>
struct RequiresExplicitCast<signed char> {
    template <class U>
    static signed char cast(U const & v)
        { return static_cast<signed char>(v); }
};

template <>
struct RequiresExplicitCast<unsigned char> {
    template <class U>
    static unsigned char cast(U const & v)
        { return static_cast<unsigned char>(v); }
};

template <>
struct RequiresExplicitCast<short> {
    template <class U>
    static short cast(U const & v)
        { return static_cast<short>(v); }
};

template <>
struct RequiresExplicitCast<unsigned short> {
    template <class U>
    static unsigned short cast(U const & v)
        { return static_cast<unsigned short>(v); }
};

template <>
struct RequiresExplicitCast<int> {
    template <class U>
    static int cast(U const & v)
        { return static_cast<int>(v); }
};

template <>
struct RequiresExplicitCast<unsigned int> {
    template <class U>
    static unsigned int cast(U const & v)
        { return static_cast<unsigned int>(v); }
};

template <>
struct RequiresExplicitCast<long> {
    template <class U>
    static long cast(U const & v)
        { return static_cast<long>(v); }
};

template <>
struct RequiresExplicitCast<unsigned long> {
    template <class U>
    static unsigned long cast(U const & v)
        { return static_cast<unsigned long>(v); }
};

} // namespace detail

/** \addtogroup DataAccessors Data Accessors

    Basic templates to encapsulate access to the data of an iterator.
    
    Data accessors are used to allow for flexible access to the data
    an interator points to. When we access the data directly, we
    are bound to what <TT>operator*</TT> returns, if this method exists at 
    all. Encapsulating access in an accessor enables a better
    decoupling of data structures and algorithms. 
    <a href="documents/DataAccessors.ps">This paper</a> contains
    a detailed description of the concept.
*/
//@{

/********************************************************/
/*                                                      */
/*                     StandardAccessor                 */
/*                                                      */
/********************************************************/

/** \brief <TT>%StandardAccessor</TT> and <TT>StandardValueAccessor</TT>.

    <TT>StandardAccessor</TT> is normally used with all relizations of 
    \ref vigra::ImageIterator. It
    simply encapsulates a call to the iterator's <TT>operator*()</TT> in both
    the read and write function. If the iterator returns its items by 
    value (such as \ref vigra::CoordinateIterator), you must use
    <TT>StandardValueAccessor</TT> instead of <TT>StandardAccessor</TT> (both have 
    the same interface, but the latter can be optimized better).

    <b>\#include</b> "<a href="accessor_8hxx-source.html">vigra/accessor.hxx</a>"
    
    Namespace: vigra
    
*/
template <class VALUETYPE>
class StandardAccessor
{
  public:
        /** the value_type
        */
    typedef VALUETYPE value_type;
    typedef typename NumericTraits<VALUETYPE>::Promote Promote;
    typedef typename NumericTraits<VALUETYPE>::RealPromote RealPromote;
    
        /** read the current data item
        */
    template <class ITERATOR>
    VALUETYPE const & operator()(ITERATOR & i) const { return *i; }
    
        /** read the data item at a distance (can be 1D or 2D or higher distance)
        */
    template <class ITERATOR, class DISTANCE>
    VALUETYPE const & operator()(ITERATOR & i, DISTANCE const & dist) const
    { 
        return i[dist]; 
    }
    
        /** Write the current data item. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
        */
    template <class V, class ITERATOR>
    void set(V const & value, ITERATOR & i) const 
    { *i = detail::RequiresExplicitCast<VALUETYPE>::cast(value); }
    
        /** Write the data item at a distance (can be 1D or 2D or higher distance).
            The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
        */
    template <class V, class ITERATOR, class DISTANCE>
    void set(V const & value, ITERATOR & i, DISTANCE const & dist) const 
    { 
        i[dist]= detail::RequiresExplicitCast<VALUETYPE>::cast(value); 
    }
};

template <class VALUETYPE>
class StandardValueAccessor
{
  public:
        /* the value_type
        */
    typedef VALUETYPE value_type;
    
        /* read the current data item
        */
    template <class ITERATOR>
    VALUETYPE operator()(ITERATOR & i) const { return *i; }
    
        /* read the data item at a distance (can be 1D or 2D or higher distance)
        */
    template <class ITERATOR, class DISTANCE>
    VALUETYPE operator()(ITERATOR & i, DISTANCE const & dist) const
    { 
        return i[dist]; 
    }
    
        /** Write the current data item. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
        */
    template <class V, class ITERATOR>
    void set(V const & value, ITERATOR & i) const { *i = detail::RequiresExplicitCast<VALUETYPE>::cast(value); }
    
        /** Write the data item at a distance (can be 1D or 2D or higher distance).
            The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
        */
    template <class V, class ITERATOR, class DISTANCE>
    void set(V const & value, ITERATOR & i, DISTANCE const & dist) const 
    { 
        i[dist]= detail::RequiresExplicitCast<VALUETYPE>::cast(value); 
    }
};

/********************************************************/
/*                                                      */
/*                StandardConstAccessor                 */
/*                                                      */
/********************************************************/

/** \brief <TT>%StandardConstAccessor</TT> and <TT>StandardConstValueAccessor</TT>.

    <TT>StandardConstAccessor</TT> is normally used with all 
    realizations of \ref vigra::ConstImageIterator. It
    simply encapsulates a call to the iterator's <TT>operator*()</TT> in its
    read function (as a const accessor, it does not have a write
    function). If the iterator returns its items by 
    value (such as the \ref vigra::CoordinateIterator), you must use
    <TT>StandardConstValueAccessor</TT> instead of <TT>StandardConstAccessor</TT> (both have 
    the same interface, but the latter can be optimized better).
    
    <b>\#include</b> "<a href="accessor_8hxx-source.html">vigra/accessor.hxx</a>"
    
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
    VALUETYPE const & operator()(ITERATOR & i) const { return *i; }
    
        /** read the data item at a distance (can be 1D or 2D or higher distance)
        */
    template <class ITERATOR, class DISTANCE>
    VALUETYPE const & operator()(ITERATOR & i, DISTANCE const & dist) const
    { 
        return i[dist]; 
    }
};

template <class VALUETYPE>
class StandardConstValueAccessor
{
  public:
    typedef VALUETYPE value_type;
    
        /* read the current data item
        */
    template <class ITERATOR>
    VALUETYPE operator()(ITERATOR & i) const { return *i; }
    
        /* read the data item at a distance (can be 1D or 2D or higher distance)
        */
    template <class ITERATOR, class DISTANCE>
    VALUETYPE operator()(ITERATOR & i, DISTANCE const & dist) const
    { 
        return i[dist]; 
    }
};

/********************************************************/
/*                                                      */
/*                 VectorComponentAccessor              */
/*                                                      */
/********************************************************/

/** \brief Accessor for one component of a vector.

    Usage:

    <b>\#include</b> "<a href="accessor_8hxx-source.html">vigra/accessor.hxx</a>"
    
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
    value_type const & operator()(ITERATOR & i) const { return (*i)[index_]; }
    
        /** read the data item at a distance (can be 1D or 2D or higher distance)
        */
    template <class ITERATOR, class DISTANCE>
    value_type const & operator()(ITERATOR & i, DISTANCE const & dist) const
    { 
        return i[dist][index_]; 
    }
    
        /** Write the current data item. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
        */
    template <class V, class ITERATOR>
    void set(V const & value, ITERATOR & i) const 
    { 
        (*i)[index_] = detail::RequiresExplicitCast<value_type>::cast(value); 
    }
    
        /** Write the data item at a distance (can be 1D or 2D or higher distance).
            The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
        */
    template <class V, class ITERATOR, class DISTANCE>
    void set(V const & value, ITERATOR & i, DISTANCE const & dist) const 
    { 
        i[dist][index_]= detail::RequiresExplicitCast<value_type>::cast(value); 
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
    
    Usage:

    <b>\#include</b> "<a href="accessor_8hxx-source.html">vigra/accessor.hxx</a>"
    
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

    /** the sequence's iterator type
    */
    typedef typename SEQUENCE::iterator iterator;
    
    /** get begin iterator for sequence at given iterator position
    */
    template <class ITERATOR>
    iterator begin(ITERATOR & i) const
    { 
        return (*i).begin(); 
    }
    
    /** get end iterator for sequence at given iterator position
    */
    template <class ITERATOR>
    iterator end(ITERATOR & i)  const
    {
         return (*i).end(); 
    }
    
    /** get begin iterator for sequence at a distance
        of given iterator position
    */
    template <class ITERATOR, class DISTANCE>
    iterator begin(ITERATOR & i, DISTANCE const & dist)  const
    { 
        return i[dist].begin(); 
    }
    
    /** get end iterator for sequence at a 2D distance
        of given iterator position
    */
    template <class ITERATOR, class DISTANCE>
    iterator end(ITERATOR & i, DISTANCE const & dist)  const
    { 
        return i[dist].end(); 
    }

    /** get size of sequence at given iterator position
    */
    template <class ITERATOR>
    int size(ITERATOR & i) const { return (*i).size(); }

    /** get size of sequence at 2D distance of given iterator position
    */
    template <class ITERATOR, class DISTANCE>
    int size(ITERATOR & i, DISTANCE const & dist) const { return i[dist].size(); }
};

/********************************************************/
/*                                                      */
/*                     VectorAccessor                   */
/*                                                      */
/********************************************************/

/** \brief Accessor for items that are STL compatible vectors.

    It encapsulates access to a vector's access functionality.
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="accessor_8hxx-source.html">vigra/accessor.hxx</a>"
    
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

        /** Read the component data at given vector index
            at given iterator position 
        */
    template <class ITERATOR>
    component_type getComponent(ITERATOR & i, int idx) const 
    { 
        return (*i)[idx]; 
    }
    
        /** Set the component data at given vector index
            at given iterator position. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
        */
    template <class V, class ITERATOR>
    void setComponent(V const & value, ITERATOR & i, int idx) const
    { 
        (*i)[idx] = detail::RequiresExplicitCast<component_type>::cast(value); 
    }
    
        /** Read the component data at given vector index
            at a distance of given iterator position
        */
    template <class ITERATOR, class DISTANCE>
    component_type getComponent(ITERATOR & i, DISTANCE const & dist, int idx) const
    { 
        return i[dist][idx]; 
    }
    
    /** Set the component data at given vector index
        at a distance of given iterator position. The type <TT>V</TT> of the passed
        in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
    */
    template <class V, class ITERATOR, class DISTANCE>
    void 
    setComponent(V const & value, ITERATOR & i, DISTANCE const & dist, int idx) const 
    { 
        i[dist][idx] = detail::RequiresExplicitCast<component_type>::cast(value); 
    }
};


/********************************************************/
/*                                                      */
/*                  InterpolatingAccessor               */
/*                                                      */
/********************************************************/

/** \brief Bilinear interpolation at non-integer positions.

    This accessor allows an image be accessed at arbitrary non-integer
    coordinates and performs an bi-linear interpolation to
    obtain a pixel value.
    It uses the given ACCESSOR (which is usually the
    accessor originally associated with the iterator)
    to access data.
    
    <b>\#include</b> "<a href="accessor_8hxx-source.html">vigra/accessor.hxx</a>"
    
    Namespace: vigra
    
    <b> Required Interface:</b>
    
    \code
    ITERATOR iter;
    ACCESSOR a;
    VALUETYPE destvalue;
    float s;
    int x, y;
    
    destvalue = s * a(iter, x, y) + s * a(iter, x, y);
    
    \endcode
*/
template <class ACCESSOR, class VALUETYPE>
class BilinearInterpolatingAccessor
{
  public:
    /** the iterators' pixel type
    */
    typedef VALUETYPE value_type;
    
    /** init from given accessor
    */
    BilinearInterpolatingAccessor(ACCESSOR a)
    : a_(a)
    {}
    
    /** Interpolate the data item at a non-integer distance.
        Ensure that no outside pixels are accessed if 
        (x, y) is near the image border (as long as
        0 <= x <= width-1, 0 <= y <= height-1).
    */
    template <class ITERATOR>
    value_type operator()(ITERATOR & i, float x, float y) const 
    { 
        int ix = x;
        int iy = y;
        float dx = x - ix;
        float dy = y - iy;
        
        value_type ret;
        
        // avoid dereferencing the iterator outside its range
        if(dx == 0.0)
        {
            if(dy == 0.0)
            {
                ret = a_(i, Diff2D(ix, iy));
            }
            else
            {
                ret = (1.0 - dy) * a_(i, Diff2D(ix, iy)) +
                  dy * a_(i, Diff2D(ix, iy + 1));
            }
        }
        else
        {
            if(dy == 0.0)
            {
                ret = (1.0 - dx) * a_(i, Diff2D(ix, iy)) + 
                  dx * a_(i, Diff2D(ix + 1, iy));
            }
            else
            {
                ret = (1.0 - dx) * (1.0 - dy) * a_(i, Diff2D(ix, iy)) +
                  dx * (1.0 - dy) * a_(i, Diff2D(ix + 1, iy)) +
                  (1.0 - dx) * dy * a_(i, Diff2D(ix, iy + 1)) +
                  dx * dy * a_(i, Diff2D(ix + 1, iy + 1));
            }
        }
            
        return ret;
    }

    /** Interpolate the data item at a non-integer distance.
        This function works as long as 0 <= x < width-1, 
        0 <= y < height-1. It is slightly faster than <TT>operator()</TT>.
    */
    template <class ITERATOR>
    value_type unchecked(ITERATOR & i, float x, float y) const 
    { 
    int ix = x;
        int iy = y;
        float dx = x - ix;
        float dy = y - iy;
        return (1.0 - dx) * (1.0 - dy) * a_(i, Diff2D(ix, iy)) +
               dx * (1.0 - dy) * a_(i, Diff2D(ix + 1, iy)) +
               (1.0 - dx) * dy * a_(i, Diff2D(ix, iy + 1)) +
               dx * dy * a_(i, Diff2D(ix + 1, iy + 1));
    }
    
  private:
    ACCESSOR a_;
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
    the cost for aggregating each pixel into a region. Somtimes, we
    need more information to calcuate this cost, for example gray value
    and local gradient magnitude. These values can be stored in two images,
    which appear as only one when we pass a <TT>MultiImageAccessor2</TT> to
    the lagorithms. Of course, the cost functor must accept a <TT>pair</TT> 
    of values for this to work. Instead of an actual image iterator, we
    pass a <a href="CoordinateIterator.html">CoordinateIterator</a> which 
    selects the right pixels form both images.
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="accessor_8hxx-source.html">vigra/accessor.hxx</a>"
    
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
    template <class DISTANCE>
    value_type operator()(DISTANCE const & d) const
    { 
        return std::make_pair(a1_(i1_, d), a2_(i2_, i.x, i.y)); 
    }
    
        /** read the data item at a distance
        */
    template <class DISTANCE1, class DISTANCE2>
    value_type operator()(DISTANCE1 const & d1, DISTANCE2 d2) const
    { 
        d2 += d1;
        return std::make_pair(a1_(i1_, d2), a2_(i2_, d2)); 
    }
    
  private:
    Iter1 i1_;
    Acc1 a1_;
    Iter2 i2_;
    Acc2 a2_;
};

//@}

} // namespace vigra

#endif // VIGRA_ACCESSOR_HXX
