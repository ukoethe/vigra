/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
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

/** @name Data Accessors

    Data accessors are used to allow for flexible access to the data
    an interator points to. When we access the data directly, we
    are bound to what #operator*# returns, if this method exists at 
    all. Encapsulating access in an accessor enables a better
    decoupling of data structures and algorithms. 
    \URL[This paper]{documents/DataAccessors.ps} contains
    a detailed description of the concept.
    
    @memo Basic templates to encapsulate access to the data of an iterator
*/
//@{

/********************************************************/
/*                                                      */
/*                     StandardAccessor                 */
/*                                                      */
/********************************************************/

/** StandardAccessor and StandardValueAccessor.
    #StandardAccessor# is normally used with all \Ref{ImageIterator}s. It
    simply encapsulates a call to the iterator's #operator*()# in both
    the read and write function. If the iterator returns its items by 
    value (such as the \Ref{CoordinateIterator}), you must use
    #StandardValueAccessor# instead of #StandardAccessor# (both have 
    the same interface, but the latter can be optimized better).

    Include-File:
    \URL[vigra/accessor.hxx]{../include/vigra/accessor.hxx}
    
*/
template <class VALUETYPE>
class StandardAccessor
{
  public:
    /** the value_type
        @memo
    */
    typedef VALUETYPE value_type;
    
    /** read the current data item
        @memo
    */
    template <class ITERATOR>
    VALUETYPE const & operator()(ITERATOR & i) const { return *i; }
    
    /** read the data item at a distance (can be 1D or 2D or higher distance)
        @memo
    */
    template <class ITERATOR, class DISTANCE>
    VALUETYPE const & operator()(ITERATOR & i, DISTANCE const & dist) const
    { 
        return i[dist]; 
    }
    
    /** write the current data item
        @memo
    */
    template <class ITERATOR>
    void set(VALUETYPE const & v, ITERATOR & i) const { *i = v; }
    
    /** write the data item at a distance (can be 1D or 2D or higher distance)
        @memo
    */
    template <class ITERATOR, class DISTANCE>
    void set(VALUETYPE const & v, ITERATOR & i, DISTANCE const & dist) const 
    { 
        i[dist]= v; 
    }
};

template <class VALUETYPE>
class StandardValueAccessor
{
  public:
    /* the value_type
        @memo
    */
    typedef VALUETYPE value_type;
    
    /* read the current data item
        @memo
    */
    template <class ITERATOR>
    VALUETYPE operator()(ITERATOR & i) const { return *i; }
    
    /* read the data item at a distance (can be 1D or 2D or higher distance)
        @memo
    */
    template <class ITERATOR, class DISTANCE>
    VALUETYPE operator()(ITERATOR & i, DISTANCE const & dist) const
    { 
        return i[dist]; 
    }
    
    /* write the current data item
        @memo
    */
    template <class ITERATOR>
    void set(VALUETYPE const & v, ITERATOR & i) const { *i = v; }
    
    /* write the data item at a distance (can be 1D or 2D or higher distance)
        @memo
    */
    template <class ITERATOR, class DISTANCE>
    void set(VALUETYPE const & v, ITERATOR & i, DISTANCE const & dist) const 
    { 
        i[dist]= v; 
    }
};

/********************************************************/
/*                                                      */
/*                StandardConstAccessor                 */
/*                                                      */
/********************************************************/

/** StandardConstAccessor and StandardConstValueAccessor.
    #StandardConstAccessor# is normally used with all \Ref{ConstImageIterator}s. It
    simply encapsulates a call to the iterator's #operator*()# in its
    read function (as a const accessor, it does not have a write
    function). If the iterator returns its items by 
    value (such as the \Ref{CoordinateIterator}), you must use
    #StandardConstValueAccessor# instead of #StandardConstAccessor# (both have 
    the same interface, but the latter can be optimized better).
    
    Include-File:
    \URL[vigra/accessor.hxx]{../include/vigra/accessor.hxx}
    
*/
template <class VALUETYPE>
class StandardConstAccessor
{
  public:
    typedef VALUETYPE value_type;
    
    /** read the current data item
        @memo
    */
    template <class ITERATOR>
    VALUETYPE const & operator()(ITERATOR & i) const { return *i; }
    
    /** read the data item at a distance (can be 1D or 2D or higher distance)
        @memo
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
        @memo
    */
    template <class ITERATOR>
    VALUETYPE operator()(ITERATOR & i) const { return *i; }
    
    /* read the data item at a distance (can be 1D or 2D or higher distance)
        @memo
    */
    template <class ITERATOR, class DISTANCE>
    VALUETYPE operator()(ITERATOR & i, DISTANCE const & dist) const
    { 
        return i[dist]; 
    }
};

/********************************************************/
/*                                                      */
/*                    SequenceAccessor                  */
/*                                                      */
/********************************************************/

/** Accessor for items that are STL compatible sequences.
    It encapsulates access to the sequences' begin() and end()
    functions.
    
    Usage:

    Include-File:
    \URL[vigra/accessor.hxx]{../include/vigra/accessor.hxx}
    
    \begin{verbatim}
    Iterator iter;
    SequenceAccessor<Iterator, Iterator::PixelType> a;
    for(SequenceAccessor<Iterator, Iterator::PixelType>::iterator i = a.begin(iter); 
        i != a.end(iter); ++i) 
    {
    *i = 10;
    }
    \end{verbatim}
*/
template <class SEQUENCE>
class SequenceAccessor
: public StandardAccessor<SEQUENCE>
{
  public:
    /** the sequence's value_type
        @memo
    */
    typedef typename SEQUENCE::value_type component_type;

    /** the sequence's iterator type
        @memo
    */
    typedef typename SEQUENCE::iterator iterator;
    
    /** get begin iterator for sequence at given iterator position
        @memo
    */
    template <class ITERATOR>
    iterator begin(ITERATOR & i) const
    { 
        return (*i).begin(); 
    }
    
    /** get end iterator for sequence at given iterator position
        @memo
    */
    template <class ITERATOR>
    iterator end(ITERATOR & i)  const
    {
         return (*i).end(); 
    }
    
    /** get begin iterator for sequence at a distance
        of given iterator position
        @memo
    */
    template <class ITERATOR, class DISTANCE>
    iterator begin(ITERATOR & i, DISTANCE const & dist)  const
    { 
        return i[dist].begin(); 
    }
    
    /** get end iterator for sequence at a 2D distance
        of given iterator position
        @memo
    */
    template <class ITERATOR, class DISTANCE>
    iterator end(ITERATOR & i, DISTANCE const & dist)  const
    { 
        return i[dist].end(); 
    }

    /** get size of sequence at given iterator position
        @memo
    */
    template <class ITERATOR>
    int size(ITERATOR & i) const { return (*i).size(); }

    /** get size of sequence at 2D distance of given iterator position
        @memo
    */
    template <class ITERATOR, class DISTANCE>
    int size(ITERATOR & i, DISTANCE const & dist) const { return i[dist].size(); }
};

/********************************************************/
/*                                                      */
/*                     VectorAccessor                   */
/*                                                      */
/********************************************************/

/** Accessor for items that are STL compatible vectors.
    It encapsulates access to a vector's access functionality.
    
    {\bf Usage:}
    
    Include-File:
    \URL[vigra/accessor.hxx]{../include/vigra/accessor.hxx}
        
    The accessor has two modes of operation:
    
    \begin{enumerate}
    \item Access the vector's iterator via the #begin()# and #end()#
    functions:
    
    \begin{verbatim}
    Iterator iter;
    VectorAccessor<Iterator, Iterator::PixelType> a;
    for(VectorAccessor<Iterator, Iterator::PixelType>::iterator i = a.begin(iter); 
        i != a.end(iter); ++i) 
    {
    *i = 10;
    }
    \end{verbatim}
    \item Access the vector's components via an index (internally calls 
    the vector's #operator[]# ):
    \begin{verbatim}
    Iterator iter;
    VectorAccessor<Iterator, Iterator::PixelType> a;
    for(int i=0; i < a.size(iter); ++i) a.setComponent(i, iter) = 10;
    \end{verbatim}
    \end{enumerate}
    
    {\bf Required Interface:}
    
    \begin{verbatim}
    VECTOR v;
    VECTOR::iterator i;
    value_type d;
    int index;
    
    d = v[index];
    v[index] = d;
    i = v.begin();
    i = v.end();
    v.size();
    \end{verbatim}
*/
template <class VECTOR>
class VectorAccessor
: public SequenceAccessor<VECTOR>
{
  public:
    /** the vector's value_type
        @memo
    */
    typedef typename VECTOR::value_type component_type;

    /** read the component data at given vector index
        at given iterator position 
        @memo
    */
    template <class ITERATOR>
    component_type getComponent(ITERATOR & i, int idx) const 
    { 
        return (*i)[idx]; 
    }
    
    /** set the component data at given vector index
        at given iterator position 
        @memo
    */
    template <class ITERATOR>
    void setComponent(component_type v, ITERATOR & i, int idx) const
    { 
        (*i)[idx] = v; 
    }
    
    /** read the component data at given vector index
        at a distance of given iterator position
        @memo
    */
    template <class ITERATOR, class DISTANCE>
    component_type getComponent(ITERATOR & i, DISTANCE const & dist, int idx) const
    { 
        return i[dist][idx]; 
    }
    
    /** set the component data at given vector index
        at a distance of given iterator position
        @memo
    */
    template <class ITERATOR, class DISTANCE>
    void 
    setComponent(component_type v, ITERATOR & i, DISTANCE const & dist, int idx) const 
    { 
        i[dist][idx] = v; 
    }
};


/********************************************************/
/*                                                      */
/*                  InterpolatingAccessor               */
/*                                                      */
/********************************************************/

/** Bilinear interpolation at non-integer positions.
    This accessor allows an image be accessed at arbitrary non-integer
    coordinates and performs an bi-linear interpolation to
    obtain a pixel value.
    It uses the given ACCESSOR (which is usually the
    accessor originally associated with the iterator)
    to access data.
    
    Include-File:
    \URL[vigra/accessor.hxx]{../include/vigra/accessor.hxx}
    
    {\bf Required Interface:}
    
    \begin{verbatim}
    ITERATOR iter;
    ACCESSOR a;
    VALUETYPE destvalue;
    float s;
    int x, y;
    
    destvalue = s * a(iter, x, y) + s * a(iter, x, y);
    
    \end{verbatim}
*/
template <class ACCESSOR, class VALUETYPE>
class BilinearInterpolatingAccessor
{
  public:
    /** the iterators' pixel type
        @memo
    */
    typedef VALUETYPE value_type;
    
    /** init from given accessor
        @memo
    */
    BilinearInterpolatingAccessor(ACCESSOR a)
    : a_(a)
    {}
    
    /** Interpolate the data item at a non-integer distance.
        Ensure that no outside pixels are accessed if 
        (x, y) is near the image border (as long as
        0 <= x <= width-1, 0 <= y <= height-1).
        @memo
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
        0 <= y < height-1. It is slightly faster than #operator()#.
        @memo
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

/** Standard Constant Accessor.
    This accessor is normally used with all \Ref{ConstImageIterator}s. It
    simply encapsulates a call to the iterator's #operator*()# in its
    read function (as a const accessor, it does not have a write
    function).
    
    Include-File:
    \URL[vigra/accessor.hxx]{../include/vigra/accessor.hxx}
    
*/

template <class Iter1, class Acc1, class Iter2, class Acc2>
class MultiImageAccessor2
{
  public:
    typedef pair<typename Acc1::value_type, typename Acc2::value_type>
            value_type;
    
    MultiImageAccessor2(Iter1 i1, Acc1 a1, Iter2 i2, Acc2 a2)
    : i1_(i1), a1_(a1), i2_(i2), a2_(a2)
    {}

    /** read the current data item
        @memo
    */
    template <class DISTANCE>
    value_type operator()(DISTANCE const & d) const
    { 
        return std::make_pair(a1_(i1_, d), a2_(i2_, i.x, i.y)); 
    }
    
    /** read the data item at a distance
        @memo
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

#endif // VIGRA_ACCESSOR_HXX
