/************************************************************************/
/*                                                                      */
/*                Copyright 2004 by Ullrich Koethe                      */
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

#ifndef VIGRA_NAVIGATOR_HXX
#define VIGRA_NAVIGATOR_HXX

namespace vigra {

/********************************************************/
/*                                                      */
/*                MultiArrayNavigator                   */
/*                                                      */
/********************************************************/

/** \brief A navigator that provides acces to the 1D subranges of an 
    n-dimensional range given by a \ref vigra::MultiIterator and an nD shape.

    Normally, the innermost loop of an iteration extends over the innermost
    dimension of a given array. Sometimes, however, it is necessary to have
    some other dimension in the inner loop. For example, instead of iterating over 
    the rows, the inner loop should extend over the columns. The class MultiArrayNavigator
    encapsulates the necessary functionality. Given an arbitrary dimensional 
    array (represented by a vigra::MultiIterator/shape pair), and the desired
    inner loop dimension <TT>d</TT>, it moves the encapsulated iterator to all possible
    starting points of 1D subsets along the given dimension (e.g. all columns). By calling 
    <TT>begin()</TT> and <TT>end()</TT>, one can then obtain an STL-compatible 1-dimensional 
    iterator for the current subset.

    The template parameters specify the embedded iterator type and its dimension.

    <b>Usage:</b>

    <b>\#include</b> "<a href="navigator_8hxx-source.html">vigra/navigator.hxx</a>"

    Namespace: vigra

    \code
    typedef vigra::MultiArray<3, int>  Array;

    Array a(Array::size_type(X, Y, Z));

    typedef vigra::MultiArrayNavigator<Array::traverser, 3> Navigator;

    for(int d=0; d<3; ++d)
    {
        // create Navigator for dimension d
        Navigator nav(a.traverser_begin(), a.shape(), d);

        // outer loop: move navigator to all starting points
        // of 1D subsets that run parallel to coordinate axis d
        for(; nav.hasMore(); ++nav)
        {
            // inner loop: linear iteration over current subset
            //             d == {0, 1, 2}: interate along {x, y, z}-axis respectively
             Navigator::iterator i = nav.begin(), end = nav.end();
            for(; i != end; ++i)
                // do something
        }
    }
    \endcode
*/
template <class MULTI_ITERATOR, unsigned int N>
class MultiArrayNavigator 
#ifndef DOXYGEN  // docygen doesn't understand this inheritance
: public MultiArrayNavigator<MULTI_ITERATOR, N-1>
#endif
{
    typedef MultiArrayNavigator<MULTI_ITERATOR, N-1> base_type;

  public:
    enum { level = N-1 };

        /** The required shape type for the given iterator type.
         */
    typedef typename MULTI_ITERATOR::multi_difference_type shape_type;

        /** The iterator type for the inner loop (result of begin() and end()).
         */
    typedef typename MULTI_ITERATOR::iterator iterator;
    
        /** Construct navigator for multi-dimensional iterator <TT>i</TT>, array shape <TT>shape</TT>
            and inner loop dimension <TT>inner_dimension</TT>.
         */
    MultiArrayNavigator(MULTI_ITERATOR const & i, shape_type const & shape, unsigned int inner_dimension)
    : base_type(i, shape, inner_dimension),
      i_(i), 
      end_(i)
    {
        if(inner_dimension != level)
            end_.template dim<level>() += shape[level];
    }
    
        /** Advance to next starting location.
         */
    void operator++()
    {
        base_type::operator++();
        if(base_type::atEnd() && i_ < end_) // this tests implicitly inner_dimension_ != level
        {
            ++i_.template dim<level>();
            if(i_ < end_)
                base_type::reset(i_);
        }
    }
    
        /** Advance to next starting location.
         */
    void operator++(int)
    {
        ++*this;
    }
    
        /** true if there are more elements.
         */
    bool hasMore() const
    {
        return this->inner_dimension_ == level ?
                    base_type::hasMore() :
                    i_ < end_; 
    }
    
        /** true if iterator is exhausted.
         */
    bool atEnd() const
    {
        return this->inner_dimension_ == level ?
                    base_type::atEnd() :
	            !( i_ < end_);
    }
    
  protected:
    void reset(MULTI_ITERATOR const & i)
    {
        end_ = i_ = i;
        if(this->inner_dimension_ != level)
            end_.template dim<level>() += this->shape_[level];
        base_type::reset(i);
    }
    
    MULTI_ITERATOR i_, end_;
};

template <class MULTI_ITERATOR>
class MultiArrayNavigator<MULTI_ITERATOR, 1>
{
  public:
    enum { level = 0 };
    typedef typename MULTI_ITERATOR::multi_difference_type shape_type;
    typedef typename MULTI_ITERATOR::iterator iterator;

    MultiArrayNavigator(MULTI_ITERATOR const & i, shape_type const & shape, unsigned int inner_dimension)
    : shape_(shape),
      inner_dimension_(inner_dimension),
      i_(i), 
      end_(i)
    {
        if(inner_dimension != level)
            end_.template dim<level>() += shape[level];
    }

    void operator++()
    {
        ++i_.template dim<level>();
    }
    
    void operator++(int)
    {
        ++*this;
    }
    
    iterator begin() const
    {
        return i_.iteratorForDimension(inner_dimension_);
    }
    
    iterator end() const
    {
        return begin() + shape_[inner_dimension_];
    }
    
    bool hasMore() const
    {
        return i_ < end_; 
    }
    
    bool atEnd() const
    {
      return !( i_ < end_);
    }
    
  protected:
    void reset(MULTI_ITERATOR const & i)
    {
        end_ = i_ = i;
        if(inner_dimension_ != level)
            end_.template dim<level>() += shape_[level];
    }
    
    shape_type shape_;
    unsigned int inner_dimension_;
    MULTI_ITERATOR i_, end_;
};

} // namespace vigra

#endif /* VIGRA_NAVIGATOR_HXX */
