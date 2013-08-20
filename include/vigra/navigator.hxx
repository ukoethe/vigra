/************************************************************************/
/*                                                                      */
/*                Copyright 2004 by Ullrich Koethe                      */
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

#ifndef VIGRA_NAVIGATOR_HXX
#define VIGRA_NAVIGATOR_HXX

namespace vigra {

/********************************************************/
/*                                                      */
/*                MultiArrayNavigator                   */
/*                                                      */
/********************************************************/

/** \brief A navigator that provides access to the 1D subranges of an
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

    <b>\#include</b> \<vigra/navigator.hxx\>

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
            //             d == {0, 1, 2}: iterate along {x, y, z}-axis respectively
            Navigator::iterator i = nav.begin(), end = nav.end();
            for(; i != end; ++i)
                // do something
        }
    }
    \endcode
*/
template <class MULTI_ITERATOR, unsigned int N>
class MultiArrayNavigator
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
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
    : base_type(i, shape, inner_dimension)
    {}

    MultiArrayNavigator(MULTI_ITERATOR const & i, shape_type const & start, shape_type const & stop, 
                        unsigned int inner_dimension)
    : base_type(i, start, stop, inner_dimension)
    {}

        /** Advance to next starting location.
         */
    void operator++()
    {
        base_type::operator++();
        if(this->point_[level-1] == this->stop_[level-1])
        {
            base_type::reset();
            ++this->point_[level];
            ++this->i_.template dim<level>();
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
        return this->point_[level] < this->stop_[level];
    }

        /** true if iterator is exhausted.
         */
    bool atEnd() const
    {
        return this->point_[level] >= this->stop_[level];
    }

  protected:
    void reset()
    {
        this->point_[level] = this->start_[level];
        this->i_.template dim<level>() -= (this->stop_[level] - this->start_[level]);
    }
};

template <class MULTI_ITERATOR>
class MultiArrayNavigator<MULTI_ITERATOR, 1>
{
  public:
    enum { level = 0 };
    typedef typename MULTI_ITERATOR::multi_difference_type shape_type;
    typedef typename MULTI_ITERATOR::iterator iterator;

    MultiArrayNavigator(MULTI_ITERATOR const & i, shape_type const & shape, unsigned int inner_dimension)
    : start_(), stop_(shape), point_(start_),
      inner_dimension_(inner_dimension),
      inner_shape_(stop_[inner_dimension] - start_[inner_dimension]),
      i_(i + start_)
    {
        if(stop_[inner_dimension] > start_[inner_dimension])
            stop_[inner_dimension] = start_[inner_dimension] + 1;
    }

    MultiArrayNavigator(MULTI_ITERATOR const & i, shape_type const & start, shape_type const & stop, 
                        unsigned int inner_dimension)
    : start_(start), stop_(stop), point_(start_),
      inner_dimension_(inner_dimension),
      inner_shape_(stop_[inner_dimension] - start_[inner_dimension]),
      i_(i + start_)
    {
        if(stop_[inner_dimension] > start_[inner_dimension])
            stop_[inner_dimension] = start_[inner_dimension] + 1;
    }

    void operator++()
    {
        ++point_[level];
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
        return begin() + inner_shape_;
    }

    bool hasMore() const
    {
        return point_[level] < stop_[level];
    }

    bool atEnd() const
    {
        return point_[level] >= stop_[level];
    }
    
    shape_type const & point() const
    {
        return point_;
    }

  protected:
    void reset()
    {
        point_[level] = start_[level];
        i_.template dim<level>() -= (stop_[level] - start_[level]);
   }

    shape_type start_, stop_, point_;
    unsigned int inner_dimension_, inner_shape_;
    MULTI_ITERATOR i_;
};

/********************************************************/
/*                                                      */
/*               MultiCoordinateNavigator               */
/*                                                      */
/********************************************************/

/** \brief A navigator that provides access to the 1D subranges of an
    n-dimensional range given by an nD shape.

    This class works similarly to \ref MultiArrayNavigator, but instead of a 
    1-dimensional iterator pair, it returns a pair of shapes whose difference
    specifies a 1-dimensional range along the desired dimension. That is, when
    the navigator refers to dimension <tt>d</tt>, the difference between
    <tt>end()</tt> and <tt>begin()</tt> is <tt>1</tt> along all dimensions
    except <tt>d</tt>.

    The template parameters specifies the dimension of the shape.

    <b>Usage:</b>

    <b>\#include</b> \<vigra/navigator.hxx\>

    Namespace: vigra

    \code
    typedef vigra::MultiArrayShape<3>::type Shape;
    typedef vigra::MultiArray<3, int>  Array;
    typedef vigra::MultiCoordinateNavigator<3> Navigator;

    Array a(Shape(X, Y, Z));

    for(int d=0; d<3; ++d)
    {
        // create Navigator for dimension d
        Navigator nav(a.shape(), d);

        // outer loop: move navigator to all starting points
        // of 1D subsets that run parallel to coordinate axis d
        for(; nav.hasMore(); ++nav)
        {
            // inner loop: linear iteration over current subset
            //             d == {0, 1, 2}: iterate along {x, y, z}-axis respectively
            Shape point = nav.begin(), end = nav.end();
            for(; point[d] != end[d]; ++point[d])
                a[point] = 5;
        }
    }
    \endcode
*/
template <unsigned int Dimensions, unsigned int N = Dimensions>
class MultiCoordinateNavigator
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
: public MultiCoordinateNavigator<Dimensions, N-1>
#endif
{
    typedef MultiCoordinateNavigator<Dimensions, N-1> base_type;

  public:
    enum { level = N-1 };

        /** The shape type for the given iterator type.
         */
    typedef typename MultiArrayShape<Dimensions>::type value_type;


        /** Construct navigator for multi-dimensional iterator <TT>i</TT>, array shape <TT>shape</TT>
            and inner loop dimension <TT>inner_dimension</TT>.
         */
    MultiCoordinateNavigator(value_type const & shape, unsigned int inner_dimension)
    : base_type(shape, inner_dimension)
    {
        this->end_[level] = (this->inner_dimension_ == level)
                                 ? 1
                                 : this->shape_[level];
    }

        /** Advance to next starting location.
         */
    void operator++()
    {
        base_type::operator++();
        if(base_type::atEnd() && this->i_[level] < this->end_[level])
        {
            ++this->i_[level];
            if(this->i_[level] < this->end_[level])
                base_type::reset();
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
        return this->inner_dimension_ == level 
                     ? base_type::hasMore() 
                     : this->i_[level] < this->end_[level];
    }

        /** true if iterator is exhausted.
         */
    bool atEnd() const
    {
        return !hasMore();
        // return this->inner_dimension_ == level 
                     // ? base_type::atEnd()
                     // : !(this->i_[level] < this->end_[level]);
    }

  protected:
    void reset()
    {
        this->i_[level] = 0;
        this->end_[level] = (this->inner_dimension_ == level)
                                 ? 1
                                 : this->shape_[level];
        base_type::reset();
    }
};

template <unsigned int Dimensions>
class MultiCoordinateNavigator<Dimensions, 1>
{
  public:
    enum { level = 0 };
    typedef typename MultiArrayShape<Dimensions>::type value_type;
 
    MultiCoordinateNavigator(value_type const & shape, unsigned int inner_dimension)
    : shape_(shape),
      inner_dimension_(inner_dimension)
    {
        end_[level] = (inner_dimension_ == level)
                         ? 1
                         : shape_[level];
    }

    void operator++()
    {
        ++i_[level];
    }

    void operator++(int)
    {
        ++*this;
    }

    value_type const & begin() const
    {
        return i_;
    }

    value_type end() const
    {
        value_type res = i_ + value_type(MultiArrayIndex(1));
        res[inner_dimension_] = shape_[inner_dimension_];
        return res;
    }

    bool hasMore() const
    {
        return i_[level] < end_[level];
    }

    bool atEnd() const
    {
      return !hasMore();
    }

  protected:
    void reset()
    {
        i_[level] = 0;
        end_[level] = (inner_dimension_ == level)
                         ? 1
                         : shape_[level];
    }

    value_type shape_, i_, end_;
    unsigned int inner_dimension_;
};

} // namespace vigra

#endif /* VIGRA_NAVIGATOR_HXX */
