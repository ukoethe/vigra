/************************************************************************/
/*                                                                      */
/*                 Copyright 2015 by Thorsten Beier                     */
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

#ifndef VIGRA_COUNTING_ITERATOR_HXX
#define VIGRA_COUNTING_ITERATOR_HXX


#include <cmath>
#include <iterator>
#include "error.hxx"
#include "tinyvector.hxx"

namespace vigra {

    /** \brief Iterator that counts upwards or downwards with a given step size.

        This iterator replicates the functionality of Python's
        well-known range-function. It is especially convenient in
        range-based for-loops. <tt>CountingIterator</tt> also works for
        floating-point counting.

        <b>Usage:</b>

        <b>\#include</b> \<vigra/counting_iterator.hxx\><br>
        Namespace: vigra

        You will normally construct instances of this iterator with
        one of the <tt>range()</tt> factory functions. There are three versions
        of this function <tt>range(end)</tt>, <tt>range(begin, end)</tt>, and
        <tt>range(begin, end, step)</tt>.
        \code
        // count upwards from 0 to 4
        for(int i: range(5))
            std::cout << i << " ";  // prints '0 1 2 3 4'

        // count upwards from 4 to 7
        for(int i: range(4, 8))
            std::cout << i << " ";  // prints '4 5 6 7'

        // count upwards from 0 to 9 with step 3
        for(int i: range(0, 9, 3))
            std::cout << i << " ";  // prints '0 3 6'

        // likewise (note upper bound)
        for(int i: range(0, 7, 3))
            std::cout << i << " ";  // prints '0 3 6'

        // count downwards from 4 to 1 with step -1
        for(int i: range(4, 0))
            std::cout << i << " ";  // prints '4 3 2 1'

        // count downwards from 8 to 2 with step -2
        for(int i: range(8, 0, -2))
            std::cout << i << " ";  // prints '8 6 4 2'
        \endcode

        Alternatively, you can create a traditional random-access iterator pair.
        The end iterator can be conveniently constructed by the begin iterator's
        <tt>end()</tt> function:
        \code
        auto iter = range(5),
             end  = iter.end();
        std::cout << std::accumulate(iter, end, 0) << std::endl; // prints '10'
        \endcode

        <tt>range()</tt> and <tt>CountingIterator</tt> also work for floating-point
        arguments. However, it is not recommended to use this variant in
        range-based for-loops or STL algorithms because the condition
        <tt>iter != end</tt> may never become false due to rounding errors, leading
        to an infinite loop. It is safer to use the condition <tt>iter < end</tt>
        or <tt>iter <= end</tt> to terminate the loop:
        \code
        auto iter = range(1.0, 1.6, 0.1),
             end  = iter.end();
        for(; iter <= end; ++iter)
            std::cout << *iter << " ";  // prints '1 1.1 1.2 1.3 1.4 1.5 1.6'
        \endcode

        Note that the termination condition is still <tt>iter <= end</tt>, even
        when the iterator counts downwards:
        \code
        auto iter = range(1.6, 1.0, -0.1),
             end  = iter.end();
        for(; iter <= end; ++iter)
            std::cout << *iter << " ";  // prints '1.6 1.5 1.4 1.3 1.2 1.1 1'
        \endcode
    */
template<class T = ptrdiff_t>
class CountingIterator
: public std::iterator<std::random_access_iterator_tag,
                       T, ptrdiff_t, T const *, T>
{
  public:
    CountingIterator()
    : begin_(0)
    , end_(0)
    , step_(1)
    {}

    CountingIterator(T begin, T end)
    : begin_(begin)
    , end_(end)
    , step_(1)
    {
        vigra_precondition(begin <= end,
            "CountingIterator(): begin must be less or equal to end.");
    }

    CountingIterator(T begin, T end, T step)
    : begin_(begin)
    , end_(end)
    , step_(step)
    {
        vigra_precondition(step != 0,
            "CountingIterator(): step must be non-zero.");
        vigra_precondition((step > 0 && begin <= end) || (step < 0 && begin >= end),
            "CountingIterator(): sign mismatch between step and (end-begin).");
    }

    CountingIterator(CountingIterator const & other, ReverseCopyTag)
    : begin_(other.end_)
    , end_(other.begin_)
    , step_(-other.step_)
    {}

  public:

    CountingIterator begin() const
    {
        return *this;
    }

    CountingIterator end() const
    {
        // since the range-based for-loop checks "iter != end",
        // (end - begin) must be a multiple of step to avoid an infinite loop
        T end = begin_ + step_*distance(begin_, end_, step_);
        return CountingIterator(end, end, step_);
    }

    bool empty() const
    {
        return (end_ - begin_)*step_ <= 0;
    }

    CountingIterator& operator++()    {begin_ += step_; return *this;} // prefix++
    CountingIterator  operator++(int) {CountingIterator tmp(*this); ++(*this); return tmp;} // postfix++
    CountingIterator& operator--()    {begin_ -= step_; return *this;} // prefix--
    CountingIterator  operator--(int) {CountingIterator tmp(*this); --(*this); return tmp;} // postfix--

    CountingIterator& operator+=(ptrdiff_t n)
    {
        begin_ += n*step_;
        return *this;
    }

    CountingIterator operator+(ptrdiff_t n) const
    {
        return CountingIterator(*this) += n;
    }

    CountingIterator& operator-=(ptrdiff_t n)
    {
        begin_ -= n*step_;
        return *this;
    }

    CountingIterator operator-(ptrdiff_t n) const
    {
        return CountingIterator(*this) -= n;
    }

    ptrdiff_t operator-(const CountingIterator& other) const
    {
        return distance(other.begin_, begin_, step_);
    }

    bool operator<(CountingIterator const & other) const
    {
        // NOTE: the more efficient '(other.begin_ - begin_)*step_ > 0'
        //       fails for unsigned arguments
        return step_ > 0
                ? begin_ < other.begin_
                : begin_ > other.begin_;
    }

    bool operator<=(CountingIterator const & other) const
    {
        return step_ > 0
                ? begin_ <= other.begin_
                : begin_ >= other.begin_;
    }

    bool operator>(CountingIterator const & other) const
    {
        return step_ > 0
                ? begin_ > other.begin_
                : begin_ < other.begin_;
    }

    bool operator>=(CountingIterator const & other) const
    {
        return step_ > 0
                ? begin_ >= other.begin_
                : begin_ <= other.begin_;
    }

    bool operator==(const CountingIterator& other) const
    {
        return begin_ == other.begin_;
    }

    bool operator!=(const CountingIterator& other) const
    {
        return begin_ != other.begin_;
    }

    T operator[](ptrdiff_t n) const {
        return begin_ + n*step_;
    }

    T operator*() const {
        return  begin_;
    }

    T const * operator->() const{
        return &begin_;
    }

  private:
    static ptrdiff_t distance(T from, T to, T step)
    {
        const double diff = (double(to) - double(from)) / double(step);
        return diff > 0.0
                 ? (ptrdiff_t)std::ceil(diff)
                 : (ptrdiff_t)std::floor(diff);
    }

    T begin_, end_, step_;
};


template <class T1, class T2, class T3>
inline CountingIterator<T1>
range(T1 begin, T2 end, T3 step)
{
    return CountingIterator<T1>(begin, end, step);
}

template <class T1, class T2>
inline CountingIterator<T1>
range(T1 begin, T2 end)
{
    return CountingIterator<T1>(begin, end, 1);
}

template <class T>
inline CountingIterator<T>
range(T end)
{
    return CountingIterator<T>(0, end, 1);
}

} // namespace vigra

#endif
