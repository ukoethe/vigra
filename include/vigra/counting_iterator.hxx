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
#include <limits>
#include <type_traits>
#include "error.hxx"
#include "tinyvector.hxx"

namespace vigra {

namespace detail  {

template <class T, bool is_float=false>
struct CountingIteratorCompare
{
    // use exact comparison for integer counting
    static bool equal(T left, T right, T /* step */)
    {
        return left == right;
    }
    static bool not_equal(T left, T right, T /* step */)
    {
        return left != right;
    }
    static bool less(T left, T right, T step)
    {
        // NOTE: the more efficient '(right - left)*step > 0'
        //       fails for unsigned arguments
        return step > 0
                ? left < right
                : left > right;
    }
    static bool less_equal(T left, T right, T step)
    {
        return step > 0
                ? left <= right
                : left >= right;
    }
    static bool greater(T left, T right, T step)
    {
        return step > 0
                ? left > right
                : left < right;
    }
    static bool greater_equal(T left, T right, T step)
    {
        return step > 0
                ? left >= right
                : left <= right;
    }
    // integer counting: if the raw distance is not divisible by step,
    // we must round upwards
    static std::ptrdiff_t distance(T from, T to, T step)
    {
        const double diff = (double(to) - double(from)) / double(step);
        return diff > 0.0
                 ? (std::ptrdiff_t)std::ceil(diff)
                 : (std::ptrdiff_t)std::floor(diff);
    }
};

template <class T>
struct CountingIteratorCompare<T, true>
{
    typedef std::numeric_limits<T> limit;

    // use comparison with tolerance for floating-point counting
    // (the natural epsilon is 0.5*step)
    static bool equal(T left, T right, T step)
    {
        return std::fabs(right-left) <= 0.5*std::fabs(step);
    }
    static bool not_equal(T left, T right, T step)
    {
        return std::fabs(right-left) > 0.5*std::fabs(step);
    }
    static bool less(T left, T right, T step)
    {
        return step > 0.0
                ? right - left > 0.5*step
                : right - left < 0.5*step;
    }
    static bool less_equal(T left, T right, T step)
    {
        return step > 0.0
                ? left - right < 0.5*step
                : left - right > 0.5*step;
    }
    static bool greater(T left, T right, T step)
    {
        return step > 0.0
                ? left - right > 0.5*step
                : left - right < 0.5*step;
    }
    static bool greater_equal(T left, T right, T step)
    {
        return step > 0.0
                ? right - left < 0.5*step
                : right - left > 0.5*step;
    }
    // floating-point counting: if the raw distance is not divisible by step,
    // we round to nearest if the difference is small, otherwise upwards
    static std::ptrdiff_t distance(T from, T to, T step)
    {
        const double diff = (double(to) - double(from)) / double(step);
        return diff > 0.0
                 ? (std::ptrdiff_t)std::ceil(diff*(1.0-2.0*limit::epsilon()))
                 : (std::ptrdiff_t)std::floor(diff*(1.0-2.0*limit::epsilon()));
    }
};

} // namespace detail

/** \addtogroup RangesAndPoints
*/
//@{

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
        arguments. As in the integer case, the upper bound is excluded from the range
        if it can be reached by an integer multiple of the step (within machine
        epsilon):
        \code
        for(auto i: range(1.0, 1.6, 0.1))  // 1.6 is excluded
            std::cout << i << " ";         // prints '1 1.1 1.2 1.3 1.4 1.5'

        for(auto i: range(1.0, 1.61, 0.1)) // 1.6 is included
            std::cout << i << " ";         // prints '1 1.1 1.2 1.3 1.4 1.5 1.6'
        \endcode

        If you use an iterator pair, you can make clear which behavior you want
        by using either <tt>iter < end</tt> or <tt>iter <= end</tt> to terminate
        the loop:
        \code
        auto iter = range(1.0, 1.6, 0.1),
             end  = iter.end();
        for(; iter < end; ++iter)       // exclude upper bound
            std::cout << *iter << " ";  // prints '1 1.1 1.2 1.3 1.4 1.5'

        iter = range(1.0, 1.6, 0.1);
        for(; iter <= end; ++iter)      // include upper bound
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
template<class T = std::ptrdiff_t>
class CountingIterator
: public std::iterator<std::random_access_iterator_tag,
                       T, std::ptrdiff_t, T const *, T>
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
        // (end - begin) must be a multiple of step to avoid infinite loops
        T end = begin_ + step_*Compare::distance(begin_, end_, step_);
        return CountingIterator(end, end, step_);
    }

    bool empty() const
    {
        return Compare::greater_equal(begin_, end_, step_);
    }

    CountingIterator& operator++()    {begin_ += step_; return *this;} // prefix++
    CountingIterator  operator++(int) {CountingIterator tmp(*this); ++(*this); return tmp;} // postfix++
    CountingIterator& operator--()    {begin_ -= step_; return *this;} // prefix--
    CountingIterator  operator--(int) {CountingIterator tmp(*this); --(*this); return tmp;} // postfix--

    CountingIterator& operator+=(std::ptrdiff_t n)
    {
        begin_ += n*step_;
        return *this;
    }

    CountingIterator operator+(std::ptrdiff_t n) const
    {
        return CountingIterator(*this) += n;
    }

    CountingIterator& operator-=(std::ptrdiff_t n)
    {
        begin_ -= n*step_;
        return *this;
    }

    CountingIterator operator-(std::ptrdiff_t n) const
    {
        return CountingIterator(*this) -= n;
    }

    std::ptrdiff_t operator-(const CountingIterator& other) const
    {
        return Compare::distance(other.begin_, begin_, step_);
    }

    bool operator<(CountingIterator const & other) const
    {
        return Compare::less(begin_, other.begin_, step_);
    }

    bool operator<=(CountingIterator const & other) const
    {
        return Compare::less_equal(begin_, other.begin_, step_);
    }

    bool operator>(CountingIterator const & other) const
    {
        return Compare::greater(begin_, other.begin_, step_);
    }

    bool operator>=(CountingIterator const & other) const
    {
        return Compare::greater_equal(begin_, other.begin_, step_);
    }

    bool operator==(const CountingIterator& other) const
    {
        return Compare::equal(begin_, other.begin_, step_);
    }

    bool operator!=(const CountingIterator& other) const
    {
        return Compare::not_equal(begin_, other.begin_, step_);
    }

    T operator[](std::ptrdiff_t n) const {
        return begin_ + n*step_;
    }

    T operator*() const {
        return  begin_;
    }

    T const * operator->() const{
        return &begin_;
    }

  private:

    typedef detail::CountingIteratorCompare<T, std::is_floating_point<T>::value> Compare;

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

//@}

} // namespace vigra

#endif
