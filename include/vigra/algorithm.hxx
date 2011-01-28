/************************************************************************/
/*                                                                      */
/*               Copyright 2010-2011 by Ullrich Koethe                  */
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

#ifndef VIGRA_ALGORITHM_HXX
#define VIGRA_ALGORITHM_HXX

#include "numerictraits.hxx"
#include <algorithm>
#include <functional>
#include <iterator>

namespace vigra {

/** \addtogroup MathFunctions
*/
//@{
    /*! Find the minimum element in a sequence.
    
        The function returns the iterator refering to the minimum element.
        This is identical to the function <tt>std::min_element()</tt>.
        
        <b>Required Interface:</b>
        
        \code
        Iterator is a standard forward iterator.
        
        bool f = *first < NumericTraits<typename std::iterator_traits<Iterator>::value_type>::max();
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator>
Iterator argMin(Iterator first, Iterator last)
{
    if(first == last)
        return last;
    Iterator best = first;
    for(++first; first != last; ++first)
        if(*first < *best)
            best = first;
    return best;
}

    /*! Find the maximum element in a sequence.
    
        The function returns the iterator refering to the maximum element.
        This is identical to the function <tt>std::max_element()</tt>.
        
        <b>Required Interface:</b>
        
        \code
        Iterator is a standard forward iterator.
        
        bool f = NumericTraits<typename std::iterator_traits<Iterator>::value_type>::min() < *first;
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator>
Iterator argMax(Iterator first, Iterator last)
{
    if(first == last)
        return last;
    Iterator best = first;
    for(++first; first != last; ++first)
        if(*best < *first)
            best = first;
    return best;
}

    /*! Find the minimum element in a sequence conforming to a condition.
    
        The function returns the iterator refering to the minimum element,
        where only elements conforming to the condition (i.e. where 
        <tt>condition(*iterator)</tt> evaluates to <tt>true</tt>) are considered.
        If no element conforms to the condition, or the sequence is empty,
        the end iterator \a last is returned.
        
        <b>Required Interface:</b>
        
        \code
        Iterator is a standard forward iterator.
        
        bool c = condition(*first);
        
        bool f = *first < NumericTraits<typename std::iterator_traits<Iterator>::value_type>::max();
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator, class UnaryFunctor>
Iterator argMinIf(Iterator first, Iterator last, UnaryFunctor condition)
{
    for(; first != last; ++first)
        if(condition(*first))
            break;
    if(first == last)
        return last;
    Iterator best = first;
    for(++first; first != last; ++first)
        if(condition(*first) && *first < *best)
            best = first;
    return best;
}

    /*! Find the maximum element in a sequence conforming to a condition.
    
        The function returns the iterator refering to the maximum element,
        where only elements conforming to the condition (i.e. where 
        <tt>condition(*iterator)</tt> evaluates to <tt>true</tt>) are considered.
        If no element conforms to the condition, or the sequence is empty,
        the end iterator \a last is returned.
        
        <b>Required Interface:</b>
        
        \code
        Iterator is a standard forward iterator.
        
        bool c = condition(*first);
        
        bool f = NumericTraits<typename std::iterator_traits<Iterator>::value_type>::min() < *first;
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator, class UnaryFunctor>
Iterator argMaxIf(Iterator first, Iterator last, UnaryFunctor condition)
{
    for(; first != last; ++first)
        if(condition(*first))
            break;
    if(first == last)
        return last;
    Iterator best = first;
    for(++first; first != last; ++first)
        if(condition(*first) && *best < *first)
            best = first;
    return best;
}

    /*! Fill an array with a sequence of numbers.
    
        The sequence starts at \a start and is incremented with \a step. Default start
        and stepsize are 0 and 1 respectively.
        
        <b> Declaration:</b>

        \code
        namespace vigra {
            template <class Iterator, class Value>
            void linearSequence(Iterator first, Iterator last, 
                          Value const & start = 0, Value const & step = 1);
        }
        \endcode
        
        <b>Required Interface:</b>
        
        \code
        Iterator is a standard forward iterator.
        
        *first = start;
        start += step;
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator, class Value>
void linearSequence(Iterator first, Iterator last, Value start, Value step)
{
    for(; first != last; ++first, start += step)
        *first = start;
}

template <class Iterator, class Value>
void linearSequence(Iterator first, Iterator last, Value start)
{
    for(; first != last; ++first, ++start)
        *first = start;
}

template <class Iterator>
void linearSequence(Iterator first, Iterator last)
{
    typedef typename std::iterator_traits<Iterator>::value_type Value;
    
    linearSequence(first, last, NumericTraits<Value>::zero());
}

namespace detail {

template <class Iterator, class Compare>
struct IndexCompare
{
    Iterator i_;
    Compare c_;
    
    IndexCompare(Iterator i, Compare c)
    : i_(i),
      c_(c)
    {}

    template <class Index>
    bool operator()(Index const & l, Index const & r) const
    {
        return c_(i_[l], i_[r]);
    }
};

} // namespace detail

    /*! Return the index permutation that would sort the input array.
    
        To actually sort an array according to the ordering thus determined, use 
        \ref applyPermutation().
        
        <b> Declarations:</b>

        \code
        namespace vigra {
            // compare using std::less
            template <class Iterator, class IndexIterator>
            void indexSort(Iterator first, Iterator last, IndexIterator index_first);

            // compare using functor Compare
            template <class Iterator, class IndexIterator, class Compare>
            void indexSort(Iterator first, Iterator last, IndexIterator index_first, Compare compare);
        }
        \endcode
        
        <b>Required Interface:</b>
        
        \code
        Iterator and IndexIterators are random access iterators.
        
        bool res = compare(first[*index_first], first[*index_first]);
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator, class IndexIterator, class Compare>
void indexSort(Iterator first, Iterator last, IndexIterator index_first, Compare c)
{
    int size = last - first;
    linearSequence(index_first, index_first+size);
    std::sort(index_first, index_first+size, 
              detail::IndexCompare<Iterator, Compare>(first, c));
}

template <class Iterator, class IndexIterator>
void indexSort(Iterator first, Iterator last, IndexIterator index_first)
{
    typedef typename std::iterator_traits<Iterator>::value_type Value;
    indexSort(first, last, index_first, std::less<Value>());
}

    /*! Sort an array according to the given index permutation.
    
        The iterators \a in and \a out may not refer to the same array, as
        this would overwrite the input prematurely.
        
        <b> Declaration:</b>

        \code
        namespace vigra {
            template <class IndexIterator, class InIterator, class OutIterator>
            void applyPermutation(IndexIterator index_first, IndexIterator index_last, 
                                  InIterator in, OutIterator out);
        }
        \endcode
        
        <b>Required Interface:</b>
        
        \code
        OutIterator and IndexIterators are forward iterators.
        InIterator is a random access iterator.
        
        *out = in[*index_first];
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class IndexIterator, class InIterator, class OutIterator>
void applyPermutation(IndexIterator index_first, IndexIterator index_last, 
                      InIterator in, OutIterator out)
{
    for(; index_first != index_last; ++index_first, ++out)
        *out = in[*index_first];
}


    /*! Compute the inverse of a given permutation.
    
        This is just another name for \ref indexSort(), referring to
        another semantics.
        
        <b> Declaration:</b>

        \code
        namespace vigra {
            template <class InIterator, class OutIterator>
            void inversePermutation(InIterator first, InIterator last, 
                                    OutIterator out);
        }
        \endcode
        
        <b>Required Interface:</b>
        
        \code
        InIterator and OutIterator are random access iterators.
        
        *out = in[*index_first];
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class InIterator, class OutIterator>
void inversePermutation(InIterator first, InIterator last, 
                        OutIterator out)
{
    indexSort(first, last, out);
}

//@}

} // namespace vigra

#endif /* VIGRA_ALGORITHM_HXX */
