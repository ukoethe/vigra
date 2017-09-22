/************************************************************************/
/*                                                                      */
/*        Copyright 2014-2015 by Ullrich Koethe and Philip Schill       */
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
#ifndef VIGRA_FILTER_ITERATOR_HXX
#define VIGRA_FILTER_ITERATOR_HXX

#include <type_traits>
#include <iterator>

#include "iteratorfacade.hxx"

namespace vigra
{

namespace detail
{
    template <typename T>
    struct is_const_pointer
    {
        static bool const value = false;
    };

    template <typename T>
    struct is_const_pointer<T const *>
    {
        static bool const value = true;
    };

    template <typename ITER>
    struct is_const_iterator
    {
        typedef typename std::iterator_traits<ITER>::pointer pointer;
        static bool const value = is_const_pointer<pointer>::value;
    };
}


/********************************************************/
/*                                                      */
/*                     FilterIterator                   */
/*                                                      */
/********************************************************/

/**
 * @brief This iterator creates a view of another iterator and skips elements that
 *        do not fulfill a given predicate.
 *
 * The iterator is compatible to an STL forward iterator as defined in the C++ standard.
 *
 * @note The equality comparison only checks, whether the iterators point to the same place. The predicate is not checked.
 */
template <typename PREDICATE, typename ITER>
class FilterIterator
: public ForwardIteratorFacade<FilterIterator<PREDICATE, ITER>,
                               typename std::iterator_traits<ITER>::value_type,
                               detail::is_const_iterator<ITER>::value>
{
public:

    typedef PREDICATE Predicate;
    typedef ITER Iter;
    typedef typename std::iterator_traits<Iter>::value_type IterValueType;
    typedef FilterIterator<Predicate, Iter> SelfType;
    typedef ForwardIteratorFacade<SelfType,
                                  IterValueType,
                                  detail::is_const_iterator<ITER>::value> Parent;
    typedef typename Parent::value_type value_type;
    typedef typename Parent::reference reference;
    typedef reference const const_reference;

    /// Construct a filter iterator with the given predicate for
    /// a base iterator range \a iter to \a end.
    FilterIterator(Predicate pred, Iter iter, Iter end = Iter())
        :
        pred_(pred),
        iter_(iter),
        end_(end)
    {
        satisfy_predicate();
    }

private:

    void satisfy_predicate()
    {
        while (iter_ != end_ && !pred_(*iter_))
            ++iter_;
    }

    void increment()
    {
        ++iter_;
        satisfy_predicate();
    }

    reference dereference() const
    {
        return *iter_;
    }

    bool equal(FilterIterator const & other) const
    {
        return iter_ == other.iter_;
    }

    Predicate pred_;
    Iter iter_;
    Iter end_;

    friend class vigra::IteratorFacadeCoreAccess;

};

template <typename PREDICATE, typename ITER>
FilterIterator<PREDICATE, ITER>
make_filter_iterator(PREDICATE pred, ITER iter, ITER end = ITER())
{
    return FilterIterator<PREDICATE, ITER>(pred, iter, end);
}



} // namespace vigra



#endif
