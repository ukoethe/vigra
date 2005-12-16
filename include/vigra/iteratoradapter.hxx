/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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


#ifndef VIGRA_ITERATORADAPTER_HXX
#define VIGRA_ITERATORADAPTER_HXX

namespace vigra {

/********************************************************/
/*                                                      */
/*                    IteratorAdaptor                   */
/*                                                      */
/********************************************************/

/*! \brief Quckly create 1-dimensional iterator adapters.

    This class supports the easy creation of 1D iterator adpaters out
    of existing iterators. To use it, you must first implement a policy class
    that defines the iterator's behavior. The policy is used to
    instantiate the IteratorAdapter template, which thus automatically
    obtains all required functions of an STL-compatible iterator.
    General information on how this works can be found on the
    <a href="http://www.boost.org/libs/utility/iterator_adaptors.htm">Boost Iterator Adaptor</a>
    page, although there are some differences in the details of the
    boost and VIGRA implementations.
    Here is an example policy class that just exports the behaviour
    of the underlying iterator:

    \code
    template <class Iterator>
    class TrivialIteratorAdaptorPolicy
    {
      public:
        // the underlying iterator
        typedef Iterator                               BaseType;

        // the adaptor's value type
        typedef typename Iterator::value_type          value_type;

        // the adaptor's difference type (result of 'iter1 - iter2',
        //                                argument of 'iter[n]')
        typedef typename Iterator::difference_type     difference_type;

        // the adaptor's reference type (result of '*iter')
        typedef typename Iterator::reference           reference;

        // the adaptor's index_reference type (result of 'iter[n]')
        typedef typename Iterator::index_reference     index_reference;

        // the adaptor's pointer type (result of 'iter.operator->()')
        typedef typename Iterator::pointer             pointer;

        // the adaptor's iterator category
        typedef typename Iterator::iterator_category   iterator_category;

        // do some additional initialization in the adaptor's constructor
        static void initialize(BaseType & d) {}

        // called by '*iter', 'iter->'
        static reference dereference(BaseType const & d)
            { return *d; }

        // called by 'iter[n]'
        static index_reference dereference(BaseType d, difference_type n)
            { return d[n]; }

        // called by 'iter1 == iter2', 'iter1 != iter2'
        static bool equal(BaseType const & d1, BaseType const & d2)
            { return d1 == d2; }

        // called by 'iter1 < iter2', 'iter1 <= iter2', 'iter1 > iter2', 'iter1 >= iter2'
        static bool less(BaseType const & d1, BaseType const & d2)
            { return d1 < d2; }

        // called by 'iter1 - iter2'
        static difference_type difference(BaseType const & d1, BaseType const & d2)
            { return d1 - d2; }

        // called by '++iter', 'iter++'
        static void increment(BaseType & d)
            { ++d; }

        // called by '--iter', 'iter--'
        static void decrement(BaseType & d)
            { --d; }

        // called by 'iter += n', 'iter -= n'
        static void advance(BaseType & d, difference_type n)
            { d += n; }
    };
    \endcode

    This policy class is used like this:

    \code
    SomeIterator iter = ...;

    vigra::IteratorAdaptor<vigra::TrivialIteratorAdaptorPolicy<SomeIterator> > iter_adaptor(iter);
    \endcode

    By changing the definition of the policy members, a wide range of
    adaptor behaviors can be achieved. If the base iterator isn't a
    random access iterator, just drop the functions that cannot be implemented.
    This simply means that some adaptor functions may not be called,
    as one would expect from an iterator that doesn't support random access.
    Note also that the <TT>BaseType</TT> needs not be an iterator -
    it can be any type that contains the information necessary for the
    adaptor to do it's work.

    <b>\#include</b> "<a href="iteratoradapter_8hxx-source.html">vigra/iteratoradapter.hxx</a>"<br>
    Namespace: vigra

*/
template <class Policy>
class IteratorAdaptor
{
  public:

    typedef typename Policy::BaseType BaseType;
    typedef typename Policy::value_type        value_type;
    typedef typename Policy::difference_type   difference_type;
    typedef typename Policy::reference         reference;
    typedef typename Policy::index_reference   index_reference;
    typedef typename Policy::pointer           pointer;
    typedef typename Policy::iterator_category iterator_category;

    IteratorAdaptor()
    : adaptee_()
    {}

        /** Construct from an instance of the policy class' BaseType
            Note that the functions of the adaptor implement the
            interface of an random access iterator as defined in the
            C++ standard, so there is no need for explicit documentation.
        */
    explicit IteratorAdaptor(BaseType const & o)
    : adaptee_(o)
    {
        Policy::initialize(adaptee_);
    }

    IteratorAdaptor(IteratorAdaptor const & o)
    : adaptee_(o.adaptee_)
    {}

    IteratorAdaptor & operator=(BaseType const & o)
    {
        if(this != &o)
        {
            adaptee_ = o;
            Policy::initialize(adaptee_);
        }
        return *this;
    }

    IteratorAdaptor & operator=(IteratorAdaptor const & o)
    {
        if(this != &o)
            adaptee_ = o.adaptee_;
        return *this;
    }

    IteratorAdaptor & operator+=(difference_type d)
    {
        Policy::advance(adaptee_, d);
        return *this;
    }

    IteratorAdaptor operator+(difference_type d) const
    {
        return IteratorAdaptor(*this) += d;
    }

    IteratorAdaptor & operator-=(difference_type d)
    {
        Policy::advance(adaptee_, -d);
        return *this;
    }

    IteratorAdaptor operator-(difference_type d) const
    {
        return IteratorAdaptor(*this) -= d;
    }

    IteratorAdaptor & operator++()
    {
        Policy::increment(adaptee_);
        return *this;
    }

    IteratorAdaptor operator++(int)
    {
        IteratorAdaptor res(*this);
        Policy::increment(adaptee_);
        return res;
    }

    IteratorAdaptor & operator--()
    {
        Policy::decrement(adaptee_);
        return *this;
    }

    IteratorAdaptor operator--(int)
    {
        IteratorAdaptor res(*this);
        Policy::decrement(adaptee_);
        return res;
    }

    bool operator==(IteratorAdaptor const & o) const
    {
        return Policy::equal(adaptee_, o.adaptee_);
    }

    bool operator!=(IteratorAdaptor const & o) const
    {
        return !Policy::equal(adaptee_, o.adaptee_);
    }

    bool operator<(IteratorAdaptor const & o) const
    {
        return Policy::less(adaptee_, o.adaptee_);
    }

    bool operator<=(IteratorAdaptor const & o) const
    {
        return !Policy::less(o.adaptee_, adaptee_);
    }

    bool operator>(IteratorAdaptor const & o) const
    {
        return Policy::less(o.adaptee_, adaptee_);
    }

    bool operator>=(IteratorAdaptor const & o) const
    {
        return !Policy::less(adaptee_, o.adaptee_);
    }

    difference_type operator-(IteratorAdaptor const & o) const
    {
        return Policy::difference(adaptee_, o.adaptee_);
    }

    reference operator*() const
    {
        return Policy::dereference(adaptee_);
    }

    index_reference operator[](difference_type d) const
    {
        return Policy::dereference(adaptee_, d);
    }

    pointer operator->() const
    {
        return &Policy::dereference(adaptee_);
    }

  protected:

    BaseType adaptee_;
};

} // namespace vigra


#endif /* VIGRA_ITERATORADAPTER_HXX */
