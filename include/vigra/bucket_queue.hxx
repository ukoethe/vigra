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


#ifndef VIGRA_BUCKET_QUEUE_HXX
#define VIGRA_BUCKET_QUEUE_HXX

#include "config.hxx"
#include "error.hxx"
#include "array_vector.hxx"
#include <queue>

namespace vigra {

/** \brief Priority queue implemented using bucket sort.

    This template implements functionality similar to <tt><a href="http://www.sgi.com/tech/stl/priority_queue.html">std::priority_queue</a></tt>,
    but uses a more efficient algorithm based on bucket sort. It can be used
    when all priorities are positive integers in a given range (typically, 0...255).
    By default, <tt>BucketQueue\<ValueType\></tt> sorts the elements in descending order,
    i.e. like in <tt>std::priority_queue</tt> the largest element has highest priority.
    An ascending queue can be specified as <tt>BucketQueue\<ValueType, true\></tt>.
    Elements with equal priorities are returned in a first-in first-out fashion.
    
    The main difference to <tt>std::priority_queue</tt> is the function <tt>push</tt>
    which explicitly takes the priority of the element to be added as a second argument.
    This allows optimization of <tt>ValueType</tt>: since the bucket uniquely
    determines an element's priority, there is no need for <tt>ValueType</tt> to
    store redundant priority information. If compatibility to <tt>std::priority_queue</tt>
    is more important, use \ref vigra::MappedBucketQueue.

    <b>\#include</b> \<vigra/bucket_queue.hxx\><br>
    Namespace: vigra
*/
template <class ValueType,
          bool Ascending = false>  // std::priority_queue is descending
class BucketQueue
{
    ArrayVector<std::queue<ValueType> > buckets_;
    std::size_t size_;
    std::ptrdiff_t top_;
    
  public:
  
    typedef ValueType value_type;
    typedef ValueType & reference;
    typedef ValueType const & const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t priority_type;
    
        /** \brief Create bucket queue with \arg bucket_count entries.
            Priorities must be integers in the range <tt>[0, ..., bucket_count-1]</tt>.
        */
    BucketQueue(size_type bucket_count = 256)
    : buckets_(bucket_count),
      size_(0), top_(0)
    {}
    
        /** \brief Number of elements in this queue.
        */
    size_type size() const
    {
        return size_;
    }
    
        /** \brief Quene contains no elements.
             Equivalent to <tt>size() == 0</tt>.
        */
    bool empty() const
    {
        return size() == 0;
    }
    
        /** \brief Maximum index (i.e. priority) allowed in this queue.
             Equivalent to <tt>bucket_count - 1</tt>.
        */
    priority_type maxIndex() const
    {
        return (priority_type)buckets_.size() - 1;
    }
    
        /** \brief Priority of the current top element.
        */
    priority_type topPriority() const
    {
        return top_;
    }
    
        /** \brief The current top element.
        */
    const_reference top() const
    {
        
        return buckets_[top_].front();
    }

        /** \brief Remove the current top element.
        */
    void pop()
    {
        --size_;
        buckets_[top_].pop();
        
        while(top_ > 0 && buckets_[top_].size() == 0)
            --top_;
    }
    
        /** \brief Insert new element \arg v with given \arg priority.
        */
    void push(value_type const & v, priority_type priority)
    {
        ++size_;
        buckets_[priority].push(v);
        
        if(priority > top_)
            top_ = priority;
    }
};

template <class ValueType> 
class BucketQueue<ValueType, true> // ascending queue
{
    ArrayVector<std::queue<ValueType> > buckets_;
    std::size_t size_;
    std::ptrdiff_t top_;
    
  public:
  
    typedef ValueType value_type;
    typedef ValueType & reference;
    typedef ValueType const & const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t priority_type;
    
    BucketQueue(size_type bucket_count = 256)
    : buckets_(bucket_count),
      size_(0), top_((priority_type)bucket_count)
    {}
    
    size_type size() const
    {
        return size_;
    }
    
    bool empty() const
    {
        return size() == 0;
    }
    
    priority_type maxIndex() const
    {
        return (priority_type)buckets_.size() - 1;
    }
    
    priority_type topPriority() const
    {
        return top_;
    }
    
    const_reference top() const
    {
        
        return buckets_[top_].front();
    }

    void pop()
    {
        --size_;
        buckets_[top_].pop();
        
        while(top_ < (priority_type)buckets_.size() && buckets_[top_].size() == 0)
            ++top_;
    }
    
    void push(value_type const & v, priority_type priority)
    {
        ++size_;
        buckets_[priority].push(v);
        
        if(priority < top_)
            top_ = priority;
    }
};

/** \brief Priority queue implemented using bucket sort (STL compatible).

    This template is compatible to <tt><a href="http://www.sgi.com/tech/stl/priority_queue.html">std::priority_queue</a></tt>,
    but uses a more efficient algorithm based on bucket sort. It us used
    like \ref vigra::BucketQueue, but has an additional <tt>PriorityFunctor</tt>
    which extracts the priority value of an element of type <tt>ValueType</tt>.
    Thus functor is called within <tt>push</tt> so that it does not need an
    extra argument.

    <b>\#include</b> \<vigra/bucket_queue.hxx\><br>
    Namespace: vigra
*/
template <class ValueType,
          class PriorityFunctor, 
          bool Ascending = false> 
class MappedBucketQueue
: public BucketQueue<ValueType, Ascending>
{
    PriorityFunctor get_priority_;
    
  public:

    typedef BucketQueue<ValueType, Ascending> BaseType;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::reference reference;
    typedef typename BaseType::const_reference const_reference;
    typedef typename BaseType::size_type size_type;
    typedef typename BaseType::priority_type priority_type;
    
        /** \brief Create a queue with \arg bucket_count entries.
            Priorities will be computed by the <tt>PriorityFunctor</tt>
            given in \arg priority (i.e. <tt>priority(v)</tt> must result in an integer,
            where <tt>v</tt> is an instance of <tt>ValueType</tt>).
        */
    MappedBucketQueue(unsigned int bucket_count = 256,
                      PriorityFunctor const & priority = PriorityFunctor())
    : BaseType(bucket_count),
      get_priority_(priority)
    {}
    
        /** \brief Insert new element \arg v.
            Its priority is calculated by <tt>priority(v)</tt>,
            where <tt>priority</tt> is an instance of the 
            <tt>PriorityFunctor</tt> passed in the constructor.
            If the priority is outside the range <tt>[0, ..., bucket_count-1]</tt>,
            it is clamped to the range borders.
        */
    void push(value_type const & v)
    {
        priority_type index = get_priority_(v);
        
        // clamp index to the allowed range
        if(index > BaseType::maxIndex())
            index = BaseType::maxIndex();
        else if (index < 0)
            index = 0;
        
        BaseType::push(v, index);
    }
};

} // namespace vigra

#endif // VIGRA_BUCKET_QUEUE_HXX
