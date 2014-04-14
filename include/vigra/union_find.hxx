/************************************************************************/
/*                                                                      */
/*               Copyright 2008-2009 by Ullrich Koethe                  */
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


#ifndef VIGRA_UNION_FIND_HXX
#define VIGRA_UNION_FIND_HXX

/*std*/
#include <map>

/*vigra*/
#include "config.hxx"
#include "error.hxx"
#include "array_vector.hxx"
#include "iteratoradapter.hxx"

namespace vigra {

namespace detail {

template <class T, class IsSigned = VigraFalseType>
struct UnionFindAccessorImpl
{
    static const T max_label = NumericTraits<T>::maxConst >> 1;
    static const T anchor_bit = ~max_label;
    
    static T max()
    {
        return max_label;
    }
    
    static T deletedAnchor()
    {
        return NumericTraits<T>::maxConst;
    }
    
    static bool isAnchor(T const & t)
    {
        return (t & anchor_bit) != 0;
    }
    
    static bool isValidAnchor(T const & t)
    {
        return isAnchor(t) && t != deletedAnchor();
    }
    
    static bool notAnchor(T const & t)
    {
        return (t & anchor_bit) == 0;
    }
    
    static T toAnchor(T const & t)
    {
        return t | anchor_bit;
    }
    
    static T fromAnchor(T const & t)
    {
        return t & max_label;
    }
};

template <class T>
struct UnionFindAccessorImpl<T, VigraTrueType>
{
    static T max()
    {
        return NumericTraits<T>::max();
    }
    
    static T deletedAnchor()
    {
        return NumericTraits<T>::min();
    }
    
    static bool isAnchor(T const & t)
    {
        return t <= 0;
    }
    
    static bool isValidAnchor(T const & t)
    {
        return isAnchor(t) && t != deletedAnchor();
    }
    
    static bool notAnchor(T const & t)
    {
        return t > 0;
    }
    
    static T toAnchor(T const & t)
    {
        return -t;
    }
    
    static T fromAnchor(T const & t)
    {
        return -t;
    }
};

template <class Array, class LabelAccessor>
class UnionFindIteratorPolicy
{
  public:
    typedef UnionFindIteratorPolicy                BaseType;
    typedef typename Array::difference_type        value_type;
    typedef typename Array::difference_type        difference_type;
    typedef value_type const &                     reference;
    typedef value_type const &                     index_reference;
    typedef value_type const *                     pointer;
    typedef typename std::forward_iterator_tag     iterator_category;

    Array const & array_;
    value_type index_;
    
    UnionFindIteratorPolicy(Array const & array, value_type index=0)
    : array_(array)
    , index_(index)
    {}
    
    static void initialize(BaseType & d) 
    {
        advanceToAnchor(d);
    }

    static reference dereference(BaseType const & d)
    { 
        return d.index_;
    }

    static bool equal(BaseType const & d1, BaseType const & d2)
    {
        return d1.index_ == d2.index_;
    }

    static bool less(BaseType const & d1, BaseType const & d2)
    {
        return d1.index_ < d2.index_;
    }

    static void increment(BaseType & d)
    {
        ++d.index_; 
        advanceToAnchor(d);
    }

    static void advanceToAnchor(BaseType & d)
    {
        while(d.index_ < (value_type)d.array_.size()-1 && 
              !LabelAccessor::isValidAnchor(d.array_[d.index_]))
        {
            ++d.index_;
        }
    }
};

} // namespace detail

template <class T>
class UnionFindArray
{
    typedef ArrayVector<T>                                             LabelArray;
    typedef typename LabelArray::difference_type                       IndexType;
    typedef detail::UnionFindAccessorImpl<T, 
                          typename NumericTraits<T>::isSigned>         LabelAccessor;
    typedef detail::UnionFindIteratorPolicy<LabelArray, LabelAccessor> IteratorPolicy;
    typedef IteratorAdaptor<IteratorPolicy>                            iterator;
    typedef iterator                                                   const_iterator;

    mutable ArrayVector<T> labels_;
    
  public:
    UnionFindArray(T next_free_index = 1)
    {
        for(T k=0; k <= next_free_index; ++k)
            labels_.push_back(LabelAccessor::toAnchor(k));
    }
    
    const_iterator begin(unsigned int start_at=0) const
    {
        return const_iterator(IteratorPolicy(labels_, start_at));
    }
    
    const_iterator end() const
    {
        return const_iterator(IteratorPolicy(labels_, labels_.size()-1));
    }
    
    T nextFreeIndex() const
    {
        return T(labels_.size() - 1);
    }
    
    T findIndex(T index) const
    {
        IndexType root = index;
        while(LabelAccessor::notAnchor(labels_[root]))
            root = (IndexType)labels_[root];
        // path compression
        while((IndexType)index != root)
        {
            T next = labels_[(IndexType)index];
            labels_[(IndexType)index] = root;
            index = next;
        }
        return (T)root;
    } 
    
    T findLabel(T index) const
    {
        return LabelAccessor::fromAnchor(labels_[findIndex(index)]);
    }
    
    void deleteIndex(T index)
    {
        labels_[findIndex(index)] = LabelAccessor::deletedAnchor();
    }
    
        // this function does not yet check for deletedIndex()
    T makeUnion(T l1, T l2)
    {
        IndexType i1 = findIndex(l1);
        IndexType i2 = findIndex(l2);
        if(i1 == i2)
        {
            return i1;
        }
        else if(i1 < i2)
        {
            labels_[i2] = i1;
            return (T)i1;
        }
        else
        {
            labels_[i1] = i2;
            return (T)i2;
        }
    }
    
    T finalizeIndex(T index)
    {
        if(index == (T)labels_.size()-1)
        {
            // indeed a new region
            vigra_invariant(index < LabelAccessor::max(),
                    "connected components: Need more labels than can be represented in the destination type.");
            // create new back entry
            labels_.push_back(LabelAccessor::toAnchor((T)labels_.size()));
        }
        else
        {
            // no new index => reset the back entry of the index array
            labels_.back() = LabelAccessor::toAnchor((T)labels_.size()-1);
        }
        return index;
    }
    
    T makeNewIndex() 
    {
        T index = LabelAccessor::fromAnchor(labels_.back());
        vigra_invariant(index < LabelAccessor::max(),
          "connected components: Need more labels than can be represented in the destination type.");
        labels_.push_back(LabelAccessor::toAnchor((T)labels_.size()));
        return index;
    }
    
    unsigned int makeContiguous()
    {
        // compress trees
        unsigned int count = 0; 
        for(IndexType i=0; i<(IndexType)(labels_.size()-1); ++i)
        {
            if(LabelAccessor::isValidAnchor(labels_[i]))
            {
                    labels_[i] = LabelAccessor::toAnchor((T)count++);
            }
            else
            {
                    labels_[i] = findIndex(i);  // path compression
            }
        }
        return count-1;   
    }
};

} // namespace vigra

#endif // VIGRA_UNION_FIND_HXX
