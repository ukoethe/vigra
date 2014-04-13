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
    
    static bool isAnchor(T const & t)
    {
        return (t & anchor_bit) != 0;
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
    
    static bool isAnchor(T const & t)
    {
        return t <= 0;
    }
    
    static bool notAnchor(T const & t)
    {
        return t > 0;
    }
    
    static T toAnchor(T const & t)
    {
        return t < 0 ? t : -t;
    }
    
    static T fromAnchor(T const & t)
    {
        return t < 0 ? -t : t;
    }
};

template <class T>
class UnionFindArray
{
    typedef typename ArrayVector<T>::difference_type IndexType;
    typedef UnionFindAccessorImpl<T, typename NumericTraits<T>::isSigned> LabelAccessor;

    mutable ArrayVector<T> labels_;
    
  public:
    UnionFindArray(T next_free_label = 1)
    {
        for(T k=0; k <= next_free_label; ++k)
            labels_.push_back(LabelAccessor::toAnchor(k));
    }
    
    T nextFreeIndex() const
    {
        return T(labels_.size() - 1);
    }
    
    T findIndex(T label) const
    {
        IndexType root = label;
        while(LabelAccessor::notAnchor(labels_[root]))
            root = (IndexType)labels_[root];
        // path compression
        while((IndexType)label != root)
        {
            T next = labels_[(IndexType)label];
            labels_[(IndexType)label] = root;
            label = next;
        }
        return (T)root;
    } 
    
    T findLabel(T label) const
    {
        return LabelAccessor::fromAnchor(labels_[findIndex(label)]);
    } 
    
    T makeUnion(T l1, T l2)
    {
        IndexType i1 = findIndex(l1);
        IndexType i2 = findIndex(l2);
        if(i1 == i2)
        {
            return i1;
        }
        if(i1 < i2)
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
            if(LabelAccessor::isAnchor(labels_[i]))
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

} // namespace detail
} // namespace vigra

#endif // VIGRA_UNION_FIND_HXX
