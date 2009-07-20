/************************************************************************/
/*                                                                      */
/*               Copyright 2008-2009 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#include "config.hxx"
#include "error.hxx"
#include "array_vector.hxx"

namespace vigra {

namespace detail {

template <class T>
class UnionFindArray
{
    typedef typename ArrayVector<T>::difference_type IndexType;
    ArrayVector<T> labels_;
    
  public:
    UnionFindArray()
    {
        labels_.push_back(0);
        labels_.push_back(1);
    }
    
    T nextFreeLabel() const
    {
        return labels_.back();
    }
    
    T find(T label) const
    {
        while(label != labels_[(IndexType)label])
            label = labels_[(IndexType)label];
        return label;
    } 
    
    T makeUnion(T l1, T l2)
    {
        l1 = find(l1);
        l2 = find(l2);
        if(l1 <= l2)
        {
            labels_[(IndexType)l2] = l1;
            return l1;
        }
        else
        {
            labels_[(IndexType)l1] = l2;
            return l2;
        }
    }
    
    T finalizeLabel(T label)
    {
        if(label == (T)labels_.size()-1)
        {
            // indeed a new region
            vigra_invariant(label < NumericTraits<T>::max(),
                    "connected components: Need more labels than can be represented in the destination type.");
            // create new back entry
            labels_.push_back(labels_.size());
        }
        else
        {
            // no new label => reset the back entry of the label array
            labels_.back() = labels_.size()-1;
        }
        return label;
    }
    
    T makeNewLabel() 
    {
        T label = labels_.back();
        vigra_invariant(label < NumericTraits<T>::max(),
          "connected components: Need more labels than can be represented in the destination type.");
        labels_.push_back(labels_.size());
        return label;
    }
    
    unsigned int makeContiguous()
    {
        // compress trees
        unsigned int count = 0; 
        for(IndexType i=0; i<(IndexType)(labels_.size()-1); ++i)
        {
            if(labels_[i] == i)
            {
                    labels_[i] = count++;
            }
            else
            {
                    labels_[i] = labels_[(IndexType)labels_[i]]; 
            }
        }
        return count-1;   
    }
    
    T operator[](T label) const
    {
        return labels_[(IndexType)label];
    }
};

} // namespace detail

} // namespace vigra

#endif // VIGRA_UNION_FIND_HXX
