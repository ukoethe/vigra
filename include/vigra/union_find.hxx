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

template <class T>
class UnionFindArray
{
    typedef typename ArrayVector<T>::difference_type IndexType;
    mutable ArrayVector<T> labels_;
    
  public:
    UnionFindArray(T next_free_label = 1)
    {
        for(T k=0; k <= next_free_label; ++k)
            labels_.push_back(k);
    }
    
    T nextFreeLabel() const
    {
        return labels_.back();
    }
    
    T find(T label) const
    {
        T root = label;
        while(root != labels_[(IndexType)root])
            root = labels_[(IndexType)root];
        // path compression
        while(label != root)
        {
            T next = labels_[(IndexType)label];
            labels_[(IndexType)label] = root;
            label = next;
        }
        return root;
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
            labels_.push_back((T)labels_.size());
        }
        else
        {
            // no new label => reset the back entry of the label array
            labels_.back() = (T)labels_.size()-1;
        }
        return label;
    }
    
    T makeNewLabel() 
    {
        T label = labels_.back();
        vigra_invariant(label < NumericTraits<T>::max(),
          "connected components: Need more labels than can be represented in the destination type.");
        labels_.push_back((T)labels_.size());
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
                    labels_[i] = (T)count++;
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



/// Disjoint set data structure with path compression.
/// \ingroup datastructures
template<class T = size_t>
class Partition {
public:
   typedef T value_type;

   Partition();
   Partition(const value_type&);

   // query
   value_type find(const value_type&) const; // without path compression
   value_type find(value_type); // with path compression
   value_type numberOfElements() const;
   value_type numberOfSets() const;
   template<class Iterator> void elementLabeling(Iterator) const;
   template<class Iterator> void representatives(Iterator) const;
   
   template<class MAP> //map from value_type to value_type
   void representativeLabeling(MAP&) const;

   // manipulation
   void reset(const value_type&);
   void merge(value_type, value_type);
   void insert(const value_type&);

private:
   std::vector<value_type> parents_;
   std::vector<value_type> ranks_;
   value_type numberOfElements_;
   value_type numberOfSets_;
};

/// Construct a partition.
template<class T>
Partition<T>::Partition()
: parents_(),
  ranks_(),
  numberOfElements_(0),
  numberOfSets_(0)
{}

/// Construct a partition.
///
/// \param size Number of distinct sets.
///
template<class T>
inline
Partition<T>::Partition
(
   const value_type& size
)
: parents_(static_cast<size_t>(size)),
  ranks_(static_cast<size_t>(size)),
  numberOfElements_(size),
  numberOfSets_(size)
{
   for(T j=0; j<size; ++j) {
      parents_[static_cast<size_t>(j)] = j;
   }
}

/// Reset a partition such that each set contains precisely one element
///
/// \param size Number of distinct sets.
///
template<class T>
inline void
Partition<T>::reset
(
   const value_type& size
)
{
   numberOfElements_ = size;
   numberOfSets_ = size;
   ranks_.resize(static_cast<size_t>(size));
   parents_.resize(static_cast<size_t>(size));
   for(T j=0; j<size; ++j) {
      ranks_[static_cast<size_t>(j)] = 0;
      parents_[static_cast<size_t>(j)] = j;
   }
}

/// Find the representative element of the set that contains the given element.
///
/// This constant function does not compress the search path.
///
/// \param element Element.
///
template<class T>
inline typename Partition<T>::value_type
Partition<T>::find
(
   const value_type& element
) const
{
   // find the root
   value_type root = element;
   while(parents_[static_cast<size_t>(root)] != root) {
      root = parents_[static_cast<size_t>(root)];
   }
   return root;
}

/// Find the representative element of the set that contains the given element.
///
/// This mutable function compresses the search path.
///
/// \param element Element.
///
template<class T>
inline typename Partition<T>::value_type
Partition<T>::find
(
   value_type element // copy to work with
)
{
   // find the root
   value_type root = element;
   while(parents_[static_cast<size_t>(root)] != root) {
      root = parents_[static_cast<size_t>(root)];
   }
   // path compression
   while(element != root) {
      value_type tmp = parents_[static_cast<size_t>(element)];
      parents_[static_cast<size_t>(element)] = root;
      element = tmp;
   }
   return root;
}

/// Merge two sets.
///
/// \param element1 Element in the first set.
/// \param element2 Element in the second set.
///
template<class T>
inline void
Partition<T>::merge
(
   value_type element1,
   value_type element2
)
{
   // merge by rank
   element1 = find(element1);
   element2 = find(element2);
   if(ranks_[static_cast<size_t>(element1)] < ranks_[static_cast<size_t>(element2)]) {
      parents_[static_cast<size_t>(element1)] = element2;
      --numberOfSets_;
   }
   else if(ranks_[static_cast<size_t>(element1)] > ranks_[static_cast<size_t>(element2)]) {
      parents_[static_cast<size_t>(element2)] = element1;
      --numberOfSets_;
   }
   else if(element1 != element2) {
      parents_[static_cast<size_t>(element2)] = element1;
      ++ranks_[static_cast<size_t>(element1)];
      --numberOfSets_;
   }
}

/// Insert new sets.
///
/// \param number Number of sets to insert.
///
template<class T>
inline void
Partition<T>::insert
(
   const value_type& number
)
{
   ranks_.insert(ranks_.end(), static_cast<size_t>(number), T(0));
   parents_.insert(parents_.end(), static_cast<size_t>(number), T(0));
   for(value_type j=numberOfElements_; j<numberOfElements_+number; ++j) {
      parents_[static_cast<size_t>(j)] = j;
   }
   numberOfElements_ += number;
   numberOfSets_ += number;
}

/// Output all elements which are set representatives.
///
/// \param it (Output) Iterator into a container.
///
template<class T>
template<class Iterator>
inline void
Partition<T>::representatives
(
   Iterator it
) const
{
   for(value_type j=0; j<numberOfElements(); ++j) {
      if(parents_[static_cast<size_t>(j)] == j) {
         *it = j;
         ++it;
      }
   }
}

/// Output a continuous labeling of the representative elements.
///
/// \param out (Output) A map that assigns each representative element to its label.
///
template<class T>
template<class MAP>
inline void
Partition<T>::representativeLabeling
(
   MAP& out
) const
{
   out.clear();
   std::vector<value_type> r(static_cast<size_t>(numberOfSets()));
   representatives(r.begin());
   for(T j=0; j<numberOfSets(); ++j) {
      out[ r[static_cast<size_t>(j)] ] = j;
   }
}

/// Output a continuous labeling of all elements.
///
/// \param out (Output) Iterator into a container in which the j-th entry is the label of the element j.
///
template<class T>
template<class Iterator>
inline void
Partition<T>::elementLabeling
(
   Iterator out
) const
{
   std::map<value_type, value_type> rl;
   representativeLabeling(rl);
   for(value_type j=0; j<numberOfElements(); ++j) {
      *out = rl[find(j)];
      ++out;
   }
}

template<class T>
inline typename Partition<T>::value_type
Partition<T>::numberOfElements() const
{
   return numberOfElements_;
}

template<class T>
inline typename Partition<T>::value_type
Partition<T>::numberOfSets() const
{
   return numberOfSets_;
}

} // namespace detail
} // namespace vigra

#endif // VIGRA_UNION_FIND_HXX
