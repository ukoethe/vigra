#pragma once
#ifndef PARTITON_PARTITION_HXX
#define PARTITON_PARTITION_HXX

#include <vector>
#include <map>
#include "macros.hxx"

#include <boost/iterator/iterator_facade.hpp>

namespace partition {

template<class T>
class Partition;



template<class T>
struct  ConstRepIter
:  public boost::iterator_facade<
      ConstRepIter<T>,
      T const,
      boost::forward_traversal_tag
   >

{
   typedef Partition<T> PartitionType;
   ConstRepIter(const PartitionType & p,const T cr)
   :  partition_(&p),
      currentRep_(cr){
   }


   ConstRepIter(){

   }
   //explicit ConstRepIter(node_base* p)


   bool isBegin()const{
      return currentRep_==partition_->firstRep();
   }
   bool isEnd()const{
      return currentRep_>partition_->lastRep();
   }

   bool equal(const ConstRepIter & other)const{
      return (this->isEnd() && other.isEnd() ) || this->currentRep_==other.currentRep_;
   }

   void increment(){
      if(partition_->jumpVec_[currentRep_].second==0){
         VIGRA_ASSERT_OP(currentRep_,==,partition_->lastRep());
         currentRep_+=1;
      }
      else{
         currentRep_+=partition_->jumpVec_[currentRep_].second;
      }
   }

   const T & dereference()const{
      return currentRep_;
   }



   const PartitionType * partition_;
   T currentRep_;
   
};

/// Disjoint set data structure with path compression.
/// \ingroup datastructures
template<class T>
class Partition {
public:
   friend class ConstRepIter<T>;
   typedef T value_type;
   typedef std::size_t SizeTType;
   Partition();
   Partition(const value_type&);

   // query
   value_type find(const value_type&) const; // without path compression
   value_type find(value_type); // with path compression
   value_type numberOfElements() const;
   value_type numberOfSets() const;
   template<class Iterator> void elementLabeling(Iterator) const;
   template<class Iterator> void representatives(Iterator) const;
   void representativeLabeling(std::map<value_type, value_type>&) const;

   // manipulation
   void reset(const value_type&);
   void merge(value_type, value_type);
   void insert(const value_type&);

   template<class ITER>
   value_type multiMerge(const value_type a,ITER begin,ITER end){
      while(begin!=end){
         this->merge(a,*begin);
         ++begin;
      }
      return this->find(a);
   }

   value_type firstRep()const{
      return firstRep_;
   }
   value_type lastRep()const{
      return lastRep_;
   }
   typedef ConstRepIter<T> const_iterator;
   const_iterator begin()const{
      return ConstRepIter<T>(*this,firstRep_);
   }
   const_iterator end()const{
      return ConstRepIter<T>(*this,lastRep_+1);
   }

   void eraseElement(const value_type & value){
      const size_t notRep=value;
      const size_t jumpMinus = jumpVec_[notRep].first;
      const size_t jumpPlus  = jumpVec_[notRep].second;
      VIGRA_ASSERT_OP(jumpMinus+jumpPlus,>=,0);

      if(jumpMinus==0){
         VIGRA_ASSERT_OP(firstRep_,==,notRep);
         const size_t nextRep = notRep+jumpPlus;
         firstRep_=nextRep;
         jumpVec_[nextRep].first=0;
      }
      else if(jumpPlus==0){
         VIGRA_ASSERT_OP(lastRep_,==,notRep);
         const size_t prevRep = notRep-jumpMinus;
         lastRep_=prevRep;
         jumpVec_[prevRep].second=0;

      }
      else{
         const size_t nextRep = notRep+jumpPlus;
         const size_t prevRep = notRep-jumpMinus;
         jumpVec_[nextRep].first+=jumpVec_[notRep].first;
         jumpVec_[prevRep].second+=jumpVec_[notRep].second;
      }   
   }

private:
   std::vector<value_type> parents_;
   std::vector<value_type> ranks_;
   std::vector< std::pair< size_t, size_t> > jumpVec_;
   value_type firstRep_;
   value_type lastRep_;
   value_type numberOfElements_;
   value_type numberOfSets_;
};

/// Construct a partition.
template<class T>
Partition<T>::Partition()
: parents_(),
  ranks_(),
  jumpVec_(),
  firstRep_(0),
  lastRep_(0),
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
: parents_(static_cast<SizeTType>(size)),
  ranks_(static_cast<SizeTType>(size)),
  jumpVec_(static_cast<SizeTType>(size)),
  firstRep_(0),
  lastRep_(static_cast<SizeTType>(size)-1),
  numberOfElements_(size),
  numberOfSets_(size)
{
   for(T j=0; j<size; ++j) {
      parents_[static_cast<SizeTType>(j)] = j;
   }

   jumpVec_.front().first=0;
   jumpVec_.front().second=1;
   for(T j=1; j<size-1;++j){
      jumpVec_[j].first =1;
      jumpVec_[j].second=1;
   }
   jumpVec_.back().first=1;
   jumpVec_.back().second=0;
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
   ranks_.resize(static_cast<SizeTType>(size));
   parents_.resize(static_cast<SizeTType>(size));
   jumpVec_.resize(static_cast<SizeTType>(size));
   firstRep_=0;
   lastRep_=static_cast<SizeTType>(size)-1;
   for(T j=0; j<size; ++j) {
      ranks_[static_cast<SizeTType>(j)] = 0;
      parents_[static_cast<SizeTType>(j)] = j;
   }

   jumpVec_.front().first=0;
   jumpVec_.front().second=1;
   for(T j=1; j<size-1;++j){
      jumpVec_[j].first =1;
      jumpVec_[j].second=1;
   }
   jumpVec_.back().first=1;
   jumpVec_.back().second=0;
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
   while(parents_[static_cast<SizeTType>(root)] != root) {
      root = parents_[static_cast<SizeTType>(root)];
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
   while(parents_[static_cast<SizeTType>(root)] != root) {
      root = parents_[static_cast<SizeTType>(root)];
   }
   // path compression
   while(element != root) {
      value_type tmp = parents_[static_cast<SizeTType>(element)];
      parents_[static_cast<SizeTType>(element)] = root;
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
   T rep,notRep;
   if(ranks_[static_cast<SizeTType>(element1)] < ranks_[static_cast<SizeTType>(element2)]) {
      parents_[static_cast<SizeTType>(element1)] = element2;
      --numberOfSets_;
      //rep=element2;
      notRep=element1;
   }
   else if(ranks_[static_cast<SizeTType>(element1)] > ranks_[static_cast<SizeTType>(element2)]) {
      parents_[static_cast<SizeTType>(element2)] = element1;
      --numberOfSets_;
      //rep=element1;
      notRep=element2;
   }
   else if(element1 != element2) {
      parents_[static_cast<SizeTType>(element2)] = element1;
      ++ranks_[static_cast<SizeTType>(element1)];
      --numberOfSets_;
      //rep=element1;
      notRep=element2;
   }

   this->eraseElement(notRep);
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
   ranks_.insert(ranks_.end(), static_cast<SizeTType>(number), T(0));
   parents_.insert(parents_.end(), static_cast<SizeTType>(number), T(0));
   for(value_type j=numberOfElements_; j<numberOfElements_+number; ++j) {
      parents_[static_cast<SizeTType>(j)] = j;
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
      if(parents_[static_cast<SizeTType>(j)] == j) {
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
inline void
Partition<T>::representativeLabeling
(
   std::map<value_type, value_type>& out
) const
{
   out.clear();
   std::vector<value_type> r(static_cast<SizeTType>(numberOfSets()));
   representatives(r.begin());
   for(T j=0; j<numberOfSets(); ++j) {
      out[ r[static_cast<SizeTType>(j)] ] = j;
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

} // namespace partition

#endif // #ifndef PARTITON_PARTITION_HXX
