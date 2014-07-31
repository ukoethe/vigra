/************************************************************************/
/*                                                                      */
/*     Copyright 2014 by Thorsten Beier and Ullrich Koethe              */
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


//! @cond Doxygen_Suppress
#pragma once
#ifndef VIGRA_RANDOM_ACCESS_SET_HXX
#define VIGRA_RANDOM_ACCESS_SET_HXX

#include <vector>
#include <algorithm>
#include <utility>

namespace vigra {
   
/// set with O(n) insert and O(1) access
///
/// \tparam Key key and value type of the set
/// \tparam Alloc allocator of the set
/// RandomAccessSet has the same interface as std::set.
/// In addition, there is operator[].
/// \warning Values in set must not be changend through the mutable iterator
/// because doing so would potentially change the order of the values
/// \\ingroup datastructures
template<class Key,class Compare=std::less<Key>,class Alloc=std::allocator<Key> >
class RandomAccessSet {
private:
   /// type of the underlying vector
   typedef std::vector<Key,Alloc> VectorType;

public:
   // typedefs
   /// key type of the set
   typedef Key key_type;
   /// value type of the set
   typedef Key ValueType;
   /// value type of the set
   typedef Key value_type;
   /// comperator
   typedef Compare key_compare;
   /// value comperator
   typedef Compare value_compare;
   /// acclocator
   typedef Alloc  allocator_type;
   /// const reference type
   typedef typename Alloc::const_reference const_reference;
   /// iterator type
   typedef typename VectorType::iterator iterator;
   /// const iterator type
   typedef typename VectorType::const_iterator const_iterator;
   /// size type
   typedef typename VectorType::size_type size_type;
   /// difference type
   typedef typename VectorType::difference_type	difference_type;
   /// const pointer type
   typedef typename VectorType::const_pointer const_pointer;
   /// const reverse iterator
   typedef typename VectorType::const_reverse_iterator	const_reverse_iterator;

   // memeber functions:
   // constructor
   RandomAccessSet(const size_t, const Compare& compare=Compare(), const Alloc& alloc=Alloc());
   RandomAccessSet(const Compare& compare=Compare(), const Alloc& alloc=Alloc());
   template <class InputIterator>
      RandomAccessSet(InputIterator, InputIterator, const Compare& compare =Compare(), const Alloc & alloc=Alloc());
   RandomAccessSet(const RandomAccessSet&);

   // operator=
   RandomAccessSet& operator=(const RandomAccessSet &);
   // operator[]
   const value_type& operator[](const size_type) const;
   // iterators
   const_iterator begin() const;
   const_iterator end() const;
   const_iterator rbegin() const;
   const_iterator rend() const;

   iterator begin();
   iterator end();
   iterator rbegin();
   iterator rend();
   bool empty() const;
   size_type size() const;
   size_type max_size() const;
   std::pair< const_iterator,bool> insert(const value_type&);
   template <class InputIterator>
      void insert(InputIterator, InputIterator);
   const_iterator insert(const_iterator , const value_type&);
   void erase(iterator position);
   size_type erase(const key_type& );
   void erase( const_iterator, const_iterator);
   void swap(RandomAccessSet&);
   void clear();
   key_compare key_comp() const;
   value_compare value_comp() const;
   const_iterator find(const key_type&) const;
   iterator find(const key_type&);
   size_type count(const key_type&) const;
   const_iterator lower_bound(const key_type&) const;
   const_iterator upper_bound(const key_type&) const;
   std::pair<const_iterator,const_iterator> equal_range(const key_type&) const;
   iterator lower_bound(const key_type&) ;
   iterator upper_bound(const key_type&) ;
   std::pair<iterator,iterator> equal_range(const key_type&) ;
   allocator_type get_allocator() const;

   // std vector functions
   void reserve(const size_t size){
       vector_.reserve(size);
   }
   size_t capacity()const{
       return vector_.capacity();
   }
   
   template<class SET>
   void assignFromSet(const SET & set){
      vector_.assign(set.begin(),set.end());
   }

private:
   std::vector<Key> vector_;
   Compare compare_;
};

/// constructor
/// \param reserveSize reserve /allocate space
/// \param compare comperator
/// \param alloc allocator
template<class Key, class Compare, class Alloc>
inline
RandomAccessSet<Key,Compare,Alloc>::RandomAccessSet
(
   const size_t reserveSize,
   const Compare& compare,
   const Alloc& alloc
)
:  vector_(alloc),
   compare_(compare) {
   vector_.reserve(reserveSize);
}

/// const access values
/// \param index index of the value in the set
/// \return value / key at the position index
template<class Key, class Compare, class Alloc>
inline const typename RandomAccessSet<Key,Compare,Alloc>::value_type&
RandomAccessSet<Key,Compare,Alloc>::operator[]
(
   const typename RandomAccessSet<Key,Compare,Alloc>::size_type  index
) const
{
   return vector_[index];
}

/// constructor
/// \param compare comperator
/// \allloc allocator
template<class Key, class Compare, class Alloc>
inline
RandomAccessSet<Key,Compare,Alloc>::RandomAccessSet
(
   const Compare& compare,
   const Alloc& alloc
)
:  vector_(alloc),
   compare_(compare)
{}

/// constructor
/// \tparam InputIterator (key/value) input iterator
/// \param beginInput
/// \param endInput
template<class Key, class Compare, class Alloc>
template <class InputIterator>
inline
RandomAccessSet<Key,Compare,Alloc>::RandomAccessSet
(
   InputIterator beginInput,
   InputIterator endInput,
   const Compare& compare,
   const Alloc& alloc
)
:  vector_(alloc),
   compare_(compare)
{
   while(beginInput!=endInput) {
      this->insert(*beginInput);
      ++beginInput;
   }
}

/// copy constructor
/// \param src other random access set
template<class Key, class Compare, class Alloc>
inline
RandomAccessSet<Key,Compare,Alloc>::RandomAccessSet
(
   const RandomAccessSet<Key,Compare,Alloc>& src
)
:  vector_(src.vector_),
   compare_(src.compare_) {
}

/// assignment operator
/// \param src other random access set
template<class Key, class Compare, class Alloc>
inline RandomAccessSet<Key,Compare,Alloc>&
RandomAccessSet<Key,Compare,Alloc>::operator=
(
   const RandomAccessSet<Key,Compare,Alloc> & src
)
{
    if(this!=&src) {
      vector_=src.vector_;
      compare_=src.compare_;
   }
   return *this;
}

/// const begin iterator
/// \returns begin iterator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::const_iterator
RandomAccessSet<Key,Compare,Alloc>::begin() const
{
   return vector_.begin();
}

/// const end iterator
/// \returns end iterator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::const_iterator
RandomAccessSet<Key,Compare,Alloc>::end() const
{
    return vector_.end();
}
/// reverse const begin iterator
/// \returns reverse begin iterator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::const_iterator
RandomAccessSet<Key,Compare,Alloc>::rbegin() const
{
   return vector_.rbegin();
}

/// reverse const end iterator
/// \param reverse end iterator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::const_iterator
RandomAccessSet<Key,Compare,Alloc>::rend() const
{
    return vector_.rend();
}

/// begin iterator
/// \param begin iterator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::iterator
RandomAccessSet<Key,Compare,Alloc>::begin()
{
   return vector_.begin();
}

/// end iterator
/// \param end iterator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::iterator
RandomAccessSet<Key,Compare,Alloc>::end()
{
    return vector_.end();
}

/// reverse  begin iterator
/// \param reverse begin iterator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::iterator
RandomAccessSet<Key,Compare,Alloc>::rbegin()
{
   return vector_.rbegin();
}

/// reverse end iterator
/// \param reverse end iterator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::iterator
RandomAccessSet<Key,Compare,Alloc>::rend()
{
    return vector_.rend();
}

/// query if the set is empty
/// \return true if empty
template<class Key, class Compare, class Alloc>
inline bool
RandomAccessSet<Key,Compare,Alloc>::empty() const
{
   return vector_.empty();
}

/// number of elements of the set
/// \returns number of elements in the set
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::size_type
RandomAccessSet<Key,Compare,Alloc>::size() const
{
   return vector_.size();
}

/// maximum size of the underlying container
/// \return the maximum size
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::size_type
RandomAccessSet<Key,Compare,Alloc>::max_size() const
{
   return vector_.max_size();
}

// modifiers
/// insert an element into the set
///
/// \param value element to insert
/// \return pair in which the first entry is an iterator pointing to inserted
/// value and the second entry is true iff the value had not already been in the
/// set
template<class Key, class Compare, class Alloc>
inline std::pair<typename RandomAccessSet<Key,Compare,Alloc>::const_iterator,bool>
RandomAccessSet<Key,Compare,Alloc>::insert
(
   const typename RandomAccessSet<Key,Compare,Alloc>::value_type& value
) {
   bool found(true);
   iterator i(lower_bound(static_cast<Key>(value)));
   if(i == end() || compare_(static_cast<Key>(value), *i)) {
      i = vector_.insert(i, static_cast<Key>(value));
      found = false;
   }
   return std::make_pair(i, !found);
}

/// insert a sequence of elements
///
/// \param first iterator to the first element
/// \param last iterator to the last element
template<class Key, class Compare, class Alloc>
template <class InputIterator>
inline void
RandomAccessSet<Key,Compare,Alloc>::insert
(
   InputIterator first,
   InputIterator last
)
{
   while(first!=last) {
      this->insert(*first);
      ++first;
   }
}

/// insert a sequence of elements with a hint for the position
///
/// \param position iterator to the position
/// \param value element to insert
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::const_iterator
RandomAccessSet<Key,Compare,Alloc>::insert
(
   typename RandomAccessSet<Key,Compare,Alloc>::const_iterator position,
   const typename RandomAccessSet<Key,Compare,Alloc>::value_type& value
)
{
   if((position == begin() || this->operator()(*(position-1),value))
   && (position == end() || this->operator()(value, *position))) {
       return vector_.insert(position, value);
   }
   return insert(value).first;
}

/// erase an element
/// \param position iterator to the position
template<class Key, class Compare, class Alloc>
inline void
RandomAccessSet<Key,Compare,Alloc>::erase
(
   typename RandomAccessSet<Key,Compare,Alloc>::iterator position
)
{
   vector_.erase(position);
}

/// erease and element
/// \param x element
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::size_type
RandomAccessSet<Key,Compare,Alloc>::erase
(
   const typename RandomAccessSet<Key,Compare,Alloc>::key_type& x
)
{
   iterator i =find(x);
   if(i!=vector_.end())
   {
      erase(i);
      return 1;
   }
   return 0;
}

/// erase a sequence of elements
/// \param first iterator to the beginning of the sequence to erase
/// \param last iterator to the end of the sequence to erase
template<class Key, class Compare, class Alloc>
inline void
RandomAccessSet<Key,Compare,Alloc>::erase
(
   const typename RandomAccessSet<Key,Compare,Alloc>::const_iterator first,
   const typename RandomAccessSet<Key,Compare,Alloc>::const_iterator last
)
{
   vector_.erase(first,last);
}

/// swap random access sets
/// \param rhs set to swap with
template<class Key, class Compare, class Alloc>
inline void
RandomAccessSet<Key,Compare,Alloc>::swap
(
   RandomAccessSet<Key,Compare,Alloc>& rhs
)
{
   vector_.swap(rhs.vector_);
   std::swap(compare_, rhs.compare_);
}

/// clear the set
///
/// erases all elements
template<class Key, class Compare, class Alloc>
inline void
RandomAccessSet<Key,Compare,Alloc>::clear()
{
   vector_.clear();
}

/// key comparator
/// \return key comparator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::key_compare
RandomAccessSet<Key,Compare,Alloc>::key_comp() const
{
   return compare_;
}

/// value comparator
/// \return value comparator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::value_compare
RandomAccessSet<Key,Compare,Alloc>::value_comp() const
{
   return compare_;
}

/// find an element
/// \param value element
/// \return const_iterator to the position where element was found or end
/// iterator if the element was not found
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::const_iterator
RandomAccessSet<Key,Compare,Alloc>::find
(
   const typename RandomAccessSet<Key,Compare,Alloc>::key_type& value
) const
{
   const_iterator i(lower_bound(value));
   if (i != end() && compare_(value, *i))
   {
       i = end();
   }
   return i;
}

/// find an element
/// \param value element
/// \return iterator to the position where the element was found or end
/// iterator if the element was not found
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::iterator
RandomAccessSet<Key,Compare,Alloc>::find
(
   const typename RandomAccessSet<Key,Compare,Alloc>::key_type& value
)
{
   iterator i(lower_bound(value));
   if (i != end() && compare_(value, *i))
   {
       i = end();
   }
   return i;
}

/// count elements
/// \param element
/// \return zero or one
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::size_type
RandomAccessSet<Key,Compare,Alloc>::count
(
   const typename RandomAccessSet<Key,Compare,Alloc>::key_type&  value
) const
{
   return find(value) != end() ? 1 : 0;
}

/// lower bound
/// \param value
/// \return iterator to lower bound
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::const_iterator
RandomAccessSet<Key,Compare,Alloc>::lower_bound
(
   const typename RandomAccessSet<Key,Compare,Alloc>::key_type& value
) const
{
   return std::lower_bound(vector_.begin(), vector_.end(), value, compare_);
}

/// lower bound
/// \param value
/// \return iterator to lower bound
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::iterator
RandomAccessSet<Key,Compare,Alloc>::lower_bound
(
   const typename RandomAccessSet<Key,Compare,Alloc>::key_type& value
)
{
   return std::lower_bound(vector_.begin(), vector_.end(), value, compare_);
}

/// upper bound
/// \param value
/// \return iterator to upper bound
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::const_iterator
RandomAccessSet<Key,Compare,Alloc>::upper_bound
(
   const typename RandomAccessSet<Key,Compare,Alloc>::key_type& value
) const
{
   return std::upper_bound(vector_.begin(), vector_.end(), value, compare_);
}

/// upper bound
/// \param value
/// \return iterator to upper bound
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::iterator
RandomAccessSet<Key,Compare,Alloc>::upper_bound
(
   const typename RandomAccessSet<Key,Compare,Alloc>::key_type& value
)
{
   return std::upper_bound(vector_.begin(), vector_.end(), value, compare_);
}

/// equal range
/// \param value
/// \return iterator pair to lower equal range
template<class Key, class Compare, class Alloc>
inline std::pair<typename RandomAccessSet<Key,Compare,Alloc>::const_iterator,typename RandomAccessSet<Key,Compare,Alloc>::const_iterator>
RandomAccessSet<Key,Compare,Alloc>::equal_range
(
   const typename RandomAccessSet<Key,Compare,Alloc>::key_type& value
) const
{
   return std::equal_range(vector_.begin(), vector_.end(), value, compare_);
}

/// equal range
/// \param value
/// \return iterator pair to lower equal range
template<class Key, class Compare, class Alloc>
inline std::pair<typename RandomAccessSet<Key,Compare,Alloc>::iterator,typename RandomAccessSet<Key,Compare,Alloc>::iterator>
RandomAccessSet<Key,Compare,Alloc>::equal_range
(
   const typename RandomAccessSet<Key,Compare,Alloc>::key_type& value
)
{
   return std::equal_range(vector_.begin(), vector_.end(), value, compare_);
}
/// allocators
/// \return allocator
template<class Key, class Compare, class Alloc>
inline typename RandomAccessSet<Key,Compare,Alloc>::allocator_type
RandomAccessSet<Key,Compare,Alloc>::get_allocator() const
{
   return vector_.get_allocator();
}

} // namespace vigra

#endif // VIGRA_RANDOM_ACCESS_SET_HXX
//! @endcond