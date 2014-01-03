#pragma once
#ifndef CGP2D_ACCESSOR_ITERATOR
#define CGP2D_ACCESSOR_ITERATOR

#include <iterator>



namespace vigra {

/// \cond HIDDEN_SYMBOLS

/// Interface that simplifies the implementation of STL-compliant
/// random access iterators.
///
/// For an object a of class A (the template parameter) with with
/// member functions
///
///    A()
///    size_t size() const
///    A::reference operator[](const size_t)
///    const A::value_type& operator[](const size_t) const
///    bool operator==(const A&) const
///
/// AccessorIterator<A, true/false> is a constant/mutable
/// STL-compliant random access iterator to the sequence
/// a[0], ..., a[a.size()-1].
///
/// We call the class A an accessor. Here's a complete example,
/// an accessor for the container std::vector
///
///    template<class T, bool isConst>
///    class VectorAccessor {
///    public:
///       typedef T value_type;
///       typedef typename meta::if_<isConst, const value_type&, value_type&>::type reference;
///       typedef typename meta::if_<isConst, const value_type*, value_type*>::type pointer;
///       typedef typename meta::if_<isConst, const std::vector<value_type>&, std::vector<value_type>&>::type vector_reference;
///       typedef typename meta::if_<isConst, const std::vector<value_type>*, std::vector<value_type>*>::type vector_pointer;
///
///       VectorAccessor(vector_pointer v = 0)
///          : vector_(v) {}
///       VectorAccessor(vector_reference v)
///          : vector_(&v) {}
///       size_t size() const
///          { return vector_ == 0 ? 0 : vector_->size(); }
///       reference operator[](const size_t j)
///          { return (*vector_)[j]; }
///       const value_type& operator[](const size_t j) const
///          { return (*vector_)[j]; }
///       template<bool isConstLocal>
///          bool operator==(const VectorAccessor<T, isConstLocal>& other) const
///             { return vector_ == other.vector_; }
///
///    private:
///       vector_pointer vector_;
///    };
///
template<class A, bool isConst = true>
class AccessorIterator {
public:
   typedef A Accessor;
   typedef typename Accessor::value_type value_type;
   typedef const value_type*  pointer;
   typedef const value_type & reference;
   typedef size_t difference_type;
   typedef std::random_access_iterator_tag iterator_category;

   AccessorIterator(const Accessor& = Accessor(), const size_t = 0);

   template<bool isConstLocal>
      bool operator==(const AccessorIterator<A, isConstLocal>&) const;
   template<bool isConstLocal>
      bool operator!=(const AccessorIterator<A, isConstLocal>&) const;
   template<bool isConstLocal>
      bool operator<=(const AccessorIterator<A, isConstLocal>&) const;
   template<bool isConstLocal>
      bool operator>=(const AccessorIterator<A, isConstLocal>&) const;
   template<bool isConstLocal>
      bool operator<(const AccessorIterator<A, isConstLocal>&) const;
   template<bool isConstLocal>
      bool operator>(const AccessorIterator<A, isConstLocal>&) const;

   pointer operator->();
   const value_type* operator->() const;
   reference operator*();
   const value_type& operator*() const;
   reference operator[](const size_t);
   const value_type& operator[](const size_t) const;

   AccessorIterator<A, isConst>& operator++();
   AccessorIterator<A, isConst> operator++(int);
   AccessorIterator<A, isConst>& operator--();
   AccessorIterator<A, isConst> operator--(int);
   AccessorIterator<A, isConst> operator+=(const size_t);
   AccessorIterator<A, isConst> operator-=(const size_t);

   AccessorIterator<A, isConst> operator+(const size_t) const;
   AccessorIterator<A, isConst> operator-(const size_t) const;
   template<bool isConstLocal>
      size_t operator-(const AccessorIterator<A, isConstLocal>&) const;

private:
   void testInvariant() const;

   Accessor accessor_;
   size_t index_;
};

// binary arithmetic operators
template<class A, bool isConst>
   AccessorIterator<A, isConst> operator+(const size_t, const AccessorIterator<A, isConst>&);

// implementation of AccessorIterator

template<class A, bool isConst>
inline
AccessorIterator<A, isConst>::AccessorIterator
(
   const typename AccessorIterator<A, isConst>::Accessor& accessor,
   const size_t index
)
:  accessor_(accessor),
   index_(index)
{}

template<class A, bool isConst>
template<bool isConstLocal>
inline bool
AccessorIterator<A, isConst>::operator==
(
   const AccessorIterator<A, isConstLocal>& it
) const
{
   cgp2d_ASSERT(it.accessor_ == accessor_);
   return it.index_ == index_;
}

template<class A, bool isConst>
template<bool isConstLocal>
inline bool
AccessorIterator<A, isConst>::operator!=
(
   const AccessorIterator<A, isConstLocal>& it
) const
{
   cgp2d_ASSERT(it.accessor_ == accessor_);
   return it.index_ != index_;
}

template<class A, bool isConst>
template<bool isConstLocal>
inline bool
AccessorIterator<A, isConst>::operator<=
(
   const AccessorIterator<A, isConstLocal>& it
) const
{
   cgp2d_ASSERT(it.accessor_ == accessor_);
   return index_ <= it.index_;
}

template<class A, bool isConst>
template<bool isConstLocal>
inline bool
AccessorIterator<A, isConst>::operator>=
(
   const AccessorIterator<A, isConstLocal>& it
) const
{
   cgp2d_ASSERT(it.accessor_ == accessor_);
   return index_ >= it.index_;
}

template<class A, bool isConst>
template<bool isConstLocal>
inline bool
AccessorIterator<A, isConst>::operator<
(
   const AccessorIterator<A, isConstLocal>& it
) const
{
   cgp2d_ASSERT(it.accessor_ == accessor_);
   return index_ < it.index_;
}

template<class A, bool isConst>
template<bool isConstLocal>
inline bool
AccessorIterator<A, isConst>::operator>
(
   const AccessorIterator<A, isConstLocal>& it
) const
{
   return index_ > it.index_;
}

template<class A, bool isConst>
inline typename AccessorIterator<A, isConst>::pointer
AccessorIterator<A, isConst>::operator->()
{
   cgp2d_ASSERT(index_ < accessor_.size());
   return &accessor_[index_]; // whether this works depends on the accessor
}

template<class A, bool isConst>
inline const typename AccessorIterator<A, isConst>::value_type*
AccessorIterator<A, isConst>::operator->() const
{
   return &accessor_[index_]; // whether this works depends on the accessor
}

template<class A, bool isConst>
inline const typename AccessorIterator<A, isConst>::value_type&
AccessorIterator<A, isConst>::operator*() const
{
   return accessor_[index_];
}

template<class A, bool isConst>
inline typename AccessorIterator<A, isConst>::reference
AccessorIterator<A, isConst>::operator*()
{
   return accessor_[index_];
}

template<class A, bool isConst>
inline const typename AccessorIterator<A, isConst>::value_type&
AccessorIterator<A, isConst>::operator[]
(
   const size_t j
) const
{
   return accessor_[index_ + j];
}

template<class A, bool isConst>
inline typename AccessorIterator<A, isConst>::reference
AccessorIterator<A, isConst>::operator[]
(
   const size_t j
)
{
   return accessor_[index_ + j];
}

template<class A, bool isConst>
inline AccessorIterator<A, isConst>&
AccessorIterator<A, isConst>::operator++()
{
   if(index_ < accessor_.size()) {
      ++index_;
   }
   testInvariant();
   return *this;
}

template<class A, bool isConst>
inline AccessorIterator<A, isConst>
AccessorIterator<A, isConst>::operator++(int)
{
   if(index_ < accessor_.size()) {
      AccessorIterator it = *this; // copy
      ++index_;
      testInvariant();
      return it;
   }
   else {
      return *this;
   }
}

template<class A, bool isConst>
inline AccessorIterator<A, isConst>&
AccessorIterator<A, isConst>::operator--()
{
   --index_;
   testInvariant();
   return *this;
}

template<class A, bool isConst>
inline AccessorIterator<A, isConst>
AccessorIterator<A, isConst>::operator--(int)
{
   cgp2d_ASSERT(index_ > 0);
   AccessorIterator it = *this; // copy
   --index_;
   testInvariant();
   return it;
}

template<class A, bool isConst>
inline AccessorIterator<A, isConst>
AccessorIterator<A, isConst>::operator+=
(
   const size_t j
)
{
   if(index_ + j <= accessor_.size()) {
      index_ += j;
   }
   else {
      index_ = accessor_.size();
   }
   testInvariant();
   return *this;
}

template<class A, bool isConst>
inline AccessorIterator<A, isConst>
AccessorIterator<A, isConst>::operator-=
(
   const size_t j
)
{
   index_ -= j;
   testInvariant();
   return *this;
}

template<class A, bool isConst>
inline void
AccessorIterator<A, isConst>::testInvariant() const
{
}

template<class A, bool isConst>
inline AccessorIterator<A, isConst>
AccessorIterator<A, isConst>::operator+
(
   const size_t j
) const
{
   AccessorIterator<A, isConst> it = *this; // copy
   it += j;
   return it;
}

template<class A, bool isConst>
inline AccessorIterator<A, isConst>
AccessorIterator<A, isConst>::operator-
(
   const size_t j
) const
{
   AccessorIterator<A, isConst> it = *this; // copy
   it -= j;
   return it;
}

template<class A, bool isConst>
template<bool isConstLocal>
inline size_t
AccessorIterator<A, isConst>::operator-
(
   const AccessorIterator<A, isConstLocal>& it
) const
{
   //gcc 4.6 bugfix
   #if __GNUC__ == 4 && __GNUC_MINOR__ >= 6
   typedef std::ptrdiff_t difference_type;
   #else
   typedef ptrdiff_t difference_type;
   #endif
   return static_cast<difference_type>(index_) - static_cast<difference_type>(it.index_);
}

// implementation of binary arithmetic operators

template<class A, bool isConst>
inline AccessorIterator<A, isConst>
operator+
(
   const size_t j,
   const AccessorIterator<A, isConst>& it
)
{
   return it + j;
}

/// \endcond

} // namespace cgp2d

#endif // #ifndef CGP2D_ACCESSOR_ITERATOR

