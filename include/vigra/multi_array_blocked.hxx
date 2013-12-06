/************************************************************************/
/*                                                                      */
/*     Copyright 2012-2014 by Ullrich Koethe and Thorben Kroeger        */
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

#ifndef VIGRA_MULTI_ARRAY_BLOCKED_HXX
#define VIGRA_MULTI_ARRAY_BLOCKED_HXX

#include "metaprogramming.hxx"
#include "multi_array.hxx"

namespace vigra {

template <class T>
struct BlockedMemory;

template <unsigned int N>
struct BlockShape;

template <>
struct BlockShape<1>
{
    static const unsigned int bits = 18;
    static const unsigned int value = 1 << bits;
    static const unsigned int mask = value -1;
};

template <>
struct BlockShape<2>
{
    static const unsigned int bits = 9;
    static const unsigned int value = 1 << bits;
    static const unsigned int mask = value -1;
};

template <>
struct BlockShape<3>
{
    static const unsigned int bits = 6;
    static const unsigned int value = 1 << bits;
    static const unsigned int mask = value -1;
    
    static std::size_t unmangle(Shape3 const & s, Shape3 & b)
    {
        typedef std::size_t UI;
        b[0] = (UI)s[0] >> bits;
        UI o = (UI)s[0] & mask;
        b[1] = (UI)s[1] >> bits;
        o += ((UI)s[1] & mask) << bits;
        b[2] = (UI)s[2] >> bits;
        o += ((UI)s[2] & mask) << (2*bits);
        return o;
        
    }
};

static int globalCount = 0;

template <class Shape>
Shape outerShape(Shape shape)
{
    static const int N = Shape::static_size;
    for(int k=0; k<N; ++k)
        shape[k] = (shape[k] + BlockShape<N>::mask) >> BlockShape<N>::bits;
    return shape;
}

template <unsigned int N, class T>
class BlockedArray
{
  public:
    typedef MultiArray<N, MultiArray<N, T> > BlockStorage;
    typedef typename BlockStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    BlockedArray(shape_type const & shape)
    : shape_(shape),
      inner_shape_(BlockShape<N>::value),
      outer_array_(outerShape(shape))
    {
        auto i = outer_array_.begin(), 
             end = outer_array_.end();
        for(; i != end; ++i)
            i->reshape(inner_shape_);
    }
    
    reference operator[](shape_type const & p)
    {
        return *ptr(p);
    }
    
    pointer ptr(shape_type const & p)
    {
        // if(!isInside(p))
            // return 0;
        
        shape_type block(SkipInitialization);
        std::size_t offset = BlockShape<N>::unmangle(p, block);
        // std::cerr << "      " << outer_array_[block].data() << " " << offset 
        // << " " << (outer_array_[block].data() + offset) << "\n";
        return outer_array_[block].data() + offset;
    }
    
    shape_type const & shape() const
    {
        return shape_;
    }
    
    bool isInside (shape_type const & p) const
    {
        for(int d=0; d<N; ++d)
            if(p[d] < 0 || p[d] >= shape_[d])
                return false;
        return true;
    }
    
    shape_type shape_, inner_shape_;
    BlockStorage outer_array_;
};

template <class T, class NEXT>
class CoupledHandle<BlockedMemory<T>, NEXT>
: public NEXT
{
public:
    typedef NEXT                                  base_type;
    typedef CoupledHandle<BlockedMemory<T>, NEXT> self_type;
    
    static const int index =                      NEXT::index + 1;    // index of this member of the chain
    static const unsigned int dimensions =        NEXT::dimensions;

    typedef BlockedArray<dimensions, T>           array_type;
    typedef BlockShape<dimensions>                block_shape;
    typedef T                                     value_type;
    typedef value_type *                          pointer;
    typedef value_type const *                    const_pointer;
    typedef value_type &                          reference;
    typedef value_type const &                    const_reference;
    typedef typename base_type::shape_type        shape_type;
    
    typedef T* (*SetPointerFct)(shape_type const & p, array_type * array);
    
    static SetPointerFct setPointer;
    
    static T* setPointerFct(shape_type const & p, array_type * array)
    {
        return array->ptr(p);
    }

    CoupledHandle()
    : base_type(),
      pointer_(), 
      array_()
    {}

    CoupledHandle(array_type & array, NEXT const & next)
    : base_type(next),
      pointer_((*setPointer)(shape_type(), &array)), 
      strides_(detail::defaultStride(shape_type(block_shape::value))),
      array_(&array)
    {}

    inline void incDim(int dim) 
    {
        base_type::incDim(dim);
        if((point()[dim] & block_shape::mask) == 0)
        {
            if(point()[dim] < shape()[dim])
                pointer_ = (*setPointer)(point(), array_);
        }
        else
        {
            pointer_ += strides_[dim];
        }
    }

    // inline void decDim(int dim) 
    // {
        // view_.unsafePtr() -= strides_[dim];
        // base_type::decDim(dim);
    // }

    inline void addDim(int dim, MultiArrayIndex d) 
    {
        base_type::addDim(dim, d);
        pointer_ = (*setPointer)(point(), array_);
        // if(dim == 0)
            // std::cerr << point() << " " << pointer_ << " " << d << "\n";
    }

    // inline void add(shape_type const & d) 
    // {
        // view_.unsafePtr() += dot(d, strides_);
        // base_type::add(d);
    // }
    
    // template<int DIMENSION>
    // inline void increment() 
    // {
        // view_.unsafePtr() += strides_[DIMENSION];
        // base_type::template increment<DIMENSION>();
    // }
    
    // template<int DIMENSION>
    // inline void decrement() 
    // {
        // view_.unsafePtr() -= strides_[DIMENSION];
        // base_type::template decrement<DIMENSION>();
    // }
    
    // // TODO: test if making the above a default case of the this hurts performance
    // template<int DIMENSION>
    // inline void increment(MultiArrayIndex offset) 
    // {
        // view_.unsafePtr() += offset*strides_[DIMENSION];
        // base_type::template increment<DIMENSION>(offset);
    // }
    
    // template<int DIMENSION>
    // inline void decrement(MultiArrayIndex offset) 
    // {
        // view_.unsafePtr() -= offset*strides_[DIMENSION];
        // base_type::template decrement<DIMENSION>(offset);
    // }
    
    // void restrictToSubarray(shape_type const & start, shape_type const & end)
    // {
        // view_.unsafePtr() += dot(start, strides_);
        // base_type::restrictToSubarray(start, end);
    // }

    // ptr access
    reference operator*()
    {
        return *pointer_;
    }

    const_reference operator*() const
    {
        return *pointer_;
    }

    pointer operator->()
    {
        return pointer_;
    }

    const_pointer operator->() const
    {
        return pointer_;
    }

    pointer ptr()
    {
        return pointer_;
    }

    const_pointer ptr() const
    {
        return pointer_;
    }

    pointer pointer_;
    shape_type strides_;
    array_type * array_;
};

template <class T, class NEXT>
typename CoupledHandle<BlockedMemory<T>, NEXT>::SetPointerFct 
CoupledHandle<BlockedMemory<T>, NEXT>::setPointer = &CoupledHandle<BlockedMemory<T>, NEXT>::setPointerFct;


#if 0
template <class T, class NEXT>
class CoupledHandle
: public NEXT
{
public:
    typedef NEXT                            base_type;
    typedef CoupledHandle<T, NEXT>          self_type;
    
    static const int index =                NEXT::index + 1;    // index of this member of the chain
    static const unsigned int dimensions =  NEXT::dimensions;

    typedef T                               value_type;
    typedef T *                             pointer;
    typedef T const *                       const_pointer;
    typedef T &                             reference;
    typedef T const &                       const_reference;
    typedef typename base_type::shape_type  shape_type;

    CoupledHandle()
    : base_type(),
      pointer_(), 
      strides_()
    {}

    CoupledHandle(const_pointer p, shape_type const & strides, NEXT const & next)
    : base_type(next),
      pointer_(const_cast<pointer>(p)), 
      strides_(strides)
    {}

    template <class Stride>
    CoupledHandle(MultiArrayView<dimensions, T, Stride> const & v, NEXT const & next)
    : base_type(next),
      pointer_(const_cast<pointer>(v.data())), 
      strides_(v.stride())
    {
        vigra_precondition(v.shape() == this->shape(), "createCoupledIterator(): shape mismatch.");
    }

    inline void incDim(int dim) 
    {
        pointer_ += strides_[dim];
        base_type::incDim(dim);
    }

    inline void decDim(int dim) 
    {
        pointer_ -= strides_[dim];
        base_type::decDim(dim);
    }

    inline void addDim(int dim, MultiArrayIndex d) 
    {
        pointer_ += d*strides_[dim];
        base_type::addDim(dim, d);
    }

    inline void add(shape_type const & d) 
    {
        pointer_ += dot(d, strides_);
        base_type::add(d);
    }
    
    template<int DIMENSION>
    inline void increment() 
    {
        pointer_ += strides_[DIMENSION];
        base_type::template increment<DIMENSION>();
    }
    
    template<int DIMENSION>
    inline void decrement() 
    {
        pointer_ -= strides_[DIMENSION];
        base_type::template decrement<DIMENSION>();
    }
    
    // TODO: test if making the above a default case of the this hurts performance
    template<int DIMENSION>
    inline void increment(MultiArrayIndex offset) 
    {
        pointer_ += offset*strides_[DIMENSION];
        base_type::template increment<DIMENSION>(offset);
    }
    
    template<int DIMENSION>
    inline void decrement(MultiArrayIndex offset) 
    {
        pointer_ -= offset*strides_[DIMENSION];
        base_type::template decrement<DIMENSION>(offset);
    }
    
    void restrictToSubarray(shape_type const & start, shape_type const & end)
    {
        pointer_ += dot(start, strides_);
        base_type::restrictToSubarray(start, end);
    }

    // ptr access
    reference operator*()
    {
        return *pointer_;
    }

    const_reference operator*() const
    {
        return *pointer_;
    }

    pointer operator->()
    {
        return pointer_;
    }

    const_pointer operator->() const
    {
        return pointer_;
    }

    pointer ptr()
    {
        return pointer_;
    }

    const_pointer ptr() const
    {
        return pointer_;
    }

    shape_type const & strides() const
    {
        return strides_;
    }

    pointer pointer_;
    shape_type strides_;
};

template <unsigned int N, class T, class TAIL>
struct ComposeCoupledHandle<N, TypeList<T, TAIL> >
{
    typedef typename ComposeCoupledHandle<N, TAIL>::type  BaseType;
    typedef typename MultiArrayShape<N>::type             shape_type;
    typedef CoupledHandle<T, BaseType>                    type;
    
    template <class S>
    type exec(MultiArrayView<N, T, S> const & m, 
              shape_type const & start, shape_type const & end,
              BaseType const & base)
    {
        return type(m.subarray(start, end).data(), m.stride(), base);
    }
    
    template <class S>
    type exec(MultiArrayView<N, T, S> const & m, BaseType const & base)
    {
        return type(m.data(), m.stride(), base);
    }
};

template <unsigned int N>
struct ComposeCoupledHandle<N, void>
{
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef CoupledHandle<shape_type, void>    type;
    
    type exec(shape_type const & shape)
    {
        return type(shape);
    }
    
    type exec(shape_type const & start, shape_type const & end)
    {
        return type(end-start);
    }
};

template <unsigned int N, class T1=void, class T2=void, class T3=void, class T4=void, class T5=void>
struct CoupledHandleType
{
    // reverse the order to get the desired index order
    typedef typename MakeTypeList<T5, T4, T3, T2, T1>::type TypeList;
    typedef typename ComposeCoupledHandle<N, TypeList>::type type;
};

template <unsigned int N, class T1, class T2, class T3, class T4, class T5>
struct CoupledHandleType<N, Multiband<T1>, T2, T3, T4, T5>
{
    // reverse the order to get the desired index order
    typedef typename MakeTypeList<T5, T4, T3, T2, Multiband<T1> >::type TypeList;
    typedef typename ComposeCoupledHandle<N-1, TypeList>::type type;
};

#endif

} // namespace vigra

#endif /* VIGRA_MULTI_ARRAY_BLOCKED_HXX */
