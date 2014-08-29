#ifndef STRIDED_PTR_RANGE_HXX_
#define STRIDED_PTR_RANGE_HXX_

#include <vigra/range/tiny_vector_size.hxx>

namespace vigra
{

namespace strided_ptr_range_detail
{

template <class T, class Shape, unsigned int N>
struct StridedPtrRangeImplCommon
{
    enum{ dimension = N };
    typedef unsigned long size_type;
    
    // enum{ n = dimension - 1 };
    
    T* begin;
    T* end;
    Shape strides;
    Shape shape;
    
    StridedPtrRangeImplCommon()
    {
        begin = 0;
        end = 0;
    }
    StridedPtrRangeImplCommon(T* begin, const Shape& strides, const Shape& shape)
      : begin(begin),
        end(begin + shape[dimension - 1] * strides[dimension - 1]),
        strides(strides),
        shape(shape)
    {}
    
    bool empty() const
    {
        return begin == end;
    }

    size_type length() const
    {
        return (end - begin) / strides[dimension - 1];
    }
    
    void pop_front(size_type i = 1)
    {
        static const int n = dimension-1;
        pop_front_dim<n>(i);
    }

    template <unsigned int n>
    void pop_front_dim(size_type i)
    {
        begin += i * strides[n];
        shape[n] -= i;
        if(n != dimension - 1)
            end += i * strides[n];
    }

    void take_front(size_type i)
    {
        static const int n = dimension - 1;
        take_front_dim<n>(i);
    }
    
    template <unsigned int n>
    void take_front_dim(size_type i)
    {
        shape[n] = i;
        if(n == dimension - 1)
            end = begin + i * strides[n];
    }
    
    //template <unsigned int n>
    void pop_back()
    {
        static const int n = dimension - 1;
        --shape[n];
        if(n == dimension - 1)
            end -= strides[n];
    }

    template <unsigned int n>
    void pop_back(size_type i)
    {
        shape[n] -= i;
        if(n == dimension - 1)
            end -= i * strides[n];
    }
    //template <unsigned int n>
    void take_back(size_type i)
    {
        static const int n = dimension - 1;
        begin += (shape[n] - i) * strides[n];
        if(n != dimension - 1)
            end += (shape[n] - i) * strides[n];
    }
    /*
    void multi_slice(const Shape& new_begin, const Shape& new_end)
    {
        auto shape_it = shape.begin();
        auto strides_it = strides.begin();
        auto new_begin_it = new_begin.begin();
        auto new_end_it = new_end.begin();
        
        for( ; shape_it != shape.end(); ++shape_it, ++new_begin_it, ++new_end_it, ++strides_it)
        {
            *shape_it = *new_end_it - *new_begin_it;
            begin += *strides_it * *new_begin_it;
        }
        end = begin + *(shape.end() - 1) * *(strides.end() - 1);
    }
    */
};

template <class T, class Shape, unsigned int N>
struct StridedPtrRangeImpl
  : StridedPtrRangeImplCommon<T, Shape, N>
{
    typedef StridedPtrRangeImplCommon<T, Shape, N> base_type;
    
    typedef StridedPtrRangeImpl<T, Shape, N - 1> value_type;
    
//    using base_type::base_type;
    StridedPtrRangeImpl()
    : base_type()
    {}

    StridedPtrRangeImpl(T* begin, const Shape& strides, const Shape& shape)
    : base_type(begin, strides, shape)
    {}

/*    
    value_type&& take_front_before(value_type&& r) const
    {
        r.take_front(this->shape[N - 1]);
        return std::move(r);
    }
*/    
    value_type front() const
    {
        return value_type(this->begin, this->strides, this->shape);
    }
    value_type back() const
    {
        return value_type(this->end - this->strides[N - 1], this->strides, this->shape);
    }
};
template <class T, class Shape>
struct StridedPtrRangeImpl<T, Shape, 1>
  : StridedPtrRangeImplCommon<T, Shape, 1>
{
    typedef StridedPtrRangeImplCommon<T, Shape, 1> base_type;
    
    typedef T& value_type;
    
    //using base_type::base_type;
    StridedPtrRangeImpl()
    : base_type()
    {}

    StridedPtrRangeImpl(T* begin, const Shape& strides, const Shape& shape)
    : base_type(begin, strides, shape)
    {}

    value_type front() const
    {
        return *this->begin;
    }
    value_type back() const
    {
        return *(this->end - 1);
    }
};
/*
template <class Shape>
Shape multiply_shape(Shape shape, size_t multiplier)
{
    for(typename Shape::iterator pos = shape.begin(); pos != shape.end(); ++pos)
        *pos *= multiplier;
    return shape;
}
*/

} // strided_ptr_range_detail

template <class T, class Shape>
class StridedPtrRange
  : strided_ptr_range_detail::StridedPtrRangeImpl<T, Shape, TinyVectorSize<Shape>::value>
{
private:
    typedef strided_ptr_range_detail::StridedPtrRangeImpl<T, Shape, TinyVectorSize<Shape>::value> base_type;
public:
    using typename base_type::value_type;
    //using base_type::dimension;
    enum { dimension = base_type::dimension };
    using typename base_type::size_type;

    //using base_type::base_type;
    StridedPtrRange()
    : base_type()
    {}

    StridedPtrRange(T* begin, const Shape& strides, const Shape& shape)
    : base_type(begin, strides, shape)
    {}

    /*
    StridedPtrRange(T* begin, T* end, const Shape& shape)
      : base_type(begin, end, strided_ptr_range_detail::multiply_shape(shape, sizeof(T)))
    {
    }
    */

    using base_type::empty;
    using base_type::length;

    using base_type::front;
    using base_type::pop_front;
    using base_type::pop_front_dim;
    using base_type::take_front;
    using base_type::take_front_dim;

    using base_type::back;
    using base_type::pop_back;
    using base_type::take_back;

    //using base_type::multi_slice;
};

}

#endif

