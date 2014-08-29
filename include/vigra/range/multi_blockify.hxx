#ifndef MULTI_BLOCKIFY_HXX_
#define MULTI_BLOCKIFY_HXX_

#include <utility>
#include <algorithm>

#include <vigra/range/multi_length.hxx>

namespace vigra
{

namespace multi_blockify_detail
{

template <class Range, unsigned int N, class Shape>
class BlockifiedMultiRangeImplCommon
{
public:
    enum{ dimension = N };
    typedef typename Range::size_type size_type;
    
    Range range;
    Shape block_shape;
    
    BlockifiedMultiRangeImplCommon(Range range, const Shape& block_shape)
      : range(std::move(range)),
        block_shape(block_shape)
    {}
    
    bool empty() const
    {
        return multi_length(range)[N - 1] == 0;
    }
    size_type length() const
    {
        size_type result = multi_length(range)[N - 1] / block_shape[N - 1];
        if(result * block_shape[N - 1] != multi_length(range)[N - 1])
            ++result;
        return result;
    }
    
    void pop_front(size_type i = 1)
    {
        static const int n = dimension - 1;
        pop_front_dim<n>(i);
    }
    
    template <unsigned int n>
    void pop_front_dim(size_type i)
    {
        //cout << "pop_front_dim " << n << endl;
        //cout << "length before: " << multi_length(range) << endl;
        range. template pop_front_dim<n>(std::min<size_type>(multi_length(range)[n], block_shape[n] * i));
        //cout << "length after: " << multi_length(range) << endl;
        //print_range(range);
    }

    template <unsigned int n>
    void take_front(size_type i = 1)
    {
        range. template take_front<n>(std::min<size_type>(multi_length(range)[n], block_shape[n] * i));
    }
    
    template <unsigned int n>
    void pop_back(size_type i = 1)
    {
        size_type last_size = multi_length(range)[n] % (block_shape[n]);
        if(last_size != 0)
        {
            range. template pop_back<n>(last_size);
            range. template pop_back<n>(block_shape[n] * (i - 1));
        }
        else
        {
            range. template pop_back<n>(block_shape[n] * i );
        }
    }

    template <unsigned int n>
    void take_back(size_type i = 1)
    {
        size_type last_size = multi_length(range)[n] % (block_shape[n]);
        if(last_size == 0)
            last_size = block_shape[n];
        range. template take_back<n>(last_size);
    }
};


template <class Range, unsigned int N, class Shape>
class BlockifiedMultiRangeImpl
  : public BlockifiedMultiRangeImplCommon<Range, N, Shape>
{
    // N >= 2
public:    
    typedef BlockifiedMultiRangeImplCommon<Range, N, Shape> base_type;

    using base_type::dimension;
    using typename base_type::size_type;
    typedef BlockifiedMultiRangeImpl<Range, N - 1, Shape> value_type;
    
    //using base_type::base_type;

    BlockifiedMultiRangeImpl(Range const & range, const Shape& block_shape)
    : base_type(range, block_shape)
    {}

    value_type front() const
    {
        value_type result(this->range, this->block_shape);
        result.range.template take_front_dim<N - 1>(std::min<size_type>(multi_length(this->range)[N - 1], this->block_shape[N - 1]));
        //cout << "front, length: " << multi_length(result.range) << endl;
        return result;
    }
    value_type back() const
    {
        value_type result(this->range, this->block_shape);
        size_type last_size = multi_length(this->range)[N - 1] % this->block_shape[N - 1];
        if(last_size == 0)
            last_size = this->block_shape[N - 1];
        result.range.template take_back<N - 1>(last_size);
        return result;
    }
};

template <class Range, class Shape>
class BlockifiedMultiRangeImpl<Range, 1, Shape>
  : public BlockifiedMultiRangeImplCommon<Range, 1, Shape>
{
public: 
    typedef BlockifiedMultiRangeImplCommon<Range, 1, Shape> base_type;
    
    using typename base_type::size_type;
    typedef Range value_type;

    //using base_type::base_type;
    BlockifiedMultiRangeImpl(Range const & range, const Shape& block_shape)
    : base_type(range, block_shape)
    {}
    
    value_type front() const
    {
        Range result = this->range;
        result.template take_front_dim<0>(std::min<size_type>(this->block_shape[0], multi_length(result)[0]));
        //cout << "end front, length: " << multi_length(result) << endl;
        return result;
    }
    value_type back() const
    {
        Range result = this->range;
        size_type last_size = multi_length(this->range)[0] % this->block_shape[0];
        if(last_size == 0)
            last_size = this->block_shape[0];
        result.template take_back<0>(last_size);
        return result;
    }
};

}

template <class Range, class Shape>
class BlockifiedMultiRange
  : private multi_blockify_detail::BlockifiedMultiRangeImpl<Range, Range::dimension, Shape>
{
private:
    typedef multi_blockify_detail::BlockifiedMultiRangeImpl<Range, Range::dimension, Shape> base_type;
public:
    BlockifiedMultiRange(Range const & range, const Shape& block_shape)
    : base_type(range, block_shape)
    {}

    using base_type::value_type;
    //using base_type::dimension;
    enum {dimension = base_type::dimension };
    using base_type::size_type;

    //using base_type::base_type;
    using base_type::empty;
    using base_type::length;
    
    using base_type::front;
    using base_type::pop_front;
    using base_type::pop_front_dim;
    using base_type::take_front;

    using base_type::back;
    using base_type::pop_back;
    using base_type::take_back;
};

template <class Range, class Shape>
BlockifiedMultiRange<Range, Shape> multi_blockify(Range range, Shape shape)
{
    return BlockifiedMultiRange<Range, Shape>(range, shape);
}

}

#endif

