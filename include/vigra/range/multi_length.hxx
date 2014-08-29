#ifndef MULTI_LENGTH_HXX_
#define MULTI_LENGTH_HXX_

#include <vigra/range/multi_range_traits.hxx>
#include <vigra/range/range.hxx>

namespace vigra
{

namespace multi_length_detail
{

template <unsigned int N>
struct multi_length_impl
{
    template <class Range, class Shape>
    static void make(const Range& r, Shape& shape)
    {
        shape[N - 1] = length(r);
        if(!r.empty())
            multi_length_impl<N - 1>::make(r.front(), shape);
        else
            multi_length_impl<N - 1>::make(shape);
    }
    template <class Shape>
    static void make(Shape& shape)
    {
        shape[N - 1] = 0;
        multi_length_impl<N - 1>::make(shape);
    }
};

template <>
struct multi_length_impl<0>
{
    template <class Range, class Shape>
    static void make(const Range&, Shape&)
    {}
    template <class Shape>
    static void make(Shape& shape)
    {}
};

}

template <class Range>
typename multi_range_traits::multi_size_type<Range>::type
multi_length(const Range& range)
{
    using namespace multi_length_detail;
    typename multi_range_traits::multi_size_type<Range>::type shape;
    multi_length_impl<Range::dimension>::make(range, shape);
    return shape;
}

}

#endif

