#ifndef MULTI_AT_HXX_
#define MULTI_AT_HXX_

#include <vigra/range/multi_range_traits.hxx>

#include <utility>

namespace vigra
{

namespace multi_at_detail
{

template <unsigned int N>
struct multi_at_impl
{
    template <class Range, class Shape>
    static typename multi_range_traits::value_type<Range>::type make(Range range, const Shape& shape)
    {
        range.pop_front(shape[N - 1]);
        return multi_at_impl<N - 1>::make(range.front(), shape);
    }
};
template <>
struct multi_at_impl<0>
{
    template <class Value, class Shape>
    static Value make(Value&& val, const Shape& shape)
    {
        //return std::forward<Value>(val);
        return val; // FIXME: not sure whether this will work if the range returns an rvalue reference
    }
};

}

template <class Range, class Shape>
auto multi_at(Range range, Shape shape) -> decltype(multi_at_detail::multi_at_impl<Range::dimension>::make(range, shape))
{
    return multi_at_detail::multi_at_impl<Range::dimension>::make(range, shape);
}

}

#endif

