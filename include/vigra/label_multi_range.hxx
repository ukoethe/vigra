#ifndef LABEL_MULTI_RANGE_HXX_
#define LABEL_MULTI_RANGE_HXX_

#include <vigra/multi_labeling.hxx>
#include <vigra/range/multi_range_traits.hxx>
#include <vigra/range/multi_range.hxx>
#include <vigra/range/multi_length.hxx>

#include <type_traits>

namespace vigra
{

namespace label_multi_range_detail
{

template <class Range>
struct RangePropertyMap
{
    Range range;
    RangePropertyMap(Range r)
      : range(std::move(r))
    {}
    typedef typename multi_range_traits::value_type<Range>::type actual_value_type;
    typedef typename std::decay<actual_value_type>::type value_type;
    // typedef typename multi_range_traits::multi_size_type<Range>::type difference_type;
    
    template <class Shape>
    actual_value_type operator[](const Shape& coordinates) const
    {
        return multi_at(range, coordinates);
    }
};

}

template <class DataRange, class LabelRange, class Equal>
typename multi_range_traits::decayed_value_type<LabelRange>::type
labelMultiRange(DataRange data_range, LabelRange label_range,
                NeighborhoodType neighborhood, Equal equal)
{
    enum{ N = DataRange::dimension };

    using namespace label_multi_range_detail;

    RangePropertyMap<DataRange> data_map(data_range);
    RangePropertyMap<LabelRange> label_map(label_range);
    
    typedef typename multi_range_traits::multi_size_type<DataRange>::type Shape;
    Shape shape = multi_length(data_range);
    GridGraph<N, undirected_tag> graph(shape, neighborhood);

    return lemon_graph::labelGraph(graph, data_map, label_map, equal);
}

}

#endif

