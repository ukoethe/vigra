#ifndef VISIT_BORDER_HXX_
#define VISIT_BORDER_HXX_

#include <vigra/range/tiny_vector_size.hxx>
#include <vigra/range/range.hxx>

#include <vigra/error.hxx>
#include <vigra/multi_shape.hxx>

using namespace vigra;

namespace visit_border_detail
{

template <class Range1, class Range2>
void assertNonEmptyEqualLength(const Range1& u, const Range2& v)
{
    vigra_assert(length(u) == length(v), "ranges do not have the same length");
    vigra_assert(!u.empty(), "empty ranges");
}

// let u and v be ranges of equal length n
// for each i in 0 ... n - 1
//     visit(u[i], v[i])
// if u and v are two adjacent hyperplanes, this corresponds to visiting all direct neighbors
template <class Range1, class Range2, class Visitor>
void forEachDirectNeighbor(Range1 u, Range2 v, Visitor visitor)
{
    assertNonEmptyEqualLength(u, v);
    auto zipped = zip(u, v);
    typedef decltype(zipped) zipped_type;
    typedef typename zipped_type::value_type pair_type;
    for_each(zip(u, v), [&visitor](pair_type p)
    {
        visitor(p.first, p.second);
    });
}
// let u and v be ranges of equal length n
// for each i in 0 ... n - 1
//   if i > 0
//     visit(u[i] , v[i - 1])    
//   visit(u[i], v[i])
//   if i < n - 1
//     visit(u[i] , v[i + 1])
// if u and v are two adjacent hyperplanes, this corresponds to visiting all indirect neighbors
template <class Range1, class Range2, class Visitor>
void forEachIndirectNeighbor(Range1 u, Range2 v, Visitor visitor)
{
    assertNonEmptyEqualLength(u, v);
    typedef typename Range1::value_type u_value_type;
    typedef typename Range2::value_type v_value_type;

    visitor(u.front(), v.front());
    Range2 v_after_begin = v;
    v_after_begin.pop_front();
    if(!v_after_begin.empty())
    {
        visitor(u.front(), v_after_begin.front());
        u.pop_front();
        Range1 u_till_before_end = u;
        u_till_before_end.pop_back();
        for_each(u_till_before_end, [&](u_value_type u_val)
        {
            Range2 v_neighbors = v;
            v_neighbors.take_front(3); // TODO
            for_each(v_neighbors, [&u_val, &visitor](v_value_type v_val)
            {
                visitor(u_val, v_val);
            });
            v.pop_front();
        });
        visitor(u.back(), v.front());
        v.pop_front();
        visitor(u.back(), v.front());
    }
}

template <unsigned int N, NeighborhoodType neighborhood, class Range1, class Range2>
struct VisitBorderImpl
{
    template <class Shape, class Visitor>
    static void make(Range1 u, Range2 v, const Shape& difference, Visitor visitor)
    {
        enum{ n = N - 1 };
        typedef typename Range1::value_type Next1;
        typedef typename Range2::value_type Next2;
        vigra_assert(-1 <= difference[n] && difference[n] <= 1, "difference out of bounds");
        auto call_recursively = [&](Next1 u_next, Next2 v_next)
        {
            VisitBorderImpl<N - 1, neighborhood, Next1, Next2>::make(u_next, v_next, difference, visitor);
        };
        if(difference[n] == 0)
        {
            if(neighborhood == DirectNeighborhood)
                forEachDirectNeighbor(u, v, call_recursively);
            else if(neighborhood == IndirectNeighborhood)
                forEachIndirectNeighbor(u, v, call_recursively);
            else
                vigra_fail("unsupported neighborhood");
        }
        else if(difference[n] == -1)
            call_recursively(u.front(), v.back());
        else // difference[n] == 1
            call_recursively(u.back(), v.front());
    }
};

template <NeighborhoodType neighborhood, class Value1, class Value2>
struct VisitBorderImpl<0, neighborhood, Value1, Value2>
{
    template <class Shape, class Visitor>
    static void make(Value1 u_val, Value2 v_val, const Shape& shape, Visitor visitor)
    {
        visitor(u_val, v_val);
    }
};

}

    
template <class Range1, class Range2, class Shape, class Visitor>
void visitBorder(Range1 u, Range2 v, const Shape& difference,
                 NeighborhoodType neighborhood, Visitor visitor)
{
    using namespace visit_border_detail;
    enum{ N = TinyVectorSize<Shape>::value };

    if(neighborhood == DirectNeighborhood)
        VisitBorderImpl<N, DirectNeighborhood, Range1, Range2>::make(u, v, difference, visitor);
    else if(neighborhood == IndirectNeighborhood)
        VisitBorderImpl<N, IndirectNeighborhood, Range1, Range2>::make(u, v, difference, visitor);
    else
        vigra_fail("unsupported neighborhood");
}

#endif

