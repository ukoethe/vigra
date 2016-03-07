/************************************************************************/
/*                                                                      */
/*               Copyright 2012-2013 by Ullrich Koethe                  */
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

#ifndef VIGRA_MULTI_LABELING_HXX
#define VIGRA_MULTI_LABELING_HXX

#include "multi_array.hxx"
#include "multi_gridgraph.hxx"
#include "union_find.hxx"
#include "any.hxx"

namespace vigra{

namespace labeling_equality{

struct Yes
{
    char d[100];
};
struct No
{
    char d[1];
};

template <size_t>
struct SizeToType;
template <>
struct SizeToType<sizeof(Yes)>
{
    typedef VigraTrueType type;
};
template <>
struct SizeToType<sizeof(No)>
{
    typedef VigraFalseType type;
};

template <class Equal>
class TakesThreeArguments
{
public:
    template <class T>
    static Yes check(typename T::WithDiffTag*);
    template <class T>
    static No check(...);

    typedef typename SizeToType<sizeof(check<Equal>(0))>::type type;
    static const unsigned int value = type::asBool;
};

template <class Equal, class Data, class Shape>
bool callEqualImpl(Equal& equal, const Data& u_data, const Data& v_data, const Shape& diff, VigraTrueType)
{
    return equal(u_data, v_data, diff);
}
template <class Equal, class Data, class Shape>
bool callEqualImpl(Equal& equal, const Data& u_data, const Data& v_data, const Shape&, VigraFalseType)
{
    return equal(u_data, v_data);
}

template< class Equal, class Data, class Shape>
bool callEqual(Equal& equal, const Data& u_data, const Data& v_data, const Shape& diff)
{
    return callEqualImpl(equal, u_data, v_data, diff, typename TakesThreeArguments<Equal>::type());
}

}


/** \addtogroup Labeling
*/
//@{

namespace lemon_graph {

template <class Graph, class T1Map, class T2Map, class Equal>
typename T2Map::value_type
labelGraph(Graph const & g,
           T1Map const & data,
           T2Map & labels,
           Equal const & equal)
{
    typedef typename Graph::NodeIt        graph_scanner;
    typedef typename Graph::OutBackArcIt  neighbor_iterator;
    typedef typename T2Map::value_type    LabelType;

    vigra::UnionFindArray<LabelType>  regions;

    // pass 1: find connected components
    for (graph_scanner node(g); node != INVALID; ++node)
    {
        typename T1Map::value_type center = data[*node];

        // define tentative label for current node
        LabelType currentIndex = regions.nextFreeIndex();

        for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
        {
            // merge regions if colors are equal
            if(equal(center, data[g.target(*arc)]))
            {
                currentIndex = regions.makeUnion(labels[g.target(*arc)], currentIndex);
            }
        }
        // set label of current node
        labels[*node] = regions.finalizeIndex(currentIndex);
    }

    LabelType count = regions.makeContiguous();

    // pass 2: make component labels contiguous
    for (graph_scanner node(g); node != INVALID; ++node)
    {
        labels[*node] = regions.findLabel(labels[*node]);
    }
    return count;
}

template <unsigned int N, class DirectedTag, class T1Map, class T2Map, class Equal>
typename T2Map::value_type
labelGraph(GridGraph<N, DirectedTag> const & g,
           T1Map const & data,
           T2Map & labels,
           Equal const & equal)
{
    typedef GridGraph<N, DirectedTag>     Graph;
    typedef typename Graph::NodeIt        graph_scanner;
    typedef typename Graph::OutBackArcIt  neighbor_iterator;
    typedef typename T2Map::value_type    LabelType;
    typedef typename Graph::shape_type    Shape;

    vigra::UnionFindArray<LabelType>  regions;

    // pass 1: find connected components
    for (graph_scanner node(g); node != INVALID; ++node)
    {
        typename T1Map::value_type center = data[*node];

        // define tentative label for current node
        LabelType currentIndex = regions.nextFreeIndex();

        for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
        {
            Shape diff = g.neighborOffset(arc.neighborIndex());
            // merge regions if colors are equal
            if(labeling_equality::callEqual(equal, center, data[g.target(*arc)], diff))
            {
                currentIndex = regions.makeUnion(labels[g.target(*arc)], currentIndex);
            }
        }
        // set label of current node
        labels[*node] = regions.finalizeIndex(currentIndex);
    }

    LabelType count = regions.makeContiguous();

    // pass 2: make component labels contiguous
    for (graph_scanner node(g); node != INVALID; ++node)
    {
        labels[*node] = regions.findLabel(labels[*node]);
    }
    return count;
}


template <class Graph, class T1Map, class T2Map, class Equal>
typename T2Map::value_type
labelGraphWithBackground(Graph const & g,
                         T1Map const & data,
                         T2Map & labels,
                         typename T1Map::value_type backgroundValue,
                         Equal const & equal)
{
    typedef typename Graph::NodeIt        graph_scanner;
    typedef typename Graph::OutBackArcIt  neighbor_iterator;
    typedef typename T2Map::value_type    LabelType;

    vigra::UnionFindArray<LabelType>  regions;

    // pass 1: find connected components
    for (graph_scanner node(g); node != INVALID; ++node)
    {
        typename T1Map::value_type center = data[*node];

        // background always gets label zero
        if(equal(center, backgroundValue))
        {
            labels[*node] = 0;
            continue;
        }

        // define tentative label for current node
        LabelType currentIndex = regions.nextFreeIndex();

        for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
        {
            // merge regions if colors are equal
            if(equal(center, data[g.target(*arc)]))
            {
                currentIndex = regions.makeUnion(labels[g.target(*arc)], currentIndex);
            }
        }
        // set label of current node
        labels[*node] = regions.finalizeIndex(currentIndex);
    }

    LabelType count = regions.makeContiguous();

    // pass 2: make component labels contiguous
    for (graph_scanner node(g); node != INVALID; ++node)
    {
        labels[*node] = regions.findLabel(labels[*node]);
    }
    return count;
}

template <unsigned int N, class DirectedTag, class T1Map, class T2Map, class Equal>
typename T2Map::value_type
labelGraphWithBackground(GridGraph<N, DirectedTag> const & g,
                         T1Map const & data,
                         T2Map & labels,
                         typename T1Map::value_type backgroundValue,
                         Equal const & equal)
{
    typedef GridGraph<N, DirectedTag>     Graph;
    typedef typename Graph::NodeIt        graph_scanner;
    typedef typename Graph::OutBackArcIt  neighbor_iterator;
    typedef typename T2Map::value_type    LabelType;
    typedef typename Graph::shape_type    Shape;

    vigra::UnionFindArray<LabelType>  regions;

    // pass 1: find connected components
    for (graph_scanner node(g); node != INVALID; ++node)
    {
        typename T1Map::value_type center = data[*node];

        // background always gets label zero
        if(labeling_equality::callEqual(equal, center, backgroundValue, Shape()))
        {
            labels[*node] = 0;
            continue;
        }

        // define tentative label for current node
        LabelType currentIndex = regions.nextFreeIndex();

        for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
        {
            // merge regions if colors are equal
            Shape diff = g.neighborOffset(arc.neighborIndex());
            if(labeling_equality::callEqual(equal, center, data[g.target(*arc)], diff))
            {
                currentIndex = regions.makeUnion(labels[g.target(*arc)], currentIndex);
            }
        }
        // set label of current node
        labels[*node] = regions.finalizeIndex(currentIndex);
    }

    LabelType count = regions.makeContiguous();

    // pass 2: make component labels contiguous
    for (graph_scanner node(g); node != INVALID; ++node)
    {
        labels[*node] = regions.findLabel(labels[*node]);
    }
    return count;
}


} // namespace lemon_graph

    /** \brief Option object for labelMultiArray().
    */
class LabelOptions
{
    Any background_value_;
    NeighborhoodType neighborhood_;

  public:

        /** \brief Set default options.
        */
    LabelOptions()
    : neighborhood_(DirectNeighborhood)
    {}

        /** \brief Choose direct or indirect neighborhood.

            Default: <tt>DirectNeighborhood</tt>
        */
    LabelOptions & neighborhood(NeighborhoodType n)
    {
        neighborhood_ = n;
        return *this;
    }

        /** \brief Query the neighborhood type (direct or indirect).
        */
    NeighborhoodType getNeighborhood() const
    {
        return neighborhood_;
    }

        /** \brief Set background value.

            If specified, labelMultiArray() will internally call
            labelMultiArrayWithBackground() with the given value
            considered as background and thus ignored. If no
            background value is set, the array gets labeled completely.
            Note that the type <tt>T</tt> must correspond to the element
            type of the data array to be labeled.

            Default: don't ignore any value.
        */
    template <class T>
    LabelOptions & ignoreBackgroundValue(T const & t)
    {
        background_value_ = t;
        return *this;
    }

        /** \brief Check if some background value is to be ignored.
        */
    bool hasBackgroundValue() const
    {
        return bool(background_value_);
    }

        /** \brief Get the background value to be ignored.

            Throws an exception if the stored background value type
            is incompatible to the data array's value type.
        */
    template <class T>
    T getBackgroundValue() const
    {
        if(background_value_.empty())
            return T();
        vigra_precondition(background_value_.template is_readable<T>(),
            "LabelOptions::getBackgroundValue<T>(): stored background value is not convertible to T.");
        return background_value_.template read<T>();
    }
};

/********************************************************/
/*                                                      */
/*                     labelMultiArray                  */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a MultiArray with arbitrary many dimensions.

    See also \ref labelMultiArrayBlockwise() for a parallel version of this algorithm.

    By specifying a background value in the \ref vigra::LabelOptions, this function
    can also realize the behavior of \ref labelMultiArrayWithBackground().

    <b> Declaration:</b>

    \code
    namespace vigra {

        // specify parameters directly
        template <unsigned int N, class T, class S1,
                                  class Label, class S2,
                  class EqualityFunctor = std::equal_to<T> >
        Label
        labelMultiArray(MultiArrayView<N, T, S1> const & data,
                        MultiArrayView<N, Label, S2> labels,
                        NeighborhoodType neighborhood = DirectNeighborhood,
                        EqualityFunctor equal = std::equal_to<T>())

        // specify parameters via LabelOptions
        template <unsigned int N, class T, class S1,
                                  class Label, class S2,
                  class Equal = std::equal<T> >
        Label
        labelMultiArray(MultiArrayView<N, T, S1> const & data,
                        MultiArrayView<N, Label, S2> labels,
                        LabelOptions const & options,
                        Equal equal = std::equal<T>());

    }
    \endcode

    Connected components are defined as regions with uniform values.
    Thus, the value type <tt>T</tt> of the input array \a data either
    must be equality comparable, or an EqualityFunctor must be
    provided that realizes the desired predicate. The destination
    array's value type <tt>Label</tt> should be large enough to hold
    the labels without overflow. Region numbers will form a consecutive
    sequence starting at <b>one</b> and ending with the region number
    returned by the function (inclusive).

    Argument \a neighborhood specifies the type of connectivity used. It can
    take the values <tt>DirectNeighborhood</tt> (which corresponds to
    4-neighborhood in 2D and 6-neighborhood in 3D, default) or
    <tt>IndirectNeighborhood</tt> (which corresponds to
    8-neighborhood in 2D and 26-neighborhood in 3D).

    Return:  the highest region label used

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_labeling.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<3,int> src(Shape3(w,h,d));
    MultiArray<3,int> dest(Shape3(w,h,d));

    // find 6-connected regions
    int max_region_label = labelMultiArray(src, dest);

    // find 26-connected regions
    max_region_label = labelMultiArray(src, dest, IndirectNeighborhood);

    // find 6-connected regions, ignore the background value 0
    // (i.e. call labelMultiArrayWithBackground() internally)
    max_region_label = labelMultiArray(src, dest,
                                       LabelOptions().neighborhood(DirectNeighborhood)
                                                     .ignoreBackgroundValue(0));
    \endcode

    <b> Required Interface:</b>

    \code
    T t;
    t == t                      // T is equality comparable

    EqualityFunctor equal;      // or a suitable functor is supplied
    equal(t, t)
    \endcode
*/
doxygen_overloaded_function(template <...> unsigned int labelMultiArray)

template <unsigned int N, class T, class S1,
                          class Label, class S2,
          class Equal>
inline Label
labelMultiArray(MultiArrayView<N, T, S1> const & data,
                MultiArrayView<N, Label, S2> labels,
                NeighborhoodType neighborhood,
                Equal equal)
{
    vigra_precondition(data.shape() == labels.shape(),
        "labelMultiArray(): shape mismatch between input and output.");

    GridGraph<N, undirected_tag> graph(data.shape(), neighborhood);
    return lemon_graph::labelGraph(graph, data, labels, equal);
}

template <unsigned int N, class T, class S1,
                          class Label, class S2>
inline Label
labelMultiArray(MultiArrayView<N, T, S1> const & data,
                MultiArrayView<N, Label, S2> labels,
                NeighborhoodType neighborhood = DirectNeighborhood)
{
    return labelMultiArray(data, labels, neighborhood, std::equal_to<T>());
}

template <unsigned int N, class T, class S1,
                          class Label, class S2>
inline Label
labelMultiArray(MultiArrayView<N, T, S1> const & data,
                MultiArrayView<N, Label, S2> labels,
                LabelOptions const & options)
{
    if(options.hasBackgroundValue())
        return labelMultiArrayWithBackground(data, labels, options.getNeighborhood(),
                                             options.template getBackgroundValue<T>());
    else
        return labelMultiArray(data, labels, options.getNeighborhood());
}

template <unsigned int N, class T, class S1,
                          class Label, class S2,
          class Equal>
inline Label
labelMultiArray(MultiArrayView<N, T, S1> const & data,
                MultiArrayView<N, Label, S2> labels,
                LabelOptions const & options,
                Equal equal)
{
    if(options.hasBackgroundValue())
        return labelMultiArrayWithBackground(data, labels, options.getNeighborhood(),
                                             options.template getBackgroundValue<T>(),
                                             equal);
    else
        return labelMultiArray(data, labels, options.getNeighborhood(), equal);
}

/********************************************************/
/*                                                      */
/*           labelMultiArrayWithBackground              */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a MultiArray with arbitrary many dimensions,
     excluding the background from labeling.

    From a user's point of view, this function is no longer needed because
    \ref labelMultiArray() can realizes the same behavior when an appropriate
    background value is specified in its \ref vigra::LabelOptions. Similarly,
    \ref labelMultiArrayBlockwise() implements a parallel version of this algorithm.

    <b> Declaration:</b>

    \code
    namespace vigra {

        template <unsigned int N, class T, class S1,
                                  class Label, class S2
                  class Equal = std::equal<T> >
        Label
        labelMultiArrayWithBackground(MultiArrayView<N, T, S1> const & data,
                                      MultiArrayView<N, Label, S2> labels,
                                      NeighborhoodType neighborhood = DirectNeighborhood,
                                      T backgroundValue = T(),
                                      Equal equal = std::equal<T>());

    }
    \endcode

    This function is the same as \ref labelMultiArray(), except for
    the additional parameter \a backgroundValue. Points in the input
    array \a data with this value (which default to zero) are ignored
    during labeling, and  their output label is automatically set to
    zero. Region numbers will be a consecutive sequence starting at
    zero (when background was present) or at one (when no background
    was present) and ending with the region number returned by the
    function (inclusive).

    Return: the number of non-background regions found (= highest region label,
    because background has label 0)

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_labeling.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<3, int> src(Shape3(w,h,d));
    MultiArray<3, int> dest(Shape3(w,h,d));

    // find 6-connected regions, ignoring background value zero (the default)
    int max_region_label = labelMultiArrayWithBackground(src, dest);

    // find 26-connected regions, using -1 as the background value
    max_region_label = labelMultiArrayWithBackground(src, dest, IndirectNeighborhood, -1);
    \endcode

    <b> Required Interface:</b>

    \code
    T t, backgroundValue;
    t == backgroundValue        // T is equality comparable

    EqualityFunctor equal;      // or a suitable functor is supplied
    equal(t, backgroundValue)
    \endcode
*/
doxygen_overloaded_function(template <...> unsigned int labelMultiArrayWithBackground)

template <unsigned int N, class T, class S1,
                          class Label, class S2,
          class Equal>
inline Label
labelMultiArrayWithBackground(MultiArrayView<N, T, S1> const & data,
                              MultiArrayView<N, Label, S2> labels,
                              NeighborhoodType neighborhood,
                              T backgroundValue,
                              Equal equal)
{
    vigra_precondition(data.shape() == labels.shape(),
        "labelMultiArrayWithBackground(): shape mismatch between input and output.");

    GridGraph<N, undirected_tag> graph(data.shape(), neighborhood);
    return lemon_graph::labelGraphWithBackground(graph, data, labels, backgroundValue, equal);
}

template <unsigned int N, class T, class S1,
                          class Label, class S2>
inline Label
labelMultiArrayWithBackground(MultiArrayView<N, T, S1> const & data,
                              MultiArrayView<N, Label, S2> labels,
                              NeighborhoodType neighborhood = DirectNeighborhood,
                              T backgroundValue = T())
{
    return labelMultiArrayWithBackground(data, labels, neighborhood, backgroundValue, std::equal_to<T>());
}

//@}

} // namespace vigra

#endif //VIGRA_MULTI_LABELING_HXX
