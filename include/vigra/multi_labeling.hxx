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

namespace vigra{

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

} // namespace lemon_graph

/********************************************************/
/*                                                      */
/*                     labelMultiArray                  */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a MultiArray with arbitrary many dimensions.

    <b> Declaration:</b>

    \code
    namespace vigra {

        template <unsigned int N, class T, class S1,
                                  class Label, class S2,
                  class EqualityFunctor = std::equal_to<T> >
        Label 
        labelMultiArray(MultiArrayView<N, T, S1> const & data,
                        MultiArrayView<N, Label, S2> labels,
                        NeighborhoodType neighborhood = DirectNeighborhood,
                        EqualityFunctor equal = std::equal_to<T>())

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

    Return:  the number of regions found (= highest region label, because labeling starts at 1)

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

/********************************************************/
/*                                                      */
/*           labelMultiArrayWithBackground              */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a MultiArray with arbitrary many dimensions,
     excluding the background from labeling.

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
