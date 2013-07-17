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


#ifndef VIGRA_MULTI_LOCALMINMAX_HXX
#define VIGRA_MULTI_LOCALMINMAX_HXX

#include <vector>
#include <functional>
#include "multi_array.hxx"
#include "localminmax.hxx"
#include "multi_gridgraph.hxx"
#include "multi_labeling.hxx"
#include "metaprogramming.hxx"

namespace vigra {

namespace boost_graph {

  // Attempt without LValue propmaps, using only the free functions
  // to access ReadablePropertyMap (input) and WritablePropertyMap (label)
template <class Graph, class T1Map, class T2Map, class Compare>
unsigned int
localMinMaxGraph(Graph const &g, 
                 T1Map const &src,
                 T2Map &dest,
                 typename property_traits<T2Map>::value_type marker,
                 typename property_traits<T1Map const>::value_type threshold,
                 Compare const &compare,
                 bool allowAtBorder = true)
{
    typedef typename graph_traits<Graph>::vertex_iterator graph_scanner;
    typedef typename graph_traits<Graph>::adjacency_iterator neighbor_iterator;

    typedef typename property_traits<T1Map const>::value_type T1;

    graph_scanner node, srcend;
    neighbor_iterator arc, nbend;

    unsigned int count = 0;
    tie(node, srcend) = vertices(g);
    for (; node != srcend; ++node) 
    {
        const T1 current = get(src, *node);

        if (!compare(current, threshold))
            continue;
          
        if(!allowAtBorder && node.atBorder())
            continue;
        
        tie(arc, nbend) = adjacent_vertices(*node, g);
        for (;arc != nbend; ++arc) 
            if (!compare(current, get(src, *arc))) 
                break;
                
        if (arc == nbend)
        {
            put(dest, *node, marker);
            ++count;
        }
    }
    return count;
}

} // namespace boost_graph

namespace lemon_graph { 

template <class Graph, class T1Map, class T2Map, class Compare>
unsigned int
localMinMaxGraph(Graph const &g, 
                 T1Map const &src,
                 T2Map &dest,
                 typename T2Map::value_type marker,
                 typename T1Map::value_type threshold,
                 Compare const &compare,
                 bool allowAtBorder = true)
{
    typedef typename Graph::NodeIt    graph_scanner;
    typedef typename Graph::OutArcIt  neighbor_iterator;

    unsigned int count = 0;
    for (graph_scanner node(g); node != INVALID; ++node) 
    {
        typename T1Map::value_type current = src[*node];

        if (!compare(current, threshold))
            continue;
          
        if(!allowAtBorder && node.atBorder())
            continue;
        
        neighbor_iterator arc(g, node);
        for (; arc != INVALID; ++arc) 
            if (!compare(current, src[g.target(*arc)])) 
                break;
                
        if (arc == INVALID)
        {
            dest[*node] = marker;
            ++count;
        }
    }
    return count;
}

template <class Graph, class T1Map, class T2Map, class Compare, class Equal>
unsigned int
extendedLocalMinMaxGraph(Graph const &g, 
                         T1Map const &src,
                         T2Map &dest,
                         typename T2Map::value_type marker,
                         typename T1Map::value_type threshold,
                         Compare const &compare,
                         Equal const &equal,
                         bool allowAtBorder = true)
{
    typename Graph::template NodeMap<unsigned int> regions(g);
    
    int max_region_label = labelGraph(g, src, regions, equal);
    
    // assume that a region is a extremum until the opposite is proved
    std::vector<unsigned char> isExtremum(max_region_label+1, (unsigned char)1);

    typedef typename Graph::NodeIt    graph_scanner;
    typedef typename Graph::OutArcIt  neighbor_iterator;

    unsigned int count = max_region_label;
    for (graph_scanner node(g); node != INVALID; ++node) 
    {
        unsigned int label = regions[*node];
        
        if(!isExtremum[label])
            continue;
        
        typename T1Map::value_type current = src[*node];

        if (!compare(current, threshold) ||
            (!allowAtBorder && node.atBorder()))
        {
            isExtremum[label] = 0;
            --count;
            continue;
        }
        
        for (neighbor_iterator arc(g, node); arc != INVALID; ++arc) 
        {
            if (label != regions[g.target(*arc)] && compare(src[g.target(*arc)], current)) 
            {
                isExtremum[label] = 0;
                --count;
                break;
            }
        }
    }
    for (graph_scanner node(g); node != INVALID; ++node) 
    {
        if(isExtremum[regions[*node]])
            dest[*node] = marker;
    }
    return count;
}

} // namespace lemon_graph

template <unsigned int N, class T1, class C1, 
                          class T2, class C2,
          class Compare,
          class EqualityFunctor>
unsigned int
localMinMax(MultiArrayView<N, T1, C1> const & src,
            MultiArrayView<N, T2, C2> dest,
            T1 threshold,
            Compare const & compare,
            EqualityFunctor const & equal,
            LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    vigra_precondition(src.shape() == dest.shape(),
        "localMinMax(): shape mismatch between input and output.");
        
    NeighborhoodType neighborhood;
    
    if(options.neigh == 0 || options.neigh == 2*N)
        neighborhood = DirectNeighborhood;
    else if(options.neigh == 1 || options.neigh == MetaPow<3, N>::value - 1)
        neighborhood = IndirectNeighborhood;
    else
        vigra_precondition(false,
            "localMinMax(): option object specifies invalid neighborhood type.");
    
    T2 marker = (T2)options.marker;
    
    GridGraph<N, undirected_tag> graph(src.shape(), neighborhood);
    if(options.allow_plateaus)
        return lemon_graph::extendedLocalMinMaxGraph(graph, src, dest, marker, threshold, 
                                            compare, equal, options.allow_at_border);
    else
        return lemon_graph::localMinMaxGraph(graph, src, dest, marker, threshold, 
                                             compare, options.allow_at_border);
}

/********************************************************/
/*                                                      */
/*                       localMinima                    */
/*                                                      */
/********************************************************/

// documentation is in localminmax.hxx
template <unsigned int N, class T1, class C1, class T2, class C2>
inline unsigned int
localMinima(MultiArrayView<N, T1, C1> const & src,
            MultiArrayView<N, T2, C2> dest,
            LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    T1 threshold = options.use_threshold
                           ? std::min(NumericTraits<T1>::max(), (T1)options.thresh)
                           : NumericTraits<T1>::max();
    return localMinMax(src, dest, threshold, std::less<T1>(), std::equal_to<T1>(), options);
}


template <unsigned int N, class T1, class S1,
                          class T2, class S2,
          class EqualityFunctor>
inline unsigned int
extendedLocalMinima(MultiArrayView<N, T1, S1> const & src,
                    MultiArrayView<N, T2, S2> dest,
                    EqualityFunctor const & equal,
                    LocalMinmaxOptions options = LocalMinmaxOptions())
{
    options.allowPlateaus();
    T1 threshold = options.use_threshold
                           ? std::min(NumericTraits<T1>::max(), (T1)options.thresh)
                           : NumericTraits<T1>::max();
    return localMinMax(src, dest, threshold, std::less<T1>(), equal, options);
}
/********************************************************/
/*                                                      */
/*                       localMaxima                    */
/*                                                      */
/********************************************************/

// documentation is in localminmax.hxx
template <unsigned int N, class T1, class C1, class T2, class C2>
inline unsigned int
localMaxima(MultiArrayView<N, T1, C1> const & src,
            MultiArrayView<N, T2, C2> dest,
            LocalMinmaxOptions const & options = LocalMinmaxOptions())
{
    T1 threshold = options.use_threshold
                           ? std::max(NumericTraits<T1>::min(), (T1)options.thresh)
                           : NumericTraits<T1>::min();
    return localMinMax(src, dest, threshold, std::greater<T1>(), std::equal_to<T1>(), options);
}


template <unsigned int N, class T1, class S1,
                          class T2, class S2,
          class EqualityFunctor>
inline unsigned int
extendedLocalMaxima(MultiArrayView<N, T1, S1> const & src,
                    MultiArrayView<N, T2, S2> dest,
                    EqualityFunctor const & equal,
                    LocalMinmaxOptions options = LocalMinmaxOptions())
{
    options.allowPlateaus();
    T1 threshold = options.use_threshold
                           ? std::max(NumericTraits<T1>::min(), (T1)options.thresh)
                           : NumericTraits<T1>::min();
    return localMinMax(src, dest, threshold, std::greater<T1>(), equal, options);
}

} // namespace vigra

#endif // VIGRA_MULTI_LOCALMINMAX_HXX
