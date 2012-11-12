/************************************************************************/
/*                                                                      */
/*     Copyright 2011-2012 Stefan Schmidt and Ullrich Koethe            */
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

/**
 * This header provides definitions of graph-related types
 * and optionally provides a gateway to popular graph libraries
 * (for now, BGL is supported).
 */

#ifndef VIGRA_GRAPH_HXX
#define VIGRA_GRAPH_HXX


#ifdef WITH_BOOST_GRAPH
#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

 // We define the namespace vigragraph to be the same as boost,
 // since our BGL-like graph interface needs to reside in that 
 // namespace to specialize their graph and property trait classes.

 // alias the vigragraph namespace with the boost namespace
 // namespace vigragraph = boost;
#define vigragraph boost


#else // WITH_BOOST_GRAPH

// emulate the BGL-style interface in our namespace
namespace vigragraph {

struct no_property {};

// tag classes were copied from boost:
// directed_category tags
struct directed_tag { };
struct undirected_tag { };
struct bidirectional_tag : public directed_tag { };

// traversal_category tags
struct incidence_graph_tag { };
struct adjacency_graph_tag { };
struct bidirectional_graph_tag : public incidence_graph_tag { };
struct vertex_list_graph_tag { };
struct edge_list_graph_tag { };
struct adjacency_matrix_tag { };

// edge_parallel_category tags
struct allow_parallel_edge_tag { };
struct disallow_parallel_edge_tag { };

// property maps:
struct readable_property_map_tag { };
struct writable_property_map_tag { };
struct read_write_property_map_tag
    : public readable_property_map_tag, 
      public writable_property_map_tag {};
struct lvalue_property_map_tag
    : public read_write_property_map_tag {};

struct vertex_index_t {}; 

struct edge_property_tag {};

template<class T1, class T2>
inline std::pair<T1 &, T2 &>
tie(T1& t1, T2& t2) 
{
    return std::pair<T1 &, T2 &>(t1, t2);
}

// graph_traits class template
template <typename G>
struct graph_traits {
    typedef typename G::vertex_descriptor      vertex_descriptor;
    typedef typename G::edge_descriptor        edge_descriptor;
    typedef typename G::adjacency_iterator     adjacency_iterator;
    typedef typename G::out_edge_iterator      out_edge_iterator;
    typedef typename G::in_edge_iterator       in_edge_iterator;
    typedef typename G::vertex_iterator        vertex_iterator;
    typedef typename G::edge_iterator          edge_iterator;

    typedef typename G::directed_category      directed_category;
    typedef typename G::edge_parallel_category edge_parallel_category;
    typedef typename G::traversal_category     traversal_category;

    typedef typename G::vertices_size_type     vertices_size_type;
    typedef typename G::edges_size_type        edges_size_type;
    typedef typename G::degree_size_type       degree_size_type;

    static inline vertex_descriptor null_vertex();
};

#if 0
    // CHECKME
#define vertex_index vertex_index_t()

// iterator_property_map necessary?

// property_map class template
template<class GRAPH, class T>
struct property_map 
{
};

// property_traits class template
template <typename PA>
struct property_traits {
    typedef typename PA::key_type key_type;
    typedef typename PA::value_type value_type; 
    typedef typename PA::reference reference;
    typedef typename PA::category   category;
};


// put_get_helper, copied from BGL:
template <class Reference, class LvaluePropertyMap>
struct put_get_helper { };

template <class PropertyMap, class Reference, class K>
inline Reference
get(const put_get_helper<Reference, PropertyMap>& pa, const K& k)
{
    Reference v = static_cast<const PropertyMap&>(pa)[k];
    return v;
}
template <class PropertyMap, class Reference, class K, class V>
inline void
put(const put_get_helper<Reference, PropertyMap>& pa, K k, const V& v)
{
    static_cast<const PropertyMap&>(pa)[k] = v;
}

#endif
} // namespace vigragraph

#endif // WITH_BOOST_GRAPH


#if 0
namespace vigragraph {

// custom get_helper for read-only property maps (by-value)

template <class ValueType, class ReadablePropertyMap>
struct get_helper { };

template <class PropertyMap, class ValueType, class K>
inline ValueType
get(const get_helper<ValueType, PropertyMap>& pa, const K& k)
{
    const ValueType v = static_cast<const PropertyMap&>(pa)[k];
    return v;
}


// ! A fallback template for adjacent_vertices() called with
// a vertex_iterator 
// (which may be specialized to be implemented more efficiently;
//  the reason is that the iterator may have more information than
//  the plain vertex_descriptor, e.g. it knows the neighborhood
//  already which otherwise needs to be reconstructed.)
template<class GRAPH>
inline
std::pair<typename graph_traits<GRAPH>::adjacency_iterator, 
          typename graph_traits<GRAPH>::adjacency_iterator >
adjacent_vertices_at_iterator(typename graph_traits<GRAPH>::vertex_iterator const &v,
                              GRAPH const &g) 
{    
    // the default implementation just derefences the iterator
    // to yield a vertex_descriptor and forwards the call
    std::cout << "FB" << std::endl;
    return adjacent_vertices(*v, g);
}

} // namespace vigragraph
#endif

#endif // VIGRA_GRAPH_HXX
