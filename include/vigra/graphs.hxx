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

#include "metaprogramming.hxx"
#include "tinyvector.hxx"

#ifdef WITH_BOOST_GRAPH

#  include <boost/tuple/tuple.hpp>
#  include <boost/graph/graph_traits.hpp>
#  include <boost/graph/properties.hpp>

#else // not WITH_BOOST_GRAPH

// emulate the BGL-style interface
namespace boost {

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


#if 1
#ifndef BOOST_TUPLE_HPP

// tie() support for std::pair, similar to Boost's one:
// (necessary because std::pair doesn't define a suitable assignment operator)
template<class T1, class T2>
class tie_adapter 
{
  public:
    tie_adapter(T1 &x, T2 &y) 
        : x_(x), y_(y)
    {}

    template<class X, class Y>
    tie_adapter & operator=(const std::pair<X, Y> &pair)
    {
        x_ = pair.first;
        y_ = pair.second;
        return *this;
    }
    
  protected:
    T1 &x_;
    T2 &y_;
};

template<class T1, class T2>
inline
tie_adapter<T1, T2>
tie(T1& t1, T2& t2) 
{
    return tie_adapter<T1, T2>(t1, t2);
}
#endif
#endif

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

    static inline vertex_descriptor null_vertex()
    {
        return vertex_descriptor(-1);
    }
};

// property_traits class template
template <typename PropMap>
struct property_traits 
{
    typedef typename PropMap::key_type    key_type;
    typedef typename PropMap::value_type  value_type; 
    typedef typename PropMap::reference   reference;
    typedef typename PropMap::category    category;
};

namespace {

vertex_index_t vertex_index;

} // anonymous namespace

} // namespace boost

#endif // WITH_BOOST_GRAPH

namespace vigra {
namespace boost_graph { 

// vigra::boost_graph contains algorithms that are compatible to the Boost Graph Library
using namespace boost;

}} // namespace vigra::boost_graph




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

#ifdef WITH_LEMON
#  include <lemon/core.h>
#else // not WITH_LEMON

// emulate the lemon interface
namespace lemon {

struct Invalid {
  public:
    bool operator==(Invalid) const { return true;  }
    bool operator!=(Invalid) const { return false; }
    bool operator< (Invalid) const { return false; }

    template <class T, int N>
    operator vigra::TinyVector<T, N>() const
    {
       return vigra::TinyVector<T, N>(-1);
    }
};

static const Invalid INVALID = Invalid();

typedef vigra::VigraTrueType   True;
typedef vigra::VigraFalseType  False;

} // namespace lemon

#endif // WITH_LEMON

namespace lemon {

template <class T>
inline bool operator==(T const & t, Invalid)
{
    return t == T(Invalid());
}

template <class T>
inline bool operator==(Invalid, T const & t)
{
    return t == T(Invalid());
}

template <class T>
inline bool operator!=(T const & t, Invalid)
{
    return t != T(Invalid());
}

template <class T>
inline bool operator!=(Invalid, T const & t)
{
    return t != T(Invalid());
}

} // namespace lemon

namespace vigra {


template<class GRAPH,class ITEM>
struct GraphItemHelper;

template<class GRAPH>
struct GraphItemHelper<GRAPH,typename GRAPH::Edge>{
    typedef typename GRAPH::index_type index_type ;
    typedef typename GRAPH::Edge Item;
    typedef typename GRAPH::EdgeIt ItemIt;


    static index_type maxItemId(const GRAPH & g){
        return g.maxEdgeId();
    }
    static index_type itemNum(const GRAPH & g){
        return g.edgeNum();
    }
    static Item itemFromId(const GRAPH & g,const index_type id){
        return g.edgeFromId(id);
    }

};

template<class GRAPH>
struct GraphItemHelper<GRAPH,typename GRAPH::Node>{
    typedef typename GRAPH::index_type index_type ;
    typedef typename GRAPH::Node Item;
    typedef typename GRAPH::NodeIt ItemIt;


    static index_type maxItemId(const GRAPH & g){
        return g.maxNodeId();
    }
    static index_type itemNum(const GRAPH & g){
        return g.nodeNum();
    }
    static Item itemFromId(const GRAPH & g,const index_type id){
        return g.nodeFromId(id);
    }
};


template<class GRAPH>
struct GraphItemHelper<GRAPH,typename GRAPH::Arc>{
    typedef typename GRAPH::index_type index_type ;
    typedef typename GRAPH::Arc Item;
    typedef typename GRAPH::ArcIt ItemIt;


    static index_type maxItemId(const GRAPH & g){
        return g.maxArcId();
    }
    static index_type itemNum(const GRAPH & g){
        return g.arcNum();
    }
    static Item itemFromId(const GRAPH & g,const index_type id){
        return g.arcFromId(id);
    }
};





namespace lemon_graph { 

// vigra::lemon_graph contains algorithms that are compatible to the LEMON graph library
using namespace lemon;

}} // namespace vigra::lemon_graph

#endif // VIGRA_GRAPH_HXX
