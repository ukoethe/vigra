/************************************************************************/
/*                                                                      */
/*     Copyright 2011-2015 by Stefan Schmidt, Philip Schill and         */
/*                            Ullrich Koethe                            */
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

#include <stdexcept>
#include <vector>
#include <map>

#include "metaprogramming.hxx"
#include "tinyvector.hxx"
#include "filter_iterator.hxx"

#ifdef WITH_BOOST_GRAPH

#  include <boost/tuple/tuple.hpp>
#  include <boost/graph/graph_traits.hpp>
#  include <boost/graph/properties.hpp>

namespace vigra {
namespace boost_graph {

// vigra::boost_graph contains algorithms that are compatible to the Boost Graph Library
using namespace boost;

}} // namespace vigra::boost_graph

#else // not WITH_BOOST_GRAPH

// emulate the BGL-style interface
namespace vigra {
namespace boost_graph {

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

}} // namespace vigra::boost_graph

#endif // WITH_BOOST_GRAPH

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



namespace vigra {
namespace detail {

template <typename INDEXTYPE>
class NodeDescriptor
{
public:
    typedef INDEXTYPE index_type;
    NodeDescriptor(lemon::Invalid = lemon::INVALID)
        : id_(-1)
    {}
    explicit NodeDescriptor(index_type const & id)
        : id_(id)
    {}
    bool operator!=(NodeDescriptor const & other) const
    {
        return id_ != other.id_;
    }
    bool operator==(NodeDescriptor const & other) const
    {
        return id_ == other.id_;
    }
    bool operator<(NodeDescriptor const & other) const
    {
        return id_ < other.id_;
    }
    index_type id() const
    {
        return id_;
    }
protected:
    index_type id_;
};

template <typename INDEXTYPE>
std::ostream & operator << (std::ostream & os, NodeDescriptor<INDEXTYPE> const & item)
{
    return os << item.id();
}

template <typename INDEXTYPE>
class ArcDescriptor
{
public:
    typedef INDEXTYPE index_type;
    ArcDescriptor(lemon::Invalid = lemon::INVALID)
        : id_(-1)
    {}
    explicit ArcDescriptor(index_type const & id)
        : id_(id)
    {}
    bool operator!=(ArcDescriptor const & other) const
    {
        return id_ != other.id_;
    }
    bool operator==(ArcDescriptor const & other) const
    {
        return id_ == other.id_;
    }
    bool operator<(ArcDescriptor const & other) const
    {
        return id_ < other.id_;
    }
    index_type id() const
    {
        return id_;
    }
protected:
    index_type id_;
};

template <typename INDEXTYPE>
std::ostream & operator << (std::ostream & os, ArcDescriptor<INDEXTYPE> const & item)
{
    return os << item.id();
}

} // namespace detail



enum ContainerTag
{
    MapTag,
    VectorTag,
    IndexVectorTag
};

/**
 * @brief The PropertyMap is used to store Node or Arc information of graphs.
 * 
 * @tparam <KEYTYPE> the key type
 * @tparam <MAPPEDTYPE> the mapped type
 * @tparam <ContainerTag = MapTag> whether to use a map or a vector as underlying storage
 * 
 * @note
 * In contrast to std::map, operator[] does not insert elements. Use insert() instead.
 * If ContainerTag == MapTag: at() and operator[] behave like std::map::at().
 * If ContainerTag == IndexVectorTag: at() behaves like std::map::at(). operator[] does not check if the key exists.
 */
template <typename KEYTYPE, typename MAPPEDTYPE, ContainerTag = MapTag >
class PropertyMap
{
public:
    typedef KEYTYPE key_type;
    typedef MAPPEDTYPE mapped_type;
    typedef std::pair<key_type const, mapped_type> value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    typedef std::map<key_type, mapped_type> Map;
    typedef typename Map::iterator iterator;
    typedef typename Map::const_iterator const_iterator;

    mapped_type & at(key_type const & k)
    {
        return map_.at(k);
    }
    mapped_type const & at(key_type const & k) const
    {
        return map_.at(k);
    }
    mapped_type & operator[](key_type const & k)
    {
        return map_.at(k);
    }
    mapped_type const & operator[](key_type const & k) const
    {
        return map_.at(k);
    }
    void insert(key_type const & k, mapped_type const & v)
    {
        map_[k] = v;
    }
    iterator begin()
    {
        return map_.begin();
    }
    const_iterator begin() const
    {
        return map_.begin();
    }
    const_iterator cbegin() const
    {
        return map_.cbegin();
    }
    iterator end()
    {
        return map_.end();
    }
    const_iterator end() const
    {
        return map_.end();
    }
    const_iterator cend() const
    {
        return map_.cend();
    }
    void clear()
    {
        map_.clear();
    }
    iterator find(key_type const & k)
    {
        return map_.find(k);
    }
    const_iterator find(key_type const & k) const
    {
        return map_.find(k);
    }
    size_t size() const
    {
        return map_.size();
    }
    size_t erase(key_type const & k)
    {
        return map_.erase(k);
    }

protected:
    Map map_;
};



namespace detail
{
    template <typename VALUE_TYPE>
    struct PMapValueSkipper
    {
    public:
        typedef VALUE_TYPE value_type;
        typedef typename value_type::first_type key_type;
        PMapValueSkipper(key_type default_key)
            :
            default_key_(default_key)
        {}
        bool operator()(value_type const & v)
        {
            return v.first != default_key_;
        }
    private:
        key_type const default_key_;
    };
}

/**
 * @brief Specialization of PropertyMap that stores the elements in a vector (size = max node id of stored elements).
 */
template <typename KEYTYPE, typename MAPPEDTYPE>
class PropertyMap<KEYTYPE, MAPPEDTYPE, VectorTag>
{
public:
    typedef KEYTYPE key_type;
    typedef MAPPEDTYPE mapped_type;
    typedef std::pair<key_type, mapped_type> value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    typedef std::vector<value_type> Map;
    typedef detail::PMapValueSkipper<value_type> ValueSkipper;
    typedef FilterIterator<ValueSkipper, typename Map::iterator> iterator;
    typedef FilterIterator<ValueSkipper, typename Map::const_iterator> const_iterator;

    PropertyMap(key_type default_key = lemon::INVALID)
        :
        num_elements_(0),
        default_key_(default_key)
    {}

    mapped_type & at(key_type const & k)
    {
#ifdef VIGRA_CHECK_BOUNDS
        if (k.id() < 0 || k.id() >= map_.size() || map_[k.id()].first == default_key_)
            throw std::out_of_range("PropertyMap::at(): Key not found.");
#endif
        return map_[k.id()].second;
    }

    mapped_type const & at(key_type const & k) const
    {
#ifdef VIGRA_CHECK_BOUNDS
        if (k.id() < 0 || k.id() >= map_.size() || map_[k.id()].first == default_key_)
            throw std::out_of_range("PropertyMap::at(): Key not found.");
#endif
        return map_[k.id()].second;
    }

    mapped_type & operator[](key_type const & k)
    {
        return map_[k.id()].second;
    }

    mapped_type const & operator[](key_type const & k) const
    {
        return map_[k.id()].second;
    }

    void insert(key_type const & k, mapped_type const & v)
    {
        if (k.id() < 0)
            throw std::out_of_range("PropertyMap::insert(): Key must not be negative.");

        if ((size_t)k.id() >= map_.size())
            map_.resize(k.id()+1, value_type(default_key_, mapped_type()));

        auto & elt = map_[k.id()];
        if (elt.first == default_key_)
            ++num_elements_;

        elt.first = k;
        elt.second = v;
    }

#define MAKE_ITER(it) make_filter_iterator(ValueSkipper(default_key_), it, map_.end())
#define MAKE_CITER(it) make_filter_iterator(ValueSkipper(default_key_), it, map_.cend())

    iterator begin()
    {
        return MAKE_ITER(map_.begin());
    }

    const_iterator begin() const
    {
        return MAKE_ITER(map_.begin());
    }

    const_iterator cbegin() const
    {
        return MAKE_CITER(map_.cbegin());
    }

    iterator end()
    {
        return MAKE_ITER(map_.end());
    }

    const_iterator end() const
    {
        return MAKE_ITER(map_.end());
    }

    const_iterator cend() const
    {
        return MAKE_CITER(map_.cend());
    }

    iterator find(key_type const & k)
    {
        if (k.id() < 0 || k.id() >= map_.size() || map_[k.id()].first == default_key_)
            return end();
        else
            return MAKE_ITER(std::next(map_.begin(), k.id()));
    }

    const_iterator find(key_type const & k) const
    {
        if (k.id() < 0 || k.id() >= map_.size() || map_[k.id()].first == default_key_)
            return end();
        else
            return MAKE_ITER(std::next(map_.begin(), k.id()));
    }

#undef MAKE_ITER
#undef MAKE_CITER

    void clear()
    {
        map_.clear();
        num_elements_ = 0;
    }

    size_t size() const
    {
        return num_elements_;
    }

    size_t erase(key_type const & k)
    {
        if (k.id() < 0 || k.id() >= map_.size() || map_[k.id()].first == default_key_)
        {
            return 0;
        }
        else
        {
            map_[k.id()].first = default_key_;
            --num_elements_;
            return 1;
        }
    }

protected:
    Map map_;
    size_t num_elements_;
    key_type default_key_;
};



/**
 * @brief
 * Specialization of PropertyMap that stores the elements in a vector (size = number of stored elements).
 * An additional index vector is needed for bookkeeping (size = max node id of stored elements).
 */
template <typename KEYTYPE, typename MAPPEDTYPE>
class PropertyMap<KEYTYPE, MAPPEDTYPE, IndexVectorTag>
{
public:
    typedef KEYTYPE key_type;
    typedef MAPPEDTYPE mapped_type;
    typedef std::pair<key_type, mapped_type> value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    typedef std::vector<value_type> Map;
    typedef typename Map::iterator iterator;
    typedef typename Map::const_iterator const_iterator;

    mapped_type & at(key_type const & k)
    {
#ifdef VIGRA_CHECK_BOUNDS
        if (indices_.at(k.id()) == -1)
            throw std::out_of_range("PropertyMap::at(): Key not found.");
#endif
        return map_[indices_[k.id()]].second;
    }

    mapped_type const & at(key_type const & k) const
    {
#ifdef VIGRA_CHECK_BOUNDS
        if (indices_.at(k.id()) == -1)
            throw std::out_of_range("PropertyMap::at(): Key not found.");
#endif
        return map_[indices_[k.id()]].second;
    }

    mapped_type & operator[](key_type const & k)
    {
        return map_[indices_[k.id()]].second;
    }

    mapped_type const & operator[](key_type const & k) const
    {
        return map_[indices_[k.id()]].second;
    }

    void insert(key_type const & k, mapped_type const & v)
    {
        if (k.id() < 0)
            throw std::out_of_range("PropertyMap::insert(): Key must not be negative.");
        
        if (k.id() >= indices_.size())
            indices_.resize(k.id()+1, -1);
        
        if (indices_[k.id()] == -1)
        {
            indices_[k.id()] = map_.size();
            map_.push_back(value_type(k, v));
        }
    }

    iterator begin()
    {
        return map_.begin();
    }

    const_iterator begin() const
    {
        return map_.begin();
    }

    const_iterator cbegin() const
    {
        return map_.cend();
    }

    iterator end()
    {
        return map_.end();
    }

    const_iterator end() const
    {
        return map_.end();
    }

    const_iterator cend() const
    {
        return map_.cend();
    }

    void clear()
    {
        map_.clear();
        indices_.clear();
    }

    iterator find(key_type const & k)
    {
        if (k.id() < 0 || k.id() >= indices_.size() || indices_[k.id()] == -1)
            return map_.end();
        else
            return std::next(map_.begin(), indices_[k.id()]);
    }

    const_iterator find(key_type const & k) const
    {
        if (k.id() < 0 || k.id() >= indices_.size() || indices_[k.id()] == -1)
            return map_.end();
        else
            return std::next(map_.begin(), indices_[k.id()]);
    }

    size_t size() const
    {
        return map_.size();
    }
    
    size_t erase(key_type const & k)
    {
        if (k.id() < 0 || k.id() >= indices_.size() || indices_[k.id()] == -1)
        {
            return 0;
        }
        else
        {
            // Erase the element from the index vector and the map.
            size_t ind = indices_[k.id()];
            indices_[k.id()] = -1;
            map_.erase(std::next(map_.begin(), ind));

            // Adjust the indices.
            for (size_t i = 0; i < indices_.size(); ++i)
                if (indices_[i] > ind)
                    --indices_[i];
            return 1;
        }
    }

protected:
    Map map_;
    std::vector<int> indices_;
};



} // namespace vigra



#endif // VIGRA_GRAPH_HXX
