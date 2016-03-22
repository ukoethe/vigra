#ifndef VIGRA_POLYTOPE_HXX
#define VIGRA_POLYTOPE_HXX

#include <set>
#include <lemon/list_graph.h>

#include "config.hxx"
#include "error.hxx"
#include "tinyvector.hxx"
#include "array_vector.hxx"
#include "linear_algebra.hxx"
#include "numerictraits.hxx"
#include "permutation.hxx"

namespace vigra {

template <unsigned int N, class T>
class Polytope
{
  public:

    enum Dimension {dimension = N};
    enum node_enum {INVALID, FACET, VERTEX};

    template <node_enum NodeType>
    struct node_type_iterator;

    typedef T                                       coordinate_type;
    typedef typename NumericTraits<T>::RealPromote  real_type;
    typedef TinyVector<T, N>                        point_type;
    typedef TinyVectorView<T, N>                    point_view_type;
    typedef typename point_type::difference_type    difference_type;
    typedef typename lemon::ListDigraph             graph_type;
    typedef typename graph_type::Node               node_type;
    typedef typename graph_type::Arc                arc_type;
    typedef typename graph_type::NodeIt             node_iterator;
    typedef typename graph_type::OutArcIt           out_arc_iterator;
    typedef typename graph_type::InArcIt            in_arc_iterator;
    typedef node_type_iterator<FACET>               facet_iterator;
    typedef node_type_iterator<VERTEX>              vertex_iterator;

    Polytope()
    : graph_()
    , type_map_(graph_)
    , vec_map_(graph_)
    , aligns_map_(graph_)
    {}

    virtual bool contains(const point_view_type & p) const = 0;

    virtual real_type nVolume() const = 0;

    virtual real_type nSurface() const = 0;

    virtual bool closed(const node_type n) const
    {
        vigra_precondition(
                type_map_[n] == FACET,
                "Polytope::closed(): Node needs do be a facet");
        return (aligns_map_[n].size() == dimension);
    }

    virtual bool closed() const
    {
        for (facet_iterator n(graph_, type_map_); n != lemon::INVALID; ++n)
        {
            if (!(this->closed(n)))
            {
                return false;
            }
        }
        return true;
    }


    virtual node_type addVertex(const point_view_type & p)
    {
        node_type ret = graph_.addNode();
        type_map_[ret] = VERTEX;
        vec_map_[ret] = p;
        aligns_map_[ret] = std::set<node_type>();
        return ret;
    }

    virtual void eraseFacet(const node_type u)
    {
        vigra_precondition(
                type_map_[u] == FACET,
                "Polytope::eraseFacet(): Node needs to be a facet");
        for (auto neighbor : aligns_map_[u])
        {
            aligns_map_[neighbor].erase(u);
        }
        graph_.erase(u);
    }

    virtual std::set<node_type> getConnected(const node_type u) const
    {
        std::set<node_type> ret;
        if (type_map_[u] == FACET)
        {
            for (out_arc_iterator a(graph_, u); a != lemon::INVALID; ++a)
            {
                ret.insert(graph_.target(a));
            }
        }
        else
        {
            for (in_arc_iterator a(graph_, u); a != lemon::INVALID; ++a)
            {
                ret.insert(graph_.source(a));
            }
        }
        return ret;
    }

    // TODO remove
    virtual ArrayVector<point_view_type> getVertices(const node_type u) const
    {
        vigra_precondition(
                type_map_[u] == FACET,
                "Polytope::getVertices(): Node must be a facet");
        ArrayVector<point_view_type> ret;
        for (out_arc_iterator a(graph_, u); a != lemon::INVALID; ++a)
        {
            ret.push_back(vec_map_[graph_.target(a)]);
        }
        return ret;
    }

    virtual ArrayVector<node_type> litFacets(const point_view_type & p) const
    {
        ArrayVector<node_type> ret;
        for (facet_iterator n(graph_, type_map_); n != lemon::INVALID; ++n)
        {
            if (distance(n, p) > 0)
            {
                ret.push_back(n);
            }
        }
        return ret;
    }

    virtual void tidyUp()
    {
        std::set<node_type> to_erase;
        for (vertex_iterator v(graph_, type_map_); v != lemon::INVALID; ++v)
        {
            vigra_assert(
                    type_map_[v] == VERTEX,
                    "Polytope::tidyUp(): vertex not a vertex");
            in_arc_iterator a(graph_, v);
            if (a == lemon::INVALID)
            {
                to_erase.insert(v);
            }
        }
        for (node_type v : to_erase)
        {
            graph_.erase(v);
        }
    }

    virtual real_type distance(const node_type u, const point_view_type & p) const
    {
        vigra_precondition(
                type_map_[u] == FACET,
                "Polytope::distance(): Node must be a facet");
        out_arc_iterator a(graph_, u);
        vigra_precondition(
                a != lemon::INVALID,
                "Polytope::distance(): Invalid facet");

        return dot(p - vec_map_[graph_.target(a)], vec_map_[u]);
    }

  protected:

    virtual node_type findNeighbor(
            const node_type u,
            const difference_type index) const
    {
        vigra_precondition(
                type_map_[u] == FACET,
                "Polytope::findNeighbor(): Node must be a facet");
        vigra_precondition(
                index < dimension,
                "Polytope::findNeighbor(): Invalid index");
        vigra_precondition(
                countOutArcs(graph_, u) == dimension,
                "Polytope::findNeighbor(): Bad facet");
        std::set<node_type> candidates;
        for (out_skip_iterator a(graph_, u, index); a != lemon::INVALID; ++a)
        {
            const node_type vertex = graph_.target(a);
            if (a.index() == 0)
            {
                candidates = getConnected(vertex);
            }
            else
            {
                std::set<node_type> to_erase;
                const auto con = getConnected(vertex);
                for (auto candidate : candidates)
                {
                    if (con.count(candidate) == 0)
                    {
                        to_erase.insert(candidate);
                    }
                }
                for (auto candidate : to_erase)
                {
                    candidates.erase(candidate);
                }
            }
        }
        candidates.erase(u);
        if (candidates.size() != 1)
        {
            vigra_assert(
                    candidates.size() == 0,
                    "Polytope::findNeighbor(): Too many aligned facets");
            return lemon::INVALID;
        }
        else
        {
            return *(candidates.begin());
        }
    }

    void assignNeighbors(const node_type u)
    {
        vigra_precondition(
                type_map_[u] == FACET,
                "Polytope::assignNeighbors(): Node must be facet");
        for (int i = 0; i < dimension; i++)
        {
            node_type neighbor = this->findNeighbor(u, i);
            if (neighbor != lemon::INVALID)
            {
                aligns_map_[u].insert(neighbor);
            }
        }
    }

    std::set<node_type> openEdge(const node_type u)
    {
        vigra_precondition(
                type_map_[u] == FACET,
                "Polytope::openEdge(): Node must be facet");
        vigra_assert(
                lemon::countOutArcs(graph_, u) == dimension,
                "Polytope::openEdge(): Got invalid facet");
        std::set<node_type> ret;
        for (int i = 0; i < dimension; i++)
        {
            if ((this->findNeighbor(u, i)) == lemon::INVALID)
            {
                for (out_skip_iterator a(graph_, u, i); a != lemon::INVALID; ++a)
                {
                    vigra_assert(
                            ret.count(graph_.target(a)) == 0,
                            "Polytope::openEdge(): Same arc twice");
                    ret.insert(graph_.target(a));
                }
                return ret;
            }
        }
        return ret;
    }

  public:

    template <node_enum NodeType>
    struct node_type_iterator : public node_type
    {
        node_type_iterator()
        {}

        node_type_iterator(lemon::Invalid i)
        : node_type(i)
        {}

        explicit node_type_iterator(
                const graph_type & graph,
                const typename graph_type::NodeMap<node_enum> & type_map)
        : graph_(graph)
        , type_map_(type_map)
        {
            graph_.first(static_cast<node_type &>(*this));
            while (type_map_[*this] != NodeType)
            {
                graph_.next(*this);
            }
        }

        node_type_iterator<NodeType> & operator++()
        {
            while (*this != lemon::INVALID)
            {
                graph_.next(*this);
                if (type_map_[*this] == NodeType)
                {
                    return *this;
                }
            }
            return *this;
        }

        const graph_type & graph_;
        const typename graph_type::NodeMap<node_enum> & type_map_;
    };

    struct out_skip_iterator : public arc_type
    {
        out_skip_iterator()
        {}

        out_skip_iterator(lemon::Invalid i)
        : arc_type(i)
        {}

        explicit out_skip_iterator(
                const graph_type & graph,
                const node_type & node,
                const difference_type skip)
        : graph_(graph)
        , skip_(skip)
        , index_(0)
        {
            graph_.firstOut(*this, node);
            if (skip_ == 0)
            {
                graph_.nextOut(*this);
            }
        }

        out_skip_iterator & operator++()
        {
            ++index_;
            graph_.nextOut(*this);
            if (index_ == skip_)
            {
                graph_.nextOut(*this);
            }
            return *this;
        }

        difference_type index() const
        {
            return index_;
        }

        const graph_type & graph_;
        const difference_type skip_;
        difference_type index_;
    };

    graph_type graph_;
    typename graph_type::NodeMap<node_enum> type_map_;
    typename graph_type::NodeMap<point_type> vec_map_;
    typename graph_type::NodeMap<std::set<node_type> > aligns_map_;
};

template <unsigned int N, class T>
class StarPolytope : public Polytope<N, T>
{
  public:

    typedef Polytope<N, T>                          base_type;
    typedef typename base_type::coordinate_type     coordinate_type;
    typedef typename base_type::real_type           real_type;
    typedef typename base_type::point_type          point_type;
    typedef typename base_type::point_view_type     point_view_type;
    typedef typename base_type::difference_type     difference_type;
    typedef typename base_type::graph_type          graph_type;
    typedef typename base_type::node_type           node_type;
    typedef typename base_type::arc_type            arc_type;
    typedef typename base_type::node_iterator       node_iterator;
    typedef typename base_type::in_arc_iterator     in_arc_iterator;
    typedef typename base_type::out_arc_iterator    out_arc_iterator;
    typedef typename base_type::out_skip_iterator   out_skip_iterator;
    typedef typename base_type::facet_iterator      facet_iterator;
    typedef typename base_type::vertex_iterator     vertex_iterator;

    using base_type::dimension;
    using base_type::graph_;
    using base_type::vec_map_;
    using base_type::type_map_;
    using base_type::aligns_map_;
    using base_type::INVALID;
    using base_type::FACET;
    using base_type::VERTEX;

  public:

    StarPolytope()
    : base_type()
    , center_(point_type())
    {}

    StarPolytope(const point_view_type & center)
    : base_type()
    , center_(center)
    {}

    StarPolytope(
            const point_view_type & a,
            const point_view_type & b,
            const point_view_type & c,
            const point_view_type & center)
    : base_type()
    , center_(center)
    {
        vigra_precondition(
                dimension == 2,
                "StarPolytope::StarPolytope(): Signature only for use in 2D");
        node_type na = this->addVertex(a);
        node_type nb = this->addVertex(b);
        node_type nc = this->addVertex(c);
        this->addFacet(nb, nc);
        this->addFacet(na, nc);
        this->addFacet(na, nb);
    }

    StarPolytope(
            const point_view_type & a,
            const point_view_type & b,
            const point_view_type & c,
            const point_view_type & d,
            const point_view_type & center)
    : base_type()
    , center_(center)
    {
        vigra_precondition(
                dimension == 3,
                "StarPolytope::StarPolytope(): Signature only for use in 3D");
        node_type na = this->addVertex(a);
        node_type nb = this->addVertex(b);
        node_type nc = this->addVertex(c);
        node_type nd = this->addVertex(d);
        this->addFacet(nb, nc, nd);
        this->addFacet(na, nc, nd);
        this->addFacet(na, nb, nd);
        this->addFacet(na, nb, nc);
    }

    virtual point_type getCenter() const 
    {
        return center_;
    }

    virtual void assignNormal(const node_type & u)
    {
        vigra_precondition(
                type_map_[u] == FACET,
                "StarPolytope::assignNormal(): Node needs to be a facet node");
        MultiArray<2, real_type> mat(dimension-1, dimension);
        out_arc_iterator a(graph_, u);
        point_view_type vertex = vec_map_[graph_.target(a)];
        ++a;
        for (int i = 0; a != lemon::INVALID; ++a, ++i)
        {
            rowVector(mat, i) = vec_map_[graph_.target(a)] - vertex;
        }
        point_view_type normal = vec_map_[u];
        for (int i = 0; i < dimension; i++)
        {
            normal[i] = 0;
        }
        for (auto permutation : permutations_)
        {
            coordinate_type val = 1;
            for (int i = 0; i < dimension - 1; i++)
            {
                val *= mat(i, permutation[i]);
            }
            val *= permutation.sign();
            normal[permutation[dimension - 1]] += val;
        }
        if (dot(normal, vertex - center_) < 0)
        {
            normal *= -1;
        }
    }

    virtual node_type addFacet(const node_type & a, const node_type & b)
    {
        vigra_precondition(
                dimension == 2,
                "StarPolytope::addFacet(): Signature only for use in 2D");
        node_type ret = graph_.addNode();
        type_map_[ret] = FACET;
        graph_.addArc(ret, a);
        graph_.addArc(ret, b);
        vigra_assert(
                lemon::countOutArcs(graph_, ret) == dimension,
                "StarPolytope::addFacet(): Invalid facet created");
        this->assignNormal(ret);
        this->assignNeighbors(ret);
        for (auto facet : aligns_map_[ret])
        {
            vigra_assert(
                    type_map_[facet] == FACET,
                    "StarPolytope::addFacet(): Node must be facet");
            this->assignNeighbors(facet);
        }
        return ret;
    }

    virtual node_type addFacet(
            const node_type & a,
            const node_type & b,
            const node_type & c)
    {
        vigra_precondition(
                dimension == 3,
                "StarPolytope::addFacet(): Signature only for use in 3D");
        node_type ret = graph_.addNode();
        type_map_[ret] = FACET;
        graph_.addArc(ret, a);
        graph_.addArc(ret, b);
        graph_.addArc(ret, c);
        vigra_assert(
                lemon::countOutArcs(graph_, ret) == dimension,
                "StarPolytope::addFacet(): Invalid facet created");
        this->assignNormal(ret);
        this->assignNeighbors(ret);
        for (auto facet : aligns_map_[ret])
        {
            vigra_assert(
                    type_map_[facet] == FACET,
                    "StarPolytope::addFacet(): Node must be facet");
            this->assignNeighbors(facet);
        }
        return ret;
    }

    virtual bool contains(const node_type & n, const point_view_type & p) const
    {
        vigra_precondition(
                type_map_[n] == FACET,
                "StarPolytope::contains(): Node needs do be a facet");
        ArrayVector<point_view_type> vertices = this->getVertices(n);
        vertices.push_back(center_);
        MultiArray<2, coordinate_type> jp_mat(dimension, dimension);
        MultiArray<2, coordinate_type> jj_mat(dimension, dimension);
        for (int j = 0; j < dimension + 1; j++)
        {
            for (int i = 0, ii = 0; i < dimension; i++, ii++)
            {
                if (i == j)
                {
                    ii++;
                }
                rowVector(jp_mat, i) = vertices[ii] - p;
                rowVector(jj_mat, i) = vertices[ii] - vertices[j];
            }
            const coordinate_type jj_det = linalg::determinant(jj_mat);
            const coordinate_type jp_det = linalg::determinant(jp_mat);
            const coordinate_type eps = std::numeric_limits<T>::epsilon() * 2;
            if (((jj_det > 0) xor (jp_det > 0)) and abs(jp_det) > eps)
            {
                return false;
            }
        }
        return true;
    }

    virtual bool contains(const point_view_type & p) const
    {
        for (facet_iterator n(graph_, type_map_); n != lemon::INVALID; ++n)
        {
            if (contains(n, p))
            {
                return true;
            }
        }
        return false;
    }

    virtual real_type nVolume(const node_type & n) const
    {
        vigra_precondition(type_map_[n] == FACET,
                           "StarPolytope::nVolume(): Node needs do be a facet");
        MultiArray<2, coordinate_type> mat(dimension, dimension);
        real_type fac = 1;
        out_arc_iterator a(graph_, n);
        for (int i = 0; i < dimension; ++i, ++a)
        {
            fac *= (i+1);
            rowVector(mat, i) = vec_map_[graph_.target(a)] - center_;
        }
        return abs(linalg::determinant(mat) / fac);
    }

    virtual real_type nVolume() const
    {
        real_type ret = 0;
        for (facet_iterator n(graph_, type_map_); n != lemon::INVALID; ++n)
        {
            ret += this->nVolume(n);
        }
        return ret;
    }

    virtual real_type nSurface(const node_type & n) const
    {
        vigra_precondition(type_map_[n] == FACET,
                           "StarPolytope::nVolume(): Node needs do be a facet");
        MultiArray<2, coordinate_type> mat(dimension, dimension);
        real_type factor = vec_map_[n].magnitude();
        out_arc_iterator a(graph_, n);
        const point_view_type vec = vec_map_[graph_.target(a)];
        ++a;
        for (int i = 1; i < dimension; ++i, ++a)
        {
            factor *= i;
            rowVector(mat, i) = vec_map_[graph_.target(a)] - vec;
        }
        rowVector(mat, 0) = vec_map_[n];
        return abs(linalg::determinant(mat)) / factor;
    }

    virtual real_type nSurface() const
    {
        real_type ret = 0;
        for (facet_iterator n(graph_, type_map_); n != lemon::INVALID; ++n)
        {
            ret += this->nSurface(n);
        }
        return ret;
    }

  protected:

    PlainChangesPermutations<N> permutations_;
    point_type center_;
};

template <unsigned int N, class T>
class ConvexPolytope : public StarPolytope<N, T>
{
  public:

    typedef StarPolytope<N, T>                      base_type;
    typedef typename base_type::coordinate_type     coordinate_type;
    typedef typename base_type::real_type           real_type;
    typedef typename base_type::point_type          point_type;
    typedef typename base_type::point_view_type     point_view_type;
    typedef typename base_type::difference_type     difference_type;
    typedef typename base_type::graph_type          graph_type;
    typedef typename base_type::node_type           node_type;
    typedef typename base_type::arc_type            arc_type;
    typedef typename base_type::node_iterator       node_iterator;
    typedef typename base_type::in_arc_iterator     in_arc_iterator;
    typedef typename base_type::out_arc_iterator    out_arc_iterator;
    typedef typename base_type::out_skip_iterator   out_skip_iterator;
    typedef typename base_type::facet_iterator      facet_iterator;
    typedef typename base_type::vertex_iterator     vertex_iterator;

    using base_type::dimension;
    using base_type::graph_;
    using base_type::vec_map_;
    using base_type::type_map_;
    using base_type::aligns_map_;
    using base_type::INVALID;
    using base_type::FACET;
    using base_type::VERTEX;

  public:

    ConvexPolytope()
    : base_type()
    {}

    ConvexPolytope(const point_view_type & center)
    : base_type(center)
    {}

    ConvexPolytope(
            const point_view_type & a,
            const point_view_type & b,
            const point_view_type & c)
    : base_type(a, b, c, (a + b + c) / 3)
    {}

    ConvexPolytope(
            const point_view_type & a,
            const point_view_type & b,
            const point_view_type & c,
            const point_view_type & d)
    : base_type(a, b, c, d, (a + b + c + d) / 4)
    {}

  protected:

    virtual void closeFacet(
            const node_type & vertex,
            const node_type & facet)
    {
        vigra_precondition(
                type_map_[vertex] == VERTEX,
                "StarPolytope::closeFacet(): Vertex needs to be a vertex node");
        vigra_precondition(
                type_map_[facet] == FACET,
                "StarPolytope::closeFacet(): Facet needs to be a facet node");
        vigra_assert(
                (this->getConnected(facet)).count(vertex) == 0,
                "StarPolytope::closeFacet(): Cannot close facet with vertex");

        while (!(this->closed(facet)))
        {
            std::set<node_type> vertices = this->openEdge(facet);
            vigra_assert(
                    vertices.size() == (dimension - 1),
                    "StarPolytope::closeFacet(): Invalid facet");
            node_type new_facet = graph_.addNode();
            type_map_[new_facet] = FACET;
            graph_.addArc(new_facet, vertex);
            for (auto n : vertices)
            {
                graph_.addArc(new_facet, n);
            }
            vigra_assert(
                    lemon::countOutArcs(graph_, new_facet) == dimension,
                    "ConvexPolytope::closeFacet(): Invalid facet created");
            this->assignNormal(new_facet);
            this->assignNeighbors(new_facet);
            for (auto neighbor : aligns_map_[new_facet])
            {
                vigra_assert(
                        type_map_[facet] == FACET,
                        "StarPolytope::addFacet(): Node must be facet");
                this->assignNeighbors(neighbor);
            }
        }
    }

  public:

    virtual bool contains(const node_type & n, const point_view_type & p) const
    {
        vigra_precondition(
                type_map_[n] == FACET,
                "ConvexPolytope::contains(): Node needs do be a facet");
        const out_arc_iterator a(graph_, n);
        const point_view_type vertex = vec_map_[graph_.target(a)];
        const point_view_type normal = vec_map_[n];
        const real_type scalar = dot(p - vertex, normal);
        return (scalar < std::numeric_limits<T>::epsilon() * 2);
    }

    virtual bool contains(const point_view_type & p) const
    {
        for (facet_iterator n(graph_, type_map_); n != lemon::INVALID; ++n)
        {
            if (!contains(n, p))
            {
                return false;
            }
        }
        return true;
    }

    virtual void addExtremeVertex(const point_view_type & p)
    {
        vigra_precondition(
                this->closed(),
                "ConvexPolytope::addExtremeVertex(): Polytope needs to be closed");
        ArrayVector<node_type> lit_facets = this->litFacets(p);
        std::set<node_type> open_facets;
        for (node_type lit_facet : lit_facets)
        {
            for (auto con : aligns_map_[lit_facet])
            {
                vigra_assert(
                        type_map_[con] == FACET,
                        "ConvexPolytope::addExtremeVertex(): facet not a facet");
                open_facets.insert(con);
            }
            open_facets.erase(lit_facet);
            this->eraseFacet(lit_facet);
        }
        this->tidyUp();
        node_type new_vertex = this->addVertex(p);
        for (auto open_facet : open_facets)
        {
            this->closeFacet(new_vertex, open_facet);
        }
    }
};

} /* namespace vigra */

#endif /* VIGRA_POLYTOPE_HXX */
