#ifndef VIGRA_POLYTOPE_HXX
#define VIGRA_POLYTOPE_HXX

#ifndef WITH_LEMON
    #error "Should only be included with flag \"WITH_LEMON\""
#endif

#include <set>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "config.hxx"
#include "error.hxx"
#include "tinyvector.hxx"
#include "array_vector.hxx"
#include "linear_algebra.hxx"
#include "numerictraits.hxx"
#include "permutation.hxx"

namespace vigra {

/** \brief Represent an n-dimensional polytope.

 \tparam N Dimension the polytope.
 \tparam T Type of the vector components of the polytope vertices.
*/
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

    /** Default constructor creates an empty polytope class.
    */
    Polytope()
    : graph_()
    , type_map_(graph_)
    , vec_map_(graph_)
    , aligns_map_(graph_)
    {}

    /** Copy constructor.
    */
    Polytope(const Polytope<N, T> & other)
    : graph_()
    , type_map_(graph_)
    , vec_map_(graph_)
    , aligns_map_(graph_)
    {
        *this = other;
    }

    /** Copy from another polytope.
    */
    virtual void operator=(const Polytope<N, T> & other)
    {
        lemon::digraphCopy(other.graph_, graph_);
        lemon::mapCopy(other.graph_, other.type_map_, type_map_);
        lemon::mapCopy(other.graph_, other.vec_map_, vec_map_);
        lemon::mapCopy(other.graph_, other.aligns_map_, aligns_map_);
    }

    virtual bool contains(const point_view_type & p) const = 0;

    virtual real_type nVolume() const = 0;

    virtual real_type nSurface() const = 0;

    /** Check if the facet aligns with other facets at each of its ridges.
    */
    virtual bool closed(const node_type n) const
    {
        vigra_assert(
                type_map_[n] == FACET,
                "Polytope::closed(): Node needs do be a facet");
        return std::find(
                aligns_map_[n].begin(),
                aligns_map_[n].end(),
                lemon::INVALID) == aligns_map_[n].end();
    }

    /** Check if the polytope has a closed surface
    */
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


    /** Add a vertex to the polytope.
    */
    virtual node_type addVertex(const point_view_type & p)
    {
        node_type ret = graph_.addNode();
        type_map_[ret] = VERTEX;
        vec_map_[ret] = p;
        for (int i = 0; i < N; i++)
        {
            aligns_map_[ret][i] = lemon::INVALID;
        }
        return ret;
    }

    /** Erase a facet.
    */
    virtual void eraseFacet(const node_type u)
    {
        vigra_assert(
                type_map_[u] == FACET,
                "Polytope::eraseFacet(): Node needs to be a facet");
        for (auto neighbor : aligns_map_[u])
        {
            if (neighbor != lemon::INVALID)
            {
                auto it = std::find(
                        aligns_map_[neighbor].begin(),
                        aligns_map_[neighbor].end(),
                        u);
                vigra_assert(
                        it != aligns_map_[neighbor].end(),
                        "Polytope::eraseFacet(): Inconsistent aligns map");
                *it = lemon::INVALID;
            }
        }
        graph_.erase(u);
    }

    /** Get the connected elements in the graph that represents the polytope.
        If a facet node is inserted, all of its vertices will be returned, if
        a vertex node is inserted, all facets having this vertex will be
        returned.
    */
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
        vigra_assert(
                type_map_[u] == FACET,
                "Polytope::getVertices(): Node must be a facet");
        ArrayVector<point_view_type> ret;
        for (out_arc_iterator a(graph_, u); a != lemon::INVALID; ++a)
        {
            ret.push_back(vec_map_[graph_.target(a)]);
        }
        return ret;
    }

    /** Get all facets whose normal has a positive scalar product with the
        vector to the given vertex.
    */
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

    /** Remove all vertices that are not part of the polytope mesh.
    */
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

    /** Get the distance between a facet and a vertex */
    virtual real_type distance(const node_type u, const point_view_type & p) const
    {
        vigra_assert(
                type_map_[u] == FACET,
                "Polytope::distance(): Node must be a facet");
        out_arc_iterator a(graph_, u);
        vigra_assert(
                a != lemon::INVALID,
                "Polytope::distance(): Invalid facet");

        return dot(p - vec_map_[graph_.target(a)], vec_map_[u]);
    }

    /** Label all elements in the array which are inside the polytope.
    */
    virtual unsigned int fill(
            MultiArrayView<N, unsigned int> & array,
            const unsigned int label,
            const point_view_type offset,
            const point_view_type scale) const
    {
        typedef MultiArrayView<N, unsigned int> array_type;

        unsigned int ret = 0;
        typename array_type::iterator it = array.begin();
        for (it = array.begin(); it != array.end(); it++)
        {
            const typename array_type::difference_type coord = it.template get<0>();
            point_type vec;
            for (unsigned int i = 0; i < vec.size(); i++)
            {
                vec[i] = coord[i]*scale[i] + offset[i];
            }
            if (this->contains(vec))
            {
                ret++;
                *it = label;
            }
        }
        return ret;
    }

    /** Label all elements in the array which are inside the polytope.
    */
    virtual unsigned int fill(
            MultiArrayView<N, unsigned int> & array,
            const unsigned int label,
            const point_view_type offset) const
    {
        vigra_assert(
                closed(),
                "Polytope::fill(): Polytope not closed.");
        typedef MultiArrayView<N, unsigned int> array_type;

        unsigned int ret = 0;
        typename array_type::iterator it = array.begin();
        for (it = array.begin(); it != array.end(); it++)
        {
            const typename array_type::difference_type coord = it.template get<0>();
            point_type vec;
            for (unsigned int i = 0; i < vec.size(); i++)
            {
                vec[i] = coord[i] + offset[i];
            }
            if (this->contains(vec))
            {
                ret++;
                *it = label;
            }
        }
        return ret;
    }

  protected:

    virtual bool isConnected(
            const node_type vertex,
            const node_type facet) const
    {
        vigra_assert(
                type_map_[vertex] == VERTEX,
                "Polytope::isConnected(): First node must be a vertex");
        vigra_assert(
                type_map_[facet] == FACET,
                "Polytope::isConnected(): Second node must be a facet");
        for (out_arc_iterator a(graph_, facet); a != lemon::INVALID; ++a)
        {
            if (graph_.target(a) == vertex)
            {
                return true;
            }
        }
        return false;
    }

    virtual node_type findNeighbor(
            const node_type u,
            const difference_type index) const
    {
        vigra_assert(
                type_map_[u] == FACET,
                "Polytope::findNeighbor(): Node must be a facet");
        vigra_assert(
                index < dimension,
                "Polytope::findNeighbor(): Invalid index");
        vigra_assert(
                countOutArcs(graph_, u) == dimension,
                "Polytope::findNeighbor(): Bad facet");
        out_skip_iterator a(graph_, u, index);
        const node_type first_vertex = graph_.target(a);
        for (node_type candidate : getConnected(first_vertex))
        {
            if (candidate != u)
            {
                out_skip_iterator b(a);
                do
                {
                    ++b;
                    if (b == lemon::INVALID)
                    {
                        return candidate;
                    }
                } while (isConnected(graph_.target(b), candidate));
            }
        }
        return lemon::INVALID;
    }

    void assignNeighbors(const node_type u)
    {
        vigra_assert(
                type_map_[u] == FACET,
                "Polytope::assignNeighbors(): Node must be facet");
        for (int i = 0; i < dimension; i++)
        {
            aligns_map_[u][i] = this->findNeighbor(u, i);
        }
    }

    void updateInvalidNeighbors(const node_type u)
    {
        vigra_assert(
                type_map_[u] == FACET,
                "Polytope::assignNeighbors(): Node must be facet");
        for (int i = 0; i < dimension; i++)
        {
            if (aligns_map_[u][i] == lemon::INVALID)
            {
                aligns_map_[u][i] = this->findNeighbor(u, i);
            }
        }
    }

    ArrayVector<node_type> openEdge(const node_type u)
    {
        vigra_assert(
                type_map_[u] == FACET,
                "Polytope::openEdge(): Node must be facet");
        vigra_assert(
                lemon::countOutArcs(graph_, u) == dimension,
                "Polytope::openEdge(): Got invalid facet");
        ArrayVector<node_type> ret;
        for (int i = 0; i < dimension; i++)
        {
            if (aligns_map_[u][i] == lemon::INVALID)
            {
                for (out_skip_iterator a(graph_, u, i); a != lemon::INVALID; ++a)
                {
                    ret.push_back(graph_.target(a));
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

        explicit node_type_iterator(
                const graph_type & graph,
                const typename graph_type::NodeMap<node_enum> & type_map)
        : graph_(graph)
        , type_map_(type_map)
        {
            graph_.first(static_cast<node_type &>(*this));
            while (*this != lemon::INVALID && type_map_[*this] != NodeType)
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

        bool operator==(lemon::Invalid i) const
        {
            return (static_cast<node_type>(*this) == i);
        }

        bool operator!=(lemon::Invalid i) const
        {
            return (static_cast<node_type>(*this) != i);
        }

        const graph_type & graph_;
        const typename graph_type::NodeMap<node_enum> & type_map_;
    };

    struct out_skip_iterator : public arc_type
    {
        out_skip_iterator()
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

        bool operator==(lemon::Invalid i) const
        {
            return (static_cast<arc_type>(*this) == i);
        }

        bool operator!=(lemon::Invalid i) const
        {
            return (static_cast<arc_type>(*this) != i);
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
    typename graph_type::NodeMap<TinyVector<node_type, N> > aligns_map_;
};

/** \brief Specialization of the polytope to polytopes which forms a star
    domain.
*/
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

    /** Constructor creates an empty StarPolytope with its center a the orign.
    */
    StarPolytope()
    : base_type()
    , center_(point_type())
    {}

    /** Copy constructor.
    */
    StarPolytope(const point_view_type & center)
    : base_type()
    , center_(center)
    {}

    /** Constructor for the 2-dimensional case taking three vertices and the
        center.
    */
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

    /** Constructor for the 3-dimensional case taking four vertices and the
        center.
    */
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

    /** Get the center of the star domain.
    */
    virtual point_type getCenter() const
    {
        return center_;
    }

    virtual void assignNormal(const node_type & u)
    {
        vigra_assert(
                type_map_[u] == FACET,
                "StarPolytope::assignNormal(): Node needs to be a facet node");
        MultiArray<2, real_type> mat(dimension-1, dimension);
        out_arc_iterator a(graph_, u);
        point_view_type vertex = vec_map_[graph_.target(a)];
        ++a;
        for (int i = 0; a != lemon::INVALID; ++a, ++i)
        {
            const point_type vec = vec_map_[graph_.target(a)] - vertex;
            std::copy(vec.begin(), vec.end(), rowVector(mat, i).begin());
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

    /** Add a facet to a 2-dimensional polytope.
    */
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
            if (facet != lemon::INVALID)
            {
                vigra_assert(
                        type_map_[facet] == FACET,
                        "StarPolytope::addFacet(): Node must be facet");
                this->updateInvalidNeighbors(facet);
            }
        }
        return ret;
    }

    /** Add a facet to a 3-dimensional polytope.
    */
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
            if (facet != lemon::INVALID)
            {
                vigra_assert(
                        type_map_[facet] == FACET,
                        "StarPolytope::addFacet(): Node must be facet");
                this->updateInvalidNeighbors(facet);
            }
        }
        return ret;
    }

    virtual void close()
    {
        vigra_precondition(
                lemon::countNodes(graph_) == dimension + 1,
                "StarPolytope::close(): Can only close for dim+1 vertices");
        // Set center of polytope
        {
            vertex_iterator v(graph_, type_map_);
            center_ = vec_map_[v];
            for (++v; v != lemon::INVALID; ++v)
            {
                center_ += vec_map_[v];
            }
            center_ /= static_cast<real_type>(dimension + 1);
        }
        // Create facets
        for (int i = 0; i < dimension + 1; i++)
        {
            node_type facet = graph_.addNode();
            type_map_[facet] = FACET;
            vertex_iterator v(graph_, type_map_);
            for (int j = 0; j < dimension; ++j, ++v)
            {
                if (i == j)
                {
                    ++v;
                }
                graph_.addArc(facet, v);
            }
            vigra_assert(
                    lemon::countOutArcs(graph_, facet) == dimension,
                    "StarPolytope::close(): Invalid facet created");
            this->assignNormal(facet);
            this->assignNeighbors(facet);
            for (auto neighbor : aligns_map_[facet])
            {
                if (neighbor != lemon::INVALID)
                {
                    vigra_assert(
                            type_map_[neighbor] == FACET,
                            "StarPolytope::close(): Node must be facet");
                    this->updateInvalidNeighbors(neighbor);
                }
            }
        }
    }

    virtual bool contains(const node_type & n, const point_view_type & p) const
    {
        vigra_assert(
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
                {
                const point_type vec = vertices[ii] - p;
                std::copy(vec.begin(), vec.end(), rowVector(jp_mat, i).begin());
                }
                {
                const point_type vec = vertices[ii] - vertices[j];
                std::copy(vec.begin(), vec.end(), rowVector(jj_mat, i).begin());
                }
            }
            const coordinate_type jj_det = linalg::determinant(jj_mat);
            const coordinate_type jp_det = linalg::determinant(jp_mat);
            const coordinate_type eps = std::numeric_limits<T>::epsilon() * 2;
            if (((jj_det > 0) != (jp_det > 0)) && abs(jp_det) > eps)
            {
                return false;
            }
        }
        return true;
    }

    /** Check if a point is inside the polytope.
    */
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
        vigra_assert(
                type_map_[n] == FACET,
                "StarPolytope::nVolume(): Node needs do be a facet");
        MultiArray<2, coordinate_type> mat(dimension, dimension);
        real_type fac = 1;
        out_arc_iterator a(graph_, n);
        for (int i = 0; i < dimension; ++i, ++a)
        {
            fac *= (i+1);
            const point_type vec = vec_map_[graph_.target(a)] - center_;
            std::copy(vec.begin(), vec.end(), rowVector(mat, i).begin());
        }
        return abs(linalg::determinant(mat) / fac);
    }

    /** Calculate the volume of the polytope.
    */
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
        vigra_assert(
                type_map_[n] == FACET,
                "StarPolytope::nVolume(): Node needs do be a facet");
        MultiArray<2, coordinate_type> mat(dimension, dimension);
        real_type factor = vec_map_[n].magnitude();
        out_arc_iterator a(graph_, n);
        const point_view_type vec = vec_map_[graph_.target(a)];
        ++a;
        for (int i = 1; i < dimension; ++i, ++a)
        {
            factor *= i;
            const point_type tmp = vec_map_[graph_.target(a)] - vec;
            std::copy(tmp.begin(), tmp.end(), rowVector(mat, i).begin());
        }
        const point_type tmp = vec_map_[n];
        std::copy(tmp.begin(), tmp.end(), rowVector(mat, 0).begin());
        return abs(linalg::determinant(mat)) / factor;
    }

    /** Calculate the surface of the polytope.
    */
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

/** Specialization of the StarPolytope to polytopes which have a convex domain.
*/
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
        vigra_assert(
                type_map_[vertex] == VERTEX,
                "ConvexPolytope::closeFacet(): Vertex needs to be a vertex node");
        vigra_assert(
                type_map_[facet] == FACET,
                "ConvexPolytope::closeFacet(): Facet needs to be a facet node");
        vigra_assert(
                (this->getConnected(facet)).count(vertex) == 0,
                "ConvexPolytope::closeFacet(): Cannot close facet with vertex");

        while (!(this->closed(facet)))
        {
            ArrayVector<node_type> vertices = this->openEdge(facet);
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
                if (neighbor != lemon::INVALID)
                {
                    vigra_assert(
                            type_map_[facet] == FACET,
                            "StarPolytope::addFacet(): Node must be facet");
                    this->updateInvalidNeighbors(neighbor);
                }
            }
        }
    }

  public:

    virtual bool contains(const node_type & n, const point_view_type & p) const
    {
        vigra_assert(
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

    /** Expand the polytope to the given point if it's outside of the current
        polytope, such that the new polytope is still convex.
    */
    virtual void addExtremeVertex(const point_view_type & p)
    {
        vigra_assert(
                this->closed(),
                "ConvexPolytope::addExtremeVertex(): Polytope needs to be closed");
        ArrayVector<node_type> lit_facets = this->litFacets(p);
        if (lit_facets.size() > 0)
        {
            std::set<node_type> open_facets;
            for (node_type lit_facet : lit_facets)
            {
                for (auto con : aligns_map_[lit_facet])
                {
                    if (con != lemon::INVALID)
                    {
                        vigra_assert(
                                type_map_[con] == FACET,
                                "ConvexPolytope::addExtremeVertex(): "
                                "facet not a facet");
                        open_facets.insert(con);
                    }
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
    }
};

} /* namespace vigra */

#endif /* VIGRA_POLYTOPE_HXX */
