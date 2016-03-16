/************************************************************************/
/*                                                                      */
/*     Copyright 2011-2012 by Stefan Schmidt and Ullrich Koethe         */
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

#ifndef VIGRA_MULTI_GRIDGRAPH_HXX
#define VIGRA_MULTI_GRIDGRAPH_HXX

#include "multi_fwd.hxx"
#include "multi_iterator.hxx"
#include "multi_array.hxx"
#include "graphs.hxx"

template <unsigned int N>
struct NeighborhoodTests;

namespace vigra {

template<unsigned int N, class T, class Stride>
inline
typename vigra::MultiArrayView<N, T, Stride>::const_reference 
get(vigra::MultiArrayView<N, T, Stride> const & pmap,
    typename vigra::MultiArrayView<N, T, Stride>::difference_type const & k)
{ 
    return pmap[k]; 
}

/** \addtogroup GraphDataStructures Graph Data Structures and Algorithms
        
        Graph algorithms and the underlying graph data structures (e.g. GridGraph and AdjacencyListGraph)
        implementing the APIs of the 
        <a href="http://www.boost.org/doc/libs/release/libs/graph/">boost::graph</a> and
        <a href="http://lemon.cs.elte.hu/">LEMON</a> libraries. 
        
        See also the \ref BoostGraphExtensions.
*/
//@{

/*
undirected edge_descriptor derived from TinyVector, 
including flag for original edge orientation, 
as necessary for source(), target() functions;
This edge_descriptor allows to directly (without adapter object)
index a MultiArrayView with one dimension more than the gridgraph, 
the last coordinate indexing the edge number
(missing edges at the border are just ignored) 

The gridgraph class is able to convert/construct these edge_descriptors
and to reconstruct the corresponding source/target nodes.

FIXME: store edge index in (*this)[0] ??
*/
template<unsigned int N>
class GridGraphArcDescriptor
    : public MultiArrayShape<N+1>::type
{
  public:
    typedef typename MultiArrayShape<N+1>::type  base_type;
    typedef typename base_type::value_type       value_type;
    typedef base_type                            edge_coord_type;
    typedef value_type                           index_type;
    typedef typename MultiArrayShape<N>::type    shape_type;
    typedef TinyVectorView<value_type, N>        vertex_descriptor_view;

    GridGraphArcDescriptor()
    : is_reversed_(false)
    {}

    GridGraphArcDescriptor(lemon::Invalid)
    : base_type(-1),
      is_reversed_(false)
    {}

    GridGraphArcDescriptor(base_type const & b, bool reversed)
    : base_type(b),
      is_reversed_(reversed)
    {}

    GridGraphArcDescriptor(shape_type const &vertex,
                           index_type edge_index,
                           bool reversed=false)
    : base_type(detail::DontInit())
    {
        set(vertex, edge_index, reversed);
    }
                                      
    void set(shape_type const &vertex, index_type edge_index, bool reversed) 
    {
        this->template subarray<0,N>() = vertex;
        (*this)[N] = edge_index;
        is_reversed_ = reversed;
    }
        
    void increment(GridGraphArcDescriptor const & diff, bool opposite=false) 
    {
        if(diff.is_reversed_)
        {
            is_reversed_ = !opposite;
            this->template subarray<0,N>() += diff.template subarray<0,N>();
        }
        else
        {
            is_reversed_ = opposite;
        }
        (*this)[N] = diff[N];
    }
        
    bool isReversed() const 
    {
        return is_reversed_;
    }
    
    vertex_descriptor_view vertexDescriptor() const
    {
        return this->template subarray<0,N>();
    }
    
    value_type edgeIndex() const
    {
        return (*this)[N];
    }

  protected:
    bool is_reversed_;
};

inline MultiArrayIndex 
gridGraphMaxDegree(unsigned int N, NeighborhoodType t)
{
    return t == DirectNeighborhood
                ? 2*N
                : pow(3.0, (int)N) - 1;
}

template <unsigned int N, NeighborhoodType>
struct GridGraphMaxDegree;

template <unsigned int N>
struct GridGraphMaxDegree<N, DirectNeighborhood>
{
    static const MultiArrayIndex value = 2*N;
};

template <unsigned int N>
struct GridGraphMaxDegree<N, IndirectNeighborhood>
{
    static const MultiArrayIndex value = MetaPow<3, N>::value - 1;
};

template <class Shape>
MultiArrayIndex 
gridGraphEdgeCount(Shape const & shape, NeighborhoodType t, bool directed)
{
    int res = 0;
    if(t == DirectNeighborhood)
    {
        for(unsigned int k=0; k<shape.size(); ++k)
            res += 2*prod(shape - Shape::unitVector(k));
    }
    else
    {
        res = prod(3*shape - Shape(2)) - prod(shape);
    }
    return directed
               ? res
               : res / 2;
}

namespace detail {

template <class Shape>
void
computeNeighborOffsets(ArrayVector<Shape> const & neighborOffsets, 
                       ArrayVector<ArrayVector<bool> > const & neighborExists,
                       ArrayVector<ArrayVector<Shape> > & incrementOffsets,
                       ArrayVector<ArrayVector<GridGraphArcDescriptor<Shape::static_size> > > & edgeDescriptorOffsets,
                       ArrayVector<ArrayVector<MultiArrayIndex> > & indices,
                       ArrayVector<ArrayVector<MultiArrayIndex> > & backIndices,
                       bool directed)
{
    typedef GridGraphArcDescriptor<Shape::static_size> EdgeDescriptor;
    
    unsigned int borderTypeCount = neighborExists.size();
    incrementOffsets.resize(borderTypeCount);
    edgeDescriptorOffsets.resize(borderTypeCount);
    indices.resize(borderTypeCount);
    backIndices.resize(borderTypeCount);
    
    for(unsigned int k=0; k<borderTypeCount; ++k)
    {
        incrementOffsets[k].clear();
        edgeDescriptorOffsets[k].clear();
        indices[k].clear();
        backIndices[k].clear();
        
        for(unsigned int j=0; j < neighborOffsets.size(); ++j)
        {
            if(neighborExists[k][j])
            {
                if(incrementOffsets[k].size() == 0)
                {
                    incrementOffsets[k].push_back(neighborOffsets[j]);
                }
                else
                {
                    incrementOffsets[k].push_back(neighborOffsets[j] - neighborOffsets[indices[k].back()]);
                }
                
                if(directed || j < neighborOffsets.size() / 2) // directed or backward edge
                {
                    edgeDescriptorOffsets[k].push_back(EdgeDescriptor(Shape(), j));
                }
                else if(edgeDescriptorOffsets[k].size() == 0 || !edgeDescriptorOffsets[k].back().isReversed()) // the first forward edge
                {
                    edgeDescriptorOffsets[k].push_back(EdgeDescriptor(neighborOffsets[j], neighborOffsets.size()-j-1, true));
                }       
                else // second or higher forward edge
                {
                    edgeDescriptorOffsets[k].push_back(EdgeDescriptor(neighborOffsets[j] - neighborOffsets[indices[k].back()], 
                                                                      neighborOffsets.size()-j-1, true));
                }
                
                indices[k].push_back(j);
                if(j < neighborOffsets.size() / 2)
                    backIndices[k].push_back(j);
            }
        }
    }
}

} // namespace detail

template<unsigned int N, bool BackEdgesOnly>
class GridGraphNeighborIterator
{
public:
    typedef typename MultiArrayShape<N>::type          shape_type;
    typedef MultiCoordinateIterator<N>                 vertex_iterator;
    typedef typename vertex_iterator::value_type       vertex_descriptor;
    typedef vertex_descriptor                          value_type;
    typedef typename vertex_iterator::pointer          pointer;
    typedef typename vertex_iterator::const_pointer    const_pointer;
    typedef typename vertex_iterator::reference        reference;
    typedef typename vertex_iterator::const_reference  const_reference;
    typedef MultiArrayIndex                            difference_type;
    typedef MultiArrayIndex                            index_type;
    typedef std::forward_iterator_tag                  iterator_category;
    
    friend struct NeighborhoodTests<N>;

    GridGraphNeighborIterator() 
    : neighborOffsets_(0),
      neighborIndices_(0),
      index_(0)
    {}

    template <class DirectedTag>
    GridGraphNeighborIterator(GridGraph<N, DirectedTag> const & g, typename GridGraph<N, DirectedTag>::Node const & v)
    : neighborOffsets_(0),
      neighborIndices_(0),
      target_(v),
      index_(0)
    {
        unsigned int nbtype = g.get_border_type(v);
        neighborOffsets_ = &(*g.neighborIncrementArray())[nbtype];
        neighborIndices_ = &(*g.neighborIndexArray(BackEdgesOnly))[nbtype];
        updateTarget();
    }

    template <class DirectedTag>
    GridGraphNeighborIterator(GridGraph<N, DirectedTag> const & g, typename GridGraph<N, DirectedTag>::NodeIt const & v)
    : neighborOffsets_(0),
      neighborIndices_(0),
      target_(v),
      index_(0)
    {
        unsigned int nbtype = g.get_border_type(v);
        neighborOffsets_ = &(*g.neighborIncrementArray())[nbtype];
        neighborIndices_ = &(*g.neighborIndexArray(BackEdgesOnly))[nbtype];
        updateTarget();
    }

    // TODO: implement a "goto-neighbor" operation
    // yielding a vertex_iterator! -> useful for 
    // watershed algo.

    GridGraphNeighborIterator & operator++()
    {
        ++index_;
        updateTarget();
        return *this;
    }

    GridGraphNeighborIterator operator++(int)
    {
        GridGraphNeighborIterator ret(*this);
        ++*this;
        return ret;
    }

    const_reference operator*() const
    {
        return target_;
    }

    const_pointer operator->() const
    {
        return &target_;
    }
    
    operator const_reference() const
    {
        return target_;
    }

    const_reference target() const
    {
        return target_;
    }

    MultiArrayIndex index() const
    {
        return index_;
    }

    MultiArrayIndex neighborIndex() const
    {
        return (*neighborIndices_)[index_];
    }

    bool operator==(GridGraphNeighborIterator const & other) const
    {
        return index_ == other.index_;
    }

    bool operator!=(GridGraphNeighborIterator const & other) const
    {
        return index_ != other.index_;
    }

    bool isValid() const
    {
        return index_ < (MultiArrayIndex)neighborIndices_->size();
    }
    
    bool atEnd() const
    {
        return index_ >= (MultiArrayIndex)neighborIndices_->size();
    }
    
    GridGraphNeighborIterator getEndIterator() const
    {
        GridGraphNeighborIterator res(*this);
        res.index_ = (MultiArrayIndex)neighborIndices_->size();
        return res;
    }

  protected:

        // for testing only
    GridGraphNeighborIterator(ArrayVector<shape_type> const & neighborOffsets,
                              ArrayVector<index_type> const & neighborIndices,
                              ArrayVector<index_type> const & backIndices,
                              vertex_descriptor source)
    : neighborOffsets_(&neighborOffsets),
      neighborIndices_(BackEdgesOnly ? &backIndices : &neighborIndices),
      target_(source),
      index_(0)
    {
        updateTarget();
    }
  
    void updateTarget()
    {
        if(isValid())
            target_ += (*neighborOffsets_)[index_];
    }

    ArrayVector<shape_type> const * neighborOffsets_;
    ArrayVector<index_type> const * neighborIndices_;
    vertex_descriptor target_;
    MultiArrayIndex index_;
};

template<unsigned int N, bool BackEdgesOnly>
class GridGraphOutEdgeIterator
{
  public:
    typedef typename MultiArrayShape<N>::type    shape_type;
    typedef MultiArrayIndex                      index_type;
    typedef GridGraphArcDescriptor<N>            arc_descriptor;
    typedef typename MultiArrayShape<N+1>::type  value_type;
    typedef value_type const &                   reference;
    typedef value_type const &                   const_reference;
    typedef value_type const *                   pointer;
    typedef value_type const *                   const_pointer;
    typedef MultiArrayIndex                      difference_type;
    typedef std::forward_iterator_tag            iterator_category;

    friend struct NeighborhoodTests<N>;
    friend class GridGraphEdgeIterator<N, BackEdgesOnly>;

    GridGraphOutEdgeIterator() 
    : neighborOffsets_(0),
      neighborIndices_(0),
      index_(0)
    {}

    template <class DirectedTag>
    GridGraphOutEdgeIterator(GridGraph<N, DirectedTag> const & g, 
                             typename GridGraph<N, DirectedTag>::NodeIt const & v,
                             bool opposite=false)
    : neighborOffsets_(0),
      neighborIndices_(0),
      edge_descriptor_(),
      index_(0)
    {
        if(v.isValid()){
            unsigned int nbtype = g.get_border_type(v);
            init(&(*g.edgeIncrementArray())[nbtype], &(*g.neighborIndexArray(BackEdgesOnly))[nbtype], *v, opposite);
        }
        else{
            index_ = (index_type)neighborIndices_->size();
        }
    }

    template <class DirectedTag>
    GridGraphOutEdgeIterator(GridGraph<N, DirectedTag> const & g, 
                             typename GridGraph<N, DirectedTag>::Node const & v,
                             bool opposite=false)
    : neighborOffsets_(0),
      neighborIndices_(0),
      edge_descriptor_(),
      index_(0)
    {
        if(isInside(g, v)){
            unsigned int nbtype = g.get_border_type(v);
            init(&(*g.edgeIncrementArray())[nbtype], &(*g.neighborIndexArray(BackEdgesOnly))[nbtype], v, opposite);
        }
        else{
            index_ = (index_type)neighborIndices_->size();
        }
    }
    
    GridGraphOutEdgeIterator & operator++()
    {
        increment(false);
        return *this;
    }

    GridGraphOutEdgeIterator  operator++(int)
    {
        GridGraphOutEdgeIterator ret(*this);
        ++*this;
        return ret;
    }

    const_reference operator*() const
    {
        return edge_descriptor_;
    }

    operator const_reference() const
    {
        return edge_descriptor_;
    }

    const_pointer operator->() const
    {
        return &edge_descriptor_;
    }

    index_type index() const
    {
        return index_;
    }

    index_type neighborIndex() const
    {
        return (*neighborIndices_)[index_];
    }

    arc_descriptor const & arcDescriptor() const
    {
        return edge_descriptor_;
    }

    bool operator==(GridGraphOutEdgeIterator const & other) const
    {
        return index_ == other.index();
    }

    bool operator!=(GridGraphOutEdgeIterator const & other) const
    {
        return index_ != other.index();
    }

    bool isValid() const 
    {
        return index_ < (index_type)neighborIndices_->size();
    }

    bool atEnd() const 
    {
        return index_ >= (index_type)neighborIndices_->size();
    }

    GridGraphOutEdgeIterator getEndIterator() const
    {
        GridGraphOutEdgeIterator res(*this);
        res.index_ = (index_type)neighborIndices_->size();
        return res;
    }

  protected:
  
        // for testing only
    GridGraphOutEdgeIterator(ArrayVector<arc_descriptor> const & neighborOffsets,
                             ArrayVector<index_type> const & neighborIndices,
                             ArrayVector<index_type> const & backIndices,
                             shape_type const & source)
    : neighborOffsets_(0),
      neighborIndices_(0),
      edge_descriptor_(),
      index_(0)
    {
        init(&neighborOffsets, BackEdgesOnly ? &backIndices : &neighborIndices, source);
    }

    void init(ArrayVector<arc_descriptor> const * neighborOffsets,
              ArrayVector<index_type> const * neighborIndices,
              shape_type const & source,
              bool opposite=false)
    {
        neighborOffsets_ = neighborOffsets;
        neighborIndices_ = neighborIndices;
        edge_descriptor_ = arc_descriptor(source, 0);
        index_ = 0;
        updateEdgeDescriptor(opposite);
    }

    void increment(bool opposite)
    {
        ++index_;
        updateEdgeDescriptor(opposite);
    }

    void updateEdgeDescriptor(bool opposite)
    {
        if(isValid())
            edge_descriptor_.increment((*neighborOffsets_)[index_], opposite);
    }
  
    ArrayVector<arc_descriptor> const * neighborOffsets_;
    ArrayVector<index_type> const * neighborIndices_;
    arc_descriptor edge_descriptor_;
    index_type index_;
};

template<unsigned int N, bool BackEdgesOnly>
class GridGraphOutArcIterator
: public GridGraphOutEdgeIterator<N, BackEdgesOnly>
{
  public:
    typedef GridGraphOutEdgeIterator<N, BackEdgesOnly>        base_type;
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef MultiArrayIndex                    index_type;
    typedef GridGraphArcDescriptor<N>          value_type;
    typedef value_type const &                 reference;
    typedef value_type const &                 const_reference;
    typedef value_type const *                 pointer;
    typedef value_type const *                 const_pointer;
    typedef MultiArrayIndex                    difference_type;
    typedef std::forward_iterator_tag          iterator_category;

    friend struct NeighborhoodTests<N>;
    friend class GridGraphEdgeIterator<N, BackEdgesOnly>;

    GridGraphOutArcIterator() 
    : base_type()
    {}

    explicit GridGraphOutArcIterator(base_type const & b) 
    : base_type(b)
    {}

    template <class DirectedTag>
    GridGraphOutArcIterator(GridGraph<N, DirectedTag> const & g, typename GridGraph<N, DirectedTag>::NodeIt const & v)
    : base_type(g, v)
    {}

    template <class DirectedTag>
    GridGraphOutArcIterator(GridGraph<N, DirectedTag> const & g, typename GridGraph<N, DirectedTag>::Node const & v)
    : base_type(g, v)
    {}

    GridGraphOutArcIterator & operator++()
    {
        base_type::operator++();
        return *this;
    }

    GridGraphOutArcIterator  operator++(int)
    {
        GridGraphOutArcIterator ret(*this);
        ++*this;
        return ret;
    }

    const_reference operator*() const
    {
        return this->edge_descriptor_;
    }

    operator const_reference() const
    {
        return this->edge_descriptor_;
    }

    const_pointer operator->() const
    {
        return &this->edge_descriptor_;
    }

    GridGraphOutArcIterator getEndIterator() const
    {
        return GridGraphOutArcIterator(base_type::getEndIterator());
    }
    
  protected:

        // for testing only
    GridGraphOutArcIterator(ArrayVector<value_type> const & neighborOffsets,
                            ArrayVector<index_type> const & neighborIndices,
                            ArrayVector<index_type> const & backIndices,
                            shape_type const & source)
    : base_type(neighborOffsets, neighborIndices, backIndices, source)
    {}
};

template<unsigned int N, bool BackEdgesOnly>
class GridGraphInArcIterator
: public GridGraphOutEdgeIterator<N, BackEdgesOnly>
{
  public:
    typedef GridGraphOutEdgeIterator<N, BackEdgesOnly>        base_type;
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef MultiArrayIndex                    index_type;
    typedef GridGraphArcDescriptor<N>          value_type;
    typedef value_type const &                 reference;
    typedef value_type const &                 const_reference;
    typedef value_type const *                 pointer;
    typedef value_type const *                 const_pointer;
    typedef MultiArrayIndex                    difference_type;
    typedef std::forward_iterator_tag          iterator_category;

    friend struct NeighborhoodTests<N>;

    GridGraphInArcIterator() 
    : base_type()
    {}

    explicit GridGraphInArcIterator(base_type const & b) 
    : base_type(b)
    {}

    template <class DirectedTag>
    GridGraphInArcIterator(GridGraph<N, DirectedTag> const & g, typename GridGraph<N, DirectedTag>::NodeIt const & v)
    : base_type(g, v, true)
    {}

    template <class DirectedTag>
    GridGraphInArcIterator(GridGraph<N, DirectedTag> const & g, typename GridGraph<N, DirectedTag>::Node const & v)
    : base_type(g, v, true)
    {}

    GridGraphInArcIterator & operator++()
    {
        base_type::increment(true);
        return *this;
    }

    GridGraphInArcIterator  operator++(int)
    {
        GridGraphInArcIterator ret(*this);
        ++*this;
        return ret;
    }

    const_reference operator*() const
    {
        return this->edge_descriptor_;
    }

    operator const_reference() const
    {
        return this->edge_descriptor_;
    }

    const_pointer operator->() const
    {
        return &this->edge_descriptor_;
    }

    GridGraphInArcIterator getEndIterator() const
    {
        return GridGraphInArcIterator(base_type::getEndIterator());
    }
};

    // Edge iterator for directed and undirected graphs. 
    // Composed of a vertex_iterator and an out_edge_iterator.
template<unsigned int N, bool BackEdgesOnly>
class GridGraphEdgeIterator
{
public:
    typedef GridGraphEdgeIterator<N, BackEdgesOnly>      self_type;
    typedef MultiCoordinateIterator<N>                   vertex_iterator;
    typedef typename vertex_iterator::value_type         vertex_descriptor;
    typedef GridGraphOutArcIterator<N, BackEdgesOnly>    out_edge_iterator;
    typedef typename MultiArrayShape<N+1>::type          edge_descriptor;
    typedef edge_descriptor                              value_type;
    typedef value_type const *                           pointer;
    typedef value_type const *                           const_pointer;
    typedef value_type const &                           reference;
    typedef value_type const &                           const_reference;
    typedef typename MultiArrayShape<N>::type            shape_type;
    typedef MultiArrayIndex                              difference_type;
    typedef MultiArrayIndex                              index_type;
    typedef std::forward_iterator_tag                    iterator_category;

    friend struct NeighborhoodTests<N>;

    GridGraphEdgeIterator() 
    : neighborOffsets_(0),
      neighborIndices_(0)
    {}
    
    template <class DirectedTag>
    GridGraphEdgeIterator(GridGraph<N, DirectedTag> const & g)
    : neighborOffsets_(g.edgeIncrementArray()),
      neighborIndices_(g.neighborIndexArray(BackEdgesOnly)),
      vertexIterator_(g),
      outEdgeIterator_(g, vertexIterator_)
    {
        if(outEdgeIterator_.atEnd()) // in a undirected graph, the first point stores no edges
        {
            ++vertexIterator_;
            if(vertexIterator_.isValid())
                outEdgeIterator_ = out_edge_iterator(g, vertexIterator_);
        }
    }



    template <class DirectedTag>
    GridGraphEdgeIterator(GridGraph<N, DirectedTag> const & g, const typename  GridGraph<N, DirectedTag>::Edge & edge)
    : neighborOffsets_(g.edgeIncrementArray()),
      neighborIndices_(g.neighborIndexArray(BackEdgesOnly)),
      vertexIterator_(g,g.u(edge)),
      outEdgeIterator_(g, vertexIterator_)
    {
        if(vertexIterator_.isValid()){
            // vigra_precondition(edge!=lemon::INVALID,"no invalid edges here");
            // vigra_precondition( allLess(*vertexIterator_,g.shape()), "fixme1");
            // vigra_precondition( allGreaterEqual(*vertexIterator_,shape_type() ), "fixme2");


            if(edge[N] >= 0  && edge[N] < g.maxUniqueDegree( ) && 
                (*( g.neighborExistsArray()))[vertexIterator_.borderType()][edge[N]] ){
                while(*outEdgeIterator_!=edge){
                    ++outEdgeIterator_;
                }
            }
            else{
                vertexIterator_ = vertexIterator_.getEndIterator();
            }

            // vigra_precondition(edge[N] >= 0 && edge[N] < g.maxUniqueDegree(),"fixme3");
            // vigra_precondition(    ,"fixme4");
            // vigra_precondition(!outEdgeIterator_.atEnd(),"fixme5");


        }
    }




    GridGraphEdgeIterator & operator++()
    {
        ++outEdgeIterator_;
        if(outEdgeIterator_.atEnd())
        {
            ++vertexIterator_;
            if(vertexIterator_.isValid())
            {
                unsigned int borderType = vertexIterator_.borderType();
                outEdgeIterator_.init(&(*neighborOffsets_)[borderType], &(*neighborIndices_)[borderType], *vertexIterator_);
            }
        }
        return *this;
    }

    GridGraphEdgeIterator  operator++(int)
    {
        GridGraphEdgeIterator ret(*this);
        ++*this;
        return ret;
    }

    const_reference operator*() const
    {
        return *outEdgeIterator_;
    }

    operator const_reference() const
    {
        return *outEdgeIterator_;
    }

    const_pointer operator->() const
    {
        return outEdgeIterator_.operator->();
    }
    
    MultiArrayIndex neighborIndex() const
    {
        return outEdgeIterator_.neighborIndex();
    }


    bool operator==(GridGraphEdgeIterator const & other) const
    {
        return (vertexIterator_ == other.vertexIterator_ && outEdgeIterator_ == other.outEdgeIterator_);
    }

    bool operator!=(GridGraphEdgeIterator const & other) const
    {
        return !operator==(other);
    }
    
    bool isValid() const
    {
        return vertexIterator_.isValid();
    }
    
    bool atEnd() const
    {
        return !isValid();
    }
    
    GridGraphEdgeIterator getEndIterator() const
    {
        GridGraphEdgeIterator ret(*this);
        ret.vertexIterator_ = vertexIterator_.getEndIterator();
        vertex_iterator lastVertex = ret.vertexIterator_ - 1;
        unsigned int borderType = lastVertex.borderType();
        ret.outEdgeIterator_.init(&(*neighborOffsets_)[borderType], &(*neighborIndices_)[borderType], *lastVertex);
        ret.outEdgeIterator_ = ret.outEdgeIterator_.getEndIterator();
        return ret;
    }

   protected:
   
        // for testing only
    GridGraphEdgeIterator(ArrayVector<ArrayVector<typename out_edge_iterator::value_type> > const & neighborOffsets,
                          ArrayVector<ArrayVector<index_type> > const & neighborIndices,
                          ArrayVector<ArrayVector<index_type> > const & backIndices,
                          shape_type const & shape)
    : neighborOffsets_(&neighborOffsets),
      neighborIndices_(BackEdgesOnly ? &backIndices : &neighborIndices),
      vertexIterator_(shape),
      outEdgeIterator_(neighborOffsets[vertexIterator_.borderType()], 
                       neighborIndices[vertexIterator_.borderType()], 
                       backIndices[vertexIterator_.borderType()], shape_type())
    {
        if(outEdgeIterator_.atEnd()) // in a undirected graph, the first point stores no edges
        {
            ++vertexIterator_;
            if(vertexIterator_.isValid())
            {
                unsigned int borderType = vertexIterator_.borderType();
                outEdgeIterator_.init(&(*neighborOffsets_)[borderType], &(*neighborIndices_)[borderType], *vertexIterator_);
            }
        }
    }

    ArrayVector<ArrayVector<typename out_edge_iterator::value_type> > const * neighborOffsets_;
    ArrayVector<ArrayVector<index_type> > const * neighborIndices_;
    vertex_iterator vertexIterator_;
    out_edge_iterator outEdgeIterator_;
};

template<unsigned int N, bool BackEdgesOnly>
class GridGraphArcIterator
: public GridGraphEdgeIterator<N, BackEdgesOnly>
{
public:
    typedef GridGraphEdgeIterator<N, BackEdgesOnly>              base_type;
    typedef GridGraphArcIterator<N, BackEdgesOnly>               self_type;
    typedef MultiCoordinateIterator<N>                           vertex_iterator;
    typedef typename vertex_iterator::value_type                 vertex_descriptor;
    typedef GridGraphOutArcIterator<N, BackEdgesOnly>            out_edge_iterator;
    typedef typename out_edge_iterator::value_type               edge_descriptor;
    typedef edge_descriptor                                      value_type;
    typedef value_type const *                                   pointer;
    typedef value_type const *                                   const_pointer;
    typedef value_type const &                                   reference;
    typedef value_type const &                                   const_reference;
    typedef typename MultiArrayShape<N>::type                    shape_type;
    typedef MultiArrayIndex                                      difference_type;
    typedef MultiArrayIndex                                      index_type;
    typedef std::forward_iterator_tag                            iterator_category;

    friend struct NeighborhoodTests<N>;

    GridGraphArcIterator() 
    : base_type()
    {}

    explicit GridGraphArcIterator(base_type const & b) 
    : base_type(b)
    {}

    template <class DirectedTag>
    GridGraphArcIterator(GridGraph<N, DirectedTag> const & g)
    : base_type(g)
    {}

    GridGraphArcIterator & operator++()
    {
        base_type::operator++();
        return *this;
    }

    GridGraphArcIterator  operator++(int)
    {
        GridGraphArcIterator ret(*this);
        ++*this;
        return ret;
    }

    const_reference operator*() const
    {
        return *(this->outEdgeIterator_);
    }

    operator const_reference() const
    {
        return *(this->outEdgeIterator_);
    }

    const_pointer operator->() const
    {
        return this->outEdgeIterator_.operator->();
    }
    
    GridGraphArcIterator getEndIterator() const
    {
        return GridGraphArcIterator(base_type::getEndIterator());
    }

  protected:
  
        // for testing only
    GridGraphArcIterator(ArrayVector<ArrayVector<value_type> > const & neighborOffsets,
                          ArrayVector<ArrayVector<index_type> > const & neighborIndices,
                          ArrayVector<ArrayVector<index_type> > const & backIndices,
                          shape_type const & shape)
    : base_type(neighborOffsets, neighborIndices, backIndices, shape)
    {}
};

template<unsigned int N>
inline bool operator==(MultiCoordinateIterator<N> const & i, lemon::Invalid)
{
    return i.atEnd();
}

template<unsigned int N>
inline bool operator!=(MultiCoordinateIterator<N> const & i, lemon::Invalid)
{
    return i.isValid();
}

template<unsigned int N>
inline bool operator==(lemon::Invalid, MultiCoordinateIterator<N> const & i)
{
    return i.atEnd();
}

template<unsigned int N>
inline bool operator!=(lemon::Invalid, MultiCoordinateIterator<N> const & i)
{
    return i.isValid();
}

#define VIGRA_LEMON_INVALID_COMPARISON(type) \
template<unsigned int N, bool BackEdgesOnly> \
inline bool operator==(type<N, BackEdgesOnly> const & i, lemon::Invalid) \
{ \
    return i.atEnd(); \
} \
template<unsigned int N, bool BackEdgesOnly> \
inline bool operator!=(type<N, BackEdgesOnly> const & i, lemon::Invalid) \
{ \
    return i.isValid(); \
} \
template<unsigned int N, bool BackEdgesOnly> \
inline bool operator==(lemon::Invalid, type<N, BackEdgesOnly> const & i) \
{ \
    return i.atEnd(); \
} \
template<unsigned int N, bool BackEdgesOnly> \
inline bool operator!=(lemon::Invalid, type<N, BackEdgesOnly> const & i) \
{ \
    return i.isValid(); \
}

VIGRA_LEMON_INVALID_COMPARISON(GridGraphNeighborIterator)
VIGRA_LEMON_INVALID_COMPARISON(GridGraphOutEdgeIterator)
VIGRA_LEMON_INVALID_COMPARISON(GridGraphOutArcIterator)
VIGRA_LEMON_INVALID_COMPARISON(GridGraphInArcIterator)
VIGRA_LEMON_INVALID_COMPARISON(GridGraphEdgeIterator)
VIGRA_LEMON_INVALID_COMPARISON(GridGraphArcIterator)

#undef VIGRA_LEMON_INVALID_COMPARISON

using boost_graph::directed_tag;
using boost_graph::undirected_tag;

namespace detail {

template <unsigned int N, class DirectedTag>
struct GridGraphBase;

template <unsigned int N>
struct GridGraphBase<N, directed_tag>
{
    template <class T>
    class ArcMap
    : public MultiArray<N+1, Multiband<T> >
    {
      public:
        typedef MultiArray<N+1, Multiband<T> >             base_type;
        typedef typename base_type::difference_type        difference_type;
        typedef typename base_type::key_type               key_type;
        typedef typename base_type::value_type             value_type; 
        typedef typename base_type::reference              reference;
        typedef typename base_type::const_reference        const_reference;
        typedef boost_graph::read_write_property_map_tag   category;
        
        typedef lemon::True                                ReferenceMapTag;
        typedef key_type                                   Key;
        typedef value_type                                 Value;
        typedef reference                                  Reference;
        typedef const_reference                            ConstReference;

        ArcMap()
        : base_type()
        {}
        
        explicit ArcMap(GridGraph<N, directed_tag> const & g)
        : base_type(g.arc_propmap_shape())
        {}
        
        ArcMap(GridGraph<N, directed_tag> const & g, T const & t)
        : base_type(g.arc_propmap_shape(), t)
        {}

        explicit ArcMap(difference_type const & shape)
        : base_type(shape)
        {}
        
        ArcMap(difference_type const & shape, T const & t)
        : base_type(shape, t)
        {}
        
        ArcMap & operator=(ArcMap const & m)
        {
            base_type::operator=(m);
            return *this;
        }
        
        ArcMap & operator=(base_type const & m)
        {
            base_type::operator=(m);
            return *this;
        }
        
        // appropriate operator[] are inherited
        
        void set(Key const & k, Value const & v)
        {
            (*this)[k] = v;
        }
    };
};

template <unsigned int N>
struct GridGraphBase<N, undirected_tag>
{
    typedef lemon::True UndirectedTag;
    
    template <class T>
    class ArcMap
    : public MultiArray<N+1, Multiband<T> >
    {
      public:
        typedef MultiArray<N+1, Multiband<T> >             base_type;
        typedef GridGraphArcDescriptor<N>                  difference_type;
        typedef difference_type                            key_type;
        typedef typename base_type::value_type             value_type; 
        typedef typename base_type::reference              reference;
        typedef typename base_type::const_reference        const_reference;
        typedef boost_graph::read_write_property_map_tag   category;
        
        typedef lemon::True                                ReferenceMapTag;
        typedef key_type                                   Key;
        typedef value_type                                 Value;
        typedef reference                                  Reference;
        typedef const_reference                            ConstReference;
        
        ArcMap()
        : base_type()
        {}
        
        explicit ArcMap(GridGraph<N, undirected_tag> const & g)
        : base_type(g.arc_propmap_shape()),
          graph_(&g)
        {}
        
        ArcMap(GridGraph<N, undirected_tag> const & g, T const & t)
        : base_type(g.arc_propmap_shape(), t),
          graph_(&g)
        {}
        
        ArcMap & operator=(ArcMap const & m)
        {
            base_type::operator=(m);
            return *this;
        }
        
        ArcMap & operator=(base_type const & m)
        {
            base_type::operator=(m);
            return *this;
        }
        
        reference operator[](difference_type const & s)
        {
            if(s.isReversed())
            {
                return base_type::operator[](graph_->directedArc(s));
            }
            else
            {
                return base_type::operator[](s);
            }
        }
        
        const_reference operator[](difference_type const & s) const
        {
            if(s.isReversed())
            {
                return base_type::operator[](graph_->directedArc(s));
            }
            else
            {
                return base_type::operator[](s);
            }
        }
        
        void set(Key const & k, Value const & v)
        {
            (*this)[k] = v;
        }
        
        GridGraph<N, undirected_tag> const * graph_;
    };
};

} // namespace detail

/********************************************************/
/*                                                      */
/*                      GridGraph                       */
/*                                                      */
/********************************************************/

/** \brief Define a grid graph in arbitrary dimensions.

    A GridGraph connects each pixel to its direct or indirect neighbors. 
    Direct neighbors are the adjacent point along the coordinate axes, 
    whereas indirect neighbors include the diagonally adjacent points. Thus, 
    direct neighbors correspond to the 4-neighborhood in 2D and the 
    6-neighborhood in 3D, whereas indirect neighbors correspond to the 
    8-neighborhood and 26-neighborhood respectively. The GridGraph class 
    extends this scheme to arbitrary dimensions. While the dimension must be 
    defined at compile time via the template parameter <tt>N</tt>, the 
    desired neighborhood can be chosen at runtime in the GridGraph's 
    constructor. The shape of the grid is also specified at runtime in terms 
    of a suitable shape object. 

    Another choice to be made at compile time is whether the graph should be directed 
    or undirected. This is defined via the <tt>DirectedTag</tt> template parameter
    which can take the values <tt>directed_tag</tt> or <tt>undirected_tag</tt> (default).
    
    The main difficulty in a grid graph is to skip the missing neighbors
    of the pixels near the grid's border. For example, the upper left pixel in a 
    2D grid has only two direct neighbors instead of the usual four. The GridGraph
    class uses a precomputed set of internal look-up tables to efficiently determine the 
    appropriate number and location of the existing neighbors. A key design decision 
    to make this fast was to give up the requirement that edge IDs are contiguous 
    integers (as in LEMON's implementation of a 2D grid graph, which became very
    complicated and slow by enforcing this requirement). Instead, edges are numbered 
    as if all nodes (including nodes at the grid's border) had the same degree.
    Since edges crossing the border are not actually present in the graph, these IDs 
    are missing, leading to gaps in the sequence of IDs. 
    
    For the sake of compatibility, the GridGraph class implements two 
    popular graph APIs: the one defined by the <a href="http://www.boost.org/doc/libs/release/libs/graph/">boost::graph</a> library and
    the one defined by the <a href="http://lemon.cs.elte.hu/">LEMON</a> library.
    Their basic principles are very similar: The graph's nodes, edges and adjacency
    structure are accessed via a set of iterators, whereas additional properties
    (like node and edge weights) are provided via <i>property maps</i> that are
    indexed with suitable node or edge descriptor objects returned by the iterators.
    
    Specifically, GridGraph implements the requirements of the following <a href="http://www.boost.org/doc/libs/release/libs/graph/doc/graph_concepts.html">concepts of the 
    boost::graph API</a>: <b>Graph, IncidenceGraph, BidirectionalGraph, AdjacencyGraph, 
    VertexListGraph, EdgeListGraph,</b> and <b>AdjacencyMatrix</b>. Likewise, it supports 
    the concepts <b>Graph</b> and <b>Digraph</b> of the LEMON API. The property maps 
    associated with a GridGraph support the boost concepts ReadablePropertyMap, 
    WritablePropertyMap, ReadWritePropertyMap, and LvaluePropertyMap as well as the 
    LEMON concepts ReadMap, WriteMap, ReadWriteMap, and ReferenceMap.
    
    VIGRA's GridGraph class is designed such that multi-dimensional coordinates
    (i.e. <tt>TinyVector<MultiArrayIndex></tt>) serve as descriptor objects, which means
    that normal <tt>MultiArray</tt>s or <tt>MultiArrayView</tt>s can serve as property
    maps in most situations. Thus, node properties like a foreground probability for 
    foreground/background segmentation can simply be stored as normal images.
    
    Since the boost::graph and LEMON APIs differ in how they call corresponding
    functionality (e.g., they use the terms <tt>vertex</tt> and <tt>node</tt> respectively
    in an exactly synonymous way), most GridGraph helper classes and functions are exposed 
    under two different names. To implement your own algorithms, you can choose the API 
    you like most. VIGRA adopts the convention that algorithms using the boost::graph 
    API go into the namespace <tt>vigra::boost_graph</tt>, while those using the 
    LEMON API are placed into the namespace <tt>vigra::lemon_graph</tt>. This helps 
    to avoid name clashes when the two APIs happen to use the same name for different 
    things. The documentation of the GridGraph members specifies which API the respective
    function or class belongs to. Please consult the documentation of these 
    libraries for tutorial introductions and full reference of the respective APIs.
    
    VIGRA adds an important new concept of its own: the back neighborhood
    and associated adjacency and edge iterators. The back neighborhood of a given vertex 
    with ID <tt>i</tt> is the set of all neighbor vertices whose ID is smaller than <tt>i</tt>.
    This concept is useful if you work on the grid graph's vertices in scan order
    and want to access just those neighbors that have already been processed. Connected
    components labeling is a typical example of an algorithm that can take advantage of this 
    possibility. In principle, a back neighborhood iterator could be implemented as
    a filter iterator that simply skips all neighbors with inappropriate IDs. However,
    in case of grid graphs it is more efficient to provide a special iterator for this purpose.
    
    <b>Usage:</b>

    <b>\#include</b> \<vigra/multi_gridgraph.hxx\> <br/>
    Namespace: vigra
    
    At present, the GridGraph class is mainly used internally to implement image 
    analysis functions for arbitrary dimensional arrays (e.g. detection of local 
    extrema, connected components labeling, watersheds, SLIC superpixels). For example,
    a dimension-independent algorithm to detect local maxima using the LEMON API
    might look like this:
    
    \code
    namespace vigra { namespace lemon_graph { 

    template <class Graph, class InputMap, class OutputMap>
    void
    localMaxima(Graph const & g, 
                InputMap const & src,
                OutputMap & local_maxima,
                typename OutputMap::value_type marker)
    {
        // define abreviations for the required iterators
        typedef typename Graph::NodeIt    graph_scanner;
        typedef typename Graph::OutArcIt  neighbor_iterator;

        // iterate over all nodes (i.e. pixels)
        for (graph_scanner node(g); node != INVALID; ++node) 
        {
            // remember the value of the current node
            typename InputMap::value_type current = src[*node];

            // iterate over all neighbors of the current node
            bool is_local_maximum = true;
            for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
            {
                // if a neighbor has larger value, the center node is not a local maximum
                if (current < src[g.target(*arc)]) 
                    is_local_maximum = false;
            }
            
            // mark the node when it is a local maximum
            if (is_local_maximum)
                local_maxima[*node] = marker;
        }
    }
    
    }} // namespace vigra::lemon_graph
    \endcode
    
    It should be noted that this algorithm will work for any LEMON-compatible graph class, 
    not just our GridGraph. When implemented in terms of the boost::graph API, the same algorithm
    looks like this:
    
    \code
    namespace vigra { namespace boost_graph {

    template <class Graph, class InputMap, class OutputMap>
    void
    localMaxima(Graph const & g, 
                InputMap const & src,
                OutputMap & local_maxima,
                typename property_traits<OutputMap>::value_type marker)
    {
        // define abreviations and variables for the required iterators
        typedef typename graph_traits<Graph>::vertex_iterator graph_scanner;
        typedef typename graph_traits<Graph>::adjacency_iterator neighbor_iterator;

        graph_scanner      node, node_end;
        neighbor_iterator  arc, arc_end;

        // iterate over all nodes (i.e. pixels)
        tie(node, node_end) = vertices(g);
        for (; node != node_end; ++node) 
        {
            // remember the value of the current node
            typename property_traits<InputMap>::value_type current = get(src, *node);

            // iterate over all neighbors of the current node
            bool is_local_maximum = true;
            tie(arc, arc_end) = adjacent_vertices(*node, g);
            for (;arc != arc_end; ++arc)
            {
                // if a neighbor has larger value, the center node is not a local maximum
                if (current < get(src, *arc))
                    is_local_maximum = false;
            }
            
            // mark the node when it is a local maximum
            if (is_local_maximum)
                put(local_maxima, *node, marker);
        }
    }

    }} // namespace vigra::boost_graph
    \endcode
    
    It can be seen that the differences between the APIs are mainly syntactic 
    (especially note that boost::graph users traits classes and free functions, 
    whereas LEMON uses nested typedefs and member functions). Either of these 
    algorithms can now serve as the backend of a local maxima detector 
    for <tt>MultiArrayViews</tt>:
    
    \code
    namespace vigra {
    
    template <unsigned int N, class T1, 
                              class T2>
    void
    localMaxima(MultiArrayView<N, T1> const & src,
                MultiArrayView<N, T2> local_maxima,
                T2 marker,
                NeighborhoodType neighborhood = IndirectNeighborhood)
    {
        vigra_precondition(src.shape() == local_maxima.shape(),
            "localMinMax(): shape mismatch between input and output.");
        
        // create a grid graph with appropriate shape and desired neighborhood
        GridGraph<N, undirected_tag> graph(src.shape(), neighborhood);
        
        // forward the call to the graph-based algorithm, using 
        // the given MultiArrayViews as property maps
        lemon_graph::localMaxima(graph, src, local_maxima, marker);
    }

    } // namespace vigra
    \endcode
    
    A slightly enhanced version of this code is actually used to implement this
    functionality in VIGRA.
*/
template<unsigned int N, class DirectedTag=undirected_tag>
class GridGraph
: public detail::GridGraphBase<N, DirectedTag>
{
public:
        /** \brief 'true' if the graph is directed (API: boost::graph)
        */
    static const bool is_directed = IsSameType<DirectedTag, directed_tag>::value;
    
    typedef detail::GridGraphBase<N, DirectedTag>   base_type;
    typedef GridGraph<N, DirectedTag>               self_type;
    
        /** \brief Dimension of the grid.
        */
    static const unsigned int dimension = N;
    
        /** \brief Shape type of the graph and a node property map.
        */
    typedef typename MultiArrayShape<N>::type       shape_type;
    
        /** \brief Shape type of an edge property map (must have one additional dimension).
        */
    typedef typename MultiArrayShape<N+1>::type     edge_propmap_shape_type;    
    
        /** \brief Type of node and edge IDs.
        */
    typedef MultiArrayIndex                         index_type;
    
        /** \brief Type to specify number of vertices (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::vertices_size_type</tt>).
        */
    typedef MultiArrayIndex                         vertices_size_type;
    
        /** \brief Type to specify number of edges (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::edges_size_type</tt>).
        */
    typedef MultiArrayIndex                         edges_size_type;
    
        /** \brief Type to specify number of neighbors (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::degree_size_type</tt>).
        */
    typedef MultiArrayIndex                         degree_size_type;

        /** \brief Iterator over the vertices of the graph (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::vertex_iterator</tt>).
        */
    typedef MultiCoordinateIterator<N>              vertex_iterator;
    
        /** \brief Iterator over the neighbor vertices of a given vertex (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::adjacency_iterator</tt>).
        */
    typedef GridGraphNeighborIterator<N>            adjacency_iterator;
    
        /** \brief Same as adjacency_iterator (API: VIGRA).
        */
    typedef adjacency_iterator                      neighbor_vertex_iterator;
    
        /** \brief Iterator over only those neighbor vertices of a given vertex 
                   that have smaller ID (API: VIGRA).
        */
    typedef GridGraphNeighborIterator<N, true>      back_neighbor_vertex_iterator;
    
        /** \brief Iterator over the incoming edges of a given vertex (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::in_edge_iterator</tt>).
        */
    typedef GridGraphInArcIterator<N>               in_edge_iterator;
    
        /** \brief Iterator over the outgoing edges of a given vertex (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::out_edge_iterator</tt>).
        */
    typedef GridGraphOutArcIterator<N>              out_edge_iterator;
    
        /** \brief Iterator over only those outgoing edges of a given vertex
                   that go to vertices with smaller IDs (API: VIGRA).
        */
    typedef GridGraphOutArcIterator<N, true>        out_back_edge_iterator;
    
        /** \brief Iterator over the edges of a graph (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::edge_iterator</tt>).
        */
    typedef GridGraphArcIterator<N, !is_directed>   edge_iterator;
    
        /** \brief The vertex descriptor (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::vertex_descriptor</tt>).
        */
    typedef shape_type                              vertex_descriptor;
    
        /** \brief The edge descriptor (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::edge_descriptor</tt>).
        */
    typedef GridGraphArcDescriptor<N>               edge_descriptor;
     
        /** \brief Is the graph directed or undirected ? (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::directed_category</tt>).
        */
    typedef DirectedTag                             directed_category;
    
        /** \brief The graph does not allow multiple edges between the same vertices 
                  (API: boost::graph, use via 
                  <tt>boost::graph_traits<Graph>::edge_parallel_category</tt>).
        */
    typedef boost_graph::disallow_parallel_edge_tag  edge_parallel_category;
    
        /** \brief The graph does not define internal property maps (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::vertex_property_type</tt>).
        */
    typedef boost_graph::no_property                 vertex_property_type; // we only support "external properties".
    // FIXME: Maybe support the vertex -> coordinate map (identity) as the only internal property map
    // and additionally the vertex_descriptor -> ID map (vertex_index = SOI).

        /** \brief Define several graph tags related to graph traversal (API: boost::graph,
             use via <tt>boost::graph_traits<Graph>::traversal_category</tt>).
        */
    struct traversal_category 
    : virtual public boost_graph::bidirectional_graph_tag,
      virtual public boost_graph::adjacency_graph_tag,
      virtual public boost_graph::vertex_list_graph_tag,
      virtual public boost_graph::edge_list_graph_tag,
      virtual public boost_graph::adjacency_matrix_tag
    {};
    
        // internal types
    typedef ArrayVector<shape_type>                      NeighborOffsetArray;
    typedef ArrayVector<NeighborOffsetArray>             RelativeNeighborOffsetsArray;
    typedef ArrayVector<ArrayVector<edge_descriptor> >   RelativeEdgeOffsetsArray;
    typedef ArrayVector<ArrayVector<MultiArrayIndex> >   IndexArray;
    typedef ArrayVector<ArrayVector<bool> >              NeighborExistsArray;
        


    ////////////////////////////////////////////////////////////////////
    
    // LEMON interface
    typedef self_type                               Graph;
    
        /** \brief The Node descriptor (API: LEMON).
        */
    typedef vertex_descriptor                       Node;
    
        /** \brief Iterator over all nodes of the graph (API: LEMON).
        */
    typedef vertex_iterator                         NodeIt;
    
        /** \brief The arc (directed edge) descriptor (API: LEMON).
        */
    typedef GridGraphArcDescriptor<N>               Arc;
    
        /** \brief Iterator over the outgoing edges of a node (API: LEMON).
        */
    typedef GridGraphOutArcIterator<N>              OutArcIt;
    
        /** \brief Iterator over only those outgoing edges of a node
            that end in a node with smaller ID (API: VIGRA).
        */
    typedef GridGraphOutArcIterator<N, true>        OutBackArcIt;
    
        /** \brief Iterator over the acrs (directed edges) of a node (API: LEMON).
        */
    typedef GridGraphArcIterator<N, false>          ArcIt;
    
        /** \brief Iterator over the incoming arcs of a node (API: LEMON).
        */
    typedef GridGraphInArcIterator<N>               InArcIt;
    
        /** \brief The edge descriptor (API: LEMON).
        */
    typedef typename MultiArrayShape<N+1>::type     Edge;
    
        /** \brief Iterator over the incident edges of a node (API: LEMON).
        */
    typedef GridGraphOutEdgeIterator<N>             IncEdgeIt;
    
        /** \brief Iterator over only those incident edges of a node that
            end in a node with smaller ID (API: VIGRA).
        */
    typedef GridGraphOutEdgeIterator<N, true>       IncBackEdgeIt;
    
        /** \brief Iterator over the edges of the graph (API: LEMON).
        */
    typedef GridGraphEdgeIterator<N, !is_directed>  EdgeIt;
    
    typedef lemon::True NodeNumTag;
    typedef lemon::True EdgeNumTag;
    typedef lemon::True ArcNumTag;
    typedef lemon::True FindEdgeTag;
    typedef lemon::True FindArcTag;
    
    ////////////////////////////////////////////////////////////////////
    
        /** \brief Type of a node property map that maps node descriptor objects 
            onto property values of type <tt>T</tt> (API: LEMON).
        */
    template <class T>
    class NodeMap
    : public MultiArray<N, T>
    {
      public:
        typedef MultiArray<N, T> base_type;
        typedef typename base_type::difference_type        difference_type;
        typedef typename base_type::key_type               key_type;
        typedef typename base_type::value_type             value_type; 
        typedef typename base_type::reference              reference;
        typedef typename base_type::const_reference        const_reference;
        typedef boost_graph::read_write_property_map_tag   category;
        
        typedef lemon::True                                ReferenceMapTag;
        typedef key_type                                   Key;
        typedef value_type                                 Value;
        typedef reference                                  Reference;
        typedef const_reference                            ConstReference;

        NodeMap()
        : base_type()
        {}
        
            /** \brief Construct property map for the given graph \a g
                (preallocates a zero-initialized entry for each node of the graph).
            */
        explicit NodeMap(GridGraph const & g)
        : base_type(g.shape())
        {}
        
            /** \brief Construct property map for the given graph \a g
                (preallocates an entry with initial value \a t for each node of the graph).
            */
        NodeMap(GridGraph const & g, T const & t)
        : base_type(g.shape(), t)
        {}

            /** \brief Construct property map for the given \a shape.
                (preallocates a zero-initialized entry for each node of a grid
                graph with this shape).
            */
        explicit NodeMap(difference_type const & shape)
        : base_type(shape)
        {}
        
            /** \brief Construct property map for the given \a shape.
                (preallocates an entry with initial value \a t for each node of a grid
                graph with this shape).
            */
        NodeMap(difference_type const & shape, T const & t)
        : base_type(shape, t)
        {}
        
        NodeMap & operator=(NodeMap const & m)
        {
            base_type::operator=(m);
            return *this;
        }
        
        NodeMap & operator=(base_type const & m)
        {
            base_type::operator=(m);
            return *this;
        }
        
        // appropriate operator[] are inherited
#ifdef DOXYGEN
            /** \brief Read/write access to the value associated with node descriptor \a key.
            */
        Value & operator[](Key const & key);

            /** \brief Read-only access to the value associated with node descriptor \a key.
            */
        Value const & operator[](Key const & key) const;
#endif // DOXYGEN
        
            /** \brief Set the property of node desctiptor \a key to value \a v.
            */
        void set(Key const & k, Value const & v)
        {
            (*this)[k] = v;
        }
    };
    
    
        /** \brief Type of an edge property map that maps edge descriptor objects 
            onto property values of type <tt>T</tt> (API: LEMON).
        */
    template <class T>
    class EdgeMap
    : public MultiArray<N+1, Multiband<T> >
    {
      public:
        typedef MultiArray<N+1, Multiband<T> >             base_type;
        typedef typename base_type::difference_type        difference_type;
        typedef typename base_type::key_type               key_type;
        typedef typename base_type::value_type             value_type; 
        typedef typename base_type::reference              reference;
        typedef typename base_type::const_reference        const_reference;
        typedef boost_graph::read_write_property_map_tag   category;
        
        typedef lemon::True                                ReferenceMapTag;
        typedef key_type                                   Key;
        typedef value_type                                 Value;
        typedef reference                                  Reference;
        typedef const_reference                            ConstReference;

        EdgeMap()
        : base_type()
        {}
        
            /** \brief Construct property map for the given graph \a g
                (preallocates a zero-initialized entry for each edge of the graph).
            */
        explicit EdgeMap(GridGraph const & g)
        : base_type(g.edge_propmap_shape())
        {}
        
            /** \brief Construct property map for the given graph \a g
                (preallocates an entry with initial value \a t for each edge of the graph).
            */
        EdgeMap(GridGraph const & g, T const & t)
        : base_type(g.edge_propmap_shape(), t)
        {}

            /** \brief Construct property map for the given \a shape
                (preallocates a zero-initialized entry for each edge of 
                a grid graph with this shape).
            */
        explicit EdgeMap(difference_type const & shape)
        : base_type(shape)
        {}
        
            /** \brief Construct property map for the given \a shape
                (preallocates an entry with initial value \a t for each edge 
                of a grid graph with this shape).
            */
        EdgeMap(difference_type const & shape, T const & t)
        : base_type(shape, t)
        {}
        
        EdgeMap & operator=(EdgeMap const & m)
        {
            base_type::operator=(m);
            return *this;
        }
        
        EdgeMap & operator=(base_type const & m)
        {
            base_type::operator=(m);
            return *this;
        }
        
        // appropriate operator[] are inherited
#ifdef DOXYGEN
            /** \brief Read/write access to the value associated with edge descriptor \a key.
            */
        Value & operator[](Key const & key);

            /** \brief Read-only access to the value associated with edge descriptor \a key.
            */
        Value const & operator[](Key const & key) const;
#endif // DOXYGEN
        
            /** \brief Set the property of edge desctiptor \a key to value \a v.
            */
        void set(Key const & k, Value const & v)
        {
            (*this)[k] = v;
        }
    };

#ifdef DOXYGEN
        /** \brief Type of an arc property map that maps arc descriptor objects 
            onto property values of type <tt>T</tt> (API: LEMON).
        */
    template <class T>
    class ArcMap
    : public MultiArray<N+1, Multiband<T> >
    {
      public:
        
            /** \brief Construct property map for the given graph \a g
                (preallocates a zero-initialized entry for each arc of the graph).
            */
        explicit ArcMap(GridGraph const & g);
        
            /** \brief Construct property map for the given graph \a g
                (preallocates an entry with initial value \a t for each arc of the graph).
            */
        ArcMap(GridGraph const & g, T const & t);

            /** \brief Construct property map for the given \a shape
                (preallocates a zero-initialized entry for each arc of 
                a grid graph with this shape).
            */
        explicit ArcMap(difference_type const & shape);
        
            /** \brief Construct property map for the given \a shape
                (preallocates an entry with initial value \a t for each arc of 
                a grid graph with this shape).
            */
        ArcMap(difference_type const & shape, T const & t);

            /** \brief Read/write access to the value associated with arc descriptor \a key.
            */
        Value & operator[](Key const & key);

            /** \brief Read-only access to the value associated with arc descriptor \a key.
            */
        Value const & operator[](Key const & key) const;
        
            /** \brief Set the property of arc desctiptor \a key to value \a v.
            */
        void set(Key const & k, Value const & v);
    };
#endif // DOXYGEN

        /** \brief Type of a property map that returns the coordinate of a given node (API: LEMON).
        */
    class IndexMap 
    {
      public:
        typedef Node                                    Key;
        typedef Node                                    Value;
        typedef Key                                     key_type;
        typedef Value                                   value_type; 
        typedef Value const &                           reference;
        typedef boost_graph::readable_property_map_tag  category;

        IndexMap()
        {}

            /** \brief Construct property map for the given graph.
            */
        explicit IndexMap(const GridGraph&)
        {}

            /** \brief Get the grid coordinate of the node descriptor \a key.
            */
        Value const & operator[](Key const & key) const 
        {
            return key;
        }
    };
        
        /** \brief Type of a property map that returns the number of incoming edges of a given node (API: LEMON, use via <tt>lemon::InDegMap<Graph></tt>).
        */
    class InDegMap
    {
      public:
      
        /// The graph type of InDegMap (works for directed and undirected graphs)
        typedef GridGraph                               Graph;
        typedef GridGraph                               Digraph;
        typedef Node                                    Key;
        typedef degree_size_type                        Value;
        typedef Key                                     key_type;
        typedef Value                                   value_type; 
        typedef Value const &                           reference;
        typedef boost_graph::readable_property_map_tag  category;
        
            /** \brief Construct property map for the given graph.
            */
        explicit InDegMap(const GridGraph& graph)
        : graph_(graph)
        {}
        
            /** \brief Get the in-degree of the node descriptor \a key.
            */
        Value operator[](const Key& key) const 
        {
            return graph_.in_degree(key);
        }
        
      protected:
      
        GridGraph const & graph_;
    };

        /** \brief Type of a property map that returns the number of outgoing edges of a given node (API: LEMON, use via <tt>lemon::OutDegMap<Graph></tt>).
        */
    class OutDegMap
    {
      public:
      
        /// The graph type of OutDegMap (works for directed and undirected graphs)
        typedef GridGraph                               Graph;
        typedef GridGraph                               Digraph;
        typedef Node                                    Key;
        typedef degree_size_type                        Value;
        typedef Key                                     key_type;
        typedef Value                                   value_type; 
        typedef Value const &                           reference;
        typedef boost_graph::readable_property_map_tag  category;
        
            /** \brief Construct property map for the given graph.
            */
        explicit OutDegMap(const GridGraph& graph)
        : graph_(graph)
        {}
        
            /** \brief Get the out-degree of the node descriptor \a key.
            */
        Value operator[](const Key& key) const 
        {
            return graph_.out_degree(key);
        }
        
      
        GridGraph const & graph_;
    };
    
    ////////////////////////////////////////////////////////////////////
    
        // dummy default constructor to satisfy adjacency_graph concept
    GridGraph()
    : max_node_id_(-1), 
      max_arc_id_(-1), 
      max_edge_id_(-1)
    {}
        
        /** \brief Construct a grid graph with given \a shape and neighborhood type \a ntype.

            The shape must have type <tt>MultiArrayShape<N>::type</tt> with the appropriate 
            dimension <tt>N</tt>. The neighborhood type can take the values 
            <tt>DirectNeighborhood</tt> to use only the axis-aligned edges (2N-neighborhood)
            and <tt>IndirectNeighborhood</tt> to use all diagonal edges as well 
            ((3<sup>N</sup>-1)-neighborhood).
        */
    GridGraph(shape_type const &shape, NeighborhoodType ntype = DirectNeighborhood) 
    : shape_(shape),
      num_vertices_(prod(shape)),
      num_edges_(gridGraphEdgeCount(shape, ntype, is_directed)), 
      max_node_id_(num_vertices_ - 1), 
      max_arc_id_(-2), 
      max_edge_id_(-2), 
      neighborhoodType_(ntype)
    {
        // populate the neighborhood tables:
        // FIXME: this might be static (but make sure that it works with multi-threading)
        detail::makeArrayNeighborhood(neighborOffsets_, neighborExists_, neighborhoodType_);
        detail::computeNeighborOffsets(neighborOffsets_, neighborExists_, incrementalOffsets_, 
                                       edgeDescriptorOffsets_, neighborIndices_, backIndices_, is_directed);
        
        // compute the neighbor offsets per neighborhood type
        // detail::makeArraySubNeighborhood(neighborhood[0], neighborExists, shape_type(1), neighborhoodIndices);
    }
    
        /** \brief Get the ID (i.e. scan-order index) for node desciptor \a v (API: LEMON).
        */
    index_type id(Node const & v) const
    {
        return detail::CoordinateToScanOrder<N>::exec(shape(), v);
    }

    index_type id(NodeIt const & v) const
    {
        return v.scanOrderIndex();
    }
    
    index_type id(neighbor_vertex_iterator const & v) const
    {
        return id(*v);
    }
    
    index_type id(back_neighbor_vertex_iterator const & v) const
    {
        return id(*v);
    }
    
        /** \brief Get node descriptor for given node ID \a i (API: LEMON).
        
            Return <tt>Node(lemon::INVALID)</tt> when the ID does not exist in this graph.
        */
    Node nodeFromId(index_type i) const
    {
        if(i < 0 || i > maxNodeId())
            return Node(lemon::INVALID);
        
        Node res(SkipInitialization);
        detail::ScanOrderToCoordinate<N>::exec(i, shape(), res);
        return res;
    }

        /** \brief Get the maximum ID of any node in this graph (API: LEMON).
        */
    index_type maxNodeId() const
    {
        return prod(shape()) - 1;
    }
    
        /** \brief Get the grid cordinate of the given node \a v (convenience function).
        */
    Node const & pos(Node const & v) const
    {
        return v;
    }
    
        /** \brief Get vertex iterator pointing to the first vertex of this graph (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::vertices(graph)</tt>,<br/>
            LEMON uses <tt>Graph::NodeIt(graph)</tt>).
        */
    vertex_iterator get_vertex_iterator() const 
    {
        return vertex_iterator(shape_);
    }

        /** \brief Get vertex iterator pointing to the given vertex (API: VIGRA).
        */
    vertex_iterator get_vertex_iterator(vertex_descriptor const & v) const 
    {
        return vertex_iterator(shape_) + v;
    }

        /** \brief Get vertex iterator pointing beyond the valid range of vertices of this graph (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::vertices(graph)</tt>,<br/>
            LEMON uses the special value <tt>lemon::INVALID</tt> instead).
        */
    vertex_iterator get_vertex_end_iterator() const 
    {
        return get_vertex_iterator().getEndIterator();
    }

        /** \brief Get an iterator pointing to the first neighbor of the given vertex (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::adjacent_vertices(v, graph)</tt>,<br/>
            LEMON uses <tt>Graph::ArcIt(g, v)</tt>).
        */
    neighbor_vertex_iterator get_neighbor_vertex_iterator(vertex_descriptor const & v) const 
    {
        return neighbor_vertex_iterator(*this, v);
    }

        /** \brief Get an iterator pointing beyond the range of neighbors of the given vertex (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::adjacent_vertices(v, graph)</tt>,<br/>
            LEMON uses the speical value <tt>lemon::INVALID</tt> instead).
        */
    neighbor_vertex_iterator get_neighbor_vertex_end_iterator(vertex_descriptor const & v) const 
    {
       return get_neighbor_vertex_iterator(v).getEndIterator();
    }

        /** \brief Get an iterator pointing to the first backward neighbor of the given vertex (API: VIGRA,<br/>
            in analogy to the boost::graph API, we also provide a free function <tt>boost::back_adjacent_vertices(v, g)</tt>,<br/>
            and the LEMON analogue is <tt>Graph::OutBackArcIt(graph, v)</tt>).
        */
    back_neighbor_vertex_iterator get_back_neighbor_vertex_iterator(vertex_descriptor const & v) const 
    {
        return back_neighbor_vertex_iterator(*this, v);
    }

        /** \brief Get an iterator pointing beyond the range of backward neighbors of the given vertex (API: VIGRA,<br/>
            in analogy to the boost::graph API, we also provide a free function <tt>boost::back_adjacent_vertices(v, g)</tt>,<br/>
            and LEMON just uses <tt>lemon::INVALID</tt> instead).
        */
    back_neighbor_vertex_iterator get_back_neighbor_vertex_end_iterator(vertex_descriptor const & v) const 
    {
       return get_back_neighbor_vertex_iterator(v).getEndIterator();
    }
    
    // --------------------------------------------------
    // support for VertexListGraph:

        /** \brief Get the number of vertices in this graph (convenience function,
            the boost::graph API provides the free function <tt>boost::num_vertices(graph)</tt>).
        */
    vertices_size_type num_vertices() const 
    {
        return num_vertices_;
    }

        /** \brief Get the number of nodes in this graph (API: LEMON).
        */
    vertices_size_type nodeNum() const 
    {
        return num_vertices();
    }

    // --------------------------------------------------
    // support for IncidenceGraph:

         /** \brief Get the ID (i.e. scan-order index in an edge property map) for the
             given edges descriptor \a e (API: LEMON).
        */
   index_type id(Edge const & e) const
    {
        return detail::CoordinateToScanOrder<N+1>::exec(edge_propmap_shape(), e);
    }
    
    index_type id(EdgeIt const & e) const
    {
        return id(*e);
    }
    
    index_type id(IncEdgeIt const & e) const
    {
        return id(*e);
    }

    index_type id(IncBackEdgeIt const & e) const
    {
        return id(*e);
    }

        /** \brief Get the edge descriptor for the given edge ID \a i (API: LEMON).
        
            Return <tt>Edge(lemon::INVALID)</tt> when the ID does not exist
            in this graph. 
        */
    Edge edgeFromId(index_type i) const
    {
        if(i < 0 || i > maxEdgeId())
            return Edge(lemon::INVALID);
        
        Edge res(SkipInitialization);
        detail::ScanOrderToCoordinate<N+1>::exec(i, edge_propmap_shape(), res);
        
        unsigned int b = detail::BorderTypeImpl<N>::exec(res.template subarray<0, N>(), shape());
        if(neighborExists_[b][res[N]])
            return res;
        else
            return Edge(lemon::INVALID);
    }
    
        /** \brief Get the maximum ID of any edge in this graph (API: LEMON).
        */
    index_type maxEdgeId() const
    {
        if(max_edge_id_ == -2) // -2 means uninitialized
            const_cast<GridGraph *>(this)->computeMaxEdgeAndArcId();
        return max_edge_id_;
    }
    
        /* Initial computation of the max_arc_id_ and max_edge_id_ (call in the constructor and
           whenever the shape changes).
        */
    void computeMaxEdgeAndArcId()
    {
        if(edgeNum() == 0)
        {
            max_arc_id_ = -1;
            max_edge_id_ = -1;
        }
        else
        {
            Node lastNode = shape() - shape_type(1);
            index_type n = neighborIndices_[get_border_type(lastNode)][0];
            Arc a(neighbor(lastNode, n), oppositeIndex(n), false);
            max_arc_id_ = detail::CoordinateToScanOrder<N+1>::exec(arc_propmap_shape(), a);
            
            if(is_directed)
            {
                max_edge_id_ = max_arc_id_;
            }
            else
            {
                Arc a(lastNode, backIndices_[get_border_type(lastNode)].back(), false);
                max_edge_id_ = detail::CoordinateToScanOrder<N+1>::exec(edge_propmap_shape(), a);
            }
        }
    }
    
        /** \brief Get the ID (i.e. scan-order index an an arc property map) for 
            the given ar \a a (API: LEMON).
        */
    index_type id(Arc const & a) const
    {
        return detail::CoordinateToScanOrder<N+1>::exec(arc_propmap_shape(), directedArc(a));
    }
    
    index_type id(ArcIt const & a) const
    {
        return id(*a);
    }
    
    index_type id(OutArcIt const & a) const
    {
        return id(*a);
    }
    
    index_type id(OutBackArcIt const & a) const
    {
        return id(*a);
    }
        
        /** \brief Get an arc descriptor for the given arc ID \a i (API: LEMON).
        
            Return <tt>Arc(lemon::INVALID)</tt> when the ID does not exist
            in this graph. 
        */
    Arc arcFromId(index_type i) const
    {
        if(i < 0 || i > maxArcId())
            return Arc(lemon::INVALID);
        
        Arc res;
        detail::ScanOrderToCoordinate<N+1>::exec(i, arc_propmap_shape(), res);
        unsigned int b = detail::BorderTypeImpl<N>::exec(res.template subarray<0, N>(), shape());
        if(neighborExists_[b][res[N]])
            return undirectedArc(res);
        else
            return Arc(lemon::INVALID);
    }
    
        /** \brief Get the maximal ID af any arc in this graph (API: LEMON).
        */
    index_type maxArcId() const
    {
        if(max_arc_id_ == -2) // -2 means uninitialized
            const_cast<GridGraph *>(this)->computeMaxEdgeAndArcId();
        return max_arc_id_;
    }
    
        /** \brief Return <tt>true</tt> when the arc is looking on the underlying
            edge in its natural (i.e. forward) direction, <tt>false</tt> otherwise (API: LEMON).
        */
    bool direction(Arc const & a) const
    {
        return !a.isReversed();
    }
    
        /** \brief Create an arc for the given edge \a e, oriented along the 
            edge's natural (<tt>forward = true</tt>) or reversed 
            (<tt>forward = false</tt>) direction (API: LEMON).
        */
    Arc direct(Edge const & e, bool forward) const
    {
        if(!is_directed || forward)
            return Arc(e, !forward);
        else
            return Arc(v(e), oppositeIndex(e[N]), true);
    }
    
        /** \brief Create an arc for the given edge \a e oriented
            so that node \a n is the starting node of the arc (API: LEMON), or
            return <tt>lemon::INVALID</tt> if the edge is not incident to this node.
        */
    Arc direct(Edge const & e, Node const & n) const
    {
        if(u(e) == n)
            return direct(e, true);
        if(v(e) == n)
            return direct(e, false);
        return Arc(lemon::INVALID);
    }
    
        /** \brief Return the opposite node of the given node \a n
            along edge \a e (API: LEMON), or return <tt>lemon::INVALID</tt>
            if the edge is not incident to this node.
        */
    Node oppositeNode(Node const & n, Edge const & e) const
    {
        Node start(u(e)), end(v(e));
        if(n == start)
            return end;
        if(n == end)
            return start;
        return Node(lemon::INVALID);
    }
    
        /** \brief Create an arc referring to the same edge as the given 
            arc \a a, but with reversed direction (API: LEMON).
        */
    Arc oppositeArc(Arc const & a) const
    {
        return is_directed
                 ? Arc(neighbor(a.vertexDescriptor(), a.edgeIndex()), oppositeIndex(a.edgeIndex()), false)
                 : Arc(a, !a.isReversed());
    }
    
        // internal function
        // transforms the arc into its directed form (i.e. a.isReversed() is 
        // guaranteed to be false in the returned arc).
    Arc directedArc(Arc const & a) const
    {
        return a.isReversed()
                 ? Arc(neighbor(a.vertexDescriptor(), a.edgeIndex()), oppositeIndex(a.edgeIndex()), false)
                 : a;
    }
    
        // internal function
        // transforms the arc into its undirected form (i.e. a.isReversed() will 
        // be true in the returned arc if this graph is undirected and the arc
        // traverses the edge backwards).
    Arc undirectedArc(Arc const & a) const
    {
        return a.edgeIndex() < maxUniqueDegree() 
                 ? a
                 : Arc(neighbor(a.vertexDescriptor(), a.edgeIndex()), oppositeIndex(a.edgeIndex()), true);
    }
    
        /** \brief Return the start node of the edge the given iterator is referring to (API: LEMON).
        */
    Node baseNode(IncEdgeIt const & e)  const
    {
        return source(e.arcDescriptor());
    }
    
        /** \brief Return the start node of the edge the given iterator is referring to (API: VIGRA).
        */
    Node baseNode(IncBackEdgeIt const & e)  const
    {
        return source(e.arcDescriptor());
    }
    
        /** \brief Return the start node of the edge the given iterator is referring to (API: LEMON).
        */
    Node baseNode(OutArcIt const & a)  const
    {
        return source(*a);
    }
    
        /** \brief Return the start node of the edge the given iterator is referring to (API: VIGRA).
        */
    Node baseNode(OutBackArcIt const & a)  const
    {
        return source(*a);
    }

        /** \brief Return the end node of the edge the given iterator is referring to (API: LEMON).
        */
    Node runningNode(IncEdgeIt const & e)  const
    {
        return target(e.arcDescriptor());
    }
    
        /** \brief Return the end node of the edge the given iterator is referring to (API: VIGRA).
        */
    Node runningNode(IncBackEdgeIt const & e)  const
    {
        return target(e.arcDescriptor());
    }
    
        /** \brief Return the end node of the edge the given iterator is referring to (API: LEMON).
        */
    Node runningNode(OutArcIt const & a)  const
    {
        return target(*a);
    }
    
        /** \brief Return the end node of the edge the given iterator is referring to (API: VIGRA).
        */
    Node runningNode(OutBackArcIt const & a)  const
    {
        return target(*a);
    }
    
        /** \brief Get the start node of the given arc \a a (API: LEMON).
        */
    Node source(Arc const & a) const 
    {
        return source_or_target(a, true);
    }

        /** \brief Get the end node of the given arc \a a (API: LEMON).
        */
    Node target(Arc const & a) const 
    {
        return source_or_target(a, false);
    }

        /** \brief Get the start node of the given edge \a e (API: LEMON,<br/>
            the boost::graph API provides the free function <tt>boost::source(e, graph)</tt>).
        */
    Node u(Edge const & e) const 
    {
        return Node(e.template subarray<0,N>());
    }

        /** \brief Get the end node of the given edge \a e (API: LEMON,<br/>
            the boost::graph API provides the free function <tt>boost::target(e, graph)</tt>).
        */
    Node v(Edge const & e) const 
    {
        return Node(e.template subarray<0,N>()) + neighborOffsets_[e[N]];
    }

        /** \brief Get an iterator pointing to the first outgoing edge of the given vertex (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::out_edges(v, graph)</tt>,<br/>
            LEMON uses <tt>Graph::OutArcIt(g, v)</tt>).
        */
    out_edge_iterator get_out_edge_iterator(vertex_descriptor const & v) const 
    {
        return out_edge_iterator(*this, v);
    }

        /** \brief Get an iterator pointing beyond the range of outgoing edges of the given vertex (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::out_edges(v, graph)</tt>,<br/>
            LEMON uses the special value <tt>lemon::INVALID</tt> instead).
        */
    out_edge_iterator get_out_edge_end_iterator(vertex_descriptor const & v) const 
    {
        return get_out_edge_iterator(v).getEndIterator();
    }

        /** \brief Get an iterator pointing to the first outgoing backward edge of the given vertex (API: VIGRA,<br/>
            in analogy to the boost::graph API, we also provide a free function <tt>boost::out_back_edges(v, g)</tt>,<br/>
            and the LEMON analogue is <tt>Graph::IncBackEdgeIt(graph, v)</tt>).
        */
    out_back_edge_iterator get_out_back_edge_iterator(vertex_descriptor const & v) const 
    {
        return out_back_edge_iterator(*this, v);
    }

        /** \brief Get an iterator pointing beyond the range of outgoing backward edges of the given vertex (API: VIGRA,<br/>
            in analogy to the boost::graph API, we also provide a free function <tt>boost::out_back_edges(v, g)</tt>,<br/>
            and LEMON uses the special value <tt>lemon::INVALID</tt> instead).
        */
    out_back_edge_iterator get_out_back_edge_end_iterator(vertex_descriptor const & v) const 
    {
        return get_out_back_edge_iterator(v).getEndIterator();
    }

        /** \brief Get an iterator pointing to the first incoming edge of the given vertex (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::in_edges(v, graph)</tt>,<br/>
            LEMON uses <tt>Graph::InArcIt(g, v)</tt>).
        */
    in_edge_iterator get_in_edge_iterator(vertex_descriptor const & v) const 
    {
        return in_edge_iterator(*this, v);
    }

        /** \brief Get an iterator pointing beyond the range of incoming edges of the given vertex (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::in_edges(v, graph)</tt>,<br/>
            LEMON uses the special value <tt>lemon::INVALID</tt> instead).
        */
    in_edge_iterator get_in_edge_end_iterator(vertex_descriptor const & v) const 
    {
        return get_in_edge_iterator(v).getEndIterator();
    }

        /** \brief Get the number of outgoing edges of the given vertex (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::out_degree(v, graph)</tt>,<br/>
            LEMON uses a special property map <tt>lemon::OutDegMap<Graph></tt>).
        */
    degree_size_type out_degree(vertex_descriptor const & v) const 
    {
        return (degree_size_type)neighborIndices_[get_border_type(v)].size();
    }
    
        /** \brief Get the number of outgoing backward edges of the given vertex (API: VIGRA).
        */
    degree_size_type back_degree(vertex_descriptor const & v) const 
    {
        return (degree_size_type)backIndices_[get_border_type(v)].size();
    }
    
        /** \brief Get the number of outgoing forward edges of the given vertex (API: VIGRA).
        */
    degree_size_type forward_degree(vertex_descriptor const & v) const 
    {
        unsigned int bt = get_border_type(v);
        return (degree_size_type)(neighborIndices_[bt].size() - backIndices_[bt].size());
    }
    
        /** \brief Get the number of incoming edges of the given vertex (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::in_degree(v, graph)</tt>,<br/>
            LEMON uses a special property map <tt>lemon::InDegMap<Graph></tt>).
        */
    degree_size_type in_degree(vertex_descriptor const & v) const 
    {
        return out_degree(v);
    }
    
        /** \brief Get the total number of edges (incoming plus outgoing) of the given vertex (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::degree(v, graph)</tt>,<br/>
            LEMON has no analogue).
        */
    degree_size_type degree(vertex_descriptor const & v) const 
    {
        return is_directed
                   ? 2*out_degree(v)
                   : out_degree(v);
    }
    
    // --------------------------------------------------
    // support for EdgeListGraph:

        /** \brief Get the number of edges in this graph (convenience function,
            boost::graph API provides the free function <tt>boost::num_edges(graph)</tt>).
        */
    edges_size_type num_edges() const 
    {
        return num_edges_;
    }

        /** \brief Get the number of edges in this graph (API: LEMON).
        */
    edges_size_type edgeNum() const 
    {
        return num_edges();
    }

        /** \brief Get the number of arc in this graph (API: LEMON).
        */
    edges_size_type arcNum() const 
    {
        return is_directed
                   ? num_edges()
                   : 2*num_edges();
    }
    
        /** \brief Get edge iterator pointing to the first edge of the graph (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::edges(graph)</tt>,<br/>
            LEMON uses <tt>Graph::EdgeIt(graph)</tt>).
        */
    edge_iterator get_edge_iterator() const 
    {
        return edge_iterator(*this);
    }

        /** \brief Get edge iterator pointing beyond the valid range of edges of this graph (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::vertices(graph)</tt>,<br/>
            LEMON uses the special value <tt>lemon::INVALID</tt> instead).
        */
    edge_iterator get_edge_end_iterator() const 
    {
        return get_edge_iterator().getEndIterator();
    }

    // --------------------------------------------------
    // support for AdjacencyMatrix concept:

        /** \brief Get a descriptor for the edge connecting vertices \a u and \a v,<br/>
            or <tt>(lemon::INVALID, false)</tt> if no such edge exists (convenience function,<br/>
            the boost::graph API provides the free function <tt>boost::edge(u, v, graph)</tt>).
        */
    std::pair<edge_descriptor, bool>
    edge(vertex_descriptor const & u, vertex_descriptor const & v) const
    {
        std::pair<edge_descriptor, bool> res(lemon::INVALID, false);

        neighbor_vertex_iterator i = get_neighbor_vertex_iterator(u),
                                 end = i.getEndIterator();
        for (; i != end; ++i) 
        {
            if (*i == v) 
            {
                res.first = make_edge_descriptor(u, i.neighborIndex());
                res.second = true;
                break;
            }
        }
        return res;
    }

        /** \brief Get a descriptor for the edge connecting vertices \a u and \a v,<br/>or <tt>lemon::INVALID</tt> if no such edge exists (API: LEMON).
        */
    Edge findEdge(Node const & u, Node const & v, Edge const & = lemon::INVALID) const 
    {
        return this->edge(u, v).first;
    }
    
        /** \brief Get a descriptor for the arc connecting vertices \a u and \a v,<br/>or <tt>lemon::INVALID</tt> if no such edge exists (API: LEMON).
        */
    Arc findArc(Node const & u, Node const & v, Arc const & = lemon::INVALID) const 
    {
        return this->edge(u, v).first;
        // std::pair<edge_descriptor, bool> res(edge(u, v));
        // return res.second
                 // ? res.first
                 // : Arc(lemon::INVALID);
    }
    
        /** \brief Create a property map that returns the coordinate of each node (API: LEMON GridGraph).
        */
    IndexMap indexMap() const 
    {
        return IndexMap();
    }

    // --------------------------------------------------
    // other helper functions:
    
    bool isDirected() const
    {
        return is_directed;
    }
    
    degree_size_type maxDegree() const
    {
        return (degree_size_type)neighborOffsets_.size();
    }

    degree_size_type maxUniqueDegree() const 
    {
         return is_directed
                    ? maxDegree()
                    : maxDegree() / 2;
    }

    shape_type const & shape() const 
    {
        return shape_;
    }

    edge_propmap_shape_type edge_propmap_shape() const 
    {
        edge_propmap_shape_type res(SkipInitialization);
        res.template subarray<0, N>() = shape_;
        res[N] = maxUniqueDegree();
        return res;
    }

    edge_propmap_shape_type arc_propmap_shape() const 
    {
        edge_propmap_shape_type res(SkipInitialization);
        res.template subarray<0, N>() = shape_;
        res[N] = maxDegree();
        return res;
    }

    unsigned int get_border_type(vertex_descriptor const & v) const
    {
        return detail::BorderTypeImpl<N>::exec(v, shape_);
    }

    unsigned int get_border_type(vertex_iterator const & v) const
    {
        return v.borderType();
    }
    
    index_type oppositeIndex(index_type neighborIndex) const
    {
        return  maxDegree() - neighborIndex - 1;
    }

        /* the given neighborIndex must be valid for the given vertex,
           otherwise this function will crash
        */
    edge_descriptor make_edge_descriptor(vertex_descriptor const & v,
                                         index_type neighborIndex) const
    {
        if(neighborIndex < maxUniqueDegree())
            return edge_descriptor(v, neighborIndex, false);
        else
            return edge_descriptor(neighbor(v, neighborIndex), oppositeIndex(neighborIndex), true);
    }
    
    shape_type const & neighborOffset(index_type neighborIndex) const
    {
        return neighborOffsets_[neighborIndex];
    }

    vertex_descriptor neighbor(vertex_descriptor const & v, index_type neighborIndex) const
    {
        return v + neighborOffsets_[neighborIndex];
    }

    vertex_descriptor 
    source_or_target(edge_descriptor const & e, bool return_source) const
    {
        // source is always the attached node (first coords) unless the
        // edge has been reversed. 
        if ((return_source && e.isReversed()) ||
            (!return_source && !e.isReversed())) 
        {
            return neighbor(e.vertexDescriptor(), e.edgeIndex());
        } 
        else 
        {
            return e.vertexDescriptor();
        }
    }
    
    NeighborOffsetArray const * neighborOffsetArray() const
    {
        return &neighborOffsets_;
    }
    
    RelativeNeighborOffsetsArray const * neighborIncrementArray() const
    {
        return &incrementalOffsets_;
    }
    
    RelativeEdgeOffsetsArray const * edgeIncrementArray() const
    {
        return &edgeDescriptorOffsets_;
    }
    
    IndexArray const * neighborIndexArray(bool backEdgesOnly) const
    {
        return backEdgesOnly 
                   ? &backIndices_
                   : &neighborIndices_;
    }

    NeighborExistsArray const * neighborExistsArray() const
    {
        return &neighborExists_;
    }

  protected:
    NeighborOffsetArray neighborOffsets_;
    NeighborExistsArray neighborExists_;
    IndexArray neighborIndices_, backIndices_;
    RelativeNeighborOffsetsArray incrementalOffsets_;
    RelativeEdgeOffsetsArray edgeDescriptorOffsets_;
    shape_type shape_;
    MultiArrayIndex num_vertices_, num_edges_, max_node_id_, max_arc_id_, max_edge_id_;
    NeighborhoodType neighborhoodType_;
};

template<unsigned int N, class DirectedTag>
inline
bool
isInside(GridGraph<N, DirectedTag> const & g,
         typename GridGraph<N, DirectedTag>::vertex_descriptor const & v) 
{
    return allLess(v, g.shape()) && allGreaterEqual(v, typename MultiArrayShape<N>::type());
}

//@}

#ifdef WITH_BOOST_GRAPH

} // namespace vigra

namespace boost {

#else

namespace boost_graph {

#endif 



template <unsigned int N, class T, class Acc>
struct property_traits<vigra::MultiArray<N, T, Acc> >
{
    typedef vigra::MultiArray<N, T, Acc>             type;
    typedef typename type::key_type                  key_type;
    typedef typename type::value_type                value_type; 
    typedef typename type::reference                 reference;
    typedef read_write_property_map_tag              category;
};

template <unsigned int N, class T, class Acc>
struct property_traits<vigra::MultiArray<N, T, Acc> const>
{
    typedef vigra::MultiArray<N, T, Acc> const       type;
    typedef typename type::key_type                  key_type;
    typedef typename type::value_type                value_type; 
    typedef typename type::const_reference           reference;
    typedef readable_property_map_tag                category;
};

template <unsigned int N, class T, class Stride>
struct property_traits<vigra::MultiArrayView<N, T, Stride> >
{
    typedef vigra::MultiArrayView<N, T, Stride>       type;
    typedef typename type::key_type                   key_type;
    typedef typename type::value_type                 value_type; 
    typedef typename type::reference                  reference;
    typedef read_write_property_map_tag               category;
};

template <unsigned int N, class T, class Stride>
struct property_traits<vigra::MultiArrayView<N, T, Stride> const>
{
    typedef vigra::MultiArrayView<N, T, Stride> const     type;
    typedef typename type::key_type                       key_type;
    typedef typename type::value_type                     value_type; 
    typedef typename type::const_reference                reference;
    typedef readable_property_map_tag                     category;
};

#ifdef WITH_BOOST_GRAPH

} // namespace boost

namespace vigra {
namespace boost_graph {

#endif 

/** \addtogroup BoostGraphExtensions GridGraph additions to namespace <tt>boost</tt>
        
        provides the required functionality to make \ref vigra::GridGraph compatible to the boost::graph library.
*/
//@{

    /** \brief Return number of outgoing edges of vertex \a v (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::degree_size_type
out_degree(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v, 
           vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.out_degree(v);
}

    /** \brief Return number of incoming edges of vertex \a v (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::degree_size_type
in_degree(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v, 
          vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.in_degree(v);
}

    /** \brief Return total number of edges (incoming and outgoing) of vertex \a v (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::degree_size_type
degree(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v, 
       vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.degree(v);
}

    /** \brief Get a (begin, end) iterator pair for the vertices of graph \a g (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
std::pair<typename vigra::GridGraph<N, DirectedTag>::vertex_iterator, 
          typename vigra::GridGraph<N, DirectedTag>::vertex_iterator>
vertices(vigra::GridGraph<N, DirectedTag> const & g) 
{
    return std::make_pair(g.get_vertex_iterator(),
                          g.get_vertex_end_iterator());    
}

    /** \brief Return the number of vertices in graph \a g (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::vertices_size_type
num_vertices(vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.num_vertices();
}

    /** \brief Get a (begin, end) iterator pair for the neighbor vertices of 
               vertex \a v in graph \a g (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
std::pair<typename vigra::GridGraph<N, DirectedTag>::adjacency_iterator, 
          typename vigra::GridGraph<N, DirectedTag>::adjacency_iterator>
adjacent_vertices(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v,
                  vigra::GridGraph<N, DirectedTag> const & g) 
{
    return std::make_pair(g.get_neighbor_vertex_iterator(v),
                          g.get_neighbor_vertex_end_iterator(v));    
}

    /** \brief Get a (begin, end) iterator pair for only thise neighbor vertices of 
               vertex \a v that have smaller ID than \a v (API: VIGRA).
    */
template<unsigned int N, class DirectedTag>
inline
std::pair<typename vigra::GridGraph<N, DirectedTag>::back_neighbor_vertex_iterator, 
          typename vigra::GridGraph<N, DirectedTag>::back_neighbor_vertex_iterator>
back_adjacent_vertices(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v,
                       vigra::GridGraph<N, DirectedTag> const & g) 
{
    return std::make_pair(g.get_back_neighbor_vertex_iterator(v),
                          g.get_back_neighbor_vertex_end_iterator(v));    
}

// adjacent_vertices variant in vigra namespace: allows to call adjacent_vertices with vertex_iterator argument
template<unsigned int N, class DirectedTag>
inline
std::pair<typename vigra::GridGraph<N, DirectedTag>::adjacency_iterator, 
          typename vigra::GridGraph<N, DirectedTag>::adjacency_iterator>
adjacent_vertices_at_iterator(typename vigra::GridGraph<N, DirectedTag>::vertex_iterator const & v,
                              vigra::GridGraph<N, DirectedTag> const & g) 
{    
    return std::make_pair(g.get_neighbor_vertex_iterator(v),
                          g.get_neighbor_vertex_end_iterator(v));    
}

    /** \brief Get a (begin, end) iterator pair for the outgoing edges of 
               vertex \a v in graph \a g (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
std::pair<typename vigra::GridGraph<N, DirectedTag>::out_edge_iterator, 
          typename vigra::GridGraph<N, DirectedTag>::out_edge_iterator>
out_edges(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v,
          vigra::GridGraph<N, DirectedTag> const & g) 
{
    return std::make_pair(g.get_out_edge_iterator(v),
                          g.get_out_edge_end_iterator(v));    
}

    /** \brief Get a (begin, end) iterator pair for only those outgoing edges of 
               vertex \a v whose ID is smaller than that of \a v (API: VIGRA).
    */
template<unsigned int N, class DirectedTag>
inline
std::pair<typename vigra::GridGraph<N, DirectedTag>::out_back_edge_iterator, 
          typename vigra::GridGraph<N, DirectedTag>::out_back_edge_iterator>
out_back_edges(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v,
               vigra::GridGraph<N, DirectedTag> const & g) 
{
    return std::make_pair(g.get_out_back_edge_iterator(v),
                          g.get_out_back_edge_end_iterator(v));    
}

    /** \brief Get a (begin, end) iterator pair for the incoming edges of 
               vertex \a v in graph \a g (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
std::pair<typename vigra::GridGraph<N, DirectedTag>::in_edge_iterator, 
          typename vigra::GridGraph<N, DirectedTag>::in_edge_iterator>
in_edges(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v,
         vigra::GridGraph<N, DirectedTag> const & g) 
{
    return std::make_pair(g.get_in_edge_iterator(v),
                          g.get_in_edge_end_iterator(v));    
}

    /** \brief Get a vertex descriptor for the start vertex of edge \a e in graph \a g (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor 
source(typename vigra::GridGraph<N, DirectedTag>::edge_descriptor const & e,
       vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.source(e);
}

    /** \brief Get a vertex descriptor for the end vertex of edge \a e in graph \a g (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor 
target(typename vigra::GridGraph<N, DirectedTag>::edge_descriptor const & e,
       vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.target(e);
}

    /** \brief Get a (begin, end) iterator pair for the edges of graph \a g (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
std::pair<typename vigra::GridGraph<N, DirectedTag>::edge_iterator, 
          typename vigra::GridGraph<N, DirectedTag>::edge_iterator>
edges(vigra::GridGraph<N, DirectedTag> const & g) 
{
    return std::make_pair(g.get_edge_iterator(), g.get_edge_end_iterator());
}

    /** \brief Return the number of edges in graph \a g (API: boost).
    */
template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::edges_size_type
num_edges(vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.num_edges();
}

// --------------------------------------------------
// support for AdjacencyMatrix concept:

// FIXME: test this
    /** \brief Return the pair (edge_descriptor, true) when an edge between vertices 
               \a u and \a v exists, or (lemon::INVALID, false) otherwise (API: boost).
    */
template<unsigned int N, class DirectedTag>
std::pair<typename vigra::GridGraph<N, DirectedTag>::edge_descriptor, bool>
edge(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & u,
     typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v,
     vigra::GridGraph<N, DirectedTag> const & g)
{
    return g.edge(u,v);
}

// provide get / put for MultiArrayViews, indexed by the 
// above-defined vertex_descriptor and edge_descriptor (in our case, a coordinate tuple):
// FIXME: place this into multi_array.hxx ?
    /** \brief Put value \a val at key \a k in the property map \a pmap (API: boost).
    */
template<unsigned int N, class T, class Stride, class U>
inline
void put(vigra::MultiArrayView<N, T, Stride> & pmap,
         typename vigra::MultiArrayView<N, T, Stride>::difference_type const & k,
         U const & val) 
{ 
    pmap[k] = val;  
}

    /** \brief Read the value at key \a k in property map \a pmap (API: boost).
    */
//template<unsigned int N, class T, class Stride>
//inline
//typename vigra::MultiArrayView<N, T, Stride>::const_reference 
//get(vigra::MultiArrayView<N, T, Stride> const & pmap,
//    typename vigra::MultiArrayView<N, T, Stride>::difference_type const & k)
//{ 
//    return pmap[k]; 
//}


//@}

}} // namespace vigra::boost_graph

namespace lemon {

template <typename GR>
class InDegMap;

    // LEMON-compatible property map for the in-degree of the nodes
template<unsigned int N, class DirectedTag>
class InDegMap<vigra::GridGraph<N, DirectedTag> >
: public vigra::GridGraph<N, DirectedTag>::InDegMap
{
  public:
    typedef vigra::GridGraph<N, DirectedTag> Graph;
    
    explicit InDegMap(const Graph& graph)
    : Graph::InDegMap(graph)
    {}
};

template <typename GR>
class OutDegMap;

    // LEMON-compatible property map for the out-degree of the nodes
template<unsigned int N, class DirectedTag>
class OutDegMap<vigra::GridGraph<N, DirectedTag> >
: public vigra::GridGraph<N, DirectedTag>::OutDegMap
{
  public:
    typedef vigra::GridGraph<N, DirectedTag> Graph;
    
    explicit OutDegMap(const Graph& graph)
    : Graph::OutDegMap(graph)
    {}
};


} // namespace lemon

namespace std {

template<unsigned int N, class DirectedTag>
ostream& operator<<(ostream& out,
                    typename vigra::GridGraph<N, DirectedTag>::vertex_iterator const & arg)
{
    out << "v" << arg.scanOrderIndex();
    return out;
}

template<unsigned int N, class DirectedTag>
ostream& operator<<(ostream& out,
                    typename vigra::GridGraph<N, DirectedTag>::adjacency_iterator const & arg)
{
    out << "nb" << arg.index();
    return out;
}

} // namespace std

#ifdef WITH_BOOST_GRAPH
namespace boost {
    using vigra::boost_graph::out_edges;
    using vigra::boost_graph::out_degree;
    using vigra::boost_graph::source;
    using vigra::boost_graph::target;
    using vigra::boost_graph::in_edges;
    using vigra::boost_graph::in_degree;
    using vigra::boost_graph::adjacent_vertices;
    using vigra::boost_graph::vertices;
    using vigra::boost_graph::edges;
    using vigra::boost_graph::edge;
    using vigra::boost_graph::num_vertices;
    using vigra::boost_graph::num_edges;
}
#endif /* WITH_BOOST_GRAPH */


#endif /* VIGRA_MULTI_GRIDGRAPH_HXX */



