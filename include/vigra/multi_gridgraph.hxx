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

#include "multi_iterator.hxx"
#include "multi_array.hxx"
#include "graphs.hxx"

template <unsigned int N>
struct NeighborhoodTests;

namespace vigra {

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

template <class Shape>
void
computeNeighborOffsetsOld(ArrayVector<Shape> const & neighborOffsets, 
                       ArrayVector<ArrayVector<bool> > const & neighborExists,
                       ArrayVector<ArrayVector<Shape> > & incrementOffsets,
                       ArrayVector<ArrayVector<GridGraphArcDescriptor<Shape::static_size> > > & edgeDescriptorOffsets,
                       ArrayVector<ArrayVector<MultiArrayIndex> > & indices,
                       bool directed, bool includeBackEdges, bool includeForwardEdges)
{
    typedef GridGraphArcDescriptor<Shape::static_size> EdgeDescriptor;
    
    unsigned int borderTypeCount = neighborExists.size();
    incrementOffsets.resize(borderTypeCount);
    edgeDescriptorOffsets.resize(borderTypeCount);
    indices.resize(borderTypeCount);
    
    for(unsigned int k=0; k<borderTypeCount; ++k)
    {
        incrementOffsets[k].clear();
        edgeDescriptorOffsets[k].clear();
        indices[k].clear();
        
        unsigned int j   = includeBackEdges
                              ? 0
                              : neighborOffsets.size() / 2,
                     end = includeForwardEdges
                              ? neighborOffsets.size()
                              : neighborOffsets.size() / 2;
        for(; j < end; ++j)
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
                
                if(directed || j < neighborOffsets.size() / 2) // directed edge or backward edge
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
            }
        }
    }
}

} // namespace detail

template<unsigned int N, bool BackEdgesOnly=false>
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
class GridGraphEdgeIterator;

template<unsigned int N, bool BackEdgesOnly=false>
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
        unsigned int nbtype = g.get_border_type(v);
        init(&(*g.edgeIncrementArray())[nbtype], &(*g.neighborIndexArray(BackEdgesOnly))[nbtype], *v, opposite);
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
        unsigned int nbtype = g.get_border_type(v);
        init(&(*g.edgeIncrementArray())[nbtype], &(*g.neighborIndexArray(BackEdgesOnly))[nbtype], v, opposite);
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

template<unsigned int N, bool BackEdgesOnly=false>
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

template<unsigned int N, bool BackEdgesOnly=false>
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

using boost::directed_tag;
using boost::undirected_tag;

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
        typedef boost::read_write_property_map_tag         category;
        
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
        typedef boost::read_write_property_map_tag         category;
        
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

    // Grid Graph class to adapt vigra MultiArrayViews to a BGL-like interface.
    //       This class only knows about
    //       - dimensions
    //       - shape
    //       - neighborhood type (DirectedNeighborhood or IndirectNeighborhood, i.e.\ including diagonal neighbors)
    //       - whether the graph is directed or undirected
template<unsigned int N, class DirectedTag>
class GridGraph
: public detail::GridGraphBase<N, DirectedTag>
{
public:
    static const bool is_directed = IsSameType<DirectedTag, directed_tag>::value;
    
    typedef detail::GridGraphBase<N, DirectedTag>   base_type;
    typedef GridGraph<N, DirectedTag>               self_type;
    typedef typename MultiArrayShape<N>::type       shape_type;
    typedef typename MultiArrayShape<N+1>::type     edge_propmap_shape_type;
    typedef MultiArrayIndex                         index_type;
    typedef MultiArrayIndex                         vertices_size_type;
    typedef MultiArrayIndex                         edges_size_type;
    typedef MultiArrayIndex                         degree_size_type;

    // Boost Graph interface
    typedef MultiCoordinateIterator<N>              vertex_iterator;
    typedef GridGraphNeighborIterator<N>            neighbor_vertex_iterator;
    typedef GridGraphNeighborIterator<N, true>      back_neighbor_vertex_iterator;
    typedef neighbor_vertex_iterator                adjacency_iterator;
    typedef GridGraphInArcIterator<N>               in_edge_iterator;
    typedef GridGraphOutArcIterator<N>              out_edge_iterator;
    typedef GridGraphOutArcIterator<N, true>        out_back_edge_iterator;
    typedef GridGraphArcIterator<N, !is_directed>   edge_iterator;

    typedef shape_type                              vertex_descriptor;
    typedef GridGraphArcDescriptor<N>               edge_descriptor;
 
    typedef DirectedTag                             directed_category;
    typedef boost::disallow_parallel_edge_tag       edge_parallel_category;
    typedef boost::no_property                      vertex_property_type; // we only support "external properties".
    // FIXME: Maybe support the vertex -> coordinate map (identity) as the only internal property map
    // and additionally the vertex_descriptor -> ID map (vertex_index = SOI).

    struct traversal_category 
    : virtual public boost::bidirectional_graph_tag,
      virtual public boost::adjacency_graph_tag,
      virtual public boost::vertex_list_graph_tag,
      virtual public boost::edge_list_graph_tag,
      virtual public boost::adjacency_matrix_tag
    {};
    
    typedef ArrayVector<shape_type>                      NeighborOffsetArray;
    typedef ArrayVector<NeighborOffsetArray>             RelativeNeighborOffsetsArray;
    typedef ArrayVector<ArrayVector<edge_descriptor> >   RelativeEdgeOffsetsArray;
    typedef ArrayVector<ArrayVector<MultiArrayIndex> >   IndexArray;

    // LEMON interface
    typedef self_type                               Graph;
    typedef vertex_descriptor                       Node;
    typedef vertex_iterator                         NodeIt;
    
    typedef GridGraphArcDescriptor<N>               Arc;
    typedef GridGraphOutArcIterator<N>              OutArcIt;
    typedef GridGraphOutArcIterator<N, true>        OutBackArcIt;
    typedef GridGraphArcIterator<N, false>          ArcIt;
    typedef GridGraphInArcIterator<N>               InArcIt;
    
    typedef typename MultiArrayShape<N+1>::type     Edge;
    typedef GridGraphOutEdgeIterator<N>             IncEdgeIt;
    typedef GridGraphOutEdgeIterator<N, true>       IncBackEdgeIt;
    typedef GridGraphEdgeIterator<N, !is_directed>  EdgeIt;
    
    typedef lemon::True NodeNumTag;
    typedef lemon::True EdgeNumTag;
    typedef lemon::True ArcNumTag;
    typedef lemon::True FindEdgeTag;
    typedef lemon::True FindArcTag;

    class IndexMap 
    {
      public:
        typedef Node                                    Key;
        typedef Node                                    Value;
        typedef Key                                     key_type;
        typedef Value                                   value_type; 
        typedef Value const &                           reference;
        typedef boost::readable_property_map_tag        category;

        IndexMap()
        {}

        IndexMap(const GridGraph&)
        {}

        Value const & operator[](Key const & key) const 
        {
            return key;
        }
    };
        
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
        typedef boost::read_write_property_map_tag         category;
        
        typedef lemon::True                                ReferenceMapTag;
        typedef key_type                                   Key;
        typedef value_type                                 Value;
        typedef reference                                  Reference;
        typedef const_reference                            ConstReference;

        NodeMap()
        : base_type()
        {}
        
        explicit NodeMap(GridGraph const & g)
        : base_type(g.shape())
        {}
        
        NodeMap(GridGraph const & g, T const & t)
        : base_type(g.shape(), t)
        {}

       explicit  NodeMap(difference_type const & shape)
        : base_type(shape)
        {}
        
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
        
        void set(Key const & k, Value const & v)
        {
            (*this)[k] = v;
        }
    };
    
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
        typedef boost::read_write_property_map_tag         category;
        
        typedef lemon::True                                ReferenceMapTag;
        typedef key_type                                   Key;
        typedef value_type                                 Value;
        typedef reference                                  Reference;
        typedef const_reference                            ConstReference;

        EdgeMap()
        : base_type()
        {}
        
        explicit EdgeMap(GridGraph const & g)
        : base_type(g.edge_propmap_shape())
        {}
        
        EdgeMap(GridGraph const & g, T const & t)
        : base_type(g.edge_propmap_shape(), t)
        {}

        explicit EdgeMap(difference_type const & shape)
        : base_type(shape)
        {}
        
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
        
        void set(Key const & k, Value const & v)
        {
            (*this)[k] = v;
        }
    };

        // dummy default constructor to satisfy adjacency_graph concept
    GridGraph()
    {}

        //! Constructor for grid graph. 
        //  @param shape                  an array of the graph's dimensions as a TinyVector
        //  @param ntype                  DirectNeighborhood for direct neighborhood (axis-aligned edges only) 
        //                                or IndirectNeighborhood for indirect neighborhood (including all diagonal edges)
    GridGraph(shape_type const &shape, NeighborhoodType ntype = DirectNeighborhood) 
    : shape_(shape),
      num_vertices_(prod(shape)),
      num_edges_(gridGraphEdgeCount(shape, ntype, is_directed)), 
      neighborhoodType_(ntype)
    {
        ArrayVector<ArrayVector<bool> > neighborExists;
        
        // populate the neighborhood tables:
        // FIXME: this might be static (but make sure that it works with multi-threading)
        detail::makeArrayNeighborhood(neighborOffsets_, neighborExists, neighborhoodType_);
        detail::computeNeighborOffsets(neighborOffsets_, neighborExists, incrementalOffsets_, 
                                       edgeDescriptorOffsets_, neighborIndices_, backIndices_, is_directed);
        
        // compute the neighbor offsets per neighborhood type
        // detail::makeArraySubNeighborhood(neighborhood[0], neighborExists, shape_type(1), neighborhoodIndices);
    }
    
        // convention: Node id equals the scan order index in an EdgeMap
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
    
    Node nodeFromId(index_type i) const
    {
        Node res(SkipInitialization);
        detail::ScanOrderToCoordinate<N>::exec(i, shape(), res);
        return res;
    }

    index_type maxNodeId() const
    {
        return prod(shape()) - 1;
    }
    
    Node const & pos(Node const & v) const
    {
        return v;
    }
    
    vertex_iterator get_vertex_iterator() const 
    {
        return vertex_iterator(shape_);
    }

    vertex_iterator get_vertex_iterator(vertex_descriptor const & v) const 
    {
        return vertex_iterator(shape_) + v;
    }

    vertex_iterator get_vertex_end_iterator() const 
    {
        return get_vertex_iterator().getEndIterator();
    }

    neighbor_vertex_iterator get_neighbor_vertex_iterator(vertex_descriptor const & v) const 
    {
        return neighbor_vertex_iterator(*this, v);
    }

    neighbor_vertex_iterator get_neighbor_vertex_iterator(vertex_iterator const & v) const 
    {
        return neighbor_vertex_iterator(*this, v);
    }

    neighbor_vertex_iterator get_neighbor_vertex_end_iterator(vertex_descriptor const & v) const 
    {
       return get_neighbor_vertex_iterator(v).getEndIterator();
    }

    neighbor_vertex_iterator get_neighbor_vertex_end_iterator(vertex_iterator const & v) const 
    {
       return get_neighbor_vertex_iterator(v).getEndIterator();
    }

    back_neighbor_vertex_iterator get_back_neighbor_vertex_iterator(vertex_descriptor const & v) const 
    {
        return back_neighbor_vertex_iterator(*this, v);
    }

    back_neighbor_vertex_iterator get_back_neighbor_vertex_iterator(vertex_iterator const & v) const 
    {
        return back_neighbor_vertex_iterator(*this, v);
    }

    back_neighbor_vertex_iterator get_back_neighbor_vertex_end_iterator(vertex_descriptor const & v) const 
    {
       return get_back_neighbor_vertex_iterator(v).getEndIterator();
    }

    back_neighbor_vertex_iterator get_back_neighbor_vertex_end_iterator(vertex_iterator const & v) const 
    {
       return get_back_neighbor_vertex_iterator(v).getEndIterator();
    }

    // --------------------------------------------------
    // support for VertexListGraph:

    vertices_size_type num_vertices() const 
    {
        return num_vertices_;
    }

    vertices_size_type nodeNum() const 
    {
        return num_vertices();
    }

    // --------------------------------------------------
    // support for IncidenceGraph:

        // convention: Edge id equals the scan order index in an EdgeMap
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

    Edge edgeFromId(index_type i) const
    {
        Edge res(SkipInitialization);
        detail::ScanOrderToCoordinate<N+1>::exec(i, edge_propmap_shape(), res);
        return res;
    }
    
    index_type maxEdgeId() const
    {
        if(is_directed)
            return maxArcId();
        if(edgeNum() == 0)
            return -1;
        Node lastNode = shape() - shape_type(1);
        Arc a(lastNode, backIndices_[get_border_type(lastNode)].back(), false);
        return detail::CoordinateToScanOrder<N+1>::exec(edge_propmap_shape(), a);
    }
    
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
    
    Arc arcFromId(index_type i) const
    {
        Arc res;
        detail::ScanOrderToCoordinate<N+1>::exec(i, arc_propmap_shape(), res);
        return undirectedArc(res);
    }
    
    index_type maxArcId() const
    {
        if(edgeNum() == 0)
            return -1;
        Node lastNode = shape() - shape_type(1);
        index_type n = neighborIndices_[get_border_type(lastNode)][0];
        Arc a(neighbor(lastNode, n), oppositeIndex(n), false);
        return detail::CoordinateToScanOrder<N+1>::exec(arc_propmap_shape(), a);
    }
    
    bool direction(Arc const & a) const
    {
        return !a.isReversed();
    }
    
    Arc direct(Edge const & e, bool forward) const
    {
        if(!is_directed || forward)
            return Arc(e, !forward);
        else
            return Arc(v(e), oppositeIndex(e[N]), true);
    }
    
    Arc direct(Edge const & e, Node const & n) const
    {
        if(u(e) == n)
            return direct(e, true);
        if(v(e) == n)
            return direct(e, false);
        return Arc(lemon::INVALID);
    }
    
    Node oppositeNode(Node const & n, Edge const & e) const
    {
        Node start(u(e)), end(v(e));
        if(n == start)
            return end;
        if(n == end)
            return start;
        return Node(lemon::INVALID);
    }
    
    Arc oppositeArc(Arc const & a) const
    {
        return is_directed
                 ? Arc(neighbor(a.vertexDescriptor(), a.edgeIndex()), oppositeIndex(a.edgeIndex()), false)
                 : Arc(a, !a.isReversed());
    }
    
    Arc directedArc(Arc const & a) const
    {
        return a.isReversed()
                 ? Arc(neighbor(a.vertexDescriptor(), a.edgeIndex()), oppositeIndex(a.edgeIndex()), false)
                 : a;
    }
    
    Arc undirectedArc(Arc const & a) const
    {
        return a.edgeIndex() < maxUniqueDegree() 
                 ? a
                 : Arc(neighbor(a.vertexDescriptor(), a.edgeIndex()), oppositeIndex(a.edgeIndex()), true);
    }
    
    Node baseNode(IncEdgeIt const & e)  const
    {
        return source(e.arcDescriptor());
    }
    
    Node runningNode(IncEdgeIt const & e)  const
    {
        return target(e.arcDescriptor());
    }
    
    Node baseNode(IncBackEdgeIt const & e)  const
    {
        return source(e.arcDescriptor());
    }
    
    Node runningNode(IncBackEdgeIt const & e)  const
    {
        return target(e.arcDescriptor());
    }
    
    Node baseNode(OutArcIt const & a)  const
    {
        return source(*a);
    }
    
    Node runningNode(OutArcIt const & a)  const
    {
        return target(*a);
    }
    
    Node baseNode(OutBackArcIt const & a)  const
    {
        return source(*a);
    }
    
    Node runningNode(OutBackArcIt const & a)  const
    {
        return target(*a);
    }
    
    vertex_descriptor source(Arc const & e) const 
    {
        return source_or_target(e, true);
    }

    vertex_descriptor target(Arc const & e) const 
    {
        return source_or_target(e, false);
    }

    vertex_descriptor u(Edge const & e) const 
    {
        return vertex_descriptor(e.template subarray<0,N>());
    }

    vertex_descriptor v(Edge const & e) const 
    {
        return vertex_descriptor(e.template subarray<0,N>()) + neighborOffsets_[e[N]];
    }

    out_edge_iterator get_out_edge_iterator(vertex_descriptor const & v) const 
    {
        return out_edge_iterator(*this, v);
    }

    out_edge_iterator get_out_edge_iterator(vertex_iterator const & v) const 
    {
        return out_edge_iterator(*this, v);
    }

    out_edge_iterator get_out_edge_end_iterator(vertex_descriptor const & v) const 
    {
        return get_out_edge_iterator(v).getEndIterator();
    }

    out_edge_iterator get_out_edge_end_iterator(vertex_iterator const & v) const 
    {
        return get_out_edge_iterator(v).getEndIterator();
    }

    out_back_edge_iterator get_out_back_edge_iterator(vertex_descriptor const & v) const 
    {
        return out_back_edge_iterator(*this, v);
    }

    out_back_edge_iterator get_out_back_edge_iterator(vertex_iterator const & v) const 
    {
        return out_back_edge_iterator(*this, v);
    }

    out_back_edge_iterator get_out_back_edge_end_iterator(vertex_descriptor const & v) const 
    {
        return get_out_back_edge_iterator(v).getEndIterator();
    }

    out_back_edge_iterator get_out_back_edge_end_iterator(vertex_iterator const & v) const 
    {
        return get_out_back_edge_iterator(v).getEndIterator();
    }

    in_edge_iterator get_in_edge_iterator(vertex_descriptor const & v) const 
    {
        return in_edge_iterator(*this, v);
    }

    in_edge_iterator get_in_edge_iterator(vertex_iterator const & v) const 
    {
        return in_edge_iterator(*this, v);
    }

    in_edge_iterator get_in_edge_end_iterator(vertex_descriptor const & v) const 
    {
        return get_in_edge_iterator(v).getEndIterator();
    }

    in_edge_iterator get_in_edge_end_iterator(vertex_iterator const & v) const 
    {
        return get_in_edge_iterator(v).getEndIterator();
    }

    degree_size_type out_degree(vertex_iterator const & v) const 
    {
        return (degree_size_type)neighborIndices_[get_border_type(v)].size();
    }

    degree_size_type out_degree(vertex_descriptor const & v) const 
    {
        return (degree_size_type)neighborIndices_[get_border_type(v)].size();
    }
    
    degree_size_type back_degree(vertex_iterator const & v) const 
    {
        return (degree_size_type)backIndices_[get_border_type(v)].size();
    }

    degree_size_type back_degree(vertex_descriptor const & v) const 
    {
        return (degree_size_type)backIndices_[get_border_type(v)].size();
    }
    
    degree_size_type forward_degree(vertex_iterator const & v) const 
    {
        unsigned int bt = get_border_type(v);
        return (degree_size_type)(neighborIndices_[bt].size() - backIndices_[bt].size());
    }

    degree_size_type forward_degree(vertex_descriptor const & v) const 
    {
        unsigned int bt = get_border_type(v);
        return (degree_size_type)(neighborIndices_[bt].size() - backIndices_[bt].size());
    }
    
    degree_size_type in_degree(vertex_iterator const & v) const 
    {
        return out_degree(v);
    }

    degree_size_type in_degree(vertex_descriptor const & v) const 
    {
        return out_degree(v);
    }
    
    degree_size_type degree(vertex_iterator const & v) const 
    {
        return degree(*v);
    }

    degree_size_type degree(vertex_descriptor const & v) const 
    {
        return is_directed
                   ? 2*out_degree(v)
                   : out_degree(v);
    }
    
    // --------------------------------------------------
    // support for EdgeListGraph:

    edges_size_type num_edges() const 
    {
        return num_edges_;
    }

    edges_size_type edgeNum() const 
    {
        return num_edges();
    }

    edges_size_type arcNum() const 
    {
        return is_directed
                   ? num_edges()
                   : 2*num_edges();
    }
    
    edge_iterator get_edge_iterator() const 
    {
        return edge_iterator(*this);
    }

    edge_iterator get_edge_end_iterator() const 
    {
        return get_edge_iterator().getEndIterator();
    }

    // --------------------------------------------------
    // support for AdjacencyMatrix concept:

    std::pair<edge_descriptor, bool>
    edge(vertex_descriptor const & u, vertex_descriptor const & v) const
    {
        std::pair<edge_descriptor, bool> res;
        res.second = false;

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

    Edge findEdge(Node const & u, Node const & v, Edge const & = lemon::INVALID) const 
    {
        std::pair<edge_descriptor, bool> res(edge(u, v));
        return res.second
                 ? res.first
                 : Edge(lemon::INVALID);
    }
    
    Arc findArc(Node const & u, Node const & v, Arc const & = lemon::INVALID) const 
    {
        std::pair<edge_descriptor, bool> res(edge(u, v));
        return res.second
                 ? res.first
                 : Arc(lemon::INVALID);
    }
    
    // --------------------------------------------------
    // other helper functions:
    
    IndexMap indexMap() const 
    {
        return IndexMap();
    }

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
    
    // //! In case of an undirected graph, the edge u->v is the same as v->u.
    // //  This function folds in the edge descriptors corresponding to "causal neighbors",
    // //  i.e. those to vertices with lower scan order index, by reversing the edge and computing
    // //  the corresponding edge descriptor.
    // //  (assuming here the neighbor-indices are in the order of noncausal, causal neighbors....
    // //   FIXME: check this again in neighborhood construction!)
    // //  (At least, the neighborhood construction is symmetrical, hence this should be OK)
    // inline 
    // edge_descriptor
    // map_to_undirected_edge(const edge_descriptor &e) const {
        // edge_descriptor res = e;
        // //assert(!res.isReversed()); 
        // if (res[N] >= halfMaxDegree()) {
            // TinyVectorView<typename edge_descriptor::value_type, N> vertex(res.data());
            // vertex += neighborCoordOffset(res[N]);
            // res[N] = maxDegree() - res[N] - 1;
            // res.setReversed(!res.isReversed());
        // }
        // return res;
    // }



  protected:
    NeighborOffsetArray neighborOffsets_;
    IndexArray neighborIndices_, backIndices_;
    RelativeNeighborOffsetsArray incrementalOffsets_;
    RelativeEdgeOffsetsArray edgeDescriptorOffsets_;
    shape_type shape_;
    MultiArrayIndex num_vertices_, num_edges_;
    NeighborhoodType neighborhoodType_;
};

} // namespace vigra

namespace boost {

template <unsigned int N, class T, class Acc>
struct property_traits<vigra::MultiArray<N, T, Acc> >
{
    typedef vigra::MultiArray<N, T, Acc>             type;
    typedef typename type::key_type                  key_type;
    typedef typename type::value_type                value_type; 
    typedef typename type::reference                 reference;
    typedef boost::read_write_property_map_tag       category;
};

template <unsigned int N, class T, class Acc>
struct property_traits<vigra::MultiArray<N, T, Acc> const>
{
    typedef vigra::MultiArray<N, T, Acc> const       type;
    typedef typename type::key_type                  key_type;
    typedef typename type::value_type                value_type; 
    typedef typename type::const_reference           reference;
    typedef boost::readable_property_map_tag         category;
};

template <unsigned int N, class T, class Stride>
struct property_traits<vigra::MultiArrayView<N, T, Stride> >
{
    typedef vigra::MultiArrayView<N, T, Stride>       type;
    typedef typename type::key_type                   key_type;
    typedef typename type::value_type                 value_type; 
    typedef typename type::reference                  reference;
    typedef boost::read_write_property_map_tag        category;
};

template <unsigned int N, class T, class Stride>
struct property_traits<vigra::MultiArrayView<N, T, Stride> const>
{
    typedef vigra::MultiArrayView<N, T, Stride> const     type;
    typedef typename type::key_type                       key_type;
    typedef typename type::value_type                     value_type; 
    typedef typename type::const_reference                reference;
    typedef boost::readable_property_map_tag              category;
};



template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::degree_size_type
out_degree(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v, 
           vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.out_degree(v);
}

template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::degree_size_type
in_degree(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v, 
          vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.in_degree(v);
}

template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::degree_size_type
degree(typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor const & v, 
       vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.degree(v);
}

template<unsigned int N, class DirectedTag>
inline
std::pair<typename vigra::GridGraph<N, DirectedTag>::vertex_iterator, 
          typename vigra::GridGraph<N, DirectedTag>::vertex_iterator>
vertices(vigra::GridGraph<N, DirectedTag> const & g) 
{
    return std::make_pair(g.get_vertex_iterator(),
                          g.get_vertex_end_iterator());    
}

template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::vertices_size_type
num_vertices(vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.num_vertices();
}

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

template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor 
source(typename vigra::GridGraph<N, DirectedTag>::edge_descriptor const & e,
       vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.source(e);
}

template<unsigned int N, class DirectedTag>
inline
typename vigra::GridGraph<N, DirectedTag>::vertex_descriptor 
target(typename vigra::GridGraph<N, DirectedTag>::edge_descriptor const & e,
       vigra::GridGraph<N, DirectedTag> const & g) 
{
    return g.target(e);
}

template<unsigned int N, class DirectedTag>
inline
std::pair<typename vigra::GridGraph<N, DirectedTag>::edge_iterator, 
          typename vigra::GridGraph<N, DirectedTag>::edge_iterator>
edges(vigra::GridGraph<N, DirectedTag> const & g) 
{
    return std::make_pair(g.get_edge_iterator(), g.get_edge_end_iterator());
}


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
template<unsigned int N, class T, class Stride, class U>
inline
void put(vigra::MultiArrayView<N, T, Stride> & pmap,
         typename vigra::MultiArrayView<N, T, Stride>::difference_type const & k,
         U const & val) 
{ 
    pmap[k] = val;  
}

template<unsigned int N, class T, class Stride>
inline
typename vigra::MultiArrayView<N, T, Stride>::const_reference 
get(vigra::MultiArrayView<N, T, Stride> const & pmap,
    typename vigra::MultiArrayView<N, T, Stride>::difference_type const & k)
{ 
    return pmap[k]; 
}

#if 0

// provide get / put for MultiArrayViews, indexed by the above-defined vertex_descriptor (in this case, a coordinate tuple):

template<class VIEW>
class MultiArrayView_property_map 
{
public:
    // typedef boost::read_write_property_map_tag category;
    typedef boost::lvalue_property_map_tag category;
    //    typedef int value_type;
    typedef typename VIEW::value_type value_type;
    typedef typename VIEW::reference reference;
    typedef typename VIEW::const_reference const_reference;
    typedef typename vigra::GridGraph<VIEW::actual_dimension> graph_type;
    typedef typename graph_type::vertex_descriptor key_type;
    MultiArrayView_property_map(const VIEW& view)
        : view_(view) { }
    template <class T2>
    inline
    reference operator[](const T2 & x) { return view_[x]; } 
    template <class T2>
    inline
    const_reference operator[](const T2 & x) const { return view_[x]; } 

protected:
    VIEW view_;
};

template<class VIEW>
inline
void put(MultiArrayView_property_map<VIEW> &pmap,
         const typename MultiArrayView_property_map<VIEW>::key_type &k,
         const typename MultiArrayView_property_map<VIEW>::value_type& val
         ) 
{ 
    pmap[k] = val;  
}

template<class VIEW>
inline
typename MultiArrayView_property_map<VIEW>::const_reference 
get(
    const MultiArrayView_property_map<VIEW> & pmap, 
    const typename MultiArrayView_property_map<VIEW>::key_type &k)
{ 
    return pmap[k]; 
}


//! use a MultiArrayView as an undirected edge property map
//  (edge_descriptor keys will be transformed by "flipping" the edges if necessary)
template<class VIEW, class GRAPH>
class MultiArrayView_undirected_edge_property_map 
{
public:
    // typedef boost::read_write_property_map_tag category;
    typedef boost::lvalue_property_map_tag category;
    //    typedef int value_type;
    typedef typename VIEW::value_type value_type;
    typedef typename VIEW::reference reference;
    typedef typename VIEW::const_reference const_reference;
    typedef GRAPH graph_type;
    typedef typename graph_type::edge_descriptor key_type;
    MultiArrayView_undirected_edge_property_map(const VIEW& view, const GRAPH& graph)
        : view_(view), graph_(graph) { }
    template <class T2>
    inline
    reference operator[](const T2 & x) { return view_[graph_.map_to_undirected_edge(x)]; } 
    template <class T2>
    inline
    const_reference operator[](const T2 & x) const { return view_[graph_.map_to_undirected_edge(x)]; } 

protected:
    VIEW view_;
    const GRAPH &graph_;
};

template<class VIEW, class GRAPH>
inline
void put(MultiArrayView_undirected_edge_property_map<VIEW, GRAPH> &pmap,
         const typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::key_type &k,
         const typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::value_type& val
         ) 
{ 
    pmap[k] = val;  
}

template<class VIEW, class GRAPH>
inline
typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::const_reference 
get(
    const MultiArrayView_undirected_edge_property_map<VIEW, GRAPH> & pmap, 
    const typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::key_type &k)
{ 
    return pmap[k]; 
}





// property map support for mapping coordinates to scan-order indices:

template<unsigned int N>
struct IDMapper {
    typedef typename vigra::GridGraph<N> graph_type;
    typedef boost::readable_property_map_tag category;
    typedef typename graph_type::index_type value_type;
    typedef typename graph_type::vertex_descriptor key_type;
    typedef const value_type& reference;


    IDMapper(const graph_type &graph) 
        : map_helper(graph.get_vertex_iterator())
    {}

    typename graph_type::vertex_iterator map_helper;
};

template<unsigned int N>
struct property_map<vigra::GridGraph<N>, boost::vertex_index_t>
{
    typedef IDMapper<N> type;
    typedef IDMapper<N> const_type;
};


template<unsigned int N>
inline
typename IDMapper<N>::value_type
get(const IDMapper<N> & mapper, 
    const typename IDMapper<N>::key_type &k)
{ 
    return (mapper.map_helper + k).scanOrderIndex();
}


template<unsigned int N>
typename boost::property_map<vigra::GridGraph<N>, boost::vertex_index_t>::type
//typename IDMapper<N>
get(boost::vertex_index_t, const vigra::GridGraph<N> &graph) {
    // return a lightweight wrapper for the CoupledIterator, which easily allows the conversion of 
    // coordinates via its += operator followed by index().
    return IDMapper<N>(graph);
}

#if 0 
// CHECK if required: also provide the direct (three-parameter) version for index lookup
template<unsigned int N>
typename vigra::GridGraph<N>::vertices_size_type
get(boost::vertex_index_t, 
    const vigra::GridGraph<N> &graph,
    const typename vigra::GridGraph<N>::vertex_descriptor &v) {
    return (IDMapper<N>(graph).map_helper + v).scanOrderIndex();
}
#endif



// TODO:
// eventually provide an edge_index property map as well?
// (edge_descriptor -> linear contiguous edge index)

#endif

} // namespace boost


#if 0

//#define VERBOSE

#include "graphs.hxx"
#include "multi_iterator.hxx"
#include "multi_iterator_coupled.hxx"
#include "mathutil.hxx"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iterator>


namespace vigra {
    // forward declare graph type:
    template<unsigned int N>
    class GridGraph;
}

#include "multi_shape.hxx"
#include "multi_iterator.hxx"
#include "multi_gridgraph_coords_neighbor_iterator.hxx"
#include "multi_gridgraph_coords_edge_descriptor.hxx"
#include "multi_gridgraph_coords_out_edge_iterator.hxx"
#include "multi_gridgraph_coords_edge_iterator_decl.hxx"


namespace vigra {

namespace detail {

template <class Shape>
void
makeArrayNeighborhood2(ArrayVector<ArrayVector<Shape> > & neighborOffsets, 
                      ArrayVector<ArrayVector<bool> > & neighborExists,
                      ArrayVector<ArrayVector<bool> > & causalNeighborExists,
                      ArrayVector<ArrayVector<bool> > & anticausalNeighborExists,
                      ArrayVector<ArrayVector<int> > & neighborIndexLookup,
                      NeighborhoodType neighborhoodType = DirectNeighborhood)
{
    enum { N = Shape::static_size };
    unsigned int size = 1 << 2*N;
    Shape strides = cumprod(Shape(MultiArrayIndex(3))) / 3; 
    
    neighborOffsets.resize(size);
    neighborOffsets[0].clear(); // [0] is the standard case of all neighbors present
    if(neighborhoodType == DirectNeighborhood)
    {
        MakeDirectArrayNeighborhood<N-1>::offsets(neighborOffsets[0]);
    }
    else
    {
        Shape point; // represents the center
        MakeIndirectArrayNeighborhood<N-1>::offsets(neighborOffsets[0], point);
    }
    
    unsigned int neighborCount = neighborOffsets[0].size(); // maximal number of neighbors

#ifdef VERBOSE    
    std::cerr << " size " << neighborCount << ": " << neighborOffsets[0] << "\n strides ";
    for(unsigned int l=0; l<neighborCount; ++l)
        std::cerr << dot(neighborOffsets[0][l], strides) << ", ";
    std::cerr << "\n\n";
#endif
    
    neighborExists.resize(size);
    causalNeighborExists.resize(size);
    anticausalNeighborExists.resize(size);
    neighborIndexLookup.resize(size);

    for(unsigned int k=0; k<size; ++k) // iterate all k neighborhood codes
    {
        if (k>0) 
            neighborOffsets[k].clear();
        neighborExists[k].clear();
        if(neighborhoodType == DirectNeighborhood)
        {
            MakeDirectArrayNeighborhood<N-1>::exists(neighborExists[k], k);
        }
        else
        {
            MakeIndirectArrayNeighborhood<N-1>::exists(neighborExists[k], k);
        }
        
        causalNeighborExists[k].resize(neighborCount);
        anticausalNeighborExists[k].resize(neighborCount);
        
        for(unsigned int l = 0; l<neighborCount; ++l)
        {
            MultiArrayIndex stride = dot(neighborOffsets[0][l], strides);
            if(stride < 0)
            {
                causalNeighborExists[k][l] = neighborExists[k][l];
                anticausalNeighborExists[k][l] = false;
            }
            else
            {
                causalNeighborExists[k][l] = false;
                anticausalNeighborExists[k][l] = neighborExists[k][l];
            }
            if (neighborExists[k][l])
                neighborIndexLookup[k].push_back(l);
            if (k>0)
                if (neighborExists[k][l])
                    neighborOffsets[k].push_back(neighborOffsets[0][l]);
        }
    }

}

template <class Shape>
void
makeArraySubNeighborhood2(const ArrayVector<Shape> & allNeighborOffsets, 
             const ArrayVector<ArrayVector<bool> > & neighborExists,
             const Shape strides,
             ArrayVector<ArrayVector<MultiArrayIndex> > & neighborIndices
             )
{
    enum { N = Shape::static_size };
    unsigned int size = 1 << 2*N;
    
    neighborIndices.resize(size);
    const unsigned int neighborCount = allNeighborOffsets.size(); // maximal number of neighbors
    
    for (unsigned int k=0; k<size; ++k)  // iterate all k neighborhood codes
    for(unsigned int l=0; l<neighborCount; ++l) 
        if (neighborExists[k][l])
        neighborIndices[k].push_back(dot(allNeighborOffsets[l], strides));
#if 0
    for (unsigned int k=0; k<size; ++k)  // iterate all k neighborhood codes
    {
    std::cerr << " NB-type " << k << ": ";
    for(unsigned int l=0; l<neighborCount; ++l) 
        if (neighborExists[k][l])
        {
        std::cerr << neighborIndices[k].back() << ", ";
        }
    std::cerr << std::endl;
    }
#endif
}

} // namespace detail


//! Grid Graph class to adapt vigra MultiArrayViews to a BGL-like interface.
//       This class only knows about
//       - dimensions
//       - shape
//       - neighborhood type (DirectedNeighborhood or IndirectNeighborhood, i.e.\ including diagonal neighbors)
template<unsigned int N>
class GridGraph
{
public:
    typedef GridGraph<N> self_type;

    typedef typename MultiArrayShape<N>::type shape_type;
    typedef typename MultiArrayShape<N+1>::type edge_propmap_shape_type;
    typedef MultiArrayIndex index_type;

    typedef detail::MultiCoordinateIterator<N> vertex_iterator;
    typedef detail::CoordsGridGraphNeighborIterator<N> neighbor_vertex_iterator;
    typedef detail::CoordsGridGraphOutEdgeIterator<N> out_edge_iterator;
    typedef detail::GridGraphEdgeIterator<self_type> edge_iterator;


    struct traversal_category : virtual public boost::incidence_graph_tag,
                                virtual public boost::adjacency_graph_tag,
                                virtual public boost::vertex_list_graph_tag,
                                virtual public boost::edge_list_graph_tag,
                                virtual public boost::adjacency_matrix_tag
                                { };

    typedef shape_type  vertex_descriptor;
    //typedef typename MultiArrayShape<N+1>::type edge_descriptor;
    typedef detail::GridGraphArcDescriptor<N> edge_descriptor;

    typedef void in_edge_iterator; // for bidirectional_graph concept, not implemented here

    typedef neighbor_vertex_iterator      adjacency_iterator; // must be a MultiPassInputIterator model

    typedef undirected_tag directed_category;
    typedef boost::disallow_parallel_edge_tag edge_parallel_category;

    typedef MultiArrayIndex     vertices_size_type;
    typedef MultiArrayIndex     edges_size_type;
    typedef MultiArrayIndex     degree_size_type;

    // we only support "external properties".
    typedef boost::no_property vertex_property_type;
    // TODO: Maybe support the vertex -> coordinate map (identity) as the only internal property map
    // and additionally the vertex_descriptor -> ID map (vertex_index = SOI).


    // dummy default constructor to satisfy adjacency_graph concept
    GridGraph()
    {}


    //! Constructor for grid graph. 
    //  @param shape                  an array of the graph's dimensions as a TinyVector
    //  @param directNeighborsOnly    true for direct neighborhood (axis-aligned edges only) 
    //                                or false for indirect neighborhood (including all diagonal edges)
    GridGraph(shape_type const &shape, borderType directNeighborsOnly = IndirectNeighborhood) 
        : shape_(shape)
    {
        // use makeArrayNeighborhood to populate neighborhood tables:
        detail::makeArrayNeighborhood(neighborhood, 
                                      neighborExists, 
                                      causalNeighborhood, 
                                      anticausalNeighborhood, 
                                      neighborIndexLookup, 
                                      directNeighborsOnly);
        
        // compute the neighbor offsets per neighborhood type
        detail::makeArraySubNeighborhood(neighborhood[0], neighborExists, shape_type(1), neighborhoodIndices);

        // compute total number of edges
        num_edges_ = 0;
        for (unsigned int i=0; i<neighborhood[0].size(); ++i) {
            size_t product = 1;
            for (unsigned int j=0; j<N; ++j) {
                product *= (neighborhood[0][i][j]==0) ? shape_[j] : shape_[j]-1;
            }
            num_edges_ += product;
        }
        num_edges_ /= 2; // because of undirectedness
    }

    inline
    vertex_iterator get_vertex_iterator() const {
        return vertex_iterator(shape_);
    }

    inline
    vertex_iterator get_vertex_end_iterator() const {
        return vertex_iterator(shape_).getEndIterator();
    }

    inline
    neighbor_vertex_iterator get_neighbor_vertex_iterator(const vertex_iterator& pos) const {
        // determine neighborhood type
        unsigned int nbtype = pos.borderType();
        // instantiate appropriate neighborhood iterator
        return neighbor_vertex_iterator(pos.point(), neighborhood[nbtype], neighborhoodIndices[nbtype], neighborIndexLookup[nbtype]);
    }

    inline
    neighbor_vertex_iterator get_neighbor_vertex_end_iterator(const vertex_iterator &pos) const {
        // determine neighborhood type
        // instantiate appropriate neighborhood iterator end
        unsigned int nbtype = pos.borderType();
        // instantiate appropriate neighborhood iterator
        return neighbor_vertex_iterator(pos.point(), neighborhood[nbtype], neighborhoodIndices[nbtype], neighborIndexLookup[nbtype], true);
    }


    // --------------------------------------------------
    // support for VertexListGraph:

    inline
    vertices_size_type
    num_vertices() const 
    {
        return prod(shape_);
    }


    // --------------------------------------------------
    // support for IncidenceGraph:

    inline
    out_edge_iterator get_out_edge_iterator(const vertex_iterator& pos) const {
        // determine neighborhood type
        unsigned int nbtype = pos.borderType();
        // instantiate appropriate neighborhood iterator
        return out_edge_iterator(pos.point(), neighborhood[nbtype], neighborhoodIndices[nbtype], neighborIndexLookup[nbtype], *this);
    }

    inline
    out_edge_iterator get_out_edge_end_iterator(const vertex_iterator &pos) const {
        // determine neighborhood type
        // instantiate appropriate neighborhood iterator end
        unsigned int nbtype = pos.borderType();
        // instantiate appropriate neighborhood iterator
        return out_edge_iterator(pos.point(), neighborhood[nbtype], neighborhoodIndices[nbtype], neighborIndexLookup[nbtype], *this, true);
    }

    inline
    const vertex_descriptor& neighborCoordOffset(int fullNeighborhoodIndex) const {
        const shape_type& neighborOffset(neighborhood[0][fullNeighborhoodIndex]);
        return neighborOffset;
    }

    inline 
    degree_size_type out_degree(const vertex_iterator &pos) const {
        // requires to fully reconstructed iterator (to accesss for neighborhood type)
        unsigned int nbtype = pos.borderType();
        return neighborhood[nbtype].size();
    }


    // --------------------------------------------------
    // support for EdgeListGraph:

    inline
    edges_size_type
    num_edges() const 
    {
        return num_edges_;
    }


    // --------------------------------------------------
    // support for AdjacencyMatrix concept:

    std::pair<edge_descriptor, bool>
    edge(const vertex_descriptor &u, const vertex_descriptor &v) const
    {
        edge_descriptor edge;
        bool found=false;

        vertex_iterator reconstructed = get_vertex_iterator();
        reconstructed += u;
        unsigned int nbtype = reconstructed.borderType();

        // check if (u-v) in neighborlist (or v-u, to save reconstruction of v!)
        shape_type diff = u-v;
        for (unsigned int i=0; i< neighborhood[nbtype].size(); ++i) {
            if (diff == neighborhood[nbtype][i]) {
                found = true;
                edge = map_to_undirected_edge(make_edge_descriptor(v, neighborIndexLookup[nbtype][i]));
                break;
            }
            else if ((-diff) == neighborhood[nbtype][i]) {
                found = true;
                // need to use edge from v to u in this case:
                edge = map_to_undirected_edge(make_edge_descriptor(u, neighborIndexLookup[nbtype][i]));
                break;
            }
        }
        return std::make_pair(edge, found);
    }


    // --------------------------------------------------
    // other helper functions:

    inline 
    degree_size_type maxDegree() const {
        // or: return max_degree;
         return neighborhood[0].size();
    }
    inline 
    degree_size_type halfMaxDegree() const {
         return maxDegree() / 2;
    }

    inline
    const shape_type& shape() const {
        return shape_;
    }

    static 
    edge_descriptor make_edge_descriptor(const vertex_descriptor &v,
                                            index_type nbindex) 
    {
        edge_descriptor res;
        TinyVectorView<typename edge_descriptor::value_type, N>(res.data()) = v;
        res[N] = nbindex;
        return res;
    }

    edge_propmap_shape_type edge_propmap_shape() const {
        edge_propmap_shape_type res;
        TinyVectorView<typename edge_propmap_shape_type::value_type, N>(res.data()) = shape_;
        // res[N] = maxDegree(); // for directed graph
        res[N] = halfMaxDegree(); // for undirected graph
        return res;
    }

    edge_propmap_shape_type directed_edge_propmap_shape() const {
        edge_propmap_shape_type res;
        TinyVectorView<typename edge_propmap_shape_type::value_type, N>(res.data()) = shape_;
        res[N] = maxDegree(); // for directed graph
        // res[N] = halfMaxDegree(); // for undirected graph
        return res;
    }

    
    //! In case of an undirected graph, the edge u->v is the same as v->u.
    //  This function folds in the edge descriptors corresponding to "causal neighbors",
    //  i.e. those to vertices with lower scan order index, by reversing the edge and computing
    //  the corresponding edge descriptor.
    //  (assuming here the neighbor-indices are in the order of noncausal, causal neighbors....
    //   FIXME: check this again in neighborhood construction!)
    //  (At least, the neighborhood construction is symmetrical, hence this should be OK)
    inline 
    edge_descriptor
    map_to_undirected_edge(const edge_descriptor &e) const {
        edge_descriptor res = e;
        //assert(!res.isReversed()); 
        if (res[N] >= halfMaxDegree()) {
            TinyVectorView<typename edge_descriptor::value_type, N> vertex(res.data());
            vertex += neighborCoordOffset(res[N]);
            res[N] = maxDegree() - res[N] - 1;
            res.setReversed(!res.isReversed());
        }
        return res;
    }



protected:
    shape_type shape_;
    size_t num_edges_;
    ArrayVector<ArrayVector<shape_type> > neighborhood;
    ArrayVector<ArrayVector<MultiArrayIndex> > neighborhoodIndices;
    ArrayVector<ArrayVector<bool> > neighborExists, causalNeighborhood, anticausalNeighborhood;
    ArrayVector<ArrayVector<int> > neighborIndexLookup;
};


    





} // namespace vigra










// Define Traits classes for BGL compatibility:
//   to obtain vertex_iterator, adjacency_iterator etc.

namespace boost {
        using namespace vigra;

        template<unsigned int N>
        inline
        std::pair<typename boost::graph_traits<vigra::GridGraph<N> >::vertex_iterator, 
                  typename boost::graph_traits<vigra::GridGraph<N> >::vertex_iterator >
        vertices(vigra::GridGraph<N> &g) 
        {
            return std::make_pair(g.get_vertex_iterator(),
                                  g.get_vertex_end_iterator());    
        }

        // const variant
        template<unsigned int N>
        inline
        std::pair<typename boost::graph_traits<vigra::GridGraph<N> >::vertex_iterator, 
                  typename boost::graph_traits<vigra::GridGraph<N> >::vertex_iterator >
        vertices(const vigra::GridGraph<N> &g) 
        {
            return std::make_pair(g.get_vertex_iterator(),
                                  g.get_vertex_end_iterator());    
        }


        
        template<unsigned int N>
        inline
        typename boost::graph_traits<vigra::GridGraph<N> >::vertices_size_type
        num_vertices(const vigra::GridGraph<N> &g) 
        {
            return g.num_vertices();
        }



        template<unsigned int N>
        inline
        std::pair<typename boost::graph_traits<vigra::GridGraph<N> >::adjacency_iterator, 
                  typename boost::graph_traits<vigra::GridGraph<N> >::adjacency_iterator >
        adjacent_vertices(typename boost::graph_traits<vigra::GridGraph<N> >::vertex_descriptor v,
                          vigra::GridGraph<N> const &g) 
        {
            // Here: need to provide a variant that converts the index vertex_descriptor
            // back into the corresponding node_iterator.
            // 
            typedef typename boost::graph_traits<vigra::GridGraph<N> >::vertex_iterator
                vertex_iterator;
            vertex_iterator reconstructed = g.get_vertex_iterator();
            reconstructed += v;
  
            return std::make_pair(g.get_neighbor_vertex_iterator(reconstructed),
                                  g.get_neighbor_vertex_end_iterator(reconstructed));    
        }


        // adjacent_vertices variant in vigra namespace: allows to call adjacent_vertices with vertex_iterator argument
        template<unsigned int N>
        inline
        std::pair<typename boost::graph_traits<vigra::GridGraph<N> >::adjacency_iterator, 
                  typename boost::graph_traits<vigra::GridGraph<N> >::adjacency_iterator >
        adjacent_vertices_at_iterator(typename boost::graph_traits<vigra::GridGraph<N> >::vertex_iterator const &v,
                                      vigra::GridGraph<N> const &g) 
        {    
            return std::make_pair(g.get_neighbor_vertex_iterator(v),
                                  g.get_neighbor_vertex_end_iterator(v));    
        }



        template<unsigned int N>
        inline
        std::pair<typename boost::graph_traits<vigra::GridGraph<N> >::out_edge_iterator, 
                  typename boost::graph_traits<vigra::GridGraph<N> >::out_edge_iterator >
        out_edges(typename boost::graph_traits<vigra::GridGraph<N> >::vertex_descriptor v,
                          vigra::GridGraph<N> const &g) 
        {
            // Here: need to provide a variant that converts the index vertex_descriptor
            // back into the corresponding node_iterator.
            // 
            typedef typename boost::graph_traits<vigra::GridGraph<N> >::vertex_iterator
                vertex_iterator;
            vertex_iterator reconstructed = g.get_vertex_iterator();
            reconstructed += v;
  
            return std::make_pair(g.get_out_edge_iterator(reconstructed),
                                  g.get_out_edge_end_iterator(reconstructed));    
        }

        template<unsigned int N>
        inline
        typename boost::graph_traits<vigra::GridGraph<N> >::vertex_descriptor 
        source_or_target(const typename boost::graph_traits<vigra::GridGraph<N> >::edge_descriptor e,
                         vigra::GridGraph<N> const &g,
                         bool return_source)
        {
            // source is always the attached node (first coords) unless the
            // edge has been reversed. 
            if ((return_source && e.isReversed()) 
                || (!return_source && !e.isReversed())) {
                typename boost::graph_traits<vigra::GridGraph<N> >::vertex_descriptor res = 
                    TinyVectorView<typename boost::graph_traits<vigra::GridGraph<N> >::edge_descriptor::value_type, N>(e.data());

                // the target is a bit more complicated, because we need the help of the graph to find the correct offset:
                res += g.neighborCoordOffset(e[N]);
                return res;
            } else {
                typename boost::graph_traits<vigra::GridGraph<N> >::vertex_descriptor res =
                    TinyVectorView<typename boost::graph_traits<vigra::GridGraph<N> >::edge_descriptor::value_type, N>(e.data());
                return res;
            }
        }


        template<unsigned int N>
        inline
        typename boost::graph_traits<vigra::GridGraph<N> >::vertex_descriptor 
        source(const typename boost::graph_traits<vigra::GridGraph<N> >::edge_descriptor e,
                          vigra::GridGraph<N> const &g) 
        {
            return source_or_target(e, g, true);
        }



        template<unsigned int N>
        inline
        typename boost::graph_traits<vigra::GridGraph<N> >::vertex_descriptor 
        target(const typename boost::graph_traits<vigra::GridGraph<N> >::edge_descriptor e,
                          vigra::GridGraph<N> const &g) 
        {
            return source_or_target(e, g, false);
        }

        template<unsigned int N>
        inline
        typename boost::graph_traits<vigra::GridGraph<N> >::degree_size_type
        out_degree(const typename boost::graph_traits<vigra::GridGraph<N> >::vertex_descriptor v,
                          vigra::GridGraph<N> const &g) 
        {
            typedef typename boost::graph_traits<vigra::GridGraph<N> >::vertex_iterator
                vertex_iterator;
            vertex_iterator reconstructed = g.get_vertex_iterator();
            reconstructed += v;
            return g.out_degree(reconstructed);
        }



        template<unsigned int N>
        inline
        std::pair<typename boost::graph_traits<vigra::GridGraph<N> >::edge_iterator, 
                  typename boost::graph_traits<vigra::GridGraph<N> >::edge_iterator >
        edges(vigra::GridGraph<N> const &g) 
        {
            typedef typename boost::graph_traits<vigra::GridGraph<N> >::edge_iterator edge_iterator;
            return std::make_pair(edge_iterator(g), edge_iterator());
        }


        template<unsigned int N>
        inline
        typename boost::graph_traits<vigra::GridGraph<N> >::edges_size_type
        num_edges(vigra::GridGraph<N> const &g) 
        {
            return g.num_edges();
        }


        // --------------------------------------------------
        // support for AdjacencyMatrix concept:

        template<unsigned int N>
        std::pair<typename boost::graph_traits<vigra::GridGraph<N> >::edge_descriptor, bool>
        edge(const typename boost::graph_traits<vigra::GridGraph<N> >::vertex_descriptor &u,
             const typename boost::graph_traits<vigra::GridGraph<N> >::vertex_descriptor &v,
             vigra::GridGraph<N> const &g)
        {
            return g.edge(u,v);
        }



        // provide get / put for MultiArrayViews, indexed by the above-defined vertex_descriptor (in this case, a coordinate tuple):

        template<class VIEW>
        class MultiArrayView_property_map 
        {
        public:
            // typedef boost::read_write_property_map_tag category;
            typedef boost::lvalue_property_map_tag category;
            //    typedef int value_type;
            typedef typename VIEW::value_type value_type;
            typedef typename VIEW::reference reference;
            typedef typename VIEW::const_reference const_reference;
            typedef typename vigra::GridGraph<VIEW::actual_dimension> graph_type;
            typedef typename graph_type::vertex_descriptor key_type;
            MultiArrayView_property_map(const VIEW& view)
                : view_(view) { }
            template <class T2>
            inline
            reference operator[](const T2 & x) { return view_[x]; } 
            template <class T2>
            inline
            const_reference operator[](const T2 & x) const { return view_[x]; } 
    
        protected:
            VIEW view_;
        };

        template<class VIEW>
        inline
        void put(MultiArrayView_property_map<VIEW> &pmap,
                 const typename MultiArrayView_property_map<VIEW>::key_type &k,
                 const typename MultiArrayView_property_map<VIEW>::value_type& val
                 ) 
        { 
            pmap[k] = val;  
        }

        template<class VIEW>
        inline
        typename MultiArrayView_property_map<VIEW>::const_reference 
        get(
            const MultiArrayView_property_map<VIEW> & pmap, 
            const typename MultiArrayView_property_map<VIEW>::key_type &k)
        { 
            return pmap[k]; 
        }


        //! use a MultiArrayView as an undirected edge property map
        //  (edge_descriptor keys will be transformed by "flipping" the edges if necessary)
        template<class VIEW, class GRAPH>
        class MultiArrayView_undirected_edge_property_map 
        {
        public:
            // typedef boost::read_write_property_map_tag category;
            typedef boost::lvalue_property_map_tag category;
            //    typedef int value_type;
            typedef typename VIEW::value_type value_type;
            typedef typename VIEW::reference reference;
            typedef typename VIEW::const_reference const_reference;
            typedef GRAPH graph_type;
            typedef typename graph_type::edge_descriptor key_type;
            MultiArrayView_undirected_edge_property_map(const VIEW& view, const GRAPH& graph)
                : view_(view), graph_(graph) { }
            template <class T2>
            inline
            reference operator[](const T2 & x) { return view_[graph_.map_to_undirected_edge(x)]; } 
            template <class T2>
            inline
            const_reference operator[](const T2 & x) const { return view_[graph_.map_to_undirected_edge(x)]; } 
    
        protected:
            VIEW view_;
            const GRAPH &graph_;
        };

        template<class VIEW, class GRAPH>
        inline
        void put(MultiArrayView_undirected_edge_property_map<VIEW, GRAPH> &pmap,
                 const typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::key_type &k,
                 const typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::value_type& val
                 ) 
        { 
            pmap[k] = val;  
        }

        template<class VIEW, class GRAPH>
        inline
        typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::const_reference 
        get(
            const MultiArrayView_undirected_edge_property_map<VIEW, GRAPH> & pmap, 
            const typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::key_type &k)
        { 
            return pmap[k]; 
        }




        
        // property map support for mapping coordinates to scan-order indices:

        template<unsigned int N>
        struct IDMapper {
            typedef typename vigra::GridGraph<N> graph_type;
            typedef boost::readable_property_map_tag category;
            typedef typename graph_type::index_type value_type;
            typedef typename graph_type::vertex_descriptor key_type;
            typedef const value_type& reference;


            IDMapper(const graph_type &graph) 
                : map_helper(graph.get_vertex_iterator())
            {}

            typename graph_type::vertex_iterator map_helper;
        };

        template<unsigned int N>
        struct property_map<vigra::GridGraph<N>, boost::vertex_index_t>
        {
            typedef IDMapper<N> type;
            typedef IDMapper<N> const_type;
        };


        template<unsigned int N>
        inline
        typename IDMapper<N>::value_type
        get(const IDMapper<N> & mapper, 
            const typename IDMapper<N>::key_type &k)
        { 
            return (mapper.map_helper + k).scanOrderIndex();
        }


        template<unsigned int N>
        typename boost::property_map<vigra::GridGraph<N>, boost::vertex_index_t>::type
        //typename IDMapper<N>
        get(boost::vertex_index_t, const vigra::GridGraph<N> &graph) {
            // return a lightweight wrapper for the CoupledIterator, which easily allows the conversion of 
            // coordinates via its += operator followed by index().
            return IDMapper<N>(graph);
        }

#if 0 
        // CHECK if required: also provide the direct (three-parameter) version for index lookup
        template<unsigned int N>
        typename vigra::GridGraph<N>::vertices_size_type
        get(boost::vertex_index_t, 
            const vigra::GridGraph<N> &graph,
            const typename vigra::GridGraph<N>::vertex_descriptor &v) {
            return (IDMapper<N>(graph).map_helper + v).scanOrderIndex();
        }
#endif



        // TODO:
        // eventually provide an edge_index property map as well?
        // (edge_descriptor -> linear contiguous edge index)


} // namespace boost




// FIXME: I'd rather like the definition of this iterator in 
// a single file with the declaration, however I didn't manage
// to resolve the circular dependency between the two templates
// in another way. Have to try some harder sometime else.        
#include "multi_gridgraph_coords_edge_iterator_defn.hxx"



namespace std {
    template<unsigned int N>
    ostream& operator<<(ostream& out,
                        const typename vigra::GridGraph<N>::vertex_iterator & arg)
    {
        out << "v" << arg.scanOrderIndex();
        return out;
    }
};

#endif


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



#endif /* VIGRA_MULTI_GRIDGRAPH_HXX */



