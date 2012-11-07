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

#ifndef VIGRA_MULTI_GRIDGRAPH_OUT_EDGE_ITERATOR_HXX
#define VIGRA_MULTI_GRIDGRAPH_OUT_EDGE_ITERATOR_HXX

#include "multi_gridgraph_edge_descriptor.hxx"
#include <cassert>

namespace vigra {

namespace detail {

template<unsigned int N>
class GridGraphOutEdgeIterator
{
  public:
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef MultiArrayIndex                    index_type;
    typedef GridGraphEdgeDescriptor<N>         value_type;
    typedef value_type const &                 reference;
    typedef value_type const &                 const_reference;
    typedef value_type cpnst *                 pointer;
    typedef value_type const *                 const_pointer;
    typedef std::forward_iterator_tag          iterator_category;

    GridGraphOutEdgeIterator() 
    : neighborOffsets_(0),
      neighborIndices_(0),
      index_(0),
      max_index_(0),
      directed_(false)
    {}

    GridGraphOutEdgeIterator(const ArrayVector<shape_type> &neighborOffsets,
                             const ArrayVector<index_type> &neighborIndices,
                             shape_type const & pos,
                             bool directed)
    : neighborOffsets_(&neighborOffsets),
      neighborIndices_(&neighborIndices),
      edge_descriptor_(pos, 0),
      pos_(pos),
      index_(0),
      max_index_((index_type)neighborOffsets.size()-1),
      directed_(directed)
    {
        updateEdgeDescriptor();
    }

    GridGraphOutEdgeIterator & operator++()
    {
        ++index_;
        updateEdgeDescriptor();
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

    index_type index() const
    {
        return index_;
    }

    index_type neighborIndex() const
    {
        return (*neighborIndices_)[index_];
    }

    bool operator==(self_type const & other) const
    {
        return index_ == other.index();
    }

    bool operator!=(self_type const & other) const
    {
        return index_ != other.index();
    }

    bool hasMore() const 
    {
        return index_ < (index_type)neighborIndices_->size();
    }

    GridGraphOutEdgeIterator getEndIterator() const
    {
        GridGraphOutEdgeIterator res(*this);
        res.index_ = (index_type)neighborIndices_->size();
        return res;
    }

  protected:
    void updateEdgeDescriptor()
    {
        if(hasMore())
        {
            index_type neighbor_index = neighborIndex();
            if(directed_ || 2*neighbor_index < max_index_)
                edge_descriptor_[N] = neighbor_index;
            else
                edge_descriptor_.set(pos_ + (*neighborOffsets_)[neighbor_index],
                                     max_index_ - neighbor_index,
                                     true);
        }
    }
  
    const ArrayVector<shape_type> *neighborOffsets_;
    const ArrayVector<index_type> *neighborIndices_;
    value_type edge_descriptor_;
    shape_type pos_;
    index_type index_, max_index_;
    bool directed_;
};

} // namespace detail

} // namespace vigra

#endif // VIGRA_MULTI_GRIDGRAPH_OUT_EDGE_ITERATOR_HXX
