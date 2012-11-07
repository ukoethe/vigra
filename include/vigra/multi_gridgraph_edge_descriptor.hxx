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

#ifndef VIGRA_MULTI_GRIDGRAPH_EDGE_DESCRIPTOR_HXX
#define VIGRA_MULTI_GRIDGRAPH_EDGE_DESCRIPTOR_HXX

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

*/

#include "multi_iterator.hxx"
#include "array_vector.hxx"

namespace vigra {

namespace detail {

template<unsigned int N>
class GridGraphEdgeDescriptor
    : public typename MultiArrayShape<N+1>::type
{
  public:
    typedef typename MultiArrayShape<N+1>::type  base_type;
    typedef typename base_type::value_type       value_type;
    typedef base_type                            edge_coord_type;
    typedef value_type                           index_type;
    typedef typename MultiArrayShape<N>::type    shape_type;
    // TODO: expose constructors

    GridGraphEdgeDescriptor()
    : is_reversed_(false)
    {}

    GridGraphEdgeDescriptor(shape_type const &vertex,
                            index_type edge_index,
                            bool reversed=false)
    : base_type(DontInit()),
      is_reversed_(reversed)
    {
        TinyVectorView<value_type, N>(this->data()) = vertex;
        *this[N] = edge_index;
    }
                                      
    void set(shape_type const &vertex, index_type edge_index, bool reversed) 
    {
        TinyVectorView<value_type, N>(this->data()) = vertex;
        *this[N] = edge_index;
        is_reversed_ = reversed;
    }
        
    bool isReversed() const 
    {
        return is_reversed_;
    }

  protected:
    bool is_reversed_;
};

} // namespace detail

} // namespace vigra

#endif // VIGRA_MULTI_GRIDGRAPH_EDGE_DESCRIPTOR_HXX
