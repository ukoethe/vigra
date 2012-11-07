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

#ifndef VIGRA_MULTI_GRIDGRAPH_VERTEX_ITERATOR_HXX
#define VIGRA_MULTI_GRIDGRAPH_VERTEX_ITERATOR_HXX

#include "multi_iterator_coupled.hxx"

namespace vigra {

namespace detail {

template<unsigned int N>
class GridGraphVertexIterator
    : public typename CoupledIteratorType<N>::type
{
  public:
    typedef typename CoupledIteratorType<N>::type  base_type;

    typedef typename base_type::shape_type         shape_type;
    typedef typename base_type::difference_type    difference_type;
    typedef GridGraphVertexIterator                iterator;
    typedef std::random_access_iterator_tag        iterator_category;

    typedef typename base_type::value_type         handle_type;
    typedef typename handle_type::value_type       value_type;
    typedef typename handle_type::reference        reference;
    typedef typename handle_type::const_reference  const_reference;
    typedef typename handle_type::pointer          pointer;

    GridGraphVertexIterator(handle_type const & handles = handle_type()) 
        : base_type(handles)
    {}

    GridGraphVertexIterator(base_type const & base) 
        : base_type(base)
    {}

    // dereferencing the iterator yields the coordinate object
    // (used as vertex_descriptor)
    reference operator*()
    {
        return this->template get<0>();
    }
    
    const_reference operator*() const
    {
        return this->template get<0>();
    }

    reference operator[](MultiArrayIndex i)
    {
        return *(GridGraphVertexIterator(*this) += i);
    }

    const_reference operator[](MultiArrayIndex i) const
    {
        return *(GridGraphVertexIterator(*this) += i);
    }

    GridGraphVertexIterator & operator++()
    {
        base_type::operator++();
        return *this;
    }
    
    GridGraphVertexIterator operator++(int)
    {
        GridGraphVertexIterator res(*this);
        ++*this;
        return res;
    }

    GridGraphVertexIterator & operator+=(MultiArrayIndex i)
    {
        base_type::operator+=(i);
        return *this;
    }

    GridGraphVertexIterator & operator+=(const shape_type &coordOffset)
    {
        base_type::operator+=(coordOffset);
        return *this;
    }

    GridGraphVertexIterator & operator--()
    {
        base_type::operator--();
        return *this;
    }

    GridGraphVertexIterator operator--(int)
    {
        GridGraphVertexIterator res(*this);
        --*this;
        return res;
    }

    GridGraphVertexIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }

    GridGraphVertexIterator & operator-=(const shape_type &coordOffset)
    {
        return operator+=(-coordOffset);
    }

    GridGraphVertexIterator getEndIterator() const
    {
        return GridGraphVertexIterator(base_type::getEndIterator());
    }

    GridGraphVertexIterator operator+(MultiArrayIndex d) const
    {
        return GridGraphVertexIterator(*this) += d;
    }

    GridGraphVertexIterator operator-(MultiArrayIndex d) const
    {
        return GridGraphVertexIterator(*this) -= d;
    }

    GridGraphVertexIterator operator+(const shape_type &coordOffset) const
    {
        return GridGraphVertexIterator(*this) += coordOffset;
    }

    GridGraphVertexIterator operator-(const shape_type &coordOffset) const
    {
        return GridGraphVertexIterator(*this) -= coordOffset;
    }
};

} // namespace detail

} // namespace vigra

#endif // VIGRA_MULTI_GRIDGRAPH_VERTEX_ITERATOR_HXX
