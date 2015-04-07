/************************************************************************/
/*                                                                      */
/*    Copyright 2003-2013 by Thorben Kroeger, Kasim Terzic,             */
/*        Christian-Dennis Rahn and Ullrich Koethe                      */
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

#ifndef VIGRA_VECTOR_DISTANCE_HXX
#define VIGRA_VECTOR_DISTANCE_HXX

#include <vector>
#include <set>
#include <functional>
#include "array_vector.hxx"
#include "multi_array.hxx"
#include "accessor.hxx"
#include "numerictraits.hxx"
#include "navigator.hxx"
#include "metaprogramming.hxx"
#include "multi_pointoperators.hxx"
#include "functorexpression.hxx"
#include "multi_distance.hxx"

#undef VECTORIAL_DIST_DEBUG 

namespace vigra
{

namespace detail
{

template <class Vector, class Value>
struct VectorialDistParabolaStackEntry
{
    double left, center, right;
    Value apex_height;
    Vector point;
    
    VectorialDistParabolaStackEntry(const Vector& vec, Value prev, double l, double c, double r)
    : left(l), center(c), right(r), apex_height(prev), point(vec)
    {}
};

#ifdef VECTORIAL_DIST_DEBUG
template <class Vector, class Value>
std::ostream& operator<<(std::ostream&o, const VectorialDistParabolaStackEntry<Vector, Value>& e) {
    o << "l=" << e.left << ", c=" << e.center << ", r=" << e.right << ", pV=" << e.apex_height << ", pVec=" << e.point;
    return o;
}
#endif

/********************************************************/
/*                                                      */
/*               vectorialDistParabola                  */
/*                                                      */
/********************************************************/

template <class VEC, class ARRAY>
inline double 
partialSquaredMagnitude(const VEC& vec, MultiArrayIndex dim, ARRAY const & pixel_pitch)
{
    //computes the squared magnitude of vec
    //considering only the first dim dimensions
    double sqMag = 0.0; 
    for(MultiArrayIndex i=0; i<=dim; ++i)
    {
        sqMag += sq(pixel_pitch[i]*vec[i]);
    }
    return sqMag;
}

template <class SrcIterator,
          class Array>
void 
vectorialDistParabola(MultiArrayIndex dimension,
                      SrcIterator is, SrcIterator iend,
                      Array const & pixel_pitch )
{
    typedef typename SrcIterator::value_type SrcType;
    typedef VectorialDistParabolaStackEntry<SrcType, double> Influence;
    
    double sigma = pixel_pitch[dimension],
           sigma2 = sq(sigma);
    double w = iend - is; //width of the scanline
    
    SrcIterator id = is;
   
    std::vector<Influence> _stack; //stack of influence parabolas
    double apex_height = partialSquaredMagnitude(*is, dimension, pixel_pitch);
    _stack.push_back(Influence(*is, apex_height, 0.0, 0.0, w));
    ++is;
    double current = 1.0;
    while(current < w)
    {
        apex_height = partialSquaredMagnitude(*is, dimension, pixel_pitch);
        Influence & s = _stack.back();
        double diff = current - s.center;
        double intersection = current + (apex_height - s.apex_height - sq(sigma*diff)) / (2.0*sigma2 * diff);
        
        if( intersection < s.left) // previous point has no influence
        {
            _stack.pop_back();
            if(_stack.empty())
                _stack.push_back(Influence(*is, apex_height, 0.0, current, w));
            else
                continue; // try new top of stack without advancing current
        }
        else if(intersection < s.right)
        {
            s.right = intersection;
            _stack.push_back(Influence(*is, apex_height, intersection, current, w));
        }
        ++is;
        ++current;
    }
   
    // Now we have the stack indicating which rows are influenced by (and therefore
    // closest to) which row. We can go through the stack and calculate the
    // distance squared for each element of the column.
    typename std::vector<Influence>::iterator it = _stack.begin();
    for(current = 0.0; current < w; ++current, ++id)
    {
        while( current >= it->right)
            ++it; 
        
        *id = it->point;
        (*id)[dimension] = it->center - current;
    }
}

template <class DestIterator, 
          class LabelIterator,
          class Array1, class Array2>
void 
boundaryVectorDistParabola(MultiArrayIndex dimension,
                           DestIterator is, DestIterator iend, 
                           LabelIterator ilabels,
                           Array1 const & pixel_pitch,
                           Array2 const & dmax,
                           bool array_border_is_active=false)
{
    double w = iend - is;
    if(w <= 0)
        return;

    typedef typename LabelIterator::value_type LabelType;
    typedef typename DestIterator::value_type DestType;
    typedef VectorialDistParabolaStackEntry<DestType, double> Influence;
    typedef std::vector<Influence> Stack;

    DestIterator id = is;
    DestType border_point = array_border_is_active
                                ? DestType(0)
                                : dmax;
    double apex_height = partialSquaredMagnitude(border_point, dimension, pixel_pitch);
    Stack _stack(1, Influence(border_point, apex_height, 0.0, -1.0, w));
    LabelType current_label = *ilabels;
    for(double begin = 0.0, current = 0.0; current <= w; ++ilabels, ++is, ++current)
    {
        DestType point = (current < w)
                             ? (current_label == *ilabels)
                                 ? *is
                                 : DestType(0)
                             : border_point;
        apex_height = partialSquaredMagnitude(point, dimension, pixel_pitch);
        while(true)
        {
            Influence & s = _stack.back();
            double diff = (current - s.center)*pixel_pitch[dimension];
            double intersection = current + (apex_height - s.apex_height - sq(diff)) / (2.0 * diff);
            
            if(intersection < s.left) // previous parabola has no influence
            {
                _stack.pop_back();
                if(_stack.empty())
                    intersection = begin; // new parabola is valid for entire present segment
                else
                    continue;  // try new top of stack without advancing to next pixel
            }
            else if(intersection < s.right)
            {
                s.right = intersection;
            }
            if(intersection < w)
                _stack.push_back(Influence(point, apex_height, intersection, current, w));
            if(current < w && current_label == *ilabels)
                break; // finished present pixel, advance to next one
                
            // label changed => finalize the current segment
            typename Stack::iterator it = _stack.begin();
            for(double c = begin; c < current; ++c, ++id)
            {
                while(c >= it->right) 
                    ++it; 
                *id = it->point;
                (*id)[dimension] = it->center - c;
            }
            if(current == w)
                break;  // stop when this was the last segment
                
            // initialize the new segment
            begin = current;
            current_label = *ilabels;
            point = *is;
            apex_height = partialSquaredMagnitude(point, dimension, pixel_pitch);
            Stack(1, Influence(DestType(0), 0.0, begin-1.0, begin-1.0, w)).swap(_stack);
            // don't advance to next pixel here, because the present pixel must also 
            // be analysed in the context of the new segment
        }
    }
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2,
          class Array>
void 
interpixelBoundaryVectorDistance(MultiArrayView<N, T1, S1> const & labels,
                                 MultiArrayView<N, T2, S2> dest, 
                                 Array const & pixelPitch)
{
    typedef typename MultiArrayShape<N>::type  Shape;
    typedef GridGraph<N>                       Graph;
    typedef typename Graph::Node               Node;
    typedef typename Graph::NodeIt             graph_scanner;
    typedef typename Graph::OutArcIt           neighbor_iterator;
    
    Graph g(labels.shape());
    for (graph_scanner node(g); node != lemon_graph::INVALID; ++node)
    {
        T1 label = labels[*node];
        double min_dist = NumericTraits<double>::max();
        Node point    = *node,
             boundary = point + Node(dest[point]),
             min_pos  = lemon::INVALID;
        T2 min_diff;

        //go to adjacent neighbour with same label as origin pixel with smallest distance
        if(labels.isInside(boundary))
        {
            for (neighbor_iterator arc(g, boundary); arc != lemon_graph::INVALID; ++arc)
            {
                if(label == labels[g.target(*arc)])
                {
                    double dist = squaredNorm(pixelPitch*(g.target(*arc) - point));
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        min_pos = g.target(*arc);
                    }
                }
            }
            if(min_pos == lemon::INVALID)
                continue;
            min_dist = NumericTraits<double>::max();
        }
        else
        {
            min_pos = clip(boundary, Shape(0), labels.shape()-Shape(1));
            min_diff = 0.5*(boundary + min_pos) - point;
            min_dist = squaredNorm(pixelPitch*min_diff);
        }
        
        //from this pixel look for the vector which points to the nearest interpixel between two label
        for (neighbor_iterator arc(g, min_pos); arc != lemon_graph::INVALID; ++arc)
        {
            if(label != labels[g.target(*arc)])
            {
                T2 diff = 0.5*(g.target(*arc) + min_pos) - point;
                double dist = squaredNorm(pixelPitch*diff);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    min_diff = diff;
                }
            }
        }
        dest[point] = min_diff;
    }
}

} // namespace detail

/** \addtogroup MultiArrayDistanceTransform
*/
//@{

template<bool PRED>
struct Error_output_pixel_type_must_be_TinyVector_of_appropriate_length 
: vigra::staticAssert::AssertBool<PRED> {};

/********************************************************/
/*                                                      */
/*               separableVectorDistance                */
/*                                                      */
/********************************************************/

    /** \brief Compute the vector distance transform of a N-dimensional binary array.

        <b> Declarations:</b>

        \code
        namespace vigra {
            template <unsigned int N, class T1, class S1,
                      class T2, class S2, class Array>
            void 
            separableVectorDistance(MultiArrayView<N, T1, S1> const & source,
                                    MultiArrayView<N, T2, S2> dest, 
                                    bool background,
                                    Array const & pixelPitch=TinyVector<double, N>(1));
        }
        \endcode

        This function works like \ref separableMultiDistance() (see there for details),
        but returns in each pixel the <i>vector</i> to the nearest background pixel 
        rather than the scalar distance. This enables much more powerful applications.

        <b> Usage:</b>

        <b>\#include</b> \<vigra/vector_distance.hxx\><br/>
        Namespace: vigra

        \code
        Shape3 shape(width, height, depth);
        MultiArray<3, unsigned char> source(shape);
        MultiArray<3, Shape3> dest(shape);
        ...

        // For each background pixel, find the vector to the nearest foreground pixel.
        separableVectorDistance(source, dest, true);
        \endcode

        \see vigra::separableMultiDistance(), vigra::boundaryVectorDistance()
    */
doxygen_overloaded_function(template <...> void separableVectorDistance)

template <unsigned int N, class T1, class S1,
          class T2, class S2, class Array>
void 
separableVectorDistance(MultiArrayView<N, T1, S1> const & source,
                        MultiArrayView<N, T2, S2> dest, 
                        bool background,
                        Array const & pixelPitch)
{
    using namespace vigra::functor;
    typedef typename MultiArrayView<N, T2, S2>::traverser Traverser;
    typedef MultiArrayNavigator<Traverser, N> Navigator;

    VIGRA_STATIC_ASSERT((Error_output_pixel_type_must_be_TinyVector_of_appropriate_length<N == T2::static_size>));
    vigra_precondition(source.shape() == dest.shape(),
        "separableVectorDistance(): shape mismatch between input and output.");
    vigra_precondition(pixelPitch.size() == N, 
        "separableVectorDistance(): pixelPitch has wrong length.");
        
    T2 maxDist(2*sum(source.shape()*pixelPitch)), rzero;
    if(background == true)
        transformMultiArray( source, dest,
                                ifThenElse( Arg1() == Param(0), Param(maxDist), Param(rzero) ));
    else
        transformMultiArray( source, dest, 
                                ifThenElse( Arg1() != Param(0), Param(maxDist), Param(rzero) ));
    
    for(int d = 0; d < N; ++d )
    {
        Navigator nav( dest.traverser_begin(), dest.shape(), d);
        for( ; nav.hasMore(); nav++ )
        {
             detail::vectorialDistParabola(d, nav.begin(), nav.end(), pixelPitch);
        }
    }
}

template <unsigned int N, class T1, class S1,
          class T2, class S2>
inline void 
separableVectorDistance(MultiArrayView<N, T1, S1> const & source,
                        MultiArrayView<N, T2, S2> dest, 
                        bool background=true)
{
    TinyVector<double, N> pixelPitch(1.0);
    separableVectorDistance(source, dest, background, pixelPitch);
}


    /** \brief Compute the vector distance transform to the implicit boundaries of a 
               multi-dimensional label array.

        <b> Declarations:</b>

        \code
        namespace vigra {
            template <unsigned int N, class T1, class S1,
                                      class T2, class S2,
                      class Array>
            void
            boundaryVectorDistance(MultiArrayView<N, T1, S1> const & labels,
                                   MultiArrayView<N, T2, S2> dest,
                                   bool array_border_is_active=false,
                                   BoundaryDistanceTag boundary=OuterBoundary,
                                   Array const & pixelPitch=TinyVector<double, N>(1));
        }
        \endcode

        This function works like \ref boundaryMultiDistance() (see there for details),
        but returns in each pixel the <i>vector</i> to the nearest boundary pixel 
        rather than the scalar distance. This enables much more powerful applications.
        Additionally, it support a <tt>pixelPitch</tt> parameter which allows to adjust
        the distance calculations for anisotropic grid resolution.

        <b> Usage:</b>

        <b>\#include</b> \<vigra/vector_distance.hxx\><br/>
        Namespace: vigra

        \code
        Shape3 shape(width, height, depth);
        MultiArray<3, UInt32> labels(shape);
        MultiArray<3, Shape3> dest(shape);
        ...

        // For each region, find the vectors to the nearest boundary pixel, including the 
        // outer border of the array.
        boundaryVectorDistance(labels, dest, true);
        \endcode

        \see vigra::boundaryMultiDistance(), vigra::separableVectorDistance()
    */
doxygen_overloaded_function(template <...> void boundaryVectorDistance)

template <unsigned int N, class T1, class S1,
                          class T2, class S2,
          class Array>
void
boundaryVectorDistance(MultiArrayView<N, T1, S1> const & labels,
                       MultiArrayView<N, T2, S2> dest,
                       bool array_border_is_active,
                       BoundaryDistanceTag boundary,
                       Array const & pixelPitch)
{
    VIGRA_STATIC_ASSERT((Error_output_pixel_type_must_be_TinyVector_of_appropriate_length<N == T2::static_size>));
    vigra_precondition(labels.shape() == dest.shape(),
        "boundaryVectorDistance(): shape mismatch between input and output.");
    vigra_precondition(pixelPitch.size() == N, 
        "boundaryVectorDistance(): pixelPitch has wrong length.");
        
    using namespace vigra::functor;
    
    if(boundary == InnerBoundary)
    {
        MultiArray<N, unsigned char> boundaries(labels.shape());
        
        markRegionBoundaries(labels, boundaries, IndirectNeighborhood);
        if(array_border_is_active)
            initMultiArrayBorder(boundaries, 1, 1);
        separableVectorDistance(boundaries, dest, true, pixelPitch);
    }
    else
    {
        if(boundary == InterpixelBoundary)
        {
            vigra_precondition(!NumericTraits<T2>::isIntegral::value,
                "boundaryVectorDistance(..., InterpixelBoundary): output pixel type must be float or double.");
        }
        
        typedef typename MultiArrayView<N, T1, S1>::const_traverser LabelIterator;
        typedef typename MultiArrayView<N, T2, S2>::traverser DestIterator;
        typedef MultiArrayNavigator<LabelIterator, N> LabelNavigator;
        typedef MultiArrayNavigator<DestIterator, N> DNavigator;
        
        T2 maxDist(2*sum(labels.shape()*pixelPitch));
        dest = maxDist;
        for( int d = 0; d < N; ++d )
        {
            LabelNavigator lnav( labels.traverser_begin(), labels.shape(), d );
            DNavigator dnav( dest.traverser_begin(), dest.shape(), d );

            for( ; dnav.hasMore(); dnav++, lnav++ )
            {
                detail::boundaryVectorDistParabola(d, dnav.begin(), dnav.end(), lnav.begin(), 
                                                   pixelPitch, maxDist, array_border_is_active);
            }
        }
        
        if(boundary == InterpixelBoundary)
        {
           detail::interpixelBoundaryVectorDistance(labels, dest, pixelPitch);
        }
    }
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2>
void
boundaryVectorDistance(MultiArrayView<N, T1, S1> const & labels,
                       MultiArrayView<N, T2, S2> dest,
                       bool array_border_is_active=false,
                       BoundaryDistanceTag boundary=OuterBoundary)
{
    TinyVector<double, N> pixelPitch(1.0);
    boundaryVectorDistance(labels, dest, array_border_is_active, boundary, pixelPitch);
}

//@}

} //-- namespace vigra

#endif        //-- VIGRA_VECTOR_DISTANCE_HXX
