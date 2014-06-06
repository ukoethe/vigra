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

#ifndef VIGRA_VECTORIAL_BOUNDARY_DISTANCE_HXX
#define VIGRA_VECTORIAL_BOUNDARY_DISTANCE_HXX

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

#undef VECTORIAL_BOUNDARY_DIST_DEBUG

namespace vigra
{

namespace detail {
    
template <class Vector, class Value>
struct VectorialBoundaryDistParabolaStackEntry
{
    double left, center, right;
    Value prevVal;
    Vector prevVector;
    bool boundary;
    bool corner1;
    bool corner2;
    
    VectorialBoundaryDistParabolaStackEntry(const Vector& vec, Value prev, double l, double c, double r, bool b, bool c1, bool c2)
    : left(l), center(c), right(r), prevVal(prev), prevVector(vec), boundary(b), corner1(c1), corner2(c2)
    {}
};

#ifdef VECTORIAL_BOUNDARY_DIST_DEBUG
template <class Vector, class Value>
std::ostream& operator<<(std::ostream&o, const VectorialBoundaryDistParabolaStackEntry<Vector, Value>& e) {
    o << "l=" << e.left << ", c=" << e.center << ", r=" << e.right
      << ", prevVal=" << e.prevVal
      << ", prevVector=" << e.prevVector
      << ", b=" << e.boundary
      << ", c1=" << e.corner1
      << ", c2=" << e.corner2;
    return o;
}
#endif

/********************************************************/
/*                                                      */
/*               vectorialBoundaryDistParabola          */
/*                                                      */
/********************************************************/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class SegmentationIterator,
          class CornersIterator>
void vectorialBoundaryDistParabola(
    MultiArrayIndex dimension,
    size_t maxDist,
    SrcIterator is, SrcIterator iend, SrcAccessor sa,
    DestIterator id, DestAccessor da,
    SegmentationIterator segIter,
    CornersIterator cornersIter1,
    CornersIterator cornersIter2
) {
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    typedef typename SrcType::value_type SrcPixelType;
    typedef VectorialBoundaryDistParabolaStackEntry<DestType, SrcPixelType> Influence;
    
    double w = iend - is; //width of the scanline
    
    #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG
    {
        std::cout << "* vectorialBoundaryDistParabola (dim=" << dimension << ", width=" << w << ")" << std::endl;
        
        std::cout << "  source line: ";
        SrcIterator is_copy(is);
        for(; is != iend; ++is) {
            std::cout << sa(is) << " ";
        }
        std::cout << std::endl;
        is = is_copy; 
        
        std::cout << "  seg line:    " << std::endl;
        SegmentationIterator segIter_copy(segIter);
        CornersIterator cornersIter1_copy(cornersIter1);
        CornersIterator cornersIter2_copy(cornersIter2);
        
        for(int i=0; i<w+1; ++i, ++cornersIter1) {
            std::cout << (int)*(cornersIter1) << " ";
        }
        std::cout << std::endl;
        std::cout << " ";
        for(int i=0; i<w; ++i, ++segIter) {
            std::cout << *segIter << " ";
        }
        std::cout << std::endl;
        for(int i=0; i<w+1; ++i, ++cornersIter2) {
            std::cout << (int)*(cornersIter2) << " ";
        }
        std::cout << std::endl;
        segIter = segIter_copy;
        cornersIter1 = cornersIter1_copy;
        cornersIter2 = cornersIter2_copy;
        
    }
    #endif
    
    std::vector<Influence> _stack; //stack of influence parabolas
    
    //variables needed for loop
    double psm = 0.0;
    SrcType srcValue;
    bool   boundary = true;
    double current = 0.0;
    double currentInterpixel = -0.5;
   
    srcValue = SrcType(static_cast<typename SrcType::value_type>(0));
    srcValue[dimension] = 0.5;
    
    psm = partialSquaredMagnitude(srcValue, dimension+1);
    _stack.push_back(Influence(srcValue, psm, currentInterpixel, currentInterpixel, w, true, *cornersIter1, *cornersIter2));
    #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG 
    std::cout << "  " << "pushed " << _stack.back() << " (before while)" << std::endl;
    #endif
    
    ++cornersIter1, ++cornersIter2;
    
    while(current < w)
    {
        srcValue = sa(is);
        currentInterpixel = current+0.5;
       
        //determine if we have an implicit boundary on this line
        if(current < w-1) {
            boundary = (*segIter != *(segIter+1));
        }
        else {
            boundary = true; //implicit boundary at image edge
        }
        
        //if we have an implicit boundary
        if(boundary) {
            srcValue = SrcType(static_cast<typename SrcType::value_type>(0));
            srcValue[dimension] = 0.5;
        }
        else if(*cornersIter1 && !*cornersIter2) {
            srcValue[dimension] = std::min(0.5, srcValue[dimension]);
            if(dimension < 1) {
                srcValue[dimension+1] = -0.5;
            }
        }
        else if(!*cornersIter1 && *cornersIter2) {
            srcValue[dimension] = std::min(0.5, srcValue[dimension]);
            if(dimension < 1) {
                srcValue[dimension+1] = 0.5;
            }
        }
        else if(*cornersIter1 && *cornersIter2) {
            //this can occur:
            //
            // AABB
            // CCCC
            // BBAA
            srcValue[dimension] = 0.5;
            srcValue[dimension+1] = 0.5;
        }
        
        SrcType xxx = srcValue;
        xxx[dimension] = 0.0;
        psm = partialSquaredMagnitude(xxx, dimension+1);
        
        Influence & s = _stack.back();
        double diff = current - s.center;
        //Bailey 2004, eq. (14)
        //Compute the intersection of the two parabolas
        double intersection = current + ( psm - s.prevVal - sq(diff) ) / (2.0* diff);
        
        #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG 
        std::cout << "  "
                  << "current=" << current
                  << ", currentInterpixel=" << currentInterpixel
                  << ", intersection=" << intersection
                  << ", srcValue=" << srcValue
                  << ", psm=" << psm
                  << " | boundary=" << boundary
                  << ", corner=" << (int)*cornersIter1
                  << std::endl;
        #endif
        
        if( intersection < s.left) // previous point has no influence
        {
            _stack.pop_back();
            #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG 
            std::cout << "  " << "popped " << _stack.back() << std::endl;
            #endif
            if(_stack.empty())
            {
                _stack.push_back(Influence(srcValue, psm, 0.0, currentInterpixel, w, boundary, *cornersIter1, *cornersIter2));
                #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG 
                std::cout << "  " << "pushed " << _stack.back() << std::endl;
                #endif
            }
            else
            {
                continue; // try new top of stack without advancing current
            }
        }
        else if(intersection < s.right)
        {
            s.right = intersection;
            _stack.push_back(Influence(srcValue, psm, intersection, currentInterpixel, w, boundary, *cornersIter1, *cornersIter2));
            #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG 
            std::cout << "  " << "pushed " << _stack.back() << " reason 2 " << std::endl;
            #endif
        }
        ++is; ++segIter, ++cornersIter1, ++cornersIter2;
        ++current;
    }
   
    #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG
    std::cout << ">>> stack for dimension=" << dimension << std::endl;
    for(int i=0; i<_stack.size(); ++i) {
        std::cout << "    stack " << i << " = " << _stack[i] << std::endl;
    }
    std::cout << "    write line as " << std::endl;
    #endif
   
    // Now we have the stack indicating which rows are influenced by (and therefore
    // closest to) which row. We can go through the stack and calculate the
    // distance squared for each element of the column.
    typename std::vector<Influence>::iterator it = _stack.begin();
    int pos = 0;
    for(current = 0.0; current < w; ++current, ++id, ++pos)
    {
        while( current >= it->right) 
            ++it; 
       
        if( it->prevVector[dimension] < 1.0 ) {
            DestType oldValue = da(id);
           
            
            double thisWidth = it->center - current; //FIXME: take sigma into account
            
            #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG 
            std::cout << "      current=" << current << ", it=" << *it << std::endl;
            std::cout << "        thisWidth=" << thisWidth << std::endl;
            std::cout << "        oldValue=" << oldValue << std::endl;
            #endif
          
            DestType newValue = it->prevVector;
            newValue[dimension] = std::min(thisWidth, oldValue[dimension]);
          
            /*
            if(dimension > 0 && std::fabs(oldValue[dimension-1]) < 1.0) {
                if(it->boundary) {
                    newValue[dimension-1]= 0.0;
                }
                else if(it->corner2) {
                    newValue[dimension-1] = 0.5;
                }
                else if(it->corner1) {
                    newValue[dimension-1] = -0.5;
                }
            }
            */
            
            #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG 
            std::cout << "        newValue=" << newValue << std::endl;
            #endif
           
            //if( partialSquaredMagnitude(newValue, dimension+1) < partialSquaredMagnitude(oldValue, dimension+1) ) {
            da.set(newValue, id);
            //}

        }
        #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG 
        std::cout << "  ==> " << da(id) << std::endl;
        #endif
    }
    #ifdef VECTORIAL_BOUNDARY_DIST_DEBUG 
    std::cout << "<<<" << std::endl;
    #endif
}
} /* namespace detail */

/********************************************************/
/*                                                      */
/*      separableMultiVectorialBoundaryDist             */
/*                                                      */
/********************************************************/

template<int N, class T>
std::pair< vigra::MultiArray<N,  vigra::TinyVector<double, N> >, vigra::MultiArray<N, unsigned char> >
separableMultiVectorialBoundaryDist(
    const vigra::MultiArrayView<N, T>& segmentation
)
{
    using namespace vigra::functor;
  
    //input: a segmenation
    typedef vigra::MultiArrayView<N, T> SegmentationType;
    typedef typename SegmentationType::const_traverser SegmentationTraverser;
    typedef MultiArrayNavigator<SegmentationTraverser, N> SegmentationNavigator;
   
    //output: a vector image/volume
    typedef vigra::TinyVector<double, N> VectorType;
    typedef vigra::MultiArrayView<N, VectorType> ResultType;
    typedef typename ResultType::traverser ResultTraverser;
    typedef MultiArrayNavigator<ResultTraverser,       N> ResultNavigator;
    
    //corners
    typedef vigra::MultiArray<N, unsigned char> CornersType;
    typedef typename CornersType::traverser CornersTraverser;
    typedef MultiArrayNavigator<CornersTraverser, N> CornersNavigator;
    
    vigra::MultiArray<N, VectorType> result(segmentation.shape());
   
    //figure out the maximal value that the components of the distance vectors
    //can obtain --> this number serves as a substitute for "infinity"
    typename VectorType::value_type dmax = 0.0;
    for( int k=0; k<N; ++k)
    {
        dmax += segmentation.shape(k);
    }
   
    //initialize the result array with "infinity" vectors
    VectorType maxDist = VectorType(std::ceil(dmax));
    initMultiArray(destMultiArrayRange(result), maxDist);
    
    //array with 'corners' marked 
    typename MultiArrayShape<N>::type cShape;
    for(int d=0; d<N; ++d) {
        cShape[d] = segmentation.shape(d)+2;
    }
    CornersType corners(cShape);
    for(MultiArrayIndex i=0; i<corners.shape(0); ++i) {
        for(MultiArrayIndex j=0; j<corners.shape(1); ++j) {
            T a1 = 0;
            T a2 = 0;
            T a3 = 0;
            T a4 = 0;
            
            MultiArrayIndex I,J;
            
            I = i;
            J = j;
            if( I >= 0 && J >= 0 && I < segmentation.shape(0) && J < segmentation.shape(1) ) {
                a1 = segmentation(I, J);
            }
            I = i;
            J = j-1;
            if( I >= 0 && J >= 0 && I < segmentation.shape(0) && J < segmentation.shape(1) ) {
                a2 = segmentation(I, J);
            }
            I = i-1;
            J = j-1;
            if( I >= 0 && J >= 0 && I < segmentation.shape(0) && J < segmentation.shape(1) ) {
                a3 = segmentation(I, J);
            }
            I = i-1;
            J = j;
            if( I >= 0 && J >= 0 && I < segmentation.shape(0) && J < segmentation.shape(1) ) {
                a4 = segmentation(I, J);
            }
            
            bool isCorner = true;
           
            if(a1 == a2 && a3 == a4) {
                isCorner = false;
            }
            if(a2 == a3 && a1 == a4) {
                isCorner = false;
            }
           
            if(isCorner) {
                corners(i,j) = 1;
            }
        }
    }

    // temporary arrays to hold the current line to enable in-place operation
    ArrayVector<VectorType> tmp;
    ArrayVector<typename SegmentationType::value_type> segmentationLine;
    ArrayVector<unsigned char> cornersLine1, cornersLine2;
    
    for( int d = 0; d < N; ++d )
    {
        ResultNavigator       resultNav      ( result.traverser_begin(),       result.shape(),       d );
        SegmentationNavigator segmentationNav( segmentation.traverser_begin(), segmentation.shape(), d );
        CornersNavigator      cornersNav(       corners.traverser_begin(),     corners.shape(),      d );
        
        tmp.resize( segmentation.shape(d) );
        segmentationLine.resize( segmentation.shape(d) );
        cornersLine1.resize( corners.shape(d) );
        cornersLine2.resize( corners.shape(d) );
        
        for( ; segmentationNav.hasMore(); resultNav++, segmentationNav++, cornersNav++ )
        {
                // first copy source to temp for maximum cache efficiency
                std::copy(resultNav.begin(), resultNav.end(), tmp.begin());
                std::copy(segmentationNav.begin(), segmentationNav.end(), segmentationLine.begin());
                std::copy(cornersNav.begin(), cornersNav.end(), cornersLine1.begin());
                if(segmentationNav.hasMore()) {
                    ++cornersNav;
                    std::copy(cornersNav.begin(), cornersNav.end(), cornersLine2.begin());
                    --cornersNav;
                }

                detail::vectorialBoundaryDistParabola(d /*dimension*/,
                    dmax,
                    tmp.begin(), tmp.end(), typename AccessorTraits<typename ResultType::value_type>::default_const_accessor(),
                    resultNav.begin(), typename AccessorTraits<typename ResultType::value_type>::default_accessor(),
                    segmentationLine.begin(),
                    cornersLine1.begin(),
                    cornersLine2.begin() 
                );
        }
        
        if(d == 0) {
            break;
        }
    }
    
    return std::make_pair(result, corners);
}

} //-- namespace vigra


#endif //-- VIGRA_VECTORIAL_BOUNDARY_DISTANCE_HXX
