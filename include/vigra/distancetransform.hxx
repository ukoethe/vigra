/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/
 
 
#ifndef VIGRA_DISTANCETRANSFORM_HXX
#define VIGRA_DISTANCETRANSFORM_HXX

#include <cmath>
#include "vigra/stdimage.hxx"

namespace vigra {

/*
 * functors to determine the distance norm 
 * these functors assume that dx and dy are positive
 * (this is OK for use in internalDistanceTransform())
 */
 
// chessboard metric
struct InternalDistanceTransformLInifinityNormFunctor
{
    float operator()(float dx, float dy) const
    {
        return (dx < dy) ? dy : dx;
    }
};

// Manhattan metric
struct InternalDistanceTransformL1NormFunctor
{
    float operator()(float dx, float dy) const
    {
        return dx + dy;
    }
};

// Euclidean metric
struct InternalDistanceTransformL2NormFunctor
{
    float operator()(float dx, float dy) const
    {
        return std::sqrt(dx*dx + dy*dy);
    }
};


template <class SrcImageIterator, class SrcAccessor,
                   class DestImageIterator, class DestAccessor,
                   class ValueType, class NormFunctor>
void
internalDistanceTransform(SrcImageIterator src_upperleft, 
                SrcImageIterator src_lowerright, SrcAccessor sa,
                DestImageIterator dest_upperleft, DestAccessor da,
                ValueType background, NormFunctor norm)
{
    int w = src_lowerright.x - src_upperleft.x;  
    int h = src_lowerright.y - src_upperleft.y;  
    
    FImage xdist(w,h), ydist(w,h);
    
    xdist = w;    // init x and
    ydist = h;    // y distances with 'large' values

    SrcImageIterator sy = src_upperleft;
    DestImageIterator ry = dest_upperleft;
    FImage::Iterator xdy = xdist.upperLeft();
    FImage::Iterator ydy = ydist.upperLeft();
    SrcImageIterator sx = sy;
    DestImageIterator rx = ry;
    FImage::Iterator xdx = xdy;
    FImage::Iterator ydx = ydy;
    
    static const Diff2D left(-1, 0);
    static const Diff2D right(1, 0);
    static const Diff2D top(0, -1);
    static const Diff2D bottom(0, 1);
        
    int x,y;
    if(sa(sx) != background)    // first pixel
    {
        *xdx = 0.0;
        *ydx = 0.0;
        da.set(0.0, rx);
    }
    else
    {
        da.set(norm(*xdx, *ydx), rx);
    }
    
    
    for(x=1, ++xdx.x, ++ydx.x, ++sx.x, ++rx.x; 
        x<w; 
        ++x, ++xdx.x, ++ydx.x, ++sx.x, ++rx.x)   // first row left to right
    {
        if(sa(sx) != background)
        {
            *xdx = 0.0;
            *ydx = 0.0;
            da.set(0.0, rx);
        }
        else
        {
            *xdx  = xdx[left] + 1.0;   // propagate x and
            *ydx  = ydx[left];         // y components of distance from left pixel
            da.set(norm(*xdx, *ydx), rx); // calculate distance from x and y components
        }
    }
    for(x=w-2, xdx.x -= 2, ydx.x -= 2, sx.x -= 2, rx.x -= 2; 
        x>=0; 
        --x, --xdx.x, --ydx.x, --sx.x, --rx.x)   // first row right to left
    {
        float d = norm(xdx[right] + 1.0, ydx[right]);
        
        if(da(rx) < d) continue;
        
        *xdx = xdx[right] + 1.0;
        *ydx = ydx[right];
        da.set(d, rx);
    }
    for(y=1, ++xdy.y, ++ydy.y, ++sy.y, ++ry.y; 
        y<h;
        ++y, ++xdy.y, ++ydy.y, ++sy.y, ++ry.y)   // top to bottom
    {
        sx = sy;
        rx = ry;
        xdx = xdy;
        ydx = ydy;
        
        if(sa(sx) != background)    // first pixel of current row
        {
            *xdx = 0.0;
            *ydx = 0.0;
            da.set(0.0, rx);
        }
        else
        {
            *xdx = xdx[top];
            *ydx = ydx[top] + 1.0;
            da.set(norm(*xdx, *ydx), rx);
        }
        
        for(x=1, ++xdx.x, ++ydx.x, ++sx.x, ++rx.x; 
            x<w; 
            ++x, ++xdx.x, ++ydx.x, ++sx.x, ++rx.x)  // current row left to right
        {
            if(sa(sx) != background)
            {
                *xdx = 0.0;
                *ydx = 0.0;
                da.set(0.0, rx);
            }
            else
            {
                float d1 = norm(xdx[left] + 1.0, ydx[left]);
                float d2 = norm(xdx[top], ydx[top] + 1.0);
                
                if(d1 < d2)
                {
                    *xdx = xdx[left] + 1.0;
                    *ydx = ydx[left];
                    da.set(d1, rx);
                }
                else
                {
                    *xdx = xdx[top];
                    *ydx = ydx[top] + 1.0;
                    da.set(d2, rx);
                }
            }
        }
        for(x=w-2, xdx.x -= 2, ydx.x -= 2, sx.x -= 2, rx.x -= 2; 
            x>=0; 
            --x, --xdx.x, --ydx.x, --sx.x, --rx.x)  // current row right to left
        {
            float d1 = norm(xdx[right] + 1.0, ydx[right]);
            
            if(da(rx) < d1) continue;
            
            *xdx = xdx[right] + 1.0;
            *ydx = ydx[right];
            da.set(d1, rx);
        }
    }
    for(y=h-2, xdy.y -= 2, ydy.y -= 2, sy.y -= 2, ry.y -= 2; 
        y>=0;
        --y, --xdy.y, --ydy.y, --sy.y, --ry.y)    // bottom to top
    {
        sx = sy;
        rx = ry;
        xdx = xdy;
        ydx = ydy;
        
        float d = norm(xdx[bottom], ydx[bottom] + 1.0);
        if(d < da(rx))    // first pixel of current row
        { 
            *xdx = xdx[bottom];
            *ydx = ydx[bottom] + 1.0;
            da.set(d, rx);
        }
            
        for(x=1, ++xdx.x, ++ydx.x, ++sx.x, ++rx.x; 
            x<w;
            ++x, ++xdx.x, ++ydx.x, ++sx.x, ++rx.x)  // current row left to right
        {
            float d1 = norm(xdx[left] + 1.0, ydx[left]);
            float d2 = norm(xdx[bottom], ydx[bottom] + 1.0);
            
            if(d1 < d2)
            {
                if(da(rx) < d1) continue;
                *xdx = xdx[left] + 1.0;
                *ydx = ydx[left];
                da.set(d1, rx);
            }
            else
            {
                if(da(rx) < d2) continue;
                *xdx = xdx[bottom];
                *ydx = ydx[bottom] + 1.0;
                da.set(d2, rx);
            }
        }
        for(x=w-2, xdx.x -= 2, ydx.x -= 2, sx.x -= 2, rx.x -= 2; 
            x>=0; 
            --x, --xdx.x, --ydx.x, --sx.x, --rx.x)  // current row right to left
        {
            float d1 = norm(xdx[right] + 1.0, ydx[right]);

            if(da(rx) < d1) continue;
            *xdx = xdx[right] + 1.0;
            *ydx = ydx[right];
            da.set(d1, rx);
        }
    }
}

/********************************************************/
/*                                                      */
/*                 distanceTransform                    */
/*                                                      */
/********************************************************/

/** @name Distance Transform
    @memo using Euclidean, Manhattan, or chessboard metrics
*/
//@{

/** For all background pixels, calculate the distance to 
    the nearest object or contour. The label of the pixels to be considered 
    background in the source image is passed in the parameter 'background'.
    Source pixels with other labels will be considered objects. In the 
    destination image, all pixels corresponding to background will be assigned 
    the their distance value, all pixels corresponding to objects will be
    assigned 0.
    
    The parameter 'norm' gives the distance norm to be used:
    
    \begin{itemize}

    \item norm == 0: use chessboard distance (L-infinity norm)
    \item norm == 1: use Manhattan distance (L1 norm)
    \item norm == 2: use Euclidean distance (L2 norm)
    
    \end{itemize}
    
    If you use the L2 norm, the destination pixels must be real valued to give
    correct results.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                           class DestImageIterator, class DestAccessor,
                           class ValueType>
        void distanceTransform(SrcImageIterator src_upperleft, 
                        SrcImageIterator src_lowerright, SrcAccessor sa,
                        DestImageIterator dest_upperleft, DestAccessor da,
                        ValueType background, int norm)
    }
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                           class DestImageIterator, class DestAccessor,
                           class ValueType>
        void distanceTransform(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest,
            ValueType background, int norm)
    }
    \end{verbatim}
    
    {\bf Usage:}
    
    Include-File:
    \URL[vigra/distancetransform.hxx]{../include/vigra/distancetransform.hxx}\\
    Namespace: vigra
    
    
    \begin{verbatim}
    
    vigra::BImage src(w,h), edges(w,h);
    vigra::FImage distance(w, h);

    // empty edge image
    edges = 0;
    ...

    // detect edges in src image (edges will be marked 1, background 0)
    vigra::differenceOfExponentialEdgeImage(srcImageRange(src), destImage(edges), 
                                     0.8, 4.0);
     
    // find distance of all pixels from nearest edge
    vigra::distanceTransform(srcImageRange(edges), destImage(distance),
                             0,                   2);
    //                       ^ background label   ^ norm (Euclidean)
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor sa;
    DestAccessor da;
    
    ValueType background;
    float distance;
    
    sa(src_upperleft) != background;
    da(dest_upperleft) < distance;
    da.set(distance, dest_upperleft);
 
    \end{verbatim}
    

    @memo
*/
template <class SrcImageIterator, class SrcAccessor,
                   class DestImageIterator, class DestAccessor,
                   class ValueType>
inline void
distanceTransform(SrcImageIterator src_upperleft, 
                SrcImageIterator src_lowerright, SrcAccessor sa,
                DestImageIterator dest_upperleft, DestAccessor da,
                ValueType background, int norm)
{
    if(norm == 1)
    {
        internalDistanceTransform(src_upperleft, src_lowerright, sa,
                                  dest_upperleft, da, background,
                                  InternalDistanceTransformL1NormFunctor());
    }
    else if(norm == 2)
    {
        internalDistanceTransform(src_upperleft, src_lowerright, sa,
                                  dest_upperleft, da, background,
                                  InternalDistanceTransformL2NormFunctor());
    }
    else
    {
        internalDistanceTransform(src_upperleft, src_lowerright, sa,
                                  dest_upperleft, da, background,
                                  InternalDistanceTransformLInifinityNormFunctor());
    }
}

template <class SrcImageIterator, class SrcAccessor,
                   class DestImageIterator, class DestAccessor,
                   class ValueType>
inline void
distanceTransform(
    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
    pair<DestImageIterator, DestAccessor> dest,
    ValueType background, int norm)
{
    distanceTransform(src.first, src.second, src.third,
                      dest.first, dest.second, background, norm);
}

//@}

} // namespace vigra

#endif // VIGRA_DISTANCETRANSFORM_HXX

