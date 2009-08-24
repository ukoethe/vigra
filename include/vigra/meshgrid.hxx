/************************************************************************/
/*                                                                      */
/*                  Copyright 2009 by Ullrich Koethe                    */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#ifndef VIGRA_MESHGRID_HXX
#define VIGRA_MESHGRID_HXX

#include "tinyvector.hxx"
#include "diff2d.hxx"

namespace vigra{
/** \addtogroup RangesAndPoints Two-dimensional Ranges and Points

    Specify a 2D position, extent, or rectangle.
*/
//@{

/********************************************************/
/*                                                      */
/*                 MeshGridAccessor                     */
/*                                                      */
/********************************************************/
/** Accessor for turning iteration over Diff2D into a mesh grid.

    The mesh grid concept is adapted from MATLAB. It is a two banded image
    (i.e. with 2D vector pixel type) whose first band contains the x-coordinates
    of the current pixel, and whose second band contains the y-coordinates.
    See \ref meshGrid() for more detailed documentation.
*/
struct MeshGridAccessor
{
        /** the value_type of a mesh grid is a 2D vector
        */
    typedef TinyVector<Diff2D::MoveX, 2> value_type;

        /** read the current data item
        */
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const
    {
        return value_type(i->x, i->y);
    }

        /** read the data item at an offset (can be 1D or 2D or higher order difference).
        */
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE const & diff) const
    {
        return value_type(i->x+diff.x, i->y+diff.y);
    }
};

/** Create a mesh grid for the specified rectangle.

    The mesh grid concept is adapted from MATLAB. It is a two banded image
    (i.e. with 2D vector pixel type) whose first band contains the x-coordinates
    of the current pixel, and whose second band contains the y-coordinates.
    If \a upperLeft is not the point (0,0), the mesh grid is translated relative to
    the pixel indices.

    <b> Declarations:</b>

    \code
    triple<Diff2D, Diff2D, MeshGridAccessor>
    meshGrid(Diff2D upperLeft, Diff2D lowerRight);

    triple<Diff2D, Diff2D, MeshGridAccessor>
    meshGrid(Rect2D const & r);

    \endcode

    <b>Usage:</b>

    \code
    #include <vigra/meshgrid.hxx>
    // create an image whose values are equal to each pixel's distance from the image center
    int width = 5, height = 7;
    int xc = width/2, yc = height/2; // the image center

    FImage dist(width, height);
    Point2D upperLeft(-xc, -yc);

    using namespace vigra::functor;
    transformImage(meshGrid(upperLeft, upperLeft+dist.size()),
                   destImage(dist),
                   norm(Arg1()));
    \endcode
*/
inline
triple<Diff2D, Diff2D, MeshGridAccessor>
meshGrid(Diff2D upperLeft, Diff2D lowerRight)
{
    return triple<Diff2D, Diff2D, MeshGridAccessor>(upperLeft, lowerRight, MeshGridAccessor());
}

inline
triple<Diff2D, Diff2D, MeshGridAccessor>
meshGrid(Rect2D const & r)
{
    return triple<Diff2D, Diff2D, MeshGridAccessor>(r.upperLeft(), r.lowerRight(), MeshGridAccessor());
}

}//namespace vigra
//@}
#endif //VIGRA_MESHGRID_HXX
