/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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


#ifndef VIGRA_LABELIMAGE_HXX
#define VIGRA_LABELIMAGE_HXX

#include <vector>
#include <functional>
#include "utilities.hxx"
#include "stdimage.hxx"
#include "union_find.hxx"
#include "sized_int.hxx"
#include "multi_shape.hxx"

namespace vigra {

/** \addtogroup Labeling Connected Components Labeling
     The 2-dimensional connected components algorithms may use either 4 or 8 connectivity.
     By means of a functor the merge criterion can be defined arbitrarily.
*/
//@{

/********************************************************/
/*                                                      */
/*                        labelImage                    */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a segmented image.

    Connected components are defined as regions with uniform pixel
    values. Thus, <TT>T1</TT> either must be
    equality comparable, or a suitable EqualityFunctor must be
    provided that realizes the desired predicate. The
    destination's value type <tt>T2</tt> should be large enough to hold the labels
    without overflow. Region numbers will be a consecutive sequence
    starting with one and ending with the region number returned by
    the function (inclusive). The parameter '<TT>eight_neighbors</TT>'
    determines whether the regions should be 4-connected or
    8-connected. 

    Return:  the number of regions found (= largest region label)
    
    See \ref labelMultiArray() for a dimension-independent implementation of 
    connected components labelling.
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class EqualityFunctor = std::equal_to<T1> >
        unsigned int
        labelImage(MultiArrayView<2, T1, S1> const & src,
                   MultiArrayView<2, T2, S2> dest,
                   bool eight_neighbors, EqualityFunctor equal = EqualityFunctor());
    }
    \endcode

    \deprecatedAPI{labelImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        unsigned int labelImage(SrcIterator upperlefts,
                                SrcIterator lowerrights, SrcAccessor sa,
                                DestIterator upperleftd, DestAccessor da,
                                bool eight_neighbors);

        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class EqualityFunctor>
        unsigned int labelImage(SrcIterator upperlefts,
                                SrcIterator lowerrights, SrcAccessor sa,
                                DestIterator upperleftd, DestAccessor da,
                                bool eight_neighbors, EqualityFunctor equal);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        unsigned int labelImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                pair<DestIterator, DestAccessor> dest,
                                bool eight_neighbors);

        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class EqualityFunctor>
        unsigned int labelImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                pair<DestIterator, DestAccessor> dest,
                                bool eight_neighbors, EqualityFunctor equal)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/labelimage.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, unsigned char> src(w,h);
    MultiArray<2, unsigned int>  labels(w,h);

    // threshold at 128
    transformImage(src, src, Threshold<int, int>(128, 256, 0, 255));

    // find 4-connected regions
    labelImage(src, labels, false);
    \endcode

    \deprecatedUsage{labelImage}
    \code
    vigra::BImage src(w,h);
    vigra::IImage labels(w,h);

    // threshold at 128
    vigra::transformImage(srcImageRange(src), destImage(src),
       vigra::Threshold<vigra::BImage::PixelType, vigra::BImage::PixelType>(
                                                    128, 256, 0, 255));

    // find 4-connected regions
    vigra::labelImage(srcImageRange(src), destImage(labels), false);
    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u == u                  // first form

    EqualityFunctor equal;      // second form
    equal(u, u)                 // second form

    int i;
    dest_accessor.set(i, dest_upperleft);
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> unsigned int labelImage)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class EqualityFunctor>
unsigned int labelImage(SrcIterator upperlefts,
                        SrcIterator lowerrights, SrcAccessor sa,
                        DestIterator upperleftd, DestAccessor da,
                        bool eight_neighbors, EqualityFunctor equal)
{
    typedef typename DestAccessor::value_type LabelType;
    
    int w = lowerrights.x - upperlefts.x;
    int h = lowerrights.y - upperlefts.y;
    int x,y,i;

    const Diff2D neighbor[] = {
        Diff2D(-1,0),  // left
        Diff2D(-1,-1), // topleft
        Diff2D(0,-1),  // top
        Diff2D(1,-1)   // topright
    };

    const int left = 0, /* unused:  topleft = 1, */ top = 2, topright = 3;
    int step = eight_neighbors ? 1 : 2;

    SrcIterator ys = upperlefts;
    DestIterator yd = upperleftd;
    
    detail::UnionFindArray<LabelType>  label;    

    // pass 1: scan image from upper left to lower right
    // to find connected components

    // Each component will be represented by a tree of pixels. Each
    // pixel contains the scan order address of its parent in the
    // tree.  In order for pass 2 to work correctly, the parent must
    // always have a smaller scan order address than the child.
    // Therefore, we can merge trees only at their roots, because the
    // root of the combined tree must have the smallest scan order
    // address among all the tree's pixels/ nodes.  The root of each
    // tree is distinguished by pointing to itself (it contains its
    // own scan order address). This condition is enforced whenever a
    // new region is found or two regions are merged


    for(y = 0; y != h; ++y, ++ys.y, ++yd.y)
    {
        SrcIterator xs = ys;
        DestIterator xd = yd;

        int endNeighbor = (y == 0) ? left : (eight_neighbors ? topright : top);

        for(x = 0; x != w; ++x, ++xs.x, ++xd.x)
        {
            int beginNeighbor = (x == 0) ? top : left;
            if(x == w-1 && endNeighbor == topright) endNeighbor = top;

            for(i=beginNeighbor; i<=endNeighbor; i+=step)
            {
                if(equal(sa(xs), sa(xs, neighbor[i])))
                {
                    LabelType neighborLabel = label.find(da(xd,neighbor[i]));

                    for(int j=i+2; j<=endNeighbor; j+=step)
                    {
                        if(equal(sa(xs), sa(xs, neighbor[j])))
                        {
                            neighborLabel = label.makeUnion(da(xd, neighbor[j]), neighborLabel);
                            break;
                        }
                    }
                    da.set(neighborLabel, xd);
                    break;
                }

            }
            if(i > endNeighbor)
            {
                da.set(label.makeNewLabel(), xd);
            }
        }
    }

    // pass 2: assign one label to each region (tree)
    // so that labels form a consecutive sequence 1, 2, ...
    unsigned int count = label.makeContiguous();    
    
    yd = upperleftd;
    for(y=0; y != h; ++y, ++yd.y)
    {
        typename DestIterator::row_iterator xd = yd.rowIterator();
        for(x = 0; x != w; ++x, ++xd)
        {
            da.set(label[da(xd)], xd);
        }
    }
    return count;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
unsigned int labelImage(SrcIterator upperlefts,
                        SrcIterator lowerrights, SrcAccessor sa,
                        DestIterator upperleftd, DestAccessor da,
                        bool eight_neighbors)
{
    return labelImage(upperlefts, lowerrights, sa,
                 upperleftd, da, eight_neighbors,
                 std::equal_to<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class EqualityFunctor>
inline unsigned int
labelImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
           pair<DestIterator, DestAccessor> dest,
           bool eight_neighbors, EqualityFunctor equal)
{
    return labelImage(src.first, src.second, src.third,
                      dest.first, dest.second, eight_neighbors, equal);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline unsigned int
labelImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
           pair<DestIterator, DestAccessor> dest,
           bool eight_neighbors)
{
    return labelImage(src.first, src.second, src.third,
                      dest.first, dest.second, eight_neighbors,
                      std::equal_to<typename SrcAccessor::value_type>());
}

template <class T1, class S1,
          class T2, class S2,
          class EqualityFunctor>
inline unsigned int
labelImage(MultiArrayView<2, T1, S1> const & src,
           MultiArrayView<2, T2, S2> dest,
           bool eight_neighbors, EqualityFunctor equal)
{
    vigra_precondition(src.shape() == dest.shape(),
        "labelImage(): shape mismatch between input and output.");
    return labelImage(srcImageRange(src),
                      destImage(dest), eight_neighbors, equal);
}

template <class T1, class S1,
          class T2, class S2>
inline unsigned int
labelImage(MultiArrayView<2, T1, S1> const & src,
           MultiArrayView<2, T2, S2> dest,
           bool eight_neighbors)
{
    return labelImage(srcImageRange(src),
                      destImage(dest), eight_neighbors,
                      std::equal_to<T1>());
}

/********************************************************/
/*                                                      */
/*             labelImageWithBackground                 */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a segmented image,
    excluding the background from labeling.

    This function works like \ref labelImage(), but considers all background pixels
    (i.e. pixels having the given '<TT>background_value</TT>') as a single region that
    is ignored when determining connected components and remains untouched in the
    destination image. Usually, you will zero-initialize the output image, so that
    the background gets label 0 (remember that actual region labels start at one).

    Return:  the number of non-background regions found (= largest region label)
    
    See \ref labelMultiArrayWithBackground() for a dimension-independent implementation
    if this algorithm.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class ValueType, 
                  class EqualityFunctor = std::equal_to<T1> >
        unsigned int 
        labelImageWithBackground(MultiArrayView<2, T1, S1> const & src,
                                 MultiArrayView<2, T2, S2> dest,
                                 bool eight_neighbors,
                                 ValueType background_value, 
                                 EqualityFunctor equal = EqualityFunctor());
    }
    \endcode

    \deprecatedAPI{labelImageWithBackground}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class ValueType>
        int labelImageWithBackground(SrcIterator upperlefts,
                       SrcIterator lowerrights, SrcAccessor sa,
                       DestIterator upperleftd, DestAccessor da,
                       bool eight_neighbors,
                       ValueType background_value );

        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class ValueType, class EqualityFunctor>
        int labelImageWithBackground(SrcIterator upperlefts,
                       SrcIterator lowerrights, SrcAccessor sa,
                       DestIterator upperleftd, DestAccessor da,
                       bool eight_neighbors,
                       ValueType background_value, EqualityFunctor equal);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class ValueType>
        int labelImageWithBackground(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                     pair<DestIterator, DestAccessor> dest,
                                     bool eight_neighbors,
                                     ValueType background_value);

        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class ValueType, class EqualityFunctor>
        int labelImageWithBackground(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                     pair<DestIterator, DestAccessor> dest,
                                     bool eight_neighbors,
                                     ValueType background_value, EqualityFunctor equal);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/labelimage.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, unsigned char> src(w,h);
    MultiArray<2, unsigned int>  labels(w,h);

    // threshold at 128
    transformImage(src, src, Threshold<int, int>(128, 256, 0, 255));

    // find 4-connected regions of foreground (= white pixels) only
    labelImageWithBackground(src, labels, false, 0);
    \endcode

    \deprecatedUsage{labelImageWithBackground}
    \code
    vigra::BImage src(w,h);
    vigra::IImage labels(w,h);

    // threshold at 128
    vigra::transformImage(srcImageRange(src), destImage(src),
        vigra::Threshold<vigra::BImage::PixelType, vigra::BImage::PixelType>(
                                                    128, 256, 0, 255));

    // find 4-connected regions of foreground (= white pixels) only
    vigra::labelImageWithBackground(srcImageRange(src), destImage(labels),
                             false, 0);
    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);
    ValueType background_value;

    u == u                  // first form
    u == background_value   // first form

    EqualityFunctor equal;      // second form
    equal(u, u)                 // second form
    equal(u, background_value)  // second form

    int i;
    dest_accessor.set(i, dest_upperleft);
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> unsigned int labelImageWithBackground)
    
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class ValueType, class EqualityFunctor>
unsigned int labelImageWithBackground(
    SrcIterator upperlefts,
    SrcIterator lowerrights, SrcAccessor sa,
    DestIterator upperleftd, DestAccessor da,
    bool eight_neighbors,
    ValueType background_value, EqualityFunctor equal)
{
    int w = lowerrights.x - upperlefts.x;
    int h = lowerrights.y - upperlefts.y;
    int x,y,i;

    const Diff2D neighbor[] = {
        Diff2D(-1,0),  // left
        Diff2D(-1,-1), // topleft
        Diff2D(0,-1),  // top
        Diff2D(1,-1)   // topright
    };

    const int left = 0, /* unused:  topleft = 1,*/ top = 2, topright = 3;
    int step = eight_neighbors ? 1 : 2;

    SrcIterator ys(upperlefts);
    SrcIterator xs(ys);
    
    // temporary image to store region labels
    typedef BasicImage<IntBiggest> TmpImage;
    TmpImage labelimage(w, h);
    TmpImage::ScanOrderIterator label = labelimage.begin();
    TmpImage::Iterator yt = labelimage.upperLeft();
    TmpImage::Iterator  xt(yt);

    // pass 1: scan image from upper left to lower right
    // find connected components

    for(y = 0; y != h; ++y, ++ys.y, ++yt.y)
    {
        xs = ys;
        xt = yt;

        int endNeighbor = (y == 0) ? left : (eight_neighbors ? topright : top);

        for(x = 0; x != w; ++x, ++xs.x, ++xt.x)
        {
            if(equal(sa(xs), background_value))
            {
                *xt = -1;
            }
            else
            {
                int beginNeighbor = (x == 0) ? top : left;
                if(x == w-1 && endNeighbor == topright) endNeighbor = top;

                for(i=beginNeighbor; i<=endNeighbor; i+=step)
                {
                    if(equal(sa(xs), sa(xs, neighbor[i])))
                    {
                        IntBiggest neighborLabel = xt[neighbor[i]];

                        for(int j=i+2; j<=endNeighbor; j+=step)
                        {
                            if(equal(sa(xs), sa(xs, neighbor[j])))
                            {
                                IntBiggest neighborLabel1 = xt[neighbor[j]];

                                if(neighborLabel != neighborLabel1)
                                {
                                    // find roots of the region trees
                                    while(neighborLabel != label[neighborLabel])
                                    {
                                        neighborLabel = label[neighborLabel];
                                    }
                                    while(neighborLabel1 != label[neighborLabel1])
                                    {
                                        neighborLabel1 = label[neighborLabel1];
                                    }

                                    // merge the trees
                                    if(neighborLabel1 < neighborLabel)
                                    {
                                        label[neighborLabel] = neighborLabel1;
                                        neighborLabel = neighborLabel1;
                                    }
                                    else if(neighborLabel < neighborLabel1)
                                    {
                                        label[neighborLabel1] = neighborLabel;
                                    }
                                }
                                break;
                            }
                        }
                        *xt = neighborLabel;
                        break;
                    }

                }
                if(i > endNeighbor)
                {
                    // new region
                    // The initial label of a new region equals the
                    // scan order address of it's first pixel.
                    // This is essential for correct operation of the algorithm.
                    *xt = x + y*w;
                }
            }
        }
    }

    // pass 2: assign contiguous labels to the regions
    DestIterator yd(upperleftd);

    int count = 0;
    i = 0;
    for(y=0; y != h; ++y, ++yd.y)
    {
        DestIterator xd(yd);
        for(x = 0; x != w; ++x, ++xd.x, ++i)
        {
            if(label[i] == -1) continue;

            if(label[i] == i)
            {
                label[i] = count++;
            }
            else
            {
                label[i] = label[label[i]];
            }
            da.set(label[i]+1, xd);
        }
    }

    return count;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class ValueType>
inline
unsigned int labelImageWithBackground(
    SrcIterator upperlefts,
    SrcIterator lowerrights, SrcAccessor sa,
    DestIterator upperleftd, DestAccessor da,
    bool eight_neighbors,
    ValueType background_value)
{
    return labelImageWithBackground(upperlefts, lowerrights, sa,
                            upperleftd, da,
                            eight_neighbors, background_value,
                            std::equal_to<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class ValueType, class EqualityFunctor>
inline unsigned int 
labelImageWithBackground(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest,
                         bool eight_neighbors,
                         ValueType background_value, EqualityFunctor equal)
{
    return labelImageWithBackground(src.first, src.second, src.third,
                                    dest.first, dest.second,
                                    eight_neighbors, background_value, equal);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class ValueType>
inline unsigned int
labelImageWithBackground(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest,
                         bool eight_neighbors,
                         ValueType background_value)
{
    return labelImageWithBackground(src.first, src.second, src.third,
                                    dest.first, dest.second,
                                    eight_neighbors, background_value,
                                    std::equal_to<typename SrcAccessor::value_type>());
}

template <class T1, class S1,
          class T2, class S2,
          class ValueType, class EqualityFunctor>
inline unsigned int 
labelImageWithBackground(MultiArrayView<2, T1, S1> const & src,
                         MultiArrayView<2, T2, S2> dest,
                         bool eight_neighbors,
                         ValueType background_value, EqualityFunctor equal)
{
    vigra_precondition(src.shape() == dest.shape(),
        "labelImageWithBackground(): shape mismatch between input and output.");
    return labelImageWithBackground(srcImageRange(src),
                                    destImage(dest),
                                    eight_neighbors, background_value, equal);
}

template <class T1, class S1,
          class T2, class S2,
          class ValueType>
inline unsigned int
labelImageWithBackground(MultiArrayView<2, T1, S1> const & src,
                         MultiArrayView<2, T2, S2> dest,
                         bool eight_neighbors,
                         ValueType background_value)
{
    vigra_precondition(src.shape() == dest.shape(),
        "labelImageWithBackground(): shape mismatch between input and output.");
    return labelImageWithBackground(srcImageRange(src),
                                    destImage(dest),
                                    eight_neighbors, background_value,
                                    std::equal_to<T1>());
}

/********************************************************/
/*                                                      */
/*            regionImageToCrackEdgeImage               */
/*                                                      */
/********************************************************/

/** \brief Transform a labeled image into a crack edge (interpixel edge) image.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2, class DestValue>
        void 
        regionImageToCrackEdgeImage(MultiArrayView<2, T1, S1> const & src,
                                    MultiArrayView<2, T2, S2> dest,
                                    DestValue edge_marker);
    }
    \endcode

    \deprecatedAPI{regionImageToCrackEdgeImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor, class DestValue>
        void regionImageToCrackEdgeImage(
                       SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                       DestIterator dul, DestAccessor da,
                       DestValue edge_marker)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor, class DestValue>
        void regionImageToCrackEdgeImage(
                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest,
                   DestValue edge_marker)
    }
    \endcode
    \deprecatedEnd

    This algorithm inserts border pixels (so called "crack edges" or "interpixel edges")
    between regions in a labeled image like this (<TT>a</TT> and
    <TT>c</TT> are the original labels, and <TT>0</TT> is the value of
    <TT>edge_marker</TT> and denotes the inserted edges):

    \code
       original image     insert zero- and one-cells

                                         a 0 c c c
          a c c                          a 0 0 0 c
          a a c               =>         a a a 0 c
          a a a                          a a a 0 0
                                         a a a a a
    \endcode

    The algorithm assumes that the original labeled image contains
    no background. Therefore, it is suitable as a post-processing
    operation of \ref labelImage() or \ref seededRegionGrowing().

    The destination image must be twice the size of the original
    (precisely, <TT>(2*w-1)</TT> by <TT>(2*h-1)</TT> pixels). The
    source value type (<TT>SrcAccessor::value-type</TT>) must be
    equality-comparable.

    <b> Usage:</b>

    <b>\#include</b> \<vigra/labelimage.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, unsigned char> src(w,h);
    MultiArray<2, unsigned int>  labels(w,h),
                                 cellgrid(2*w-1, 2*h-1);

    // threshold at 128
    transformImage(src, src, Threshold<int, int>(128, 256, 0, 255));

    // find 4-connected regions
    labelImage(src, labels, false);

    // create cell grid image, mark edges with 0
    regionImageToCrackEdgeImage(labels, cellgrid, 0);
    \endcode

    \deprecatedUsage{regionImageToCrackEdgeImage}
    \code
    vigra::BImage src(w,h);
    vigra::IImage labels(w,h);
    vigra::IImage cellgrid(2*w-1, 2*h-1);

    // threshold at 128
    vigra::transformImage(srcImageRange(src), destImage(src),
       vigra::Threshold<vigra::BImage::PixelType, vigra::BImage::PixelType>(
                                                    128, 256, 0, 255));

    // find 4-connected regions
    vigra::labelImage(srcImageRange(src), destImage(labels), false);

    // create cell grid image, mark edges with 0
    vigra::regionImageToCrackEdgeImage(srcImageRange(labels), destImage(cellgrid), 0);
    \endcode
    <b> Required Interface:</b>
    \code
    ImageIterator src_upperleft, src_lowerright;
    ImageIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u != u

    DestValue edge_marker;
    dest_accessor.set(edge_marker, dest_upperleft);
    \endcode
    \deprecatedEnd

    <b> Preconditions:</b>

    The destination image must have twice the size of the source:
    \code
    w_dest = 2 * w_src - 1
    h_dest = 2 * h_src - 1
    \endcode
*/
doxygen_overloaded_function(template <...> void regionImageToCrackEdgeImage)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
void regionImageToCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
               DestIterator dul, DestAccessor da,
               DestValue edge_marker)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    int x,y;

    const Diff2D right(1,0);
    const Diff2D left(-1,0);
    const Diff2D bottomright(1,1);
    const Diff2D bottom(0,1);
    const Diff2D top(0,-1);

    SrcIterator iy = sul;
    DestIterator dy = dul;

    for(y=0; y<h-1; ++y, ++iy.y, dy.y+=2)
    {
        SrcIterator ix = iy;
        DestIterator dx = dy;

        for(x=0; x<w-1; ++x, ++ix.x, dx.x+=2)
        {
            da.set(sa(ix), dx);
            da.set(sa(ix), dx, bottomright);

            if(sa(ix, right) != sa(ix))
            {
                da.set(edge_marker, dx, right);
            }
            else
            {
                da.set(sa(ix), dx, right);
            }
            if(sa(ix, bottom) != sa(ix))
            {
                da.set(edge_marker, dx, bottom);
            }
            else
            {
                da.set(sa(ix), dx, bottom);
            }

        }

        da.set(sa(ix), dx);
        if(sa(ix, bottom) != sa(ix))
        {
            da.set(edge_marker, dx, bottom);
        }
        else
        {
            da.set(sa(ix), dx, bottom);
        }
    }

    SrcIterator ix = iy;
    DestIterator dx = dy;

    for(x=0; x<w-1; ++x, ++ix.x, dx.x+=2)
    {
        da.set(sa(ix), dx);
        if(sa(ix, right) != sa(ix))
        {
            da.set(edge_marker, dx, right);
        }
        else
        {
            da.set(sa(ix), dx, right);
        }
    }
    da.set(sa(ix), dx);

    dy = dul + Diff2D(1,1);

    const Diff2D dist[] = {right, top, left, bottom };
    // find missing 0-cells
    for(y=0; y<h-1; ++y, dy.y+=2)
    {
        DestIterator dx = dy;

        for(x=0; x<w-1; ++x, dx.x+=2)
        {
            int i;
            for(i=0; i<4; ++i)
            {
                if(da(dx, dist[i]) == edge_marker) break;
            }

            if(i < 4) da.set(edge_marker, dx);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void 
regionImageToCrackEdgeImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest,
                            DestValue edge_marker)
{
    regionImageToCrackEdgeImage(src.first, src.second, src.third,
                                dest.first, dest.second,
                                edge_marker);
}

template <class T1, class S1,
          class T2, class S2, class DestValue>
inline void 
regionImageToCrackEdgeImage(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, T2, S2> dest,
                            DestValue edge_marker)
{
    vigra_precondition(2*src.shape()-Shape2(1) == dest.shape(),
        "regionImageToCrackEdgeImage(): shape mismatch between input and output.");
    regionImageToCrackEdgeImage(srcImageRange(src),
                                destImage(dest),
                                edge_marker);
}

/********************************************************/
/*                                                      */
/*                regionImageToEdgeImage                */
/*                                                      */
/********************************************************/

/** \brief Transform a labeled image into an edge image.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2, class DestValue>
        void 
        regionImageToEdgeImage(MultiArrayView<2, T1, S1> const & src,
                               MultiArrayView<2, T2, S2> dest,
                               DestValue edge_marker);
    }
    \endcode

    \deprecatedAPI{regionImageToEdgeImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor, class DestValue>
        void regionImageToEdgeImage(
                       SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                       DestIterator dul, DestAccessor da,
                       DestValue edge_marker)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor, class DestValue>
        void regionImageToEdgeImage(
                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest,
                   DestValue edge_marker)
    }
    \endcode
    \deprecatedEnd

    This algorithm marks all pixels with the given <TT>edge_marker</TT>
    which belong to a different region (label) than their right or lower
    neighbors:

    \code
       original image                     edges
                                 (assuming edge_marker == 1)

          a c c                            1 1 *
          a a c               =>           * 1 1
          a a a                            * * *
    \endcode

    The non-edge pixels of the destination image will not be touched.
    The source value type <TT>T1</TT> must be
    equality-comparable.

    <b> Usage:</b>

    <b>\#include</b> \<vigra/labelimage.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, unsigned char> src(w,h),
                                 edges(w,h);
    MultiArray<2, unsigned int>  labels(w,h);

    edges = 255;  // init background (non-edge) to 255

    // threshold at 128
    transformImage(src, src, Threshold<int, int>(128, 256, 0, 255));

    // find 4-connected regions
    labelImage(src, labels, false);

    // create edge image, mark edges with 0
    regionImageToEdgeImage(labels, edges, 0);
    \endcode

    \deprecatedUsage{regionImageToEdgeImage}
    \code
    vigra::BImage src(w,h);
    vigra::IImage labels(w,h);
    vigra::IImage edges(w, h);
    edges = 255;  // init background (non-edge) to 255

    // threshold at 128
    vigra::transformImage(srcImageRange(src), destImage(src),
      vigra::Threshold<vigra::BImage::PixelType, vigra::BImage::PixelType>(
                                                    128, 256, 0, 255));

    // find 4-connected regions
    vigra::labelImage(srcImageRange(src), destImage(labels), false);

    // create edge image, mark edges with 0
    vigra::regionImageToEdgeImage(srcImageRange(labels), destImage(edges), 0);
    \endcode
    <b> Required Interface:</b>
    \code
    ImageIterator src_upperleft, src_lowerright;
    ImageIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u != u

    DestValue edge_marker;
    dest_accessor.set(edge_marker, dest_upperleft);
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void regionImageToEdgeImage)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
void regionImageToEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
               DestIterator dul, DestAccessor da,
               DestValue edge_marker)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    int x,y;

    const Diff2D right(1,0);
    const Diff2D left(-1,0);
    const Diff2D bottomright(1,1);
    const Diff2D bottom(0,1);
    const Diff2D top(0,-1);

    SrcIterator iy = sul;
    DestIterator dy = dul;

    for(y=0; y<h-1; ++y, ++iy.y, ++dy.y)
    {
        SrcIterator ix = iy;
        DestIterator dx = dy;

        for(x=0; x<w-1; ++x, ++ix.x, ++dx.x)
        {
            if(sa(ix, right) != sa(ix))
            {
                da.set(edge_marker, dx);
            }
            if(sa(ix, bottom) != sa(ix))
            {
                da.set(edge_marker, dx);
            }
        }

        if(sa(ix, bottom) != sa(ix))
        {
            da.set(edge_marker, dx);
        }
    }

    SrcIterator ix = iy;
    DestIterator dx = dy;

    for(x=0; x<w-1; ++x, ++ix.x, ++dx.x)
    {
        if(sa(ix, right) != sa(ix))
        {
            da.set(edge_marker, dx);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
inline void 
regionImageToEdgeImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                       pair<DestIterator, DestAccessor> dest,
                       DestValue edge_marker)
{
    regionImageToEdgeImage(src.first, src.second, src.third,
                           dest.first, dest.second,
                           edge_marker);
}

template <class T1, class S1,
          class T2, class S2, class DestValue>
inline void 
regionImageToEdgeImage(MultiArrayView<2, T1, S1> const & src,
                       MultiArrayView<2, T2, S2> dest,
                       DestValue edge_marker)
{
    vigra_precondition(src.shape() == dest.shape(),
        "regionImageToEdgeImage(): shape mismatch between input and output.");
    regionImageToEdgeImage(srcImageRange(src),
                           destImage(dest),
                           edge_marker);
}

//@}

} // namespace vigra

#endif // VIGRA_LABELIMAGE_HXX
