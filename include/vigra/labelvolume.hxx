/************************************************************************/
/*                                                                      */
/*     Copyright 2006-2007 by F. Heinrich, B. Seppke, Ullrich Koethe    */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#ifndef VIGRA_LABELVOLUME_HXX
#define VIGRA_LABELVOLUME_HXX


#include "vigra/voxelneighborhood.hxx"
#include "vigra/multi_array.hxx"

namespace vigra{

/** \addtogroup Labeling Connected Components Labeling
     The connected components algorithm may use either 6 or 26 connectivity.
     By means of a functor the merge criterium can be defined arbitrarily.
*/
//@{

/********************************************************/
/*                                                      */
/*                        labelVolume                   */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a segmented volume.

    Connected components are defined as regions with uniform voxel
    values. Thus, <TT>SrcAccessor::value_type</TT> either must be
    equality comparable (first form), or an EqualityFunctor must be
    provided that realizes the desired predicate (second form). The
    destination's value type should be large enough to hold the labels
    without overflow. Region numbers will be a consecutive sequence
    starting with one and ending with the region number returned by
    the function (inclusive).

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {

        template <class SrcIterator, class SrcAccessor,class SrcShape,
                  class DestIterator, class DestAccessor,
                  class Neighborhood3D>
        unsigned int labelVolume(SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                                 DestIterator d_Iter, DestAccessor da,
                                 Neighborhood3D neighborhood3D);

        template <class SrcIterator, class SrcAccessor,class SrcShape,
                          class DestIterator, class DestAccessor,
                          class Neighborhood3D, class EqualityFunctor>
        unsigned int labelVolume(SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                                 DestIterator d_Iter, DestAccessor da,
                                 Neighborhood3D neighborhood3D, EqualityFunctor equal);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {

        template <class SrcIterator, class SrcAccessor,class SrcShape,
                  class DestIterator, class DestAccessor,
                  class Neighborhood3D>
        unsigned int labelVolume(triple<SrcIterator, SrcShape, SrcAccessor> src,
                                 pair<DestIterator, DestAccessor> dest,
                                 Neighborhood3D neighborhood3D);

        template <class SrcIterator, class SrcAccessor,class SrcShape,
                 class DestIterator, class DestAccessor,
                 class Neighborhood3D, class EqualityFunctor>
        unsigned int labelVolume(triple<SrcIterator, SrcShape, SrcAccessor> src,
                                 pair<DestIterator, DestAccessor> dest,
                                 Neighborhood3D neighborhood3D, EqualityFunctor equal);

    }
    \endcode
    
    use with 3D-Six-Neighborhood:
    \code
    namespace vigra {    
    
        template <class SrcIterator, class SrcAccessor,class SrcShape,
                  class DestIterator, class DestAccessor>
        unsigned int labelVolumeSix(triple<SrcIterator, SrcShape, SrcAccessor> src,
                                    pair<DestIterator, DestAccessor> dest);
                                    
    }
    \endcode

    Return:  the number of regions found (= largest region label)

    <b> Usage:</b>

    <b>\#include</b> "<a href="labelvolume_8hxx-source.html">vigra/labelvolume.hxx</a>"<br>
    Namespace: vigra

    \code
    typedef vigra::MultiArray<3,int> IntVolume;
    IntVolume src(IntVolume::difference_type(w,h,d));
    IntVolume dest(IntVolume::difference_type(w,h,d));
    
    // find 6-connected regions
    int max_region_label = vigra::labelVolumeSix(srcMultiArrayRange(src), destMultiArray(dest));

    // find 26-connected regions
    int max_region_label = vigra::labelVolume(srcMultiArrayRange(src), destMultiArray(dest), NeighborCode3DTwentySix());
    \endcode

    <b> Required Interface:</b>

    \code
    SrcIterator src_begin;
    SrcShape shape;
    DestIterator dest_begin;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_begin);

    u == u                      // first form

    EqualityFunctor equal;      // second form
    equal(u, u)                 // second form

    int i;
    dest_accessor.set(i, dest_begin);
    \endcode

*/
template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor>
unsigned int labelVolumeSix(triple<SrcIterator, SrcShape, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest)
{
    return labelVolume(src.first, src.second, src.third, dest.first, dest.second, NeighborCode3DSix(), std::equal_to<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D>
unsigned int labelVolume(SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                         DestIterator d_Iter, DestAccessor da,
                         Neighborhood3D neighborhood3D)
{
        return labelVolume(s_Iter, srcShape, sa, d_Iter, da, neighborhood3D, std::equal_to<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D>
unsigned int labelVolume(triple<SrcIterator, SrcShape, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest,
                         Neighborhood3D neighborhood3D)
{
    return labelVolume(src.first, src.second, src.third, dest.first, dest.second, neighborhood3D, std::equal_to<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D, class EqualityFunctor>
unsigned int labelVolume(triple<SrcIterator, SrcShape, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest,
                         Neighborhood3D neighborhood3D, EqualityFunctor equal)
{
    return labelVolume(src.first, src.second, src.third, dest.first, dest.second, neighborhood3D, equal);
}

template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D, class EqualityFunctor>
unsigned int labelVolume(SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                         DestIterator d_Iter, DestAccessor da,
                         Neighborhood3D, EqualityFunctor equal)
{

    //basically needed for iteration and border-checks
    int w = srcShape[0], h = srcShape[1], d = srcShape[2];
    int x,y,z, i;
        
    //declare and define Iterators for all three dims at src
    SrcIterator zs = s_Iter;
    SrcIterator ys(zs);
    SrcIterator xs(ys);
        
    // temporary image to store region labels
    typedef vigra::MultiArray<3,int> LabelVolume;
    typedef LabelVolume::traverser LabelTraverser;
    
    LabelVolume labelvolume(srcShape);
        
    //Declare traversers for all three dims at target
    LabelTraverser zt = labelvolume.traverser_begin();
    LabelTraverser yt(zt);
    LabelTraverser xt(yt);

    // Kovalevsky's clever idea to use
    // image iterator and scan order iterator simultaneously
    // memory order indicates label
    LabelVolume::iterator label = labelvolume.begin();
    
    // initialize the neighborhood traversers

#if 0
    NeighborOffsetTraverser<Neighborhood3D> nc(Neighborhood3D::CausalFirst);
    NeighborOffsetTraverser<Neighborhood3D> nce(Neighborhood3D::CausalLast);
#endif /* #if 0 */

    NeighborOffsetCirculator<Neighborhood3D> nc(Neighborhood3D::CausalFirst);
    NeighborOffsetCirculator<Neighborhood3D> nce(Neighborhood3D::CausalLast);
    ++nce;
    // pass 1: scan image from upper left front to lower right back
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
    i=0;
    for(z = 0; z != d; ++z, ++zs.dim2(), ++zt.dim2())
    {
        ys = zs;
        yt = zt;

        for(y = 0; y != h; ++y, ++ys.dim1(), ++yt.dim1())
        {
            xs = ys;
            xt = yt;

            for(x = 0; x != w; ++x, ++xs.dim0(), ++xt.dim0(), ++i)
            {
                *xt = i; // default: new region    

                //queck whether there is a special borde threatment to be used or not
                AtVolumeBorder atBorder = isAtVolumeBorderCausal(x,y,z,w,h,z);
                    
                //We are not at the border!
                if(atBorder == NotAtBorder)
                {
                    

#if 0
                    nc = NeighborOffsetTraverser<Neighborhood3D>(Neighborhood3D::CausalFirst);
#endif /* #if 0 */

                    nc = NeighborOffsetCirculator<Neighborhood3D>(Neighborhood3D::CausalFirst);
                
                    do
                    {            
                        // if colors are equal
                        if(equal(*xs, xs[*nc]))
                        {
                            int neighborLabel = xt[*nc];

                            // find the root label of a label tree
                            while(neighborLabel != label[neighborLabel])
                            {
                                neighborLabel = label[neighborLabel];
                            }

                            if(neighborLabel < *xt) // always keep the smallest among the possible neighbor labels
                            {
                                label[*xt] = neighborLabel;
                                *xt = neighborLabel;
                            }
                            else
                            {
                                label[neighborLabel] = *xt;
                            }
                        }
                        ++nc;
                    }while(nc!=nce);
                }
                //we are at a border - handle this!!
                else
                {

#if 0
                    nc = NeighborOffsetTraverser<Neighborhood3D>(Neighborhood3D::nearBorderDirectionsCausal(atBorder,0));
#endif /* #if 0 */

                    nc = NeighborOffsetCirculator<Neighborhood3D>(Neighborhood3D::nearBorderDirectionsCausal(atBorder,0));
                    int j=0;
                    while(nc.direction() != Neighborhood3D::Error)
                    {
                        //   colors equal???
                        if(equal(*xs, xs[*nc]))
                        {
                            int neighborLabel = xt[*nc];

                            // find the root label of a label tree
                            while(neighborLabel != label[neighborLabel])
                            {
                                neighborLabel = label[neighborLabel];
                            }

                            if(neighborLabel < *xt) // always keep the smallest among the possible neighbor labels
                            {
                                label[*xt] = neighborLabel;
                                *xt = neighborLabel;
                            }
                            else
                            {
                                label[neighborLabel] = *xt;
                            }
                        }
                        nc.turnTo(Neighborhood3D::nearBorderDirectionsCausal(atBorder,++j));
                    }
                }
            }
        }
    }

    // pass 2: assign one label to each region (tree)
    // so that labels form a consecutive sequence 1, 2, ...
    DestIterator zd = d_Iter;

    unsigned int count = 0; 

    i= 0;

    for(z=0; z != d; ++z, ++zd.dim2())
    {
        DestIterator yd(zd);

        for(y=0; y != h; ++y, ++yd.dim1())
        {
            DestIterator xd(yd);

            for(x = 0; x != w; ++x, ++xd.dim0(), ++i)
            {
                if(label[i] == i)
                {
                        label[i] = count++;
                }
                else
                {
                        label[i] = label[label[i]]; // compress trees
                }
                da.set(label[i]+1, xd);
            }
        }
    }
    return count;
}


/********************************************************/
/*                                                      */
/*               labelVolumeWithBackground              */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a segmented volume,
     excluding the background from labeling.

    Connected components are defined as regions with uniform voxel
    values. Thus, <TT>SrcAccessor::value_type</TT> either must be
    equality comparable (first form), or an EqualityFunctor must be
    provided that realizes the desired predicate (second form). All
    voxel equal to the given '<TT>background_value</TT>' are ignored
    when determining connected components and remain untouched in the
    destination volume.

    The destination's value type should be large enough to hold the
    labels without overflow. Region numbers will be a consecutive
    sequence starting with one and ending with the region number
    returned by the function (inclusive).

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {

        template <class SrcIterator, class SrcAccessor,class SrcShape,
                          class DestIterator, class DestAccessor,
                          class Neighborhood3D, class ValueType>
        unsigned int labelVolumeWithBackground(    SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                                                          DestIterator d_Iter, DestAccessor da,
                                                          Neighborhood3D neighborhood3D, ValueType background_value);

        template <class SrcIterator, class SrcAccessor,class SrcShape,
                          class DestIterator, class DestAccessor,
                          class Neighborhood3D, class ValueType, class EqualityFunctor>
        unsigned int labelVolumeWithBackground(    SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                                                            DestIterator d_Iter, DestAccessor da,
                                                          Neighborhood3D neighborhood3D, ValueType background_value,
                                                            EqualityFunctor equal);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {

        template <class SrcIterator, class SrcAccessor,class SrcShape,
                          class DestIterator, class DestAccessor,
                          class Neighborhood3D, class ValueType>
        unsigned int labelVolumeWithBackground(    triple<SrcIterator, SrcShape, SrcAccessor> src,
                                                          pair<DestIterator, DestAccessor> dest,
                                                          Neighborhood3D neighborhood3D, ValueType background_value);

        template <class SrcIterator, class SrcAccessor,class SrcShape,
                          class DestIterator, class DestAccessor,
                          class Neighborhood3D, class ValueType, class EqualityFunctor>
        unsigned int labelVolumeWithBackground(    triple<SrcIterator, SrcShape, SrcAccessor> src,
                                                        pair<DestIterator, DestAccessor> dest,
                                                        Neighborhood3D neighborhood3D, ValueType background_value,
                                                        EqualityFunctor equal);

    }
    \endcode

    Return:  the number of regions found (= largest region label)

    <b> Usage:</b>

    <b>\#include</b> "<a href="labelvolume_8hxx-source.html">vigra/labelvolume.hxx</a>"<br>
    Namespace: vigra

    \code
    typedef vigra::MultiArray<3,int> IntVolume;
    IntVolume src(IntVolume::difference_type(w,h,d));
    IntVolume dest(IntVolume::difference_type(w,h,d));

    // find 6-connected regions
    int max_region_label = vigra::labelVolumeWithBackground(
    srcMultiArrayRange(src), destMultiArray(dest), NeighborCode3DSix(), 0);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcIterator src_begin;
    SrcShape shape;
    DestIterator dest_begin;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_begin);

    u == u                      // first form

    EqualityFunctor equal;      // second form
    equal(u, u)                 // second form

    int i;
    dest_accessor.set(i, dest_begin);
    \endcode

*/
template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D,
          class ValueType>
unsigned int labelVolumeWithBackground(SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                                       DestIterator d_Iter, DestAccessor da,
                                       Neighborhood3D neighborhood3D, ValueType backgroundValue)
{
    return labelVolumeWithBackground(s_Iter, srcShape, sa, d_Iter, da, neighborhood3D, backgroundValue, std::equal_to<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D,
          class ValueType>
unsigned int labelVolumeWithBackground(triple<SrcIterator, SrcShape, SrcAccessor> src,
                                       pair<DestIterator, DestAccessor> dest,
                                       Neighborhood3D neighborhood3D, ValueType backgroundValue)
{
    return labelVolumeWithBackground(src.first, src.second, src.third, dest.first, dest.second, neighborhood3D, backgroundValue, std::equal_to<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D,
          class ValueType, class EqualityFunctor>
unsigned int labelVolumeWithBackground(triple<SrcIterator, SrcShape, SrcAccessor> src,
                                       pair<DestIterator, DestAccessor> dest,
                                       Neighborhood3D neighborhood3D, ValueType backgroundValue, EqualityFunctor equal)
{
    return labelVolumeWithBackground(src.first, src.second, src.third, dest.first, dest.second, neighborhood3D, backgroundValue, equal);
}
 
template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D,
          class ValueType, class EqualityFunctor>
unsigned int labelVolumeWithBackground(SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                                       DestIterator d_Iter, DestAccessor da,
                                       Neighborhood3D,
                                       ValueType backgroundValue, EqualityFunctor equal)
{

    //basically needed for iteration and border-checks
    int w = srcShape[0], h = srcShape[1], d = srcShape[2];
    int x,y,z, i;
        
    //declare and define Iterators for all three dims at src
    SrcIterator zs = s_Iter;
    SrcIterator ys(zs);
    SrcIterator xs(ys);
        
    // temporary image to store region labels
    typedef vigra::MultiArray<3,int> LabelVolume;
    typedef LabelVolume::traverser LabelTraverser;
    
    LabelVolume labelvolume(srcShape);
        
    //Declare traversers for all three dims at target
    LabelTraverser zt = labelvolume.traverser_begin();
    LabelTraverser yt(zt);
    LabelTraverser xt(yt);

    // Kovalevsky's clever idea to use
    // image iterator and scan order iterator simultaneously
    // memory order indicates label
    LabelVolume::iterator label = labelvolume.begin();
    
    // initialize the neighborhood traversers

#if 0
    NeighborOffsetTraverser<Neighborhood3D> nc(Neighborhood3D::CausalFirst);
    NeighborOffsetTraverser<Neighborhood3D> nce(Neighborhood3D::CausalLast);
#endif /* #if 0 */

    NeighborOffsetCirculator<Neighborhood3D> nc(Neighborhood3D::CausalFirst);
    NeighborOffsetCirculator<Neighborhood3D> nce(Neighborhood3D::CausalLast);
    ++nce;
    // pass 1: scan image from upper left front to lower right back
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
    i=0;
    int backgroundAddress=-1;
    for(z = 0; z != d; ++z, ++zs.dim2(), ++zt.dim2())
    {
        ys = zs;
        yt = zt;

        for(y = 0; y != h; ++y, ++ys.dim1(), ++yt.dim1())
        {
            xs = ys;
            xt = yt;

            for(x = 0; x != w; ++x, ++xs.dim0(), ++xt.dim0(), ++i)
            {
                if(equal(*xs, backgroundValue)){
                        if(backgroundAddress==-1) backgroundAddress = i;
                        label[i]=backgroundAddress;
                }
                else{
                    *xt = i; // default: new region    

                    //queck whether there is a special borde threatment to be used or not
                    AtVolumeBorder atBorder = isAtVolumeBorderCausal(x,y,z,w,h,z);
                        
                    //We are not at the border!
                    if(atBorder == NotAtBorder)
                    {
                                

#if 0
                        nc = NeighborOffsetTraverser<Neighborhood3D>(Neighborhood3D::CausalFirst);
#endif /* #if 0 */

                        nc = NeighborOffsetCirculator<Neighborhood3D>(Neighborhood3D::CausalFirst);
                    
                        do
                        {            
                            // if colors are equal
                            if(equal(*xs, xs[*nc]))
                            {
                                int neighborLabel = xt[*nc];

                                // find the root label of a label tree
                                while(neighborLabel != label[neighborLabel])
                                {
                                    neighborLabel = label[neighborLabel];
                                }

                                if(neighborLabel < *xt) // always keep the smallest among the possible neighbor labels
                                {
                                    label[*xt] = neighborLabel;
                                    *xt = neighborLabel;
                                }
                                else
                                {
                                    label[neighborLabel] = *xt;
                                }
                            }
                            ++nc;
                        }while(nc!=nce);
                    }
                    //we are at a border - handle this!!
                    else
                    {

#if 0
                        nc = NeighborOffsetTraverser<Neighborhood3D>(Neighborhood3D::nearBorderDirectionsCausal(atBorder,0));
#endif /* #if 0 */

                        nc = NeighborOffsetCirculator<Neighborhood3D>(Neighborhood3D::nearBorderDirectionsCausal(atBorder,0));
                        int j=0;
                        while(nc.direction() != Neighborhood3D::Error)
                        {
                            //   colors equal???
                            if(equal(*xs, xs[*nc]))
                            {
                                int neighborLabel = xt[*nc];

                                // find the root label of a label tree
                                while(neighborLabel != label[neighborLabel])
                                {
                                    neighborLabel = label[neighborLabel];
                                }

                                if(neighborLabel < *xt) // always keep the smallest among the possible neighbor labels
                                {
                                    label[*xt] = neighborLabel;
                                    *xt = neighborLabel;
                                }
                                else
                                {
                                    label[neighborLabel] = *xt;
                                }
                            }
                            nc.turnTo(Neighborhood3D::nearBorderDirectionsCausal(atBorder,++j));
                        }
                    }
                }
            }
        }
    }

    // pass 2: assign one label to each region (tree)
    // so that labels form a consecutive sequence 1, 2, ...
    DestIterator zd = d_Iter;

    unsigned int count = 0; 

    i= 0;

    for(z=0; z != d; ++z, ++zd.dim2())
    {
        DestIterator yd(zd);

        for(y=0; y != h; ++y, ++yd.dim1())
        {
            DestIterator xd(yd);

            for(x = 0; x != w; ++x, ++xd.dim0(), ++i)
            {
                if(label[i] == backgroundAddress){
                    label[i]=0;
                    continue;
                }
                
                if(label[i] == i)
                {
                    label[i] = count++;
                }
                else
                {
                    label[i] = label[label[i]]; // compress trees
                }
                da.set(label[i]+1, xd);
            }
        }
    }
    return count;
}

//@}

} //end of namespace vigra

#endif //VIGRA_LABELVOLUME_HXX
