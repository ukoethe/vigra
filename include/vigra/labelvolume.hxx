/************************************************************************/
/*                                                                      */
/*     Copyright 2006-2007 by F. Heinrich, B. Seppke, Ullrich Koethe    */
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

#ifndef VIGRA_LABELVOLUME_HXX
#define VIGRA_LABELVOLUME_HXX


#include "voxelneighborhood.hxx"
#include "multi_array.hxx"
#include "union_find.hxx"

namespace vigra{

/** \addtogroup Labeling Connected Components Labeling
     The 3-dimensional connected components algorithms may use either 6 or 26 connectivity.
     By means of a functor the merge criterion can be defined arbitrarily.
*/
//@{

/********************************************************/
/*                                                      */
/*                        labelVolume                   */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a segmented volume.

    Connected components are defined as regions with uniform voxel
    values. Thus, <TT>T1</TT> either must be equality comparable, 
    or an EqualityFunctor must be provided explicitly that realizes 
    the desired equivalence predicate. The destination's value type 
    <tt>T2</tt> should be large enough to hold the labels
    without overflow. Region numbers will be a consecutive sequence
    starting with one and ending with the region number returned by
    the function (inclusive).

    Return:  the number of regions found (= largest region label)
    
    See \ref labelMultiArray() for a dimension-independent implementation of 
    connected components labelling.

    <b> Declarations:</b>

    pass 3D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class Neighborhood3D, 
                  class EqualityFunctor = std::equal_to<T1> >
        unsigned int
        labelVolume(MultiArrayView<3, T1, S1> const & source,
                    MultiArrayView<3, T2, S2> dest,
                    Neighborhood3D neighborhood3D,
                    EqualityFunctor equal  = EqualityFunctor());

    }
    \endcode

    \deprecatedAPI{labelVolume}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
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
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    \deprecatedEnd    

    <b> Usage:</b>

    <b>\#include</b> \<vigra/labelvolume.hxx\><br>
    Namespace: vigra

    \code
    typedef MultiArray<3,int> IntVolume;
    IntVolume src(Shape3(w,h,d));
    IntVolume dest(Shape3(w,h,d));
    
    // find 6-connected regions
    int max_region_label = labelVolumeSix(src, dest);

    // find 26-connected regions
    int max_region_label = labelVolume(src, dest, NeighborCode3DTwentySix());
    \endcode

    \deprecatedUsage{labelVolume}
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
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> unsigned int labelVolume)


template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D, class EqualityFunctor>
unsigned int labelVolume(SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                         DestIterator d_Iter, DestAccessor da,
                         Neighborhood3D, EqualityFunctor equal)
{
    typedef typename DestAccessor::value_type LabelType;
    
    //basically needed for iteration and border-checks
    int w = srcShape[0], h = srcShape[1], d = srcShape[2];
    int x,y,z;       
        
    // temporary image to store region labels
    UnionFindArray<LabelType>  label;
        
    //Declare traversers for all three dims at target
    SrcIterator zs = s_Iter;
    DestIterator zd = d_Iter;

    // initialize the neighborhood traversers
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
    for(z = 0; z != d; ++z, ++zs.dim2(), ++zd.dim2())
    {
        SrcIterator ys(zs);
        DestIterator yd(zd);

        for(y = 0; y != h; ++y, ++ys.dim1(), ++yd.dim1())
        {
            SrcIterator xs(ys);
            DestIterator xd(yd);

            for(x = 0; x != w; ++x, ++xs.dim0(), ++xd.dim0())
            {
                LabelType currentIndex = label.nextFreeIndex();

                //check whether there is a special border treatment to be used or not
                AtVolumeBorder atBorder = isAtVolumeBorderCausal(x,y,z,w,h,d);

                //We are not at the border!
                if(atBorder == NotAtBorder)
                {
                    NeighborOffsetCirculator<Neighborhood3D> nc(Neighborhood3D::CausalFirst);
                
                    do
                    {            
                        // if colors are equal
                        if(equal(sa(xs), sa(xs, *nc)))
                        {
                            currentIndex = label.makeUnion(da(xd,*nc), currentIndex);
                        }
                        ++nc;
                    }
                    while(nc!=nce);
                }               
                else //we are at a border - handle this!!
                {
                    NeighborOffsetCirculator<Neighborhood3D> nc(Neighborhood3D::nearBorderDirectionsCausal(atBorder,0));
                    int j=0;
                    while(nc.direction() != Neighborhood3D::Error)
                    {
                        int dummy = x+(*nc)[0];  // prevents an apparently incorrect optimization in gcc 4.8
                        if (dummy<0)
                        {  
                            std::cerr << "internal error " << dummy << std::endl;
                        }
                        //   colors equal???
                        if(equal(sa(xs), sa(xs, *nc)))
                        {
                            currentIndex = label.makeUnion(da(xd,*nc), currentIndex);
                        }
                        nc.turnTo(Neighborhood3D::nearBorderDirectionsCausal(atBorder,++j));
                    }
                }
                da.set(label.finalizeIndex(currentIndex), xd);
            }
        }
    }
    
    LabelType count = label.makeContiguous();

    // pass 2: assign one label to each region (tree)
    // so that labels form a consecutive sequence 1, 2, ...
    zd = d_Iter;
    for(z=0; z != d; ++z, ++zd.dim2())
    {
        DestIterator yd(zd);

        for(y=0; y != h; ++y, ++yd.dim1())
        {
            DestIterator xd(yd);

            for(x = 0; x != w; ++x, ++xd.dim0())
            {
                da.set(label.findLabel(da(xd)), xd);
            }
        }
    }
    return count;
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

template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor,
          class Neighborhood3D>
unsigned int labelVolume(triple<SrcIterator, SrcShape, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest,
                         Neighborhood3D neighborhood3D)
{
    return labelVolume(src.first, src.second, src.third, dest.first, dest.second, neighborhood3D, std::equal_to<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood3D, class EqualityFunctor>
unsigned int labelVolume(triple<SrcIterator, SrcShape, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest,
                         Neighborhood3D neighborhood3D, EqualityFunctor equal)
{
    return labelVolume(src.first, src.second, src.third, dest.first, dest.second, neighborhood3D, equal);
}

template <class T1, class S1, 
          class T2, class S2,
          class Neighborhood3D, class EqualityFunctor>
inline unsigned int
labelVolume(MultiArrayView<3, T1, S1> const & source,
            MultiArrayView<3, T2, S2> dest,
            Neighborhood3D neighborhood3D,
            EqualityFunctor equal)
{
    vigra_precondition(source.shape() == dest.shape(),
        "labelVolume(): shape mismatch between input and output.");
    return labelVolume(srcMultiArrayRange(source), destMultiArray(dest), neighborhood3D, equal);
}

template <class T1, class S1, 
          class T2, class S2,
          class Neighborhood3D>
inline unsigned int
labelVolume(MultiArrayView<3, T1, S1> const & source,
            MultiArrayView<3, T2, S2> dest,
            Neighborhood3D neighborhood3D)
{
    vigra_precondition(source.shape() == dest.shape(),
        "labelVolume(): shape mismatch between input and output.");
    return labelVolume(srcMultiArrayRange(source), destMultiArray(dest), neighborhood3D, std::equal_to<T1>());
}

/********************************************************/
/*                                                      */
/*                    labelVolumeSix                    */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a segmented volume
     using the 6-neighborhood.
     
     See \ref labelVolume() for detailed documentation.

*/
template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor>
unsigned int labelVolumeSix(triple<SrcIterator, SrcShape, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest)
{
    return labelVolume(src.first, src.second, src.third, dest.first, dest.second, NeighborCode3DSix(), std::equal_to<typename SrcAccessor::value_type>());
}

template <class T1, class S1,
          class T2, class S2>
unsigned int labelVolumeSix(MultiArrayView<3, T1, S1> const & source,
                            MultiArrayView<3, T2, S2> dest)
{
    return labelVolume(srcMultiArrayRange(source), destMultiArray(dest), 
                       NeighborCode3DSix(), std::equal_to<T1>());
}

/********************************************************/
/*                                                      */
/*               labelVolumeWithBackground              */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a segmented volume,
     excluding the background from labeling.

    This function works like \ref labelVolume(), but considers all background voxels
    (i.e. voxels having the given '<TT>background_value</TT>') as a single region that
    is ignored when determining connected components and remains untouched in the
    destination array. Usually, you will zero-initialize the output array, so that
    the background gets label 0 (remember that actual region labels start at one).

    Return:  the number of regions found (= largest region label)

    See \ref labelMultiArrayWithBackground() for a dimension-independent implementation
    if this algorithm.

    <b> Declarations:</b>

    pass 3D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class Neighborhood3D,
                  class ValueType,
                  class EqualityFunctor = std::equalt_to<T1> >
        unsigned int
        labelVolumeWithBackground(MultiArrayView<3, T1, S1> const & source,
                                  MultiArrayView<3, T2, S2> dest,
                                  Neighborhood3D neighborhood3D,
                                  ValueType backgroundValue,
                                  EqualityFunctor equal = EqualityFunctor());
    }
    \endcode

    \deprecatedAPI{labelVolumeWithBackground}
    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
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
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/labelvolume.hxx\><br>
    Namespace: vigra

    \code
    typedef vigra::MultiArray<3,int> IntVolume;
    
    IntVolume src(Shape3(w,h,d));
    IntVolume dest(Shape3(w,h,d));

    // find 6-connected regions
    int max_region_label = labelVolumeWithBackground(src, dest, NeighborCode3DSix(), 0);
    \endcode

    \deprecatedUsage{labelVolumeWithBackground}
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
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> unsigned int labelVolumeWithBackground)

template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D,
          class ValueType, class EqualityFunctor>
unsigned int labelVolumeWithBackground(SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                                       DestIterator d_Iter, DestAccessor da,
                                       Neighborhood3D,
                                       ValueType backgroundValue, EqualityFunctor equal)
{
    typedef typename DestAccessor::value_type LabelType;
    
    //basically needed for iteration and border-checks
    int w = srcShape[0], h = srcShape[1], d = srcShape[2];
    int x,y,z;       
        
    // temporary image to store region labels
    UnionFindArray<LabelType>  label;
        
    //Declare traversers for all three dims at target
    SrcIterator zs = s_Iter;
    DestIterator zd = d_Iter;

    // initialize the neighborhood traversers
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
    for(z = 0; z != d; ++z, ++zs.dim2(), ++zd.dim2())
    {
        SrcIterator ys(zs);
        DestIterator yd(zd);

        for(y = 0; y != h; ++y, ++ys.dim1(), ++yd.dim1())
        {
            SrcIterator xs(ys);
            DestIterator xd(yd);

            for(x = 0; x != w; ++x, ++xs.dim0(), ++xd.dim0())
            {
                if(equal(sa(xs), backgroundValue))
                {
                    //da.set(label.getIndex(0), xd);
                    da.set(0, xd);
                    continue;
                }

                LabelType currentIndex = label.nextFreeIndex();

                //check whether there is a special border treatment to be used or not
                AtVolumeBorder atBorder = isAtVolumeBorderCausal(x,y,z,w,h,d);
                    
                //We are not at the border!
                if(atBorder == NotAtBorder)
                {
                    NeighborOffsetCirculator<Neighborhood3D> nc(Neighborhood3D::CausalFirst);
                
                    do
                    {            
                        // if colors are equal
                        if(equal(sa(xs), sa(xs, *nc)))
                        {
                            currentIndex = label.makeUnion(da(xd,*nc), currentIndex);
                        }
                        ++nc;
                    }
                    while(nc!=nce);
                }               
                else //we are at a border - handle this!!
                {
                    NeighborOffsetCirculator<Neighborhood3D> nc(Neighborhood3D::nearBorderDirectionsCausal(atBorder,0));
                    int j=0;
                    while(nc.direction() != Neighborhood3D::Error)
                    {
                        int dummy = x+(*nc)[0];  // prevents an apparently incorrect optimization in gcc 4.8
                        if (dummy<0)
                        {  
                            std::cerr << "internal error " << dummy << std::endl;
                        }
                        //   colors equal???
                        if(equal(sa(xs), sa(xs, *nc)))
                        {
                            currentIndex = label.makeUnion(da(xd,*nc), currentIndex);
                        }
                        nc.turnTo(Neighborhood3D::nearBorderDirectionsCausal(atBorder,++j));
                    }
                }
                da.set(label.finalizeIndex(currentIndex), xd);
            }
        }
    }
    
    LabelType count = label.makeContiguous();

    // pass 2: assign one label to each region (tree)
    // so that labels form a consecutive sequence 1, 2, ...
    zd = d_Iter;
    for(z=0; z != d; ++z, ++zd.dim2())
    {
        DestIterator yd(zd);

        for(y=0; y != h; ++y, ++yd.dim1())
        {
            DestIterator xd(yd);

            for(x = 0; x != w; ++x, ++xd.dim0())
            {
                da.set(label.findLabel(da(xd)), xd);
            }
        }
    }
    return count;
}

template <class SrcIterator, class SrcAccessor,class SrcShape,
          class DestIterator, class DestAccessor,
          class Neighborhood3D,
          class ValueType>
inline unsigned int 
labelVolumeWithBackground(SrcIterator s_Iter, SrcShape srcShape, SrcAccessor sa,
                          DestIterator d_Iter, DestAccessor da,
                          Neighborhood3D neighborhood3D, ValueType backgroundValue)
{
    return labelVolumeWithBackground(s_Iter, srcShape, sa, d_Iter, da, neighborhood3D, backgroundValue, std::equal_to<typename SrcAccessor::value_type>());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood3D,
          class ValueType,
          class EqualityFunctor>
inline unsigned int
labelVolumeWithBackground(triple<SrcIterator, SrcShape, SrcAccessor> src,
                          pair<DestIterator, DestAccessor> dest,
                          Neighborhood3D neighborhood3D, ValueType backgroundValue, EqualityFunctor equal)
{
    return labelVolumeWithBackground(src.first, src.second, src.third, dest.first, dest.second, neighborhood3D, backgroundValue, equal);
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood3D,
          class ValueType>
inline unsigned int
labelVolumeWithBackground(triple<SrcIterator, SrcShape, SrcAccessor> src,
                          pair<DestIterator, DestAccessor> dest,
                          Neighborhood3D neighborhood3D, ValueType backgroundValue)
{
    return labelVolumeWithBackground(src.first, src.second, src.third, dest.first, dest.second, 
                                     neighborhood3D, backgroundValue, std::equal_to<typename SrcAccessor::value_type>());
}

template <class T1, class S1,
          class T2, class S2,
          class Neighborhood3D,
          class ValueType,
          class EqualityFunctor>
inline unsigned int
labelVolumeWithBackground(MultiArrayView<3, T1, S1> const & source,
                          MultiArrayView<3, T2, S2> dest,
                          Neighborhood3D neighborhood3D,
                          ValueType backgroundValue,
                          EqualityFunctor equal)
{
    vigra_precondition(source.shape() == dest.shape(),
        "labelVolumeWithBackground(): shape mismatch between input and output.");
    return labelVolumeWithBackground(srcMultiArrayRange(source), destMultiArray(dest), 
                                     neighborhood3D, backgroundValue, equal);
}

template <class T1, class S1,
          class T2, class S2,
          class Neighborhood3D,
          class ValueType>
inline unsigned int
labelVolumeWithBackground(MultiArrayView<3, T1, S1> const & source,
                          MultiArrayView<3, T2, S2> dest,
                          Neighborhood3D neighborhood3D,
                          ValueType backgroundValue)
{
    vigra_precondition(source.shape() == dest.shape(),
        "labelVolumeWithBackground(): shape mismatch between input and output.");
    return labelVolumeWithBackground(srcMultiArrayRange(source), destMultiArray(dest), 
                                     neighborhood3D, backgroundValue,
                                     std::equal_to<T1>());
}

//@}

} //end of namespace vigra

#endif //VIGRA_LABELVOLUME_HXX
