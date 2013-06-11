/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyanalysis_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/localminmax.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/watersheds.hxx>
#include <vigra/seededregiongrowing.hxx>
#include <vigra/labelvolume.hxx>
#include <vigra/watersheds3d.hxx>
#include <vigra/seededregiongrowing3d.hxx>
#include <vigra/multi_watersheds.hxx>
#include <vigra/convolution.hxx>
#include <vigra/multi_convolution.hxx>
#include <vigra/slic.hxx>

#include <string>
#include <cmath>

#include "tws.hxx"

namespace python = boost::python;

namespace vigra
{

template < class PixelType >
NumpyAnyArray 
pythonLabelImage(NumpyArray<2, Singleband<PixelType> > image,
                 int neighborhood = 4, 
                 NumpyArray<2, Singleband<npy_uint32> > res = NumpyArray<2, Singleband<npy_uint32> >())
{
    vigra_precondition(neighborhood == 4 || neighborhood == 8,
        "labelImage(): neighborhood must be 4 or 8.");

    std::string description("connected components, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
            "labelImage(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        switch (neighborhood)
        {
            case 4:
            {
                labelImage(srcImageRange(image), destImage(res), false);
                break;
            }
            case 8:
            {
                labelImage(srcImageRange(image), destImage(res), true);
                break;
            }
        }
    }
    
    return res;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pyLabelImage, pythonLabelImage)

template < class PixelType >
NumpyAnyArray 
pythonLabelImageWithBackground(NumpyArray<2, Singleband<PixelType> > image,
                               int neighborhood = 4,
                               PixelType background_value = 0,
                               NumpyArray<2, Singleband<npy_uint32> > res = NumpyArray<2, Singleband<npy_uint32> >())
{
    vigra_precondition(neighborhood == 4 || neighborhood == 8,
        "labelImageWithBackground(): neighborhood must be 4 or 8.");

    std::string description("connected components with background, neighborhood=");
    description += asString(neighborhood)+ ", bglabel=" + asString(background_value);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
        "labelImageWithBackground(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        switch (neighborhood)
        {
            case 4:
            {
                labelImageWithBackground(srcImageRange(image),
                    destImage(res), false, background_value);
                break;
            }
            case 8:
            {
                labelImageWithBackground(srcImageRange(image),
                    destImage(res), true, background_value);
                break;
            }
        }
    }
    return res;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pyLabelImageWithBackground, pythonLabelImageWithBackground)

template < class VoxelType >
NumpyAnyArray 
pythonLabelVolume(NumpyArray<3, Singleband<VoxelType> > volume, 
                  int neighborhood=6,
                  NumpyArray<3, Singleband<npy_uint32> > res = NumpyArray<3, Singleband<npy_uint32> >())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
        "labelVolume(): neighborhood must be 6 or 26.");
    
    std::string description("connected components, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description), 
            "labelVolume(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        switch (neighborhood)
        {
            case 6:
            {
                labelVolume(srcMultiArrayRange(volume),
                    destMultiArray(res), NeighborCode3DSix());
                break;
            }
            case 26:
            {
                labelVolume(srcMultiArrayRange(volume),
                    destMultiArray(res), NeighborCode3DTwentySix());
                break;
            }
        }
    }
    return res;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pyLabelVolume, pythonLabelVolume)

template < class VoxelType >
NumpyAnyArray 
pythonLabelVolumeWithBackground(NumpyArray<3, Singleband<VoxelType> > volume, 
                                int neighborhood=6,
                                VoxelType background_value = 0,
                                NumpyArray<3, Singleband<npy_uint32> > res = NumpyArray<3, Singleband<npy_uint32> >())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
        "labelVolumeWithBackground(): neighborhood must be 6 or 26.");
    
    std::string description("connected components with background, neighborhood=");
    description += asString(neighborhood)+ ", bglabel=" + asString(background_value);
    
    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description), 
        "labelVolumeWithBackground(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        switch (neighborhood)
        {
            case 6:
            {
                labelVolumeWithBackground(srcMultiArrayRange(volume),
                    destMultiArray(res), NeighborCode3DSix(),
                    background_value);
                break;
            }
            case 26:
            {
                labelVolumeWithBackground(srcMultiArrayRange(volume),
                    destMultiArray(res), NeighborCode3DTwentySix(),
                    background_value);
                break;
            }
        }
    }
    return res;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pyLabelVolumeWithBackground, pythonLabelVolumeWithBackground)

/*********************************************************************************/

// FIXME: support output of label images from localMinim/Maxima functions

template < class PixelType >
NumpyAnyArray 
pythonLocalMinima2D(NumpyArray<2, Singleband<PixelType> > image,
                    PixelType marker = NumericTraits<PixelType>::one(),
                    int neighborhood = 8,
                    NumpyArray<2, Singleband<PixelType> > res = NumpyArray<2, Singleband<PixelType> >())
{
    vigra_precondition(neighborhood == 4 || neighborhood == 8,
        "localMinima(): neighborhood must be 4 or 8.");

    std::string description("local minima, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
            "localMinima(): Output array has wrong shape.");

    {
            PyAllowThreads _pythread;
        switch (neighborhood)
        {
            case 4:
            {
                localMinima(srcImageRange(image), destImage(res), marker,
                    FourNeighborCode());
                break;
            }
            case 8:
            {
                localMinima(srcImageRange(image), destImage(res), marker,
                    EightNeighborCode());
                break;
            }
        }
    }
    
    return res;
}

template<class PixelType>
NumpyAnyArray 
pythonLocalMinima3D(NumpyArray<3, Singleband<PixelType> > volume,
                    PixelType marker = NumericTraits<PixelType>::one(), 
                    int neighborhood = 6, 
                    NumpyArray<3, Singleband<PixelType> > res = NumpyArray<3, Singleband<PixelType> >())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
            "localMinima(): neighborhood must be 6 or 26.");

    std::string description("local minima, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description), 
            "localMinima(): Output array has wrong shape.");
            
    switch (neighborhood)
    {
        case 6:
        {
            localMinima3D(srcMultiArrayRange(volume), destMultiArray(res), marker,
                    NeighborCode3DSix());
            break;
        }
        case 26:
        {
            localMinima3D(srcMultiArrayRange(volume), destMultiArray(res), marker,
                    NeighborCode3DTwentySix());
            break;
        }
    }

    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonExtendedLocalMinima2D(NumpyArray<2, Singleband<PixelType> > image,
                            PixelType marker = NumericTraits<PixelType>::one(),
                            int neighborhood = 8,
                            NumpyArray<2, Singleband<PixelType> > res = NumpyArray<2, Singleband<PixelType> >())
{
    vigra_precondition(neighborhood == 4 || neighborhood == 8,
        "extendedLocalMinima(): neighborhood must be 4 or 8.");

    std::string description("extended local minima, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
        "extendedLocalMinima(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        switch (neighborhood)
        {
            case 4:
            {
                extendedLocalMinima(srcImageRange(image), destImage(res),
                    marker, FourNeighborCode());
                break;
            }
            case 8:
            {
                extendedLocalMinima(srcImageRange(image), destImage(res),
                    marker, EightNeighborCode());
                break;
            }
        }
    }
    return res;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pyExtendedLocalMinima2D, pythonExtendedLocalMinima2D)

template<class PixelType>
NumpyAnyArray 
pythonExtendedLocalMinima3D(NumpyArray<3, Singleband<PixelType> > volume, 
                            PixelType marker = NumericTraits<PixelType>::one(), 
                            int neighborhood = 6,
                            NumpyArray<3, Singleband<PixelType> > res = NumpyArray<3, Singleband<PixelType> >())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
            "extendedLocalMinima(): neighborhood must be 6 or 26.");

    std::string description("extended local minima, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description), 
            "extendedLocalMinima(): Output array has wrong shape.");
    switch (neighborhood)
    {
        case 6:
        {
            extendedLocalMinima3D(srcMultiArrayRange(volume), destMultiArray(res),
                marker, NeighborCode3DSix());
            break;
        }
        case 26:
        {
            extendedLocalMinima3D(srcMultiArrayRange(volume), destMultiArray(res),
                marker, NeighborCode3DTwentySix());
            break;
        }
    }

    return res;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pyExtendedLocalMinima3D, pythonExtendedLocalMinima3D)

template < class PixelType >
NumpyAnyArray 
pythonLocalMaxima2D(NumpyArray<2, Singleband<PixelType> > image,
                    PixelType marker = NumericTraits<PixelType>::one(),
                    int neighborhood = 8,
                    NumpyArray<2, Singleband<PixelType> > res = NumpyArray<2, Singleband<PixelType> >())
{
    vigra_precondition(neighborhood == 4 || neighborhood == 8,
        "localMaxima(): neighborhood must be 4 or 8.");

    std::string description("local maxima, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
            "localMaxima(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        switch (neighborhood)
        {
            case 4:
            {
                localMaxima(srcImageRange(image), destImage(res), marker,
                    FourNeighborCode());
                break;
            }
            case 8:
            {
                localMaxima(srcImageRange(image), destImage(res), marker,
                    EightNeighborCode());
                break;
            }
        }
    }
    
    return res;
}

template<class PixelType>
NumpyAnyArray 
pythonLocalMaxima3D(NumpyArray<3, Singleband<PixelType> > volume,
                    PixelType marker = NumericTraits<PixelType>::one(), 
                    int neighborhood = 6, 
                    NumpyArray<3, Singleband<PixelType> > res = NumpyArray<3, Singleband<PixelType> >())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
            "localMaxima(): neighborhood must be 6 or 26.");

    std::string description("local maxima, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description), 
            "localMaxima(): Output array has wrong shape.");
    switch (neighborhood)
    {
        case 6:
        {
            localMaxima3D(srcMultiArrayRange(volume), destMultiArray(res), marker,
                NeighborCode3DSix());
            break;
        }
        case 26:
        {
            localMaxima3D(srcMultiArrayRange(volume), destMultiArray(res), marker,
                NeighborCode3DTwentySix());
            break;
        }
    }

    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonExtendedLocalMaxima2D(NumpyArray<2, Singleband<PixelType> > image,
                            PixelType marker = NumericTraits<PixelType>::one(),
                            int neighborhood = 8,
                            NumpyArray<2, Singleband<PixelType> > res = NumpyArray<2, Singleband<PixelType> >())
{
    vigra_precondition(neighborhood == 4 || neighborhood == 8,
        "extendedLocalMaxima(): neighborhood must be 4 or 8.");

    std::string description("extended local maxima, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
            "extendedLocalMaxima(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        switch (neighborhood)
        {
            case 4:
            {
                extendedLocalMaxima(srcImageRange(image), destImage(res),
                    marker, FourNeighborCode());
                break;
            }
            case 8:
            {
                extendedLocalMaxima(srcImageRange(image), destImage(res),
                    marker, EightNeighborCode());
                break;
            }
        }
    }
    return res;
}

template<class PixelType>
NumpyAnyArray 
pythonExtendedLocalMaxima3D(NumpyArray<3, Singleband<PixelType> > volume, 
                            PixelType marker = NumericTraits<PixelType>::one(), 
                            int neighborhood = 6,
                            NumpyArray<3, Singleband<PixelType> > res = NumpyArray<3, Singleband<PixelType> >())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
            "extendedLocalMaxima(): neighborhood must be 6 or 26.");

    std::string description("extended local maxima, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description), 
            "extendedLocalMaxima(): Output array has wrong shape.");
    switch (neighborhood)
    {
        case 6:
        {
            extendedLocalMaxima3D(srcMultiArrayRange(volume), destMultiArray(res),
                                  marker, NeighborCode3DSix());
            break;
        }
        case 26:
        {
            extendedLocalMaxima3D(srcMultiArrayRange(volume), destMultiArray(res),
                                  marker, NeighborCode3DTwentySix());
            break;
        }
    }

    return res;
}

/*************************************************************************/

#if 0
template < class PixelType >
python::tuple 
pythonWatersheds2DOld(NumpyArray<2, Singleband<PixelType> > image,
                   int neighborhood = 4,
                   NumpyArray<2, Singleband<npy_uint32> > seeds = python::object(),
                   std::string method = "RegionGrowing", 
                   SRGType srgType = CompleteGrow, 
                   PixelType max_cost = 0.0, 
                   NumpyArray<2, Singleband<npy_uint32> > res = NumpyArray<2, Singleband<npy_uint32> >())
{
    vigra_precondition(neighborhood == 4 || neighborhood == 8,
           "watersheds2D(): neighborhood must be 4 or 8.");

    method = tolower(method);
    
    bool haveSeeds = seeds.hasData();
    unsigned int maxRegionLabel = 0;
    
    if(method == "")
        method = "regiongrowing";
    
    if(method == "regiongrowing")
    {
        seeds.reshapeIfEmpty(image.shape(), 
                "watersheds(): Seed array has wrong shape.");
        
        if(!haveSeeds)
        {
            MultiArray<2, UInt8> minima(image.shape());
            localMinima(srcImageRange(image), destImage(minima), 1, EightNeighborCode());
            maxRegionLabel = labelImageWithBackground(srcImageRange(minima), destImage(seeds), true, 0);
        }
        else
        {
            FindMinMax< npy_uint32 > minmax;
            inspectImage(srcImageRange(seeds), minmax);
            maxRegionLabel = minmax.max;
        }
           
        res.reshapeIfEmpty(image.shape(), "watersheds(): Output array has wrong shape.");

        ArrayOfRegionStatistics< SeedRgDirectValueFunctor< PixelType > > stats(maxRegionLabel);
        if(neighborhood == 4)
        {
            seededRegionGrowing(srcImageRange(image), srcImage(seeds), destImage(res), 
                                stats, srgType, FourNeighborCode(), max_cost);
        }
        else
        {
            seededRegionGrowing(srcImageRange(image), srcImage(seeds), destImage(res), 
                                stats, srgType, EightNeighborCode(), max_cost);
        }
    }
    else if(method == "unionfind")
    {
        vigra_precondition(!haveSeeds,
           "watersheds(): UnionFind does not support seed images.");
        vigra_precondition(srgType == CompleteGrow,
           "watersheds(): UnionFind only supports 'CompleteGrow' mode.");
           
        res.reshapeIfEmpty(image.shape(), "watersheds(): Output array has wrong shape.");
        
        if(neighborhood == 4)
        {
            maxRegionLabel = watershedsUnionFind(srcImageRange(image), destImage(res),
                                        FourNeighborCode());
        }
        else
        {
            maxRegionLabel = watershedsUnionFind(srcImageRange(image), destImage(res),
                                        EightNeighborCode());
        }
    }
    else
    {
        vigra_precondition(false, "watersheds(): Unknown watershed method requested.");
    }

    return python::make_tuple(res, maxRegionLabel);
}
#endif

template < class PixelType >
python::tuple 
pythonWatersheds2D(NumpyArray<2, Singleband<PixelType> > image,
                   int neighborhood = 4,
                   NumpyArray<2, Singleband<npy_uint32> > seeds = NumpyArray<2, Singleband<npy_uint32> >(),
                   std::string method = "", 
                   SRGType srgType = CompleteGrow, 
                   PixelType max_cost = 0.0, 
                   NumpyArray<2, Singleband<npy_uint32> > res = NumpyArray<2, Singleband<npy_uint32> >())
{
    vigra_precondition(neighborhood == 4 || neighborhood == 8,
           "watersheds2D(): neighborhood must be 4 or 8.");

    method = tolower(method);
    if(method == "")
    {
        if(IsSameType<PixelType, npy_uint8>::value)
            method = "turbo";
        else
            method = "regiongrowing";
    }
    
    std::string description("watershed labeling, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
            "watersheds(): Output array has wrong shape.");
    
    WatershedOptions options;
    options.srgType(srgType);
    
    if(max_cost > 0.0)
    {
        vigra_precondition(method != "unionfind",
           "watersheds(): UnionFind does not support a cost threshold.");
        options.stopAtThreshold(max_cost);
    }
    
    if(seeds.hasData())
    {
        vigra_precondition(method != "unionfind",
           "watersheds(): UnionFind does not support seed images.");
        res = seeds;
    }
    else
    {
        if(method == "turbo")
            options.seedOptions(SeedOptions().extendedMinima());
        else
            options.seedOptions(SeedOptions().minima());
    }
    
    if(method == "turbo")
    {
        vigra_precondition((IsSameType<PixelType, npy_uint8>::value),
           "watersheds(): Turbo method only works for uint8 images.");
        options.turboAlgorithm();
        method = "regiongrowing";
    }
    
    npy_uint32 maxRegionLabel = 0;
    if(method == "regiongrowing")
    {
        PyAllowThreads _pythread;
        if(neighborhood == 4)
        {
            maxRegionLabel = watershedsRegionGrowing(srcImageRange(image), destImage(res), 
                                    FourNeighborCode(), options);
        }
        else
        {
            maxRegionLabel = watershedsRegionGrowing(srcImageRange(image), destImage(res), 
                                    EightNeighborCode(), options);
        }
    }
    else if(method == "unionfind")
    {
        vigra_precondition(srgType == CompleteGrow,
           "watersheds(): UnionFind only supports 'CompleteGrow' mode.");
           
        PyAllowThreads _pythread;
        if(neighborhood == 4)
        {
            maxRegionLabel = watershedsUnionFind(srcImageRange(image), destImage(res),
                                        FourNeighborCode());
        }
        else
        {
            maxRegionLabel = watershedsUnionFind(srcImageRange(image), destImage(res),
                                        EightNeighborCode());
        }
    }
    else
    {
        vigra_precondition(false, "watersheds(): Unknown watershed method requested.");
    }

    return python::make_tuple(res, maxRegionLabel);
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pywatersheds2D, pythonWatersheds2D)

template <unsigned int N, class PixelType >
python::tuple 
pythonWatershedsNew(NumpyArray<N, Singleband<PixelType> > image,
                    int neighborhood = 0,
                    NumpyArray<N, Singleband<npy_uint32> > seeds = NumpyArray<N, Singleband<npy_uint32> >(),
                    std::string method = "", 
                    SRGType srgType = CompleteGrow, 
                    PixelType max_cost = 0.0, 
                    NumpyArray<N, Singleband<npy_uint32> > res = NumpyArray<N, Singleband<npy_uint32> >())
{
    method = tolower(method);
    if(method == "" || method == "turbo")
    {
        method = "regiongrowing";
    }
    
    std::string description("watershed labeling, neighborhood=");
    description += asString(neighborhood);
    
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
            "watersheds(): Output array has wrong shape.");
    
    WatershedOptions options;
    options.srgType(srgType);
    
    if(method == "regiongrowing")
    {
        options.regionGrowing();
    }
    else if(method == "unionfind")
    {
        options.unionFind();
    }
    else
    {
        vigra_precondition(false, "watersheds(): Unknown watershed method requested.");
    }
    
    if(max_cost > 0.0)
    {
        vigra_precondition(method != "unionfind",
           "watersheds(): UnionFind does not support a cost threshold.");
        options.stopAtThreshold(max_cost);
    }
    
    if(seeds.hasData())
    {
        vigra_precondition(method != "unionfind",
           "watersheds(): UnionFind does not support seed images.");
        res = seeds;
    }
    else
    {
        options.seedOptions(SeedOptions().extendedMinima());
    }
    
    NeighborhoodType n = (neighborhood == 0)
                             ? DirectNeighborhood
                             : IndirectNeighborhood;    
    npy_uint32 maxRegionLabel = 0;
    {
        PyAllowThreads _pythread;
        maxRegionLabel = watershedsMultiArray(image, res, n, options);
    }

    return python::make_tuple(res, maxRegionLabel);
}

template <class PixelType >
python::tuple 
pythonWatersheds2DNew(NumpyArray<2, Singleband<PixelType> > image,
                      int neighborhood = 4,
                      NumpyArray<2, Singleband<npy_uint32> > seeds = NumpyArray<2, Singleband<npy_uint32> >(),
                      std::string method = "", 
                      SRGType srgType = CompleteGrow, 
                      PixelType max_cost = 0.0, 
                      NumpyArray<2, Singleband<npy_uint32> > res = NumpyArray<2, Singleband<npy_uint32> >())
{
    vigra_precondition(neighborhood == 4 || neighborhood == 8,
           "watersheds2D(): neighborhood must be 4 or 8.");
    neighborhood = (neighborhood == 4)
                        ? 0
                        : 1;
    return pythonWatershedsNew(image, neighborhood, seeds, method, srgType, max_cost, res);
}

template <class PixelType >
python::tuple 
pythonWatersheds3DNew(NumpyArray<3, Singleband<PixelType> > image,
                      int neighborhood = 6,
                      NumpyArray<3, Singleband<npy_uint32> > seeds = NumpyArray<3, Singleband<npy_uint32> >(),
                      std::string method = "", 
                      SRGType srgType = CompleteGrow, 
                      PixelType max_cost = 0.0, 
                      NumpyArray<3, Singleband<npy_uint32> > res = NumpyArray<3, Singleband<npy_uint32> >())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
           "watersheds3D(): neighborhood must be 6 or 26.");
    neighborhood = (neighborhood == 6)
                        ? 0
                        : 1;
    return pythonWatershedsNew(image, neighborhood, seeds, method, srgType, max_cost, res);
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pywatersheds2DNew, pythonWatersheds2DNew)
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pywatersheds3DNew, pythonWatersheds3DNew)

template < class PixelType >
python::tuple 
pythonWatersheds3D(NumpyArray<3, Singleband<PixelType> > image,
                   int neighborhood = 6,
                   NumpyArray<3, Singleband<npy_uint32> > seeds = NumpyArray<3, Singleband<npy_uint32> >(),
                   std::string method = "RegionGrowing", 
                   SRGType srgType = CompleteGrow, 
                   PixelType max_cost = 0.0, 
                   NumpyArray<3, Singleband<npy_uint32> > res = NumpyArray<3,Singleband<npy_uint32> >())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
           "watersheds3D(): neighborhood must be 6 or 26.");

    method = tolower(method);
    
    bool haveSeeds = seeds.hasData();
    unsigned int maxRegionLabel;
    
    if(method == "")
    {
        if(IsSameType<PixelType, npy_uint8>::value)
            method = "turbo";
        else
            method = "regiongrowing";
    }
    
    if(method == "turbo")
    {
        vigra_precondition((Or<typename IsSameType<PixelType, npy_uint8>::type,
                               typename IsSameType<PixelType, float>::type>::value),
           "watersheds3D(): Turbo algorithm requires input dtype = uint8 or dtype = float.");
        vigra_precondition(neighborhood == 6,
           "watersheds3D(): Turbo algorithm requires neighborhood = 6.");
        vigra_precondition(srgType == CompleteGrow,
           "watersheds3D(): Turbo algorithm requires termination = CompleteGrow.");
        vigra_precondition(max_cost == 0,
           "watersheds3D(): Turbo algorithm doesn't support 'max_cost'.");
    }
    
    if(method == "regiongrowing" || method == "turbo")
    {
        std::string description("watershed seeds");
        
        seeds.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
                "watersheds(): Seed array has wrong shape.");
        
        if(!haveSeeds)
        {
            PyAllowThreads _pythread;
            maxRegionLabel = 0;
            
            MultiArray<3, npy_uint32> minima(seeds.shape());
            
            if (neighborhood ==6)
            {
                extendedLocalMinima3D(srcMultiArrayRange(image), destMultiArray(minima),
                                      (npy_uint32)1, NeighborCode3DSix());
                maxRegionLabel = labelVolumeWithBackground(srcMultiArrayRange(minima),
                                          destMultiArray(seeds), NeighborCode3DSix(),
                                          (npy_uint32)0);
            }
            else
            {
                extendedLocalMinima3D(srcMultiArrayRange(image), destMultiArray(minima),
                                      (npy_uint32)1, NeighborCode3DTwentySix());
                maxRegionLabel = labelVolumeWithBackground(srcMultiArrayRange(minima),
                                          destMultiArray(seeds), NeighborCode3DTwentySix(),
                                          (npy_uint32)0);
            }
        }
        else
        {
            PyAllowThreads _pythread;
            FindMinMax< npy_uint32 > minmax;
            inspectMultiArray(srcMultiArrayRange(seeds), minmax);
            maxRegionLabel = minmax.max;
        }
           
        description = "watershed labeling, neighborhood=";
        description += asString(neighborhood);
        
        res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
                "watersheds(): Output array has wrong shape.");

        PyAllowThreads _pythread;
        ArrayOfRegionStatistics< SeedRgDirectValueFunctor< PixelType > > stats(maxRegionLabel);
        if(neighborhood == 6)
        {
            if(method == "turbo")
            {
                res = seeds;
                
                TWS<PixelType>::exec(image, res);
            }
            else
            {
                seededRegionGrowing3D(srcMultiArrayRange(image), srcMultiArray(seeds), 
                                      destMultiArray(res), 
                                      stats, srgType, NeighborCode3DSix(), max_cost);
            }
        }
        else
        {
            seededRegionGrowing3D(srcMultiArrayRange(image), srcMultiArray(seeds), 
                                  destMultiArray(res), 
                                  stats, srgType, NeighborCode3DTwentySix(), max_cost);
        }
    }
    else if(method == "unionfind")
    {
        vigra_precondition(!haveSeeds,
           "watersheds(): UnionFind does not support seed images.");
        vigra_precondition(srgType == CompleteGrow,
           "watersheds(): UnionFind only supports 'CompleteGrow' mode.");
           
        std::string description("watershed labeling, neighborhood=");
        description += asString(neighborhood);
        
        res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description), 
                "watersheds(): Output array has wrong shape.");
        
        PyAllowThreads _pythread;
        if(neighborhood == 6)
        {
            maxRegionLabel = watersheds3DSix(srcMultiArrayRange(image), destMultiArray(res));
        }
        else
        {
            maxRegionLabel = watersheds3DTwentySix(srcMultiArrayRange(image), destMultiArray(res));
        }
    }
    else
    {
        vigra_precondition(false, "watersheds(): Unknown watershed method requested.");
    }

    return python::make_tuple(res, maxRegionLabel);
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pywatersheds3D, pythonWatersheds3D)

template <unsigned int N, class PixelType >
python::tuple 
pythonSlic(NumpyArray<N, PixelType > array,
           double intensityScaling,
           unsigned int seedDistance,
           unsigned int minSize = 0,            // choose minSize automatically 
           unsigned int iterations = 10, 
           NumpyArray<N, Singleband<npy_uint32> > res = NumpyArray<N, Singleband<npy_uint32> >())
{
    typedef typename detail::ResolveMultiband<PixelType>::type ValueType;
    typedef typename NormTraits<ValueType>::NormType TmpType;
    typedef typename NumpyArray<N, ValueType >::view_type View;
    
    std::string description("Slic superpixels");
    
    res.reshapeIfEmpty(array.taggedShape().setChannelDescription(description), 
            "slicSuperpixels(): Output array has wrong shape.");
    
    npy_uint32 maxRegionLabel = 0;
    {
        PyAllowThreads _pythread;
        
        MultiArray<N, TmpType> gradMag(array.shape());
        
        // the original code uses the symmetric difference instead of a Gaussian gradient
        gaussianGradientMagnitude(array, gradMag, 1.0);
        // search radius of 1 is also used in the original code
        generateSlicSeeds(gradMag, res, seedDistance, 1);
        
        maxRegionLabel = slicSuperpixels(array, res, intensityScaling, seedDistance,
                                         SlicOptions().iterations(iterations)
                                                      .minSize(minSize));
    }

    return python::make_tuple(res, maxRegionLabel);
}

template <class PixelType >
python::tuple 
pythonSlic2D(NumpyArray<2, PixelType > image,
             double intensityScaling,
             unsigned int seedDistance,
             unsigned int minSize = 0,            // choose minSize automatically 
             unsigned int iterations = 10, 
             NumpyArray<2, Singleband<npy_uint32> > res = NumpyArray<2, Singleband<npy_uint32> >())
{
    return pythonSlic(image, intensityScaling, seedDistance, minSize, iterations, res);
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pySlic2D, pythonSlic2D)

template <class PixelType >
python::tuple 
pythonSlic3D(NumpyArray<3, PixelType > image,
             double intensityScaling,
             unsigned int seedDistance,
             unsigned int minSize = 0,            // choose minSize automatically 
             unsigned int iterations = 10, 
             NumpyArray<3, Singleband<npy_uint32> > res = NumpyArray<3, Singleband<npy_uint32> >())
{
    return pythonSlic(image, intensityScaling, seedDistance, minSize, iterations, res);
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pySlic3D, pythonSlic3D)

void defineSegmentation()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);

    multidef("labelImage", pyLabelImage<npy_uint8, float>(),
        (arg("image"), 
        arg("neighborhood") = 4,
        arg("out")=python::object()),
        "Find the connected components of a segmented image. Parameter 'neighborhood' specifies "
        "the pixel neighborhood to be used and can be 4 (default) or 8.\n\n"
        "For details see labelImage_ in the vigra C++ documentation.\n");

    multidef("labelImageWithBackground", pyLabelImageWithBackground<npy_uint8, float>(),
        (arg("image"), 
        arg("neighborhood") = 4,
        arg("background_value") = 0,
        arg("out")=python::object()),
        "Find the connected components of a segmented image, excluding the "
        "background from labeling, where the background is the set of all pixels with "
        "the given 'background_value'. Parameter 'neighborhood' specifies "
        "the pixel neighborhood to be used and can be 4 (default) or 8.\n\n"
        "For details see labelImageWithBackground_ in the vigra C++ documentation.\n");

    multidef("labelVolume", pyLabelVolume<npy_uint8, npy_uint32, float>(),
        (arg("volume"), 
        arg("neighborhood")=6,
        arg("out")=python::object()),
        "Find the connected components of a segmented volume. Parameter 'neighborhood' specifies "
        "the pixel neighborhood to be used and can be 6 (default) or 26.\n"
        "\n"
        "For details see labelVolume_ in the vigra C++ documentation.\n");

    multidef("labelVolumeWithBackground", pyLabelVolumeWithBackground<npy_uint8, npy_uint32, float>(),
        (arg("volume"), 
         arg("neighborhood")=6, 
         arg("background_value")=0,
         arg("out")=python::object()),
        "Find the connected components of a segmented volume, excluding the "
        "background from labeling, where the background is the set of all pixels with "
        "the given 'background_value'. Parameter 'neighborhood' specifies "
        "the pixel neighborhood to be used and can be 6 (default) or 26.\n"
        "\n"
        "For details see labelVolumeWithBackground_ in the vigra C++ documentation.\n");
    
    /******************************************************************************/
    
    def("localMinima",
        registerConverters(&pythonLocalMinima2D<float>),
        (arg("image"), 
         arg("marker")=1.0, 
         arg("neighborhood") = 8,
         arg("out")=python::object()),
        "Find local minima in an image and mark them with the given 'marker'. Parameter "
        "'neighborhood' specifies the pixel neighborhood to be used and can be "
        "4 or 8 (default).\n\n"
        "For details see localMinima_ in the vigra C++ documentation.\n");

    def("localMinima3D",
            registerConverters(&pythonLocalMinima3D<float> ),
            (arg("volume"), arg("marker") = 1.0, arg("neighborhood") = 6, arg(
                    "out") = python::object()),
            "Find local minima in a volume and mark them with the given 'marker'. Parameter "
                "'neighborhood' specifies the pixel neighborhood to be used and can be "
                "6 or 26 (default).\n\n"
                "For details see localMinima3D_ in the vigra C++ documentation.\n");

    // def("extendedLocalMinima",
        // registerConverters(&pythonExtendedLocalMinima2D<float>),
        // (arg("image"), 
         // arg("marker")=1.0, 
         // arg("neighborhood") = 8,
         // arg("out")=python::object()),
        // "Find local minima and minimal plateaus in an image and mark them with "
        // "the given 'marker'. Parameter 'neighborhood' specifies the pixel "
        // "neighborhood to be used and can be 4 or 8 (default).\n\n"
        // "For details see extendedLocalMinima_ in the vigra C++ documentation.\n"
        // );

    multidef("extendedLocalMinima",
        pyExtendedLocalMinima2D<npy_uint8, float>(),
        (arg("image"), 
         arg("marker")=1.0, 
         arg("neighborhood") = 8,
         arg("out")=python::object()),
        "Find local minima and minimal plateaus in an image and mark them with "
        "the given 'marker'. Parameter 'neighborhood' specifies the pixel "
        "neighborhood to be used and can be 4 or 8 (default).\n\n"
        "For details see extendedLocalMinima_ in the vigra C++ documentation.\n"
        );

    multidef("extendedLocalMinima3D",
        pyExtendedLocalMinima3D<float, npy_uint8>(), 
        (arg("volume"), arg("marker") = 1, arg("neighborhood") = 6, arg("out") = python::object()),
        "Find local minima and minimal plateaus in a volume and mark them with "
        "the given 'marker'. Parameter 'neighborhood' specifies the pixel "
        "neighborhood to be used and can be 6(default) or 26 .\n\n"
        "For details see extendedLocalMinima3D_ in the vigra C++ documentation.\n");

    def("localMaxima",
        registerConverters(&pythonLocalMaxima2D<float>),
        (arg("image"), 
         arg("marker")=1.0, 
         arg("neighborhood") = 8,
         arg("out")=python::object()),
        "Find local maxima in an image and mark them with the given 'marker'. Parameter "
        "'neighborhood' specifies the pixel neighborhood to be used and can be "
        "4 or 8 (default).\n\n"
        "For details see localMaxima_ in the vigra C++ documentation.\n");

    def("localMaxima3D", registerConverters(&pythonLocalMaxima3D<float> ), 
         (arg("volume"), arg("marker") = 1.0, arg("neighborhood") = 6, arg("out") = python::object()),
            "Find local maxima and maximal plateaus in a volume and mark them with "
            "the given 'marker'. Parameter 'neighborhood' specifies the pixel "
            "neighborhood to be used and can be 6(default) or 26 .\n\n"
            "For details see localMaxima3D_ in the vigra C++ documentation.\n");

    def("extendedLocalMaxima",
        registerConverters(&pythonExtendedLocalMaxima2D<float>),
        (arg("image"), 
         arg("marker")=1.0, 
         arg("neighborhood") = 8,
         arg("out")=python::object()),
        "Find local maxima and maximal plateaus in an image and mark them with "
        "the given 'marker'. Parameter 'neighborhood' specifies the pixel "
        "neighborhood to be used and can be 4 or 8 (default).\n\n"
        "For details see extendedLocalMaxima_ in the vigra C++ documentation.\n");

    def("extendedLocalMaxima3D", 
        registerConverters(&pythonExtendedLocalMaxima3D<float> ), 
        (arg("volume"), arg("marker") = 1.0, arg("neighborhood") = 6, arg("out") = python::object()),
        "Find local maxima and maximal plateaus in a volume and mark them with "
        "the given 'marker'. Parameter 'neighborhood' specifies the pixel "
        "neighborhood to be used and can be 6 (default) or 26 .\n\n"
        "For details see extendedLocalMaxima3D_ in the vigra C++ documentation.\n");

    /*************************************************************************/

    enum_<vigra::SRGType>("SRGType")
        .value("CompleteGrow", vigra::CompleteGrow)
        .value("KeepContours", vigra::KeepContours)
        .value("StopAtThreshold", vigra::StopAtThreshold)
        ;

    /*  FIXME: int64 is unsupported by the C++ code (hard-coded int) */
    multidef("watersheds", pywatersheds2D< npy_uint8, float >(),
      (arg("image"), 
       arg("neighborhood") = 4, 
       arg("seeds")=python::object(), 
       arg("method")="",
       arg("terminate")=CompleteGrow,
       arg("max_cost")=0,
       arg("out")=python::object()),
        "Compute the watersheds of a 2D image.\n"
        "\n"
        "   watersheds(image, neighborhood=4, seeds = None, methods = 'RegionGrowing', \n"
        "              terminate=CompleteGrow, threshold=0, out = None) -> (labelimage, max_ragion_label)\n"
        "\n"
        "Parameters:\n\n"
        " image:\n"
        "    the image or volume containing the boundary indicator values "
        "    (high values = high edgeness, dtype=numpy.uint8 or numpy.float32).\n"
        " neighborhood:\n"
        "    the pixel neighborhood to be used. Feasible values depend on the "
        "    dimension and method:\n\n"
        "      2-dimensional data:\n"
        "        4 (default) or 8.\n"
        "      3-dimensional data:\n"
        "        6 (default) or 26\n\n"
        " seeds:\n"
        "    a label image specifying region seeds, only supported by methods 'RegionGrowing' and 'Turbo'"
        "    (with dtype=numpy.uint32).\n" 
        " method:\n"
        "    the algorithm to be used for watershed computation. Possible values:\n\n"
        "      'Turbo':\n"
        "        (default if input dtype == uint8) use fastSeededRegionGrowing_ or tws() respectively\n"
        "      'RegionGrowing':\n"
        "        (default if input dtype != uint8) use seededRegionGrowing_ or seededRegionGrowing3D_ respectively\n"
        "      'UnionFind:\n"
        "        use watershedsUnionFind_ or watersheds3D_ respectively\n\n"
        " terminate:\n"
        "    when to stop growing. Possible values:\n\n"
        "      CompleteGrow:\n"
        "        (default) grow until all pixels are assigned to a region\n"
        "      KeepCountours:\n"
        "        keep a 1-pixel wide contour between all regions, only supported "
        "        by method 'RegionGrowing'\n"
        "      StopAtThreshold:\n"
        "        stop when the boundary indicator values exceed the threshold given by "
        "        parameter 'max_cost', only supported by method 'RegionGrowing'\n"
        "      KeepCountours | StopAtThreshold:\n"
        "        keep 1-pixel wide contour and stop at given 'max_cost', only "
        "        supported by method 'RegionGrowing'\n\n"
        " max_cost:\n"
        "    terminate growing when boundary indicator exceeds this value (ignored when "
        "    'terminate' is not StopAtThreshold or method is not 'RegionGrowing')\n" 
        " out:\n"
        "    the label image (with dtype=numpy.uint32) to be filled by the algorithm. "
        "    It will be allocated by the watershed function if not provided)\n\n"
         "The function returns a Python tuple (labelImage, maxRegionLabel)\n\n"
         );

    multidef("watersheds", pywatersheds3D< npy_uint8, float >(),
      (arg("volume"), 
       arg("neighborhood") = 6, 
       arg("seeds")=python::object(), 
       arg("method")="",
       arg("terminate")=CompleteGrow,
       arg("max_cost")=0,
       arg("out")=python::object()),
       "Likewise, compute watersheds of a volume.\n");

    multidef("watershedsNew", pywatersheds2DNew< npy_uint8, float >(),
      (arg("image"), 
       arg("neighborhood") = 4, 
       arg("seeds")=python::object(), 
       arg("method")="",
       arg("terminate")=CompleteGrow,
       arg("max_cost")=0,
       arg("out")=python::object()),
       "graph-based watershed");

    multidef("watershedsNew", pywatersheds3DNew< npy_uint8, float >(),
      (arg("image"), 
       arg("neighborhood") = 6, 
       arg("seeds")=python::object(), 
       arg("method")="",
       arg("terminate")=CompleteGrow,
       arg("max_cost")=0,
       arg("out")=python::object()),
       "graph-based watershed");

    multidef("slicSuperpixels", pySlic2D< TinyVector<float, 3>, Singleband<float> >(),
      (arg("image"), 
       arg("intensityScaling"), 
       arg("seedDistance"), 
       arg("minSize")=0,
       arg("iterations")=10,
       arg("out")=python::object()),
        "Compute Slic superpixels for a 2D image.\n\n"

        "Parameters:\n\n"
        " image:\n"
        "    The 2D-image on which the superpixels will be calculated. Accepts single- and threeband images. \n\n"
        " intensityScaling:\n"
        "    Scale (divide) color/intensity difference by this parameter before comparing to spatial distance. \n\n"
        " seedDistance:\n"
        "    specify the radius of the window around each seed in which the algorithm looks for potential members of the corresponding superpixel"
        " thus limiting the superpixel size. The grid spacing for seed placement is determined by this parameter.\n\n"
        " minSize:\n"
        "    Minimum size for superpixels. By default the algorithm merges all regions smaller than a quarter of the average superpixel size.\n\n"
        " iterations:\n"
        "    Specify number of iterations. The default is 10."
        " out:\n"
        "    The label image (with dtype=numpy.uint32) to be filled by the algorithm. "
        "    It will be allocated by the slicSuperpixels function if not provided)\n\n"
        "The function returns a Python tuple (labelImage, maxRegionLabel)\n\n");

    multidef("slicSuperpixels", pySlic3D< TinyVector<float, 3>, Singleband<float> >(),
      (arg("image"), 
       arg("intensityScaling"), 
       arg("seedDistance"), 
       arg("minSize")=0,
       arg("iterations")=10,
       arg("out")=python::object()),
       "Likewise compute Slic superpixels for a 3D volume, either single- or threeband.\n");
}

void defineEdgedetection();
void defineInterestpoints();
void defineAccumulators();

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(analysis)
{
    import_vigranumpy();
    defineSegmentation();
    defineEdgedetection();
    defineInterestpoints();
    defineAccumulators();
}
