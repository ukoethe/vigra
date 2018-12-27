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
#include <vigra/blockwise_watersheds.hxx>
#include <vigra/seededregiongrowing.hxx>
#include <vigra/labelvolume.hxx>
#include <vigra/watersheds3d.hxx>
#include <vigra/seededregiongrowing3d.hxx>
#include <vigra/multi_watersheds.hxx>
#include <vigra/convolution.hxx>
#include <vigra/multi_convolution.hxx>
#include <vigra/slic.hxx>
#include <vigra/seg_to_seeds.hxx>
#include <vigra/multi_pointoperators.hxx>

#include <string>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

#include <boost/python/stl_iterator.hpp>

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

template < class VoxelType, unsigned int ndim >
NumpyAnyArray
pythonLabelMultiArray(NumpyArray<ndim, Singleband<VoxelType> > volume,
                      // std::string neighborhood="",
                      python::object neighborspec=python::object(),
                      NumpyArray<ndim, Singleband<npy_uint32> > res = NumpyArray<ndim, Singleband<npy_uint32> >())
{
    std::string neighborhood;
    if(neighborspec == python::object())
    {
        neighborhood = "direct";
    }
    else if(python::extract<int>(neighborspec).check())
    {
        int n = python::extract<int>(neighborspec)();
        if(n == 2*ndim || n == 0)
        {
            neighborhood = "direct";
        }
        else if(n == std::pow(3, ndim) - 1)
        {
            neighborhood = "indirect";
        }
    }
    else if(python::extract<std::string>(neighborspec).check())
    {
        neighborhood = tolower(python::extract<std::string>(neighborspec)());
        if (neighborhood == "")
        {
            neighborhood = "direct";
        }
    }

    vigra_precondition(neighborhood == "direct" || neighborhood == "indirect",
        "labelMultiArray(): neighborhood must be 'direct' or 'indirect' or '' (defaulting "
        "to 'direct') or the appropriate number of neighbors (4 or 8 in 2D, 6 or 26 in 3D).");

    std::string description("connected components, neighborhood=" + neighborhood);

    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description),
        "labelMultiArray(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        if (neighborhood == "direct")
        {
            labelMultiArray(volume, res, DirectNeighborhood);
        }
        else
        {
            labelMultiArray(volume, res, IndirectNeighborhood);
        }
    }
    return res;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR_NDIM(pyLabelMultiArray, pythonLabelMultiArray)

template < class VoxelType, unsigned int ndim >
NumpyAnyArray
pythonLabelMultiArrayWithBackground(NumpyArray<ndim, Singleband<VoxelType> > volume,
                                // std::string neighborhood="",
                                python::object neighborspec=python::object(),
                                VoxelType background_value = 0,
                                NumpyArray<ndim, Singleband<npy_uint32> > res = NumpyArray<ndim, Singleband<npy_uint32> >())
{
    std::string neighborhood;
    if(neighborspec == python::object())
    {
        neighborhood = "direct";
    }
    else if(python::extract<int>(neighborspec).check())
    {
        int n = python::extract<int>(neighborspec)();
        if(n == 2*ndim || n == 0)
        {
            neighborhood = "direct";
        }
        else if(n == std::pow(3, ndim) - 1)
        {
            neighborhood = "indirect";
        }
    }
    else if(python::extract<std::string>(neighborspec).check())
    {
        neighborhood = tolower(python::extract<std::string>(neighborspec)());
        if (neighborhood == "")
        {
            neighborhood = "direct";
        }
    }

    vigra_precondition(neighborhood == "direct" || neighborhood == "indirect",
        "labelMultiArrayWithBackground(): neighborhood must be 'direct' or 'indirect' or '' (defaulting "
        "to 'direct') or the appropriate number of neighbors (4 or 8 in 2D, 6 or 26 in 3D).");

    std::string description("connected components with background, neighborhood=");
    description += neighborhood + ", bglabel=" + asString(background_value);

    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description),
        "labelMultiArrayWithBackground(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        if (neighborhood == "direct")
        {
            labelMultiArrayWithBackground(volume,
                res, DirectNeighborhood,
                background_value);
        }
        else
        {
            labelMultiArrayWithBackground(volume,
                res, IndirectNeighborhood,
                background_value);
        }
    }
    return res;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR_NDIM(pyLabelMultiArrayWithBackground, pythonLabelMultiArrayWithBackground)

/*********************************************************************************/

// FIXME: support output of label images from localMinim/Maxima functions

template < class PixelType >
NumpyAnyArray
pythonLocalMinima2D(NumpyArray<2, Singleband<PixelType> > image,
                    PixelType marker = NumericTraits<PixelType>::one(),
                    int neighborhood = 8,
                    bool allowAtBorder = false,
                    bool allowPlateaus = false,
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
        localMinima(image, res,
            LocalMinmaxOptions()
                .neighborhood(neighborhood)
                .allowAtBorder(allowAtBorder)
                .markWith(marker)
                .allowPlateaus(allowPlateaus));
    }

    return res;
}

template<class PixelType>
NumpyAnyArray
pythonLocalMinima3D(NumpyArray<3, Singleband<PixelType> > volume,
                    PixelType marker = NumericTraits<PixelType>::one(),
                    int neighborhood = 6,
                    bool allowAtBorder = false,
                    bool allowPlateaus = false,
                    NumpyArray<3, Singleband<PixelType> > res = NumpyArray<3, Singleband<PixelType> >())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
            "localMinima(): neighborhood must be 6 or 26.");

    std::string description("local minima, neighborhood=");
    description += asString(neighborhood);

    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description),
            "localMinima(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        localMinima(volume, res,
            LocalMinmaxOptions()
                .neighborhood(neighborhood)
                .allowAtBorder(allowAtBorder)
                .markWith(marker)
                .allowPlateaus(allowPlateaus));
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
                    bool allowAtBorder = false,
                    bool allowPlateaus = false,
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
        localMaxima(image, res,
            LocalMinmaxOptions()
                .neighborhood(neighborhood)
                .allowAtBorder(allowAtBorder)
                .markWith(marker)
                .allowPlateaus(allowPlateaus));
    }

    return res;
}

template<class PixelType>
NumpyAnyArray
pythonLocalMaxima3D(NumpyArray<3, Singleband<PixelType> > volume,
                    PixelType marker = NumericTraits<PixelType>::one(),
                    int neighborhood = 6,
                    bool allowAtBorder = false,
                    bool allowPlateaus = false,
                    NumpyArray<3, Singleband<PixelType> > res = NumpyArray<3, Singleband<PixelType> >())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
            "localMaxima(): neighborhood must be 6 or 26.");

    std::string description("local maxima, neighborhood=");
    description += asString(neighborhood);

    res.reshapeIfEmpty(volume.taggedShape().setChannelDescription(description),
            "localMaxima(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        localMaxima(volume, res,
            LocalMinmaxOptions()
                .neighborhood(neighborhood)
                .allowAtBorder(allowAtBorder)
                .markWith(marker)
                .allowPlateaus(allowPlateaus));
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

template <class PixelType, int N>
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

VIGRA_PYTHON_MULTITYPE_FUNCTOR_NDIM(pySlic, pythonSlic)

template<unsigned int DIM>
NumpyAnyArray  pythonShrinkLabels(
    NumpyArray<DIM,npy_uint32> labels,
    const size_t shrinkNpixels,
    NumpyArray<DIM,Singleband<npy_uint32> > out = NumpyArray<DIM,Singleband<npy_uint32> >()
){
    out.reshapeIfEmpty(labels.shape());
    shrinkLabels(labels,shrinkNpixels,out);
    return out;
}


template<class T>
vigra::NumpyAnyArray pySizeFilterSegInplace(vigra::NumpyArray<3, T>  seg, const vigra::UInt32 maxLabel, const vigra::UInt32 sizeLimit, bool checkAtBorder=false){


    std::vector<bool > atBorder(maxLabel+1, false);

    if (! checkAtBorder){
        for(std::ptrdiff_t z=0;z<seg.shape(2); ++z)
        for(std::ptrdiff_t y=0;y<seg.shape(1); ++y){
            atBorder[seg(0,y,z)] = true;
            atBorder[seg(seg.shape(0)-1,y,z)] = true;
        }

        for(std::ptrdiff_t z=0;z<seg.shape(2); ++z)
        for(std::ptrdiff_t x=0;x<seg.shape(0); ++x){
            atBorder[seg(x,0,z)] = true;
            atBorder[seg(x,seg.shape(1)-1,z)] = true;
        }

        for(std::ptrdiff_t y=0;y<seg.shape(1); ++y)
        for(std::ptrdiff_t x=0;x<seg.shape(0); ++x){
            atBorder[seg(x,y,0)] = true;
            atBorder[seg(x,y,seg.shape(2)-1)] = true;
        }
    }



    std::vector<size_t > counts(maxLabel+1,0);

    for(auto iter = seg.begin(); iter!=seg.end(); ++iter){
        counts[*iter] += 1;
    }



    for(auto iter = seg.begin(); iter!=seg.end(); ++iter){
        const auto l = *iter;
        const auto c = counts[l];
        if(c<sizeLimit && atBorder[l] == false){
            *iter = 0;
        }
    }

    return seg;
}


template<unsigned int DIM>
python::tuple  pyUnionFindWatershedsBlockwise(
    NumpyArray<DIM,float> data,
    TinyVector<Int64, DIM> blockShape,
    NumpyArray<DIM, UInt32 > out
){
    out.reshapeIfEmpty(data.shape());
    UInt64 nSeg =  unionFindWatershedsBlockwise(data, out,
                                                BlockwiseLabelOptions().neighborhood(DirectNeighborhood)
                                                                       .blockShape(blockShape));
    return python::make_tuple(out, nSeg);
}

/** \brief Map all values in src to new values using the given mapping (a dict).
 *  See python docstring for details.
*/
template <unsigned int NDIM, class SrcVoxelType, class DestVoxelType>
NumpyAnyArray
pythonApplyMapping(NumpyArray<NDIM, Singleband<SrcVoxelType> > src,
                   python::dict mapping,
                   bool allow_incomplete_mapping = false,
                   NumpyArray<NDIM, Singleband<DestVoxelType> > res = NumpyArray<NDIM, Singleband<SrcVoxelType> >())
{
    using namespace boost::python;

    res.reshapeIfEmpty(src.taggedShape(), "applyMapping(): Output array has wrong shape.");

    // Copy dict into a c++ unordered_map of ints,
    // which is ~10x faster than using a Python dict
    typedef std::unordered_map<SrcVoxelType, DestVoxelType> labelmap_t;
    labelmap_t labelmap(2*len(mapping)); // Using 2*N buckets seems to speed things up by 10%

    typedef stl_input_iterator<tuple> dict_iter_t;

#if PY_MAJOR_VERSION < 3
    dict_iter_t map_iter = mapping.iteritems();
#else
    dict_iter_t map_iter = mapping.items();
#endif

    for (; map_iter != dict_iter_t(); ++map_iter)
    {
        object key = (*map_iter)[0];
        object value = (*map_iter)[1];
        labelmap[extract<SrcVoxelType>(key)] = extract<DestVoxelType>(value);
    }

    // Enforce const capture in the lambda below.
    labelmap_t const & _labelmap = labelmap;

    {
        std::unique_ptr<PyAllowThreads> pythread_ptr(new PyAllowThreads);

        transformMultiArray(src, res,
            [&_labelmap, allow_incomplete_mapping, &pythread_ptr](SrcVoxelType px) -> DestVoxelType {
                typename labelmap_t::const_iterator iter = _labelmap.find(px);

                if (iter != _labelmap.end())
                {
                    return iter->second;
                }

                if (allow_incomplete_mapping)
                {
                    // Key is missing. Return the original value.
                    return static_cast<DestVoxelType>(px);
                }

                // Reclaim the GIL before setting the error string.
                pythread_ptr.reset();

                std::ostringstream err_msg;
                err_msg << "Key not found in mapping: " << +px;
                PyErr_SetString( PyExc_KeyError, err_msg.str().c_str() );
                python::throw_error_already_set();

                return 0; // unreachable line
            });
    }

    return res;
}

// Unfortunately, can't use this macro because the template args uses TWO dtypes
//VIGRA_PYTHON_MULTITYPE_FUNCTOR_NDIM(pyApplyMapping, pythonApplyMapping)


/** \brief Find unique values in the given array.
*/
template <class VoxelType, unsigned int NDIM>
NumpyAnyArray
pythonUnique(NumpyArray<NDIM, Singleband<VoxelType> > src, bool sort=true)
{
    std::unordered_set<VoxelType> labelset;
    auto f = [&labelset](VoxelType px) { labelset.insert(px); };
    inspectMultiArray(src, f);

    NumpyArray<1, VoxelType> result;
    result.reshape( Shape1(labelset.size()) );
    std::copy( labelset.begin(), labelset.end(), result.begin() );

    if (sort)
    {
        std::sort( result.begin(), result.end() );
    }
    return result;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR_NDIM(pyUnique, pythonUnique)


/** \brief Relabel an array such that all labels are consecutive
 * (i.e. there are no gaps in the label values used by the array)
 * See python docstring below for details.
*/
template <unsigned int NDIM, class SrcVoxelType, class DestVoxelType>
boost::python::tuple
pythonRelabelConsecutive(NumpyArray<NDIM, Singleband<SrcVoxelType> > src,
                         DestVoxelType start_label = 1,
                         bool keep_zeros = true,
                         NumpyArray<NDIM, Singleband<DestVoxelType> > res = NumpyArray<NDIM, Singleband<SrcVoxelType> >())
{
    using namespace boost::python;
    res.reshapeIfEmpty(src.taggedShape(), "relabelConsecutive(): Output array has wrong shape.");

    std::unordered_map<SrcVoxelType, DestVoxelType> labelmap;
    if (keep_zeros)
    {
        vigra_precondition(!keep_zeros || start_label > 0,
            "relabelConsecutive(): start_label must be non-zero if using keep_zeros=True");

        // pre-initialize the mapping to keep zeros unchanged
        labelmap[0] = 0;
    }

    {
        PyAllowThreads _pythread;

        transformMultiArray(src, res,
            [&](SrcVoxelType px) -> DestVoxelType {
                auto iter = labelmap.find(px);
                if (iter != labelmap.end())
                {
                    return iter->second;
                }
                // We haven't seen this label yet.
                // Create a new entry in the hash table.
                DestVoxelType newlabel = labelmap.size() - int(keep_zeros) + start_label;
                labelmap[px] = newlabel;
                return newlabel;
            });
    }

    // Convert labelmap to dict
    dict labelmap_dict;
    for (auto old_new_pair : labelmap)
    {
        labelmap_dict[old_new_pair.first] = old_new_pair.second;
    }

    DestVoxelType max_label = labelmap.size() - int(keep_zeros) - 1 + start_label;
    return make_tuple(res, max_label, labelmap_dict);
}

// Unfortunately, can't use this macro because the template args uses TWO dtypes
//VIGRA_PYTHON_MULTITYPE_FUNCTOR_NDIM(pyRelabelConsecutive, pythonRelabelConsecutive)


void defineSegmentation()
{
    using namespace python;

    docstring_options doc_options(true, true, false);


    python::def("unionFindWatershed3D",
        registerConverters(&pyUnionFindWatershedsBlockwise<3>),
        (
            python::arg("image"),
            python::arg("blockShape"),
            python::arg("out") = python::object()
        )
    );

    python::def("segToSeeds", registerConverters(pythonShrinkLabels<2>),
        (
            python::arg("image"),
            python::arg("shrinkN"),
            python::arg("out")=python::object()
        ),
        "shrink / ungrow a labeling / segmentation"
    );

    python::def("segToSeeds", registerConverters(pythonShrinkLabels<3>),
        (
            python::arg("image"),
            python::arg("shrinkN"),
            python::arg("out")=python::object()
        ),
        "shrink / ungrow a labeling / segmentation"
    );

    multidef("labelImage",
        pyLabelMultiArray<2, 2, npy_uint8, npy_uint32, float>().installFallback(),
        (arg("image"),
         arg("neighborhood") = 4,
         arg("out")=python::object()),
        "Find the connected components of a segmented image. Parameter 'neighborhood' specifies "
        "the pixel neighborhood to be used and can be 4 (default) or 8.\n\n"
        "For details see labelMultiArray_ in the vigra C++ documentation.\n");

    multidef("labelImageWithBackground",
        pyLabelMultiArrayWithBackground<2, 2, npy_uint8, npy_uint32, float>().installFallback(),
        (arg("image"),
         arg("neighborhood") = 4,
         arg("background_value") = 0,
         arg("out")=python::object()),
        "Find the connected components of a segmented image, excluding the "
        "background from labeling, where the background is the set of all pixels with "
        "the given 'background_value'. Parameter 'neighborhood' specifies "
        "the pixel neighborhood to be used and can be 4 (default) or 8.\n\n"
        "For details see labelMultiArrayWithBackground_ in the vigra C++ documentation.\n");

    multidef("labelVolume",
        pyLabelMultiArray<3, 3, npy_uint8, npy_uint32, float>().installFallback(),
        (arg("volume"),
         arg("neighborhood")=6,
         arg("out")=python::object()),
        "Find the connected components of a segmented volume. Parameter 'neighborhood' specifies "
        "the pixel neighborhood to be used and can be 6 (default) or 26.\n"
        "\n"
        "For details see labelMultiArray_ in the vigra C++ documentation.\n");

    multidef("labelVolumeWithBackground",
        pyLabelMultiArrayWithBackground<3, 3, npy_uint8, npy_uint32, float>().installFallback(),
        (arg("volume"),
         arg("neighborhood")=6,
         arg("background_value")=0,
         arg("out")=python::object()),
        "Find the connected components of a segmented volume, excluding the "
        "background from labeling, where the background is the set of all pixels with "
        "the given 'background_value'. Parameter 'neighborhood' specifies "
        "the pixel neighborhood to be used and can be 6 (default) or 26.\n"
        "\n"
        "For details see labelMultiArrayWithBackground_ in the vigra C++ documentation.\n");

    multidef("labelMultiArray",
        pyLabelMultiArray<2, 5, npy_uint8, npy_uint32, float>().installFallback(),
        (arg("array"),
         arg("neighborhood")="",
         arg("out")=python::object()),
        "Find the connected components of a segmented multi-dimensional array\n"
        "(supported dimensions: 2 to 5).\n"
        "Parameter 'neighborhood' specifies the pixel neighborhood to be used\n"
        "and can be 'direct' (default) or 'indirect' or the exact number of\n"
        "neighbors (2D: 4 or 8, 3D: 6 or 26, 4D: 8 or 80, 5D: 10 or 242).\n"
        "\n"
        "For details see labelMultiArray_ in the vigra C++ documentation.\n");

    multidef("labelMultiArrayWithBackground",
        pyLabelMultiArrayWithBackground<2, 5, npy_uint8, npy_uint32, float>().installFallback(),
        (arg("array"),
         arg("neighborhood")="",
         arg("background_value")=0,
         arg("out")=python::object()),
        "Find the connected components of a segmented multi-dimensional array\n"
        "(supported dimensions: 2 to 5), excluding the background from labeling,\n"
        "where background is the set of all pixels with the given 'background_value'.\n"
        "Parameter 'neighborhood' specifies the pixel neighborhood to be used\n"
        "and can be 'direct' (default) or 'indirect' or the exact number of\n"
        "neighbors (2D: 4 or 8, 3D: 6 or 26, 4D: 8 or 80, 5D: 10 or 242).\n"
        "\n"
        "For details see labelMultiArrayWithBackground_ in the vigra C++ documentation.\n");

    def("sizeFilterSegInplace",registerConverters(&pySizeFilterSegInplace<UInt32>),
        (arg("seg"),
         arg("maxLabel"),
         arg("sizeLimit"),
         arg("checkAtBorder") = false),
        "replace every occurrence of each number in the array 'seg' with zeros if this number"
        " occures less than 'sizeLimit' times in the array. If 'checkAtBorder' is false (default) "
        "segments that touch the border of the array will not be changed.\n"
        "'maxLabel' is the maximum label in seg\n"
    );

    /******************************************************************************/

    def("localMinima",
        registerConverters(&pythonLocalMinima2D<float>),
        (arg("image"),
         arg("marker")=1.0,
         arg("neighborhood") = 8,
         arg("allowAtBorder") = false,
         arg("allowPlateaus") = false,
         arg("out")=python::object()),
        "Find local minima in an image and mark them with the given 'marker'. Parameter "
        "'neighborhood' specifies the pixel neighborhood to be used and can be "
        "4 or 8 (default). \n"
        "If 'allowAtBorder' is true local minima at the image border will be detected.\n"
        "If 'allowPlateaus' is true regions of constant gray value whose neighbors are all higher than the value of the region will be detected."
        "\n\n"
        "For details see localMinima_ in the vigra C++ documentation.\n");

    def("localMinima3D",
        registerConverters(&pythonLocalMinima3D<float> ),
        (arg("volume"),
         arg("marker") = 1.0,
         arg("neighborhood") = 6,
         arg("allowAtBorder") = false,
         arg("allowPlateaus") = false,
         arg("out") = python::object()),
        "Find local minima in a volume and mark them with the given 'marker'. Parameter "
        "6 (default) or 26.\n"
        "If 'allowAtBorder' is set to 'True' local minima at the volume border will be detected.\n"
        "If 'allowPlateaus' is set to 'True' regions of constant gray value whose neighbors are all higher than the value of the region will be detected."
        "\n\n"
        "For details see localMinima_ in the vigra C++ documentation.\n");

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
        pyExtendedLocalMinima2D<npy_uint8, float>().installFallback(),
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
        pyExtendedLocalMinima3D<float, npy_uint8>().installFallback(),
        (arg("volume"),
         arg("marker") = 1,
         arg("neighborhood") = 6,
         arg("out") = python::object()),
        "Find local minima and minimal plateaus in a volume and mark them with "
        "the given 'marker'. Parameter 'neighborhood' specifies the pixel "
        "neighborhood to be used and can be 6(default) or 26 .\n\n"
        "For details see extendedLocalMinima3D_ in the vigra C++ documentation.\n");

    def("localMaxima",
        registerConverters(&pythonLocalMaxima2D<float>),
        (arg("image"),
         arg("marker")=1.0,
         arg("neighborhood") = 8,
         arg("allowAtBorder") = false,
         arg("allowPlateaus") = false,
         arg("out")=python::object()),
        "Find local maxima in an image and mark them with the given 'marker'. Parameter "
        "'neighborhood' specifies the pixel neighborhood to be used and can be "
        "4 or 8 (default).\n"
        "If 'allowAtBorder' is set to 'True' local maxima at image border will be detected.\n"
        "If 'allowPlateaus' is set to 'True' regions of constant gray value whose neighbors are all lower than the value of the region will be detected."
        "\n\n"
        "For details see localMaxima_ in the vigra C++ documentation.\n");

    def("localMaxima3D", registerConverters(&pythonLocalMaxima3D<float> ),
        (arg("volume"),
         arg("marker") = 1.0,
         arg("neighborhood") = 6,
         arg("allowAtBorder") = false,
         arg("allowPlateaus") = false,
         arg("out") = python::object()),
        "Find local maxima and maximal plateaus in a volume and mark them with "
        "the given 'marker'. Parameter 'neighborhood' specifies the pixel "
        "neighborhood to be used and can be 6(default) or 26.\n"
        "If 'allowAtBorder' is set to 'True' local maxima at the volume border will be detected.\n"
        "If 'allowPlateaus' is set to 'True' regions of constant gray value whose neighbors are all lower than the value of the region will be detected."
        "\n\n"
        "For details see localMaxima_ in the vigra C++ documentation.\n");

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
    multidef("watersheds",
        pywatersheds2D< npy_uint8, float >().installFallback().noPythonSignature(),
        (arg("image"),
         arg("neighborhood") = 4,
         arg("seeds")=python::object(),
         arg("method")="",
         arg("terminate")=CompleteGrow,
         arg("max_cost")=0,
         arg("out")=python::object()),
        "\n"
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
        "        (default if input dtype == uint8) use fastSeededRegionGrowing() (in 2D) or tws() (in 3D)\n"
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

    multidef("watersheds",
        pywatersheds3D< npy_uint8, float >().noPythonSignature(),
        (arg("volume"),
         arg("neighborhood") = 6,
         arg("seeds")=python::object(),
         arg("method")="",
         arg("terminate")=CompleteGrow,
         arg("max_cost")=0,
         arg("out")=python::object()),
        "Likewise, compute watersheds of a volume.\n");

    multidef("watershedsNew",
        pywatersheds2DNew< npy_uint8, float >().installFallback(),
        (arg("image"),
         arg("neighborhood") = 4,
         arg("seeds")=python::object(),
         arg("method")="",
         arg("terminate")=CompleteGrow,
         arg("max_cost")=0,
         arg("out")=python::object()),
        "graph-based watershed");

    multidef("watershedsNew",
        pywatersheds3DNew< npy_uint8, float >(),
        (arg("image"),
         arg("neighborhood") = 6,
         arg("seeds")=python::object(),
         arg("method")="",
         arg("terminate")=CompleteGrow,
         arg("max_cost")=0,
         arg("out")=python::object()),
       "graph-based watershed");

    multidef("slicSuperpixels",
        pySlic<2, 3, TinyVector<float, 3>, Singleband<float> >().installFallback(),
        (arg("array"),
         arg("intensityScaling"),
         arg("seedDistance"),
         arg("minSize")=0,
         arg("iterations")=10,
         arg("out")=python::object()),
        "\n"
        "Compute Slic superpixels for a 2D or 3D image.\n"
        "\n"
        "Parameters:\n\n"
        " array:\n"
        "    The array on which the superpixels will be calculated. Accepts single- and\n"
        "    threeband images/volumes. \n"
        "\n"
        " intensityScaling:\n"
        "    Scale (divide) color/intensity difference by this parameter before comparing\n"
        "    to spatial distance. \n"
        "\n"
        " seedDistance:\n"
        "    specify the radius of the window around each seed in which the algorithm looks\n"
        "    for potential members of the corresponding superpixel thus limiting the\n"
        "    superpixel size. The grid spacing for seed placement is determined by this parameter.\n"
        "\n"
        " minSize:\n"
        "    Minimum size for superpixels. By default the algorithm merges all regions smaller\n"
        "    than a quarter of the average superpixel size.\n"
        "\n"
        " iterations:\n"
        "    Specify number of iterations. The default is 10.\n"
        "\n"
        " out:\n"
        "    The label image (with dtype=numpy.uint32) to be filled by the algorithm. "
        "    It will be allocated by the slicSuperpixels function if not provided)\n"
        "\n"
        "The function returns a Python tuple (labelImage, maxRegionLabel)\n"
        "\n");

    multidef("unique",
        pyUnique<1,5,npy_uint8, npy_uint32, npy_uint64, npy_int64>(),
        (arg("arr"), arg("sort")=true),
        "Find unique values in the given label array.\n"
        "If ``sort`` is True, then the output is sorted.\n"
        "Much faster then ``numpy.unique()``.\n");

    //-- 3D relabelConsecutive
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<3, npy_uint64, npy_uint32>),
        (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()),
        "Relabel the given label image to have consecutive label values.\n"
        "Note: The relative order between label values will not necessarily be preserved.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "labels: ndarray\n"
        "start_label: The lowest label of the output array.\n"
        "keep_zeros: Don't relabel zero-valued items.\n"
        "out: ndarray to hold the data. If None, it will be allocated for you.\n"
        "     A combination of uint64 labels and uint32 'out' is permitted.\n"
        "\n"
        "Returns a tuple of ``(newlabels, maxlabel, mapping)``, where:\n"
        "``maxlabel`` is the maximum label of the new labels, and\n"
        "``mapping`` is a dict showing how the old labels were converted to the new label values.\n"
        "\n"
        "Note: As with other vigra functions, you should provide accurate axistags for optimal performance.\n");
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<3, npy_uint64, npy_uint64>), (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<3, npy_uint32, npy_uint32>), (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<3, npy_uint8, npy_uint8>),   (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));

    //-- 2D relabelConsecutive
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<2, npy_uint64, npy_uint32>), (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<2, npy_uint64, npy_uint64>), (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<2, npy_uint32, npy_uint32>), (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<2, npy_uint8, npy_uint8>),   (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));

    //-- 1D relabelConsecutive
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<1, npy_uint64, npy_uint32>), (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<1, npy_uint64, npy_uint64>), (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<1, npy_uint32, npy_uint32>), (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));
    def("relabelConsecutive", registerConverters(&pythonRelabelConsecutive<1, npy_uint8, npy_uint8>),   (arg("labels"), arg("start_label")=1, arg("keep_zeros")=true, arg("out")=python::object()));


    // Lots of overloads here to allow mapping between arrays of different dtypes.
    // -- 3D
    // 64 --> 32
    def("applyMapping", registerConverters(&pythonApplyMapping<3, npy_uint64, npy_uint32>),
        (arg("labels"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()),
        "Map all values in `labels` to new values using the given mapping (a dict).\n"
        "Useful for maps with large values, for which a numpy index array would need too much RAM.\n"
        "To relabel in-place, set `out=labels`.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "labels: ndarray\n"
        "mapping: dict of ``{old_label : new_label}``\n"
        "allow_incomplete_mapping: If True, then any voxel values in the original data that are missing\n"
        "                          from the mapping dict will be copied (and casted) into the output.\n"
        "                          Otherwise, an ``IndexError`` will be raised if the map is incomplete\n"
        "                          for the input data.\n"
        "out: ndarray to hold the data. If None, it will be allocated for you.\n"
        "     The dtype of ``out`` is allowed to be smaller (or bigger) than the dtype of ``labels``.\n"
        "\n"
        "Note: As with other vigra functions, you should provide accurate axistags for optimal performance.\n");

    // 8 <--> 32
    def("applyMapping", registerConverters(&pythonApplyMapping<3, npy_uint8, npy_uint32>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<3, npy_uint32, npy_uint8>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // 32 <--> 64
    def("applyMapping", registerConverters(&pythonApplyMapping<3, npy_uint32, npy_uint64>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    //def("applyMapping", registerConverters(&pythonApplyMapping<3, npy_uint64, npy_uint32>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // 8 <--> 64
    def("applyMapping", registerConverters(&pythonApplyMapping<3, npy_uint8, npy_uint64>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<3, npy_uint64, npy_uint8>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // Cases for same input/output dtypes must come last, so they are chosen by default!
    // 8 <--> 8, 32 <--> 32, 64 <--> 64
    def("applyMapping", registerConverters(&pythonApplyMapping<3, npy_uint64, npy_uint64>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<3, npy_uint32, npy_uint32>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<3, npy_uint8, npy_uint8>),   (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // -- 2D
    // 8 <--> 32
    def("applyMapping", registerConverters(&pythonApplyMapping<2, npy_uint8, npy_uint32>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<2, npy_uint32, npy_uint8>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // 32 <--> 64
    def("applyMapping", registerConverters(&pythonApplyMapping<2, npy_uint32, npy_uint64>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<2, npy_uint64, npy_uint32>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // 8 <--> 64
    def("applyMapping", registerConverters(&pythonApplyMapping<2, npy_uint8, npy_uint64>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<2, npy_uint64, npy_uint8>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // Cases for same input/output dtypes must come last, so they are chosen by default!
    // 8 <--> 8, 32 <--> 32, 64 <--> 64
    def("applyMapping", registerConverters(&pythonApplyMapping<2, npy_uint32, npy_uint32>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<2, npy_uint64, npy_uint64>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<2, npy_uint8, npy_uint8>),   (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // -- 1D

    // 8 <--> 32
    def("applyMapping", registerConverters(&pythonApplyMapping<1, npy_uint8, npy_uint32>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<1, npy_uint32, npy_uint8>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // 32 <--> 64
    def("applyMapping", registerConverters(&pythonApplyMapping<1, npy_uint32, npy_uint64>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<1, npy_uint64, npy_uint32>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // 8 <--> 64
    def("applyMapping", registerConverters(&pythonApplyMapping<1, npy_uint8, npy_uint64>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<1, npy_uint64, npy_uint8>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));

    // Cases for same input/output dtypes must come last, so they are chosen by default!
    // 8 <--> 8, 32 <--> 32, 64 <--> 64
    def("applyMapping", registerConverters(&pythonApplyMapping<1, npy_uint32, npy_uint32>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<1, npy_uint64, npy_uint64>), (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
    def("applyMapping", registerConverters(&pythonApplyMapping<1, npy_uint8, npy_uint8>),   (arg("src"), arg("mapping"), arg("allow_incomplete_mapping")=false, arg("out")=python::object()));
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
