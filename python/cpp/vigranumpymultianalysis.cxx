/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpycmodule_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx> 
#include <vigra/numpy_array_converters.hxx>
#include <vigra/localminmax.hxx>
#include <vigra/labelvolume.hxx>
#include <vigra/multi_distance.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/watersheds3d.hxx>
#include <vigra/seededregiongrowing3d.hxx>

#include <cmath>

namespace python = boost::python;

namespace vigra
{

template < class PixelType >
NumpyAnyArray pythonWatersheds3D(NumpyArray<3, Singleband<PixelType> > volume,
                           int neighborhood = 6,
                           NumpyArray<3, Singleband<Int32> > res=python::object())
{
    vigra_precondition(neighborhood == 6 || neighborhood == 26,
       "watersheds3D(volume, neighborhood): neighborhood must be 6 or 26.");
    res.reshapeIfEmpty(volume.shape(), "Watersheds(): Output array has wrong shape.");  
    switch (neighborhood)
    {
        case 6:
        {
            watersheds3DSix(srcMultiArrayRange(volume), destMultiArray(res));
            break;
        }
        case 26:
        {
            watersheds3DTwentySix(srcMultiArrayRange(volume), destMultiArray(res));
            break;
        }
    }
    return res;
}
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pywatersheds3D, pythonWatersheds3D)









template < class VoxelType >
NumpyAnyArray pythonSeededRegionGrowing3D(NumpyArray<3, Singleband<VoxelType> > volume,
                                    SRGType keepContours,
                                    NumpyArray<3, Singleband<Int32> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "seededRegionGrowing3D(): Output array has wrong shape.");    
    extendedLocalMinima(srcMultiArrayRange(volume), destMultiArray(res));
    int max_region_label =
        labelVolumeWithBackground(srcMultiArrayRange(res),
        destMultiArray(res), NeighborCode3DSix(), 0);
    ArrayOfRegionStatistics< SeedRgDirectValueFunctor< VoxelType > >
        stats(max_region_label);

    seededRegionGrowing3D(srcMultiArrayRange(volume),
        srcMultiArray(res), 
        destMultiArray(res), stats,
        keepContours);
    return res;
}

template < class VoxelType >
NumpyAnyArray pythonSeededRegionGrowingSeeded3D(NumpyArray<3, Singleband<VoxelType> > volume,
                                          NumpyArray<3, Singleband<Int32> > seeds,
                                          SRGType srgType=CompleteGrow,
                                          NumpyArray<3, Singleband<Int32> > res=python::object())
{
    vigra_precondition(volume.shape() == seeds.shape(),
         "seededRegionGrowingSeeded3D(): magnitude and seed volumes must have the same "
         "size.");
    res.reshapeIfEmpty(seeds.shape(), "seededRegionGrowingSeeded3D(): Output array has wrong shape.");    

    FindMinMax< Int32 > minmax;
    inspectMultiArray(srcMultiArrayRange(seeds), minmax);
    // create a statistics functor for region growing
    ArrayOfRegionStatistics< SeedRgDirectValueFunctor< VoxelType > >
        stats(minmax.max);
    seededRegionGrowing3D(srcMultiArrayRange(volume), srcMultiArray(seeds),
        destMultiArray(res), stats, srgType);

    return res;
}

template < class VoxelType >
NumpyAnyArray pythonLabelVolume3D(NumpyArray<3, Singleband<VoxelType> > volume, 
                            int neighborhood=6,
                            NumpyArray<3, Singleband<VoxelType> > res=python::object())
{
    
    res.reshapeIfEmpty(volume.shape(), "labelVolume3D(): Output array has wrong shape.");
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
    return res;
}

template < class VoxelType >
NumpyAnyArray pythonLabelVolumeWithBackground3D(NumpyArray<3, Singleband<VoxelType> > volume, 
                                          int neighborhood=6,
                                          Int32 background_value = 0,
                                          NumpyArray<3, Singleband<VoxelType> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "labelVolumeWithBackground3D(): Output array has wrong shape.");
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
    return res;
}

template < class VoxelType >
NumpyAnyArray pythonDistanceTransform3D(NumpyArray<3, Singleband<VoxelType> > volume, 
                                  bool background,
                                  NumpyArray<3, Singleband<VoxelType> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "distanceTransform3D(): Output array has wrong shape.");
    
    separableMultiDistance(srcMultiArrayRange(volume),
        destMultiArray(res), background);
    return res;
}




namespace { void DUMMY_FUNCTION(int, int, int, int, int) {} }

void defineMultiAnalysisFunctions()
{
    using namespace python;

    multidef("watersheds", pywatersheds3D< UInt8, float >(),
      (arg("volume"), arg("neighborhood")=6, arg("out")=python::object()));

/*
    def("seededRegionGrowing3D",
        registerConverters(&seededRegionGrowing3D<float>),
        (arg("volume"),
        arg("srgType")=CompleteGrow,
        arg("out")=python::object()));
*/
    def("seededRegionGrowingSeeded3D",
        registerConverters(&pythonSeededRegionGrowingSeeded3D<float>),
        (arg("volume"), arg("seedVolume"),
        arg("srgType")=CompleteGrow,
        arg("out")=python::object()));


    def("labelVolume3D",
        registerConverters(&pythonLabelVolume3D<Int32>),
        (arg("volume"), 
        arg("neighborhood3D")=6,
        arg("out")=python::object()));

    def("labelVolumeWithBackground3D",
        registerConverters(&pythonLabelVolumeWithBackground3D<Int32>),
        (arg("volume"), 
         arg("neighborhood")=6, 
         arg("background_value")=0,
         arg("out")=python::object()));

    def("distanceTransform3D",
        registerConverters(&pythonDistanceTransform3D<Int32>),
        (arg("array"), arg("background"),
         arg("out")=python::object()));
}

} // namespace vigra
