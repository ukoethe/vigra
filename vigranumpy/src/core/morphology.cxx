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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyfilters_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/flatmorphology.hxx>
#include <vigra/multi_morphology.hxx>
#include <vigra/distancetransform.hxx>
#include <vigra/multi_distance.hxx>

namespace python = boost::python;

namespace vigra
{

template < class PixelType >
NumpyAnyArray 
pythonDiscRankOrderFilter(NumpyArray<3, Multiband<PixelType> > image,
                          int radius, float rank, 
                          NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((rank >= 0.0) && (rank <= 1.0),
        "Rank must be in the range 0.0 <= rank <= 1.0");
    vigra_precondition(radius >= 0, "Radius must be >= 0.");

    res.reshapeIfEmpty(image.taggedShape(),
            "discRankOrderFilter(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    for(int k=0; k<image.shape(2); ++k)
    { 
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        discRankOrderFilter(srcImageRange(bimage,StandardValueAccessor<UInt8>()), 
                            destImage(bres), radius, rank);
    }
    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonDiscRankOrderFilterWithMask(NumpyArray<3, Multiband<PixelType> > image,
                                  NumpyArray<3, Multiband<PixelType> > mask,
                                  int radius, float rank,
                                  NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((rank >= 0.0) && (rank <= 1.0),
        "Rank must be in the range 0.0 <= rank <= 1.0");
    vigra_precondition(radius >= 0, "Radius must be >= 0.");
    vigra_precondition(mask.shape(2)==1 || mask.shape(2)==image.shape(2),
               "discRankOrderFilterWithMask(): mask image must either have 1 channel or as many as the input image");
    vigra_precondition(mask.shape(0)==image.shape(0) && mask.shape(1)==image.shape(1),
               "discRankOrderFilterWithMaks(): mask dimensions must be same as image dimensions");

    res.reshapeIfEmpty(image.taggedShape(),
            "discRankOrderFilterWithMask(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    for(int k=0; k<image.shape(2); ++k)
    { 
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bmask = mask.bindOuter(mask.shape(2)==1?0:k);
        discRankOrderFilterWithMask(srcImageRange(bimage,StandardValueAccessor<UInt8>()), 
                                    srcImage(bmask),
                                    destImage(bres), radius, rank);
    }

    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonDiscErosion(NumpyArray<3, Multiband<PixelType> > image,
                  int radius, 
                  NumpyArray<3, Multiband<PixelType> > res)
{
    return pythonDiscRankOrderFilter(image, radius, 0.0f, res);
}

template < class PixelType >
NumpyAnyArray 
pythonDiscDilation(NumpyArray<3, Multiband<PixelType> > image,
                   int radius, 
                   NumpyArray<3, Multiband<PixelType> > res)
{
    return pythonDiscRankOrderFilter(image, radius, 1.0f, res);
}

template < class PixelType >
NumpyAnyArray 
pythonDiscMedian(NumpyArray<3, Multiband<PixelType> > image,
                 int radius, 
                 NumpyArray<3, Multiband<PixelType> > res)
{
    return pythonDiscRankOrderFilter(image, radius, 0.5f, res);
}

template < class PixelType >
NumpyAnyArray 
pythonDiscOpening(NumpyArray<3, Multiband<PixelType> > image, 
                  int radius, 
                  NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition(radius >= 0, "Radius must be >=0.");

    res.reshapeIfEmpty(image.taggedShape(),
            "discOpening(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    MultiArray<2,PixelType> tmp(MultiArrayShape<2>::type(image.shape(0), image.shape(1)));

    for(int k=0; k<image.shape(2); ++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        discErosion(srcImageRange(bimage), destImage(tmp), radius);
        discDilation(srcImageRange(tmp), destImage(bres), radius);
    }
    return res;
}

template < class PixelType >
NumpyAnyArray 
pythonDiscClosing(NumpyArray<3, Multiband<PixelType> > image, 
                  int radius, 
                  NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition(radius >= 0, "Radius must be >=0.");

    res.reshapeIfEmpty(image.taggedShape(),
            "discClosing(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    MultiArray<2,PixelType> tmp(MultiArrayShape<2>::type(image.shape(0), image.shape(1)));

    for(int k=0; k<image.shape(2); ++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        discDilation(srcImageRange(bimage), destImage(tmp), radius);
        discErosion(srcImageRange(tmp), destImage(bres), radius);
    }
    return res;
}

template < int dim, class PixelType >
NumpyAnyArray 
pythonMultiBinaryErosion(NumpyArray<dim, Multiband<PixelType> > array,
                         double radius, 
                         NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(array.taggedShape(),
            "multiBinaryErosion(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    for(int k=0; k<array.shape(dim-1); ++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> barray = array.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        multiBinaryErosion(srcMultiArrayRange(barray), destMultiArray(bres), radius);
    }
    return res;
}

template < int dim, class PixelType >
NumpyAnyArray 
pythonMultiBinaryDilation(NumpyArray<dim, Multiband<PixelType> > array,
                          double radius, 
                          NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(array.taggedShape(),
            "multiBinaryDilation(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    for(int k=0; k<array.shape(dim-1); ++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> barray = array.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        multiBinaryDilation(srcMultiArrayRange(barray), destMultiArray(bres), radius);
    }
    return res;
}

template <int dim, class PixelType >
NumpyAnyArray 
pythonMultiBinaryOpening(NumpyArray<dim, Multiband<PixelType> > array, 
                         double radius, 
                         NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(array.taggedShape(),
            "multiBinaryOpening(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    MultiArray<dim-1,PixelType> tmp(typename MultiArrayShape<dim-1>::type(array.shape().begin()));

    for(int k=0; k<array.shape(dim-1); ++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> barray = array.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        multiBinaryErosion(srcMultiArrayRange(barray), destMultiArray(tmp), radius);
        multiBinaryDilation(srcMultiArrayRange(tmp), destMultiArray(bres), radius);
    }
    return res;
}

template <int dim, class PixelType >
NumpyAnyArray 
pythonMultiBinaryClosing(NumpyArray<dim, Multiband<PixelType> > array, 
                         double radius, 
                         NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(array.taggedShape(),
            "multiBinaryOpening(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    MultiArray<dim-1,PixelType> tmp(typename MultiArrayShape<dim-1>::type(array.shape().begin()));

    for(int k=0; k<array.shape(dim-1); ++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> barray = array.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        multiBinaryDilation(srcMultiArrayRange(barray), destMultiArray(tmp), radius);
        multiBinaryErosion(srcMultiArrayRange(tmp), destMultiArray(bres), radius);
    }
    return res;
}

template < int dim , class PixelType>
NumpyAnyArray 
pythonMultiGrayscaleErosion(NumpyArray<dim, Multiband<PixelType> > array, 
                            double sigma, 
                            NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(array.taggedShape(),
            "multiGrayscaleErosion(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    for(int k=0; k<array.shape(dim-1); ++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> barray = array.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        multiGrayscaleErosion(srcMultiArrayRange(barray), destMultiArray(bres), sigma);
    }
    return res;
}
template < int dim, class PixelType >
NumpyAnyArray 
pythonMultiGrayscaleDilation(NumpyArray<dim, Multiband<PixelType> > array, 
                             double sigma, 
                             NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(array.taggedShape(),
            "multiGrayscaleDilation(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    for(int k=0; k<array.shape(dim-1); ++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> barray = array.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        multiGrayscaleDilation(srcMultiArrayRange(barray), destMultiArray(bres), sigma);
    }
    return res;
}

template <int dim, class PixelType>
NumpyAnyArray 
pythonMultiGrayscaleOpening(NumpyArray<dim, Multiband<PixelType> > array, 
                            double sigma, 
                            NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(array.taggedShape(),
            "multiGrayscaleOpening(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    MultiArray<dim-1,PixelType> tmp(typename MultiArrayShape<dim-1>::type(array.shape().begin()));

    for(int k=0; k<array.shape(dim-1); ++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> barray = array.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        multiGrayscaleErosion(srcMultiArrayRange(barray), destMultiArray(tmp), sigma);
        multiGrayscaleDilation(srcMultiArrayRange(tmp), destMultiArray(bres), sigma);
    }
    return res;
}

template <int dim, class PixelType>
NumpyAnyArray 
pythonMultiGrayscaleClosing(NumpyArray<dim, Multiband<PixelType> > array, 
                            double sigma, 
                            NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(array.taggedShape(),
            "multiGrayscaleClosing(): Output image has wrong dimensions");

    PyAllowThreads _pythread;
    MultiArray<dim-1,PixelType> tmp(typename MultiArrayShape<dim-1>::type(array.shape().begin()));

    for(int k=0; k<array.shape(dim-1); ++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> barray = array.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        multiGrayscaleDilation(srcMultiArrayRange(barray), destMultiArray(tmp), sigma);
        multiGrayscaleErosion(srcMultiArrayRange(tmp), destMultiArray(bres), sigma);
    }
    return res;
}

namespace detail {

template <class PixelType>
struct IsBackgroundAccessor
{
    typedef bool value_type;
    
    template <class Iterator>
    value_type operator()(Iterator const & i) const
    {
        return *i == NumericTraits<PixelType>::zero();
    }
};

} // namespace detail

template < class PixelType, typename DestPixelType >
NumpyAnyArray 
pythonDistanceTransform2D(NumpyArray<2, Singleband<PixelType> > image,
                          bool background, 
                          int norm,
                          NumpyArray<2, Singleband<DestPixelType> > res = python::object())
{
    res.reshapeIfEmpty(image.taggedShape(), 
            "distanceTransform2D(): Output array has wrong shape.");
    
    PyAllowThreads _pythread;
    if(background)
    {
        distanceTransform(srcImageRange(image), destImage(res), 
                          NumericTraits<PixelType>::zero(), norm);
    }
    else
    {
        distanceTransform(srcImageRange(image, detail::IsBackgroundAccessor<PixelType>()), 
                          destImage(res), false, norm);
    }
    return res;
}

template < class VoxelType >
NumpyAnyArray 
pythonDistanceTransform3D(NumpyArray<3, Singleband<VoxelType> > volume, 
                          bool background, 
                          NumpyArray<3, Singleband<VoxelType> > res=python::object())
{
    res.reshapeIfEmpty(volume.taggedShape(), 
            "distanceTransform3D(): Output array has wrong shape.");
    
    PyAllowThreads _pythread;
    separableMultiDistance(srcMultiArrayRange(volume), destMultiArray(res), background);
    return res;
}

void defineMorphology()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);

    def("discRankOrderFilter",
        registerConverters(&pythonDiscRankOrderFilter<float>),
        (arg("image"), arg("radius"), arg("rank"), arg("out")=object()),
        "Apply rank order filter with disc structuring function to a float image.\n\n"
        "The pixel values of the source image  must be in the range 0...255. Radius must be >= 0. "
        "Rank must be in the range 0.0 <= rank <= 1.0. The filter acts as a minimum filter if rank = 0.0, as a median "
        "if rank = 0.5, and as a maximum filter if rank = 1.0. "
        "This function also works for multiband images, it is then executed on every band.\n"
        "\n" 
        "For details see discRankOrderFilter_ in the C++ documentation.\n"
       );

    def("discRankOrderFilter",
        registerConverters(&pythonDiscRankOrderFilter<UInt8>),
        (arg("image"), arg("radius"), arg("rank"), arg("out")=object()),
        "Likewise for a uint8 image.\n");

    def("discRankOrderFilterWithMask",
        registerConverters(&pythonDiscRankOrderFilterWithMask<float>),
        (arg("image"), arg("mask"), arg("radius"), arg("rank"), arg("out")=object()),
        "Apply rank order filter with disc structuring function to a float image using a mask.\n"
        "\n"
        "The pixel values of the source image must be in the range 0...255. Radius must be >= 0."
        "Rank must be in the range 0.0 <= rank <= 1.0. The filter acts as a minimum filter if rank = 0.0,"
        "as a median if rank = 0.5, and as a maximum filter if rank = 1.0.\n"
        "\n"
        "The mask is only applied to the input image, i.e. the function generates an output "
        "wherever the current disc contains at least one pixel with non-zero mask value. "
        "Source pixels with mask value zero are ignored during the calculation of "
        "the rank order.\n\n"
        "This function also works for multiband images, it is then executed on every band. "
        "If the mask has only one band, it is used for every image band. If the mask has "
        "the same number of bands, as the image the bands are used for the corresponding image bands.\n\n"
        "For details see discRankOrderFilterWithMask_ in the C++ documentation.\n"
        );
    
    def("discRankOrderFilterWithMask",
        registerConverters(&pythonDiscRankOrderFilterWithMask<UInt8>),
        (arg("image"), arg("mask"), arg("radius"), arg("rank"), arg("out")=object()),
        "Likewise for a uint8 image.\n");

    def("discErosion",
        registerConverters(&pythonDiscErosion<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply erosion (minimum) filter with disc of given radius to image.\n\n"
        "This is an abbreviation for the rank order filter with rank = 0.0. "
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "See discErosion_ in the C++ documentation for more information.\n"
        );

    def("discDilation",
        registerConverters(&pythonDiscDilation<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply dilation (maximum) filter with disc of given radius to image.\n\n"
        "This is an abbreviation for the rank order filter with rank = 1.0. "
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "See discDilation_ in the C++ documentation for more information.\n"
       );

    def("discMedian",
        registerConverters(&pythonDiscMedian<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply median filter with disc of given radius to image.\n\n"
        "This is an abbreviation for the rank order filter with rank = 0.5. "
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "See discMedian_ in the C++ documentation for more information.\n"
        );

    def("discOpening",
        registerConverters(&pythonDiscOpening<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply a opening filter with disc of given radius to image.\n\n"
        "This is an abbreviation for applying an erosion and a dilation filter in sequence. "
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "See discRankOrderFilter_ in the C++ documentation for more information.\n"
       );

    def("discClosing",
        registerConverters(&pythonDiscClosing<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply a closing filter with disc of given radius to image.\n\n"
        "This is an abbreviation for applying a dilation and an erosion  filter in sequence. "
        "This function also works for multiband images, it is then executed on every band.\n\n"
        "See discRankOrderFilter_ in the C++ documentation for more information.\n"
       );

    def("multiBinaryErosion",
        registerConverters(&pythonMultiBinaryErosion<4, UInt8>),
        (arg("volume"), arg("radius"), arg("out")=object()),
       "Binary erosion on a 3D scalar or multiband uint8 array.\n"
       "\n"
       "This function applies a flat circular erosion operator with a given radius. "
       "The operation is isotropic. The input is a uint8 or boolean multi-dimensional array "
       "where non-zero pixels represent foreground and zero pixels represent background. "
       "This function also works for multiband arrays, it is then executed on every band.\n"
       "\n"
       "For details see multiBinaryErosion_ in the C++ documentation.\n"
        );
        
    def("multiBinaryErosion",
        registerConverters(&pythonMultiBinaryErosion<4, bool>),
        (arg("volume"), arg("radius"), arg("out")=object()),
        "Likewise for a bool array.\n");

    def("multiBinaryDilation",
        registerConverters(&pythonMultiBinaryDilation<4, UInt8>),
        (arg("volume"), arg("radius"), arg("out")=object()),
       "Binary dilation on a 3D scalar or multiband uint8 array.\n"
       "\n"
       "This function applies a flat circular dilation operator with a given radius. "
       "The operation is isotropic. The input is a uint8 or boolean multi-dimensional array "
       "where non-zero pixels represent foreground and zero pixels represent background. "
       "This function also works for multiband arrays, it is then executed on every band.\n"
       "\n"
       "For details see multiBinaryDilation_ in the C++ documentation.\n"
       );
       
    def("multiBinaryDilation",
        registerConverters(&pythonMultiBinaryDilation<4, bool>),
        (arg("volume"), arg("radius"), arg("out")=object()),
        "Likewise for bool arrays.\n");
    
    def("multiBinaryOpening",
        registerConverters(&pythonMultiBinaryOpening<4, UInt8>),
        (arg("volume"), arg("radius"), arg("out")=object()),
        "Binary opening on a 3D scalar or multiband uint8 array.\n"
        "\n"
        "This function applies a flat circular opening operator (sequential erosion "
        "and dilation) with a given radius. The operation is isotropic. "
        "The input is a uint8 or boolean multi-dimensional array where non-zero pixels represent "
        "foreground and zero pixels represent background. "
        "This function also works for multiband arrays, it is then executed on every band.\n"
        "\n"
        "For details see vigra C++ documentation (multiBinaryDilation_ and multiBinaryErosion_).\n"
        );
        
    def("multiBinaryOpening",
        registerConverters(&pythonMultiBinaryOpening<4, bool>),
        (arg("volume"), arg("radius"), arg("out")=object()),
        "Likewise for a bool array.\n");
        
    def("multiBinaryClosing",
        registerConverters(&pythonMultiBinaryClosing<4, UInt8>),
        (arg("volume"), arg("radius"), arg("out")=object()),
        "Binary closing on a 3D scalar or multiband uint8 array.\n"
        "\n"
        "This function applies a flat circular opening operator (sequential dilation "
        "and erosion) with a given radius. The operation is isotropic. "
        "The input is a uint8 or boolean multi-dimensional array where non-zero pixels represent "
        "foreground and zero pixels represent background. "
        "This function also works for multiband arrays, it is then executed on every band.\n"
        "\n"
        "For details see vigra C++ documentation (multiBinaryDilation_ and multiBinaryErosion_).\n"
        );
        
    def("multiBinaryClosing",
        registerConverters(&pythonMultiBinaryClosing<4, bool>),
        (arg("volume"), arg("radius"), arg("out")=object()),
        "Likewise for a bool array.\n");
    
    def("multiGrayscaleErosion",
        registerConverters(&pythonMultiGrayscaleErosion<4,UInt8>),
        (arg("volume"), arg("sigma"), arg("out")=object()),
        "Parabolic grayscale erosion on a 3D scalar or multiband uint8 array.\n"
        "\n"
        "This function applies a parabolic erosion operator with a given spread 'sigma' on a grayscale array. "
        "The operation is isotropic. The input is a grayscale multi-dimensional array. "
        "This function also works for multiband arrays, it is then executed on every band.\n"
        "\n"
        "For details see multiGrayscaleErosion_ in the C++ documentation.\n"
        );
                
    def("multiGrayscaleErosion",
        registerConverters(&pythonMultiGrayscaleErosion<4,float>),
        (arg("volume"), arg("sigma"), arg("out")=object()),
        "Likewise for a 3D float array.\n");

    def("multiGrayscaleErosion",
        registerConverters(&pythonMultiGrayscaleErosion<3,UInt8>),
        (arg("image"), arg("sigma"), arg("out")=object()),
        "Likewise for a 2D uint8 array.\n");
    
    def("multiGrayscaleErosion",
        registerConverters(&pythonMultiGrayscaleErosion<3,float>),
        (arg("image"), arg("sigma"), arg("out")=object()),
        "Likewise for a 2D float array.\n");

    def("multiGrayscaleDilation",
        registerConverters(&pythonMultiGrayscaleDilation<4,UInt8>),
        (arg("volume"), arg("sigma"), arg("out")=object()),
        "Parabolic grayscale dilation on multi-dimensional arrays.\n"
        "\n"
        "This function applies a parabolic dilation operator with a given spread 'sigma' on a grayscale array. "
        "The operation is isotropic. The input is a grayscale multi-dimensional array. "
        "This function also works for multiband arrays, it is then executed on every band.\n"
        "\n"
        "For details see multiGrayscaleDilation_ in the C++ documentation.\n"
        );
        
    def("multiGrayscaleDilation",
        registerConverters(&pythonMultiGrayscaleDilation<4,float>),
        (arg("volume"), arg("sigma"), arg("out")=object()),
        "Likewise for a 3D float array.\n");

    def("multiGrayscaleDilation",
        registerConverters(&pythonMultiGrayscaleDilation<3,UInt8>),
        (arg("image"), arg("sigma"), arg("out")=object()),
        "Likewise for a 2D uint8 array.\n");

    def("multiGrayscaleDilation",
        registerConverters(&pythonMultiGrayscaleDilation<3,float>),
        (arg("image"), arg("sigma"), arg("out")=object()));

    def("multiGrayscaleOpening",
        registerConverters(&pythonMultiGrayscaleOpening<4,UInt8>),
        (arg("volume"), arg("sigma"), arg("out")=object()),
        "Parabolic grayscale opening on multi-dimensional arrays.\n"
        "\n"
        "This function applies a parabolic opening (sequential erosion and dilation) "
        "operator with a given spread 'sigma' on a grayscale array. "
        "The operation is isotropic. The input is a grayscale multi-dimensional array. "
        "This function also works for multiband arrays, it is then executed on every band.\n"
        "\n"
        "For details see multiGrayscaleDilation_ and multiGrayscaleErosion_ in the C++ documentation.\n"
        );

    def("multiGrayscaleOpening",
        registerConverters(&pythonMultiGrayscaleOpening<4,float>),
        (arg("volume"), arg("sigma"), arg("out")=object()),
        "Likewise for a 3D float array.\n");

    def("multiGrayscaleOpening",
        registerConverters(&pythonMultiGrayscaleOpening<3,UInt8>),
        (arg("image"), arg("sigma"), arg("out")=object()),
        "Likewise for a 2D uint8 array.\n");

    def("multiGrayscaleOpening",
        registerConverters(&pythonMultiGrayscaleOpening<3,float>),
        (arg("image"), arg("sigma"), arg("out")=object()));

    def("multiGrayscaleClosing",
        registerConverters(&pythonMultiGrayscaleClosing<4,UInt8>),
        (arg("volume"), arg("sigma"), arg("out")=object()),
        "Parabolic grayscale closing on multi-dimensional arrays.\n"
        "\n"
        "This function applies a parabolic closing (sequential dilation and erosion) "
        "operator with a given spread 'sigma' on a grayscale array. "
        "The operation is isotropic. The input is a grayscale multi-dimensional array. "
        "This function also works for multiband arrays, it is then executed on every band.\n"
        "\n"
        "For details see multiGrayscaleDilation_ and multiGrayscaleErosion_ in the C++ documentation.\n"
        );
    
    def("multiGrayscaleClosing",
        registerConverters(&pythonMultiGrayscaleClosing<4,float>),
        (arg("volume"), arg("sigma"), arg("out")=object()),
        "Likewise for a 3D float array.\n");

    def("multiGrayscaleClosing",
        registerConverters(&pythonMultiGrayscaleClosing<3,UInt8>),
        (arg("image"), arg("sigma"), arg("out")=object()),
        "Likewise for a 2D uint8 array.\n");

    def("multiGrayscaleClosing",
        registerConverters(&pythonMultiGrayscaleClosing<3,float>),
        (arg("image"), arg("sigma"), arg("out")=object()));

    def("distanceTransform2D",
        registerConverters(&pythonDistanceTransform2D<float, float>),
        (arg("image"), 
         arg("background")=true, 
         arg("norm")=2,
         arg("out")=python::object()),
        "Compute the distance transform of a 2D scalar float image.\n"
        "All pixels with a value of 0.0 are considered to be background pixels,\n"
        "while all pixels with a nonzero value are considered to be foreground pixels.\n"
        "The parameter 'background' is a Boolean scalar that specifies whether to\n"
        "compute the distance of all background pixels to the nearest foreground pixels\n"
        "(if it is 'True', default) or vice versa (if it is 'False').\n"
        "Hence in the destination image, for background==True all background pixels\n"
        "will be assigned their distance value, while all foreground pixels will be assigned 0.\n"
        "For background==False, it is exactly the other way around.\n\n"
        "The 'norm' parameter gives the distance norm to use\n"
        "(0: infinity norm, 1: L1 norm, 2: Euclidean norm).\n\n"
        "For details see distanceTransform_ in the vigra C++ documentation.\n");

        def("distanceTransform2D",
        registerConverters(&pythonDistanceTransform2D<UInt8,float>),
        (arg("image"), 
         arg("background")=true, 
         arg("norm")=2,
         arg("out")=python::object()),
        "Likewise for a 2D uint8 input array.\n");

    def("distanceTransform3D",
        registerConverters(&pythonDistanceTransform3D<float>),
        (arg("array"), 
         arg("background") = true, 
         arg("out")=python::object()),
        "Compute the Euclidean distance transform of a 3D scalar float volume.\n"
        "All voxels with a value of 0.0 are considered to be background voxels,\n"
        "while all voxels with a nonzero value are considered to be foreground voxels.\n"
        "The parameter 'background' is a Boolean scalar that specifies whether to\n"
        "compute the distance of all background voxels to the nearest foreground voxel\n"
        "(if it is 'True', default) or vice versa (if it is 'False').\n"
        "Hence in the destination volume, for background==True all background voxels\n"
        "will be assigned their distance value, while all foreground voxels will be assigned 0.\n"
        "For background==False, it is exactly the other way around.\n"
        "\n"
        "For more details see separableMultiDistance_ in the vigra C++ documentation.\n");
}

} // namespace vigra
