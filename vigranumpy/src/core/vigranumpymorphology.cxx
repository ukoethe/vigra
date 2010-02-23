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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpymorphology_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/flatmorphology.hxx>
#include <vigra/multi_morphology.hxx>

#include <cmath>

namespace python = boost::python;

namespace vigra
{

template < class PixelType >
NumpyAnyArray pythonDiscRankOrderFilter(NumpyArray<3, Multiband<PixelType> > image,
                                        int radius, float rank, NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((rank >= 0.0) && (rank <= 1.0),
        "Rank must be in the range 0.0 <= rank <= 1.0");
    vigra_precondition(radius >= 0, "Radius must be >= 0.");

    res.reshapeIfEmpty(image.shape(),"discRankOrderFilter(): Output image has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    { 
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        discRankOrderFilter(srcImageRange(bimage,StandardValueAccessor<UInt8>()), destImage(bres), radius, rank);
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonDiscErosion(NumpyArray<3, Multiband<PixelType> > image,int radius, NumpyArray<3, Multiband<PixelType> > res)
{
    return pythonDiscRankOrderFilter(image,radius,0.0f,res);
}

template < class PixelType >
NumpyAnyArray pythonDiscDilation(NumpyArray<3, Multiband<PixelType> > image,int radius, NumpyArray<3, Multiband<PixelType> > res)
{
    return pythonDiscRankOrderFilter(image,radius,1.0f,res);
}

template < class PixelType >
NumpyAnyArray pythonDiscMedian(NumpyArray<3, Multiband<PixelType> > image,int radius, NumpyArray<3, Multiband<PixelType> > res)
{
    return pythonDiscRankOrderFilter(image,radius,0.5f,res);
}

template < class PixelType >
NumpyAnyArray pythonDiscOpening(NumpyArray<3, Multiband<PixelType> > image, int radius, NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition(radius >= 0, "Radius must be >=0.");

    res.reshapeIfEmpty(image.shape(),"discOpening(): Output image has wrong dimensions");

    MultiArray<2,PixelType> tmp(MultiArrayShape<2>::type(image.shape(0),image.shape(1)));

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        discErosion(srcImageRange(bimage), destImage(tmp), radius);
        discDilation(srcImageRange(tmp), destImage(bres), radius);
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonDiscClosing(NumpyArray<3, Multiband<PixelType> > image, int radius, NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition(radius >= 0, "Radius must be >=0.");

    res.reshapeIfEmpty(image.shape(),"discClosing(): Output image has wrong dimensions");

    MultiArray<2,PixelType> tmp(MultiArrayShape<2>::type(image.shape(0),image.shape(1)));

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        discDilation(srcImageRange(bimage), destImage(tmp), radius);
        discErosion(srcImageRange(tmp), destImage(bres), radius);
    }
    return res;
}

template < int dim >
NumpyAnyArray pythonMultiBinaryErosion(NumpyArray<dim, Multiband<bool> > image,int radius, NumpyArray<dim, Multiband<bool> > res)
{
    res.reshapeIfEmpty(image.shape(),"multiBinaryErosion(): Output image has wrong dimensions");

    for(int k=0;k<image.shape(dim-1);++k)
    {
        MultiArrayView<dim-1, bool, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<dim-1, bool, StridedArrayTag> bres=res.bindOuter(k);
        multiBinaryErosion(srcMultiArrayRange(bimage),destMultiArray(bres), radius);
    }
    return res;
}
template < int dim >
NumpyAnyArray pythonMultiBinaryDilation(NumpyArray<dim, Multiband<bool> > image,int radius, NumpyArray<dim, Multiband<bool> > res)
{
    res.reshapeIfEmpty(image.shape(),"multiBinaryDilation(): Output image has wrong dimensions");

    for(int k=0;k<image.shape(dim-1);++k)
    {
        MultiArrayView<dim-1, bool, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<dim-1, bool, StridedArrayTag> bres=res.bindOuter(k);
        multiBinaryDilation(srcMultiArrayRange(bimage),destMultiArray(bres), radius);
    }
    return res;
}

template <int dim>
NumpyAnyArray pythonMultiBinaryOpening(NumpyArray<dim, Multiband<bool> > image, int radius, NumpyArray<dim, Multiband<bool> > res)
{
    res.reshapeIfEmpty(image.shape(),"multiBinaryOpening(): Output image has wrong dimensions");

    TinyVector<UInt32,dim-1> size;
    for(int ii=0;ii<dim-1;++ii)
        size[ii]=image.shape(ii);
    MultiArray<dim-1,bool> tmp(size);

    for(int k=0;k<image.shape(dim-1);++k)
    {
        MultiArrayView<dim-1, bool, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<dim-1, bool, StridedArrayTag> bres = res.bindOuter(k);
        multiBinaryErosion(srcMultiArrayRange(bimage),destMultiArray(tmp), radius);
        multiBinaryDilation(srcMultiArrayRange(tmp),destMultiArray(bres), radius);
    }
    return res;
}

template <int dim>
NumpyAnyArray pythonMultiBinaryClosing(NumpyArray<dim, Multiband<bool> > image, int radius, NumpyArray<dim, Multiband<bool> > res)
{
    res.reshapeIfEmpty(image.shape(),"multiBinaryOpening(): Output image has wrong dimensions");

    TinyVector<UInt32,dim-1> size;
    for(int ii=0;ii<dim-1;++ii)
        size[ii]=image.shape(ii);
    MultiArray<dim-1,bool> tmp(size);

    for(int k=0;k<image.shape(dim-1);++k)
    {
        MultiArrayView<dim-1, bool, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<dim-1, bool, StridedArrayTag> bres = res.bindOuter(k);
        multiBinaryDilation(srcMultiArrayRange(bimage),destMultiArray(tmp), radius);
        multiBinaryErosion(srcMultiArrayRange(tmp),destMultiArray(bres), radius);
    }
    return res;
}

template < int dim , class PixelType>
NumpyAnyArray pythonMultiGrayscaleErosion(NumpyArray<dim, Multiband<PixelType> > image, double sigma, NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(image.shape(),"multiGrayscaleErosion(): Output image has wrong dimensions");

    for(int k=0;k<image.shape(dim-1);++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres=res.bindOuter(k);
        multiGrayscaleErosion(srcMultiArrayRange(bimage),destMultiArray(bres), sigma);
    }
    return res;
}
template < int dim, class PixelType >
NumpyAnyArray pythonMultiGrayscaleDilation(NumpyArray<dim, Multiband<PixelType> > image, double sigma, NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(image.shape(),"multiGrayscaleDilation(): Output image has wrong dimensions");

    for(int k=0;k<image.shape(dim-1);++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres=res.bindOuter(k);
        multiGrayscaleDilation(srcMultiArrayRange(bimage),destMultiArray(bres), sigma);
    }
    return res;
}

template <int dim, class PixelType>
NumpyAnyArray pythonMultiGrayscaleOpening(NumpyArray<dim, Multiband<PixelType> > image, int radius, NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(image.shape(),"multiGrayscaleOpening(): Output image has wrong dimensions");

    TinyVector<UInt32,dim-1> size;
    for(int ii=0;ii<dim-1;++ii)
    {
        size[ii]=image.shape(ii);
    }
    MultiArray<dim-1,PixelType> tmp(size);

    for(int k=0;k<image.shape(dim-1);++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        multiGrayscaleErosion(srcMultiArrayRange(bimage),destMultiArray(tmp), radius);
        multiGrayscaleDilation(srcMultiArrayRange(tmp),destMultiArray(bres), radius);
    }
    return res;
}

template <int dim, class PixelType>
NumpyAnyArray pythonMultiGrayscaleClosing(NumpyArray<dim, Multiband<PixelType> > image, int radius, NumpyArray<dim, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(image.shape(),"multiGrayscaleClosing(): Output image has wrong dimensions");

    TinyVector<UInt32,dim-1> size;
    for(int ii=0;ii<dim-1;++ii)
        size[ii]=image.shape(ii);
    MultiArray<dim-1,PixelType> tmp(size);

    for(int k=0;k<image.shape(dim-1);++k)
    {
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        multiGrayscaleDilation(srcMultiArrayRange(bimage),destMultiArray(tmp), radius);
        multiGrayscaleErosion(srcMultiArrayRange(tmp),destMultiArray(bres), radius);
    }
    return res;
}
template < class PixelType >
NumpyAnyArray pythonDiscRankOrderFilterWithMask(NumpyArray<3, Multiband<PixelType> > image,
                                                NumpyArray<3, Multiband<PixelType> > mask,
                                                int radius, float rank,
                                                NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((rank >= 0.0) && (rank <= 1.0),
        "Rank must be in the range 0.0 <= rank <= 1.0");
    vigra_precondition(radius >= 0, "Radius must be >= 0.");

    res.reshapeIfEmpty(image.shape(),"discRankOrderFilterWithMask(): Output image has wrong dimensions");
    vigra_precondition(mask.shape(2)==1 || mask.shape(2)==image.shape(2),
                       "discRankOrderFilterWithMask(): mask image must either have 1 channel or as many as the input image");
    vigra_precondition(mask.shape(0)==image.shape(0) && mask.shape(1)==image.shape(1),
                       "discRankOrderFilterWithMaks(): mask dimensions must be same as image dimensions");

    for(int k=0;k<image.shape(2);++k)
    { 
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bmask = mask.bindOuter(mask.shape(2)==1?0:k);
        discRankOrderFilterWithMask(srcImageRange(bimage,StandardValueAccessor<UInt8>()), srcImage(bmask),
                                    destImage(bres), radius, rank);
    }

    return res;
}


void defineMorphology()
{
    using namespace python;

    def("discRankOrderFilter",
        registerConverters(&pythonDiscRankOrderFilter<UInt8>),
        (arg("image"), arg("radius"), arg("rank"), arg("out")=object()),
        "Apply rank order filter with disc structuring function to the image.\n\n"
        "The pixel values of the source image  must be in the range 0...255. Radius must be >= 0.\n"
        "Rank must be in the range 0.0 <= rank <= 1.0. The filter acts as a minimum filter if rank = 0.0, as a median\n"
        "if rank = 0.5, and as a maximum filter if rank = 1.0.\n"
        "\n" 
        "For details see vigra documentation discRankOrderFilter_ .\n"
        "This function also works for multiband images, it is then executed on every band.\n"
       );
    def("discRankOrderFilter",
        registerConverters(&pythonDiscRankOrderFilter<bool>),
        (arg("image"), arg("radius"), arg("rank"), arg("out")=object()));

    def("discErosion",
        registerConverters(&pythonDiscErosion<bool>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply erosion (minimum) filter with disc of given radius to image.\n\n"
        "This is an abbreviation for the rank order filter with rank = 0.0. See discErosion_ for more information.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
        );
    def("discErosion",
        registerConverters(&pythonDiscErosion<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()));

    def("discDilation",
        registerConverters(&pythonDiscDilation<bool>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply dilation (maximum) filter with disc of given radius to image.\n\n"
        "This is an abbreviation for the rank order filter with rank = 1.0. See discDilation_ for more information.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
       );
    def("discDilation",
        registerConverters(&pythonDiscDilation<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()));

    def("discMedian",
        registerConverters(&pythonDiscMedian<bool>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply median filter with disc of given radius to image.\n\n"
        "This is an abbreviation for the rank order filter with rank = 0.5. See discMedian_ for more information.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
        );
    def("discMedian",
        registerConverters(&pythonDiscMedian<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()));

    def("discOpening",
        registerConverters(&pythonDiscOpening<bool>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply a opening filter with disc of given radius to image.\n\n"
        "This is an abbreviation for applying an erosion and a dilation filter in sequence. See discRankOrderFilter_ for more information\n"
        "This function also works for multiband images, it is then executed on every band.\n"
       );
    def("discOpening",
        registerConverters(&pythonDiscOpening<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()));

    def("discClosing",
        registerConverters(&pythonDiscClosing<bool>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply a closing filter with disc of given radius to image.\n\n"
        "This is an abbreviation for applying a dilation and an erosion  filter in sequence.\n"
        "See discRankOrderFilter_ for more information\n"
        "This function also works for multiband images, it is then executed on every band.\n"
       );
    def("discClosing",
        registerConverters(&pythonDiscClosing<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()));

    def("multiBinaryErosion",
        registerConverters(&pythonMultiBinaryErosion<4>),
        (arg("image"), arg("radius"), arg("out")=object()),
       "Binary erosion on multi-dimensional arrays.\n"
       "\n"
       "This function applies a flat circular erosion operator with a given radius. The operation is isotropic.\n"
       "The input is a binary multi-dimensional array where non-zero pixels represent foreground and zero pixels represent background.\n"
       "\n"
       "For details see vigra documentation multiBinaryErosion_.\n"
       "This function also works for multiband arrays, it is then executed on every band.\n"
        );
    def("multiBinaryDilation",
        registerConverters(&pythonMultiBinaryDilation<4>),
        (arg("image"), arg("radius"), arg("out")=object()),
       "Binary dilation on multi-dimensional arrays.\n"
       "\n"
       "This function applies a flat circular dilation operator with a given radius. The operation is isotropic.\n"
       "The input is a binary multi-dimensional array where non-zero pixels represent foreground and zero pixels represent background.\n"
       "\n"
       "For details see vigra documentation multiBinaryDilation_.\n"
       "This function also works for multiband arrays, it is then executed on every band.\n");
    
    def("multiBinaryOpening",
        registerConverters(&pythonMultiBinaryOpening<4>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Binary opening on multi-dimensional arrays.\n"
        "\n"
        "This function applies a flat circular opening operator (sequential erosion and dilation) with a given radius.\n"
        "The operation is isotropic.\n"
        "The input is a binary multi-dimensional array where non-zero pixels represent foreground and zero pixels represent background.\n"
        "\n"
        "For details see vigra documentation (multiBinaryDilation_ and multiBinaryErosion_).\n"
        "This function also works for multiband arrays, it is then executed on every band.\n");
    def("multiBinaryClosing",
        registerConverters(&pythonMultiBinaryClosing<4>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Binary closing on multi-dimensional arrays.\n"
        "\n"
        "This function applies a flat circular opening operator (sequential dilation and erosion) with a given radius.\n"
        "The operation is isotropic.\n"
        "The input is a binary multi-dimensional array where non-zero pixels represent foreground and zero pixels represent background.\n"
        "\n"
        "For details see vigra documentation (multiBinaryDilation_ and multiBinaryErosion_).\n"
        "This function also works for multiband arrays, it is then executed on every band.\n");
    
    def("multiGrayscaleErosion",
        registerConverters(&pythonMultiGrayscaleErosion<4,UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Parabolic grayscale erosion on multi-dimensional arrays.\n"
        "\n"
        "This function applies a parabolic erosion operator with a given spread (sigma) on a grayscale array.\n"
        "The operation is isotropic. The input is a grayscale multi-dimensional array.\n"
        "\n"
        "For details see vigra documentation multiGrayscaleErosion_.\n"
        "This function also works for multiband arrays, it is then executed on every band.\n");
    def("multiGrayscaleDilation",
        registerConverters(&pythonMultiGrayscaleDilation<4,UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Parabolic grayscale dilation on multi-dimensional arrays.\n"
        "\n"
        "This function applies a parabolic dilation operator with a given spread (sigma) on a grayscale array.\n"
        "The operation is isotropic. The input is a grayscale multi-dimensional array.\n"
        "\n"
        "For details see vigra documentation multiGrayscaleDilation_.\n"
        "This function also works for multiband arrays, it is then executed on every band.\n");
    def("multiGrayscaleErosion",
        registerConverters(&pythonMultiGrayscaleErosion<4,float>),
        (arg("image"), arg("sigma"), arg("out")=object()));
    def("multiGrayscaleDilation",
        registerConverters(&pythonMultiGrayscaleDilation<4,float>),
        (arg("image"), arg("sigma"), arg("out")=object()));

    def("multiGrayscaleOpening",
        registerConverters(&pythonMultiGrayscaleOpening<4,UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Parabolic grayscale opening on multi-dimensional arrays.\n"
        "\n"
        "This function applies a parabolic opening (sequencial erosion and dilation) operator with a given spread (sigma)"
        "on a grayscale array.\n"
        "The operation is isotropic. The input is a grayscale multi-dimensional array.\n"
        "\n"
        "For details see vigra documentation multiGrayscaleDilation_ and multiGrayscaleErosion_.\n"
        "This function also works for multiband arrays, it is then executed on every band.\n");
    def("multiGrayscaleClosing",
        registerConverters(&pythonMultiGrayscaleClosing<4,UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Parabolic grayscale closing on multi-dimensional arrays.\n"
        "\n"
        "This function applies a parabolic closing (sequencial dilation and erosion) operator with a given spread (sigma)"
        "on a grayscale array.\n"
        "The operation is isotropic. The input is a grayscale multi-dimensional array.\n"
        "\n"
        "For details see vigra documentation multiGrayscaleDilation_ and multiGrayscaleErosion_.\n"
        "This function also works for multiband arrays, it is then executed on every band.\n");
    def("multiGrayscaleOpening",
        registerConverters(&pythonMultiGrayscaleOpening<4,float>),
        (arg("image"), arg("sigma"), arg("out")=object()));
    def("multiGrayscaleClosing",
        registerConverters(&pythonMultiGrayscaleClosing<4,float>),
        (arg("image"), arg("sigma"), arg("out")=object()));

    def("discRankOrderFilterWithMask",
        registerConverters(&pythonDiscRankOrderFilterWithMask<float>),
        (arg("image"), arg("mask"), arg("radius"), arg("rank"), arg("out")=object()),
        "Apply rank order filter with disc structuring function to the image using a mask.\n"
        "\n"
        "The pixel values of the source image must be in the range 0...255. Radius must be >= 0.\n"
        "Rank must be in the range 0.0 <= rank <= 1.0. The filter acts as a minimum filter if rank = 0.0,"
        "as a median if rank = 0.5, and as a maximum filter if rank = 1.0.\n"
        "\n"
        "The mask is only applied to th input image, i.e. the function generates an output wherever the current disc contains\n"
        "at least one pixel with mask value 'true'. Source pixels with mask value 'false' are ignored during the calculation of\n"
        "the rank order.\n"
        "If the mask has only one band, it is used for every image band. If the mask has the same number of bands, as the image\n"
        "the bands are used for the according images.\n"
        "For details see vigra documentation discRankOrderFilterWithMask_.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
        );
    
    def("discRankOrderFilterWithMask",
        registerConverters(&pythonDiscRankOrderFilterWithMask<UInt8>),
        (arg("image"), arg("mask"), arg("radius"), arg("rank"), arg("out")=object()));

}

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(morphology)
{
    import_vigranumpy();
    defineMorphology();
}
