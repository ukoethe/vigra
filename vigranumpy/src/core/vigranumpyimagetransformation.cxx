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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpycore_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/affinegeometry.hxx>
#include <vigra/basicgeometry.hxx>
#include <vigra/resizeimage.hxx>
#include <vigra/splines.hxx>
#include <vigra/flatmorphology.hxx>
#include <vigra/noise_normalization.hxx>
#include <vigra/multi_morphology.hxx>
#include <vigra/multi_resize.hxx>
#include <vigra/distancetransform.hxx>
#include <vigra/symmetry.hxx>

#include <cmath>

#include "tensors.hxx"
namespace python = boost::python;

namespace vigra
{
    
template < class PixelType >
NumpyAnyArray pythonResampleImage(NumpyArray<3, Multiband<PixelType> > image, double factor, NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((image.shape(0) > 1) && (image.shape(1) > 1),
        "The input image must have a size of at least 2x2.");
    int height,width;
    if(factor < 1.0)
    {
        width = (int) std::ceil(factor * image.shape(0));
        height= (int) std::ceil(factor * image.shape(1));
    }
    else
    {
        width = (int) std::ceil(factor * image.shape(0));
        height= (int) std::ceil(factor * image.shape(1));
    }
    res.reshapeIfEmpty(MultiArrayShape<3>::type(width,height,image.shape(2)),"resampleImage(): Output images has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        resampleImage(srcImageRange(bimage), destImage(bres), factor);
    }

    return res;
}

enum RotationDirection
{
    ROTATE_CW,
    ROTATE_CCW,
    UPSIDE_DOWN
};

template < class PixelType>
NumpyAnyArray pythonFixedRotateImage(NumpyArray<3, Multiband<PixelType> > image, RotationDirection dir, NumpyArray<3,Multiband<PixelType> > res)
{
    int degree=0;
    switch(dir)
    {
    case ROTATE_CW:
        degree=270;
        break;
    case ROTATE_CCW:
        degree=90;
        break;
    case UPSIDE_DOWN:
        degree=180;
        break;
    }
    //make the output
    if(degree % 180 == 0)
        res.reshapeIfEmpty(image.shape(),"rotateImageSimple(): Output images has wrong dimensions");
    else
        res.reshapeIfEmpty(MultiArrayShape<3>::type(image.shape(1),image.shape(0),image.shape(2)),"rotateImage(): Output image has wrong dimensions");
    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        rotateImage(srcImageRange(bimage),destImage(bres),degree);
    }
    return res;
}

template < class PixelType>
NumpyAnyArray pythonFreeRotateImageDegree(NumpyArray<3, Multiband<PixelType> > image, double degree,RotationDirection dir, int splineOrder,NumpyArray<3,Multiband<PixelType> > res)
{
    return pythonFreeRotateImageRadiant(image,degree*M_PI/180.0,dir,splineOrder,res);
}

template < class PixelType>
NumpyAnyArray pythonFreeRotateImageRadiant(NumpyArray<3, Multiband<PixelType> > image, double radiant, RotationDirection dir, int splineOrder,NumpyArray<3,Multiband<PixelType> > res)
{
    //reshape, if empty. Otherwise accept res dimensions
    if(!res.hasData())
        res.reshapeIfEmpty(image.shape(),"");
    vigra_precondition(res.shape(2)==image.shape(2),"rotateImageRadiant(): number of channels of image and result have to be equal");
    //res.init(NumericTraits< PixelType>::zero());
    if(dir==ROTATE_CW)
        radiant=-radiant;

    //Define the transformation
    linalg::TemporaryMatrix< double >  transform=translationMatrix2D(TinyVector<double,2>(res.shape(0)/2.0,res.shape(1)/2.0))*
        rotationMatrix2DRadians(radiant,TinyVector<double,2>(0.0,0.0))*
        translationMatrix2D(TinyVector<double,2>(-image.shape(0)/2.0,-image.shape(1)/2.0));

    
    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        switch (splineOrder)
        {
        case 0:
            {
                SplineImageView< 0, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
                break;
            }
        case 1:
            {
                SplineImageView< 1, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
                break;
            }
        case 2:
            {
                SplineImageView< 2, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
            }
            break;
        case 3:
            { 
                SplineImageView< 3, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
                break;
            }
        case 4:
            {
                SplineImageView< 4, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
                break;
            }
        case 5:
            {
                SplineImageView< 5, PixelType > spline(srcImageRange(bimage));
                affineWarpImage(spline,destImageRange(bres),transform);
                break;
            }
        default:
            PyErr_SetString(PyExc_ValueError, "Spline order not supported.");
            python::throw_error_already_set();
        }
    }
    return res;
}

template < class PixelType>
NumpyAnyArray pythonResizeImageNoInterpolation(NumpyArray<3, Multiband<PixelType> > image, 
                                               python::object destSize,
                                               NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((image.shape(0) > 1) && (image.shape(1) > 1),
        "The input image must have a size of at least 2x2.");
    vigra_precondition((destSize!=python::object() && !res.hasData()) || (destSize==python::object() && res.hasData()),
                       "destSize or out has to be given, but only one of them");
    MultiArrayShape<2>::type size;
    if(!res.hasData())
    {
        size=python::extract<MultiArrayShape<2>::type>(destSize)();
    }
    else
    {
        size[0]=res.shape(0);
        size[1]=res.shape(1);
    }
    
    res.reshapeIfEmpty( MultiArrayShape<3>::type(size[0],size[1],image.shape(2)),
                        "Output image has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        resizeImageNoInterpolation(srcImageRange(bimage),destImageRange(bres));
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonResizeImageLinearInterpolation(NumpyArray<3, Multiband<PixelType> > image,
                                                   python::object destSize,
                                                   NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((image.shape(0) > 1) && (image.shape(1) > 1),
        "The input image must have a size of at least 2x2.");

    vigra_precondition((destSize!=python::object() && !res.hasData()) || (destSize==python::object() && res.hasData()),
                       "destSize or out has to be given, but only one of them");
    MultiArrayShape<2>::type size;
    if(!res.hasData())
    {
        size=python::extract<MultiArrayShape<2>::type>(destSize)();
    }
    else
    {
        size[0]=res.shape(0);
        size[1]=res.shape(1);
    }
    
    res.reshapeIfEmpty( MultiArrayShape<3>::type(size[0],size[1],image.shape(2)),
                        "Output image has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        resizeImageLinearInterpolation(srcImageRange(bimage), destImageRange(bres));
    }
    return res;
}

template < class PixelType, int dim >
NumpyAnyArray pythonResizeImageSplineInterpolation(NumpyArray<dim, Multiband<PixelType> > image,
                                                   python::object destSize,
                                                   int splineOrder=3, NumpyArray<dim, Multiband<PixelType> > res=python::object())
{
    vigra_precondition((image.shape(0) > 1) && (image.shape(1) > 1),
        "The input image must have a size of at least 2x2.");
    vigra_precondition((destSize!=python::object() && !res.hasData()) || (destSize==python::object() && res.hasData()),
                       "destSize or out has to be given, but only one of them");

    TinyVector<UInt32,dim> out_shape;
    if(!res.hasData())
    {
        //TinyVector<UInt32,dim-1> size;
        //size=python::extract<TinyVector<UInt32,dim-1> >(destSize)();
        typedef typename MultiArrayShape<dim-1>::type shape;
        shape size;
        size=python::extract<shape>(destSize)();
        for(int ii=0;ii<dim-1;++ii)
            out_shape[ii]=size[ii];
        out_shape[dim-1]=image.shape(dim-1);
    }
    else
    {
        for(int ii=0;ii<dim;++ii)
            out_shape[ii]=res.shape(ii);
    }
    
    res.reshapeIfEmpty(out_shape, "Output image has wrong dimensions");

    for(int k=0;k<image.shape(dim-1);++k)
    {
        
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<dim-1, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        switch (splineOrder)
        {
        case 1:
            {
                BSpline< 1, double > spline;
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(bimage),
                                               destMultiArrayRange(bres), spline);
                break;
            }
        case 2:
        {
            BSpline< 2, double > spline;
            resizeMultiArraySplineInterpolation(srcMultiArrayRange(bimage),
                destMultiArrayRange(bres), spline);
            break;
        }
        case 3:
        {
            BSpline< 3, double > spline;
            resizeMultiArraySplineInterpolation(srcMultiArrayRange(bimage),
                destMultiArrayRange(bres), spline);
            break;
        }
        case 4:
        {
            BSpline< 4, double > spline;
            resizeMultiArraySplineInterpolation(srcMultiArrayRange(bimage),
                destMultiArrayRange(bres), spline);
            break;
        }
        case 5:
        {
            BSpline< 5, double > spline;
            resizeMultiArraySplineInterpolation(srcMultiArrayRange(bimage),
                destMultiArrayRange(bres), spline);
            break;
        }
        default:
        {
            PyErr_SetString(PyExc_ValueError, "Spline order not supported.");
            python::throw_error_already_set();
        }
        }
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonResizeImageCatmullRomInterpolation(NumpyArray<3, Multiband<PixelType> > image,
                                                       python::object destSize,
                                                       NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((image.shape(0) > 3) && (image.shape(1) > 3),
        "The input image must have a size of at least 4x4.");
    vigra_precondition((destSize!=python::object() && !res.hasData()) || (destSize==python::object() && res.hasData()),
                       "destSize or out has to be given, but only one of them");
    MultiArrayShape<2>::type size;
    if(!res.hasData())
    {
        size=python::extract<MultiArrayShape<2>::type>(destSize)();
    }
    else
    {
        size[0]=res.shape(0);
        size[1]=res.shape(1);
    }
    
    vigra_precondition((size[0] > 1) && (size[1] > 1),
        "The destination image must have a size of at least 2x2.");
    res.reshapeIfEmpty( MultiArrayShape<3>::type(size[0],size[1],image.shape(2)),
                        "Output image has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);

        resizeImageCatmullRomInterpolation(srcImageRange(bimage),destImageRange(bres));
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonResizeImageCoscotInterpolation(NumpyArray<3, Multiband<PixelType> > image,
                                                   python::object destSize,
                                                   NumpyArray<3, Multiband<PixelType> > res)
{
    vigra_precondition((image.shape(0) > 3) && (image.shape(1) > 3),
        "The input image must have a size of at least 4x4.");
    vigra_precondition((destSize!=python::object() && !res.hasData()) || (destSize==python::object() && res.hasData()),
                       "destSize or out has to be given, but only one of them");
    MultiArrayShape<2>::type size;
    if(!res.hasData())
    {
        size=python::extract<MultiArrayShape<2>::type>(destSize)();
    }
    else
    {
        size[0]=res.shape(0);
        size[1]=res.shape(1);
    }
    
    vigra_precondition((size[0] > 1) && (size[1] > 1),
        "The destination image must have a size of at least 2x2.");
    res.reshapeIfEmpty( MultiArrayShape<3>::type(size[0],size[1],image.shape(2)),
                        "Output image has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        resizeImageCoscotInterpolation(srcImageRange(bimage),destImageRange(bres));
    }
    return res;
}



void defineImageTransformations()
{
    using namespace python;

    enum_<RotationDirection>("RotationDirection")
        .value("CLOCKWISE",ROTATE_CW)
        .value("COUNTER_CLOCKWISE",ROTATE_CCW)
        .value("UPSIDE_DOWN",UPSIDE_DOWN);
    def("rotateImageRadiant",
        registerConverters(&pythonFreeRotateImageRadiant<float>),
        (arg("image"), arg("radiant"), arg("direction")=ROTATE_CW, arg("splineOrder")=0, arg("out")=object()),
        "Rotate an image by an arbitrary angle using splines for interpolation around its center.\n"
        "\n"
        "The angle may be given in radiant (parameter radiant).\n"
        "The parameter splineOrder indicates the order of the splines used for interpolation.\n"
        "If the \"out\" parameter is given, the image is cropped for it's dimensions. If the \"out\"\n"
        "parameter is not given, an output image with the same dimensions as the input image is created.\n" 
        );
    def("rotateImageDegree",
        registerConverters(&pythonFreeRotateImageDegree<float>),
        (arg("image"), arg("degree"), arg("direction")=ROTATE_CW, arg("splineOrder")=0, arg("out")=object()),
        "Rotate an image by an arbitrary angle using splines for interpolation around its center.\n"
        "\n"
        "The angle may be given in degree (parameter degree).\n"
        "The parameter splineOrder indicates the order of the splines used for interpolation.\n"
        "If the \"out\" parameter is given, the image is cropped for it's dimensions. If the \"out\"\n"
        "parameter is not given, an output image with the same dimensions as the input image is created.\n" 
        );
    def("rotateImageSimple",
        registerConverters(&pythonFixedRotateImage<float>),
        (arg("image"), arg("orientation")=ROTATE_CW,arg("out")=object()),
        "Rotate an image by a multiple of 90 degrees.\n"
        "\n"
        "The \"orientation\" parameter (which must be one of CLOCKWISE, COUNTER_CLOCKWISE and UPSIDE_DOWN\n"
        "indicates the rotation direction. The \"out\" parameter must, if given, have the according dimensions.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
        );

//    def("rotateImageAboutCenter",
//       &DUMMY_FUNCTION,               // also multiband
//        (arg("image"), arg("degrees"), arg("center"), arg("splineOrder")));

    def("resampleImage",
        registerConverters(&pythonResampleImage<float>),               // also multiband
        (arg("image"), arg("factor"),arg("out")=object()),
        "Resample an image by the factor \"factor\""
        "\n"
        "For details, see the vigra documentation, resampleImage_.\n"
        "The \"out\" parameter must have, if given, the according dimensions.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
        );

    def("resizeImageNoInterpolation",
        registerConverters(&pythonResizeImageNoInterpolation<float>),               // also multiband
        (arg("image"), arg("destSize")=object(), arg("out")=object()),
        "Resize image by repeating the nearest pixel values.\n"
        "\n"
        "For details see vigra documentation resizeImageNoInterpolation_.\n"
        "The dimensions of the output image is taken either from \"destSize\" or \"out\".\n"
        "If both are given, they must agree on the output dimensions.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
        );


    def("resizeImageLinearInterpolation",
        registerConverters(&pythonResizeImageLinearInterpolation<float>),               // also multiband>
        (arg("image"), arg("destSize")=object(), arg("out")=object()),
        "Resize image using linear interpolation.\n"
        "The function uses the standard separable bilinear interpolation algorithm to obtain a good compromise between quality and speed.\n"
        "\n" 
        "For details see vigra documentation resizeImageLinearInterpolation_.\n"
        "The dimensions of the output image is taken either from \"destSize\" or \"out\".\n"
        "If both are given, they must agree on the output dimensions.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
        );

    def("resizeImageSplineInterpolation",
        registerConverters(&pythonResizeImageSplineInterpolation<float,3>),               // also multiband
        (arg("image"), arg("destSize")=object(), arg("splineOrder") = 3, arg("out") = object()),
        "Resize image using B-spline interpolation.\n"
        "\n"
        "For details see vigra documentation resizeImageSplineInterpolation_.\n"
        "The spline order is given in the parameter \"splineOrder\".\n"
        "The dimensions of the output image is taken either from \"destSize\" or \"out\".\n"
        "If both are given, they must agree on the output dimensions.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
       );

    def("resizeImageCatmullRomInterpolation",
        registerConverters(&pythonResizeImageCatmullRomInterpolation<float>),               // also multiband
        (arg("image"), arg("destSize")=object(), arg("out")=object()),
        "Resize image using the Catmull/Rom interpolation function.\n"
        "\n" 
        "For details see vigra documentation resizeImageCatmullRomInterpolation_.\n"
        "The dimensions of the output image is taken either from \"destSize\" or \"out\".\n"
        "If both are given, they must agree on the output dimensions.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
       );

    def("resizeImageCoscotInterpolation",
        registerConverters(&pythonResizeImageCoscotInterpolation<float>),               // also multiband
        (arg("image"), arg("destSize")=object(), arg("out")=object()),
        "Resize image using the Coscot interpolation function.\n" 
        "\n" 
        "For details see vigra documentation resizeImageCoscotInterpolation_.\n"
        "The dimensions of the output image is taken either from \"destSize\" or \"out\".\n"
        "If both are given, they must agree on the output dimensions.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
        );


}
} // namespace vigra
