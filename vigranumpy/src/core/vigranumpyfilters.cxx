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
#include <vigra/affinegeometry.hxx>
#include <vigra/basicgeometry.hxx>
#include <vigra/resizeimage.hxx>
#include <vigra/splines.hxx>
#include <vigra/flatmorphology.hxx>
#include <vigra/noise_normalization.hxx>
#include <vigra/multi_morphology.hxx>
#include <vigra/multi_resize.hxx>

#include <cmath>

#include "tensors.hxx"
namespace python = boost::python;

namespace vigra
{
namespace { void DUMMY_FUNCTION(int, int, int, int, int, int, int, int, int, int) {} }

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



//TODO
NumpyAnyArray vectorToArray(std::vector< TinyVector< double, 2 > > & result)
{
    NumpyArray<2,double> res(MultiArrayShape<2>::type(result.size(),2));

    for(size_t ii=0;ii<result.size();++ii)
    {
        res(ii,0)=result[ii][0];
        res(ii,1)=result[ii][1];
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonNoiseVarianceEstimation(NumpyArray<2, Singleband<PixelType> > image,
                                            bool useGradient=true, unsigned int windowRadius=6,
                                            unsigned int clusterCount=10, double averagingQuantile=0.8,
                                            double noiseEstimationQuantile=1.5, double noiseVarianceInitialGuess=10.0,
                                            NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    NoiseNormalizationOptions noiseNormalizationOptions;
    noiseNormalizationOptions.useGradient(useGradient).windowRadius(windowRadius).clusterCount(clusterCount).averagingQuantile(averagingQuantile).noiseEstimationQuantile(noiseEstimationQuantile).noiseVarianceInitialGuess(noiseVarianceInitialGuess);
    std::vector< TinyVector< double, 2 > > result;

    noiseVarianceEstimation(srcImageRange(image), result, noiseNormalizationOptions);

    return vectorToArray(result);
}

template < class PixelType >
NumpyAnyArray pythonNoiseVarianceClustering(NumpyArray<2, Singleband<PixelType> > image,
                                            bool useGradient=true, unsigned int windowRadius=6,
                                            unsigned int clusterCount=10, double averagingQuantile=0.8,
                                            double noiseEstimationQuantile=1.5, double noiseVarianceInitialGuess=10.0,
                                            NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    NoiseNormalizationOptions noiseNormalizationOptions;
    noiseNormalizationOptions.useGradient(useGradient).windowRadius(windowRadius).clusterCount(clusterCount).averagingQuantile(averagingQuantile).noiseEstimationQuantile(noiseEstimationQuantile).noiseVarianceInitialGuess(noiseVarianceInitialGuess);
    std::vector< TinyVector< double, 2 > > result;
    noiseVarianceClustering(srcImageRange(image), result,
        noiseNormalizationOptions);
    return vectorToArray(result);
}

template < class PixelType >
NumpyAnyArray pythonNonparametricNoiseNormalization(NumpyArray<3, Multiband<PixelType> > image,
                                                    bool useGradient=true, unsigned int windowRadius=6,
                                                    unsigned int clusterCount=10, double averagingQuantile=0.8,
                                                    double noiseEstimationQuantile=1.5, double noiseVarianceInitialGuess=10.0,
                                                    NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    NoiseNormalizationOptions noiseNormalizationOptions;
    noiseNormalizationOptions.useGradient(useGradient).windowRadius(windowRadius).clusterCount(clusterCount).averagingQuantile(averagingQuantile).noiseEstimationQuantile(noiseEstimationQuantile).noiseVarianceInitialGuess(noiseVarianceInitialGuess);
    
    res.reshapeIfEmpty(image.shape(),"nonparametricNoiseNormalization(): Output images has wrong dimensions");
    
    for(int k=0;k<image.shape(2);++k)
    {
        nonparametricNoiseNormalization(srcImageRange(image),
                                        destImage(res), noiseNormalizationOptions);
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonQuadraticNoiseNormalizationEstimated( NumpyArray<3, Multiband<PixelType> > image,
                                                          bool useGradient=true, unsigned int windowRadius=6,
                                                          unsigned int clusterCount=10, double averagingQuantile=0.8,
                                                          double noiseEstimationQuantile=1.5, double noiseVarianceInitialGuess=10.0,
                                                          NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    NoiseNormalizationOptions noiseNormalizationOptions;
    noiseNormalizationOptions.useGradient(useGradient).windowRadius(windowRadius).clusterCount(clusterCount).averagingQuantile(averagingQuantile).noiseEstimationQuantile(noiseEstimationQuantile).noiseVarianceInitialGuess(noiseVarianceInitialGuess);

    res.reshapeIfEmpty(image.shape(),"quadraticNoiseNormalizationEstimated(): Output images has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        quadraticNoiseNormalization(srcImageRange(bimage),
                                    destImage(bres), noiseNormalizationOptions);
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonLinearNoiseNormalizationEstimated(NumpyArray<3, Multiband<PixelType> > image,
                                                      bool useGradient=true, unsigned int windowRadius=6,
                                                      unsigned int clusterCount=10, double averagingQuantile=0.8,
                                                      double noiseEstimationQuantile=1.5, double noiseVarianceInitialGuess=10.0,
                                                      NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    NoiseNormalizationOptions noiseNormalizationOptions;
    noiseNormalizationOptions.useGradient(useGradient).windowRadius(windowRadius).clusterCount(clusterCount).averagingQuantile(averagingQuantile).noiseEstimationQuantile(noiseEstimationQuantile).noiseVarianceInitialGuess(noiseVarianceInitialGuess);
    
    res.reshapeIfEmpty(image.shape(),"linearNoiseNormalizationEstimated(): Output images has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        linearNoiseNormalization(srcImageRange(bimage),
                                 destImage(bres), noiseNormalizationOptions);
    }
    return res;
}


template < class PixelType >
NumpyAnyArray pythonQuadraticNoiseNormalization(NumpyArray<3, Multiband<PixelType> > image,
                                                double a0, double a1, double a2,
                                                NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    res.reshapeIfEmpty(image.shape(),"quadraticNoiseNormalization(): Output images has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        quadraticNoiseNormalization(srcImageRange(bimage), destImage(bres),a0, a1, a2);
    }
    
    return res;
}

template < class PixelType >
NumpyAnyArray pythonLinearNoiseNormalization(NumpyArray<3, Multiband<PixelType> > image,
                                             double a0, double a1, NumpyArray<3, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(image.shape(),"linearNoiseNormalization(): Output images has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        linearNoiseNormalization(srcImageRange(bimage), destImage(bres),a0, a1);
    } 
    return res;
}


template < class PixelType>
NumpyAnyArray pythonEigHessian2d(NumpyArray<2,Singleband<PixelType> > image, double scale,
                                 NumpyArray<3,Singleband<PixelType>  > res)
{
    //reshape if empty
    res.reshapeIfEmpty(MultiArrayShape<3>::type(image.shape(0),image.shape(1),2));

    MultiArray<3,PixelType> hessian(MultiArrayShape<3>::type(image.shape(0),image.shape(1),2));

    hessianOfGaussian(image, hessian, scale);

    eigenValuesPerPixel(hessian,res);

    return res;
}

template<class PixelType>
NumpyAnyArray pythonEigStructureTensor2d(NumpyArray<2,Singleband<PixelType> > image, double innerScale,double outerScale,
                                         NumpyArray<3,Singleband<PixelType> > res)
{
    //reshape output
    res.reshapeIfEmpty(MultiArrayShape<3>::type(image.shape(0),image.shape(1),2));

    MultiArray<3,PixelType> st(MultiArrayShape<3>::type(image.shape(0),image.shape(1), 3));

    structureTensor(image,st,innerScale,outerScale);

    eigenValuesPerPixel(st,res);
    return res;
}

template<class PixelType>
NumpyAnyArray pythonGaussianSmooth2d(NumpyArray<2,Singleband<PixelType> > image, double scale,
                                     NumpyArray<2,Singleband<PixelType> > res)
{
    //reshape output
    res.reshapeIfEmpty(image.shape());

    gaussianSmoothMultiArray(srcMultiArrayRange(image), destMultiArray(res), scale);
   

    return res;
}


void defineFilters2D()
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
        "This is an abbreviation for the rank order filter with rank = 0.0. See discRankOrderFilter_ for more information.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
        );
    def("discErosion",
        registerConverters(&pythonDiscErosion<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()));

    def("discDilation",
        registerConverters(&pythonDiscDilation<bool>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply dilation (maximum) filter with disc of given radius to image.\n\n"
        "This is an abbreviation for the rank order filter with rank = 1.0. See discRankOrderFilter_ for more information.\n"
        "This function also works for multiband images, it is then executed on every band.\n"
       );
    def("discDilation",
        registerConverters(&pythonDiscDilation<UInt8>),
        (arg("image"), arg("radius"), arg("out")=object()));

    def("discMedian",
        registerConverters(&pythonDiscMedian<bool>),
        (arg("image"), arg("radius"), arg("out")=object()),
        "Apply median filter with disc of given radius to image.\n\n"
        "This is an abbreviation for the rank order filter with rank = 0.5. See discRankOrderFilter_ for more information.\n"
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
        "See discRankOrderFilter_ in the vigra C++ documentation for more information\n"
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
        "For details see vigra documentation (multiBanariyDilation_ and multiBinaryErosion_).\n"
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
        "For details see vigra documentation (multiBanariyDilation_ and multiBinaryErosion_).\n"
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

    def("noiseVarianceEstimation",
        registerConverters(&pythonNoiseVarianceEstimation<float>),
        (arg("image"), arg("useGradient")=true, arg("windowRadius")=6,
        arg("clusterCount")=10, arg("averagingQuantile")=0.8,
        arg("noiseEstimationQuantile")=1.5,
        arg("noiseVarianceInitialGuess")=10.0, arg("out")=object()),
        "Determine the noise variance as a function of the image intensity.\n"
        "\n"
        "Returns an array with the means in the first column and the variances in the second column.\n"
        "Since the length of the resulting array is not known beforhand, it can not be written into an preallocated array\n"
        "(the \"out\" argument in most other vigra python functions.\n"
        "For details see the vigra documentation noiseVarianceEstimation_."
        );

    def("noiseVarianceClustering",
        registerConverters(&pythonNoiseVarianceClustering<float>),
        (arg("image"), arg("useGradient")=true, arg("windowRadius")=6,
        arg("clusterCount")=10, arg("averagingQuantile")=0.8,
        arg("noiseEstimationQuantile")=1.5,
        arg("noiseVarianceInitialGuess")=10.0, arg("out")=object()),
        "Determine the noise variance as a function of the image intensity and cluster the results.\n"
        "This operator first calls noiseVarianceEstimation() to obtain a sequence of intensity/variance pairs,\n"
        "which are then clustered using the median cut algorithm. Then the cluster centers (i.e. average variance vs. average intensity)\n"
        "are determined and returned in the result sequence.\n"
        "\n"
        "Since the length of the resulting array is not known beforhand, it can not be written into an preallocated array\n"
        "(the \"out\" argument in most other vigra python functions.\n"
        "For details see the vigra documentation noiseVarianceClustering_."
        );

    def("nonparametricNoiseNormalization",
        registerConverters(&pythonNonparametricNoiseNormalization<float>),    // also multiband
        (arg("image"), arg("useGradient")=true, arg("windowRadius")=6,
        arg("clusterCount")=10, arg("averagingQuantile")=0.8,
        arg("noiseEstimationQuantile")=1.5,
        arg("noiseVarianceInitialGuess")=10.0, arg("out")=object()),
        "Noise normalization by means of an estimated non-parametric noise model.\n"
        "For details see the vigra documentation nonparametricNoiseNormalization_.");

    def("quadraticNoiseNormalizationEstimated",
        registerConverters(&pythonQuadraticNoiseNormalizationEstimated<float>),    // also multiband
        (arg("image"), arg("useGradient")=true, arg("windowRadius")=6,
        arg("clusterCount")=10, arg("averagingQuantile")=0.8,
        arg("noiseEstimationQuantile")=1.5,
        arg("noiseVarianceInitialGuess")=10.0, arg("out")=object()));

    def("linearNoiseNormalizationEstimated",
        registerConverters(&pythonLinearNoiseNormalizationEstimated<float>),    // also multiband
        (arg("image"), arg("useGradient")=true, arg("windowRadius")=6,
        arg("clusterCount")=10, arg("averagingQuantile")=0.8,
        arg("noiseEstimationQuantile")=1.5,
        arg("noiseVarianceInitialGuess")=10.0, arg("out")=object()));

    def("quadraticNoiseNormalization",
        registerConverters(&pythonQuadraticNoiseNormalization<float>),    // also multiband
        (arg("image"), arg("a0"), arg("a1"), arg("a2"), arg("out")=object()),
        "Noise normalization by means of an estimated quadratic noise model.\n"
        "For details see the vigra documentation."
        );

    def("linearNoiseNormalization",
        registerConverters(&pythonLinearNoiseNormalization<float>),    // also multiband
        (arg("image"), arg("a0"), arg("a1"), arg("out")=object()),
        "Noise normalization by means of an estimated linear noise model.\n"
        "For details see the vigra documentaton."
       );
    def("eigHessian2d",
        registerConverters(&pythonEigHessian2d<float>),
        (arg("image"), arg("scale"), arg("out")=object()),
        "computes the eigenvalues of the Hessian in each pixel of an image.\n"
        "\n"
        "The 'scale' parameter gives the scale of 2nd derivative filter.\n"
       );
    def("eigStructureTensor2d",
        registerConverters(&pythonEigStructureTensor2d<float>),
        (arg("image"), arg("scale"), arg("out")=object()),
        "computes the eigenvalues of the structure tensor in each pixel of an image.\n"
        "\n"
        "The innerScale parameter gives the scale of the gradient filter\n"
        "The outerScale parameter gives the outer scale of the smoothing filter\n"
       );
    def("gaussianSmooth2d",
        registerConverters(&pythonGaussianSmooth2d<float>),
        (arg("image"), arg("scale"), arg("out")=object()),
        "computes the gradient magnitude in each pixel of an image.\n"
        "\n"
        "The scale parameter gives the scale of the derivative filter"
       );
}

} // namespace vigra

