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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyfourier_PyArray_API
#include <vigra/config.hxx>
#include <Python.h>
#include <boost/python.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/fftw3.hxx>
#include <vigra/multi_fft.hxx>
#include <vigra/gaborfilter.hxx>
#include <iostream>

namespace vigra {

template <>
struct NumpyArrayValuetypeTraits<FFTWComplex<float> >
{
    static bool isValuetypeCompatible(PyArrayObject const * obj) /* obj must not be NULL */
    {
        return PyArray_EquivTypenums(NPY_CFLOAT, PyArray_DESCR((PyObject *)obj)->type_num) &&
               PyArray_ITEMSIZE((PyObject *)obj) == sizeof(FFTWComplex<float>);
    }

    static NPY_TYPES const typeCode = NPY_CFLOAT;

    static std::string typeName()
    {
        return "complex64";
    }

    static std::string typeNameImpex()
    {
        return "";
    }

    static PyObject * typeObject()
    {
        return PyArray_TypeObjectFromType(NPY_CFLOAT);
    }
};

template <class T>
NumpyAnyArray
pythonCreateGaborFilter(typename MultiArrayView<2,T>::difference_type shape,
                        double orientation, 
                        double centerFrequency,
                        double angularSigma, 
                        double radialSigma, 
                        NumpyArray<2,Singleband<T> > res)
{
    res.reshapeIfEmpty(TaggedShape(shape, NumpyAnyArray::defaultAxistags(3)).toFrequencyDomain(), 
            "createGaborFilter(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        createGaborFilter(destImageRange(res),
                          orientation, centerFrequency, angularSigma, radialSigma);
    }
    return res;
}

// template <unsigned int N, int SIGN>
// NumpyAnyArray
// pythonFourierTransform(NumpyArray<N, Multiband<FFTWComplex<float> > > in, 
                       // NumpyArray<N, Multiband<FFTWComplex<float> > > res)
// {
    // res.reshapeIfEmpty(in.shape(), in.strideOrdering(),
        // "fourierTransform(): Output array must have the same shape and stride ordering as input array.", true);

    // for(MultiArrayIndex k=0; k<in.shape(N-1); ++k)
    // {
        // MultiArrayView<N-1, FFTWComplex<float> , StridedArrayTag> bin = in.bindOuter(k).permuteStridesDescending();
        // MultiArrayView<N-1, FFTWComplex<float> , StridedArrayTag> bres = res.bindOuter(k).permuteStridesDescending();

        // TinyVector<int, N-1> bshape(bin.shape()), itotal(bin.shape()), ototal(bres.shape());
        // float norm = (float)bshape[0];
        // for(int j=1; j<(int)N-1; ++j)
        // {
            // itotal[j] = bin.stride(j-1) / bin.stride(j);
            // ototal[j] = bres.stride(j-1) / bres.stride(j);
            // norm *= (float)bshape[j];
        // }

        // fftwf_plan plan = fftwf_plan_many_dft(N-1, bshape.begin(), 1,
                                            // (fftwf_complex*)bin.data(), itotal.begin(),
                                            // bin.stride(N-2), 0,
                                            // (fftwf_complex*)bres.data(), ototal.begin(),
                                            // bres.stride(N-2), 0,
                                            // SIGN, FFTW_ESTIMATE);
        // vigra_postcondition(plan != 0, "fourierTransform(): Unable to create plan.");
        // fftwf_execute(plan);
        // fftwf_destroy_plan(plan);
        // if(SIGN == FFTW_BACKWARD)
        // {
            // bres *= FFTWComplex<float>(1.0f / norm);
        // }
    // }
    // return res;
// }

template <unsigned int N, int SIGN>
NumpyAnyArray
pythonFourierTransform(NumpyArray<N, Multiband<FFTWComplex<float> > > in, 
                       NumpyArray<N, Multiband<FFTWComplex<float> > > res)
{
    res.reshapeIfEmpty(in.taggedShape().toFrequencyDomain(SIGN == FFTW_FORWARD
                                                               ? 1
                                                               : -1),
            "fourierTransform(): Output has wrong shape.");

    {
        PyAllowThreads _pythread;
        FFTWPlan<N-1, float> plan(in.bindOuter(0), res.bindOuter(0), SIGN);
        
        for(MultiArrayIndex k=0; k<in.shape(N-1); ++k)
        {
            plan.execute(in.bindOuter(k), res.bindOuter(k));
        }
    }
    return res;
}

// NumpyAnyArray
// pythonFourierTransformR2C(NumpyAnyArray in, NumpyAnyArray res)
// {
    // switch(in.spatialDimensions())
    // {
      // case 2:
      // {
        // NumpyArray<3, Multiband<FFTWComplex<float> > > inc(in, true), out(res, false);
        // return pythonFourierTransform<3, FFTW_FORWARD>(inc, out);
      // }
      // case 3:
      // {
        // NumpyArray<4, Multiband<FFTWComplex<float> > > inc(in, true), out(res, false);
        // return pythonFourierTransform<4, FFTW_FORWARD>(inc, out);
      // }
      // default:
        // vigra_fail("fourierTransform(): "
                   // "Can only handle 2 or 3 spatial dimensions.");
    // }
    // return res; // this will never be reached
// }

// FIXME: implement the correct R2C transform (is already provided in multi_fft.hxx)
// FIXME: numpy arrays are not allocated with fftw_malloc() - will this cause alignment
//        problems (sudden crashes)?
template <unsigned int N>
NumpyAnyArray
pythonFourierTransformR2C(NumpyArray<N, Multiband<float> > in, 
                          NumpyArray<N, Multiband<FFTWComplex<float> > > res)
{
    res.reshapeIfEmpty(in.taggedShape().toFrequencyDomain(),
            "fourierTransformR2C(): Output has wrong shape.");
        
    {
        PyAllowThreads _pythread;
        // static_cast<typename NumpyArray<N, Multiband<FFTWComplex<float> > >::view_type &>(res) = 
            // static_cast<typename NumpyArray<N, Multiband<float> >::view_type const &>(in);
        res = in;
        FFTWPlan<N-1, float> plan(res.bindOuter(0), res.bindOuter(0), FFTW_FORWARD);
        
        for(MultiArrayIndex k=0; k<res.shape(N-1); ++k)
        {
            plan.execute(res.bindOuter(k), res.bindOuter(k));
        }
    }
    return res;
}

} // namespace vigra

using namespace vigra;

BOOST_PYTHON_MODULE_INIT(fourier)
{
    using boost::python::arg;
    using boost::python::docstring_options;
    using boost::python::object;
    using boost::python::def;

    import_vigranumpy();

    docstring_options doc_options(true, true, false);

    def("fourierTransform", registerConverters(&pythonFourierTransformR2C<2>),
        (arg("image"), arg("out") = object()),
        "Perform 2-dimensional Fourier transformation of a scalar float32 image."
        "If the input array has multiple channels, each channel is transformed separately.\n"
        );

    def("fourierTransform", registerConverters(&pythonFourierTransformR2C<3>),
        (arg("volume"), arg("out") = object()),
        "Likewise for a 3D float32 volume.\n"
        );

    def("fourierTransform", registerConverters(&pythonFourierTransform<3, FFTW_FORWARD>),
        (arg("image"), arg("out") = object()),
        "Likewise for a 2D complex64 image.\n");

    def("fourierTransform", registerConverters(&pythonFourierTransform<4, FFTW_FORWARD>),
        (arg("volume"), arg("out") = object()),
        "Likewise for a 3D complex64 volume.\n");

    def("fourierTransformInverse", registerConverters(&pythonFourierTransform<3, FFTW_BACKWARD>),
        (arg("image"), arg("out") = object()),
        "Perform 2-dimensional inverse Fourier transformation of a complex64 array."
        "If the input array has multiple channels, each channel is transformed separately.\n");

    def("fourierTransformInverse", registerConverters(&pythonFourierTransform<4, FFTW_BACKWARD>),
        (arg("volume"), arg("out") = object()),
        "Likewise for a 3D complex128 volume.\n");

    def("createGaborFilter", registerConverters(&pythonCreateGaborFilter<float>),
        (arg("shape"),arg("orientation"),arg("centerFrequency"),arg("angularSigma"),arg("radialSigma"), arg("out") = object()),
        "Create a 2-dimensional gabor filter in frequency space.");

    def("radialGaborSigma", &radialGaborSigma,
        "Calculate sensible radial sigma for given parameters.");

    def("angularGaborSigma", &angularGaborSigma,
        "Calculate sensible angular sigma for given parameters.");
}
