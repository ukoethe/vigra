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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpytest_PyArray_API
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/signature.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <iostream>

namespace python = boost::python;

namespace vigra {

python::tuple testAny(NumpyAnyArray array)
{
    NumpyAnyArray copy(array, true);
    return python::make_tuple(array.shape(), copy, python::object(), python::object());
}

template <unsigned int N, class T, class Stride>
python::tuple test(NumpyArray<N, T, Stride> const & array)
{
    NumpyAnyArray anyarray(array);
    
    NumpyArray<N, T> copy(array, true);
    vigra_postcondition(copy.pyObject()->ob_refcnt == 1, 
          "freshly created NumpyArray<N, T> has reference count > 1.");

    NumpyArray<N, T> same_shape(array.shape());
    same_shape = array;

    NumpyArray<N, T> same_shape_and_tags;
    same_shape_and_tags.reshapeIfEmpty(array.taggedShape());
    same_shape_and_tags = anyarray;

    return python::make_tuple(anyarray.shape(), copy, 
                               same_shape, same_shape_and_tags);
}

template <unsigned int N, class T, class Stride>
typename MultiArrayShape<N>::type 
testView(MultiArrayView<N, T, Stride> array)
{
    array[typename MultiArrayShape<N>::type()] = NumericTraits<T>::one();
    return array.shape();
}

#if 0 // FIXME: temporarily disabled
// (right now, a compile-only test:)
void testMakeReference()
{
    MultiArray<2, vigra::UInt8> cpp_memory(MultiArrayShape<2>::type(100, 100));

    NumpyArray<2, npy_uint8, vigra::UnstridedArrayTag> python_view;
    python_view.makeReference(cpp_memory);
}
#endif

template <unsigned int N, class T>
python::tuple 
checkTaggedShape(NumpyArray<N, T> in)
{
    NumpyArray<3, Multiband<float> > res1(in);
    
    NumpyArray<3, Multiband<float> > res2;
    res2.reshapeIfEmpty(in.taggedShape().setChannelDescription("res2"), "res2 failed");
    
    NumpyArray<3, Multiband<float> > res3;
    res3.reshapeIfEmpty(in.taggedShape().setChannelCount(1).setChannelDescription("res3"), 
                        "res3 failed");
    
    NumpyArray<3, Multiband<float> > res4;
    res4.reshapeIfEmpty(in.taggedShape().setChannelCount(3).setChannelDescription("res4"), 
                        "res4 failed");
    
    NumpyArray<2, Singleband<float> > res5;
    res5.reshapeIfEmpty(in.taggedShape().setChannelDescription("res5"), "res5 failed");
    
    NumpyArray<2, RGBValue<float> > res6;
    res6.reshapeIfEmpty(in.taggedShape().setChannelDescription("res6"), 
                        "res6 failed");
    
    return python::make_tuple(res1, res2, res3, res4, res5, res6); 
}

} // namespace vigra

using namespace boost::python;
using namespace vigra;

BOOST_PYTHON_MODULE_INIT(vigranumpytest)
{
    import_vigranumpy();
    
    def("testAny", &testAny);

    def("testArray2Strided", registerConverters(&test<2, float, StridedArrayTag>));
    def("testImageSinglebandStrided", registerConverters(&test<2, Singleband<float>, StridedArrayTag>));
    def("testImageRGBStrided", registerConverters(&test<2, RGBValue<float>, StridedArrayTag >));
    def("testImageVector2Strided", registerConverters(&test<2, TinyVector<float, 2>, StridedArrayTag >));
    def("testImageMultibandStrided", registerConverters(&test<3, Multiband<float>, StridedArrayTag >));
    
    def("testArray2Unstrided", registerConverters(&test<2, float, UnstridedArrayTag>));
    def("testImageSinglebandUnstrided", registerConverters(&test<2, Singleband<float>, UnstridedArrayTag>));
    def("testImageRGBUnstrided", registerConverters(&test<2, RGBValue<float>, UnstridedArrayTag >));
    def("testImageVector2Unstrided", registerConverters(&test<2, TinyVector<float, 2>, UnstridedArrayTag >));
    def("testImageMultibandUnstrided", registerConverters(&test<3, Multiband<float>, UnstridedArrayTag >));

    def("testArray3Strided", registerConverters(&test<3, float, StridedArrayTag>));
    def("testVolumeSinglebandStrided", registerConverters(&test<3, Singleband<float>, StridedArrayTag>));
    def("testVolumeRGBStrided", registerConverters(&test<3, RGBValue<float>, StridedArrayTag >));
    def("testVolumeVector2Strided", registerConverters(&test<3, TinyVector<float, 2>, StridedArrayTag >));
    def("testVolumeMultibandStrided", registerConverters(&test<4, Multiband<float>, StridedArrayTag >));
    
    def("testArray3Unstrided", registerConverters(&test<3, float, UnstridedArrayTag>));
    def("testVolumeSinglebandUnstrided", registerConverters(&test<3, Singleband<float>, UnstridedArrayTag>));
    def("testVolumeRGBUnstrided", registerConverters(&test<3, RGBValue<float>, UnstridedArrayTag >));
    def("testVolumeVector2Unstrided", registerConverters(&test<3, TinyVector<float, 2>, UnstridedArrayTag >));
    def("testVolumeMultibandUnstrided", registerConverters(&test<4, Multiband<float>, UnstridedArrayTag >));

    def("testArray4Strided", registerConverters(&test<4, float, StridedArrayTag>));
    def("testArray4Unstrided", registerConverters(&test<4, float, UnstridedArrayTag>));

    def("viewArray2Strided", registerConverters(&testView<2, float, StridedArrayTag>));
    def("viewImageRGBStrided", registerConverters(&testView<2, RGBValue<float>, StridedArrayTag >));
    def("viewImageVector2Strided", registerConverters(&testView<2, TinyVector<float, 2>, StridedArrayTag >));
    
    def("viewArray2Unstrided", registerConverters(&testView<2, float, UnstridedArrayTag>));
    def("viewImageRGBUnstrided", registerConverters(&testView<2, RGBValue<float>, UnstridedArrayTag >));
    def("viewImageVector2Unstrided", registerConverters(&testView<2, TinyVector<float, 2>, UnstridedArrayTag >));

    def("viewArray3Strided", registerConverters(&testView<3, float, StridedArrayTag>));
    def("viewVolumeRGBStrided", registerConverters(&testView<3, RGBValue<float>, StridedArrayTag >));
    def("viewVolumeVector2Strided", registerConverters(&testView<3, TinyVector<float, 2>, StridedArrayTag >));
    
    def("viewArray3Unstrided", registerConverters(&testView<3, float, UnstridedArrayTag>));
    def("viewVolumeRGBUnstrided", registerConverters(&testView<3, RGBValue<float>, UnstridedArrayTag >));
    def("viewVolumeVector2Unstrided", registerConverters(&testView<3, TinyVector<float, 2>, UnstridedArrayTag >));    

    def("viewArray4Strided", registerConverters(&testView<4, float, StridedArrayTag>));
    def("viewArray4Unstrided", registerConverters(&testView<4, float, UnstridedArrayTag>));
    
    def("checkTaggedShapeMultiband", registerConverters(&checkTaggedShape<3, Multiband<float> >));
    def("checkTaggedShapeSingleband", registerConverters(&checkTaggedShape<2, Singleband<float> >));
}
