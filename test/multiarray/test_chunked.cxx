/************************************************************************/
/*                                                                      */
/*     Copyright 2013-2014 by Ullrich Koethe                            */
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

#include <stdio.h>

#include "vigra/unittest.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_array_chunked.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/multi_math.hxx"
#include "vigra/algorithm.hxx"
#include "vigra/random.hxx"
#include "vigra/timing.hxx"
//#include "marray.hxx"

using namespace vigra;
using namespace vigra::functor;

#define shouldEqualIndexing(N, a, b) \
{ \
    MultiCoordinateIterator<N> cccccccc(a.shape()), ccccccccend(cccccccc.getEndIterator()); \
    for(; cccccccc != ccccccccend; ++cccccccc) \
        if(a[*cccccccc] != b[*cccccccc]) \
            shouldEqual(a[*cccccccc], b[*cccccccc]); \
}

template <class Array>
class ChunkedMultiArrayTest
{
public:

    // typedef typename vigra::detail::ResolveMultiband<T>::type   scalar_type;
    typedef typename Array::value_type T;
    typedef MultiArray <3, T> plain;
    // typedef MultiArray <3, T> array3_type;
    // typedef MultiArrayView<3, Multiband<scalar_type> >  MultibandView3;
    // typedef typename array3_type::view_type       array3_view_type;
    // typedef typename array3_type::actual_stride   array3_stride;
    // typedef typename array2_type::difference_type difference2_type;
    // typedef typename array3_type::difference_type difference3_type;
    
    Shape3 shape, chunk_shape;
    VIGRA_UNIQUE_PTR<Array> array;
    plain ref;
    int cache_max;

    ChunkedMultiArrayTest ()
        : shape (20, 21, 22), chunk_shape (8),
          ref(shape),
          cache_max(9)
    {
        linearSequence(ref.begin(), ref.end());
        createArray(array);
        linearSequence(array->begin(), array->end());
    }
    
    void createArray(VIGRA_UNIQUE_PTR<ChunkedArrayFull<3, T> > & a)
    {
        a.reset(new ChunkedArrayFull<3, T>(shape));
    }
    
    void createArray(VIGRA_UNIQUE_PTR<ChunkedArrayLazy<3, T> > & a)
    {
        a.reset(new ChunkedArrayLazy<3, T>(shape, chunk_shape));
    }
    
    void createArray(VIGRA_UNIQUE_PTR<ChunkedArrayCompressed<3, T> > & a)
    {
        a.reset(new ChunkedArrayCompressed<3, T>(shape, cache_max, LZ4, chunk_shape));
    }
    
    void createArray(VIGRA_UNIQUE_PTR<ChunkedArrayHDF5<3, T> > & a)
    {
        HDF5File hdf5_file("chunked_test.h5", HDF5File::New);
        a.reset(new ChunkedArrayHDF5<3, T>(hdf5_file, "test", shape, cache_max, 1, chunk_shape));
    }
    
    void createArray(VIGRA_UNIQUE_PTR<ChunkedArrayTmpFile<3, T> > & a)
    {
        a.reset(new ChunkedArrayTmpFile<3, T>(shape, cache_max, "", chunk_shape));
    }

    // void testHasData ()
    // {
        // should(array3.hasData());
        // array3_type empty;
        // should(!empty.hasData());
        // array3_type empty2(Shape3(2,1,0));
        // should(!empty2.hasData());
    // }

    void test_assignment()
    {
        MultiArrayView <3, T, ChunkedArrayTag> v(shape);
        should(!v.hasData());    
        
        array->viewSubarray(Shape3(), v);
        should(v.hasData());        
        
        MultiArrayView <3, T, ChunkedArrayTag> vc;
        should(!vc.hasData());
        
        vc = v;
        should(vc.hasData());
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
        
        vc = T(7);
        std::vector<T> v7ref(vc.size(), T(7));
        should(vc.hasData());
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), v7ref.begin());
        shouldEqualSequence(v.begin(), v.end(), v7ref.begin());
        
        vc = ref;
        should(vc.hasData());
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
        
        MultiArrayView <3, T, ChunkedArrayTag> vs(Shape3(4));
        should(!vs.hasData());    
        array->viewSubarray(Shape3(), vs);
        should(vs.hasData());
        
        try
        {
            vc = vs;
            failTest("shape mismatch in assignment failed to throw exception");
        }
        catch(PreconditionViolation & e)
        {
            std::string expected("\nPrecondition violation!\nMultiArrayView::operator=(): shape mismatch.\n"),
                        actual(e.what());
            shouldEqual(actual.substr(0, expected.size()), expected);
        }
        
        vc += T(1);
        ref += T(1);
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
        
        vc += v;
        ref *= T(2);
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
         
        vc += T(42);
        ref += T(42);
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
        
        vc -= T(42);
        ref -= T(42);
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
       
        ref /= T(2);
        vc -= ref;
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
        
        vc *= v;
        ref *= ref;
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
        
        vc *= T(4);
        ref *= T(4);
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
         
        vc /= T(4);
        ref /= T(4);
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
       
        using namespace multi_math;
        ref = sqrt(ref);
        vc /= ref;
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
    }
    
    void testEquality ()
    {
        shouldEqual(array->shape(), ref.shape());
        
        shouldEqualSequence(array->begin(), array->end(), ref.begin());
        
        should(*array == ref);
        should(*array != ref.subarray(Shape3(1),ref.shape()));
        
        ++ref[ref.size()-1];
        should(*array != ref);
        
        // should(array3 != array3.subarray(Shape(1,1,1), Shape(2,2,2)));
        // should(array3.subarray(Shape(0,0,0), Shape(10,1,1)) != array3.subarray(Shape(0,1,0), Shape(10,2,1)));

        // array3_type a(array3.shape());
        // linearSequence(a.begin(), a.end());
        // should(a == array3);

        // for(unsigned int k=0; k<10; ++k)
            // array3(k,0,0) += 10;
        // should(array3.subarray(Shape(0,0,0), Shape(10,1,1)) == array3.subarray(Shape(0,1,0), Shape(10,2,1)));

        // MultibandView3 channel_view(a.multiband());
        // shouldEqual(a.shape(), channel_view.shape());
        // shouldEqual(a.data(), channel_view.data());
    }
    
    // // bindInner tests
    // void test_bindInner ()
    // {
        // TinyVector <int, 2> inner_indices (2, 5);
        // MultiArrayView <1, scalar_type, StridedArrayTag>
            // array = array3.bindInner(inner_indices);
        // shouldEqual ((array [TinyVector <int, 1> (0)]), 52);
        // shouldEqual ((array [TinyVector <int, 1> (1)]), 152);
        // shouldEqual ((array (0)), 52);
        // shouldEqual ((array (1)), 152);
    // }
    
    // // bindOuter tests
    // void test_bindOuter ()
    // {
        // TinyVector <int, 2> outer_indices (2, 5);
        // MultiArrayView <1, scalar_type, array3_stride>
            // array = array3.bindOuter(outer_indices);
        // shouldEqual ((array [TinyVector <int, 1> (0)]), 520);
        // shouldEqual ((array [TinyVector <int, 1> (1)]), 521);
        // shouldEqual ((array (0)), 520);
        // shouldEqual ((array (1)), 521);
    // }

    // bindAt tests
    void test_bindAt ()
    {
        MultiArrayView <2, T, ChunkedArrayTag> v = array->bindAt (1, 4);
        MultiArrayView <2, T, StridedArrayTag> vr = ref.bindAt (1, 4);
        
        shouldEqual(v.shape(), vr.shape());
        should(v == vr);
        shouldEqualSequence(v.begin(), v.end(), vr.begin());
        shouldEqualIndexing(2, v, vr);
        
        MultiArrayView <2, T, ChunkedArrayTag> vt = v.transpose();
        MultiArrayView <2, T, StridedArrayTag> vtr = vr.transpose();
        
        shouldEqual(vt.shape(), vtr.shape());
        should(vt == vtr);
        shouldEqualSequence(vt.begin(), vt.end(), vtr.begin());
        shouldEqualIndexing(2, vt, vtr);
        
        MultiArrayView <1, T, ChunkedArrayTag> v1 = v.bindAt (0, 11);
        MultiArrayView <1, T, StridedArrayTag> v1r = vr.bindAt (0, 11);
        
        shouldEqual(v1.shape(), v1r.shape());
        should(v1 == v1r);
        shouldEqualSequence(v1.begin(), v1.end(), v1r.begin());
        shouldEqualIndexing(1, v1, v1r);
        
        MultiArrayView <1, T, ChunkedArrayTag> v1t = v1.transpose();
        
        shouldEqual(v1t.shape(), v1r.shape());
        should(v1t == v1r);
        shouldEqualSequence(v1t.begin(), v1t.end(), v1r.begin());
        shouldEqualIndexing(1, v1t, v1r);
    }
    
    // // bind tests
    // void test_bind ()
    // {
        // MultiArrayView <2, scalar_type, array3_stride>
            // array = array3.template bind<1>(4);
        // shouldEqual ((array [TinyVector <int, 2> (0, 0)]), 40);
        // shouldEqual ((array [TinyVector <int, 2> (1, 0)]), 41);
        // shouldEqual ((array [TinyVector <int, 2> (0, 1)]), 140);
        // shouldEqual ((array (0, 0)), 40);
        // shouldEqual ((array (1, 0)), 41);
        // shouldEqual ((array (0, 1)), 140);
    // }

    // void test_bind0 ()
    // {
        // MultiArrayView <2, scalar_type, StridedArrayTag>
            // array = array3.template bind <0> (4);
        // shouldEqual ((array [TinyVector <int, 2> (0, 0)]), 4);
        // shouldEqual ((array [TinyVector <int, 2> (1, 0)]), 14);
        // shouldEqual ((array [TinyVector <int, 2> (0, 1)]), 104);
        // shouldEqual ((array (0, 0)), 4);
        // shouldEqual ((array (1, 0)), 14);
        // shouldEqual ((array (0, 1)), 104);
    // }

    // void test_singletonDimension ()
    // {
        // MultiArrayView <4, scalar_type, UnstridedArrayTag> a0 = array3.insertSingletonDimension(0);
        // shouldEqual ((a0 [TinyVector <int, 4> (0, 4, 0, 0)]), 4);
        // shouldEqual ((a0 [TinyVector <int, 4> (0, 4, 1, 0)]), 14);
        // shouldEqual ((a0 [TinyVector <int, 4> (0, 4, 0, 1)]), 104);

        // MultiArrayView <4, scalar_type, array3_stride> a1 = array3.insertSingletonDimension(1);
        // shouldEqual ((a1 [TinyVector <int, 4> (4, 0, 0, 0)]), 4);
        // shouldEqual ((a1 [TinyVector <int, 4> (4, 0, 1, 0)]), 14);
        // shouldEqual ((a1 [TinyVector <int, 4> (4, 0, 0, 1)]), 104);

        // MultiArrayView <4, scalar_type, array3_stride> a3 = array3.insertSingletonDimension(3);
        // shouldEqual ((a3 [TinyVector <int, 4> (4, 0, 0, 0)]), 4);
        // shouldEqual ((a3 [TinyVector <int, 4> (4, 1, 0, 0)]), 14);
        // shouldEqual ((a3 [TinyVector <int, 4> (4, 0, 1, 0)]), 104);
    // }

    // subarray and diagonal tests
    void test_subarray ()
    {
        {
            Shape3 start, stop(ref.shape());  // whole array
            MultiArrayView <3, T, ChunkedArrayTag> v(stop-start);
            array->viewSubarray(start, v);
            
            MultiArrayView <3, T, ChunkedArrayTag> vt(v.transpose());

            ChunkedSubarrayCopy <3, T> c(stop-start);
            array->copySubarray(start, c);
            
            MultiArrayView <3, T, StridedArrayTag> vr = ref.subarray(start, stop);
            MultiArrayView <3, T, StridedArrayTag> vtr = vr.transpose();
            
            shouldEqual(v.shape(), vr.shape());
            should(v == vr);
            shouldEqualSequence(v.begin(), v.end(), vr.begin());
            shouldEqualIndexing(3, v, vr);
            
            shouldEqual(vt.shape(), vtr.shape());
            should(vt == vtr);
            shouldEqualSequence(vt.begin(), vt.end(), vtr.begin());
            shouldEqualIndexing(3, vt, vtr);
            
            shouldEqual(c.shape(), vr.shape());
            should(c == vr);
            shouldEqualSequence(c.begin(), c.end(), vr.begin());
            shouldEqualIndexing(3, c, vr);
        }
        
        {
            Shape3 start(3,2,1), stop(4,5,6);  // single chunk
            MultiArrayView <3, T, ChunkedArrayTag> v(stop-start);
            array->viewSubarray(start, v);
            
            MultiArrayView <3, T, ChunkedArrayTag> vt(v.transpose());

            ChunkedSubarrayCopy <3, T> c(stop-start);
            array->copySubarray(start, c);
            
            MultiArrayView <3, T, StridedArrayTag> vr = ref.subarray(start, stop);
            MultiArrayView <3, T, StridedArrayTag> vtr = vr.transpose();
            
            shouldEqual(v.shape(), vr.shape());
            should(v == vr);
            shouldEqualSequence(v.begin(), v.end(), vr.begin());
            shouldEqualIndexing(3, v, vr);
            
            shouldEqual(vt.shape(), vtr.shape());
            should(vt == vtr);
            shouldEqualSequence(vt.begin(), vt.end(), vtr.begin());
            shouldEqualIndexing(3, vt, vtr);
            
            shouldEqual(c.shape(), vr.shape());
            should(c == vr);
            shouldEqualSequence(c.begin(), c.end(), vr.begin());
            shouldEqualIndexing(3, c, vr);
        }
        
        {
            Shape3 start(7,6,5), stop(9,10,11); // across chunk borders
            MultiArrayView <3, T, ChunkedArrayTag> v(stop-start);
            array->viewSubarray(start, v);
            
            ChunkedSubarrayCopy <3, T> c(stop-start);
            array->copySubarray(start, c);
            
            MultiArrayView <3, T, StridedArrayTag> vr = ref.subarray(start, stop);
            
            shouldEqual(v.shape(), vr.shape());
            should(v == vr);
            shouldEqualSequence(v.begin(), v.end(), vr.begin());
            shouldEqualIndexing(3, v, vr);
            
            shouldEqual(c.shape(), vr.shape());
            should(c == vr);
            shouldEqualSequence(c.begin(), c.end(), vr.begin());
            shouldEqualIndexing(3, c, vr);
        }
        
        // MultiArray<1, scalar_type> diagRef(Shape1(10));
        // linearSequence(diagRef.begin(), diagRef.end(), 0, 111);

        // MultiArrayView <1, scalar_type, StridedArrayTag> diag = array3.diagonal();
        // shouldEqual(diag.shape(0), 10);
        // shouldEqualSequence(diagRef.begin(), diagRef.end(), diag.begin());

        // typedef difference3_type Shape;
        
        // Shape offset (1,1,1);
        // Shape size (5,5,5);
        // MultiArrayView <3, scalar_type> array = array3.subarray (offset, size);
        // shouldEqual (array [Shape (0,0,0)], 111);
        // shouldEqual (array [Shape (5,2,1)], 236);
        // shouldEqual (array (0,0,0), 111);
        // shouldEqual (array (5,2,1), 236);

        // shouldEqual(array.shape(), array3.subarray(offset, Shape3(-5)).shape());
        // shouldEqualSequence(array.begin(), array.end(), array3.subarray(offset, Shape3(-5)).begin());

        // shouldEqual(array.shape(), array3.subarray(Shape3(-9), size).shape());
        // shouldEqualSequence(array.begin(), array.end(), array3.subarray(offset, size).begin());

        // shouldEqual(array.shape(), array3.subarray(Shape3(-9), Shape3(-5)).shape());
        // shouldEqualSequence(array.begin(), array.end(), array3.subarray(Shape3(-9), Shape3(-5)).begin());
        
        // // test swap
        // array3.subarray(Shape(0,0,0), Shape(10,10,1)).swapData( 
             // array3.subarray(Shape(0,0,1), Shape(10,10,2)));
        // for(int k=0; k<100; ++k)
            // shouldEqual(array3[k], k+100);
        // for(int k=100; k<200; ++k)
            // shouldEqual(array3[k], k-100);

    }
        
    // // stridearray tests
    // void test_stridearray ()
    // {
        // difference3_type stride (2, 2, 2);
        // difference3_type traverser (0, 0, 0);
        // MultiArrayView <3, scalar_type, StridedArrayTag>
            // array = array3.stridearray (stride);
        // shouldEqual (array [traverser], 0);
        // traverser [0] += 1;
        // shouldEqual (array [traverser], 2);
        // traverser [1] += 1;
        // shouldEqual (array [traverser], 22);
        // traverser [2] += 1;
        // shouldEqual (array [traverser], 222);
        // shouldEqual (array (1,1,1), 222);
    // }

    // void testIsUnstrided()
    // {
        // typedef difference3_type Shape;

        // should(array3.isUnstrided());
        // should(array3.isUnstrided(0));
        // should(array3.isUnstrided(1));
        // should(array3.isUnstrided(2));
        // should(array3.bindOuter(0).isUnstrided());
        // should(!array3.bindInner(0).isUnstrided());
        // should(!array3.bindAt(1, 0).isUnstrided());
        // should(array3.bindAt(1, 0).isUnstrided(0));
        // should(!array3.subarray(Shape(), array3.shape()-Shape(1)).isUnstrided());
        // should(!array3.subarray(Shape(), array3.shape()-Shape(1)).isUnstrided(1));
        // should(array3.subarray(Shape(), array3.shape()-Shape(1)).isUnstrided(0));
        // should(!array3.subarray(Shape(), array3.shape()-Shape(0,2,2)).isUnstrided());
        // should(array3.subarray(Shape(), array3.shape()-Shape(0,2,2)).isUnstrided(1));
        // should(array3.subarray(Shape(), array3.shape()-Shape(0,2,2)).isUnstrided(0));
    // }

    // void testIsStrided()
    // {
        // // for MultiArray<3, Multiband<T> >

        // should(!array3.isUnstrided());
        // should(array3.permuteStridesAscending().isUnstrided());
        // should(!array3.isUnstrided(0));
        // should(!array3.isUnstrided(1));
        // should(!array3.isUnstrided(2));
        // should(array3.bindInner(Shape2(0,0)).isUnstrided());
    // }

    // // permute and transpose tests
    // void testPermute ()
    // {
        // typedef MultiArrayView <3, scalar_type, StridedArrayTag> transposed_view;
        // array3.reshape(difference3_type(3,5,7));
        // for(int k=0; k<array3.size(); ++k)
            // array3[k] = k;

        // array3_type ref(difference3_type(array3.shape(2), array3.shape(0), array3.shape(1)));
        // for(int k=0; k<array3.shape(0); ++k)
            // for(int l=0; l<array3.shape(1); ++l)
                // for(int m=0; m<array3.shape(2); ++m)
                    // ref(m,k,l) = array3(k,l,m);

        // MultiArrayView <3, scalar_type, StridedArrayTag>
                // parray = array3.transpose (difference3_type (2, 0, 1));        
        // shouldEqual(ref.shape(), parray.shape());
        // should(ref == parray);
        
        // if(vigra::detail::ResolveMultiband<T>::value) // array is Multiband<T>
        // {
            // shouldEqual(array3.strideOrdering(), Shape3(1,2,0));
            
            // array3_type ref_ascending(difference3_type(array3.shape(2), array3.shape(0), array3.shape(1)));
            // for(int k=0; k<array3.shape(0); ++k)
                // for(int l=0; l<array3.shape(1); ++l)
                    // for(int m=0; m<array3.shape(2); ++m)
                        // ref_ascending(m,k,l) = array3(k,l,m);

            // transposed_view parray_ascending = array3.permuteStridesAscending();
            // shouldEqual(ref_ascending.shape(), parray_ascending.shape());
            // should(ref_ascending == parray_ascending);
            
            // array3_type ref_descending(difference3_type(array3.shape(1), array3.shape(0), array3.shape(2)));
            // for(int k=0; k<array3.shape(0); ++k)
                // for(int l=0; l<array3.shape(1); ++l)
                    // for(int m=0; m<array3.shape(2); ++m)
                        // ref_descending(l,k,m) = array3(k,l,m);

            // transposed_view parray_descending = array3.permuteStridesDescending();
            // shouldEqual(ref_descending.shape(), parray_descending.shape());
            // should(ref_descending == parray_descending);
        // }
        // else
        // {
            // shouldEqual(array3.strideOrdering(), Shape3(0,1,2));
            // transposed_view parray_ascending = parray.permuteStridesAscending();        
            // shouldEqual(array3.shape(), parray_ascending.shape());
            // should(array3 == parray_ascending);

            // array3_type ref_descending(difference3_type(array3.shape(2), array3.shape(1), array3.shape(0)));
            // for(int k=0; k<array3.shape(0); ++k)
                // for(int l=0; l<array3.shape(1); ++l)
                    // for(int m=0; m<array3.shape(2); ++m)
                        // ref_descending(m,l,k) = array3(k,l,m);

            // MultiArrayView <3, scalar_type, StridedArrayTag>
                    // parray_descending = array3.permuteStridesDescending();        
            // shouldEqual(ref_descending.shape(), parray_descending.shape());
            // should(ref_descending == parray_descending);
        // }

        // array2_type ref2(difference2_type(array3.shape(1), array3.shape(0)));
        // for(int k=0; k<array3.shape(0); ++k)
            // for(int l=0; l<array3.shape(1); ++l)
                // ref2(l, k) = array3(k, l, 0);

        // MultiArrayView <2, scalar_type, StridedArrayTag>
                // array2 = array3.bindOuter(0).transpose ();
        // shouldEqual(ref2.shape(), array2.shape());
        // should(ref2 == array2);

        // try {
            // array3.transpose (difference3_type (2, 0, 0));   
            // failTest("no exception thrown");
        // }
        // catch(vigra::ContractViolation & c)
        // {
            // std::string expected("\nPrecondition violation!\nMultiArrayView::transpose(): every dimension must occur exactly once");
            // std::string message(c.what());
            // should(0 == expected.compare(message.substr(0,expected.size())));
        // }
    // }

    // void testMethods ()
    // {
        // shouldEqual(array3.squaredNorm(), 332833500);
        
        // shouldEqual(array3.norm(), std::sqrt(332833500.0));
        // shouldEqual(array3.norm(0), 999.0);
        // shouldEqual(array3.norm(1), 499500.0);
        // shouldEqualTolerance(array3.norm(2, false), std::sqrt(332833500.0), 1e-14);
        
        // difference3_type first(0,0,0), last(1,1,1);
        // shouldEqual(array3.subarray(first, last).norm(), 0.0);
        // shouldEqual(array3.subarray(first, last).norm(0), 0.0);
        // shouldEqual(array3.subarray(first, last).norm(1), 0.0);
        // shouldEqual(array3.subarray(first, last).norm(2, false), 0.0);

        // shouldEqual(array3.squaredNorm(), squaredNorm(array3));
        // shouldEqual(array3.norm(), vigra::norm(array3));

        // should(array3.any());
        // should(!array3.subarray(first, last).any());
        // should(!array3.all());
        // should(array3.subarray(last, array3.shape()).all());

        // shouldEqual(array3.template sum<int>(), 499500);
        // shouldEqual(array3.subarray(Shape3(1,1,1),Shape3(3,3,2)).template product<int>(), 183521184);

        // Shape3 reducedShape(1, 1, array3.shape(2));
        // array3_type reducedSums(reducedShape);
        // array3.sum(reducedSums);
        // int res = 4950;
        // for(int k=0; k<reducedShape[2]; ++k, res += 10000)
            // shouldEqual(reducedSums(0,0,k), res);

        // scalar_type minimum, maximum;
        // array3.minmax(&minimum, &maximum);
        // shouldEqual(minimum, 0);
        // shouldEqual(maximum, array3.size()-1);

        // double mean, variance;
        // array3.meanVariance(&mean, &variance);
        // shouldEqual(mean, 499.5);
        // shouldEqual(variance, 83333.25);
    // }
    
    // void testScanOrderAccess()
    // {
        // shouldEqual(array3.size(), 1000);
        // for(int k=0; k< array3.size(); ++k)
        // {
            // shouldEqual(array3[k], k);
            // shouldEqual(array3[array3.scanOrderIndexToCoordinate(k)], k);
            // shouldEqual(array3.coordinateToScanOrderIndex(array3.scanOrderIndexToCoordinate(k)), k);
        // }
            
        // MultiArrayView <2, scalar_type, StridedArrayTag> array = array3.bindInner(2);
        // shouldEqual(array.size(), 100);
        // for(int k=0; k< array.size(); ++k)
            // shouldEqual(array[k], 10*k+2);
            
        // MultiArrayView <2, scalar_type, array3_stride>
            // subarray = array3.bindOuter(2).subarray(Shape2(1,0), Shape2(10,9));
        // shouldEqual(subarray.size(), 81);
        // for(int k=0, l=200; k< subarray.size(); ++k, ++l)
        // {
            // if(k%9 == 0)
                // ++l;
            // shouldEqual(subarray[k], l);
        // }
    // }
    
    // void testAssignmentAndReset()
    // {
        // typedef Shape3 Shape;
        // typename array3_type::view_type array;
        // array = array3;
        // should(array3 == array);
        // try {
            // array = array3.subarray(Shape(0,0,0), Shape(10,1,1));
            // failTest("no exception thrown");
        // }
        // catch(vigra::ContractViolation & c)
        // {
            // std::string expected("\nPrecondition violation!\nMultiArrayView::operator=(MultiArrayView const &): shape mismatch");
            // std::string message(c.what());
            // should(0 == expected.compare(message.substr(0,expected.size())));
        // }
        // MultiArrayView <3, scalar_type, array3_stride> subarray = array3.subarray(Shape(0,0,0), Shape(10,1,1));
        // subarray = array3.subarray(Shape(0,1,0), Shape(10,2,1)); // should overwrite the data
        // for(unsigned int k=0; k<10; ++k)
            // shouldEqual(array3(k,0,0), array3(k,1,0));
        // subarray += array3.subarray(Shape(0,1,0), Shape(10,2,1)); // should overwrite the data
        // for(unsigned int k=0; k<10; ++k)
            // shouldEqual(array3(k,0,0), 2.0*array3(k,1,0));
        // subarray -= array3.subarray(Shape(0,1,0), Shape(10,2,1)); // should overwrite the data
        // for(unsigned int k=0; k<10; ++k)
            // shouldEqual(array3(k,0,0), array3(k,1,0));
        // subarray *= array3.subarray(Shape(0,1,0), Shape(10,2,1)); // should overwrite the data
        // for(unsigned int k=0; k<10; ++k)
            // shouldEqual(array3(k,0,0), sq(array3(k,1,0)));
        // subarray /= array3.subarray(Shape(0,1,0), Shape(10,2,1)); // should overwrite the data
        // for(unsigned int k=0; k<10; ++k)
            // shouldEqual(array3(k,0,0), array3(k,1,0));
            
        // subarray += 1; // should overwrite the data
        // for(unsigned int k=0; k<10; ++k)
            // shouldEqual(array3(k,0,0), array3(k,1,0) + 1);
        // subarray -= 2; // should overwrite the data
        // for(unsigned int k=0; k<10; ++k)
            // shouldEqual(array3(k,0,0), array3(k,1,0) - 1);
        // subarray *= 2; // should overwrite the data
        // for(unsigned int k=0; k<10; ++k)
            // shouldEqual(array3(k,0,0), 2*(array3(k,1,0) - 1));
        // subarray /= 2; // should overwrite the data
        // for(unsigned int k=0; k<10; ++k)
            // shouldEqual(array3(k,0,0), array3(k,1,0) - 1);

        // // test assignment from UInt8 => double (reproduces compiler crash in VisualStudio 2010)
        // MultiArray<2, UInt8> d1(Shape2(50, 100));
        // MultiArray<2, double> d2(d1.shape());
        // d2 = d1;
    // }
};

#if 0

class MultiArrayTest
{
public:
    typedef unsigned char scalar_type;
    typedef MultiArray <1, scalar_type> array1_t;
    typedef array1_t::difference_type shape1_t;
    typedef MultiArray <3, scalar_type> array3_t;
    typedef array3_t::difference_type shape3_t;
    typedef array3_t::traverser traverser3_t;
    typedef traverser3_t::next_type traverser2_t;
    typedef traverser2_t::next_type traverser1_t;
    typedef MultiArrayView<3, scalar_type> array_view3_t;
    typedef array_view3_t::iterator iterator3_t;
    
    shape3_t s;
    array3_t a3;
    
    MultiArrayTest()
    : s(shape3_t(2,3,5)),
      a3(s, 1)
    {
    }
    
    void test_default_ctor ()
    {
        array1_t a;
        shouldEqual (a.shape (0), 0);
    }

    void test_first_ctor ()
    {
        using namespace multi_math;
        array1_t a (shape1_t(2));
        should (a.shape (0) == 2);
        should (a.size() == 2);
        should(all(a == 0));

        array1_t b(4);
        should (b.shape (0) == 4);
        should (b.size() == 4);
        should(all(b == 0));

        typedef MultiArray <2, unsigned char> array2_t;
        array2_t a2 (Shape2(2, 4));
        should (a2.shape (0) == 2);
        should (a2.shape (1) == 4);
        should (a2.width() == 2);
        should (a2.height() == 4);
        should (a2.size() == 8);
        should(all(a2 == 0));

        array2_t b2(5,3);
        should (b2.shape (0) == 5);
        should (b2.shape (1) == 3);
        should (b2.width() == 5);
        should (b2.height() == 3);
        should (b2.size() == 15);
        should(all(b2 == 0));
    }

    void test_second_ctor ()
    {
        array1_t a (shape1_t(2), 1);
        shouldEqual (a.shape (0), 2);
        shouldEqual (a(0), 1);
        shouldEqual (a(1), 1);
        
        a.init(2);
        shouldEqual(a(0), 2);
        shouldEqual(a(1), 2);
        should(a == array1_t(shape1_t(2), 4).init(2));
    }

    void test_assignment ()
    {
        array1_t a (shape1_t(2), 1);
        array1_t b;
        b = a;
        shouldEqual (b.shape (0), 2);
    }

    void test_copy_construction ()
    {
        array1_t a (shape1_t(2), 1);
        array1_t b (a);
        shouldEqual (b.shape (0), 2);
    }

    void testShape()
    {
        shouldEqual(a3.shape(0), 2);
        shouldEqual(a3.shape(1), 3);
        shouldEqual(a3.shape(2), 5);
        should(a3.isInside(shape3_t(1,2,3)));
        should(!a3.isInside(shape3_t(1,23,3)));
    }

    void test_iterator ()
    {
        // test scan-order navigation
        array_view3_t av = a3;
        iterator3_t i1 = av.begin();
        iterator3_t i2 = av.begin();
        iterator3_t iend = av.end();
        iterator3_t i3;
        MultiCoordinateIterator<3> c(av.shape()),
                                   cend = c.getEndIterator();

        should(i1.isValid() && !i1.atEnd());
        should(!iend.isValid() && iend.atEnd());
        should(iend.getEndIterator() == iend);
        
        shouldEqual(i1.point(), *c);
        shouldEqual((i1+0).point(), *(c+0));
        shouldEqual((i1+1).point(), *(c+1));
        shouldEqual((i1+2).point(), *(c+2));
        shouldEqual((i1+3).point(), *(c+3));
        shouldEqual((i1+6).point(), *(c+6));
        shouldEqual((i1+7).point(), *(c+7));
        shouldEqual((i1+9).point(), *(c+9));

        shouldEqual((i1+0).point(), c[0]);
        shouldEqual((i1+1).point(), c[1]);
        shouldEqual((i1+2).point(), c[2]);
        shouldEqual((i1+3).point(), c[3]);
        shouldEqual((i1+6).point(), c[6]);
        shouldEqual((i1+7).point(), c[7]);
        shouldEqual((i1+9).point(), c[9]);

        shouldEqual((iend-1).point(), *(cend-1));
        shouldEqual((iend-2).point(), *(cend-2));
        shouldEqual((iend-3).point(), *(cend-3));
        shouldEqual((iend-7).point(), *(cend-7));
        shouldEqual((iend-8).point(), *(cend-8));
        shouldEqual((iend-10).point(), *(cend-10));

        shouldEqual(&i1[0], &a3(0,0,0));
        shouldEqual(&i1[1], &a3(1,0,0));
        shouldEqual(&i1[2], &a3(0,1,0));
        shouldEqual(&i1[3], &a3(1,1,0));
        shouldEqual(&i1[6], &a3(0,0,1));
        shouldEqual(&i1[7], &a3(1,0,1));
        shouldEqual(&i1[9], &a3(1,1,1));

        shouldEqual(&*(i1+0), &a3(0,0,0));
        shouldEqual(&*(i1+1), &a3(1,0,0));
        shouldEqual(&*(i1+2), &a3(0,1,0));
        shouldEqual(&*(i1+3), &a3(1,1,0));
        shouldEqual(&*(i1+6), &a3(0,0,1));
        shouldEqual(&*(i1+7), &a3(1,0,1));
        shouldEqual(&*(i1+9), &a3(1,1,1));

        shouldEqual(&*(i1+shape3_t(0,0,0)), &a3(0,0,0));
        shouldEqual(&*(i1+shape3_t(1,0,0)), &a3(1,0,0));
        shouldEqual(&*(i1+shape3_t(0,1,0)), &a3(0,1,0));
        shouldEqual(&*(i1+shape3_t(1,1,0)), &a3(1,1,0));
        shouldEqual(&*(i1+shape3_t(0,0,1)), &a3(0,0,1));
        shouldEqual(&*(i1+shape3_t(1,0,1)), &a3(1,0,1));
        shouldEqual(&*(i1+shape3_t(1,1,1)), &a3(1,1,1));

        shouldEqual(&*(iend-1), &a3(1,2,4));
        shouldEqual(&*(iend-2), &a3(0,2,4));
        shouldEqual(&*(iend-3), &a3(1,1,4));
        shouldEqual(&*(iend-7), &a3(1,2,3));
        shouldEqual(&*(iend-8), &a3(0,2,3));
        shouldEqual(&*(iend-10), &a3(0,1,3));
        
        i3 = iend;
        --i3;
        should(i3.isValid() && !i3.atEnd());
        should(i3.getEndIterator() == iend);
        shouldEqual(&*iend, &*(i3+1));
        shouldEqual(&*i3, &a3(1,2,4));
        --i3;
        should(i3.isValid() && !i3.atEnd());
        should(i3.getEndIterator() == iend);
        shouldEqual(&*i3, &a3(0,2,4));
        --i3;
        should(i3.isValid() && !i3.atEnd());
        should(i3.getEndIterator() == iend);
        shouldEqual(&*i3, &a3(1,1,4));
        i3 -= 4;
        should(i3.isValid() && !i3.atEnd());
        should(i3.getEndIterator() == iend);
        shouldEqual(&*i3, &a3(1,2,3));
        --i3;
        should(i3.isValid() && !i3.atEnd());
        should(i3.getEndIterator() == iend);
        shouldEqual(&*i3, &a3(0,2,3));
        --i3;
        --i3;
        should(i3.isValid() && !i3.atEnd());
        should(i3.getEndIterator() == iend);
        shouldEqual(&*i3, &a3(0,1,3));

        i3 = iend-1;
        shouldEqual(&*(i3-shape3_t(0,0,0)), &a3(1,2,4));
        shouldEqual(&*(i3-shape3_t(1,0,0)), &a3(0,2,4));
        shouldEqual(&*(i3-shape3_t(0,1,0)), &a3(1,1,4));
        shouldEqual(&*(i3-shape3_t(1,1,0)), &a3(0,1,4));
        shouldEqual(&*(i3-shape3_t(0,0,1)), &a3(1,2,3));
        shouldEqual(&*(i3-shape3_t(1,0,1)), &a3(0,2,3));
        shouldEqual(&*(i3-shape3_t(1,1,1)), &a3(0,1,3));

        shouldEqual(&iend[-1], &a3(1,2,4));
        shouldEqual(&iend[-2], &a3(0,2,4));
        shouldEqual(&iend[-3], &a3(1,1,4));
        shouldEqual(&iend[-7], &a3(1,2,3));
        shouldEqual(&iend[-8], &a3(0,2,3));
        shouldEqual(&iend[-10], &a3(0,1,3));

        i3 = i1;
        i3 += shape3_t(0,0,1);
        shouldEqual(i3.index(), 6);
        shouldEqual(i3.point(), shape3_t(0,0,1));
        i3 -= shape3_t(0,0,1);
        shouldEqual(i3.index(), 0);
        shouldEqual(i3.point(), shape3_t(0,0,0));
        should(i3 == i1);

        unsigned int count = 0;
        shape3_t p;
        i3 = av.begin();
        iterator3_t i4(av.begin()); 
        iterator3_t i5(av.begin()); 
        iterator3_t i6(av.begin()); 

        // iterate over the third dimension
        for (p[2]=0, i3.resetDim(2), i4.setDim(2, 0), i5.dim<2>() = 0, i6.resetDim(2); 
                i3.point(2) != s[2]; 
                i3.incDim(2), i4.addDim(2, 1), ++i5.dim<2>(), i6.dim<2>() += 1, ++p[2]) 
        {
            for (p[1]=0, i3.resetDim(1), i4.setDim(1, 0), i5.dim<1>() = 0, i6.resetDim(1); 
                    i3.point(1) != s[1]; 
                    i3.incDim(1), i4.addDim(1, 1), ++i5.dim<1>(), i6.dim<1>() += 1, ++p[1]) 
            {
                for (p[0]=0, i3.resetDim(0), i4.setDim(0, 0), i5.dim<0>() = 0, i6.resetDim(0); 
                        i3.point(0) != s[0]; 
                        i3.incDim(0), i4.addDim(0, 1), ++i5.dim<0>(), i6.dim<0>() += 1, ++p[0], ++i1, ++c, i2 += 1, ++count)
                {
                    shouldEqual(&*i1, &a3[p]);
                    shouldEqual(&*i2, &a3[p]);
                    shouldEqual(&*i3, &a3[p]);
                    shouldEqual(&*i4, &a3[p]);
                    shouldEqual(&*i5, &a3[p]);
                    shouldEqual(&*i6, &a3[p]);
                    shouldEqual(i1.operator->(), &a3[p]);
                    shouldEqual(i2.operator->(), &a3[p]);
                    shouldEqual(*c, p);

                    shouldEqual(i1.point(), p);
                    shouldEqual(i2.point(), p);
                    shouldEqual(i3.point(), p);
                    shouldEqual(i4.point(), p);
                    shouldEqual(i5.point(), p);
                    shouldEqual(i6.point(), p);

                    shouldEqual(i1.index(), count);
                    shouldEqual(i2.index(), count);
                    shouldEqual(i3.index(), count);
                    shouldEqual(i4.index(), count);
                    shouldEqual(i5.index(), count);
                    shouldEqual(i6.index(), count);

                    should(i1 != iend);
                    should(!(i1 == iend));
                    should(i1 < iend);
                    should(i1 <= iend);
                    should(!(i1 > iend));
                    should(!(i1 >= iend));

                    should(i5.dim<2>() == p[2]);
                    should(i5.dim<1>() == p[1]);
                    should(i5.dim<0>() == p[0]);
                    should(i5.dim<2>() != s[2]);
                    should(i5.dim<1>() != s[1]);
                    should(i5.dim<0>() != s[0]);
                    should(i5.dim<2>() <= p[2]);
                    should(i5.dim<1>() <= p[1]);
                    should(i5.dim<0>() <= p[0]);
                    should(i5.dim<2>() < s[2]);
                    should(i5.dim<1>() < s[1]);
                    should(i5.dim<0>() < s[0]);
                    should(i5.dim<2>() >= 0);
                    should(i5.dim<1>() >= 0);
                    should(i5.dim<0>() >= 0);
                    shouldNot(i5.dim<2>() > s[2]);
                    shouldNot(i5.dim<1>() > s[1]);
                    shouldNot(i5.dim<0>() > s[0]);

                    shouldEqual(iend - i1, av.size() - count);

                    bool atBorder = p[0] == 0 || p[0] == s[0]-1 || p[1] == 0 || p[1] == s[1]-1 ||
                                    p[2] == 0 || p[2] == s[2]-1;
                    if(!atBorder)
                    {
                        should(!i1.atBorder());
                        should(!i2.atBorder());
                    }
                    else
                    {
                        should(i1.atBorder());
                        should(i2.atBorder());
                    }
                }
            }
        }

        should(c == cend);
        should(i1 == iend);
        should(!(i1 != iend));
        should(!(i1 < iend));
        should(i1 <= iend);
        should(!(i1 > iend));
        should(i1 >= iend);

        should(i2 == iend);
        should(!(i2 != iend));
        should(!(i2 < iend));
        should(i2 <= iend);
        should(!(i2 > iend));
        should(i2 >= iend);

        shouldEqual(iend - i1, 0);
        shouldEqual(iend - i2, 0);
        shouldEqual (count, av.size());

        --i1;
        i2 -= 1;
        shouldEqual(&*i1, &a3(1,2,4));
        shouldEqual(&*i2, &a3(1,2,4));
    }

    void test_const_iterator()
    {
        typedef array_view3_t::const_iterator iterator;

        // test const scan-order navigation
        array_view3_t av = a3;
        iterator i1 = const_cast<array_view3_t const &>(av).begin();
        iterator i2 = const_cast<array_view3_t const &>(av).begin();
        iterator iend = const_cast<array_view3_t const &>(av).end();

        shouldEqual(&i1[0], &a3(0,0,0));
        shouldEqual(&i1[1], &a3(1,0,0));
        shouldEqual(&i1[2], &a3(0,1,0));
        shouldEqual(&i1[3], &a3(1,1,0));
        shouldEqual(&i1[6], &a3(0,0,1));
        shouldEqual(&i1[7], &a3(1,0,1));
        shouldEqual(&i1[9], &a3(1,1,1));

        shouldEqual(&*(i1+0), &a3(0,0,0));
        shouldEqual(&*(i1+1), &a3(1,0,0));
        shouldEqual(&*(i1+2), &a3(0,1,0));
        shouldEqual(&*(i1+3), &a3(1,1,0));
        shouldEqual(&*(i1+6), &a3(0,0,1));
        shouldEqual(&*(i1+7), &a3(1,0,1));
        shouldEqual(&*(i1+9), &a3(1,1,1));

        shouldEqual(&*(iend-1), &a3(1,2,4));
        shouldEqual(&*(iend-2), &a3(0,2,4));
        shouldEqual(&*(iend-3), &a3(1,1,4));
        shouldEqual(&*(iend-7), &a3(1,2,3));
        shouldEqual(&*(iend-8), &a3(0,2,3));
        shouldEqual(&*(iend-10), &a3(0,1,3));

        shouldEqual(&iend[-1], &a3(1,2,4));
        shouldEqual(&iend[-2], &a3(0,2,4));
        shouldEqual(&iend[-3], &a3(1,1,4));
        shouldEqual(&iend[-7], &a3(1,2,3));
        shouldEqual(&iend[-8], &a3(0,2,3));
        shouldEqual(&iend[-10], &a3(0,1,3));

        unsigned int count = 0;
        shape3_t p;

        // iterate over the third dimension
        for (p[2]=0; p[2] != s[2]; ++p[2]) 
        {
            for (p[1]=0; p[1] != s[1]; ++p[1]) 
            {
                for (p[0]=0; p[0] != s[0]; ++p[0], ++i1, i2 += 1, ++count)
                {
                    shouldEqual(&*i1, &a3[p]);
                    shouldEqual(&*i2, &a3[p]);
                    shouldEqual(i1.operator->(), &a3[p]);
                    shouldEqual(i2.operator->(), &a3[p]);
                    shouldEqual(i1.point(), p);
                    shouldEqual(i2.point(), p);
                    shouldEqual(i1.index(), count);
                    shouldEqual(i2.index(), count);

                    should(i1 != iend);
                    should(!(i1 == iend));
                    should(i1 < iend);
                    should(i1 <= iend);
                    should(!(i1 > iend));
                    should(!(i1 >= iend));

                    shouldEqual(iend - i1, av.size() - count);

                    bool atBorder = p[0] == 0 || p[0] == s[0]-1 || p[1] == 0 || p[1] == s[1]-1 ||
                                    p[2] == 0 || p[2] == s[2]-1;
                    if(!atBorder)
                    {
                        should(!i1.atBorder());
                        should(!i2.atBorder());
                    }
                    else
                    {
                        should(i1.atBorder());
                        should(i2.atBorder());
                    }
                }
            }
        }

        should(i1 == iend);
        should(!(i1 != iend));
        should(!(i1 < iend));
        should(i1 <= iend);
        should(!(i1 > iend));
        should(i1 >= iend);

        should(i2 == iend);
        should(!(i2 != iend));
        should(!(i2 < iend));
        should(i2 <= iend);
        should(!(i2 > iend));
        should(i2 >= iend);

        shouldEqual(iend - i1, 0);
        shouldEqual(iend - i2, 0);
        shouldEqual (count, av.size());

        --i1;
        i2 -= 1;
        shouldEqual(&*i1, &a3(1,2,4));
        shouldEqual(&*i2, &a3(1,2,4));
    }

    void test_coupled_iterator ()
    {
        // test scan-order navigation
        typedef CoupledIteratorType<3, scalar_type>::type Iterator;
        Iterator i1 = createCoupledIterator(a3);
        Iterator iend = i1.getEndIterator();
        Iterator i2 = i1;
        Iterator i3;

        shouldEqual(&(*(*i1)), &a3(0,0,0));
        shouldEqual(&get<1>(i1[0]), &a3(0,0,0));
        shouldEqual(&get<1>(i1[1]), &a3(1,0,0));
        shouldEqual(&get<1>(i1[2]), &a3(0,1,0));
        shouldEqual(&get<1>(i1[3]), &a3(1,1,0));
        shouldEqual(&get<1>(i1[6]), &a3(0,0,1));
        shouldEqual(&get<1>(i1[7]), &a3(1,0,1));
        shouldEqual(&get<1>(i1[9]), &a3(1,1,1));

        shouldEqual(&get<1>(*(i1+0)), &a3(0,0,0));
        shouldEqual(&get<1>(*(i1+1)), &a3(1,0,0));
        shouldEqual(&get<1>(*(i1+2)), &a3(0,1,0));
        shouldEqual(&get<1>(*(i1+3)), &a3(1,1,0));
        shouldEqual(&get<1>(*(i1+6)), &a3(0,0,1));
        shouldEqual(&get<1>(*(i1+7)), &a3(1,0,1));
        shouldEqual(&get<1>(*(i1+9)), &a3(1,1,1));

        shouldEqual(&get<1>(*(i1+shape3_t(0,0,0))), &a3(0,0,0));
        shouldEqual(&get<1>(*(i1+shape3_t(1,0,0))), &a3(1,0,0));
        shouldEqual(&get<1>(*(i1+shape3_t(0,1,0))), &a3(0,1,0));
        shouldEqual(&get<1>(*(i1+shape3_t(1,1,0))), &a3(1,1,0));
        shouldEqual(&get<1>(*(i1+shape3_t(0,0,1))), &a3(0,0,1));
        shouldEqual(&get<1>(*(i1+shape3_t(1,0,1))), &a3(1,0,1));
        shouldEqual(&get<1>(*(i1+shape3_t(1,1,1))), &a3(1,1,1));

        shouldEqual(&get<1>(*(iend-1)), &a3(1,2,4));
        shouldEqual(&get<1>(*(iend-2)), &a3(0,2,4));
        shouldEqual(&get<1>(*(iend-3)), &a3(1,1,4));
        shouldEqual(&get<1>(*(iend-7)), &a3(1,2,3));
        shouldEqual(&get<1>(*(iend-8)), &a3(0,2,3));
        shouldEqual(&get<1>(*(iend-10)), &a3(0,1,3));

        i3 = iend-1;
        shouldEqual(&get<1>(*(i3-shape3_t(0,0,0))), &a3(1,2,4));
        shouldEqual(&get<1>(*(i3-shape3_t(1,0,0))), &a3(0,2,4));
        shouldEqual(&get<1>(*(i3-shape3_t(0,1,0))), &a3(1,1,4));
        shouldEqual(&get<1>(*(i3-shape3_t(1,1,0))), &a3(0,1,4));
        shouldEqual(&get<1>(*(i3-shape3_t(0,0,1))), &a3(1,2,3));
        shouldEqual(&get<1>(*(i3-shape3_t(1,0,1))), &a3(0,2,3));
        shouldEqual(&get<1>(*(i3-shape3_t(1,1,1))), &a3(0,1,3));

        shouldEqual(&get<1>(iend[-1]), &a3(1,2,4));
        shouldEqual(&get<1>(iend[-2]), &a3(0,2,4));
        shouldEqual(&get<1>(iend[-3]), &a3(1,1,4));
        shouldEqual(&get<1>(iend[-7]), &a3(1,2,3));
        shouldEqual(&get<1>(iend[-8]), &a3(0,2,3));
        shouldEqual(&get<1>(iend[-10]), &a3(0,1,3));

        i3 = i1;
        i3 += shape3_t(0,0,1);
        shouldEqual(i3.scanOrderIndex(), 6);
        shouldEqual(i3.point(), shape3_t(0,0,1));
        i3 -= shape3_t(0,0,1);
        shouldEqual(i3.scanOrderIndex(), 0);
        shouldEqual(i3.point(), shape3_t(0,0,0));
        should(i3 == i1);

        unsigned int count = 0;
        shape3_t p;

        // iterate over the third dimension
        for (p[2]=0; p[2] != s[2]; ++p[2]) 
        {
            for (p[1]=0; p[1] != s[1]; ++p[1]) 
            {
                for (p[0]=0; p[0] != s[0]; ++p[0], ++i1, i2 += 1, ++count)
                {
                    shouldEqual(&get<1>(*i1), &a3[p]);
                    shouldEqual(&get<1>(*i2), &a3[p]);
                    shouldEqual(&i1.get<1>(), &a3[p]);
                    shouldEqual(&i2.get<1>(), &a3[p]);
                    //shouldEqual(i1.operator->(), &a3[p]);
                    //shouldEqual(i2.operator->(), &a3[p]);
                    shouldEqual(i1.point(), p);
                    shouldEqual(i2.point(), p);
                    shouldEqual(i1.get<0>(), p);
                    shouldEqual(i2.get<0>(), p);
                    shouldEqual(i1.scanOrderIndex(), count);
                    shouldEqual(i2.scanOrderIndex(), count);

                    should(i1 != iend);
                    should(!(i1 == iend));
                    should(i1 < iend);
                    should(i1 <= iend);
                    should(!(i1 > iend));
                    should(!(i1 >= iend));

                    shouldEqual(iend - i1, a3.size() - count);

                    bool atBorder = p[0] == 0 || p[0] == s[0]-1 || p[1] == 0 || p[1] == s[1]-1 ||
                                    p[2] == 0 || p[2] == s[2]-1;
                    if(!atBorder)
                    {
                        should(!i1.atBorder());
                        should(!i2.atBorder());
                    }
                    else
                    {
                        should(i1.atBorder());
                        should(i2.atBorder());
                    }
                }
            }
        }

        should(i1 == iend);
        should(!(i1 != iend));
        should(!(i1 < iend));
        should(i1 <= iend);
        should(!(i1 > iend));
        should(i1 >= iend);

        should(i2 == iend);
        should(!(i2 != iend));
        should(!(i2 < iend));
        should(i2 <= iend);
        should(!(i2 > iend));
        should(i2 >= iend);

        shouldEqual(iend - i1, 0);
        shouldEqual(iend - i2, 0);
        shouldEqual (count, a3.size());

        --i1;
        i2 -= 1;
        shouldEqual(&get<1>(*i1), &a3(1,2,4));
        shouldEqual(&get<1>(*i2), &a3(1,2,4));

        int idata[] = { 10, 11, 12, 13 };
        double ddata[] = { 20.0, 21.0, 22.0, 23.0 };

        typedef CoupledIteratorType<1>::type Iterator0;
        typedef CoupledIteratorType<1, int, double>::type Iterator1;
        MultiArrayView<1, int> vi(Shape1(4), idata);
        MultiArrayView<1, double> vd(Shape1(4), ddata);

        Iterator0 i0 = createCoupledIterator(Shape1(4));
        Iterator1 it = createCoupledIterator(vi, vd),
                  end = it.getEndIterator();

        count = 0;
        for(; it < end; ++it, ++i0, ++count)
        {
            shouldEqual(i0.get<0>(), Shape1(count));
            shouldEqual(it.get<0>(), Shape1(count));
            shouldEqual(it.get<1>(), count+10);
            shouldEqual(it.get<2>(), count+20.0);
        }
        shouldEqual(count, 4);

        // test multiband
        MultiArrayView<3, scalar_type, StridedArrayTag> at = a3.transpose();

        typedef CoupledIteratorType<3, Multiband<scalar_type>, scalar_type, scalar_type>::type MultibandIterator;
        MultibandIterator im = createCoupledIterator(MultiArrayView<3, Multiband<scalar_type>, StridedArrayTag>(at),
                                                     at.bindOuter(0), at.bindOuter(1));
        MultibandIterator imend = im.getEndIterator();
        count = 0;
        for(; im < imend; ++im, ++count)
        {
            shouldEqual(im.get<1>().shape(), Shape1(2));
            shouldEqual(&(im.get<1>()[0]), &(im.get<2>()));
            shouldEqual(&(im.get<1>()[1]), &(im.get<3>()));
        }
        shouldEqual(count, 15);
    }

    void test_traverser ()
    {
        // test hierarchical navigation and 
        traverser3_t i3_f = a3.traverser_begin ();
        traverser3_t i3_l = a3.traverser_end ();
        array3_t::iterator seqi = a3.begin();

        int countx = 0, county = 0, countz = 0;

        // iterate over the third dimension
        for (int z=0; i3_f != i3_l; ++i3_f, ++z) 
        {
            traverser2_t i2_f = i3_f.begin ();
            traverser2_t i2_l = i3_f.end ();
            // iterate over the second dimension
            for (int y=0; i2_f != i2_l; ++i2_f, ++y) 
            {
                traverser1_t i1_f = i2_f.begin ();
                traverser1_t i1_l = i2_f.end ();
                // iterate over the first dimension
                for (int x=0; i1_f != i1_l; ++i1_f, ++x, ++seqi)
                {
                    ++countx;
                    shouldEqual(&*i1_f, &a3(x,y,z));
                    shouldEqual(&*seqi, &a3(x,y,z));
                }
                ++county;
            }
            ++countz;
        }

        shouldEqual (countx, 30);
        shouldEqual (county, 15);
        shouldEqual (countz, 5);
        shouldEqual (seqi, a3.end());
        
        // test direct navigation
        traverser3_t i3 = a3.traverser_begin();
        shouldEqual(&*i3, &a3[shape3_t(0,0,0)]);

        i3.dim<2>()++;
        i3.dim<1>()++;
        i3.dim<0>()++;        
        shouldEqual(&*i3, &a3[shape3_t(1,1,1)]);

        i3.dim<2>()+= 3;
        i3.dim<1>()+= 2;
        i3.dim<0>()+= 1;        
        shouldEqual(&*i3, &a3[shape3_t(2,3,4)]);
        shouldEqual(&i3[shape3_t(-2,-3,-4)], &a3[shape3_t(0,0,0)]);
        shouldEqual(&*(i3-shape3_t(2,3,4)), &a3[shape3_t(0,0,0)]);

        i3.dim<2>()--;
        i3.dim<1>()--;
        i3.dim<0>()--;        
        shouldEqual(&*i3, &a3[shape3_t(1,2,3)]);

        i3.dim<2>()-= 3;
        i3.dim<1>()-= 2;
        i3.dim<0>()-= 1;        
        shouldEqual(&*i3, &a3[shape3_t(0,0,0)]);

        shouldEqual(&i3[shape3_t(2,3,4)], &a3[shape3_t(2,3,4)]);
        shouldEqual(&*(i3+shape3_t(2,3,4)), &i3[shape3_t(2,3,4)]);
    }

    void test_const_traverser ()
    {
        typedef array3_t::const_traverser traverser;

        // test hierarchical navigation
        traverser i3_f = const_cast<array3_t const &>(a3).traverser_begin ();
        traverser i3_l = const_cast<array3_t const &>(a3).traverser_end ();

        unsigned int countx = 0, county = 0, countz = 0;

        // iterate over the third dimension
        for (; i3_f != i3_l; ++i3_f) {
            traverser::next_type i2_f = i3_f.begin ();
            traverser::next_type i2_l = i3_f.end ();
            // iterate over the second dimension
            for (; i2_f != i2_l; ++i2_f) {
                traverser::next_type::next_type i1_f = i2_f.begin ();
                traverser::next_type::next_type i1_l = i2_f.end ();
                // iterate over the first dimension
                for (; i1_f != i1_l; ++i1_f)
                    ++countx;
                ++county;
            }
            ++countz;
        }

        shouldEqual (countx, 30u);
        shouldEqual (county, 15u);
        shouldEqual (countz, 5u);
        
        // test direct navigation
        traverser i3 = const_cast<array3_t const &>(a3).traverser_begin();
        shouldEqual(&*i3, &a3[shape3_t(0,0,0)]);

        i3.dim<2>()++;
        i3.dim<1>()++;
        i3.dim<0>()++;        
        shouldEqual(&*i3, &a3[shape3_t(1,1,1)]);

        i3.dim<2>()+= 3;
        i3.dim<1>()+= 2;
        i3.dim<0>()+= 1;        
        shouldEqual(&*i3, &a3[shape3_t(2,3,4)]);
        shouldEqual(&i3[shape3_t(-2,-3,-4)], &a3[shape3_t(0,0,0)]);
        shouldEqual(&*(i3-shape3_t(2,3,4)), &a3[shape3_t(0,0,0)]);

        i3.dim<2>()--;
        i3.dim<1>()--;
        i3.dim<0>()--;        
        shouldEqual(&*i3, &a3[shape3_t(1,2,3)]);

        i3.dim<2>()-= 3;
        i3.dim<1>()-= 2;
        i3.dim<0>()-= 1;        
        shouldEqual(&*i3, &a3[shape3_t(0,0,0)]);

        shouldEqual(&i3[shape3_t(2,3,4)], &a3[shape3_t(2,3,4)]);
        shouldEqual(&*(i3+shape3_t(2,3,4)), &i3[shape3_t(2,3,4)]);
    }

    void test_bindOuter ()
    {
        MultiArrayView <2, unsigned char> ba = a3.bindOuter(TinyVector<int, 1>(2));

        shouldEqual (ba.shape (0), 2);
        shouldEqual (ba.shape (1), 3);
    }


    void test_bindInner ()
    {
        MultiArrayView <2, unsigned char, StridedArrayTag>
            fa = a3.bindInner(TinyVector <int, 1>(1));

        shouldEqual (fa.shape (0), 3);
        shouldEqual (fa.shape (1), 5);
    }

    void test_bindInnerAll ()
    {
        MultiArrayView <0, unsigned char, StridedArrayTag>
            fa = a3.bindInner(shape3_t(1,1,1));

        shouldEqual (fa.shape (0), 1);
        shouldEqual (fa[shape1_t(shape1_t::value_type(0))], 1.0);
    }

    void test_bindAt ()
    {
        MultiArrayView <2, unsigned char, StridedArrayTag>
            fa = a3.bindAt (1, 1);

        shouldEqual (fa.shape (0), 2);
        shouldEqual (fa.shape (1), 5);
    }

    void test_bind ()
    {
        MultiArrayView <2, unsigned char, StridedArrayTag>
            fa = a3.bind <0> (1);

        shouldEqual (fa.shape (0), 3);
        shouldEqual (fa.shape (1), 5);
    }

    void test_reshape ()
    {
        shouldEqual(a3(0,0), 1);
        a3.reshape(a3.shape());
        shouldEqual(a3(0,0), 0);
    
        a3.reshape (shape3_t(20,30,50));

        shouldEqual (a3.shape (0), 20);
        shouldEqual (a3.shape (1), 30);
        shouldEqual (a3.shape (2), 50);
        shouldEqual(a3(0,0), 0);
    }
    
    void test_copy_int_float()
    {
        MultiArray<2, float> a(MultiArrayShape<2>::type(2,2));
        a.init(3);
        MultiArray<2, unsigned short> b(a.shape());
        b = a;
    }

    void test_subarray ()
    {
        a3.reshape (shape3_t(20,30,50));

        MultiArrayView <3, unsigned char> st = 
            a3.subarray (shape3_t(1,3,3),shape3_t(5,3,5));

        shouldEqual (st.shape (0), 4);
        shouldEqual (st.shape (1), 0);
        shouldEqual (st.shape (2), 2);
    }

    void test_stridearray ()
    {
        typedef MultiArrayView <3, unsigned char, StridedArrayTag>
                strided_array_t;
        
        a3.reshape (shape3_t(20,30,50));
        
        strided_array_t st = a3.stridearray (shape3_t(3,4,6));

        shouldEqual (st.shape (0), 6);
        shouldEqual (st.shape (1), 7);
        shouldEqual (st.shape (2), 8);
        
        // test hierarchical navigation
        typedef strided_array_t::traverser Traverser3;
        Traverser3 i3_f = st.traverser_begin ();
        Traverser3 i3_l = st.traverser_end ();

        unsigned int countx = 0, county = 0, countz = 0;

        // iterate over the third dimension
        for (; i3_f != i3_l; ++i3_f) {
            typedef Traverser3::next_type Traverser2;
            Traverser2 i2_f = i3_f.begin ();
            Traverser2 i2_l = i3_f.end ();
            // iterate over the second dimension
            for (; i2_f != i2_l; ++i2_f) {
                typedef Traverser2::next_type Traverser1;
                Traverser1 i1_f = i2_f.begin ();
                Traverser1 i1_l = i2_f.end ();
                // iterate over the first dimension
                for (; i1_f != i1_l; ++i1_f)
                    ++countx;
                ++county;
            }
            ++countz;
        }

        shouldEqual (countx, 336u);
        shouldEqual (county, 56u);
        shouldEqual (countz, 8u);
  
        // test direct navigation
        strided_array_t::traverser i = st.traverser_begin();        
        shouldEqual(&*i, &st[shape3_t(0,0,0)]);
        shouldEqual(&*i, &a3[shape3_t(0,0,0)]);
        
        i.dim<2>()++;
        i.dim<1>()++;
        i.dim<0>()++;
        shouldEqual(&*i, &st[shape3_t(1,1,1)]);
        shouldEqual(&*i, &a3[shape3_t(3,4,6)]);
        
        i.dim<2>()+= 3;
        i.dim<1>()+= 2;
        i.dim<0>()+= 1;
        shouldEqual(&*i, &st[shape3_t(2,3,4)]);
        shouldEqual(&*i, &a3[shape3_t(6,12,24)]);
        shouldEqual(&i[shape3_t(-2,-3,-4)], &a3[shape3_t(0,0,0)]);
        shouldEqual(&*(i-shape3_t(2,3,4)), &a3[shape3_t(0,0,0)]);
        
        i.dim<2>()--;
        i.dim<1>()--;
        i.dim<0>()--;
        shouldEqual(&*i, &st[shape3_t(1,2,3)]);
        shouldEqual(&*i, &a3[shape3_t(3,8,18)]);
        
        i.dim<2>()-= 3;
        i.dim<1>()-= 2;
        i.dim<0>()-= 1;
        shouldEqual(&*i, &st[shape3_t(0,0,0)]);
        shouldEqual(&*i, &a3[shape3_t(0,0,0)]);
        shouldEqual(&i[shape3_t(2,3,4)], &a3[shape3_t(6,12,24)]);
        shouldEqual(&*(i+shape3_t(2,3,4)), &a3[shape3_t(6,12,24)]);
    }

    void test_expandElements()
    {
        using namespace multi_math;

        MultiArray<3, TinyVector<int, 3> > a(Shape3(4,3,2));
        a.init(TinyVector<int, 3>(1,2,3));

        MultiArrayView<4, int, StridedArrayTag> ex = a.expandElements(0);
        MultiArrayView<4, int, StridedArrayTag>::iterator i = ex.begin();
        while(i != ex.end())
        {
            shouldEqual(*i, 1); ++i;
            shouldEqual(*i, 2); ++i;
            shouldEqual(*i, 3); ++i;
        }

        MultiArrayView<4, int, StridedArrayTag> ex2 = a.expandElements(3);
        i = ex2.begin();
        for(int k=0; k < a.size(); ++i, ++k)
            shouldEqual(*i, 1);
        for(int k=0; k < a.size(); ++i, ++k)
            shouldEqual(*i, 2);
        for(int k=0; k < a.size(); ++i, ++k)
            shouldEqual(*i, 3);

        MultiArray<3, bool> b = (a.bindElementChannel(0) == 1);
        should(b.all());
        b = (a.bindElementChannel(1) == 2);
        should(b.all());
        b = (a.bindElementChannel(2) == 3);
        should(b.all());
    }
};

struct MultiArrayPointoperatorsTest
{

    typedef float PixelType;
    typedef MultiArray<3,PixelType> Image3D;
    typedef MultiArrayView<3,PixelType> View3D;
    typedef Image3D::difference_type Size3;
    typedef MultiArray<1,PixelType> Image1D;
    typedef Image1D::difference_type Size1;

    Image3D img;

    MultiArrayPointoperatorsTest()
    : img(Size3(5,4,3))
    {
        int i;
        PixelType c = 0.1f;
        for(i=0; i<img.elementCount(); ++i, ++c)
            img.data()[i] = c;
    }

    void testInit()
    {
        Image3D res(img.shape());
        const Image3D::value_type ini = 1.1f;
        should(res.shape() == Size3(5,4,3));

        initMultiArray(destMultiArrayRange(res), ini);

        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), ini);

        using namespace multi_math;
        should(all(res == ini));

        initMultiArray(res, 2.2f);
        should(all(res == 2.2f));

        res = 3.3f;
        should(all(res == 3.3f));

        res.init(4.4f);
        should(all(res == 4.4f));
    }

    void testCopy()
    {
        Image3D res(img.shape(), 1.0), res1(img.shape(), 1.0);
        
        copyMultiArray(srcMultiArrayRange(img), destMultiArray(res));
        copyMultiArray(img, res1);

        should(img == res);
        should(img == res1);
    }

    void testCopyOuterExpansion()
    {
        Image3D res(img.shape());

        copyMultiArray(img.subarray(Size3(0,0,0), Size3(5,1,1)), res);
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(x,0,0));
    }

    void testCopyInnerExpansion()
    {
        Image3D res(img.shape());

        copyMultiArray(img.subarray(Size3(0,0,0), Size3(1,1,3)), res);
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(0,0,z));
    }

    void testTransform()
    {
        Image3D res(img.shape()), res1(img.shape());
        transformMultiArray(srcMultiArrayRange(img), destMultiArray(res),
                            Arg1() + Arg1());
        transformMultiArray(img, res1, Arg1() + Arg1());
        
        using namespace multi_math;
        should(all(2.0*img == res));
        should(all(2.0*img == res1));
    }

    void testTransformOuterExpand()
    {
        Image3D res(img.shape());
        transformMultiArray(img.subarray(Size3(0,0,0), Size3(5,1,1)), res,
                            Arg1() + Arg1());
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,0,0));
    }

    void testTransformInnerExpand()
    {
        Image3D res(img.shape());

        transformMultiArray(img.subarray(Size3(0,0,0), Size3(1,1,3)), res,
                            Arg1() + Arg1());
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(0,0,z));
    }

    void testTransformOuterReduce()
    {
        Image3D res(Size3(5,1,1));

        transformMultiArray(img, res, reduceFunctor(Arg1() + Arg2(), 0.0));
        
        int x,y,z;
        for(x=0; x<img.shape(0); ++x)
        {
            double sum = 0.0;
            for(y=0; y<img.shape(1); ++y)
                for(z=0; z<img.shape(2); ++z)
                    sum += img(x,y,z);
            shouldEqual(res(x,0,0), sum);
        }
        
        Image1D res1(Size1(5));
        MultiArrayView<3,PixelType> res3 = res1.insertSingletonDimension(1).insertSingletonDimension(2);
        transformMultiArray(img, res3, FindSum<PixelType>());
        shouldEqualSequenceTolerance(res1.data(), res1.data()+5, res.data(), 1e-7);       
    }

    void testTransformInnerReduce()
    {
        Image3D res(Size3(1,1,3));
        
        transformMultiArray(img, res, reduceFunctor(Arg1() + Arg2(), 0.0));
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
        {
            double sum = 0.0;
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    sum += img(x,y,z);
            shouldEqual(res(0,0,z), sum);
        }
        
        Image1D res1(Size1(3));
        MultiArrayView<3,PixelType> res3 = res1.insertSingletonDimension(0).insertSingletonDimension(0);
        transformMultiArray(img, res3, FindSum<PixelType>());
        shouldEqualSequenceTolerance(res1.data(), res1.data()+3, res.data(), 1e-6);       
    }

    void testCombine2()
    {
        Image3D res(img.shape()), res1(img.shape());
        
        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArray(img), 
                              destMultiArray(res),
                              Arg1() + Arg2());
        combineTwoMultiArrays(img, img, res1, Arg1() + Arg2());
        
        using namespace multi_math;
        should(all(2.0*img == res));
        should(all(2.0*img == res1));
    }

    void testCombine2OuterExpand()
    {
        Image3D res(img.shape());
        
        combineTwoMultiArrays(img.subarray(Size3(0,0,0), Size3(5,1,1)), img, res,
                              Arg1() + Param(2.0)*Arg2());       
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,y,z) + img(x,0,0));

        combineTwoMultiArrays(img, img.subarray(Size3(0,0,0), Size3(5,1,1)), res,
                              Arg1() + Param(2.0)*Arg2());       
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(x,y,z) + 2.0*img(x,0,0));

        View3D view = img.subarray(Size3(0,0,0), Size3(5,1,1));
        combineTwoMultiArrays(srcMultiArrayRange(view), srcMultiArrayRange(view), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 3.0*img(x,0,0));
    }

    void testCombine2InnerExpand()
    {
        Image3D res(img.shape());
        
        View3D view = img.subarray(Size3(0,0,0), Size3(1,1,3));
        combineTwoMultiArrays(view, img, res,
                              Arg1() + Param(2.0)*Arg2());       
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,y,z) + img(0,0,z));

        combineTwoMultiArrays(img, view, res,
                              Arg1() + Param(2.0)*Arg2());       
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(x,y,z) + 2.0*img(0,0,z));

        combineTwoMultiArrays(srcMultiArrayRange(view), srcMultiArrayRange(view), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 3.0*img(0,0,z));
    }

    void testCombine2OuterReduce()
    {
        Image3D res(Size3(5,1,1));
        
        combineTwoMultiArrays(img, img, res,
                              reduceFunctor(Arg1() + Arg2() + Arg3(), 0.0));
        
        int x,y,z;
        for(x=0; x<img.shape(0); ++x)
        {
            double sum = 0.0;
            for(y=0; y<img.shape(1); ++y)
                for(z=0; z<img.shape(2); ++z)
                    sum += img(x,y,z);
            shouldEqual(res(x,0,0), 2.0*sum);
        }
    }

    void testCombine2InnerReduce()
    {
        Image3D res(Size3(1,1,3));
        
        combineTwoMultiArrays(img, img, res,
                              reduceFunctor(Arg1() + Arg2() + Arg3(), 0.0));
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
        {
            double sum = 0.0;
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    sum += img(x,y,z);
            shouldEqual(res(0,0,z), 2.0*sum);
        }
    }

    void testCombine3()
    {
        Image3D res(img.shape()), res1(img.shape());
        
        combineThreeMultiArrays(srcMultiArrayRange(img), 
                                srcMultiArray(img), srcMultiArray(img), 
                                destMultiArray(res),
                                Arg1() + Arg2() + Arg3());
        combineThreeMultiArrays(img, img, img, res1,
                                Arg1() + Arg2() + Arg3());

        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                {
                    shouldEqual(res(x,y,z), 3.0*img(x,y,z));
                    shouldEqual(res1(x,y,z), 3.0*img(x,y,z));
                }
    }
    
    void testInitMultiArrayBorder(){
        typedef vigra::MultiArray<1,int> IntLine;
        typedef vigra::MultiArray<2,int> IntImage;
        typedef vigra::MultiArray<3,int> IntVolume;
        
        const int desired_vol[] ={  0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,

                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,

                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 5, 5, 0, 0,
                                    0, 0, 5, 5, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,

                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 5, 5, 0, 0,
                                    0, 0, 5, 5, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,

                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,

                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0};

        const int desired_img[] ={  0, 0, 0, 0, 0, 0,
                                    0, 5, 5, 5, 5, 0,
                                    0, 5, 5, 5, 5, 0,
                                    0, 5, 5, 5, 5, 0,
                                    0, 5, 5, 5, 5, 0,
                                    0, 0, 0, 0, 0, 0};

        const int desired_lin[] ={  0, 0, 0, 5, 0, 0, 0 };

        const int desired_vol2[] ={  0, 0,
                                     0, 0,

                                     0, 0, 
                                     0, 0};

        IntVolume vol(IntVolume::difference_type(6,6,6));
        
        for(IntVolume::iterator iter=vol.begin(); iter!=vol.end(); ++iter)
            *iter=5;
        initMultiArrayBorder(destMultiArrayRange(vol),2,0);
        shouldEqualSequence(vol.begin(), vol.end(), desired_vol);

        IntImage img(IntImage::difference_type(6,6));
        
        for(IntImage::iterator iter=img.begin(); iter!=img.end(); ++iter)
            *iter=5;
        initMultiArrayBorder(destMultiArrayRange(img),1,0);
        shouldEqualSequence(img.begin(), img.end(), desired_img);

        IntLine lin(IntLine::difference_type(7));
        
        for(IntLine::iterator iter=lin.begin(); iter!=lin.end(); ++iter)
            *iter=5;
        initMultiArrayBorder(destMultiArrayRange(lin),3,0);
        shouldEqualSequence(lin.begin(), lin.end(), desired_lin);

        IntVolume vol2(IntVolume::difference_type(2,2,2));
        
        for(IntVolume::iterator iter=vol2.begin(); iter!=vol2.end(); ++iter)
            *iter=5;
        initMultiArrayBorder(vol2, 9, 0);
        shouldEqualSequence(vol2.begin(), vol2.end(), desired_vol2);

    }

    void testInspect()
    {
        vigra::FindMinMax<PixelType> minmax;

        inspectMultiArray(img, minmax);

        shouldEqual(minmax.count, img.size());
        shouldEqual(minmax.min, 0.1f);
        shouldEqual(minmax.max, 59.1f);

        vigra::MultiArray<3, unsigned char> labels(img.shape());
        labels.subarray(Shape3(1,0,0), img.shape()-Shape3(1,0,0)) = 1;

        vigra::ArrayOfRegionStatistics<vigra::FindMinMax<PixelType> > stats(1);

        inspectTwoMultiArrays(img, labels, stats);

        shouldEqual(stats[0].count, 24);
        shouldEqual(stats[0].min, 0.1f);
        shouldEqual(stats[0].max, 59.1f);
        shouldEqual(stats[1].count, 36);
        shouldEqual(stats[1].min, 1.1f);
        shouldEqual(stats[1].max, 58.1f);
    }
    
    void testTensorUtilities()
    {
        MultiArrayShape<2>::type shape(3,4);
        int size = shape[0]*shape[1];
        
        MultiArray<2, TinyVector<double, 2> > vector(shape), rvector(shape);
        MultiArray<2, TinyVector<double, 3> > tensor1(shape), tensor2(shape), rtensor(shape);
        MultiArray<2, double > trace(shape), rtrace(shape);
        MultiArray<2, double > determinant(shape), rdet(shape);
        
        for(int k=0; k<size; ++k)
        {
            for(int l=0; l<2; ++l)
                vector[k][l] = randomMT19937().uniform();
            for(int l=0; l<3; ++l)
                tensor1[k][l] = randomMT19937().uniform();
            rdet[k] = tensor1[k][0]*tensor1[k][2] - sq(tensor1[k][1]);
        }
        
        vectorToTensor(srcImageRange(vector), destImage(rtensor));
        vectorToTensorMultiArray(srcMultiArrayRange(vector), destMultiArray(tensor2));
        shouldEqualSequence(tensor2.data(), tensor2.data()+size, rtensor.data());
        tensor2.init(TinyVector<double, 3>());
        vectorToTensorMultiArray(vector, tensor2);
        shouldEqualSequence(tensor2.data(), tensor2.data()+size, rtensor.data());
                
        tensorTrace(srcImageRange(tensor1), destImage(rtrace));
        tensorTraceMultiArray(srcMultiArrayRange(tensor1), destMultiArray(trace));
        shouldEqualSequence(trace.data(), trace.data()+size, rtrace.data());
        trace = 0;
        tensorTraceMultiArray(tensor1, trace);
        shouldEqualSequence(trace.data(), trace.data()+size, rtrace.data());
                
        tensorDeterminantMultiArray(srcMultiArrayRange(tensor1), destMultiArray(determinant));
        shouldEqualSequence(determinant.data(), determinant.data()+size, rdet.data());
        determinant = 0;
        tensorDeterminantMultiArray(tensor1, determinant);
        shouldEqualSequence(determinant.data(), determinant.data()+size, rdet.data());
                
        determinant = 1000.0;
        tensorDeterminantMultiArray(srcMultiArrayRange(tensor2), destMultiArray(determinant));
        shouldEqualTolerance(norm(determinant), 0.0, 1e-14);

        tensorEigenRepresentation(srcImageRange(tensor1), destImage(rtensor));
        tensorEigenvaluesMultiArray(srcMultiArrayRange(tensor1), destMultiArray(vector));
        shouldEqualSequenceTolerance(vector.begin(), vector.end(), rtensor.begin(), (TinyVector<double, 2>(1e-14)));

        vector = TinyVector<double, 2>();
        tensorEigenvaluesMultiArray(tensor1, vector);
        shouldEqualSequenceTolerance(vector.begin(), vector.end(), rtensor.begin(), (TinyVector<double, 2>(1e-14)));
    }
};

#endif

template <class T>
class ChunkedMultiArraySpeedTest
{
public:

    typedef ChunkedArray<3, T>                                      BaseArrayType;
    typedef ChunkedArrayFull<3, T>                                  FullArrayType;
    typedef ChunkedArrayTmpFile<3, T>                               TmpFileArrayType;
    typedef ChunkedArrayLazy<3, T>                                  AllocArrayType;
    typedef ChunkedArrayCompressed<3, T>                            CompressedArrayType;
    typedef ChunkedArrayHDF5<3, T>                                  HDF5ArrayType;
    typedef typename CoupledHandleType<3, ChunkedMemory<T> >::type  P1;
    typedef typename P1::base_type                                  P0;
    typedef CoupledScanOrderIterator<3, P1>                IteratorType;

    Shape3 s, outer_shape;
    
    ChunkedMultiArraySpeedTest()
    : s(200, 201, 202),
      outer_shape(chunkArrayShape(s, Shape3(6), Shape3(63)))
    {
        std::cerr << "chunked multi array test for type " << typeid(T).name() << ": \n";        
    }
    
    void testBaselineSpeed()
    {
        using namespace vigra::multi_math;
        using namespace vigra::functor;

        std::string t;
        T count = 1;
        USETICTOC;
        MultiArray<3, T> u(s);
        typename MultiArray<3, T>::iterator i   = u.begin(),
                                            end = i.getEndIterator();
        count = 1;
        for(; i != end; ++i, ++count)
            *i = count;
            
        i   = u.begin();
        count = 1;
        TIC;
        for(; i != end; ++i, ++count)
            if(count != *i)
            {
                shouldEqual(*i, count);
            }
        t = TOCS;
        std::cerr << "    contiguous array (baseline): " << t << "\n";
    }
    
    template <class Array>
    void testSpeedImpl(Array & v)
    {
        using namespace vigra::multi_math;
        using namespace vigra::functor;

        std::string t;
        T count = 1;
        USETICTOC;

        // IteratorType bi(P1(v, P0(s))),
                     // end = bi.getEndIterator();
        typename Array::iterator bi  = v.begin(),
                                 end = v.end();

        count = 1;
        int start = 0;
        // int stop = size;
        TIC;
        for(bi.setDim(2,start); bi.coord(2) < s[2]; bi.incDim(2))
            for(bi.setDim(1,start); bi.coord(1) < s[1]; bi.incDim(1))
                for(bi.setDim(0,start); bi.coord(0) < s[0]; bi.incDim(0), ++count)
                {
                    bi.template get<1>() = count;
                }
        t = TOCS;
        std::cerr << "    chunked iterator create and init: " << t << "\n";
        
        // bi = IteratorType(P1(v, P0(s)));
        bi = v.begin();
        count = 1;
        TIC;
        // for(bi.setDim(2,start); bi.coord(2) < s[2]; bi.incDim(2))
            // for(bi.setDim(1,start); bi.coord(1) < s[1]; bi.incDim(1))
                // for(bi.setDim(0,start); bi.coord(0) < s[0]; bi.incDim(0), ++count)
                // {
                    // if(bi.template get<1>() != count)
                        // shouldEqual(bi.template get<1>(), count);
                // }
        for(; bi != end; ++bi, ++count)
        {
            if(bi.template get<1>() != count)
            {
                std::cerr << bi.point() << "\n";
                shouldEqual(bi.template get<1>(), count);
            }
        }
        t = TOCS;
        std::cerr << "    chunked iterator read: " << t << "\n";
        
    }
    
    void testSpeedFullArray()
    {
        std::cerr << "    backend: full array\n";
        FullArrayType v(s);
        testSpeedImpl(v);
    }
    
    void testSpeedTmpFileAllCached()
    {
        int cache_max = prod(outer_shape);
        std::cerr << "    backend: tmp file, max cache: " << cache_max << "\n";
        TmpFileArrayType v(s, cache_max);
        testSpeedImpl(v);
    }
    
    void testSpeedTmpFileSliceCached()
    {
        int cache_max = outer_shape[0]*outer_shape[1];
        std::cerr << "    backend: tmp file, max cache: " << cache_max << "\n";
        TmpFileArrayType v(s, cache_max);
        testSpeedImpl(v);
    }
    
    void testSpeedTmpFileRowCached()
    {
        int cache_max = outer_shape[0];
        std::cerr << "    backend: tmp file, max cache: " << cache_max << "\n";
        TmpFileArrayType v(s, cache_max);
        testSpeedImpl(v);
    }
    
    void testSpeedAllocArray()
    {
        std::cerr << "    backend: alloc\n";
        AllocArrayType v(s);
        testSpeedImpl(v);
    }
    
    void testSpeedCompressedArrayLZ4()
    {
        int cache_max = outer_shape[0]*outer_shape[1];
        std::cerr << "    backend: LZ4, max cache: " << cache_max << "\n";
        CompressedArrayType v(s, cache_max, LZ4);
        testSpeedImpl(v);
    }
    
    void testSpeedCompressedArrayZlib()
    {
        int cache_max = outer_shape[0]*outer_shape[1];
        std::cerr << "    backend: ZLIB, max cache: " << cache_max << "\n";
        CompressedArrayType v(s, cache_max, ZLIB);
        testSpeedImpl(v);
    }
    
    void testSpeedCompressedArrayZlibFast()
    {
        int cache_max = outer_shape[0]*outer_shape[1];
        std::cerr << "    backend: ZLIB_FAST, max cache: " << cache_max << "\n";
        CompressedArrayType v(s, cache_max, ZLIB_FAST);
        testSpeedImpl(v);
    }
    
    void testSpeedHDF5Array()
    {
        int cache_max = outer_shape[0]*outer_shape[1];
        int compression = 1;
        std::cerr << "    backend: HDF5 (compression: " << compression << "), max cache: " << cache_max << "\n";
        
        std::string filename = "test_chunked_" + asString(sizeof(T)) + ".h5";
        HDF5File file(filename, HDF5File::New);
        HDF5ArrayType v(file, "test", s, cache_max, compression);
        testSpeedImpl(v);
    }

    template <class Array>
    void testSubarraySpeedImpl(Array & v)
    {
        using namespace vigra::multi_math;
        using namespace vigra::functor;

        std::string t;
        T count = 1;
        USETICTOC;

        typename Array::iterator bi  = v.begin(),
                                 end = v.end();

        count = 1;
        int start = 0;
        // int stop = size;
        TIC;
        for(bi.setDim(2,start); bi.coord(2) < s[2]; bi.incDim(2))
            for(bi.setDim(1,start); bi.coord(1) < s[1]; bi.incDim(1))
                for(bi.setDim(0,start); bi.coord(0) < s[0]; bi.incDim(0), ++count)
                {
                    bi.template get<1>() = count;
                }
        t = TOCS;
        std::cerr << "    chunked iterator create and init: " << t << "\n";
        
        // {
            // Shape3 start(63, 63, 63), stop(66, 66, 66);
            // typename Array::SubarrayType sub(stop - start);
            // v.copySubarray(Shape3(63, 63, 63), sub, true);
            // for(auto i = sub.begin(), end = sub.end(); i != end; ++i)
            // {
                // std::cerr << (UInt64)*i << " ";
                // if(i.point()[0] == 2)
                    // std::cerr << "\n";
            // }
            // std::cerr << "\n";
            
            // typename Array::ViewType sub2(stop - start);
            // v.viewSubarray(start, sub2);
            // MultiCoordinateIterator<3> c(stop - start),
                                       // cend(c.getEndIterator());
            // for(; c != cend; ++c)
            // {
                // std::cerr << (UInt64)sub2[c] << " ";
                // if(c.point()[0] == 2)
                    // std::cerr << "\n";
            // }
            // std::cerr << "\n";
            
            // MultiArrayView<2, T, ChunkedArrayTag> sub3(sub2.bindAt(1, 1));
            // MultiCoordinateIterator<2> cc(sub3.shape_),
                                       // ccend(cc.getEndIterator());
            // for(; cc != ccend; ++cc)
            // {
                // std::cerr << (UInt64)sub3[cc] << " ";
                // if(cc.point()[0] == 2)
                    // std::cerr << "\n";
            // }
            // std::cerr << "\n";
            
            // for(auto i = sub3.begin(), end = sub3.end(); i != end; ++i)
            // {
                // std::cerr << (UInt64)*i << " ";
                // if(i.point()[0] == 2)
                    // std::cerr << "\n";
            // }
            // std::cerr << "\n";
            
            // MultiArrayView<2, T, ChunkedArrayTag> 
                // sub4(sub3.subarray(Shape2(1,1), Shape2(3,3)));
            // for(auto i = sub4.begin(), end = sub4.end(); i != end; ++i)
            // {
                // std::cerr << (UInt64)*i << " ";
                // if(i.point()[0] == 1)
                    // std::cerr << "\n";
            // }
        // }

        MultiArrayView<3, T, ChunkedArrayTag> sub(v.shape());
        v.viewSubarray(Shape3(), sub);
        ChunkedSubarrayCopy<3, T> copy(v.shape());
        v.copySubarray(Shape3(), copy);
        // count = 1;
        // MultiCoordinateIterator<3> c(v.shape()), cend = c.getEndIterator();
        // TIC;
        // for(; c != cend; ++c, ++count)
        // {
            // if(sub4[*c] != count)
            // {
                // std::cerr << c << "\n";
                // shouldEqual(sub4[*c], count);
            // }
        // }
        // t = TOCS;
        // std::cerr << "    chunked subarray read: " << t << "\n";

        MultiArrayView<3, T, ChunkedArrayTag> trans(sub.transpose());
        MultiArrayView<3, T> transCopy(copy.transpose());
        MultiCoordinateIterator<3> c(trans.shape()), cend = c.getEndIterator();
        TIC;
        for(; c != cend; ++c)
        {
            if(trans[*c] != transCopy[*c])
            {
                std::cerr << c << "\n";
                shouldEqual(trans[*c], transCopy[*c]);
            }
        }
        t = TOCS;
        std::cerr << "    chunked subarray read: " << t << "\n";
        
    }
    
    void testSubarraySpeedAllocArray()
    {
        std::cerr << "    backend: alloc\n";
        AllocArrayType v(s);
        testSubarraySpeedImpl(v);
    }
    
    static void testMultiThreadedRun(BaseArrayType * v, int startIndex, int d, int * go)
    {
        while(*go == 0)
            threading::this_thread::yield();
            
        Shape3 s = v->shape();
        int sliceSize = s[0]*s[1];
        
        IteratorType bi(P1(*v, P0(s)));
        T count = 1 + startIndex*sliceSize;
        for(bi.setDim(2,startIndex); bi.coord(2) < s[2]; bi.addDim(2, d), count += (d-1)*sliceSize)
            for(bi.setDim(1,0); bi.coord(1) < s[1]; bi.incDim(1))
                for(bi.setDim(0,0); bi.coord(0) < s[0]; bi.incDim(0), ++count)
                {
                    bi.template get<1>() = count;
                }
    }

    template <class Array>
    void testMultiThreadedImpl(Array & v)
    {
        int go = 0;
        
        threading::thread t1(testMultiThreadedRun, &v, 0, 4, &go);
        threading::thread t2(testMultiThreadedRun, &v, 1, 4, &go);
        threading::thread t3(testMultiThreadedRun, &v, 2, 4, &go);
        threading::thread t4(testMultiThreadedRun, &v, 3, 4, &go);
     
        go = 1;
     
        t4.join();
        t3.join();
        t2.join();
        t1.join();
        
        T count = 1;

        IteratorType bi(P1(v, P0(s)));
        for(bi.setDim(2,0); bi.coord(2) < s[2]; bi.incDim(2))
            for(bi.setDim(1,0); bi.coord(1) < s[1]; bi.incDim(1))
                for(bi.setDim(0,0); bi.coord(0) < s[0]; bi.incDim(0), ++count)
                {
                    if(bi.template get<1>() != count)
                        shouldEqual(bi.template get<1>(), count);
                }
    }
    
    void testMultiThreaded()
    {
        int cache_max = outer_shape[0]*outer_shape[1];
        std::cerr << "    multi-threaded, backend: tmp file, max cache: " << cache_max << "\n";
        TmpFileArrayType v(s, cache_max);
        testMultiThreadedImpl(v);
    }
    
    void testMultiThreadedAlloc()
    {
        std::cerr << "    multi-threaded, backend: alloc\n";
        AllocArrayType v(s);
        testMultiThreadedImpl(v);
    }
    
    void testMultiThreadedCompressed()
    {
        int cache_max = outer_shape[0]*outer_shape[1];
        std::cerr << "    multi-threaded, backend: LZ4, max cache: " << cache_max << "\n";
        CompressedArrayType v(s, cache_max, LZ4);
        testMultiThreadedImpl(v);
    }
    
    void testMultiThreadedHDF5()
    {
        int cache_max = outer_shape[0]*outer_shape[1];
        std::cerr << "    multi-threaded, backend: HDF5, max cache: " << cache_max << "\n";
        
        std::string filename = "test_chunked_" + asString(sizeof(T)) + ".h5";
        HDF5File file(filename, HDF5File::New);
        
        HDF5ArrayType v(file, "test", s, cache_max);
        testMultiThreadedImpl(v);
    }

#if 0
    void testBasicArithmetic()
    {
        using namespace vigra::multi_math;
        using namespace vigra::functor;

        // test all overload variants
        r1 = b*d;
        rv = bv*dv;
        r4.transpose() = b.transpose()*d.transpose();
        r5 = b.transpose()*d.transpose();
        combineTwoMultiArrays(srcMultiArrayRange(b), srcMultiArray(d), destMultiArray(r2),
                              Arg1()*Arg2());
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());
        shouldEqualSequence(r4.begin(), r4.end(), r2.begin());
        shouldEqualSequence(r2.begin(), r2.end(), r5.transpose().begin());

        r1 = b*cv;
        rv = cv*b;
        combineTwoMultiArrays(srcMultiArrayRange(b), srcMultiArray(c), destMultiArray(r2),
                              Arg1()*Arg2());
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = 2.0*d;
        rv = 2.0*dv;
        r4.transpose() = 2.0*d.transpose();
        r5 = 2.0*d.transpose();
        transformMultiArray(srcMultiArrayRange(d), destMultiArray(r2),
                            Param(2.0)*Arg1());
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());
        shouldEqualSequence(r4.begin(), r4.end(), r2.begin());
        shouldEqualSequence(r2.begin(), r2.end(), r5.transpose().begin());

        r1 = d*4.0;
        rv = dv*4.0;
        transformMultiArray(srcMultiArrayRange(d), destMultiArray(r2),
                            Arg1()*Param(4.0));
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = (d*2.0)+b;
        rv = (dv*2.0)+bv;
        combineTwoMultiArrays(srcMultiArrayRange(b), srcMultiArray(d), destMultiArray(r2),
                              (Arg2()*Param(2.0))+Arg1());
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = b/(d*2.0);
        rv = bv/(dv*2.0);
        combineTwoMultiArrays(srcMultiArrayRange(b), srcMultiArray(d), destMultiArray(r2),
                              Arg1()/(Arg2()*Param(2.0)));
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = (2.0*b)-(d*2.0);
        rv = (2.0*bv)-(dv*2.0);
        combineTwoMultiArrays(srcMultiArrayRange(b), srcMultiArray(d), destMultiArray(r2),
                              (Param(2.0)*Arg1())-(Arg2()*Param(2.0)));
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = 2.0*(d*2.0);
        rv = 2.0*(dv*2.0);
        transformMultiArray(srcMultiArrayRange(d), destMultiArray(r2),
                            Param(2.0)*(Arg1()*Param(2.0)));
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = (d*3.0)*6.0;
        rv = (dv*3.0)*6.0;
        transformMultiArray(srcMultiArrayRange(d), destMultiArray(r2),
                            (Param(3.0)*Arg1())*Param(6.0));
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());


        r1 = sqrt(b);
        rv = sqrt(bv);
        r4.transpose() = sqrt(b.transpose());
        r5 = sqrt(b.transpose());
        transformMultiArray(srcMultiArrayRange(b), destMultiArray(r2),
                            sqrt(Arg1()));
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());
        shouldEqualSequence(r4.begin(), r4.end(), r2.begin());
        shouldEqualSequence(r2.begin(), r2.end(), r5.transpose().begin());

        r1 = sqrt(-d);
        rv = sqrt(-dv);
        transformMultiArray(srcMultiArrayRange(d), destMultiArray(r2),
                            sqrt(-Arg1()));
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = sqrt(b)*d;
        rv = sqrt(bv)*dv;
        combineTwoMultiArrays(srcMultiArrayRange(b), srcMultiArray(d), destMultiArray(r2),
                              sqrt(Arg1())*Arg2());
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = d*sqrt(b);
        rv = dv*sqrt(bv);
        combineTwoMultiArrays(srcMultiArrayRange(b), srcMultiArray(d), destMultiArray(r2),
                              Arg2()*sqrt(Arg1()));
        shouldEqualSequenceTolerance(r1.begin(), r1.end(), r2.begin(), 1e-15);
        shouldEqualSequenceTolerance(rv.begin(), rv.end(), r2.begin(), 1e-15);

        r1 = sqrt(b)*(d*2.0);
        rv = sqrt(bv)*(dv*2.0);
        combineTwoMultiArrays(srcMultiArrayRange(b), srcMultiArray(d), destMultiArray(r2),
                              sqrt(Arg1())*(Arg2()*Param(2.0)));
        shouldEqualSequenceTolerance(r1.begin(), r1.end(), r2.begin(), 1e-15);
        shouldEqualSequenceTolerance(rv.begin(), rv.end(), r2.begin(), 1e-15);

        r1 = (d*2.0)*sqrt(b);
        rv = (dv*2.0)*sqrt(bv);
        combineTwoMultiArrays(srcMultiArrayRange(b), srcMultiArray(d), destMultiArray(r2),
                              (Arg2()*Param(2.0))*sqrt(Arg1()));
        shouldEqualSequenceTolerance(r1.begin(), r1.end(), r2.begin(), 1e-15);
        shouldEqualSequenceTolerance(rv.begin(), rv.end(), r2.begin(), 1e-15);


        r1 = b*(-c);
        rv = bv*(-cv);
        combineTwoMultiArrays(srcMultiArrayRange(b), srcMultiArray(c), destMultiArray(r2),
                              Arg1()*(-Arg2()));
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = b*c*d;
        rv = bv*cv*dv;
        combineThreeMultiArrays(srcMultiArrayRange(b), srcMultiArray(c), srcMultiArray(d), destMultiArray(r2),
                                Arg1()*Arg2()*Arg3());
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = sqrt(b)*c*d*b;
        rv = sqrt(bv)*cv*dv*bv;
        combineThreeMultiArrays(srcMultiArrayRange(b), srcMultiArray(c), srcMultiArray(d), destMultiArray(r2),
                                sqrt(Arg1())*Arg2()*Arg3()*Arg1());
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());

        r1 = (b*c)*2.0*(abs(d)*b);
        rv = (bv*cv)*2.0*(abs(dv)*bv);
        combineThreeMultiArrays(srcMultiArrayRange(b), srcMultiArray(c), srcMultiArray(d), destMultiArray(r2),
                                (Arg1()*Arg2())*Param(2.0)*(abs(Arg3())*Arg1()));
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
        shouldEqualSequence(rv.begin(), rv.end(), r2.begin());
    
        shouldEqual(sum(b+0.5, 0.0), 300.0);
        shouldEqual(sum<double>(b+0.5), 300.0);
        shouldEqual(sum<double>(b), 288.0);

        shouldEqual(product(b.subarray(Shape3(1,0,0), Shape3(2,2,2))+0.5, 1.0), 3024.0);
        shouldEqual(product<double>(b.subarray(Shape3(1,0,0), Shape3(2,2,2))+0.5), 3024.0);
        shouldEqual(product<double>(MultiArray<3, float>(shape3)), 0.0);

        should(all(b > 0.0));
        should(!all(b > 10.0));

        should(any(b > 10.0));
        should(!any(b > 100.0));

        try 
        {
            array3_type fails(shape3_type(2, 7, 6));
            r1 = b * fails;
            failTest("shape mismatch exception not thrown");
        }
        catch(PreconditionViolation & e)
        {
            std::string expected("\nPrecondition violation!\nmulti_math: shape mismatch in expression.\n"),
                        actual(e.what());
            shouldEqual(actual.substr(0, expected.size()), expected);
        }
    }

#define VIGRA_TEST_UNARY_FUNCTION(FCT, RHS) \
        r1 = FCT(RHS); \
        for(int k=0; k<r2.size(); ++k) \
            r2[k] = FCT(RHS[k]); \
        shouldEqualSequence(r2.begin(), r2.end(), r1.begin()); 

#define VIGRA_TEST_BINARY_FUNCTION(FCT, MFCT, V1, V2) \
        r1 = MFCT(V1, V2); \
        for(int k=0; k<r2.size(); ++k) \
            r2[k] = FCT(V1[k], V2[k]); \
        shouldEqualSequence(r2.begin(), r2.end(), r1.begin()); 

#define VIGRA_TEST_BINARY_OPERATOR(OP, V1, V2) \
        r1 = V1 OP V2; \
        for(int k=0; k<r2.size(); ++k) \
            r2[k] = V1[k] OP V2[k]; \
        shouldEqualSequence(r2.begin(), r2.end(), r1.begin()); 

    void testAllFunctions()
    {
        using namespace vigra::multi_math;

        MultiArray<3, unsigned int> i(r2.shape()), j(r2.shape());
        linearSequence(i.begin(), i.end());
        j.init(2);

        VIGRA_TEST_UNARY_FUNCTION(-, b)
        VIGRA_TEST_UNARY_FUNCTION(!, i)
        VIGRA_TEST_UNARY_FUNCTION(~, i)

        VIGRA_TEST_UNARY_FUNCTION(abs, d)
        VIGRA_TEST_UNARY_FUNCTION(erf, b)
        VIGRA_TEST_UNARY_FUNCTION(even, i)
        VIGRA_TEST_UNARY_FUNCTION(odd, i)
        VIGRA_TEST_UNARY_FUNCTION(sign, d)
        VIGRA_TEST_UNARY_FUNCTION(signi, b)
        VIGRA_TEST_UNARY_FUNCTION(sq, b)
        VIGRA_TEST_UNARY_FUNCTION(norm, d)
        VIGRA_TEST_UNARY_FUNCTION(squaredNorm, d)
        VIGRA_TEST_UNARY_FUNCTION(vigra::multi_math::round, b)
        VIGRA_TEST_UNARY_FUNCTION(roundi, b)
        VIGRA_TEST_UNARY_FUNCTION(sqrti, i)
        VIGRA_TEST_UNARY_FUNCTION(sin_pi, b)
        VIGRA_TEST_UNARY_FUNCTION(cos_pi, b)
        VIGRA_TEST_UNARY_FUNCTION(vigra::multi_math::gamma, b)
        VIGRA_TEST_UNARY_FUNCTION(loggamma, b)
        VIGRA_TEST_UNARY_FUNCTION(sqrt, b)
        VIGRA_TEST_UNARY_FUNCTION(exp, b)
        VIGRA_TEST_UNARY_FUNCTION(log, b)
        VIGRA_TEST_UNARY_FUNCTION(log10, b)
        VIGRA_TEST_UNARY_FUNCTION(sin, b)
        VIGRA_TEST_UNARY_FUNCTION(asin, a)
        VIGRA_TEST_UNARY_FUNCTION(cos, b)
        VIGRA_TEST_UNARY_FUNCTION(acos, a)
        VIGRA_TEST_UNARY_FUNCTION(tan, b)
        VIGRA_TEST_UNARY_FUNCTION(atan, b)
        VIGRA_TEST_UNARY_FUNCTION(floor, b)
        VIGRA_TEST_UNARY_FUNCTION(ceil, b)

        VIGRA_TEST_BINARY_FUNCTION(atan2, atan2, b, c)
        VIGRA_TEST_BINARY_FUNCTION(pow, pow, b, c)
        VIGRA_TEST_BINARY_FUNCTION(fmod, fmod, b, c)
        VIGRA_TEST_BINARY_FUNCTION(std::min, minimum, b, c)
        VIGRA_TEST_BINARY_FUNCTION(std::max, maximum, b, c)
        VIGRA_TEST_BINARY_FUNCTION(std::min, multi_math::min, b, c)
        VIGRA_TEST_BINARY_FUNCTION(std::max, multi_math::max, b, c)

        VIGRA_TEST_BINARY_OPERATOR(+, b, c)
        VIGRA_TEST_BINARY_OPERATOR(-, b, c)
        VIGRA_TEST_BINARY_OPERATOR(*, b, c)
        VIGRA_TEST_BINARY_OPERATOR(/, b, c)
        VIGRA_TEST_BINARY_OPERATOR(%, i, j)
        VIGRA_TEST_BINARY_OPERATOR(&&, i, j)
        VIGRA_TEST_BINARY_OPERATOR(||, i, j)
        VIGRA_TEST_BINARY_OPERATOR(==, i, j)
        VIGRA_TEST_BINARY_OPERATOR(!=, i, j)
        VIGRA_TEST_BINARY_OPERATOR(<, i, j)
        VIGRA_TEST_BINARY_OPERATOR(<=, i, j)
        VIGRA_TEST_BINARY_OPERATOR(>, i, j)
        VIGRA_TEST_BINARY_OPERATOR(>=, i, j)
        VIGRA_TEST_BINARY_OPERATOR(<<, i, j)
        VIGRA_TEST_BINARY_OPERATOR(>>, i, j)
        VIGRA_TEST_BINARY_OPERATOR(&, i, j)
        VIGRA_TEST_BINARY_OPERATOR(|, i, j)
        VIGRA_TEST_BINARY_OPERATOR(^, i, j)
    }

#undef VIGRA_TEST_UNARY_FUNCTION
#undef VIGRA_TEST_BINARY_FUNCTION
#undef VIGRA_TEST_BINARY_OPERATOR

    void testMixedExpressions()
    {
        using namespace vigra::multi_math;

        r1 = 2.0 + b;
        r2 = 2 + b;
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());

        MultiArray<3, int> i(r2.shape());

        i = 2 + c;
        r1 = i;
        r2 = 2 + c;
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());
    }

    void testComputedAssignment()
    {
        using namespace vigra::multi_math;
        shouldEqual(r1.size(), 0);
        r1 += sq(b);
        rv.init(0.0);
        rv += sq(b);
        r2 = sq(b);
        shouldEqual(r1.size(), 24);
        shouldEqualSequence(r2.begin(), r2.end(), r1.begin());
        shouldEqualSequence(r2.begin(), r2.end(), rv.begin());

        r1 *= sq(b);
        rv *= sq(b);
        r2 = r2 * sq(b);
        shouldEqualSequence(r2.begin(), r2.end(), r1.begin());
        shouldEqualSequence(r2.begin(), r2.end(), rv.begin());

        r1 -= sq(b);
        rv -= sq(b);
        r2 = r2 - sq(b);
        shouldEqualSequence(r2.begin(), r2.end(), r1.begin());
        shouldEqualSequence(r2.begin(), r2.end(), rv.begin());

        r1 /= sq(b);
        rv /= sq(b);
        r2 = r2 / sq(b);
        shouldEqualSequence(r2.begin(), r2.end(), r1.begin());
        shouldEqualSequence(r2.begin(), r2.end(), rv.begin());
    }

    void testNonscalarValues()
    {
        using namespace vigra::multi_math;
        using namespace vigra::functor;

        typedef std::complex<double> C;
        MultiArray<3, C> q(r2.shape()), qr(q.shape()), qref(q.shape());
        linearSequence(q.begin(), q.end(), 1.0);

        qr = -q;
        transformMultiArray(srcMultiArrayRange(q), destMultiArray(qref),
                            -Arg1());
        shouldEqualSequence(qref.begin(), qref.end(), qr.begin());
        shouldEqual(qr[0], C(-1.0));

        qr = C(2.0, -2.0)*q;
        transformMultiArray(srcMultiArrayRange(q), destMultiArray(qref),
                            Param(C(2.0, -2.0))*Arg1());
        shouldEqualSequence(qref.begin(), qref.end(), qr.begin());
        shouldEqual(qr[1], C(4.0, -4.0));

        qr = 2.0*q;
        transformMultiArray(srcMultiArrayRange(q), destMultiArray(qref),
                            Param(2.0)*Arg1());
        shouldEqualSequence(qref.begin(), qref.end(), qr.begin());
        shouldEqual(qr[1], C(4.0));

        typedef RGBValue<int> R;
        MultiArray<3, R> r(r2.shape()), rr(q.shape()), rref(q.shape());
        linearSequence(r.begin(), r.end(), 1);

        rr = sq(r);
        transformMultiArray(srcMultiArrayRange(r), destMultiArray(rref),
                            sq(Arg1()));
        shouldEqualSequence(rref.begin(), rref.end(), rr.begin());
        shouldEqual(rr[1], R(4));

        rr = R(2)*r;
        transformMultiArray(srcMultiArrayRange(r), destMultiArray(rref),
                            Param(R(2))*Arg1());
        shouldEqualSequence(rref.begin(), rref.end(), rr.begin());
        shouldEqual(rr[0], R(2));

        rr = 4.0*r;
        transformMultiArray(srcMultiArrayRange(r), destMultiArray(rref),
                            Param(4.0)*Arg1());
        shouldEqualSequence(rref.begin(), rref.end(), rr.begin());
        shouldEqual(rr[0], R(4));
    }

    void testExpandMode()
    {
        using namespace vigra::multi_math;
        using namespace vigra::functor;

        array3_type s(shape3_type(1,1,1));
        s.init(3.0);

        r1 = s*d;
        transformMultiArray(srcMultiArrayRange(d), destMultiArray(r2),
                            Param(3.0)*Arg1());
        shouldEqualSequence(r1.begin(), r1.end(), r2.begin());

        MultiArray<1, double> ss(Shape1(d.shape(1)));
        linearSequence(ss.begin(), ss.end(), 1.0);

        r1 = ss.insertSingletonDimension(0).insertSingletonDimension(2) + d;

        for(int z=0; z<d.shape(2); ++z)
            for(int y=0; y<d.shape(1); ++y)
                for(int x=0; x<d.shape(0); ++x)
                    shouldEqual(d(x,y,z)+ss(y), r1(x,y,z));
    }

    void testComplex()
    {
        using namespace vigra::multi_math;
        MultiArray<3, std::complex<double> > ac(a.shape());
        ac.init(std::complex<double>(3.0, 4.0));

        MultiArray<3, std::complex<double> > bc = conj(ac);
        for(int z=0; z<bc.shape(2); ++z)
            for(int y=0; y<bc.shape(1); ++y)
                for(int x=0; x<bc.shape(0); ++x)
                    shouldEqual(bc(x,y,z), std::complex<double>(3.0, -4.0));

        bc = ac + ac;
        for(int z=0; z<bc.shape(2); ++z)
            for(int y=0; y<bc.shape(1); ++y)
                for(int x=0; x<bc.shape(0); ++x)
                    shouldEqual(bc(x,y,z), std::complex<double>(6.0, 8.0));

        a = real(ac);
        for(int z=0; z<a.shape(2); ++z)
            for(int y=0; y<a.shape(1); ++y)
                for(int x=0; x<a.shape(0); ++x)
                    shouldEqual(a(x,y,z), 3.0);

        a = imag(ac);
        for(int z=0; z<a.shape(2); ++z)
            for(int y=0; y<a.shape(1); ++y)
                for(int x=0; x<a.shape(0); ++x)
                    shouldEqual(a(x,y,z), 4.0);

        a = abs(ac);
        for(int z=0; z<a.shape(2); ++z)
            for(int y=0; y<a.shape(1); ++y)
                for(int x=0; x<a.shape(0); ++x)
                    shouldEqual(a(x,y,z), 5.0);

        a = arg(ac);
        for(int z=0; z<a.shape(2); ++z)
            for(int y=0; y<a.shape(1); ++y)
                for(int x=0; x<a.shape(0); ++x)
                    shouldEqualTolerance(a(x,y,z), std::atan2(4.0, 3.0), 1e-16);
    }
#endif
};

#if 0
struct MultiArrayTestSuite
: public vigra::test_suite
{
    MultiArrayTestSuite()
    : vigra::test_suite("MultiArrayTestSuite")
    {
        add( testCase( &MultiArrayTest::test_default_ctor ) );
        add( testCase( &MultiArrayTest::test_first_ctor ) );
        add( testCase( &MultiArrayTest::test_second_ctor ) );
        add( testCase( &MultiArrayTest::test_assignment ) );
        add( testCase( &MultiArrayTest::test_copy_construction ) );
        add( testCase( &MultiArrayTest::testShape ) );
        add( testCase( &MultiArrayTest::test_iterator ) );
        add( testCase( &MultiArrayTest::test_const_iterator ) );
        add( testCase( &MultiArrayTest::test_coupled_iterator ) );
        add( testCase( &MultiArrayTest::test_traverser ) );
        add( testCase( &MultiArrayTest::test_const_traverser ) );
        add( testCase( &MultiArrayTest::test_bindOuter ) );
        add( testCase( &MultiArrayTest::test_bindInner ) );
        add( testCase( &MultiArrayTest::test_bindInnerAll ) );
        add( testCase( &MultiArrayTest::test_bindAt ) );
        add( testCase( &MultiArrayTest::test_bind ) );
        add( testCase( &MultiArrayTest::test_reshape) );
        add( testCase( &MultiArrayTest::test_subarray ) );
        add( testCase( &MultiArrayTest::test_stridearray ) );
        add( testCase( &MultiArrayTest::test_copy_int_float ) );
        add( testCase( &MultiArrayTest::test_expandElements ) );

        add( testCase( &MultiImpexTest::testImpex ) );
#if defined(HasTIFF)
        add( testCase( &MultiImpexTest::testMultipageTIFF ) );
#endif
    }
};

struct MultiArrayDataTestSuite
: public vigra::test_suite
{
    MultiArrayDataTestSuite()
    : vigra::test_suite("MultiArrayDataTestSuite")
    {
        {
            typedef int T;
            add( testCase( &MultiArrayDataTest<T>::testHasData ) );
            add( testCase( &MultiArrayDataTest<T>::testEquality ) );
            add( testCase( &MultiArrayDataTest<T>::test_subarray ) );
            add( testCase( &MultiArrayDataTest<T>::test_stridearray ) );
            add( testCase( &MultiArrayDataTest<T>::test_bindOuter ) );
            add( testCase( &MultiArrayDataTest<T>::test_bindInner ) );
            add( testCase( &MultiArrayDataTest<T>::test_bindAt ) );
            add( testCase( &MultiArrayDataTest<T>::test_bind ) );
            add( testCase( &MultiArrayDataTest<T>::test_bind0 ) );
            add( testCase( &MultiArrayDataTest<T>::testIsUnstrided ) );
            add( testCase( &MultiArrayDataTest<T>::test_singletonDimension ) );
            add( testCase( &MultiArrayDataTest<T>::testPermute ) );
            add( testCase( &MultiArrayDataTest<T>::testMethods ) );
            add( testCase( &MultiArrayDataTest<T>::testScanOrderAccess ) );
            add( testCase( &MultiArrayDataTest<T>::testAssignmentAndReset ) );
            add( testCase( &MultiArrayNavigatorTest::testNavigator ) );
            add( testCase( &MultiArrayNavigatorTest::testCoordinateNavigator ) );
        }
        {
            typedef Multiband<int> T;
            add( testCase( &MultiArrayDataTest<T>::testHasData ) );
            add( testCase( &MultiArrayDataTest<T>::testEquality ) );
            add( testCase( &MultiArrayDataTest<T>::test_subarray ) );
            add( testCase( &MultiArrayDataTest<T>::test_stridearray ) );
            add( testCase( &MultiArrayDataTest<T>::test_bindOuter ) );
            add( testCase( &MultiArrayDataTest<T>::test_bindInner ) );
            add( testCase( &MultiArrayDataTest<T>::test_bindAt ) );
            add( testCase( &MultiArrayDataTest<T>::test_bind ) );
            add( testCase( &MultiArrayDataTest<T>::test_bind0 ) );
            add( testCase( &MultiArrayDataTest<T>::testIsStrided ) );
            add( testCase( &MultiArrayDataTest<T>::test_singletonDimension ) );
            add( testCase( &MultiArrayDataTest<T>::testPermute ) );
            add( testCase( &MultiArrayDataTest<T>::testMethods ) );
            add( testCase( &MultiArrayDataTest<T>::testScanOrderAccess ) );
            add( testCase( &MultiArrayDataTest<T>::testAssignmentAndReset ) );
        }
    }
};

#endif

struct ChunkedMultiArraySpeedTestTestSuite
: public vigra::test_suite
{
  ChunkedMultiArraySpeedTestTestSuite()
    : vigra::test_suite("ChunkedMultiArraySpeedTestTestSuite")
    {
        {
            typedef ChunkedArrayFull<3, float> T;
            add( testCase( &ChunkedMultiArrayTest<T>::testEquality ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_assignment ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_bindAt ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_subarray ) );
        }

        {
            typedef ChunkedArrayLazy<3, float> T;
            add( testCase( &ChunkedMultiArrayTest<T>::testEquality ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_assignment ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_bindAt ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_subarray ) );
        }

        {
            typedef ChunkedArrayCompressed<3, float> T;
            add( testCase( &ChunkedMultiArrayTest<T>::testEquality ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_assignment ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_bindAt ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_subarray ) );
        }

        {
            typedef ChunkedArrayHDF5<3, float> T;
            add( testCase( &ChunkedMultiArrayTest<T>::testEquality ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_assignment ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_bindAt ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_subarray ) );
        }

        {
            typedef ChunkedArrayTmpFile<3, float> T;
            add( testCase( &ChunkedMultiArrayTest<T>::testEquality ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_assignment ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_bindAt ) );
            add( testCase( &ChunkedMultiArrayTest<T>::test_subarray ) );
        }

            //add( testCase( &MultiArrayPointoperatorsTest::testInit ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testCopy ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testCopyOuterExpansion ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testCopyInnerExpansion ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testTransform ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testTransformOuterExpand ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testTransformInnerExpand ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testTransformOuterReduce ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testTransformInnerReduce ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testCombine2 ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testCombine2OuterExpand ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testCombine2InnerExpand ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testCombine2OuterReduce ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testCombine2InnerReduce ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testCombine3 ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testInitMultiArrayBorder ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testInspect ) );
        //add( testCase( &MultiArrayPointoperatorsTest::testTensorUtilities ) );

        add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testBaselineSpeed ) );
        add( testCase( &ChunkedMultiArraySpeedTest<float>::testBaselineSpeed ) );
        add( testCase( &ChunkedMultiArraySpeedTest<double>::testBaselineSpeed ) );
        
        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testSpeedFullArray ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testSpeedFullArray ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testSpeedFullArray ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testSpeedAllocArray ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testSpeedAllocArray ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testSpeedAllocArray ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testSpeedTmpFileAllCached ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testSpeedTmpFileAllCached ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testSpeedTmpFileAllCached ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testSpeedTmpFileSliceCached ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testSpeedTmpFileSliceCached ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testSpeedTmpFileSliceCached ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testSpeedTmpFileRowCached ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testSpeedTmpFileRowCached ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testSpeedTmpFileRowCached ) );


        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testSpeedCompressedArrayLZ4 ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testSpeedCompressedArrayLZ4 ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testSpeedCompressedArrayLZ4 ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testSpeedCompressedArrayZlib ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testSpeedCompressedArrayZlib ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testSpeedCompressedArrayZlib ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testSpeedCompressedArrayZlibFast ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testSpeedCompressedArrayZlibFast ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testSpeedCompressedArrayZlibFast ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testSpeedHDF5Array ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testSpeedHDF5Array ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testSpeedHDF5Array ) );
        
        add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testSubarraySpeedAllocArray ) );
        add( testCase( &ChunkedMultiArraySpeedTest<float>::testSubarraySpeedAllocArray ) );
        add( testCase( &ChunkedMultiArraySpeedTest<double>::testSubarraySpeedAllocArray ) );
        
        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testMultiThreaded ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testMultiThreaded ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testMultiThreaded ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testMultiThreadedAlloc ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testMultiThreadedAlloc ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testMultiThreadedAlloc ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testMultiThreadedCompressed ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testMultiThreadedCompressed ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testMultiThreadedCompressed ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<UInt8>::testMultiThreadedHDF5 ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<float>::testMultiThreadedHDF5 ) );
        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testMultiThreadedHDF5 ) );

        // add( testCase( &ChunkedMultiArraySpeedTest<double>::testSpeed2 ) );
        //add( testCase( &MultiMathTest::testBasicArithmetic ) );
        //add( testCase( &MultiMathTest::testExpandMode ) );
        //add( testCase( &MultiMathTest::testAllFunctions ) );
        //add( testCase( &MultiMathTest::testComputedAssignment ) );
        //add( testCase( &MultiMathTest::testNonscalarValues ) );
        //add( testCase( &MultiMathTest::testMixedExpressions ) );
        //add( testCase( &MultiMathTest::testComplex ) );
    }
}; // struct MultiArrayPointOperatorsTestSuite


int main(int argc, char ** argv)
{
    int failed = 0;

    // run the multi-array point operator test suite
    ChunkedMultiArraySpeedTestTestSuite test0;
    failed += test0.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test0.report() << std::endl;

    //// run the multi-array testsuite
    //MultiArrayTestSuite test1;
    //int failed = test1.run(vigra::testsToBeExecuted(argc, argv));
    //std::cout << test1.report() << std::endl;

    //// run the multi-array data-testsuite
    //MultiArrayDataTestSuite test1a;
    //failed += test1a.run(vigra::testsToBeExecuted(argc, argv));
    //std::cout << test1a.report() << std::endl;
    
    return (failed != 0);
}


