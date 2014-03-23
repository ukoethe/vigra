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
#ifdef HasHDF5
#include "vigra/multi_array_chunked_hdf5.hxx"
#endif
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
    typedef MultiArray <3, T> PlainArray;
    typedef ChunkedArray<3, T> BaseArray;
    typedef VIGRA_UNIQUE_PTR<BaseArray> ArrayPtr;
    typedef typename BaseArray::iterator Iterator;
    
    Shape3 shape, chunk_shape;
    ArrayPtr array;
    PlainArray ref;

    ChunkedMultiArrayTest ()
        : shape(20,21,22),
          chunk_shape(8),
          ref(shape)
    {
        linearSequence(ref.begin(), ref.end());
        array = createArray(shape, chunk_shape, (Array *)0);
        linearSequence(array->begin(), array->end());
    }
    
    static ArrayPtr createArray(Shape3 const & shape, 
                                Shape3 const & chunk_shape,
                                ChunkedArrayFull<3, T> *)
    {
        return ArrayPtr(new ChunkedArrayFull<3, T>(shape));
    }
    
    static ArrayPtr createArray(Shape3 const & shape, 
                                Shape3 const & chunk_shape,
                                ChunkedArrayLazy<3, T> *)
    {
        return ArrayPtr(new ChunkedArrayLazy<3, T>(shape, chunk_shape));
    }
    
    static ArrayPtr createArray(Shape3 const & shape, 
                                Shape3 const & chunk_shape,
                                ChunkedArrayCompressed<3, T> *)
    {
        return ArrayPtr(new ChunkedArrayCompressed<3, T>(LZ4, shape, chunk_shape));
    }
    
#ifdef HasHDF5
    static ArrayPtr createArray(Shape3 const & shape, 
                                Shape3 const & chunk_shape,
                                ChunkedArrayHDF5<3, T> *)
    {
        HDF5File hdf5_file("chunked_test.h5", HDF5File::New);
        return ArrayPtr(new ChunkedArrayHDF5<3, T>(hdf5_file, "test", HDF5File::New, 
                                                   shape, chunk_shape));
    }
#endif
    
    static ArrayPtr createArray(Shape3 const & shape, 
                                Shape3 const & chunk_shape,
                                ChunkedArrayTmpFile<3, T> *)
    {
        return ArrayPtr(new ChunkedArrayTmpFile<3, T>(shape, chunk_shape, -1, ""));
    }
    
    void test_construction ()
    {
        should(array->isInside(Shape3(1,2,3)));
        should(!array->isInside(Shape3(1,23,3)));
        should(!array->isInside(Shape3(1,2,-3)));
        
        shouldEqual(array->shape(), ref.shape());
        shouldEqual(array->shape(0), ref.shape(0));
        shouldEqual(array->shape(1), ref.shape(1));
        shouldEqual(array->shape(2), ref.shape(2));
        
        if(IsSameType<Array, ChunkedArrayFull<3, T> >::value)
        {
            shouldEqual(array->chunkArrayShape(), Shape3(1));
        }
        else
        {
            shouldEqual(array->chunkArrayShape(), Shape3(3));
        }
        
        shouldEqualSequence(array->begin(), array->end(), ref.begin());
        
        should(*array == ref);
        should(*array != ref.subarray(Shape3(1),ref.shape()));
        
        shouldEqual(array->getItem(Shape3(1,8,17)), ref[Shape3(1,8,17)]);
        
        ref[ref.size()-1] = ref[ref.size()-1] + T(1);
        should(*array != ref);
        array->setItem(ref.shape()-Shape3(1), ref[ref.size()-1]);
        should(*array == ref);

        // FIXME: test copy construction?
        
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

    void test_assignment()
    {
        MultiArrayView <3, T, ChunkedArrayTag> v;
        should(!v.hasData());    
        
        v = array->subarray(Shape3(), ref.shape());
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
        
        MultiArrayView <3, T, ChunkedArrayTag> vs(array->subarray(Shape3(), Shape3(4)));
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
       
        vc /= PlainArray(ref.shape(), T(1));
        shouldEqual(vc.shape(), ref.shape());
        shouldEqualSequence(vc.begin(), vc.end(), ref.begin());
    }

    void test_bindAt ()
    {
        MultiArrayView <2, T, ChunkedArrayTag> v = array->bindAt (1, 4);
        MultiArrayView <2, T, ChunkedArrayTag> vv = array->template bind<1>(4);
        MultiArrayView <2, T, StridedArrayTag> vr = ref.bindAt (1, 4);
        
        shouldEqual(v.shape(), vr.shape());
        shouldEqual(vv.shape(), vr.shape());
        should(v == vr);
        should(vv == vr);
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
    
    void test_bindInner ()
    {
        MultiArrayView <2, T, ChunkedArrayTag> v = array->bindInner(2);
        MultiArrayView <2, T, StridedArrayTag> vr = ref.bindInner(2);
        shouldEqual(v.shape(), vr.shape());
        should(v == vr);
        
        TinyVector <int, 2> inner_indices (2, 5);
        MultiArrayView <1, T, ChunkedArrayTag> v1 = array->bindInner(inner_indices);
        MultiArrayView <1, T, StridedArrayTag> v1r = ref.bindInner(inner_indices);
        shouldEqual(v1.shape(), v1r.shape());
        should(v1 == v1r);
        
        MultiArrayView <1, T, ChunkedArrayTag> v21 = v.bindInner(5);
        shouldEqual(v21.shape(), v1r.shape());
        should(v21 == v1r);
    }
    
    void test_bindOuter ()
    {
        MultiArrayView <2, T, ChunkedArrayTag> v = array->bindOuter(2);
        MultiArrayView <2, T, StridedArrayTag> vr = ref.bindOuter(2);
        shouldEqual(v.shape(), vr.shape());
        should(v == vr);
        
        TinyVector <int, 2> inner_indices (5, 2);
        MultiArrayView <1, T, ChunkedArrayTag> v1 = array->bindOuter(inner_indices);
        MultiArrayView <1, T, StridedArrayTag> v1r = ref.bindOuter(inner_indices);
        shouldEqual(v1.shape(), v1r.shape());
        should(v1 == v1r);
        
        MultiArrayView <1, T, ChunkedArrayTag> v21 = v.bindOuter(5);
        shouldEqual(v21.shape(), v1r.shape());
        should(v21 == v1r);
    }

    void test_subarray ()
    {
        {
            Shape3 start, stop(ref.shape());  // whole array
            MultiArrayView <3, T, ChunkedArrayTag> v(array->subarray(start, stop));
            MultiArrayView <3, T, ChunkedArrayTag> vt(v.transpose());

            MultiArray <3, T> c(stop-start);
            array->checkoutSubarray(start, c);
            
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
            MultiArrayView <3, T, ChunkedArrayTag> v(array->subarray(start, stop));
            MultiArrayView <3, T, ChunkedArrayTag> vt(v.transpose());

            MultiArray <3, T> c(stop-start);
            array->checkoutSubarray(start, c);
            
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
            MultiArrayView <3, T, ChunkedArrayTag> v(array->subarray(start, stop));
            MultiArrayView <3, T, ChunkedArrayTag> vt(v.transpose());
            
            MultiArray <3, T> c(stop-start);
            array->checkoutSubarray(start, c);
            
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
    }
    
    void test_iterator ()
    {
        Shape3 s(ref.shape());
        typedef typename ChunkedArray<3, T>::iterator Iterator;
        MultiArrayView <3, T, ChunkedArrayTag> v(array->subarray(Shape3(), s));
        Iterator i1 = array->begin();
        Iterator iend = array->end();
        MultiCoordinateIterator<3> c(s),
                                   cend = c.getEndIterator();

        should(i1.isValid() && !i1.atEnd());
        should(!iend.isValid() && iend.atEnd());
        should(iend.getEndIterator() == iend);
        
        shouldEqual(i1.point(), *c);
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

        shouldEqual(&i1[0], &v[Shape3(0,0,0)]);
        shouldEqual(&i1[1], &v[Shape3(1,0,0)]);
        shouldEqual(&i1[s[0]], &v[Shape3(0,1,0)]);
        shouldEqual(&i1[s[0]*9+1], &v[Shape3(1,9,0)]);
        shouldEqual(&i1[s[0]*s[1]], &v[Shape3(0,0,1)]);
        shouldEqual(&i1[s[0]*s[1]*9+1], &v[Shape3(1,0,9)]);
        shouldEqual(&i1[(s[0]+1)*s[1]], &v[Shape3(1,1,1)]);

        shouldEqual(&*(i1+0), &v[Shape3(0,0,0)]);
        shouldEqual(&*(i1+1), &v[Shape3(1,0,0)]);
        shouldEqual(&*(i1+s[0]), &v[Shape3(0,1,0)]);
        shouldEqual(&*(i1+s[0]*9+1), &v[Shape3(1,9,0)]);
        shouldEqual(&*(i1+s[0]*s[1]), &v[Shape3(0,0,1)]);
        shouldEqual(&*(i1+s[0]*s[1]*9+1), &v[Shape3(1,0,9)]);
        shouldEqual(&*(i1+(s[0]+1)*s[1]), &v[Shape3(1,1,1)]);

        shouldEqual(&*(i1+Shape3(0,0,0)), &v[Shape3(0,0,0)]);
        shouldEqual(&*(i1+Shape3(1,0,0)), &v[Shape3(1,0,0)]);
        shouldEqual(&*(i1+Shape3(0,1,0)), &v[Shape3(0,1,0)]);
        shouldEqual(&*(i1+Shape3(1,11,0)), &v[Shape3(1,11,0)]);
        shouldEqual(&*(i1+Shape3(0,0,1)), &v[Shape3(0,0,1)]);
        shouldEqual(&*(i1+Shape3(1,0,11)), &v[Shape3(1,0,11)]);
        shouldEqual(&*(i1+Shape3(1,1,1)), &v[Shape3(1,1,1)]);

        shouldEqual(&*(iend-1),  &v[Shape3(19,20,21)]);
        shouldEqual(&*(iend-2),  &v[Shape3(18,20,21)]);
        shouldEqual(&*(iend-10), &v[Shape3(10,20,21)]);
        shouldEqual(&*(iend-s[0]-1),  &v[Shape3(19,19,21)]);

        shouldEqual(&iend[-1], &v[Shape3(19,20,21)]);
        shouldEqual(&iend[-2], &v[Shape3(18,20,21)]);
        shouldEqual(&iend[-10], &v[Shape3(10,20,21)]);
        shouldEqual(&iend[-s[0]-1], &v[Shape3(19,19,21)]);
        
        Iterator i2;
        i2 = iend;
        should(i2 == iend);
        should(!i2.isValid() && i2.atEnd());
        --i2;
        should(i2.isValid() && !i2.atEnd());
        should(i2.getEndIterator() == iend);
        shouldEqual(i2.point(), Shape3(19,20,21));
        shouldEqual(&*i2, &v[Shape3(19,20,21)]);
        for(int k=0; k<20; ++k)
            --i2;
        should(i2.isValid() && !i2.atEnd());
        should(i2.getEndIterator() == iend);
        shouldEqual(i2.point(), Shape3(19,19,21));
        shouldEqual(&*i2, &v[Shape3(19,19,21)]);
        for(int k=0; k<420; ++k)
            --i2;
        should(i2.isValid() && !i2.atEnd());
        should(i2.getEndIterator() == iend);
        shouldEqual(i2.point(), Shape3(19,19,20));
        shouldEqual(&*i2, &v[Shape3(19,19,20)]);

        i2 = iend-1;
        shouldEqual(&*(i2-Shape3(0,0,0)), &v[Shape3(19,20,21)]);
        shouldEqual(&*(i2-Shape3(1,0,0)), &v[Shape3(18,20,21)]);
        shouldEqual(&*(i2-Shape3(0,1,0)), &v[Shape3(19,19,21)]);
        shouldEqual(&*(i2-Shape3(9,1,0)), &v[Shape3(10,19,21)]);
        shouldEqual(&*(i2-Shape3(0,0,1)), &v[Shape3(19,20,20)]);
        shouldEqual(&*(i2-Shape3(9,0,1)), &v[Shape3(10,20,20)]);
        shouldEqual(&*(i2-Shape3(9,9,1)), &v[Shape3(10,11,20)]);
        shouldEqual(&*(i2-Shape3(9,9,9)), &v[Shape3(10,11,12)]);

        unsigned int count = 0;
        Shape3 p;
        i2 = array->begin();
        Iterator i3 = array->begin();
        Iterator i4 = array->begin();
        Iterator i5 = array->begin();
        Iterator i6 = array->begin();

        for (p[2]=0, i3.resetDim(2), i4.setDim(2, 0), i5.template dim<2>() = 0, i6.resetDim(2); 
                i3.point(2) != s[2]; 
                i3.incDim(2), i4.addDim(2, 1), ++i5.template dim<2>(), i6.template dim<2>() += 1, ++p[2]) 
        {
            for (p[1]=0, i3.resetDim(1), i4.setDim(1, 0), i5.template dim<1>() = 0, i6.resetDim(1); 
                    i3.point(1) != s[1]; 
                    i3.incDim(1), i4.addDim(1, 1), ++i5.template dim<1>(), i6.template dim<1>() += 1, ++p[1]) 
            {
                for (p[0]=0, i3.resetDim(0), i4.setDim(0, 0), i5.template dim<0>() = 0, i6.resetDim(0); 
                        i3.point(0) != s[0]; 
                        i3.incDim(0), i4.addDim(0, 1), ++i5.template dim<0>(), i6.template dim<0>() += 1, ++p[0], ++i1, ++c, i2 += 1, ++count)
                {
                    shouldEqual(&*i1, &v[p]);
                    shouldEqual(&*i2, &v[p]);
                    shouldEqual(&*i3, &v[p]);
                    shouldEqual(&*i4, &v[p]);
                    shouldEqual(&*i5, &v[p]);
                    shouldEqual(&*i6, &v[p]);
                    shouldEqual(i1.operator->(), &v[p]);
                    shouldEqual(i2.operator->(), &v[p]);
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

                    should(i5.template dim<2>() == p[2]);
                    should(i5.template dim<1>() == p[1]);
                    should(i5.template dim<0>() == p[0]);
                    should(i5.template dim<2>() != s[2]);
                    should(i5.template dim<1>() != s[1]);
                    should(i5.template dim<0>() != s[0]);
                    should(i5.template dim<2>() <= p[2]);
                    should(i5.template dim<1>() <= p[1]);
                    should(i5.template dim<0>() <= p[0]);
                    should(i5.template dim<2>() < s[2]);
                    should(i5.template dim<1>() < s[1]);
                    should(i5.template dim<0>() < s[0]);
                    should(i5.template dim<2>() >= 0);
                    should(i5.template dim<1>() >= 0);
                    should(i5.template dim<0>() >= 0);
                    shouldNot(i5.template dim<2>() > s[2]);
                    shouldNot(i5.template dim<1>() > s[1]);
                    shouldNot(i5.template dim<0>() > s[0]);

                    shouldEqual(iend - i1, v.size() - count);

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
        shouldEqual (count, v.size());

        --i1;
        i2 -= 1;
        shouldEqual(&*i1, &v[Shape3(19,20,21)]);
        shouldEqual(&*i2, &v[Shape3(19,20,21)]);
    }
    
    static void testMultiThreadedRun(BaseArray * v, int startIndex, int d, int * go)
    {
        while(*go == 0)
            threading::this_thread::yield();
            
        Shape3 s = v->shape();
        int sliceSize = s[0]*s[1];
        
        Iterator bi(v->begin());
        T count(startIndex*sliceSize), start((d-1)*sliceSize), inc(1);
        for(bi.setDim(2,startIndex); bi.coord(2) < s[2]; bi.addDim(2, d), count += start)
            for(bi.setDim(1,0); bi.coord(1) < s[1]; bi.incDim(1))
                for(bi.setDim(0,0); bi.coord(0) < s[0]; bi.incDim(0), count += inc)
                {
                    *bi = count;
                }
    }

    void testMultiThreaded()
    {
        array.reset(0); // close the file if backend is HDF5
        ArrayPtr a = createArray(Shape3(200, 201, 202), Shape3(), (Array *)0);
    
        int go = 0;
        
        threading::thread t1(testMultiThreadedRun, a.get(), 0, 4, &go);
        threading::thread t2(testMultiThreadedRun, a.get(), 1, 4, &go);
        threading::thread t3(testMultiThreadedRun, a.get(), 2, 4, &go);
        threading::thread t4(testMultiThreadedRun, a.get(), 3, 4, &go);
     
        go = 1;
     
        t4.join();
        t3.join();
        t2.join();
        t1.join();
        
        PlainArray ref(a->shape());
        linearSequence(ref.begin(), ref.end());
        
        shouldEqualSequence(a->begin(), a->end(), ref.begin());
    }
        
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

    // void test_expandElements()
    // {
        // using namespace multi_math;

        // MultiArray<3, TinyVector<int, 3> > a(Shape3(4,3,2));
        // a.init(TinyVector<int, 3>(1,2,3));

        // MultiArrayView<4, int, StridedArrayTag> ex = a.expandElements(0);
        // MultiArrayView<4, int, StridedArrayTag>::iterator i = ex.begin();
        // while(i != ex.end())
        // {
            // shouldEqual(*i, 1); ++i;
            // shouldEqual(*i, 2); ++i;
            // shouldEqual(*i, 3); ++i;
        // }

        // MultiArrayView<4, int, StridedArrayTag> ex2 = a.expandElements(3);
        // i = ex2.begin();
        // for(int k=0; k < a.size(); ++i, ++k)
            // shouldEqual(*i, 1);
        // for(int k=0; k < a.size(); ++i, ++k)
            // shouldEqual(*i, 2);
        // for(int k=0; k < a.size(); ++i, ++k)
            // shouldEqual(*i, 3);

        // MultiArray<3, bool> b = (a.bindElementChannel(0) == 1);
        // should(b.all());
        // b = (a.bindElementChannel(1) == 2);
        // should(b.all());
        // b = (a.bindElementChannel(2) == 3);
        // should(b.all());
    // }
};

// struct MultiArrayPointoperatorsTest
// {

    // typedef float PixelType;
    // typedef MultiArray<3,PixelType> Image3D;
    // typedef MultiArrayView<3,PixelType> View3D;
    // typedef Image3D::difference_type Size3;
    // typedef MultiArray<1,PixelType> Image1D;
    // typedef Image1D::difference_type Size1;

    // Image3D img;

    // MultiArrayPointoperatorsTest()
    // : img(Size3(5,4,3))
    // {
        // int i;
        // PixelType c = 0.1f;
        // for(i=0; i<img.elementCount(); ++i, ++c)
            // img.data()[i] = c;
    // }

    // void testInit()
    // {
        // Image3D res(img.shape());
        // const Image3D::value_type ini = 1.1f;
        // should(res.shape() == Size3(5,4,3));

        // initMultiArray(destMultiArrayRange(res), ini);

        // int x,y,z;
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), ini);

        // using namespace multi_math;
        // should(all(res == ini));

        // initMultiArray(res, 2.2f);
        // should(all(res == 2.2f));

        // res = 3.3f;
        // should(all(res == 3.3f));

        // res.init(4.4f);
        // should(all(res == 4.4f));
    // }

    // void testCopy()
    // {
        // Image3D res(img.shape(), 1.0), res1(img.shape(), 1.0);
        
        // copyMultiArray(srcMultiArrayRange(img), destMultiArray(res));
        // copyMultiArray(img, res1);

        // should(img == res);
        // should(img == res1);
    // }

    // void testCopyOuterExpansion()
    // {
        // Image3D res(img.shape());

        // copyMultiArray(img.subarray(Size3(0,0,0), Size3(5,1,1)), res);
        
        // int x,y,z;
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), img(x,0,0));
    // }

    // void testCopyInnerExpansion()
    // {
        // Image3D res(img.shape());

        // copyMultiArray(img.subarray(Size3(0,0,0), Size3(1,1,3)), res);
        
        // int x,y,z;
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), img(0,0,z));
    // }

    // void testTransform()
    // {
        // Image3D res(img.shape()), res1(img.shape());
        // transformMultiArray(srcMultiArrayRange(img), destMultiArray(res),
                            // Arg1() + Arg1());
        // transformMultiArray(img, res1, Arg1() + Arg1());
        
        // using namespace multi_math;
        // should(all(2.0*img == res));
        // should(all(2.0*img == res1));
    // }

    // void testTransformOuterExpand()
    // {
        // Image3D res(img.shape());
        // transformMultiArray(img.subarray(Size3(0,0,0), Size3(5,1,1)), res,
                            // Arg1() + Arg1());
        
        // int x,y,z;
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), 2.0*img(x,0,0));
    // }

    // void testTransformInnerExpand()
    // {
        // Image3D res(img.shape());

        // transformMultiArray(img.subarray(Size3(0,0,0), Size3(1,1,3)), res,
                            // Arg1() + Arg1());
        
        // int x,y,z;
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), 2.0*img(0,0,z));
    // }

    // void testTransformOuterReduce()
    // {
        // Image3D res(Size3(5,1,1));

        // transformMultiArray(img, res, reduceFunctor(Arg1() + Arg2(), 0.0));
        
        // int x,y,z;
        // for(x=0; x<img.shape(0); ++x)
        // {
            // double sum = 0.0;
            // for(y=0; y<img.shape(1); ++y)
                // for(z=0; z<img.shape(2); ++z)
                    // sum += img(x,y,z);
            // shouldEqual(res(x,0,0), sum);
        // }
        
        // Image1D res1(Size1(5));
        // MultiArrayView<3,PixelType> res3 = res1.insertSingletonDimension(1).insertSingletonDimension(2);
        // transformMultiArray(img, res3, FindSum<PixelType>());
        // shouldEqualSequenceTolerance(res1.data(), res1.data()+5, res.data(), 1e-7);       
    // }

    // void testTransformInnerReduce()
    // {
        // Image3D res(Size3(1,1,3));
        
        // transformMultiArray(img, res, reduceFunctor(Arg1() + Arg2(), 0.0));
        
        // int x,y,z;
        // for(z=0; z<img.shape(2); ++z)
        // {
            // double sum = 0.0;
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // sum += img(x,y,z);
            // shouldEqual(res(0,0,z), sum);
        // }
        
        // Image1D res1(Size1(3));
        // MultiArrayView<3,PixelType> res3 = res1.insertSingletonDimension(0).insertSingletonDimension(0);
        // transformMultiArray(img, res3, FindSum<PixelType>());
        // shouldEqualSequenceTolerance(res1.data(), res1.data()+3, res.data(), 1e-6);       
    // }

    // void testCombine2()
    // {
        // Image3D res(img.shape()), res1(img.shape());
        
        // combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArray(img), 
                              // destMultiArray(res),
                              // Arg1() + Arg2());
        // combineTwoMultiArrays(img, img, res1, Arg1() + Arg2());
        
        // using namespace multi_math;
        // should(all(2.0*img == res));
        // should(all(2.0*img == res1));
    // }

    // void testCombine2OuterExpand()
    // {
        // Image3D res(img.shape());
        
        // combineTwoMultiArrays(img.subarray(Size3(0,0,0), Size3(5,1,1)), img, res,
                              // Arg1() + Param(2.0)*Arg2());       
        // int x,y,z;
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), 2.0*img(x,y,z) + img(x,0,0));

        // combineTwoMultiArrays(img, img.subarray(Size3(0,0,0), Size3(5,1,1)), res,
                              // Arg1() + Param(2.0)*Arg2());       
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), img(x,y,z) + 2.0*img(x,0,0));

        // View3D view = img.subarray(Size3(0,0,0), Size3(5,1,1));
        // combineTwoMultiArrays(srcMultiArrayRange(view), srcMultiArrayRange(view), 
                              // destMultiArrayRange(res),
                              // Arg1() + Param(2.0)*Arg2());       
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), 3.0*img(x,0,0));
    // }

    // void testCombine2InnerExpand()
    // {
        // Image3D res(img.shape());
        
        // View3D view = img.subarray(Size3(0,0,0), Size3(1,1,3));
        // combineTwoMultiArrays(view, img, res,
                              // Arg1() + Param(2.0)*Arg2());       
        // int x,y,z;
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), 2.0*img(x,y,z) + img(0,0,z));

        // combineTwoMultiArrays(img, view, res,
                              // Arg1() + Param(2.0)*Arg2());       
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), img(x,y,z) + 2.0*img(0,0,z));

        // combineTwoMultiArrays(srcMultiArrayRange(view), srcMultiArrayRange(view), 
                              // destMultiArrayRange(res),
                              // Arg1() + Param(2.0)*Arg2());       
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // shouldEqual(res(x,y,z), 3.0*img(0,0,z));
    // }

    // void testCombine2OuterReduce()
    // {
        // Image3D res(Size3(5,1,1));
        
        // combineTwoMultiArrays(img, img, res,
                              // reduceFunctor(Arg1() + Arg2() + Arg3(), 0.0));
        
        // int x,y,z;
        // for(x=0; x<img.shape(0); ++x)
        // {
            // double sum = 0.0;
            // for(y=0; y<img.shape(1); ++y)
                // for(z=0; z<img.shape(2); ++z)
                    // sum += img(x,y,z);
            // shouldEqual(res(x,0,0), 2.0*sum);
        // }
    // }

    // void testCombine2InnerReduce()
    // {
        // Image3D res(Size3(1,1,3));
        
        // combineTwoMultiArrays(img, img, res,
                              // reduceFunctor(Arg1() + Arg2() + Arg3(), 0.0));
        
        // int x,y,z;
        // for(z=0; z<img.shape(2); ++z)
        // {
            // double sum = 0.0;
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                    // sum += img(x,y,z);
            // shouldEqual(res(0,0,z), 2.0*sum);
        // }
    // }

    // void testCombine3()
    // {
        // Image3D res(img.shape()), res1(img.shape());
        
        // combineThreeMultiArrays(srcMultiArrayRange(img), 
                                // srcMultiArray(img), srcMultiArray(img), 
                                // destMultiArray(res),
                                // Arg1() + Arg2() + Arg3());
        // combineThreeMultiArrays(img, img, img, res1,
                                // Arg1() + Arg2() + Arg3());

        // int x,y,z;
        // for(z=0; z<img.shape(2); ++z)
            // for(y=0; y<img.shape(1); ++y)
                // for(x=0; x<img.shape(0); ++x)
                // {
                    // shouldEqual(res(x,y,z), 3.0*img(x,y,z));
                    // shouldEqual(res1(x,y,z), 3.0*img(x,y,z));
                // }
    // }
    
    // void testInitMultiArrayBorder(){
        // typedef vigra::MultiArray<1,int> IntLine;
        // typedef vigra::MultiArray<2,int> IntImage;
        // typedef vigra::MultiArray<3,int> IntVolume;
        
        // const int desired_vol[] ={  0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,

                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,

                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 5, 5, 0, 0,
                                    // 0, 0, 5, 5, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,

                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 5, 5, 0, 0,
                                    // 0, 0, 5, 5, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,

                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,

                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0,
                                    // 0, 0, 0, 0, 0, 0};

        // const int desired_img[] ={  0, 0, 0, 0, 0, 0,
                                    // 0, 5, 5, 5, 5, 0,
                                    // 0, 5, 5, 5, 5, 0,
                                    // 0, 5, 5, 5, 5, 0,
                                    // 0, 5, 5, 5, 5, 0,
                                    // 0, 0, 0, 0, 0, 0};

        // const int desired_lin[] ={  0, 0, 0, 5, 0, 0, 0 };

        // const int desired_vol2[] ={  0, 0,
                                     // 0, 0,

                                     // 0, 0, 
                                     // 0, 0};

        // IntVolume vol(IntVolume::difference_type(6,6,6));
        
        // for(IntVolume::iterator iter=vol.begin(); iter!=vol.end(); ++iter)
            // *iter=5;
        // initMultiArrayBorder(destMultiArrayRange(vol),2,0);
        // shouldEqualSequence(vol.begin(), vol.end(), desired_vol);

        // IntImage img(IntImage::difference_type(6,6));
        
        // for(IntImage::iterator iter=img.begin(); iter!=img.end(); ++iter)
            // *iter=5;
        // initMultiArrayBorder(destMultiArrayRange(img),1,0);
        // shouldEqualSequence(img.begin(), img.end(), desired_img);

        // IntLine lin(IntLine::difference_type(7));
        
        // for(IntLine::iterator iter=lin.begin(); iter!=lin.end(); ++iter)
            // *iter=5;
        // initMultiArrayBorder(destMultiArrayRange(lin),3,0);
        // shouldEqualSequence(lin.begin(), lin.end(), desired_lin);

        // IntVolume vol2(IntVolume::difference_type(2,2,2));
        
        // for(IntVolume::iterator iter=vol2.begin(); iter!=vol2.end(); ++iter)
            // *iter=5;
        // initMultiArrayBorder(vol2, 9, 0);
        // shouldEqualSequence(vol2.begin(), vol2.end(), desired_vol2);

    // }

    // void testInspect()
    // {
        // vigra::FindMinMax<PixelType> minmax;

        // inspectMultiArray(img, minmax);

        // shouldEqual(minmax.count, img.size());
        // shouldEqual(minmax.min, 0.1f);
        // shouldEqual(minmax.max, 59.1f);

        // vigra::MultiArray<3, unsigned char> labels(img.shape());
        // labels.subarray(Shape3(1,0,0), img.shape()-Shape3(1,0,0)) = 1;

        // vigra::ArrayOfRegionStatistics<vigra::FindMinMax<PixelType> > stats(1);

        // inspectTwoMultiArrays(img, labels, stats);

        // shouldEqual(stats[0].count, 24);
        // shouldEqual(stats[0].min, 0.1f);
        // shouldEqual(stats[0].max, 59.1f);
        // shouldEqual(stats[1].count, 36);
        // shouldEqual(stats[1].min, 1.1f);
        // shouldEqual(stats[1].max, 58.1f);
    // }
    
    // void testTensorUtilities()
    // {
        // MultiArrayShape<2>::type shape(3,4);
        // int size = shape[0]*shape[1];
        
        // MultiArray<2, TinyVector<double, 2> > vector(shape), rvector(shape);
        // MultiArray<2, TinyVector<double, 3> > tensor1(shape), tensor2(shape), rtensor(shape);
        // MultiArray<2, double > trace(shape), rtrace(shape);
        // MultiArray<2, double > determinant(shape), rdet(shape);
        
        // for(int k=0; k<size; ++k)
        // {
            // for(int l=0; l<2; ++l)
                // vector[k][l] = randomMT19937().uniform();
            // for(int l=0; l<3; ++l)
                // tensor1[k][l] = randomMT19937().uniform();
            // rdet[k] = tensor1[k][0]*tensor1[k][2] - sq(tensor1[k][1]);
        // }
        
        // vectorToTensor(srcImageRange(vector), destImage(rtensor));
        // vectorToTensorMultiArray(srcMultiArrayRange(vector), destMultiArray(tensor2));
        // shouldEqualSequence(tensor2.data(), tensor2.data()+size, rtensor.data());
        // tensor2.init(TinyVector<double, 3>());
        // vectorToTensorMultiArray(vector, tensor2);
        // shouldEqualSequence(tensor2.data(), tensor2.data()+size, rtensor.data());
                
        // tensorTrace(srcImageRange(tensor1), destImage(rtrace));
        // tensorTraceMultiArray(srcMultiArrayRange(tensor1), destMultiArray(trace));
        // shouldEqualSequence(trace.data(), trace.data()+size, rtrace.data());
        // trace = 0;
        // tensorTraceMultiArray(tensor1, trace);
        // shouldEqualSequence(trace.data(), trace.data()+size, rtrace.data());
                
        // tensorDeterminantMultiArray(srcMultiArrayRange(tensor1), destMultiArray(determinant));
        // shouldEqualSequence(determinant.data(), determinant.data()+size, rdet.data());
        // determinant = 0;
        // tensorDeterminantMultiArray(tensor1, determinant);
        // shouldEqualSequence(determinant.data(), determinant.data()+size, rdet.data());
                
        // determinant = 1000.0;
        // tensorDeterminantMultiArray(srcMultiArrayRange(tensor2), destMultiArray(determinant));
        // shouldEqualTolerance(norm(determinant), 0.0, 1e-14);

        // tensorEigenRepresentation(srcImageRange(tensor1), destImage(rtensor));
        // tensorEigenvaluesMultiArray(srcMultiArrayRange(tensor1), destMultiArray(vector));
        // shouldEqualSequenceTolerance(vector.begin(), vector.end(), rtensor.begin(), (TinyVector<double, 2>(1e-14)));

        // vector = TinyVector<double, 2>();
        // tensorEigenvaluesMultiArray(tensor1, vector);
        // shouldEqualSequenceTolerance(vector.begin(), vector.end(), rtensor.begin(), (TinyVector<double, 2>(1e-14)));
    // }
// };

template <class Array>
class ChunkedMultiArraySpeedTest
{
public:

    typedef typename Array::value_type T;
    typedef ChunkedArray<3, T> BaseArray;
    typedef VIGRA_UNIQUE_PTR<BaseArray> ArrayPtr;
    typedef typename BaseArray::iterator Iterator;
    
    Shape3 shape;
    ArrayPtr array;

    ChunkedMultiArraySpeedTest ()
        : shape(200, 201, 202)
    {
        array = createArray(shape, (Array *)0);
        linearSequence(array->begin(), array->end());
        std::cerr << "chunked multi array test for type " << typeid(Array).name() << ": \n";        
    }
    
    static ArrayPtr createArray(Shape3 const & shape, 
                                ChunkedArrayFull<3, T> *)
    {
        return ArrayPtr(new ChunkedArrayFull<3, T>(shape));
    }
    
    static ArrayPtr createArray(Shape3 const & shape, 
                                ChunkedArrayLazy<3, T> *)
    {
        return ArrayPtr(new ChunkedArrayLazy<3, T>(shape));
    }
    
    static ArrayPtr createArray(Shape3 const & shape, 
                                ChunkedArrayCompressed<3, T> *)
    {
        return ArrayPtr(new ChunkedArrayCompressed<3, T>(LZ4, shape));
    }
    
#ifdef HasHDF5
    static ArrayPtr createArray(Shape3 const & shape, 
                                ChunkedArrayHDF5<3, T> *)
    {
        HDF5File hdf5_file("chunked_test.h5", HDF5File::New);
        return ArrayPtr(new ChunkedArrayHDF5<3, T>(hdf5_file, "test", HDF5File::New, 
                                                   ZLIB_NONE, shape));
    }
#endif
    
    static ArrayPtr createArray(Shape3 const & shape, 
                                ChunkedArrayTmpFile<3, T> *)
    {
        return ArrayPtr(new ChunkedArrayTmpFile<3, T>(shape));
    }
    
    void testBaselineSpeed()
    {
        std::cerr << "############ chunked iterator speed #############\n";
        ChunkedArrayFull<3, T> a(shape);
        linearSequence(a.begin(), a.end());

        typename ChunkedArrayFull<3, T>::iterator i   = a.begin(),
                                                  end = a.end();
        USETICTOC;
        T count = 0;
        TIC;
        for(; i != end; ++i, ++count)
            if(count != *i)
            {
                shouldEqual(*i, count);
            }
        std::string t = TOCS;
        std::cerr << "    baseline:  " << t << "\n";
    }
    
    void testIteratorSpeed()
    {
        Iterator i   = array->begin(),
                 end = array->end();
        USETICTOC;
        T count = 0;
        TIC;
        for(; i != end; ++i, ++count)
        {
            if(count != *i)
            {
                shouldEqual(*i, count);
            }
        }
        std::string t = TOCS;
        std::cerr << "    read time: " << t << " (cache: " << array->cacheSize() << ")\n";
    }
    
    void testNestedLoopSpeed()
    {
        Iterator i   = array->begin(),
                 end = array->end();
        USETICTOC;
        T count = 0;
        TIC;
        for(i.setDim(2,0); i.coord(2) < shape[2]; i.incDim(2))
            for(i.setDim(1,0); i.coord(1) < shape[1]; i.incDim(1))
                for(i.setDim(0,0); i.coord(0) < shape[0]; i.incDim(0), ++count)
                {
                    if(count != *i)
                    {
                        shouldEqual(*i, count);
                    }
                }
        std::string t = TOCS;
        std::cerr << "    loop time: " << t << " (cache: " << array->cacheSize() << ")\n";
    }
    
    void testIteratorSpeed_LargeCache()
    {
        array = createArray(shape, (Array *)0);
        array->setCacheMaxSize(prod(array->chunkArrayShape()));
        linearSequence(array->begin(), array->end());
        testIteratorSpeed();
    }
    
    void testIndexingBaselineSpeed()
    {
        std::cerr << "################## indexing speed ####################\n";
        ChunkedArrayFull<3, T> a(shape);
        linearSequence(a.begin(), a.end());

        MultiCoordinateIterator<3> i(shape),
                                   end = i.getEndIterator();
        USETICTOC;
        T count = 0;
        TIC;
        for(; i != end; ++i, ++count)
            if(count != a[*i])
            {
                shouldEqual(a[*i], count);
            }
        std::string t = TOCS;
        std::cerr << "    baseline:  " << t << "\n";
    }
    
    void testIndexingSpeed()
    {
        MultiArrayView<3, T, ChunkedArrayTag> sub(array->subarray(Shape3(), shape));
        MultiCoordinateIterator<3> i(shape),
                                   end = i.getEndIterator();
        USETICTOC;
        T count = 0;
        TIC;
        for(; i != end; ++i, ++count)
        {
            if(count != sub[*i])
            {
                shouldEqual(sub[*i], count);
            }
        }
        std::string t = TOCS;
        std::cerr << "    indexing:  " << t << " (cache: " << array->cacheSize() << ")\n";
    }
};

struct ChunkedMultiArrayTestSuite
: public vigra::test_suite
{
    template <class Array>
    void testImpl()
    {
        add( testCase( &ChunkedMultiArrayTest<Array>::test_construction ) );
        add( testCase( &ChunkedMultiArrayTest<Array>::test_assignment ) );
        add( testCase( &ChunkedMultiArrayTest<Array>::test_bindAt ) );
        add( testCase( &ChunkedMultiArrayTest<Array>::test_bindInner ) );
        add( testCase( &ChunkedMultiArrayTest<Array>::test_bindOuter ) );
        add( testCase( &ChunkedMultiArrayTest<Array>::test_subarray ) );
        add( testCase( &ChunkedMultiArrayTest<Array>::test_iterator ) );
        add( testCase( &ChunkedMultiArrayTest<Array>::testMultiThreaded ) );
    }
    
    template <class T>
    void testSpeedImpl()
    {
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayFull<3, T> >::testBaselineSpeed )));
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayFull<3, T> >::testIteratorSpeed )));
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayLazy<3, T> >::testNestedLoopSpeed )));
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayLazy<3, T> >::testIteratorSpeed )));
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayCompressed<3, T> >::testIteratorSpeed )));
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayTmpFile<3, T> >::testIteratorSpeed )));
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayTmpFile<3, T> >::testIteratorSpeed_LargeCache )));
#ifdef HasHDF5
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayHDF5<3, T> >::testIteratorSpeed )));
#endif
    }
    
    template <class T>
    void testIndexingSpeedImpl()
    {
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayFull<3, T> >::testIndexingBaselineSpeed )));
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayFull<3, T> >::testIndexingSpeed )));
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayLazy<3, T> >::testIndexingSpeed )));
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayCompressed<3, T> >::testIndexingSpeed )));
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayTmpFile<3, T> >::testIndexingSpeed )));
#ifdef HasHDF5
        add( testCase( (&ChunkedMultiArraySpeedTest<ChunkedArrayHDF5<3, T> >::testIndexingSpeed )));
#endif
    }
    
    ChunkedMultiArrayTestSuite()
    : vigra::test_suite("ChunkedMultiArrayTestSuite")
    {
        testImpl<ChunkedArrayFull<3, float> >();
        testImpl<ChunkedArrayLazy<3, float> >();
        testImpl<ChunkedArrayCompressed<3, float> >();
        testImpl<ChunkedArrayTmpFile<3, float> >();
#ifdef HasHDF5
        testImpl<ChunkedArrayHDF5<3, float> >();
#endif
        
        testImpl<ChunkedArrayFull<3, TinyVector<float, 3> > >();
        testImpl<ChunkedArrayLazy<3, TinyVector<float, 3> > >();
        testImpl<ChunkedArrayCompressed<3, TinyVector<float, 3> > >();
        testImpl<ChunkedArrayTmpFile<3, TinyVector<float, 3> > >();
#ifdef HasHDF5
        testImpl<ChunkedArrayHDF5<3, TinyVector<float, 3> > >();
#endif
        
        testSpeedImpl<unsigned char>();
        testSpeedImpl<float>();
        testSpeedImpl<double>();
        
        testIndexingSpeedImpl<unsigned char>();
        testIndexingSpeedImpl<float>();
        testIndexingSpeedImpl<double>();
        
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
    }
};


int main(int argc, char ** argv)
{
    int failed = 0;

    ChunkedMultiArrayTestSuite test0;
    failed += test0.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test0.report() << std::endl;
    
    return (failed != 0);
}


