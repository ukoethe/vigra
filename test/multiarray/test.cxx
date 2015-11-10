/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
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

#include "vigra/unittest.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_iterator_coupled.hxx"
#include "vigra/multi_hierarchical_iterator.hxx"
#include "vigra/multi_impex.hxx"
#include "vigra/basicimageview.hxx"
#include "vigra/navigator.hxx"
#include "vigra/multi_pointoperators.hxx"
#include "vigra/tensorutilities.hxx"
#include "vigra/multi_tensorutilities.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/multi_math.hxx"
#include "vigra/algorithm.hxx"
#include "vigra/random.hxx"
#include "vigra/timing.hxx"
//#include "marray.hxx"

using namespace vigra;
using namespace vigra::functor;

template <class T>
class MultiArrayDataTest
{
public:

    // test a three-dimensional image of size 10*10*10.
    // the image is filled from start to beginning (memory-wise) with
    // ascending numbers.

    typedef typename vigra::detail::ResolveMultiband<T>::type   scalar_type;
    typedef MultiArray <2, scalar_type> array2_type;
    typedef MultiArray <3, T> array3_type;
    typedef MultiArrayView<3, Multiband<scalar_type> >  MultibandView3;
    typedef typename array3_type::view_type       array3_view_type;
    typedef typename array3_type::actual_stride   array3_stride;
    typedef typename array2_type::difference_type difference2_type;
    typedef typename array3_type::difference_type difference3_type;
    
    difference3_type shape3;
    array3_type array3;

    MultiArrayDataTest ()
        : shape3 (10, 10, 10), array3 (shape3, 1)
    {
        // initialize the array to the test data
        for (int i = 0; i < 1000; ++i)
            array3[i] = i;
    }

    void testHasData ()
    {
        should(array3.hasData());
        array3_type empty;
        should(!empty.hasData());
        array3_type empty2(Shape3(2,1,0));
        should(!empty2.hasData());
    }

    void testEquality ()
    {
        typedef difference3_type Shape;
        should(array3 == array3);
        should(array3 != array3.subarray(Shape(1,1,1), Shape(2,2,2)));
        should(array3.subarray(Shape(0,0,0), Shape(10,1,1)) != array3.subarray(Shape(0,1,0), Shape(10,2,1)));

        array3_type a(array3.shape());
        linearSequence(a.begin(), a.end());
        should(a == array3);

        for(unsigned int k=0; k<10; ++k)
            array3(k,0,0) += 10;
        should(array3.subarray(Shape(0,0,0), Shape(10,1,1)) == array3.subarray(Shape(0,1,0), Shape(10,2,1)));

        MultibandView3 channel_view(a.multiband());
        shouldEqual(a.shape(), channel_view.shape());
        shouldEqual(a.data(), channel_view.data());
    }
    
    // bindInner tests
    void test_bindInner ()
    {
        TinyVector <int, 2> inner_indices (2, 5);
        MultiArrayView <1, scalar_type, StridedArrayTag>
            array = array3.bindInner(inner_indices);
        shouldEqual ((array [TinyVector <int, 1> (0)]), 52);
        shouldEqual ((array [TinyVector <int, 1> (1)]), 152);
        shouldEqual ((array (0)), 52);
        shouldEqual ((array (1)), 152);
    }
    
    // bindOuter tests
    void test_bindOuter ()
    {
        TinyVector <int, 2> outer_indices (2, 5);
        MultiArrayView <1, scalar_type, array3_stride>
            array = array3.bindOuter(outer_indices);
        shouldEqual ((array [TinyVector <int, 1> (0)]), 520);
        shouldEqual ((array [TinyVector <int, 1> (1)]), 521);
        shouldEqual ((array (0)), 520);
        shouldEqual ((array (1)), 521);
    }

    // bindAt tests
    void test_bindAt ()
    {
        MultiArrayView <2, scalar_type, StridedArrayTag>
            array = array3.bindAt (1, 4);
        shouldEqual ((array [TinyVector <int, 2> (0, 0)]), 40);
        shouldEqual ((array [TinyVector <int, 2> (1, 0)]), 41);
        shouldEqual ((array [TinyVector <int, 2> (0, 1)]), 140);
        shouldEqual ((array (0, 0)), 40);
        shouldEqual ((array (1, 0)), 41);
        shouldEqual ((array (0, 1)), 140);
    }
    
    // bind tests
    void test_bind ()
    {
        MultiArrayView <2, scalar_type, array3_stride>
            array = array3.template bind<1>(4);
        shouldEqual ((array [TinyVector <int, 2> (0, 0)]), 40);
        shouldEqual ((array [TinyVector <int, 2> (1, 0)]), 41);
        shouldEqual ((array [TinyVector <int, 2> (0, 1)]), 140);
        shouldEqual ((array (0, 0)), 40);
        shouldEqual ((array (1, 0)), 41);
        shouldEqual ((array (0, 1)), 140);
    }

    void test_bind0 ()
    {
        MultiArrayView <2, scalar_type, StridedArrayTag>
            array = array3.template bind <0> (4);
        shouldEqual ((array [TinyVector <int, 2> (0, 0)]), 4);
        shouldEqual ((array [TinyVector <int, 2> (1, 0)]), 14);
        shouldEqual ((array [TinyVector <int, 2> (0, 1)]), 104);
        shouldEqual ((array (0, 0)), 4);
        shouldEqual ((array (1, 0)), 14);
        shouldEqual ((array (0, 1)), 104);
    }

    void test_singletonDimension ()
    {
        MultiArrayView <4, scalar_type, UnstridedArrayTag> a0 = array3.insertSingletonDimension(0);
        shouldEqual ((a0 [TinyVector <int, 4> (0, 4, 0, 0)]), 4);
        shouldEqual ((a0 [TinyVector <int, 4> (0, 4, 1, 0)]), 14);
        shouldEqual ((a0 [TinyVector <int, 4> (0, 4, 0, 1)]), 104);

        MultiArrayView <4, scalar_type, array3_stride> a1 = array3.insertSingletonDimension(1);
        shouldEqual ((a1 [TinyVector <int, 4> (4, 0, 0, 0)]), 4);
        shouldEqual ((a1 [TinyVector <int, 4> (4, 0, 1, 0)]), 14);
        shouldEqual ((a1 [TinyVector <int, 4> (4, 0, 0, 1)]), 104);

        MultiArrayView <4, scalar_type, array3_stride> a3 = array3.insertSingletonDimension(3);
        shouldEqual ((a3 [TinyVector <int, 4> (4, 0, 0, 0)]), 4);
        shouldEqual ((a3 [TinyVector <int, 4> (4, 1, 0, 0)]), 14);
        shouldEqual ((a3 [TinyVector <int, 4> (4, 0, 1, 0)]), 104);
    }

    // subarray and diagonal tests
    void test_subarray ()
    {
        MultiArray<1, scalar_type> diagRef(Shape1(10));
        linearSequence(diagRef.begin(), diagRef.end(), 0, 111);

        MultiArrayView <1, scalar_type, StridedArrayTag> diag = array3.diagonal();
        shouldEqual(diag.shape(0), 10);
        shouldEqualSequence(diagRef.begin(), diagRef.end(), diag.begin());

        typedef difference3_type Shape;
        
        Shape offset (1,1,1);
        Shape size (5,5,5);
        MultiArrayView <3, scalar_type> array = array3.subarray (offset, size);
        shouldEqual (array [Shape (0,0,0)], 111);
        shouldEqual (array [Shape (5,2,1)], 236);
        shouldEqual (array (0,0,0), 111);
        shouldEqual (array (5,2,1), 236);

        shouldEqual(array.shape(), array3.subarray(offset, Shape3(-5)).shape());
        shouldEqualSequence(array.begin(), array.end(), array3.subarray(offset, Shape3(-5)).begin());

        shouldEqual(array.shape(), array3.subarray(Shape3(-9), size).shape());
        shouldEqualSequence(array.begin(), array.end(), array3.subarray(offset, size).begin());

        shouldEqual(array.shape(), array3.subarray(Shape3(-9), Shape3(-5)).shape());
        shouldEqualSequence(array.begin(), array.end(), array3.subarray(Shape3(-9), Shape3(-5)).begin());
        
        // test swap
        array3.subarray(Shape(0,0,0), Shape(10,10,1)).swapData( 
             array3.subarray(Shape(0,0,1), Shape(10,10,2)));
        for(int k=0; k<100; ++k)
            shouldEqual(array3[k], k+100);
        for(int k=100; k<200; ++k)
            shouldEqual(array3[k], k-100);
            
        MultiArrayView <3, scalar_type> varray1 = array;
        MultiArrayView <3, scalar_type> varray2 = array3.subarray (Shape(0,1,0), Shape(3,4,5));
        shouldEqual (varray1.shape(), Shape(4,4,4));
        shouldEqual (varray2.shape(), Shape(3,3,5));
        varray2.swap(varray1);
        shouldEqual (varray2.shape(), Shape(4,4,4));
        shouldEqual (varray1.shape(), Shape(3,3,5));
        shouldEqual (varray2.data(), array.data());
        shouldEqual (varray1.data(), &array3(0,1,0));
    }
        
    // stridearray tests
    void test_stridearray ()
    {
        difference3_type stride (2, 2, 2);
        difference3_type traverser (0, 0, 0);
        MultiArrayView <3, scalar_type, StridedArrayTag>
            array = array3.stridearray (stride);
        shouldEqual (array [traverser], 0);
        traverser [0] += 1;
        shouldEqual (array [traverser], 2);
        traverser [1] += 1;
        shouldEqual (array [traverser], 22);
        traverser [2] += 1;
        shouldEqual (array [traverser], 222);
        shouldEqual (array (1,1,1), 222);
    }

    void testIsUnstrided()
    {
        typedef difference3_type Shape;

        should(array3.isUnstrided());
        should(array3.isUnstrided(0));
        should(array3.isUnstrided(1));
        should(array3.isUnstrided(2));
        should(array3.bindOuter(0).isUnstrided());
        should(!array3.bindInner(0).isUnstrided());
        should(!array3.bindAt(1, 0).isUnstrided());
        should(array3.bindAt(1, 0).isUnstrided(0));
        should(!array3.subarray(Shape(), array3.shape()-Shape(1)).isUnstrided());
        should(!array3.subarray(Shape(), array3.shape()-Shape(1)).isUnstrided(1));
        should(array3.subarray(Shape(), array3.shape()-Shape(1)).isUnstrided(0));
        should(!array3.subarray(Shape(), array3.shape()-Shape(0,2,2)).isUnstrided());
        should(array3.subarray(Shape(), array3.shape()-Shape(0,2,2)).isUnstrided(1));
        should(array3.subarray(Shape(), array3.shape()-Shape(0,2,2)).isUnstrided(0));
    }

    void testIsStrided()
    {
        // for MultiArray<3, Multiband<T> >

        should(!array3.isUnstrided());
        should(array3.permuteStridesAscending().isUnstrided());
        should(!array3.isUnstrided(0));
        should(!array3.isUnstrided(1));
        should(!array3.isUnstrided(2));
        should(array3.bindInner(Shape2(0,0)).isUnstrided());
    }

    // permute and transpose tests
    void testPermute ()
    {
        typedef MultiArrayView <3, scalar_type, StridedArrayTag> transposed_view;
        array3.reshape(difference3_type(3,5,7));
        for(int k=0; k<array3.size(); ++k)
            array3[k] = k;

        array3_type ref(difference3_type(array3.shape(2), array3.shape(0), array3.shape(1)));
        for(int k=0; k<array3.shape(0); ++k)
            for(int l=0; l<array3.shape(1); ++l)
                for(int m=0; m<array3.shape(2); ++m)
                    ref(m,k,l) = array3(k,l,m);

        MultiArrayView <3, scalar_type, StridedArrayTag>
                parray = array3.transpose (difference3_type (2, 0, 1));        
        shouldEqual(ref.shape(), parray.shape());
        should(ref == parray);
        
        if(vigra::detail::ResolveMultiband<T>::value) // array is Multiband<T>
        {
            shouldEqual(array3.strideOrdering(), Shape3(1,2,0));
            
            array3_type ref_ascending(difference3_type(array3.shape(2), array3.shape(0), array3.shape(1)));
            for(int k=0; k<array3.shape(0); ++k)
                for(int l=0; l<array3.shape(1); ++l)
                    for(int m=0; m<array3.shape(2); ++m)
                        ref_ascending(m,k,l) = array3(k,l,m);

            transposed_view parray_ascending = array3.permuteStridesAscending();
            shouldEqual(ref_ascending.shape(), parray_ascending.shape());
            should(ref_ascending == parray_ascending);
            
            array3_type ref_descending(difference3_type(array3.shape(1), array3.shape(0), array3.shape(2)));
            for(int k=0; k<array3.shape(0); ++k)
                for(int l=0; l<array3.shape(1); ++l)
                    for(int m=0; m<array3.shape(2); ++m)
                        ref_descending(l,k,m) = array3(k,l,m);

            transposed_view parray_descending = array3.permuteStridesDescending();
            shouldEqual(ref_descending.shape(), parray_descending.shape());
            should(ref_descending == parray_descending);
        }
        else
        {
            shouldEqual(array3.strideOrdering(), Shape3(0,1,2));
            transposed_view parray_ascending = parray.permuteStridesAscending();        
            shouldEqual(array3.shape(), parray_ascending.shape());
            should(array3 == parray_ascending);

            array3_type ref_descending(difference3_type(array3.shape(2), array3.shape(1), array3.shape(0)));
            for(int k=0; k<array3.shape(0); ++k)
                for(int l=0; l<array3.shape(1); ++l)
                    for(int m=0; m<array3.shape(2); ++m)
                        ref_descending(m,l,k) = array3(k,l,m);

            MultiArrayView <3, scalar_type, StridedArrayTag>
                    parray_descending = array3.permuteStridesDescending();        
            shouldEqual(ref_descending.shape(), parray_descending.shape());
            should(ref_descending == parray_descending);
        }

        array2_type ref2(difference2_type(array3.shape(1), array3.shape(0)));
        for(int k=0; k<array3.shape(0); ++k)
            for(int l=0; l<array3.shape(1); ++l)
                ref2(l, k) = array3(k, l, 0);

        MultiArrayView <2, scalar_type, StridedArrayTag>
                array2 = array3.bindOuter(0).transpose ();
        shouldEqual(ref2.shape(), array2.shape());
        should(ref2 == array2);

        try {
            array3.transpose (difference3_type (2, 0, 0));   
            failTest("no exception thrown");
        }
        catch(vigra::ContractViolation & c)
        {
            std::string expected("\nPrecondition violation!\nMultiArrayView::transpose(): every dimension must occur exactly once");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
    }

    void testMethods ()
    {
        shouldEqual(array3.squaredNorm(), 332833500);
        
        shouldEqual(array3.norm(), std::sqrt(332833500.0));
        shouldEqual(array3.norm(0), 999.0);
        shouldEqual(array3.norm(1), 499500.0);
        shouldEqualTolerance(array3.norm(2, false), std::sqrt(332833500.0), 1e-14);
        
        difference3_type first(0,0,0), last(1,1,1);
        shouldEqual(array3.subarray(first, last).norm(), 0.0);
        shouldEqual(array3.subarray(first, last).norm(0), 0.0);
        shouldEqual(array3.subarray(first, last).norm(1), 0.0);
        shouldEqual(array3.subarray(first, last).norm(2, false), 0.0);

        shouldEqual(array3.squaredNorm(), squaredNorm(array3));
        shouldEqual(array3.norm(), vigra::norm(array3));

        should(array3.any());
        should(!array3.subarray(first, last).any());
        should(!array3.all());
        should(array3.subarray(last, array3.shape()).all());

        shouldEqual(array3.template sum<int>(), 499500);
        shouldEqual(array3.subarray(Shape3(1,1,1),Shape3(3,3,2)).template product<int>(), 183521184);

        Shape3 reducedShape(1, 1, array3.shape(2));
        array3_type reducedSums(reducedShape);
        array3.sum(reducedSums);
        int res = 4950;
        for(int k=0; k<reducedShape[2]; ++k, res += 10000)
            shouldEqual(reducedSums(0,0,k), res);

        scalar_type minimum, maximum;
        array3.minmax(&minimum, &maximum);
        shouldEqual(minimum, 0);
        shouldEqual(maximum, array3.size()-1);

        double mean, variance;
        array3.meanVariance(&mean, &variance);
        shouldEqual(mean, 499.5);
        shouldEqual(variance, 83333.25);
    }
    
    void testScanOrderAccess()
    {
        shouldEqual(array3.size(), 1000);
        for(int k=0; k< array3.size(); ++k)
        {
            shouldEqual(array3[k], k);
            shouldEqual(array3[array3.scanOrderIndexToCoordinate(k)], k);
            shouldEqual(array3.coordinateToScanOrderIndex(array3.scanOrderIndexToCoordinate(k)), k);
        }
            
        MultiArrayView <2, scalar_type, StridedArrayTag> array = array3.bindInner(2);
        shouldEqual(array.size(), 100);
        for(int k=0; k< array.size(); ++k)
            shouldEqual(array[k], 10*k+2);
            
        MultiArrayView <2, scalar_type, array3_stride>
            subarray = array3.bindOuter(2).subarray(Shape2(1,0), Shape2(10,9));
        shouldEqual(subarray.size(), 81);
        for(int k=0, l=200; k< subarray.size(); ++k, ++l)
        {
            if(k%9 == 0)
                ++l;
            shouldEqual(subarray[k], l);
        }
    }
    
    void testAssignmentAndReset()
    {
        typedef Shape3 Shape;
        typename array3_type::view_type array;
        array = array3;
        should(array3 == array);
        try {
            array = array3.subarray(Shape(0,0,0), Shape(10,1,1));
            failTest("no exception thrown");
        }
        catch(vigra::ContractViolation & c)
        {
            std::string expected("\nPrecondition violation!\nMultiArrayView::operator=(MultiArrayView const &): shape mismatch");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
	array.reset();
	array = array3.subarray(Shape(0,0,0), Shape(10,1,1)); // possible after reset.
        MultiArrayView <3, scalar_type, array3_stride> subarray = array3.subarray(Shape(0,0,0), Shape(10,1,1));
	should(subarray == array);
        subarray = array3.subarray(Shape(0,1,0), Shape(10,2,1)); // should overwrite the data
        for(unsigned int k=0; k<10; ++k)
            shouldEqual(array3(k,0,0), array3(k,1,0));
        subarray += array3.subarray(Shape(0,1,0), Shape(10,2,1)); // should overwrite the data
        for(unsigned int k=0; k<10; ++k)
            shouldEqual(array3(k,0,0), 2.0*array3(k,1,0));
        subarray -= array3.subarray(Shape(0,1,0), Shape(10,2,1)); // should overwrite the data
        for(unsigned int k=0; k<10; ++k)
            shouldEqual(array3(k,0,0), array3(k,1,0));
        subarray *= array3.subarray(Shape(0,1,0), Shape(10,2,1)); // should overwrite the data
        for(unsigned int k=0; k<10; ++k)
            shouldEqual(array3(k,0,0), sq(array3(k,1,0)));
        subarray /= array3.subarray(Shape(0,1,0), Shape(10,2,1)); // should overwrite the data
        for(unsigned int k=0; k<10; ++k)
            shouldEqual(array3(k,0,0), array3(k,1,0));
            
        subarray += 1; // should overwrite the data
        for(unsigned int k=0; k<10; ++k)
            shouldEqual(array3(k,0,0), array3(k,1,0) + 1);
        subarray -= 2; // should overwrite the data
        for(unsigned int k=0; k<10; ++k)
            shouldEqual(array3(k,0,0), array3(k,1,0) - 1);
        subarray *= 2; // should overwrite the data
        for(unsigned int k=0; k<10; ++k)
            shouldEqual(array3(k,0,0), 2*(array3(k,1,0) - 1));
        subarray /= 2; // should overwrite the data
        for(unsigned int k=0; k<10; ++k)
            shouldEqual(array3(k,0,0), array3(k,1,0) - 1);

        // test assignment from UInt8 => double (reproduces compiler crash in VisualStudio 2010)
        MultiArray<2, UInt8> d1(Shape2(50, 100));
        MultiArray<2, double> d2(d1.shape());
        d2 = d1;
    }
};


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

        array1_t b(Shape1(5), LinearSequence);
        shouldEqual (b.shape (0), 5);
        int ref[] = { 0, 1, 2, 3, 4 };
        shouldEqualSequence(b.begin(), b.end(), ref);
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

        shouldEqual(&i1[Shape3(0,0,0)], &a3(0,0,0));
        shouldEqual(&i1[Shape3(1,0,0)], &a3(1,0,0));
        shouldEqual(&i1[Shape3(0,1,0)], &a3(0,1,0));
        shouldEqual(&i1[Shape3(1,1,0)], &a3(1,1,0));
        shouldEqual(&i1[Shape3(0,0,1)], &a3(0,0,1));
        shouldEqual(&i1[Shape3(1,0,1)], &a3(1,0,1));
        shouldEqual(&i1[Shape3(1,1,1)], &a3(1,1,1));

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

        shouldEqual(&i3[-shape3_t(0,0,0)], &a3(1,2,4));
        shouldEqual(&i3[-shape3_t(1,0,0)], &a3(0,2,4));
        shouldEqual(&i3[-shape3_t(0,1,0)], &a3(1,1,4));
        shouldEqual(&i3[-shape3_t(1,1,0)], &a3(0,1,4));
        shouldEqual(&i3[-shape3_t(0,0,1)], &a3(1,2,3));
        shouldEqual(&i3[-shape3_t(1,0,1)], &a3(0,2,3));
        shouldEqual(&i3[-shape3_t(1,1,1)], &a3(0,1,3));
        
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

        shouldEqual(&get<1>(i1[shape3_t(0,0,0)]), &a3(0,0,0));
        shouldEqual(&get<1>(i1[shape3_t(1,0,0)]), &a3(1,0,0));
        shouldEqual(&get<1>(i1[shape3_t(0,1,0)]), &a3(0,1,0));
        shouldEqual(&get<1>(i1[shape3_t(1,1,0)]), &a3(1,1,0));
        shouldEqual(&get<1>(i1[shape3_t(0,0,1)]), &a3(0,0,1));
        shouldEqual(&get<1>(i1[shape3_t(1,0,1)]), &a3(1,0,1));
        shouldEqual(&get<1>(i1[shape3_t(1,1,1)]), &a3(1,1,1));

        shouldEqual(&get<1>(*(iend-1)), &a3(1,2,4));
        shouldEqual(&get<1>(*(iend-2)), &a3(0,2,4));
        shouldEqual(&get<1>(*(iend-3)), &a3(1,1,4));
        shouldEqual(&get<1>(*(iend-7)), &a3(1,2,3));
        shouldEqual(&get<1>(*(iend-8)), &a3(0,2,3));
        shouldEqual(&get<1>(*(iend-10)), &a3(0,1,3));

        should(cast<1>(*i1).arrayView() == a3);

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
                    shouldEqual(&get<1>(i1), &a3[p]);
                    shouldEqual(&get<1>(i2), &a3[p]);
                    shouldEqual(&i1.get<1>(), &a3[p]);
                    shouldEqual(&i2.get<1>(), &a3[p]);
                    //shouldEqual(i1.operator->(), &a3[p]);
                    //shouldEqual(i2.operator->(), &a3[p]);
                    shouldEqual(i1.point(), p);
                    shouldEqual(i2.point(), p);
                    shouldEqual(i1.get<0>(), p);
                    shouldEqual(i2.get<0>(), p);
                    shouldEqual(get<0>(i1), p);
                    shouldEqual(get<0>(i2), p);
                    shouldEqual(i1.scanOrderIndex(), count);
                    shouldEqual(i2.scanOrderIndex(), count);

                    should(cast<1>(*i1).arrayView() == a3);

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
        Iterator1 iz = zip(vi.begin(), vd.begin());

        count = 0;
        for(; count < 4; ++iz, ++it, ++i0, ++count)
        {
            shouldEqual(i0.get<0>(), Shape1(count));
            shouldEqual(it.get<0>(), Shape1(count));
            shouldEqual(it.get<1>(), count+10);
            shouldEqual(it.get<2>(), count+20.0);
            shouldEqual(iz.get<0>(), Shape1(count));
            shouldEqual(iz.get<1>(), count+10);
            shouldEqual(iz.get<2>(), count+20.0);
        }
        shouldEqual(it, end);
        shouldEqual(iz, end);

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
            should(cast<1>(*im).arrayView() == at);
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

    void test_hierarchical ()
    {
        // test hierarchical navigation and 
        auto i3_f = createHierarchicalIterator(a3),
             i3_l = i3_f.getEndIterator();
        auto i3_c = zip(i3_f, i3_f);
        array3_t::iterator seqi = a3.begin();

        int countx = 0, county = 0, countz = 0;

        // iterate over the third dimension
        for (int z=0; i3_f != i3_l; ++i3_f, ++i3_c, ++z) 
        {
            auto i2_f = i3_f.begin ();
            auto i2_c = i3_c.begin ();
            auto i2_l = i3_f.end ();
            // iterate over the second dimension
            for (int y=0; i2_f != i2_l; ++i2_f, ++i2_c, ++y) 
            {
                auto i1_f = i2_f.begin ();
                auto i1_c = i2_c.begin ();
                auto i1_l = i2_f.end ();
                // iterate over the first dimension
                for (int x=0; i1_f != i1_l; ++i1_f, ++i1_c, ++x, ++seqi)
                {
                    ++countx;
                    shouldEqual(&*seqi, &a3(x,y,z));
                    shouldEqual(i1_f.point(), Shape3(x,y,z));
                    shouldEqual(&*i1_f, &a3(x,y,z));
                    shouldEqual(i1_c->get<0>(), Shape3(x,y,z));
                    shouldEqual(&i1_c->get<1>(), &a3(x,y,z));
                    shouldEqual(&i1_c->get<2>(), &a3(x,y,z));
                    shouldEqual(&get<1>(*i1_c), &a3(x,y,z));
                    shouldEqual(&get<2>(*i1_c), &a3(x,y,z));
                }
                ++county;
            }
            ++countz;
        }

        shouldEqual (countx, 30);
        shouldEqual (county, 15);
        shouldEqual (countz, 5);
        shouldEqual (seqi, a3.end());
        //
        //// test direct navigation
        //traverser3_t i3 = a3.traverser_begin();
        //shouldEqual(&*i3, &a3[shape3_t(0,0,0)]);

        //i3.dim<2>()++;
        //i3.dim<1>()++;
        //i3.dim<0>()++;        
        //shouldEqual(&*i3, &a3[shape3_t(1,1,1)]);

        //i3.dim<2>()+= 3;
        //i3.dim<1>()+= 2;
        //i3.dim<0>()+= 1;        
        //shouldEqual(&*i3, &a3[shape3_t(2,3,4)]);
        //shouldEqual(&i3[shape3_t(-2,-3,-4)], &a3[shape3_t(0,0,0)]);
        //shouldEqual(&*(i3-shape3_t(2,3,4)), &a3[shape3_t(0,0,0)]);

        //i3.dim<2>()--;
        //i3.dim<1>()--;
        //i3.dim<0>()--;        
        //shouldEqual(&*i3, &a3[shape3_t(1,2,3)]);

        //i3.dim<2>()-= 3;
        //i3.dim<1>()-= 2;
        //i3.dim<0>()-= 1;        
        //shouldEqual(&*i3, &a3[shape3_t(0,0,0)]);

        //shouldEqual(&i3[shape3_t(2,3,4)], &a3[shape3_t(2,3,4)]);
        //shouldEqual(&*(i3+shape3_t(2,3,4)), &i3[shape3_t(2,3,4)]);
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

        shouldEqual (st.shape (0), 7);
        shouldEqual (st.shape (1), 8);
        shouldEqual (st.shape (2), 9);
        
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

        shouldEqual (countx, 504u);
        shouldEqual (county, 72u);
        shouldEqual (countz, 9u);
  
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

class MultiArrayNavigatorTest
{
public:

    typedef unsigned int scalar_type;
    typedef MultiArray <3, scalar_type> array3_type;
    typedef MultiArrayView <3, scalar_type> array3_view_type;
    typedef array3_type::difference_type difference3_type;
    
    difference3_type shape3;
    array3_type array3;

    MultiArrayNavigatorTest ()
        : shape3 (4, 3, 2), array3 (shape3)
    {
        // initialize the array to the test data
        for (unsigned int i = 0; i < 24; ++i)
            array3.data () [i] = i;
    }

    void testNavigator ()
    {
        unsigned char expected[][24] = 
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
            {0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11, 12, 16, 20, 13, 17, 21, 14, 18, 22, 15, 19, 23},
            {0, 12, 1, 13, 2, 14, 3, 15, 4, 16, 5, 17, 6, 18, 7, 19, 8, 20, 9, 21, 10, 22, 11, 23}};
        typedef MultiArrayNavigator<array3_type::traverser, 3> Navigator;
        for(int d=0; d<3; ++d)
        {
            Navigator nav(array3.traverser_begin(), array3.shape(), d);
            int k = 0;
            for(; nav.hasMore(); ++nav)
            {
                Navigator::iterator i = nav.begin(), end = nav.end();
                for(; i != end; ++i, ++k)
                    shouldEqual(*i, expected[d][k]);
            }
        }
        
        Shape3 start(1, 1, 0), stop(3, 3, 2);
        unsigned char sexpected[][8] = 
            {{5,  6,  9, 10, 17, 18, 21, 22},
            {5,  9,  6, 10, 17, 21, 18, 22},
            {5, 17,  6, 18,  9, 21, 10, 22}};
        for(int d=0; d<3; ++d)
        {
            Navigator nav(array3.traverser_begin(), start, stop, d);
            int k = 0;
            for(; nav.hasMore(); ++nav)
            {
                Navigator::iterator i = nav.begin(), end = nav.end();
                for(; i != end; ++i, ++k)
                    shouldEqual(*i, sexpected[d][k]);
            }
        }
    }

    void testCoordinateNavigator ()
    {
        unsigned char expected[][24] = 
            {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
            {0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11, 12, 16, 20, 13, 17, 21, 14, 18, 22, 15, 19, 23},
            {0, 12, 1, 13, 2, 14, 3, 15, 4, 16, 5, 17, 6, 18, 7, 19, 8, 20, 9, 21, 10, 22, 11, 23}};
        typedef MultiCoordinateNavigator<3> Navigator;
        for(int d=0; d<3; ++d)
        {
            Navigator nav(array3.shape(), d);
            int k = 0;
            for(; nav.hasMore(); ++nav)
            {
                Navigator::value_type i = nav.begin(), end = nav.end();
                for(; i[d] != end[d]; ++i[d], ++k)
                    shouldEqual(array3[i], expected[d][k]);
            }
        }
    }    
};

struct MultiImpexTest
{
    typedef MultiArray<3, unsigned char> Array;
    typedef Array::difference_type Shape;
    typedef Array::traverser Traverser;
    
    Array array;
    
    MultiImpexTest()
    : array(Shape(2,3,4))
    {
        int value = 1;
        
        Traverser i3 = array.traverser_begin ();

        for (; i3 != array.traverser_end(); ++i3, ++value) 
        {
            typedef Traverser::next_type Traverser2;
            Traverser2 i2 = i3.begin ();
            for (; i2 != i3.end(); ++i2) 
            {
                typedef Traverser2::next_type Traverser1;
                Traverser1 i1 = i2.begin ();
                for (; i1 != i2.end(); ++i1)
                {
                    *i1 = value;
                }
            }
        }
    }
    
    void testImpex()
    {
#if defined(HasTIFF)
        const char * ext1 = ".tif";
#else
        const char * ext1 = ".xv";
#endif
        exportVolume(array, VolumeExportInfo("test", ext1));
        
        Array result;
        
        VolumeImportInfo import_info("test", ext1);
        shouldEqual(Shape(2,3,4), import_info.shape());

        result.reshape(import_info.shape());
        importVolume(import_info, result);
        
        shouldEqual(result(0,1,0), 1);
        shouldEqual(result(0,1,1), 2);
        shouldEqual(result(0,1,2), 3);
        shouldEqual(result(0,1,3), 4);

#if defined(HasPNG)
        const char * ext2 = ".png";
#else
        const char * ext2 = ".pnm";
#endif
        exportVolume(array, std::string("impex/test"), std::string(ext2));
        
        importVolume(result, std::string("impex/test"), std::string(ext2));
        
        shouldEqual(result.shape(), Shape(2,3,4));
        shouldEqual(result(0,1,0), 1);
        shouldEqual(result(0,1,1), 2);
        shouldEqual(result(0,1,2), 3);
        shouldEqual(result(0,1,3), 4);

#ifdef _MSC_VER
        exportVolume(array, VolumeExportInfo("impex\\test", ext2));
        
        importVolume(result, std::string("impex\\test"), std::string(ext2));
        
        shouldEqual(result.shape(), Shape(2,3,4));
        shouldEqual(result(0,1,0), 1);
        shouldEqual(result(0,1,1), 2);
        shouldEqual(result(0,1,2), 3);
        shouldEqual(result(0,1,3), 4);
#endif // _MSC_VER
    }

#if defined(HasTIFF)
    void testMultipageTIFF()
    {
        exportVolume(array, VolumeExportInfo("multipage.tif"));

        VolumeImportInfo info("multipage.tif");
        shouldEqual(Shape(2,3,4), info.shape());

        Array result(info.shape());
        importVolume(info, result);
        shouldEqual(result(0,1,0), 1);
        shouldEqual(result(0,1,1), 2);
        shouldEqual(result(0,1,2), 3);
        shouldEqual(result(0,1,3), 4);
    }
#endif

};

template <class IMAGE>
struct ImageViewTest
{
    typedef typename IMAGE::value_type value_type;
    typedef MultiArray<2, value_type> MA;
    typedef BasicImageView<value_type> Image;
    static typename Image::value_type data[];

    ImageViewTest()
    : ma(TinyVector<int, 2>(3,3)),
      img(makeBasicImageView(ma))
    {
        typename Image::Accessor acc = img.accessor();
        typename Image::iterator i = img.begin();

        acc.set(data[0], i);
        ++i;
        acc.set(data[1], i);
        ++i;
        acc.set(data[2], i);
        ++i;
        acc.set(data[3], i);
        ++i;
        acc.set(data[4], i);
        ++i;
        acc.set(data[5], i);
        ++i;
        acc.set(data[6], i);
        ++i;
        acc.set(data[7], i);
        ++i;
        acc.set(data[8], i);
        ++i;
        should(i == img.end());
    }

    template <class Iterator>
    void scanImage(Iterator ul, Iterator lr)
    {
        Iterator y = ul;
        Iterator x = ul;
        typename Image::Accessor acc = img.accessor();

        should(acc(x) == data[0]);
        ++x.x;
        should(acc(x) == data[1]);
        ++x.x;
        should(acc(x) == data[2]);
        ++x.x;
        should(x.x == lr.x);

        ++y.y;
        x = y;
        should(acc(x) == data[3]);
        ++x.x;
        should(acc(x) == data[4]);
        ++x.x;
        should(acc(x) == data[5]);
        ++y.y;
        x = y;
        should(acc(x) == data[6]);
        ++x.x;
        should(acc(x) == data[7]);
        ++x.x;
        should(acc(x) == data[8]);
        ++y.y;
        should(y.y == lr.y);

        y = ul;
        should(acc(y, vigra::Diff2D(1,1)) == data[4]);
    }

    template <class Iterator>
    void scanRows(Iterator r1, Iterator r2, Iterator r3, int w)
    {
        Iterator end = r1 + w;
        typename Image::Accessor acc = img.accessor();

        should(acc(r1) == data[0]);
        ++r1;
        should(acc(r1) == data[1]);
        ++r1;
        should(acc(r1) == data[2]);
        ++r1;
        should(r1 == end);

        end = r2 + w;
        should(acc(r2) == data[3]);
        ++r2;
        should(acc(r2) == data[4]);
        ++r2;
        should(acc(r2) == data[5]);
        ++r2;
        should(r2 == end);

        end = r3 + w;
        should(acc(r3) == data[6]);
        ++r3;
        should(acc(r3) == data[7]);
        ++r3;
        should(acc(r3) == data[8]);
        ++r3;
        should(r3 == end);
    }

    template <class Iterator>
    void scanColumns(Iterator c1, Iterator c2, Iterator c3, int h)
    {
        Iterator end = c1 + h;
        typename Image::Accessor acc = img.accessor();

        should(acc(c1) == data[0]);
        ++c1;
        should(acc(c1) == data[3]);
        ++c1;
        should(acc(c1) == data[6]);
        ++c1;
        should(c1 == end);

        end = c2 + h;
        should(acc(c2) == data[1]);
        ++c2;
        should(acc(c2) == data[4]);
        ++c2;
        should(acc(c2) == data[7]);
        ++c2;
        should(c2 == end);

        end = c3 + h;
        should(acc(c3) == data[2]);
        ++c3;
        should(acc(c3) == data[5]);
        ++c3;
        should(acc(c3) == data[8]);
        ++c3;
        should(c3 == end);
    }

    void testStridedImageView()
    {
        // create stride MultiArrayView
        typename MA::difference_type
            start(0,0), end(2,2);
        MA roi = ma.subarray(start, end);

        // inspect both the MultiArrayView and the corresponding BasicImageView
        vigra::FindSum<typename Image::value_type> sum1, sum2;
        vigra::inspectMultiArray(srcMultiArrayRange(roi), sum1);
        vigra::inspectImage(srcImageRange(makeBasicImageView(roi)), sum2);

        shouldEqual(sum1.sum(), sum2.sum());
        shouldEqual(data[0] + data[1] + data[3] + data[4], sum2.sum());
    }

    void testBasicImageIterator()
    {
        typename Image::Iterator ul = img.upperLeft();
        typename Image::Iterator lr = img.lowerRight();

        scanImage(ul, lr);
        scanRows(ul.rowIterator(), (ul+Diff2D(0,1)).rowIterator(),
                 (ul+Diff2D(0,2)).rowIterator(), img.width());
        scanColumns(ul.columnIterator(), (ul+Diff2D(1,0)).columnIterator(),
                 (ul+Diff2D(2,0)).columnIterator(), img.height());

        typename Image::Accessor acc = img.accessor();
        typename Image::iterator i = img.begin();
        should(acc(i, 4) == data[4]);
    }

    void testImageIterator()
    {
        vigra::ImageIterator<typename Image::value_type>
            ul(img.begin(), img.width());
        vigra::ImageIterator<typename Image::value_type> lr = ul + img.size();
        scanImage(ul, lr);
        scanRows(ul.rowIterator(), (ul+Diff2D(0,1)).rowIterator(),
                 (ul+Diff2D(0,2)).rowIterator(), img.width());
        scanColumns(ul.columnIterator(), (ul+Diff2D(1,0)).columnIterator(),
                 (ul+Diff2D(2,0)).columnIterator(), img.height());
    }

    void copyImage()
    {
        typedef typename Image::value_type Value;

        Image img1(img);
        typename Image::iterator i = img.begin();
        typename Image::iterator i1 = img1.begin();
        typename Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i, ++i1)
        {
            should(acc(i) == acc(i1));
        }

        img.init(NumericTraits<Value>::zero());
        for(; i != img.end(); ++i)
        {
            should(acc(i) == Value(NumericTraits<Value>::zero()));
        }
        img(1,1) = Value(200);
        img1 = img;
        i = img.begin();
        i1 = img1.begin();

        for(; i != img.end(); ++i, ++i1)
        {
            should(acc(i) == acc(i1));
        }
    }

    MA ma;
    Image img;
};

typedef ImageViewTest<vigra::BImage> BImageViewTest;

template <>
unsigned char BImageViewTest::data[] = {1,2,3,4,5,6,7,8,9};

typedef ImageViewTest<vigra::DImage> DImageViewTest;

template <>
double DImageViewTest::data[] = {1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9};

typedef ImageViewTest<vigra::BRGBImage> BRGBImageViewTest;
typedef vigra::RGBValue<unsigned char> BRGB;
template <>
BRGB BRGBImageViewTest::data[] = {
    BRGB(1,1,1),
    BRGB(2,2,2),
    BRGB(3,3,3),
    BRGB(4,4,4),
    BRGB(5,5,5),
    BRGB(6,6,6),
    BRGB(7,7,7),
    BRGB(8,8,8),
    BRGB(9,9,9)
};

typedef ImageViewTest<vigra::FRGBImage> FRGBImageViewTest;
typedef vigra::RGBValue<float> FRGB;
template <>
FRGB FRGBImageViewTest::data[] = {
    FRGB(1.1f, 1.1f, 1.1f),
    FRGB(2.2f, 2.2f, 2.2f),
    FRGB(3.3f, 3.3f, 3.3f),
    FRGB(4.4f, 4.4f, 4.4f),
    FRGB(5.5f, 5.5f, 5.5f),
    FRGB(6.6f, 6.6f, 6.6f),
    FRGB(7.7f, 7.7f, 7.7f),
    FRGB(8.8f, 8.8f, 8.8f),
    FRGB(9.9f, 9.9f, 9.9f)
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

        const int desired_asym_img[] ={  0, 0, 0, 0, 0, 0,
                                         0, 0, 5, 5, 5, 0,
                                         0, 0, 5, 5, 5, 0,
                                         0, 0, 5, 5, 5, 0,
                                         0, 0, 0, 0, 0, 0,
                                         0, 0, 0, 0, 0, 0};
        img = 5;
        initMultiArrayBorder(img, Shape2(2,1), Shape2(1,2), 0);
        shouldEqualSequence(img.begin(), img.end(), desired_asym_img);
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

class MultiMathTest
{
public:

    typedef double scalar_type;
    typedef MultiArray <3, scalar_type> array3_type;
    typedef MultiArrayView <3, scalar_type> view_type;
    typedef array3_type::difference_type shape3_type;
    
    shape3_type shape3;
    array3_type r1, r2, r3, r4, r5, a, b, c, d;
    view_type rv, bv, cv, dv;

    MultiMathTest ()
        : shape3 (4, 3, 2), 
          r2(shape3),
          r3(shape3),
          r4(shape3),
          r5(Shape3(2,3,4)),
          a(shape3),
          b(shape3),
          c(shape3),
          d(shape3),
          rv(r3),
          bv(b),
          cv(c),
          dv(d)
    {
        // initialize the array to the test data
        for (unsigned int i = 0; i < 24; ++i)
        {
            a[i] = std::exp(-0.2*i);
            b[i] = 0.5 + i;
            c[i] = 100.0 + i;
            d[i] = -(i+1.0);
        }
    }

    void testSpeed()
    {
        using namespace vigra::multi_math;
        using namespace vigra::functor;
        const int size = 200;
        Shape3 s(size, size, size);
        array3_type u(s),
                    v(s),
                    w(s);
        std::string t;
        USETICTOC;
        std::cerr << "multi_math speed test: \n";
#if 0
        marray::Marray<double> mu(s.begin(), s.end()),
                               mv(s.begin(), s.end()),
                               mw(s.begin(), s.end());
        
        TIC;
        marray::Marray<double>::iterator iu = mu.begin(), end = mu.end(),
                                         iv = mv.begin(), iw = mw.begin();
        for(; iu != end; ++iu, ++iv, ++iw)
            *iw = *iu * *iv;
        t = TOCS;
        std::cerr << "    marray iterator: " << t << "\n";
        TIC;
        mw = mu*mv;
        t = TOCS;
        std::cerr << "    marray expression: " << t << "\n";
#endif
        //TIC;
        //typedef array3_type::view_type View;
        //View::iterator wi = ((View &)w).begin(), wend = wi.getEndIterator(),
        //               ui = ((View &)u).begin(), vi = ((View &)v).begin();
        //for(; wi != wend; ++wi, ++ui, ++vi)
        //            *wi = *ui * *vi;
        //t = TOCS;
        //std::cerr << "    StridedScanOrderIterator: " << t << "\n";
        TIC;
        typedef CoupledIteratorType<3, scalar_type, scalar_type, scalar_type>::type CI;
        CI i = createCoupledIterator(w, u, v), end = i.getEndIterator();
        for(; i != end; ++i)
                    i.get<1>() = i.get<2>() * i.get<3>();
        t = TOCS;
        std::cerr << "    CoupledScanOrderIterator: " << t << "\n";
        //TIC;
        //w = u*v;
        //t = TOCS;
        //std::cerr << "    multi_math expression: " << t << "\n";
        //TIC;
        //w.transpose() = u.transpose()*v.transpose();
        //t = TOCS;
        //std::cerr << "    transposed multi_math expression: " << t << "\n";
        //TIC;
        //combineTwoMultiArrays(srcMultiArrayRange(u), srcMultiArray(v), destMultiArray(w),
        //                      Arg1()*Arg2());
        //t = TOCS;
        //std::cerr << "    lambda expression: " << t << "\n";
        //TIC;
        //combineTwoMultiArrays(srcMultiArrayRange(u.transpose()), srcMultiArray(v.transpose()), destMultiArray(w),
        //                      Arg1()*Arg2());
        //t = TOCS;
        //std::cerr << "    transposed lambda expression: " << t << "\n";
        //TIC;
        //for(int z=0; z<size; ++z)
        //    for(int y=0; y<size; ++y)
        //        for(int x=0; x<size; ++x)
        //            w(x,y,z) = u(x,y,z) * v(x,y,z);
        //t = TOCS;
        //std::cerr << "    explicit loops: " << t << "\n";
        i = createCoupledIterator(w, u, v);
        TIC;
        for(i.dim<2>() = 0; i.dim<2>() <size; ++i.dim<2>())
            for(i.dim<1>() = 0; i.dim<1>()<size; ++i.dim<1>())
                for(i.dim<0>() = 0; i.dim<0>()<size; ++i.dim<0>())
                    i.get<1>() = i.get<2>() * i.get<3>();
        t = TOCS;
        std::cerr << "    coupled iterator explicit template loops: " << t << "\n";
        i = createCoupledIterator(w, u, v);
        TIC;
        for(i.resetDim(2); i.coord(2) < size; i.incDim(2))
            for(i.resetDim(1); i.coord(1) < size; i.incDim(1))
                for(i.resetDim(0); i.coord(0) < size; i.incDim(0))
                    i.get<1>() = i.get<2>() * i.get<3>();
        t = TOCS;
        std::cerr << "    coupled iterator explicit runtime loops: " << t << "\n";
    }

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

};


struct ImageViewTestSuite
: public vigra::test_suite
{
    ImageViewTestSuite()
    : vigra::test_suite("ImageViewTestSuite")
    {
		add( testCase( &BImageViewTest::testStridedImageView));
        add( testCase( &BImageViewTest::testBasicImageIterator));
        add( testCase( &BImageViewTest::testImageIterator));
        add( testCase( &BImageViewTest::copyImage));
		add( testCase( &DImageViewTest::testStridedImageView));
        add( testCase( &DImageViewTest::testBasicImageIterator));
        add( testCase( &DImageViewTest::testImageIterator));
        add( testCase( &DImageViewTest::copyImage));
		add( testCase( &BRGBImageViewTest::testStridedImageView));
        add( testCase( &BRGBImageViewTest::testBasicImageIterator));
        add( testCase( &BRGBImageViewTest::testImageIterator));
        add( testCase( &BRGBImageViewTest::copyImage));
		add( testCase( &FRGBImageViewTest::testStridedImageView));
        add( testCase( &FRGBImageViewTest::testBasicImageIterator));
        add( testCase( &FRGBImageViewTest::testImageIterator));
        add( testCase( &FRGBImageViewTest::copyImage));
    }
};

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
        add( testCase( &MultiArrayTest::test_hierarchical ) );
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

struct MultiArrayPointOperatorsTestSuite
: public vigra::test_suite
{
  MultiArrayPointOperatorsTestSuite()
    : vigra::test_suite("MultiArrayPointOperatorsTestSuite")
    {
        add( testCase( &MultiArrayPointoperatorsTest::testInit ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCopy ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCopyOuterExpansion ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCopyInnerExpansion ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTransform ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTransformOuterExpand ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTransformInnerExpand ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTransformOuterReduce ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTransformInnerReduce ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine2 ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine2OuterExpand ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine2InnerExpand ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine2OuterReduce ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine2InnerReduce ) );
        add( testCase( &MultiArrayPointoperatorsTest::testCombine3 ) );
        add( testCase( &MultiArrayPointoperatorsTest::testInitMultiArrayBorder ) );
        add( testCase( &MultiArrayPointoperatorsTest::testInspect ) );
        add( testCase( &MultiArrayPointoperatorsTest::testTensorUtilities ) );

        add( testCase( &MultiMathTest::testSpeed ) );
        add( testCase( &MultiMathTest::testBasicArithmetic ) );
        add( testCase( &MultiMathTest::testExpandMode ) );
        add( testCase( &MultiMathTest::testAllFunctions ) );
        add( testCase( &MultiMathTest::testComputedAssignment ) );
        add( testCase( &MultiMathTest::testNonscalarValues ) );
        add( testCase( &MultiMathTest::testMixedExpressions ) );
        add( testCase( &MultiMathTest::testComplex ) );
    }
}; // struct MultiArrayPointOperatorsTestSuite


int main(int argc, char ** argv)
{
    // run the multi-array testsuite
    MultiArrayTestSuite test1;
    int failed = test1.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test1.report() << std::endl;

    // run the multi-array data-testsuite
    MultiArrayDataTestSuite test1a;
    failed += test1a.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test1a.report() << std::endl;

    // run the image testsuite
    ImageViewTestSuite test2;
    failed += test2.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test2.report() << std::endl;

    // run the multi-array point operator test suite
    MultiArrayPointOperatorsTestSuite test3;
    failed += test3.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test3.report() << std::endl;
    
    return (failed != 0);
}


