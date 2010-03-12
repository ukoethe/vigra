/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
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

#include "unittest.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_impex.hxx"
#include "vigra/basicimageview.hxx"
#include "vigra/navigator.hxx"
#include "vigra/multi_pointoperators.hxx"
#include "vigra/functorexpression.hxx"

using namespace vigra;
using namespace vigra::functor;

class MultiArrayDataTest
{
public:

    // test a three-dimensional image of size 10*10*10.
    // the image is filled from start to beginning (memory-wise) with
    // ascending numbers.

    typedef int scalar_type;
    typedef MultiArray <2, scalar_type> array2_type;
    typedef MultiArray <3, scalar_type> array3_type;
    typedef MultiArrayView <3, scalar_type> array3_view_type;
    typedef array2_type::difference_type difference2_type;
    typedef array3_type::difference_type difference3_type;
    
    difference3_type shape3;
    array3_type array3;

    MultiArrayDataTest ()
        : shape3 (10, 10, 10), array3 (shape3, 1)
    {
        // initialize the array to the test data
        for (int i = 0; i < 1000; ++i)
            array3.data () [i] = i;
    }

    void testHasData ()
    {
        should(array3.hasData());
        array3_type empty;
        should(!empty.hasData());
    }

    void testEquality ()
    {
        typedef difference3_type Shape;
        should(array3 == array3);
        should(array3 != array3.subarray(Shape(1,1,1), Shape(2,2,2)));
        should(array3.subarray(Shape(0,0,0), Shape(10,1,1)) != array3.subarray(Shape(0,1,0), Shape(10,2,1)));
        for(unsigned int k=0; k<10; ++k)
            array3(k,0,0) += 10;
        should(array3.subarray(Shape(0,0,0), Shape(10,1,1)) == array3.subarray(Shape(0,1,0), Shape(10,2,1)));
    }
    
    // bindInner tests
    void test_bindInner ()
    {
        TinyVector <int, 2> inner_indices (2, 5);
        MultiArrayView <1, scalar_type, StridedArrayTag>
            array = array3.bindInner <2> (inner_indices);
        shouldEqual ((array [TinyVector <int, 1> (0)]), 52);
        shouldEqual ((array [TinyVector <int, 1> (1)]), 152);
        shouldEqual ((array (0)), 52);
        shouldEqual ((array (1)), 152);
    }
    
    // bindOuter tests
    void test_bindOuter ()
    {
        TinyVector <int, 2> outer_indices (2, 5);
        MultiArrayView <1, scalar_type, UnstridedArrayTag>
            array = array3.bindOuter <2> (outer_indices);
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
        MultiArrayView <2, scalar_type, UnstridedArrayTag>
            array = array3.bind <1> (4);
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
            array = array3.bind <0> (4);
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

        MultiArrayView <4, scalar_type, UnstridedArrayTag> a1 = array3.insertSingletonDimension(1);
        shouldEqual ((a1 [TinyVector <int, 4> (4, 0, 0, 0)]), 4);
        shouldEqual ((a1 [TinyVector <int, 4> (4, 0, 1, 0)]), 14);
        shouldEqual ((a1 [TinyVector <int, 4> (4, 0, 0, 1)]), 104);

        MultiArrayView <4, scalar_type, UnstridedArrayTag> a3 = array3.insertSingletonDimension(3);
        shouldEqual ((a3 [TinyVector <int, 4> (4, 0, 0, 0)]), 4);
        shouldEqual ((a3 [TinyVector <int, 4> (4, 1, 0, 0)]), 14);
        shouldEqual ((a3 [TinyVector <int, 4> (4, 0, 1, 0)]), 104);
    }

    // subarray tests
    void test_subarray ()
    {
        typedef difference3_type Shape;
        
        Shape offset (1,1,1);
        Shape size (5,5,5);
        MultiArrayView <3, scalar_type, UnstridedArrayTag>
            array = array3.subarray (offset, size);
        shouldEqual (array [Shape (0,0,0)], 111);
        shouldEqual (array [Shape (5,2,1)], 236);
        shouldEqual (array (0,0,0), 111);
        shouldEqual (array (5,2,1), 236);
        
        // test swap
        array3.subarray(Shape(0,0,0), Shape(10,10,1)).swapData( 
             array3.subarray(Shape(0,0,1), Shape(10,10,2)));
        for(unsigned int k=0; k<100; ++k)
            shouldEqual(array3[k], k+100);
        for(unsigned int k=100; k<200; ++k)
            shouldEqual(array3[k], k-100);
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

    // permute and transpose tests
    void testPermute ()
    {   
        array3.reshape(difference3_type(3,5,7));
        for(int k=0; k<array3.size(); ++k)
            array3[k] = k;

        array3_type ref(difference3_type(array3.shape(2), array3.shape(0), array3.shape(1)));
        for(int k=0; k<array3.shape(0); ++k)
            for(int l=0; l<array3.shape(1); ++l)
                for(int m=0; m<array3.shape(2); ++m)
                    ref(m,k,l) = array3(k,l,m);

        MultiArrayView <3, scalar_type, StridedArrayTag>
                parray = array3.permuteDimensions (difference3_type (2, 0, 1));        
        shouldEqual(ref.shape(), parray.shape());
        should(ref == parray);

        MultiArrayView <3, scalar_type, StridedArrayTag>
                parray_ascending = parray.permuteStridesAscending();        
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

        array2_type ref2(difference2_type(array3.shape(1), array3.shape(0)));
        for(int k=0; k<array3.shape(0); ++k)
            for(int l=0; l<array3.shape(1); ++l)
                ref2(l, k) = array3(k, l, 0);

        MultiArrayView <2, scalar_type, StridedArrayTag>
                array2 = array3.bindOuter(0).transpose ();
        shouldEqual(ref2.shape(), array2.shape());
        should(ref2 == array2);

        try {
            array3.permuteDimensions (difference3_type (2, 0, 0));   
            failTest("no exception thrown");
        }
        catch(vigra::ContractViolation & c)
        {
            std::string expected("\nPrecondition violation!\nMultiArrayView::permuteDimensions(): every dimension must occur exactly once");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
    }

    void testNorm ()
    {
        shouldEqual(array3.squaredNorm(), 332833500);
        
        shouldEqual(array3.norm(), std::sqrt(332833500.0));
        shouldEqual(array3.norm(0), 999.0);
        shouldEqual(array3.norm(1), 499500.0);
        shouldEqualTolerance(array3.norm(2, false), std::sqrt(332833500.0), 1e-14);
        
        difference3_type first(0,0,0), last(0,1,1);
        shouldEqual(array3.subarray(first, last).norm(), 0.0);
        shouldEqual(array3.subarray(first, last).norm(0), 0.0);
        shouldEqual(array3.subarray(first, last).norm(1), 0.0);
        shouldEqual(array3.subarray(first, last).norm(2, false), 0.0);

        shouldEqual(array3.squaredNorm(), squaredNorm(array3));
        shouldEqual(array3.norm(), vigra::norm(array3));
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
            
        typedef MultiArrayView <2, scalar_type, UnstridedArrayTag>::difference_type Shape;
        MultiArrayView <2, scalar_type, UnstridedArrayTag>
            subarray = array3.bindOuter(2).subarray(Shape(1,0), Shape(10,9));
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
        typedef MultiArrayView <3, scalar_type, UnstridedArrayTag>::difference_type Shape;
        MultiArrayView <3, scalar_type, UnstridedArrayTag> array;
        array = array3;
        should(array3 == array);
        try {
            array = array3.subarray(Shape(0,0,0), Shape(10,1,1));
            failTest("no exception thrown");
        }
        catch(vigra::ContractViolation & c)
        {
            std::string expected("\nPrecondition violation!\nMultiArrayView::operator=(MultiArrayView const &) size mismatch");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
        MultiArrayView <3, scalar_type, UnstridedArrayTag> subarray = array3.subarray(Shape(0,0,0), Shape(10,1,1));
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
    }
};


class MultiArrayTest
{
public:
    typedef MultiArray <1, unsigned char> array1_t;
    typedef array1_t::difference_type shape1_t;
    typedef MultiArray <3, unsigned char> array3_t;
    typedef array3_t::difference_type shape3_t;
    typedef array3_t::traverser iter3_t;
    typedef iter3_t::next_type iter2_t;
    typedef iter2_t::next_type iter1_t;
    
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
        TinyVector <int, 1> s;
        s[0] = 2;
        typedef MultiArray <1, unsigned char> array1_t;
        array1_t a (s);
        should (a.shape (0) == 2);
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

    void test_iter ()
    {
        // test hierarchical navigation and 
        iter3_t i3_f = a3.traverser_begin ();
        iter3_t i3_l = a3.traverser_end ();
        array3_t::iterator seqi = a3.begin();

        unsigned int countx = 0, county = 0, countz = 0;

        // iterate over the third dimension
        for (int z=0; i3_f != i3_l; ++i3_f, ++z) 
        {
            iter2_t i2_f = i3_f.begin ();
            iter2_t i2_l = i3_f.end ();
            // iterate over the second dimension
            for (int y=0; i2_f != i2_l; ++i2_f, ++y) 
            {
                iter1_t i1_f = i2_f.begin ();
                iter1_t i1_l = i2_f.end ();
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
        iter3_t i3 = a3.traverser_begin();
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

    void test_const_iter ()
    {
        // test hierarchical navigation
        iter3_t i3_f = a3.traverser_begin ();
        iter3_t i3_l = a3.traverser_end ();

        unsigned int countx = 0, county = 0, countz = 0;

        // iterate over the third dimension
        for (; i3_f != i3_l; ++i3_f) {
            iter2_t i2_f = i3_f.begin ();
            iter2_t i2_l = i3_f.end ();
            // iterate over the second dimension
            for (; i2_f != i2_l; ++i2_f) {
                iter1_t i1_f = i2_f.begin ();
                iter1_t i1_l = i2_f.end ();
                // iterate over the first dimension
                for (; i1_f != i1_l; ++i1_f)
                    ++countx;
                ++county;
            }
            ++countz;
        }

        shouldEqual (countx, 30);
        shouldEqual (county, 15);
        shouldEqual (countz, 5);
        
        // test direct navigation
        iter3_t i3 = a3.traverser_begin();
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
        MultiArrayView <2, unsigned char> ba = 
            a3.bindOuter<1> (TinyVector<int, 1>(2));

        shouldEqual (ba.shape (0), 2);
        shouldEqual (ba.shape (1), 3);
    }


    void test_bindInner ()
    {
        MultiArrayView <2, unsigned char, StridedArrayTag>
            fa = a3.bindInner <1> (TinyVector <int, 1>(1));

        shouldEqual (fa.shape (0), 3);
        shouldEqual (fa.shape (1), 5);
    }

    void test_bindInnerAll ()
    {
        MultiArrayView <0, unsigned char, StridedArrayTag>
            fa = a3.bindInner <3> (shape3_t(1,1,1));

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

        shouldEqual (countx, 336);
        shouldEqual (county, 56);
        shouldEqual (countz, 8);
  
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

    // bindInner tests
    void testNavigator ()
    {
        int expected[][24] = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
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
        result.reshape(import_info.shape());
        importVolume(import_info, result);
        
        shouldEqual(result.shape(), Shape(2,3,4));
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

#ifdef _WIN32
        exportVolume(array, VolumeExportInfo("impex\\test", ext2));
        
        importVolume(result, std::string("impex\\test"), std::string(ext2));
        
        shouldEqual(result.shape(), Shape(2,3,4));
        shouldEqual(result(0,1,0), 1);
        shouldEqual(result(0,1,1), 2);
        shouldEqual(result(0,1,2), 3);
        shouldEqual(result(0,1,3), 4);
#endif // _WIN32
    }
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
    }

    void testCopy()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        copyMultiArray(srcMultiArrayRange(img), destMultiArray(res));
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(x,y,z));
    }

    void testCopyOuterExpansion()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(5,1,1));
        
        copyMultiArray(srcMultiArrayRange(view), destMultiArrayRange(res));
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(x,0,0));
    }

    void testCopyInnerExpansion()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(1,1,3));
        
        copyMultiArray(srcMultiArrayRange(view), destMultiArrayRange(res));
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(0,0,z));
    }

    void testTransform()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        transformMultiArray(srcMultiArrayRange(img), destMultiArray(res),
                            Arg1() + Arg1());
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,y,z));
    }

    void testTransformOuterExpand()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(5,1,1));
        
        transformMultiArray(srcMultiArrayRange(view), destMultiArrayRange(res),
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
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(1,1,3));
        
        transformMultiArray(srcMultiArrayRange(view), destMultiArrayRange(res),
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
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        transformMultiArray(srcMultiArrayRange(img), destMultiArrayRange(res),
                            reduceFunctor(Arg1() + Arg2(), 0.0));
        
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
        transformMultiArray(srcMultiArrayRange(img), 
                            destMultiArrayRange(res3),
                            FindSum<PixelType>());
        shouldEqualSequenceTolerance(res1.data(), res1.data()+5, res.data(), 1e-7);       
    }

    void testTransformInnerReduce()
    {
        Image3D res(Size3(1,1,3));
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        transformMultiArray(srcMultiArrayRange(img), destMultiArrayRange(res),
                            reduceFunctor(Arg1() + Arg2(), 0.0));
        
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
        transformMultiArray(srcMultiArrayRange(img), 
                            destMultiArrayRange(res3),
                            FindSum<PixelType>());
        shouldEqualSequenceTolerance(res1.data(), res1.data()+3, res.data(), 1e-6);       
    }

    void testCombine2()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArray(img), 
                              destMultiArray(res),
                              Arg1() + Arg2());
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,y,z));
    }

    void testCombine2OuterExpand()
    {
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(5,1,1));
        combineTwoMultiArrays(srcMultiArrayRange(view), srcMultiArrayRange(img), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,y,z) + img(x,0,0));

        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArrayRange(view), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), img(x,y,z) + 2.0*img(x,0,0));

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
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        View3D view = img.subarray(Size3(0,0,0), Size3(1,1,3));
        combineTwoMultiArrays(srcMultiArrayRange(view), srcMultiArrayRange(img), 
                              destMultiArrayRange(res),
                              Arg1() + Param(2.0)*Arg2());       
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 2.0*img(x,y,z) + img(0,0,z));

        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArrayRange(view), 
                              destMultiArrayRange(res),
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
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArrayRange(img), 
                              destMultiArrayRange(res),
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
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        combineTwoMultiArrays(srcMultiArrayRange(img), srcMultiArrayRange(img), 
                              destMultiArrayRange(res),
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
        Image3D res(img.shape());
        
        initMultiArray(destMultiArrayRange(res), 0.0);
        
        combineThreeMultiArrays(srcMultiArrayRange(img), 
                                srcMultiArray(img), srcMultiArray(img), 
                                destMultiArray(res),
                                Arg1() + Arg2() + Arg3());
        
        int x,y,z;
        for(z=0; z<img.shape(2); ++z)
            for(y=0; y<img.shape(1); ++y)
                for(x=0; x<img.shape(0); ++x)
                    shouldEqual(res(x,y,z), 3.0*img(x,y,z));
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
        initMultiArrayBorder(destMultiArrayRange(vol2),9,0);
        shouldEqualSequence(vol2.begin(), vol2.end(), desired_vol2);

    }
};


struct ImageViewTestSuite
: public vigra::test_suite
{
    ImageViewTestSuite()
    : vigra::test_suite("ImageViewTestSuite")
    {
        add( testCase( &BImageViewTest::testBasicImageIterator));
        add( testCase( &BImageViewTest::testImageIterator));
        add( testCase( &BImageViewTest::copyImage));
        add( testCase( &DImageViewTest::testBasicImageIterator));
        add( testCase( &DImageViewTest::testImageIterator));
        add( testCase( &DImageViewTest::copyImage));
        add( testCase( &BRGBImageViewTest::testBasicImageIterator));
        add( testCase( &BRGBImageViewTest::testImageIterator));
        add( testCase( &BRGBImageViewTest::copyImage));
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
        add( testCase( &MultiArrayTest::test_iter ) );
        add( testCase( &MultiArrayTest::test_const_iter ) );
        add( testCase( &MultiArrayTest::test_bindOuter ) );
        add( testCase( &MultiArrayTest::test_bindInner ) );
        add( testCase( &MultiArrayTest::test_bindInnerAll ) );
        add( testCase( &MultiArrayTest::test_bindAt ) );
        add( testCase( &MultiArrayTest::test_bind ) );
        add( testCase( &MultiArrayTest::test_reshape) );
        add( testCase( &MultiArrayTest::test_subarray ) );
        add( testCase( &MultiArrayTest::test_stridearray ) );
        add( testCase( &MultiArrayTest::test_copy_int_float ) );

        add( testCase( &MultiImpexTest::testImpex ) );
	}
};

struct MultiArrayDataTestSuite
: public vigra::test_suite
{
    MultiArrayDataTestSuite()
    : vigra::test_suite("MultiArrayDataTestSuite")
    {
        add( testCase( &MultiArrayDataTest::testHasData ) );
        add( testCase( &MultiArrayDataTest::testEquality ) );
        add( testCase( &MultiArrayDataTest::test_subarray ) );
        add( testCase( &MultiArrayDataTest::test_stridearray ) );
        add( testCase( &MultiArrayDataTest::test_bindOuter ) );
        add( testCase( &MultiArrayDataTest::test_bindInner ) );
        add( testCase( &MultiArrayDataTest::test_bindAt ) );
        add( testCase( &MultiArrayDataTest::test_bind ) );
        add( testCase( &MultiArrayDataTest::test_bind0 ) );
        add( testCase( &MultiArrayDataTest::test_singletonDimension ) );
        add( testCase( &MultiArrayDataTest::testPermute ) );
        add( testCase( &MultiArrayDataTest::testNorm ) );
        add( testCase( &MultiArrayDataTest::testScanOrderAccess ) );
        add( testCase( &MultiArrayDataTest::testAssignmentAndReset ) );
        add( testCase( &MultiArrayNavigatorTest::testNavigator ) );
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


