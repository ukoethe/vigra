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

#include <iostream>
#include <iterator>
#include "unittest.hxx"
#include "vigra/accessor.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/copyimage.hxx"
#include "vigra/sized_int.hxx"

using namespace vigra;

struct ArrayVectorTest
{
    typedef int value_type;
    typedef ArrayVector<value_type> Vector;
    typedef Vector::iterator Iterator;
    typedef Vector::const_iterator ConstIterator;
    typedef StandardValueAccessor<value_type> Accessor;
    
    Vector vector_;

    ArrayVectorTest()
    {}
    
    void testAccessor()
    {
        vector_.resize(3, 0);
        Iterator i = vector_.begin();
        ConstIterator ci = const_cast<Vector const &>(vector_).begin();
        
        StandardValueAccessor<value_type> sva;
        StandardConstValueAccessor<value_type> scva;
        sva.set(3, i);
        shouldEqual(vector_[0], 3);
        shouldEqual(sva(i), 3);
        shouldEqual(scva(i), 3);
        shouldEqual(scva(ci), 3);
        sva.set(4, i+1);
        shouldEqual(vector_[1], 4);
        shouldEqual(sva(i+1), 4);
        shouldEqual(scva(i+1), 4);
        shouldEqual(scva(ci+1), 4);
        sva.set(5, i, 2);
        shouldEqual(vector_[2], 5);
        shouldEqual(sva(i, 2), 5);
        shouldEqual(scva(i, 2), 5);
        shouldEqual(scva(ci, 2), 5);
        
        StandardAccessor<value_type> sa;
        StandardConstAccessor<value_type> sca;
        sa.set(6, i);
        shouldEqual(vector_[0], 6);
        shouldEqual(sa(i), 6);
        shouldEqual(sca(i), 6);
        shouldEqual(sca(ci), 6);
        sa.set(7, i+1);
        shouldEqual(vector_[1], 7);
        shouldEqual(sa(i+1), 7);
        shouldEqual(sca(i+1), 7);
        shouldEqual(sca(ci+1), 7);
        sa.set(8, i, 2);
        shouldEqual(vector_[2], 8);
        shouldEqual(sa(i, 2), 8);
        shouldEqual(sca(i, 2), 8);
        shouldEqual(sca(ci, 2), 8);
        
        Vector varray[] = { vector_, vector_, vector_ };
        Vector  * v = varray;
        Vector const * cv = varray;
        int k;
        
        VectorComponentAccessor<Vector> vca(0);
        for(k = 0; k<2; ++k, vca.setIndex(k))
        {
            shouldEqual(vca(v), 6 + k);
            shouldEqual(vca(v + 1), 6 + k);
            shouldEqual(vca(v, 2), 6 + k);
            shouldEqual(vca(cv), 6 + k);
            shouldEqual(vca(cv + 1), 6 + k);
            shouldEqual(vca(cv, 2), 6 + k);
            vca.set(3 + k, v);
            vca.set(3 + k, v + 1);
            vca.set(3 + k, v, 2);
            shouldEqual(v[0][k], 3 + k);
            shouldEqual(v[1][k], 3 + k);
            shouldEqual(v[2][k], 3 + k);
        }
        
        VectorComponentValueAccessor<Vector> vcva(0);
        for(k = 0; k<2; ++k, vcva.setIndex(k))
        {
            shouldEqual(vcva(v), 3 + k);
            shouldEqual(vcva(v + 1), 3 + k);
            shouldEqual(vcva(v, 2), 3 + k);
            shouldEqual(vcva(cv), 3 + k);
            shouldEqual(vcva(cv + 1), 3 + k);
            shouldEqual(vcva(cv, 2), 3 + k);
            vcva.set(6 + k, v);
            vcva.set(6 + k, v + 1);
            vcva.set(6 + k, v, 2);
            shouldEqual(v[0][k], 6 + k);
            shouldEqual(v[1][k], 6 + k);
            shouldEqual(v[2][k], 6 + k);
        }

        VectorAccessor<Vector> va;
        VectorElementAccessor<VectorAccessor<Vector> > vea(0, va);
        for(k = 0; k<2; ++k, vea.setIndex(k))
        {
            shouldEqual(vea(v), 6 + k);
            shouldEqual(vea(v + 1), 6 + k);
            shouldEqual(vea(v, 2), 6 + k);
            shouldEqual(vea(cv), 6 + k);
            shouldEqual(vea(cv + 1), 6 + k);
            shouldEqual(vea(cv, 2), 6 + k);
            vea.set(3 + k, v);
            vea.set(3 + k, v + 1);
            vea.set(3 + k, v, 2);
            shouldEqual(v[0][k], 3 + k);
            shouldEqual(v[1][k], 3 + k);
            shouldEqual(v[2][k], 3 + k);
        }

        for(k = 0; k<2; ++k)
        {
            shouldEqual(va.getComponent(v, k), 3 + k);
            shouldEqual(va.getComponent(v + 1, k), 3 + k);
            shouldEqual(va.getComponent(v, 2, k), 3 + k);
            shouldEqual(va.getComponent(cv, k), 3 + k);
            shouldEqual(va.getComponent(cv + 1, k), 3 + k);
            shouldEqual(va.getComponent(cv, 2, k), 3 + k);
            va.setComponent(6 + k, v, k);
            va.setComponent(6 + k, v + 1, k);
            va.setComponent(6 + k, v, 2, k);
            shouldEqual(v[0][k], 6 + k);
            shouldEqual(v[1][k], 6 + k);
            shouldEqual(v[2][k], 6 + k);
        }

        SequenceAccessor<Vector> sqa;
        SequenceAccessor<const Vector> sqca;
        shouldEqual(sqa.size(v), 3);
        shouldEqual(sqa.size(v + 1), 3);
        shouldEqual(sqa.size(v, 2), 3);
        shouldEqual(sqca.size(cv), 3);
        shouldEqual(sqca.size(cv + 1), 3);
        shouldEqual(sqca.size(cv, 2), 3);
        should(sqa.end(v) == v[0].end());
        should(sqa.end(v + 1) == v[1].end());
        should(sqa.end(v, 2) == v[2].end());
        should(sqca.end(cv) == cv[0].end());
        should(sqca.end(cv + 1) == cv[1].end());
        should(sqca.end(cv, 2) == cv[2].end());
        for(k = 0; k<2; ++k)
        {
            shouldEqual(sqa.begin(v)[k], 6 + k);
            shouldEqual(sqa.begin(v + 1)[k], 6 + k);
            shouldEqual(sqa.begin(v, 2)[k], 6 + k);
            shouldEqual(sqca.begin(cv)[k], 6 + k);
            shouldEqual(sqca.begin(cv + 1)[k], 6 + k);
            shouldEqual(sqca.begin(cv, 2)[k], 6 + k);
            sqa.begin(v)[k] = 3 + k;
            sqa.begin(v + 1)[k] = 3 + k;
            sqa.begin(v, 2)[k] = 3 + k;
            shouldEqual(v[0][k], 3 + k);
            shouldEqual(v[1][k], 3 + k);
            shouldEqual(v[2][k], 3 + k);
        }
    }

    void testBackInsertion()
    {
        static value_type data[] = { 0, 1, 2, 3, 4 };
        
        shouldEqual(vector_.size(), 0);
        
        Accessor a;
        copyLine(data, data + 5, a, std::back_inserter(vector_), a);
        
        shouldEqual(vector_.size(), 5);
        shouldEqualSequence(vector_.begin(), vector_.end(), data);
    }
};

struct SizedIntTest
{
    void testSizedInt()
    {
        shouldEqual(sizeof(Int8), 1);
        shouldEqual(sizeof(Int16), 2);
        shouldEqual(sizeof(Int32), 4);
        shouldEqual(sizeof(UInt8), 1);
        shouldEqual(sizeof(UInt16), 2);
        shouldEqual(sizeof(UInt32), 4);
        should(sizeof(IntBiggest) >= 4);
        should(sizeof(UIntBiggest) >= 4);
    }
};

struct UtilitiesTestSuite
: public vigra::test_suite
{
    UtilitiesTestSuite()
    : vigra::test_suite("UtilitiesTestSuite")
    {
        add( testCase( &ArrayVectorTest::testAccessor));
        add( testCase( &ArrayVectorTest::testBackInsertion));
        add( testCase( &SizedIntTest::testSizedInt));
    }
};

int main(int argc, char ** argv)
{
    UtilitiesTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

