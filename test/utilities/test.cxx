#include <iostream>
#include <iterator>
#include "unittest.hxx"
#include "vigra/accessor.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/copyimage.hxx"

using namespace vigra;

struct ArrayVectorTest
{
    typedef int value_type;
    typedef ArrayVector<value_type> Vector;
    typedef Vector::iterator Iterator;
    typedef StandardValueAccessor<value_type> Accessor;
    
    Vector vector_;

    ArrayVectorTest()
    {}
    
    void testAccessor()
    {
        vector_.resize(3, 0);
        Iterator i = vector_.begin();
        
        StandardValueAccessor<value_type> sva;
        StandardConstValueAccessor<value_type> scva;
        sva.set(3, i);
        shouldEqual(vector_[0], 3);
        shouldEqual(sva(i), 3);
        shouldEqual(scva(i), 3);
        sva.set(4, i+1);
        shouldEqual(vector_[1], 4);
        shouldEqual(sva(i+1), 4);
        shouldEqual(scva(i+1), 4);
        sva.set(5, i, 2);
        shouldEqual(vector_[2], 5);
        shouldEqual(sva(i, 2), 5);
        shouldEqual(scva(i, 2), 5);
        
        StandardAccessor<value_type> sa;
        StandardConstAccessor<value_type> sca;
        sa.set(6, i);
        shouldEqual(vector_[0], 6);
        shouldEqual(sa(i), 6);
        shouldEqual(sca(i), 6);
        sa.set(7, i+1);
        shouldEqual(vector_[1], 7);
        shouldEqual(sa(i+1), 7);
        shouldEqual(sca(i+1), 7);
        sa.set(8, i, 2);
        shouldEqual(vector_[2], 8);
        shouldEqual(sa(i, 2), 8);
        shouldEqual(sca(i, 2), 8);
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

struct ArrayVectorTestSuite
: public vigra::test_suite
{
    ArrayVectorTestSuite()
    : vigra::test_suite("ArrayVectorTestSuite")
    {
        add( testCase( &ArrayVectorTest::testAccessor));
        add( testCase( &ArrayVectorTest::testBackInsertion));
    }
};

int main()
{
    ArrayVectorTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

