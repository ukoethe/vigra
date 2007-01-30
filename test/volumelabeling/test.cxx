#include <iostream>
#include <functional>
#include <cmath>
#include "unittest.hxx"

#include "vigra/labelvolume.hxx"

using namespace vigra;

struct LabelingTest
{
    typedef vigra::MultiArray<3,int> IntVolume;
    typedef vigra::MultiArray<3,double> DoubleVolume;

    LabelingTest()
    : vol1(IntVolume::difference_type(4,4,4)),vol2(IntVolume::difference_type(4,4,4)),vol3(IntVolume::difference_type(5,5,5)),vol4(DoubleVolume::difference_type(5,5,5)),vol5(DoubleVolume::difference_type(5,5,5))
    {
        static const int in1[] = { 0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0,
                                   0, 0, 0, 0,    0, 1, 1, 0,    0, 1, 1, 0,    0, 0, 0, 0,
                                   0, 0, 0, 0,    0, 1, 1, 0,    0, 1, 1, 0,    0, 0, 0, 0,
                                   0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0};

        IntVolume::iterator i = vol1.begin();
        IntVolume::iterator end = vol1.end();
        const int * p = in1;

        for(; i != end; ++i, ++p)
        {
            *i=*p;
        }

        static const int in2[] = { 0, 1, 0, 1,    1, 0, 1, 0,    0, 1, 0, 1,    1, 0, 1, 0,
                                   1, 0, 1, 0,    0, 1, 0, 1,    1, 0, 1, 0,    0, 1, 0, 1,
                                   0, 1, 0, 1,    1, 0, 1, 0,    0, 1, 0, 1,    1, 0, 1, 0,
                                   1, 0, 1, 0,    0, 1, 0, 1,    1, 0, 1, 0,    0, 1, 0, 1};

        i = vol2.begin();
        end = vol2.end();
        p = in2;

        for(; i != end; ++i, ++p)
        {
            *i=*p;
        }

                        
        static const int in3[] = { 0, 1, 0, 0, 0,    1, 1, 0, 0, 0,    0, 0, 0, 0, 0,    0, 0, 0, 0, 0,    0, 0, 0, 0, 0,
                                   1, 1, 1, 0, 0,    1, 1, 1, 0, 0,    1, 1, 1, 0, 0,    0, 0, 0, 0, 0,    0, 0, 0, 0, 0,
                                   0, 0, 1, 0, 0,    0, 0, 1, 0, 0,    1, 1, 1, 0, 0,    0, 0, 0, 0, 0,    0, 0, 0, 0, 0,
                                   1, 1, 1, 1, 0,    1, 1, 1, 1, 0,    1, 1, 1, 1, 0,    1, 1, 1, 1, 0,    0, 0, 0, 0, 0,
                                   0, 0, 0, 1, 0,    0, 0, 0, 1, 0,    0, 0, 0, 1, 0,    1, 1, 1, 1, 0,    0, 0, 0, 0, 0};

        i = vol3.begin();
        end = vol3.end();
        p = in3;

        for(; i != end; ++i, ++p)
        {
            *i=*p;
        }

                static const double in4[] = { 1.0, 0.0, 0.0, 0.0, 1.0,    0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,    1.0, 0.0, 0.0, 0.0, 1.0,
                                              0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 1.0, 0.0, 1.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 1.0, 0.0, 1.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,
                                              0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,    1.0, 1.0, 1.0, 1.0, 1.0,    0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,
                                              0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 1.0, 0.0, 1.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 1.0, 0.0, 1.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,
                                              1.0, 0.0, 0.0, 0.0, 1.0,    0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,    1.0, 0.0, 0.0, 0.0, 1.0};

        DoubleVolume::iterator id = vol4.begin();
        DoubleVolume::iterator endd = vol4.end();
        const double * pd = in4;

        for(; id != endd; ++id, ++pd)
        {
            *id=*pd;
        }

        static const double in5[] = { 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,
                                      2.0, 2.0, 0.0, 2.0, 2.0,    2.0, 1.0, 0.0, 1.0, 2.0,    2.0, 2.0, 0.0, 2.0, 2.0,    2.0, 1.0, 0.0, 1.0, 2.0,    2.0, 2.0, 0.0, 2.0, 2.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 2.0, 2.0, 2.0, 0.0,    0.0, 2.0, 1.0, 2.0, 0.0,    0.0, 2.0, 2.0, 2.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,
                                      2.0, 2.0, 0.0, 2.0, 2.0,    2.0, 1.0, 0.0, 1.0, 2.0,    2.0, 2.0, 0.0, 2.0, 2.0,    2.0, 1.0, 0.0, 1.0, 2.0,    2.0, 2.0, 0.0, 2.0, 2.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0};

        id = vol5.begin();
        endd = vol5.end();
        pd = in5;

        for(; id != endd; ++id, ++pd)
        {
            *id=*pd;
        }

    }

    void labelingSixTest1()
    {
        IntVolume res(vol1);

        should(2 == labelVolumeSix(srcMultiArrayRange(vol1), destMultiArray(res)));

        IntVolume::iterator i1 = vol1.begin();
        IntVolume::iterator i1end = vol1.end();
        IntVolume::iterator i2 = res.begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should( *i1 == (*i2 - 1.0) );
        }
    }

        void labelingSixTest2()
    {
        IntVolume res(vol2);

        should(64 == labelVolume(srcMultiArrayRange(vol2), destMultiArray(res), NeighborCode3DSix()));

        IntVolume::iterator i2 = res.begin();
        IntVolume::iterator i2end = res.end();
        int address = 0;

        for(; i2 != i2end; ++i2, ++address)
        {
            should( *i2 == address+1 );
        }
    }

        void labelingSixTest3()
    {
        IntVolume res(vol3);

        should(5 == labelVolume(srcMultiArrayRange(vol3), destMultiArray(res), NeighborCode3DSix()));

        static const int out3[] = { 1, 2, 3, 3, 3,    2, 2, 3, 3, 3,    3, 3, 3, 3, 3,    3, 3, 3, 3, 3,    3, 3, 3, 3, 3,
                                    2, 2, 2, 3, 3,    2, 2, 2, 3, 3,    2, 2, 2, 3, 3,    3, 3, 3, 3, 3,    3, 3, 3, 3, 3,
                                    4, 4, 2, 3, 3,    4, 4, 2, 3, 3,    2, 2, 2, 3, 3,    3, 3, 3, 3, 3,    3, 3, 3, 3, 3,
                                    2, 2, 2, 2, 3,    2, 2, 2, 2, 3,    2, 2, 2, 2, 3,    2, 2, 2, 2, 3,    3, 3, 3, 3, 3,
                                    5, 5, 5, 2, 3,    5, 5, 5, 2, 3,    5, 5, 5, 2, 3,    2, 2, 2, 2, 3,    3, 3, 3, 3, 3};

        IntVolume::iterator i2 = res.begin();
        IntVolume::iterator i2end = res.end();
        const int * p = out3;

        for(; i2 != i2end; ++i2, ++p)
        {
            should( *i2 == *p );
        }
    }

    void labelingSixTest4()
    {
        IntVolume res(vol4);

        should(18 == labelVolume(srcMultiArrayRange(vol4), destMultiArray(res), NeighborCode3DSix()));

        static const int out4[] = { 1, 2, 2, 2, 3,    2, 2, 2, 2, 2,    2, 2, 4, 2, 2,    2, 2, 2, 2, 2,    5, 2, 2, 2, 6,
                                    2, 2, 2, 2, 2,    2, 7, 2, 8, 2,    2, 2, 4, 2, 2,    2, 9, 2,10, 2,    2, 2, 2, 2, 2,
                                    2, 2, 4, 2, 2,    2, 2, 4, 2, 2,    4, 4, 4, 4, 4,    2, 2, 4, 2, 2,    2, 2, 4, 2, 2,
                                    2, 2, 2, 2, 2,    2,11, 2,12, 2,    2, 2, 4, 2, 2,    2,13, 2,14, 2,    2, 2, 2, 2, 2,
                                   15, 2, 2, 2,16,    2, 2, 2, 2, 2,    2, 2, 4, 2, 2,    2, 2, 2, 2, 2,   17, 2, 2, 2,18};

        IntVolume::iterator i2 = res.begin();
        IntVolume::iterator i2end = res.end();
        const int * p = out4;

        for(; i2 != i2end; ++i2, ++p)
        {
            should( *i2 == *p );
        }
    }

    void labelingSixWithBackgroundTest1()
    {
        IntVolume res(vol5);

        should(4 == labelVolumeWithBackground(srcMultiArrayRange(vol5), destMultiArray(res), NeighborCode3DSix(), 0));

        static const int out5[] = { 0, 0, 0, 0, 0,    0, 1, 1, 1, 0,    0, 1, 1, 1, 0,    0, 1, 1, 1, 0,    0, 0, 0, 0, 0,
                                    2, 2, 0, 2, 2,    2, 1, 0, 1, 2,    2, 2, 0, 2, 2,    2, 1, 0, 1, 2,    2, 2, 0, 2, 2,
                                    0, 0, 0, 0, 0,    0, 2, 2, 2, 0,    0, 2, 3, 2, 0,    0, 2, 2, 2, 0,    0, 0, 0, 0, 0,
                                    2, 2, 0, 2, 2,    2, 4, 0, 4, 2,    2, 2, 0, 2, 2,    2, 4, 0, 4, 2,    2, 2, 0, 2, 2,
                                    0, 0, 0, 0, 0,    0, 4, 4, 4, 0,    0, 4, 4, 4, 0,    0, 4, 4, 4, 0,    0, 0, 0, 0, 0};

        IntVolume::iterator i2 = res.begin();
        IntVolume::iterator i2end = res.end();
        const int * p = out5;

        for(; i2 != i2end; ++i2, ++p)
        {
            should( *i2 == *p );
        }
    }


    void labelingTwentySixTest1()
    {
        IntVolume res(vol1);

        should(2 == labelVolume(srcMultiArrayRange(vol1), destMultiArray(res), NeighborCode3DTwentySix()));

        IntVolume::iterator i1 = vol1.begin();
        IntVolume::iterator i1end = vol1.end();
        IntVolume::iterator i2 = res.begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should( *i1 == (*i2 - 1.0) );
        }
    }

    void labelingTwentySixTest2()
    {
        IntVolume res(vol2);

        should(2 == labelVolume(srcMultiArrayRange(vol2), destMultiArray(res), NeighborCode3DTwentySix()));

        IntVolume::iterator i1 = vol2.begin();
        IntVolume::iterator i1end = vol2.end();
        IntVolume::iterator i2 = res.begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should( *i1 == (*i2 - 1.0) );
        }
    }

    void labelingTwentySixTest3()
    {
        IntVolume res(vol4);

        should(2 == labelVolume(srcMultiArrayRange(vol4), destMultiArray(res), NeighborCode3DTwentySix()));

        DoubleVolume::iterator i1 = vol4.begin();
        DoubleVolume::iterator i1end = vol4.end();
        IntVolume::iterator i2 = res.begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should( *i1 == 2-*i2 );
        }
    }

    void labelingTwentySixWithBackgroundTest1()
    {
        IntVolume res(vol5);

        should(2 == labelVolumeWithBackground(srcMultiArrayRange(vol5), destMultiArray(res), NeighborCode3DTwentySix(), 0));

        DoubleVolume::iterator i1 = vol5.begin();
        DoubleVolume::iterator i1end = vol5.end();
        IntVolume::iterator i2 = res.begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should( *i1 == *i2 );
        }
    }

    IntVolume vol1, vol2, vol3;
    DoubleVolume vol4, vol5;
};



struct SimpleAnalysisTestSuite
: public vigra::test_suite
{
    SimpleAnalysisTestSuite()
    : vigra::test_suite("SimpleAnalysisTestSuite")
    {
        add( testCase( &LabelingTest::labelingSixTest1));
        add( testCase( &LabelingTest::labelingSixTest2));
        add( testCase( &LabelingTest::labelingSixTest3));
        add( testCase( &LabelingTest::labelingSixTest4));
        add( testCase( &LabelingTest::labelingSixWithBackgroundTest1));
        add( testCase( &LabelingTest::labelingTwentySixTest1));
        add( testCase( &LabelingTest::labelingTwentySixTest2));
        add( testCase( &LabelingTest::labelingTwentySixTest3));
        add( testCase( &LabelingTest::labelingTwentySixWithBackgroundTest1));
    }
};

int main()
{
    SimpleAnalysisTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;
    return (failed != 0);
}

