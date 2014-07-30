#include <functional>
#include "vigra/unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/openmp_vigra.h"
#include "vigra/multi_array.hxx"
#include "vigra/error.hxx"

using namespace std;
using namespace vigra;

//#define PARA_CHECK;
//#ifdef OPENMP

class ParallelTestFunctor {
    int num_threads;
public:
    ParallelTestFunctor() :
            num_threads(0) {
    }

    //The copy constructor is initialized per thread
    ParallelTestFunctor(ParallelTestFunctor const&) :
            num_threads(omp_get_num_threads()) {
    }

    template<class T>
    T operator()(T const& arg1) const {
        vigra_precondition(num_threads > 0,
                "The current team should not empty");
        return std::cos(arg1) + num_threads;
    }

    template<class T>
    T operator()(T const& arg1, T const& arg2) const {
        vigra_precondition(num_threads > 0,
                "The current team should not empty");
        return arg1 + arg2 + num_threads;
    }

    template<class T>
    T operator()(T const& arg1, T const& arg2, T const& arg3) const {
        vigra_precondition(num_threads > 0,
                "The current team should not empty");
        return arg1 + arg2 + arg3 + num_threads;
    }
};

struct OpenMPWrapperTest {
    typedef vigra::DImage Image;

    OpenMPWrapperTest() :
            img(10000, 10000), mask(Shape2(10000, 10000)) {
        Image::Accessor acc = img.accessor();
        Image::ScanOrderIterator i = img.begin();

        for (int j = 1; i < img.end(); ++i, ++j) {
            acc.set(j, i);
        }

        //First pixel and last pixel are not affected
        mask.init(1);
        mask.begin()[0] = 0;
        mask.begin()[10000 * 10000 - 1] = 0;
    }

    void copyImageTest() {
        Image img1(10000, 10000);

        vigra::omp::copyImage(srcImageRange(img), destImage(img1));

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        for (; i != img.end(); ++i, ++i1) {
            should(acc(i) == acc(i1));
        }
    }

    void combineTwoImagesTest()
    {
        Image img1(10000, 10000);

        vigra::omp::combineTwoImages(srcImageRange(img), srcImage(img), destImage(img1), ParallelTestFunctor());

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        int num_threads = acc(i1) - 2*acc(i);

        for(; i != img.end(); ++i, ++i1)
        {
            should( 2.0*acc(i) + num_threads == acc(i1) );
        }

#ifdef OPENMP
        should( num_threads > 1 );
#endif
    }

    void combineTwoImagesIfTest()
    {
        Image img1(10000, 10000);
        img1 = 27.0; //Could be parallel this !

        vigra::omp::combineTwoImagesIf(srcImageRange(img), srcImage(img), maskImage(mask), destImage(img1), ParallelTestFunctor());

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator i1end = img1.end();
        Image::Accessor acc = img.accessor();

        //First pixel and last pixel are untouched
        should(acc(i1) == 27.0);

        i++;
        i1++;
        i1end--;

        should(acc(i1end) == 27.0);

        int num_threads = acc(i1) - 2*acc(i);

        for(; i1 != i1end; ++i, ++i1)
        {
            should( 2.0*acc(i) + num_threads == acc(i1) );
        }

#ifdef OPENMP
        should( num_threads > 1 );
#endif
    }

    void combineThreeImagesTest()
    {
        Image img1(10000, 10000);
        Image img2(10000, 10000);

        std::plus<Image::value_type> add;

        vigra::omp::combineTwoImages(srcImageRange(img), srcImage(img), destImage(img1), add);

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i, ++i1)
        {
            should(2.0*acc(i) == acc(i1));
        }

        vigra::omp::combineThreeImages(srcImageRange(img), srcImage(img), srcImage(img1), destImage(img2), ParallelTestFunctor());

        i = img.begin();
        Image::ScanOrderIterator i2 = img2.begin();

        int num_threads = acc(i2) - 4*acc(i);

        for(; i != img.end(); ++i, ++i2)
        {
            should(4.0*acc(i) + num_threads == acc(i2));
        }

#ifdef OPENMP
        should( num_threads > 1 );
#endif
    }

    void transformImageTest()
    {
        Image img1(10000, 10000);

        vigra::omp::transformImage(srcImageRange(img), destImage(img1),
                ParallelTestFunctor());

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        int num_threads = acc(i1) - cos(acc(i));

        for(; i != img.end(); ++i, ++i1)
        {
            should(cos(acc(i)) + num_threads == acc(i1));
        }

#ifdef OPENMP
        should( num_threads > 1 );
#endif
    }

    void transformImageIfTest()
    {
        Image img1(10000, 10000);
        img1 = 27.0; //should be parallel?

        vigra::omp::transformImageIf(srcImageRange(img), maskImage(mask), destImage(img1), ParallelTestFunctor());

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator i1end = img1.end();
        Image::Accessor acc = img.accessor();

        //First pixel and last pixel are untouched
        should(acc(i1) == 27.0);

        i++;
        i1++;
        i1end--;

        should(acc(i1end) == 27.0);

        int num_threads = acc(i1) - cos(acc(i));

        for(; i1 != i1end; ++i, ++i1)
        {
            should( cos(acc(i)) + num_threads == acc(i1));
        }

#ifdef OPENMP
        should( num_threads > 1 );
#endif

    }

    Image img;
    vigra::MultiArray<2, unsigned char> mask;
};

struct DistanceTransformTest {
    typedef vigra::DImage Image;

    DistanceTransformTest() :
            img(7, 7) {
        static const double in[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0 };

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator end = img.end();
        Image::Accessor acc = img.accessor();
        const double * p = in;

        for (; i != end; ++i, ++p) {
            acc.set(*p, i);
        }
    }

    void distanceTransformL1Test() {
        Image res(img);
        Image res1(img);

        vigra::omp::distanceTransform(srcImageRange(img), destImage(res), 0.0, 1);

        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x, y;

        for (y = 0; y < 7; ++y) {
            for (x = 0; x < 7; ++x) {
                double dist = acc(i, vigra::Diff2D(x, y));
                double dist1 = std::abs(2.0 - x) + std::abs(2.0 - y);
                double dist2 = std::abs(5.0 - x) + std::abs(5.0 - y);
                double desired = (dist1 < dist2) ? dist1 : dist2;

                shouldEqualTolerance(dist, desired, 1e-14);
            }
        }
    }

    void distanceTransformL2Test() {

        Image res(img);

        vigra::omp::distanceTransform(srcImageRange(img), destImage(res), 0.0, 2);

        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x, y;

        for (y = 0; y < 7; ++y) {
            for (x = 0; x < 7; ++x) {
                double dist = acc(i, vigra::Diff2D(x, y));
                double dist1 = VIGRA_CSTD::sqrt(
                        (2.0 - x) * (2.0 - x) + (2.0 - y) * (2.0 - y));
                double dist2 = VIGRA_CSTD::sqrt(
                        (5.0 - x) * (5.0 - x) + (5.0 - y) * (5.0 - y));
                double desired = (dist1 < dist2) ? dist1 : dist2;

                shouldEqualTolerance(dist, desired, 1e-7);
            }
        }
    }

    Image img;
};

struct OpenMPWrapperTestSuite: public vigra::test_suite {
    OpenMPWrapperTestSuite() :
            vigra::test_suite("OpenMPWrapperTestSuite") {

        add( testCase( &OpenMPWrapperTest::combineTwoImagesTest));
        add( testCase( &OpenMPWrapperTest::combineTwoImagesIfTest));
        add( testCase( &OpenMPWrapperTest::combineThreeImagesTest));
        add( testCase( &OpenMPWrapperTest::transformImageTest));
        add( testCase( &OpenMPWrapperTest::transformImageIfTest));

//        add( testCase( &OpenMPWrapperTest::copyImageTest));
//        add( testCase( &DistanceTransformTest::distanceTransformL1Test));
//        add( testCase( &DistanceTransformTest::distanceTransformL2Test));
    }
};

int main(int argc, char ** argv) {

    OpenMPWrapperTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

//#endif //OpenMP
