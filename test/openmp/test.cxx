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

#include <functional>
#include "vigra/unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/openmp_vigra.h"
#include "vigra/multi_array.hxx"
#include "vigra/error.hxx"
#include "vigra/timing.hxx"

using namespace std;
using namespace vigra;

#define IMAGE_WIDTH 1000
#define IMAGE_HEIGHT 1000
#define DEFAULT_PIXEL_VALUE -27.0


class ParallelTestFunctor 
{
    int num_threads, thread_id;
public:
    ParallelTestFunctor() 
    : num_threads(1)
    , thread_id(0)
    {
    }

#ifdef OPENMP
    //This constructor is initialized per thread
    ParallelTestFunctor(ParallelTestFunctor const&) 
    : num_threads(omp_get_num_threads())
    , thread_id(omp_get_thread_num())
    {}
#endif

    template<class T>
    TinyVector<T, 3> operator()(T const& arg1) const 
    {
        return TinyVector<T, 3>(std::cos(arg1), num_threads, thread_id);
    }

    template<class T>
    TinyVector<T, 3> operator()(T const& arg1, T const& arg2) const 
    {
        return TinyVector<T, 3>(arg1 + arg2, num_threads, thread_id);
    }

    template<class T>
    TinyVector<T, 3> operator()(T const& arg1, T const& arg2, T const& arg3) const
    {
        return TinyVector<T, 3>(arg1 + arg2 + arg3, num_threads, thread_id);
    }
};

struct OpenMPWrapperTest {
    typedef vigra::DImage Image;

    OpenMPWrapperTest() 
    : img(IMAGE_WIDTH, IMAGE_HEIGHT)
    , mask(Shape2(IMAGE_WIDTH, IMAGE_HEIGHT))
#ifdef OPENMP
    , max_threads(omp_get_num_procs() / 2)
#else
    , max_threads(1)
#endif
    {
#ifdef OPENMP
        omp_set_num_threads(max_threads);
#endif
        Image::Accessor acc = img.accessor();
        Image::ScanOrderIterator i = img.begin();

        for (int j = 1; i < img.end(); ++i, ++j) {
            acc.set(j, i);
        }

        //First pixel and last pixel are not affected
        mask.init(1);
        mask.begin()[0] = 0;
        mask.begin()[IMAGE_WIDTH * IMAGE_HEIGHT - 1] = 0;
    }

    void copyImageTest() 
    {
        Image img1(IMAGE_WIDTH, IMAGE_HEIGHT);

        int num_threads = vigra::omp::copyImage(srcImageRange(img), destImage(img1));

        should( num_threads <= max_threads);
#ifdef OPENMP
        should( num_threads > 1);
#endif

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        for (; i != img.end(); ++i, ++i1) {
            should(acc(i) == acc(i1));
        }
    }

    void combineTwoImagesTest()
    {
        typedef MultiArray<2, TinyVector<double, 3> > Image2;
        Image2 img1(Shape2(IMAGE_WIDTH, IMAGE_HEIGHT));

        int num_threads = vigra::omp::combineTwoImages(srcImageRange(img), srcImage(img), destImage(img1), ParallelTestFunctor());

        should( num_threads <= max_threads);
#ifdef OPENMP
        should( num_threads > 1);
#endif

        Image::ScanOrderIterator i = img.begin();
        Image2::iterator i1 = img1.begin();
        Image2::iterator i1end = img1.end();

        std::vector<int> count_per_thread(num_threads, 0);

        for(; i1 != i1end; ++i, ++i1)
        {
            shouldEqual( (*i)*2.0, (*i1)[0]);
            shouldEqual( (double)num_threads, (*i1)[1]);
            ++count_per_thread[int((*i1)[2])];
        }

        int active_threads = 0, count = 0, 
            m = NumericTraits<int>::max(), M = NumericTraits<int>::min();
        for(int k=0; k<num_threads; ++k)
        {
            if(count_per_thread[k] > 0)
                ++active_threads;
            count += count_per_thread[k];
            m = std::min(m, count_per_thread[k]);
            M = std::max(M, count_per_thread[k]);
        }

        shouldEqual(count, img1.size());
        std::cerr << "combineTwoImages() active threads: " << active_threads << " of " << num_threads << "\n";
        std::cerr << "                       load range: " << m << " to " << M << "\n";
    }

    void combineTwoImagesIfTest()
    {
        typedef MultiArray<2, TinyVector<double, 3> > Image2;
        Image2 img1(Shape2(IMAGE_WIDTH, IMAGE_HEIGHT), TinyVector<double, 3>(DEFAULT_PIXEL_VALUE));

        int num_threads = vigra::omp::combineTwoImagesIf(srcImageRange(img), srcImage(img), maskImage(mask), destImage(img1), ParallelTestFunctor());

        should( num_threads <= max_threads);
#ifdef OPENMP
        should( num_threads > 1);
#endif

        Image::ScanOrderIterator i = img.begin();
        Image2::iterator i1 = img1.begin();
        Image2::iterator i1end = img1.end();

        //First pixel and last pixel are preserved
        shouldEqual((*i1)[0], DEFAULT_PIXEL_VALUE);

        i++;
        i1++;
        i1end--;

        shouldEqual((*i1end)[0], DEFAULT_PIXEL_VALUE);

        std::vector<int> count_per_thread(num_threads, 0);
        count_per_thread[0] += 2;

        for(; i1 != i1end; ++i, ++i1)
        {
            shouldEqual( (*i)*2.0, (*i1)[0]);
            shouldEqual( (double)num_threads, (*i1)[1]);
            ++count_per_thread[int((*i1)[2])];
        }

        int active_threads = 0, count = 0, 
            m = NumericTraits<int>::max(), M = NumericTraits<int>::min();
        for(int k=0; k<num_threads; ++k)
        {
            if(count_per_thread[k] > 0)
                ++active_threads;
            count += count_per_thread[k];
            m = std::min(m, count_per_thread[k]);
            M = std::max(M, count_per_thread[k]);
        }

        shouldEqual(count, img1.size());
        std::cerr << "combineTwoImagesIf() active threads: " << active_threads << " of " << num_threads << "\n";
        std::cerr << "                         load range: " << m << " to " << M << "\n";
    }

    void combineThreeImagesTest()
    {
        typedef MultiArray<2, TinyVector<double, 3> > Image2;
        Image2 img1(Shape2(IMAGE_WIDTH, IMAGE_HEIGHT), TinyVector<double, 3>(DEFAULT_PIXEL_VALUE));

        int num_threads = vigra::omp::combineThreeImages(srcImageRange(img), srcImage(img), srcImage(img), destImage(img1), ParallelTestFunctor());

        should( num_threads <= max_threads);
#ifdef OPENMP
        should( num_threads > 1);
#endif
        Image::ScanOrderIterator i = img.begin();
        Image2::iterator i1 = img1.begin();

        std::vector<int> count_per_thread(num_threads, 0);

        for(; i != img.end(); ++i, ++i1)
        {
            shouldEqual( (*i)*3.0, (*i1)[0]);
            shouldEqual( (double)num_threads, (*i1)[1]);
            ++count_per_thread[int((*i1)[2])];
        }

        int active_threads = 0, count = 0, 
            m = NumericTraits<int>::max(), M = NumericTraits<int>::min();
        for(int k=0; k<num_threads; ++k)
        {
            if(count_per_thread[k] > 0)
                ++active_threads;
            count += count_per_thread[k];
            m = std::min(m, count_per_thread[k]);
            M = std::max(M, count_per_thread[k]);
        }

        shouldEqual(count, img1.size());
        std::cerr << "combineThreeImages() active threads: " << active_threads << " of " << num_threads << "\n";
        std::cerr << "                         load range: " << m << " to " << M << "\n";
    }

    void transformImageTest()
    {
        typedef MultiArray<2, TinyVector<double, 3> > Image2;
        Image2 img1(Shape2(IMAGE_WIDTH, IMAGE_HEIGHT));

        int num_threads = vigra::omp::transformImage(srcImageRange(img), destImage(img1), ParallelTestFunctor());

        should( num_threads <= max_threads);
#ifdef OPENMP
        should( num_threads > 1);
#endif

        Image::ScanOrderIterator i = img.begin();
        Image2::iterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        std::vector<int> count_per_thread(num_threads, 0);

        for(; i != img.end(); ++i, ++i1)
        {
            shouldEqual( cos(*i), (*i1)[0]);
            shouldEqual( (double)num_threads, (*i1)[1]);
            ++count_per_thread[int((*i1)[2])];
        }

        int active_threads = 0, count = 0, 
            m = NumericTraits<int>::max(), M = NumericTraits<int>::min();
        for(int k=0; k<num_threads; ++k)
        {
            if(count_per_thread[k] > 0)
                ++active_threads;
            count += count_per_thread[k];
            m = std::min(m, count_per_thread[k]);
            M = std::max(M, count_per_thread[k]);
        }

        shouldEqual(count, img1.size());
        std::cerr << "transformImage() active threads: " << active_threads << " of " << num_threads << "\n";
        std::cerr << "                     load range: " << m << " to " << M << "\n";
    }

    void transformImageIfTest()
    {
        typedef MultiArray<2, TinyVector<double, 3> > Image2;
        Image2 img1(Shape2(IMAGE_WIDTH, IMAGE_HEIGHT), TinyVector<double, 3>(DEFAULT_PIXEL_VALUE));

        int num_threads = vigra::omp::transformImageIf(srcImageRange(img), maskImage(mask), destImage(img1), ParallelTestFunctor());

        should( num_threads <= max_threads);
#ifdef OPENMP
        should( num_threads > 1);
#endif

        Image::ScanOrderIterator i = img.begin();
        Image2::iterator i1 = img1.begin();
        Image2::iterator i1end = img1.end();

        //First pixel and last pixel are preserved
        shouldEqual((*i1)[0], DEFAULT_PIXEL_VALUE);

        i++;
        i1++;
        i1end--;

        shouldEqual((*i1end)[0], DEFAULT_PIXEL_VALUE);

        std::vector<int> count_per_thread(num_threads, 0);
        count_per_thread[0] += 2;

        for(; i1 != i1end; ++i, ++i1)
        {
            shouldEqual( cos(*i), (*i1)[0]);
            shouldEqual( (double)num_threads, (*i1)[1]);
            ++count_per_thread[int((*i1)[2])];
        }

        int active_threads = 0, count = 0, 
            m = NumericTraits<int>::max(), M = NumericTraits<int>::min();
        for(int k=0; k<num_threads; ++k)
        {
            if(count_per_thread[k] > 0)
                ++active_threads;
            count += count_per_thread[k];
            m = std::min(m, count_per_thread[k]);
            M = std::max(M, count_per_thread[k]);
        }

        shouldEqual(count, img1.size());
        std::cerr << "transformImageIf() active threads: " << active_threads << " of " << num_threads << "\n";
        std::cerr << "                       load range: " << m << " to " << M << "\n";
    }

    Image img;
    vigra::MultiArray<2, unsigned char> mask;
    int max_threads;
};

struct DistanceTransformTest {
    typedef vigra::DImage Image;

    DistanceTransformTest() 
    : img(7, 7) 
#ifdef OPENMP
    , max_threads(omp_get_num_procs() / 2)
#else
    , max_threads(1)
#endif
    {
#ifdef OPENMP
        omp_set_num_threads(max_threads);
#endif
        static const double in[] = {
        		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        		0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

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

        int num_threads = vigra::omp::distanceTransform(srcImageRange(img), destImage(res), 0.0, 1);

        should( num_threads <= max_threads);
#ifdef OPENMP
        should( num_threads > 1);
#endif
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

        int num_threads = vigra::omp::distanceTransform(srcImageRange(img), destImage(res), 0.0, 2);

        should( num_threads <= max_threads);
#ifdef OPENMP
        should( num_threads > 1);
#endif
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
    int max_threads;
};

struct OpenMPWrapperTestSuite: public vigra::test_suite {
    OpenMPWrapperTestSuite() :
            vigra::test_suite("OpenMPWrapperTestSuite") {

        add( testCase( &OpenMPWrapperTest::copyImageTest));
        add( testCase( &OpenMPWrapperTest::combineTwoImagesTest));
        add( testCase( &OpenMPWrapperTest::combineTwoImagesIfTest));
        add( testCase( &OpenMPWrapperTest::combineThreeImagesTest));
        add( testCase( &OpenMPWrapperTest::transformImageTest));
        add( testCase( &OpenMPWrapperTest::transformImageIfTest));

        add( testCase( &DistanceTransformTest::distanceTransformL1Test));
        add( testCase( &DistanceTransformTest::distanceTransformL2Test));
    }
};

int main(int argc, char ** argv) {

    OpenMPWrapperTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
