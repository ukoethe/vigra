/************************************************************************/
/*                                                                      */
/*       Copyright 2007 by F. Heinrich, B. Seppke, Ullrich Koethe       */
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

#include <iostream>
#include <functional>
#include <cmath>
#include <list>
#include <vigra/unittest.hxx>

#include <vigra/multi_distance.hxx>
#include <vigra/distancetransform.hxx>
#include <vigra/eccentricitytransform.hxx>
#include <vigra/impex.hxx>

#include "test_data.hxx"

using namespace vigra;

struct MultiDistanceTest
{
    typedef vigra::MultiArray<3,int> IntVolume;
    typedef vigra::MultiArray<3,double> DoubleVolume; 
    typedef vigra::MultiArray<2,double> Double2DArray;
    typedef vigra::DImage Image;
    typedef vigra::MultiArrayView<2,Image::value_type> ImageView;
    typedef vigra::TinyVector<int,3> IntVec;

#if 1
    enum { WIDTH    =   15,  // 
           HEIGHT   =   15,  // Volume-Dimensionen
           DEPTH    =   15}; //
#else
    enum { WIDTH    =   4,  // 
           HEIGHT   =   4,  // Volume-Dimensionen
           DEPTH    =   1}; //
#endif

    std::list<std::list<IntVec> > pointslists;
    std::vector<Image> images;
    Double2DArray img2;
    DoubleVolume volume;
    IntVolume shouldVol;

    MultiDistanceTest()
    : images(3, Image(7,7)), img2(Double2DArray::difference_type(7,1)),
      volume(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH)),
      shouldVol(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH))
    {
        std::list<IntVec> temp;
        temp.push_back(IntVec(      0,        0,       0));
        temp.push_back(IntVec(WIDTH-1,        0,       0));
        temp.push_back(IntVec(      0, HEIGHT-1,       0));
        temp.push_back(IntVec(WIDTH-1, HEIGHT-1,       0));
        temp.push_back(IntVec(      0,        0, DEPTH-1));
        temp.push_back(IntVec(WIDTH-1,        0, DEPTH-1));
        temp.push_back(IntVec(      0, HEIGHT-1, DEPTH-1));
        temp.push_back(IntVec(WIDTH-1, HEIGHT-1, DEPTH-1));
        pointslists.push_back(temp);

        temp.clear();
        temp.push_back(IntVec(      0, HEIGHT/2, DEPTH/2));
        temp.push_back(IntVec(WIDTH/2, HEIGHT/2,       0));
        temp.push_back(IntVec(WIDTH/2,        0, DEPTH/2));
        temp.push_back(IntVec(WIDTH-1, HEIGHT/2, DEPTH/2));
        temp.push_back(IntVec(WIDTH/2, HEIGHT/2, DEPTH-1));
        temp.push_back(IntVec(WIDTH/2, HEIGHT-1, DEPTH/2));
        pointslists.push_back(temp);
    

        static const double in[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        {
            Image::ScanOrderIterator i = images[0].begin();
            Image::ScanOrderIterator end = images[0].end();
            Image::Accessor acc = images[0].accessor();
            const double * p = in;

            for(; i != end; ++i, ++p)
            {
                acc.set(*p, i);
            }
        }
        
        static const unsigned char in2[] = {
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1};

        {
            Image::ScanOrderIterator i = images[1].begin();
            Image::ScanOrderIterator end = images[1].end();
            Image::Accessor acc = images[1].accessor();
            const unsigned char * p = in2;

            for(; i != end; ++i, ++p)
            {
                acc.set(*p, i);
            }
        }

        static const unsigned char in3[] = {
            1, 1, 1, 1, 1, 1, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 1, 1, 1, 1, 1, 1};

        {
            Image::ScanOrderIterator i = images[2].begin();
            Image::ScanOrderIterator end = images[2].end();
            Image::Accessor acc = images[2].accessor();
            const unsigned char * p = in3;

            for(; i != end; ++i, ++p)
            {
                acc.set(*p, i);
            }
        }

        static const double in2d[] = {0, 0, 0, 1, 0, 0, 0};
        const double * p=in2d;
        for(Double2DArray::iterator iter=img2.begin(); iter!=img2.end(); ++iter, ++p){
            *iter=*p;
        }
    }

    void testDistanceVolumes()
    {    
        DoubleVolume desired(volume);
        for(std::list<std::list<IntVec> >::iterator list_iter=pointslists.begin(); 
                                          list_iter!=pointslists.end(); ++list_iter)
        {
            IntVec temp;
            for(int z=0; z<DEPTH; ++z)
                for(int y=0; y<HEIGHT; ++y)
                    for(int x=0; x<WIDTH; ++x){
                        temp = IntVec(x,y,z);
                        int tempVal=10000000;
                        for(std::list<IntVec>::iterator iter=(*list_iter).begin(); iter!=(*list_iter).end(); ++iter){
                            if((temp-*iter).squaredMagnitude()<tempVal){
                                tempVal = (temp-*iter).squaredMagnitude();
                            }
                        }
                        desired(x,y,z)=tempVal;
                    }

            for(DoubleVolume::iterator vol_iter = volume.begin(); vol_iter != volume.end(); ++vol_iter)
                *vol_iter=0;
            for(std::list<IntVec>::iterator iter=(*list_iter).begin(); iter!=(*list_iter).end(); ++iter){
                *(volume.traverser_begin()+*iter)=1;
            }

            separableMultiDistSquared(volume, volume, true);
            shouldEqualSequence(volume.begin(),volume.end(),desired.begin());
        }

        typedef MultiArrayShape<3>::type Shape;
        MultiArrayView<3, double> vol(Shape(12,10,35), volume_data);
        
        MultiArray<3, double> res(vol.shape());

        separableMultiDistSquared(srcMultiArrayRange(vol), destMultiArray(res), false);
                
        shouldEqualSequence(res.data(), res.data()+res.elementCount(), ref_dist2);
    }

    void testDistanceAxesPermutation()
    {
        typedef MultiArrayShape<3>::type Shape;
        MultiArrayView<3, double> vol(Shape(12,10,35), volume_data);
        
        MultiArray<3, double> res1(vol.shape()), res2(vol.shape());
        MultiArrayView<3, double, StridedArrayTag> pvol(vol.transpose()), pres2(res2.transpose());
        
        separableMultiDistSquared(vol, res1, true);
        separableMultiDistSquared(pvol, pres2, true);
                
        shouldEqualSequence(res1.data(), res1.data()+res1.elementCount(), res2.data());
        
        separableMultiDistSquared(vol, res1, false);
        separableMultiDistSquared(pvol, pres2, false);
                
        shouldEqualSequence(res1.data(), res1.data()+res1.elementCount(), res2.data());
    }

    void testDistanceVolumesAnisoptopic()
    {    
        double epsilon = 1e-14;
        TinyVector<double, 3> pixelPitch(1.2, 1.0, 2.4);
        
        DoubleVolume desired(volume);
        for(std::list<std::list<IntVec> >::iterator list_iter=pointslists.begin(); 
                                          list_iter!=pointslists.end(); ++list_iter){

            for(DoubleVolume::iterator vol_iter = volume.begin(); vol_iter != volume.end(); ++vol_iter)
                *vol_iter=0;
            for(std::list<IntVec>::iterator iter=(*list_iter).begin(); iter!=(*list_iter).end(); ++iter)
                *(volume.traverser_begin()+*iter)=1;

            IntVec temp;
            for(int z=0; z<DEPTH; ++z)
                for(int y=0; y<HEIGHT; ++y)
                {
                    for(int x=0; x<WIDTH; ++x)
                    {
                        temp = IntVec(x,y,z);
                        double tempVal=10000000.0;
                        for(std::list<IntVec>::iterator iter=(*list_iter).begin(); iter!=(*list_iter).end(); ++iter){
                            double squaredMag = (pixelPitch*(temp-*iter)).squaredMagnitude();
                            if(squaredMag<tempVal){
                                tempVal = squaredMag;
                            }
                        }
                        desired(x,y,z)=tempVal;
                    }
                }


            separableMultiDistSquared(volume, volume, true, pixelPitch);
            shouldEqualSequenceTolerance(volume.begin(),volume.end(),desired.begin(), epsilon);
        }
    }

    void distanceTransform2DCompare()
    {
        for(unsigned int k=0; k<images.size(); ++k)
        {
            Image res(images[k]);
            ImageView img_array(ImageView::difference_type(images[k].width(), images[k].height()), &images[k](0,0));

            distanceTransform(srcImageRange(images[k]), destImage(res), 0.0, 2);

            separableMultiDistance(img_array, img_array, true);

            Image::Iterator i = res.upperLeft();
            Image::Accessor acc = res.accessor();

            int x,y;

            for(y=0; y<7; ++y)
            {
                for(x=0; x<7; ++x)
                {
                    double dist_old = acc(i, vigra::Diff2D(x,y));

                    shouldEqualTolerance(dist_old, img_array(x,y), 1e-7);
                }
            }
        }
    }

    void distanceTest1D()
    {
        vigra::MultiArray<2,double> res(img2);
        
        static const int desired[] = {3, 2, 1, 0, 1, 2, 3};
        separableMultiDistance(img2, res, true);
        shouldEqualSequence(res.begin(), res.end(), desired);
    }
};

struct BoundaryMultiDistanceTest
{
    typedef vigra::MultiArray<3,int> IntVolume;
    typedef vigra::MultiArray<3,double> DoubleVolume; 
    typedef vigra::MultiArray<2,double> Double2DArray;
    typedef vigra::MultiArray<1,double> Double1DArray;
    typedef vigra::DImage Image;
    typedef vigra::MultiArrayView<2,Image::value_type> ImageView;
    typedef vigra::TinyVector<int,3> IntVec;

#if 1
    enum { WIDTH    =   50,  // 
           HEIGHT   =   50,  // Volume-Dimensionen
           DEPTH    =   1}; //
#else
    enum { WIDTH    =   4,  // 
           HEIGHT   =   4,  // Volume-Dimensionen
           DEPTH    =   1}; //
#endif

    std::list<std::list<IntVec> > pointslists;
    std::vector<Image> images;
    Double1DArray img2;
    DoubleVolume volume;
    IntVolume shouldVol;

    BoundaryMultiDistanceTest()
    : images(3, Image(7,7)), img2(Shape1(7)),
      volume(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH)),
      shouldVol(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH))
    {
        std::list<IntVec> temp;
        temp.push_back(IntVec(      0,        0,       0));
        temp.push_back(IntVec(WIDTH-1,        0,       0));
        temp.push_back(IntVec(      0, HEIGHT-1,       0));
        temp.push_back(IntVec(WIDTH-1, HEIGHT-1,       0));
        temp.push_back(IntVec(      0,        0, DEPTH-1));
        temp.push_back(IntVec(WIDTH-1,        0, DEPTH-1));
        temp.push_back(IntVec(      0, HEIGHT-1, DEPTH-1));
        temp.push_back(IntVec(WIDTH-1, HEIGHT-1, DEPTH-1));
        pointslists.push_back(temp);

        temp.clear();
        temp.push_back(IntVec(      0, HEIGHT/2, DEPTH/2));
        temp.push_back(IntVec(WIDTH/2, HEIGHT/2,       0));
        temp.push_back(IntVec(WIDTH/2,        0, DEPTH/2));
        temp.push_back(IntVec(WIDTH-1, HEIGHT/2, DEPTH/2));
        temp.push_back(IntVec(WIDTH/2, HEIGHT/2, DEPTH-1));
        temp.push_back(IntVec(WIDTH/2, HEIGHT-1, DEPTH/2));
        pointslists.push_back(temp);
    

        static const double in[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        {
            Image::ScanOrderIterator i = images[0].begin();
            Image::ScanOrderIterator end = images[0].end();
            Image::Accessor acc = images[0].accessor();
            const double * p = in;

            for(; i != end; ++i, ++p)
            {
                acc.set(*p, i);
            }
        }
        
        static const unsigned char in2[] = {
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1};

        {
            Image::ScanOrderIterator i = images[1].begin();
            Image::ScanOrderIterator end = images[1].end();
            Image::Accessor acc = images[1].accessor();
            const unsigned char * p = in2;

            for(; i != end; ++i, ++p)
            {
                acc.set(*p, i);
            }
        }

        static const unsigned char in3[] = {
            1, 1, 1, 1, 1, 1, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 0, 0, 0, 0, 0, 1,
            1, 1, 1, 1, 1, 1, 1};

        {
            Image::ScanOrderIterator i = images[2].begin();
            Image::ScanOrderIterator end = images[2].end();
            Image::Accessor acc = images[2].accessor();
            const unsigned char * p = in3;

            for(; i != end; ++i, ++p)
            {
                acc.set(*p, i);
            }
        }

        static const double in2d[] = {0, 0, 0, 1, 0, 0, 0};
        const double * p=in2d;
        for(Double1DArray::iterator iter=img2.begin(); iter!=img2.end(); ++iter, ++p){
            *iter=*p;
        }
    }

    void testDistanceVolumes()
    {    
        MultiArrayView<2, double> vol(Shape2(50,50), bndMltDst_data);    
        MultiArray<2, double> res(vol.shape());

        boundaryMultiDistance(vol, res);
        shouldEqualSequenceTolerance(res.begin(), res.end(), bndMltDst_ref, 1e-6);

        boundaryMultiDistance(vol, res, true);
        shouldEqualSequenceTolerance(res.begin(), res.end(), bndMltDstArrayBorder_ref, 1e-6);
    }

    void distanceTest1D()
    {
        Double1DArray res(img2);
        
        static const float desired[] = {2.5, 1.5, 0.5, 0.5, 0.5, 1.5, 2.5};
        boundaryMultiDistance(img2, res);

        shouldEqualSequence(res.begin(), res.end(), desired);
        res = 0.0;

        static const float desired2[] = {0.5, 1.5, 0.5, 0.5, 0.5, 1.5, 0.5};
        boundaryMultiDistance(img2, res, true);

        shouldEqualSequence(res.begin(), res.end(), desired2);
    }
};

struct EccentricityTest
{
    void testEccentricityCenters()
    {
        typedef Shape2 Point;
        typedef Shape3 Point3;
        {
            MultiArray<2, int> labels(Shape2(4,2), 1);
            labels.subarray(Shape2(2,0), Shape2(4,2)) = 2;

            ArrayVector<Point> centers;
            MultiArray<2, float> distances(labels.shape());
            eccentricityTransformOnLabels(labels, distances, centers);

            shouldEqual(centers.size(), 3);
            shouldEqual(centers[1], Point(1,1));
            shouldEqual(centers[2], Point(3,1));

            float ref[] = {1.41421f, 1, 1.41421f, 1, 1, 0, 1, 0 };
            shouldEqualSequenceTolerance(distances.begin(), distances.end(), ref, 1e-5f);
        }
        {
            int image_data_small[100] = {
                1, 2, 3, 3, 3, 3, 2, 2, 4, 4,
                2, 2, 2, 3, 3, 3, 2, 2, 2, 4,
                5, 2, 2, 2, 2, 2, 2, 6, 2, 2,
                5, 5, 2, 2, 2, 2, 6, 6, 2, 2,
                5, 5, 5, 5, 6, 2, 6, 2, 2, 2,
                5, 5, 5, 5, 6, 6, 6, 2, 2, 7,
                5, 5, 2, 2, 6, 6, 2, 2, 2, 7,
                5, 5, 2, 2, 2, 2, 2, 2, 2, 7,
                5, 5, 2, 2, 2, 2, 2, 7, 7, 7,
                5, 5, 2, 2, 7, 7, 7, 7, 7, 7
            };
            MultiArrayView<2, int> labels(Shape2(10,10), image_data_small);

            ArrayVector<Point> centers;
            MultiArray<2, float> distances(labels.shape());
            eccentricityTransformOnLabels(labels, distances, centers);

            shouldEqual(centers.size(), 8);

            Point centers_ref[] = {
                Point(0, 0),
                Point(8, 2),
                Point(3, 1),
                Point(9, 1),
                Point(1, 5),
                Point(6, 4),
                Point(8, 8)
            };
            shouldEqualSequence(centers.begin()+1, centers.end(), centers_ref);

            float dist_ref[] = {
                0.000000f, 8.656855f, 1.414214f, 1.000000f, 1.414214f, 2.414214f, 2.828427f, 2.414214f, 1.414214f, 1.000000f,
                9.242640f, 8.242640f, 7.242641f, 0.000000f, 1.000000f, 2.000000f, 2.414214f, 1.414214f, 1.000000f, 0.000000f,
                3.414214f, 7.828427f, 6.828427f, 5.828427f, 4.828427f, 3.828427f, 2.828427f, 2.414214f, 0.000000f, 1.000000f,
                2.414214f, 2.000000f, 7.242640f, 6.242640f, 5.242640f, 4.242640f, 1.000000f, 1.414214f, 1.000000f, 1.414214f,
                1.414214f, 1.000000f, 1.414214f, 2.414214f, 2.828427f, 5.242640f, 0.000000f, 2.414214f, 2.000000f, 2.414214f,
                1.000000f, 0.000000f, 1.000000f, 2.000000f, 2.414214f, 1.414214f, 1.000000f, 3.414214f, 3.000000f, 3.414214f,
                1.414214f, 1.000000f, 9.656855f, 8.656855f, 2.828427f, 2.414214f, 4.828427f, 4.414214f, 4.000000f, 2.414214f,
                2.414214f, 2.000000f, 9.242641f, 8.242641f, 7.242641f, 6.242641f, 5.828427f, 5.414214f, 5.000000f, 1.414214f,
                3.414214f, 3.000000f, 9.656855f, 8.656855f, 7.656855f, 7.242641f, 6.828427f, 1.000000f, 0.000000f, 1.000000f,
                4.414214f, 4.000000f, 10.071068f, 9.071068f, 4.414214f, 3.414214f, 2.414214f, 1.414214f, 1.000000f, 1.414214f,
            };
            shouldEqualSequenceTolerance(distances.begin(), distances.end(), dist_ref, 1e-5f);
        }
        {
            MultiArrayView<2, unsigned int> labels(Shape2(100, 100), eccTrafo_data);

            ArrayVector<Point> centers;
            MultiArray<2, float> distances(labels.shape());
            eccentricityTransformOnLabels(labels, distances, centers);

            shouldEqual(centers.size(), 98);

            Point centers_ref[98];
            for (int i=0; i<98; ++i) {
                centers_ref[i] = Point(eccTrafo_centers[2*i], eccTrafo_centers[2*i+1]);
            }
            shouldEqualSequence(centers.begin(), centers.end(), centers_ref);

            shouldEqualSequenceTolerance(distances.begin(), distances.end(), eccTrafo_ref, 1e-5f);
        }
        {
            MultiArrayView<3, unsigned int> labels(Shape3(40, 40, 40), eccTrafo_volume);

            ArrayVector<Point3> centers;
            MultiArray<3, float> distances(labels.shape());
            eccentricityTransformOnLabels(labels, distances, centers);

            shouldEqual(centers.size(), 221);

            Point3 centers_ref[221];
            for (int i=0; i<221; ++i) {
                centers_ref[i] = Point3(eccTrafo_volume_centers[3*i], eccTrafo_volume_centers[3*i+1], eccTrafo_volume_centers[3*i+2]);
            }
            shouldEqualSequence(centers.begin(), centers.end(), centers_ref);

            shouldEqualSequenceTolerance(distances.begin(), distances.end(), eccTrafo_volume_ref, 1e-5f);
        }
    }
};


struct DistanceTransformTestSuite
: public vigra::test_suite
{
    DistanceTransformTestSuite()
    : vigra::test_suite("DistanceTransformTestSuite")
    {
        add( testCase( &MultiDistanceTest::testDistanceVolumes));
        add( testCase( &MultiDistanceTest::testDistanceAxesPermutation));
        add( testCase( &MultiDistanceTest::testDistanceVolumesAnisoptopic));
        add( testCase( &MultiDistanceTest::distanceTransform2DCompare));
        add( testCase( &MultiDistanceTest::distanceTest1D));
        add( testCase( &BoundaryMultiDistanceTest::distanceTest1D));
        add( testCase( &BoundaryMultiDistanceTest::testDistanceVolumes));
        add( testCase( &EccentricityTest::testEccentricityCenters));
    }
};

int main(int argc, char ** argv)
{
    DistanceTransformTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;
    return (failed != 0);
}

