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
#include <vigra/vector_distance.hxx>


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
    typedef vigra::MultiArray<3,vigra::TinyVector<int,3> > IntVecVolume;
    typedef vigra::MultiArray<3,vigra::TinyVector<double,3> > DoubleVecVolume;
    typedef vigra::MultiArray<2,vigra::TinyVector<double,2> > DoubleVecImage;


#if 1
    enum { WIDTH    =   15,  // 
           HEIGHT   =   15,  // Volume-Dimensionen
           DEPTH    =   15}; //
#else
    enum { WIDTH    =   4,  // 
           HEIGHT   =   4,  // Volume-Dimensionen
           DEPTH    =   1}; //
#endif

    std::vector<std::vector<IntVec> > pointslists;
    std::vector<Image> images;
    Double2DArray img2;
    DoubleVolume volume;
    IntVolume shouldVol;

    MultiDistanceTest()
    : images(3, Image(7,7)), img2(Double2DArray::difference_type(7,1)),
      volume(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH)),
      shouldVol(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH))
    {
        std::vector<IntVec> temp;
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
        DoubleVolume dt(volume.shape()), desired(volume.shape());
        DoubleVecVolume vecDesired(volume.shape()); 
        for(unsigned k = 0; k<pointslists.size(); ++k)
        {
            DoubleVolume::iterator i = desired.begin();
            for(; i.isValid(); ++i)
            {
                UInt64 minDist = NumericTraits<UInt64>::max();
                int nearest = -1;
                for(unsigned j=0; j<pointslists[k].size(); ++j)
                {
                    UInt64 dist = squaredNorm(pointslists[k][j] - i.point());
                    if(dist < minDist)
                    {
                        minDist = dist;
                        nearest = j;
                    }
                }
                *i = minDist;
                vecDesired[i.point()] = pointslists[k][nearest] - i.point();
            }

            volume = 0.0;
            for(unsigned j=0; j<pointslists[k].size(); ++j)
                volume[pointslists[k][j]] = 1;

            separableMultiDistSquared(volume, dt, true);
            shouldEqualSequence(dt.begin(), dt.end(), desired.begin());

            {
                //test vectorial distance
                using functor::Arg1;
                DoubleVecVolume vecVolume(volume.shape()); 
                separableVectorDistance(volume, vecVolume, true);
                DoubleVolume distVolume(volume.shape());
                transformMultiArray(vecVolume, distVolume, squaredNorm(Arg1()));
                shouldEqualSequence(distVolume.begin(), distVolume.end(), desired.begin());
                // FIXME: this test fails because the nearest point may be ambiguous
                //shouldEqualSequence(vecVolume.begin(), vecVolume.end(), vecDesired.begin());
            }
        }

        typedef MultiArrayShape<3>::type Shape;
        MultiArrayView<3, double> vol(Shape(12,10,35), volume_data);
        
        MultiArray<3, double> res(vol.shape());

        separableMultiDistSquared(vol, res, false);
                
        shouldEqualSequence(res.data(), res.data()+res.elementCount(), ref_dist2);

        {
            //test vectorial distance
            using functor::Arg1;
            DoubleVecVolume vecVolume(vol.shape()); 
            separableVectorDistance(vol, vecVolume, false);
            DoubleVolume distVolume(vol.shape());
            transformMultiArray(vecVolume, distVolume, squaredNorm(Arg1()));
            shouldEqualSequence(distVolume.begin(), distVolume.end(), ref_dist2);
        }
    }

    void testDistanceAxesPermutation()
    {
        using namespace vigra::functor;
        typedef MultiArrayShape<3>::type Shape;
        MultiArrayView<3, double> vol(Shape(12,10,35), volume_data);
        
        MultiArray<3, double> res1(vol.shape()), res2(vol.shape());
        DoubleVecVolume vecVolume(reverse(vol.shape())); 
        
        separableMultiDistSquared(vol, res1, true);
        separableMultiDistSquared(vol.transpose(), res2.transpose(), true);
        shouldEqualSequence(res1.data(), res1.data()+res1.elementCount(), res2.data());

        res2 = 0.0;
        separableVectorDistance(vol.transpose(), vecVolume, true);
        transformMultiArray(vecVolume, res2.transpose(), squaredNorm(Arg1()));
        shouldEqualSequence(res1.data(), res1.data()+res1.elementCount(), res2.data());

        separableMultiDistSquared(vol, res1, false);
        separableMultiDistSquared(vol.transpose(), res2.transpose(), false);
        shouldEqualSequence(res1.data(), res1.data()+res1.elementCount(), res2.data());

        res2 = 0.0;
        separableVectorDistance(vol.transpose(), vecVolume, false);
        transformMultiArray(vecVolume, res2.transpose(), squaredNorm(Arg1()));
        shouldEqualSequence(res1.data(), res1.data()+res1.elementCount(), res2.data());
     }

    void testDistanceVolumesAnisotropic()
    {    
        double epsilon = 1e-14;
        TinyVector<double, 3> pixelPitch(1.2, 1.0, 2.4);
        
        DoubleVolume res(volume.shape()), desired(volume.shape());
        for(unsigned k = 0; k<pointslists.size(); ++k)
        {
            DoubleVolume::iterator i = desired.begin();
            for(; i.isValid(); ++i)
            {
                double minDist = NumericTraits<double>::max();
                int nearest = -1;
                for(unsigned j=0; j<pointslists[k].size(); ++j)
                {
                    double dist = squaredNorm(pixelPitch*(pointslists[k][j] - i.point()));
                    if(dist < minDist)
                    {
                        minDist = dist;
                        nearest = j;
                    }
                }
                *i = minDist;
                //vecDesired[i.point()] = pointslists[k][nearest] - i.point();
            }

            volume = 0.0;
            for(unsigned j=0; j<pointslists[k].size(); ++j)
                volume[pointslists[k][j]] = 1;

            separableMultiDistSquared(volume, res, true, pixelPitch);
            shouldEqualSequenceTolerance(res.begin(), res.end(), desired.begin(), epsilon);

            {
                //test vectorial distance
                using namespace functor;
                DoubleVecVolume vecVolume(volume.shape()); 
                separableVectorDistance(volume, vecVolume, true, pixelPitch);
                DoubleVolume distVolume(volume.shape());
                transformMultiArray(vecVolume, distVolume, squaredNorm(Param(pixelPitch)*Arg1()));
                shouldEqualSequenceTolerance(distVolume.begin(), distVolume.end(), desired.begin(), epsilon);
            }
        }
    }

    void distanceTransform2DCompare()
    {
        for(unsigned int k=0; k<images.size(); ++k)
        {
            Image res_old(images[k]);
            ImageView img_array(Shape2(images[k].width(), images[k].height()), &images[k](0,0));
            MultiArray<2, Image::value_type> res_new(img_array.shape());

            distanceTransform(srcImageRange(images[k]), destImage(res_old), 0.0, 2);

            separableMultiDistance(img_array, res_new, true);
            shouldEqualSequenceTolerance(res_new.begin(), res_new.end(), res_old.data(), 1e-7);

            DoubleVecImage vec_image(img_array.shape());
            separableVectorDistance(img_array, vec_image, true);
            MultiArray<2, double> img_array_dist(img_array.shape());
            using namespace functor;
            transformMultiArray(vec_image, img_array_dist, norm(Arg1()));
            shouldEqualSequenceTolerance(img_array_dist.begin(), img_array_dist.end(), res_old.data(), 1e-7);
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

    std::vector<std::vector<IntVec> > pointslists;
    std::vector<Image> images;
    Double1DArray img2;
    DoubleVolume volume;
    IntVolume shouldVol;

    BoundaryMultiDistanceTest()
    : images(3, Image(7,7)), img2(Shape1(7)),
      volume(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH)),
      shouldVol(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH))
    {
        std::vector<IntVec> temp;
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
        using namespace multi_math;
        MultiArrayView<2, double> vol(Shape2(50,50), bndMltDst_data);    
        MultiArray<2, double> res(vol.shape());

        boundaryMultiDistance(vol, res);
        shouldEqualSequenceTolerance(res.begin(), res.end(), bndMltDst_ref, 1e-6);

        boundaryMultiDistance(vol, res, true);
        shouldEqualSequenceTolerance(res.begin(), res.end(), bndMltDstArrayBorder_ref, 1e-6);

        
        MultiArray<2, double> res2(vol.shape());
        MultiArray<2, TinyVector<double, 2> > res_vec(vol.shape());

        boundaryMultiDistance(vol, res, false, InterpixelBoundary);
        boundaryVectorDistance(vol, res_vec, false, InterpixelBoundary);
        res2 = norm(res_vec);
        shouldEqualSequenceTolerance(res.begin(), res.end(), res2.begin(), 0.25); // FIXME: check this -- 0.25 is a lot
            
        boundaryMultiDistance(vol, res, false, OuterBoundary);
        boundaryVectorDistance(vol, res_vec, false, OuterBoundary);
        res2 = norm(res_vec);
        shouldEqualSequenceTolerance(res.begin(), res.end(), res2.begin(), 1e-15);
            
        boundaryMultiDistance(vol, res, false, InnerBoundary);
        boundaryVectorDistance(vol, res_vec, false, InnerBoundary);
        res2 = norm(res_vec);
        shouldEqualSequenceTolerance(res.begin(), res.end(), res2.begin(), 1e-15);

        boundaryMultiDistance(vol, res, true, InterpixelBoundary);
        boundaryVectorDistance(vol, res_vec, true, InterpixelBoundary);
        res2 = norm(res_vec);
        shouldEqualSequenceTolerance(res.begin(), res.end(), res2.begin(), 0.25); // FIXME: check this -- 0.25 is a lot

        boundaryMultiDistance(vol, res, true, OuterBoundary);
        boundaryVectorDistance(vol, res_vec, true, OuterBoundary);
        res2 = norm(res_vec);
        shouldEqualSequenceTolerance(res.begin(), res.end(), res2.begin(), 1e-15);
            
        boundaryMultiDistance(vol, res, true, InnerBoundary);
        boundaryVectorDistance(vol, res_vec, true, InnerBoundary);
        res2 = norm(res_vec);
        shouldEqualSequenceTolerance(res.begin(), res.end(), res2.begin(), 1e-15);
            
            // FIXME: add tests for alternative boundary definitions
    }

    void distanceTest1D()
    {
        {
            // OuterBoundary
            Double1DArray res(img2.shape());
        
            static const float desired[] = {3, 2, 1, 1, 1, 2, 3};
            boundaryMultiDistance(img2, res, false, OuterBoundary);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
        {
            //InterpixelBoundary
            Double1DArray res(img2.shape());
        
            static const float desired[] = {2.5, 1.5, 0.5, 0.5, 0.5, 1.5, 2.5};
            boundaryMultiDistance(img2, res);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
        {
            // InnerBoundary
            Double1DArray res(img2.shape());
        
            static const float desired[] = {2, 1, 0, 0, 0, 1, 2};
            boundaryMultiDistance(img2, res, false, InnerBoundary);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
        {
            //OuterBoundary and image border
            Double1DArray res(img2.shape());
        
            static const float desired[] = {1, 2, 1, 1, 1, 2, 1};
            boundaryMultiDistance(img2, res, true, OuterBoundary);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
        {
            //InterpixelBoundary and image border
            Double1DArray res(img2.shape());
        
            static const float desired[] = {0.5, 1.5, 0.5, 0.5, 0.5, 1.5, 0.5};
            boundaryMultiDistance(img2, res, true);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
        {
            //InnerBoundary and image border
            Double1DArray res(img2.shape());
        
            static const float desired[] = {0, 1, 0, 0, 0, 1, 0};
            boundaryMultiDistance(img2, res, true, InnerBoundary);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
    }

    void vectorDistanceTest1D()
    {
        typedef TinyVector<double, 1> P;
        {
            // OuterBoundary
            MultiArray<1, P> res(img2.shape());
        
            static const P desired[] = {P(3), P(2), P(1), P(1), P(-1), P(-2), P(-3) };
            boundaryVectorDistance(img2, res, false, OuterBoundary);
            //for(int k=0; k<7; ++k)
            //    std::cerr << res[k] << " ";
            //std::cerr << "\n";
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
        {
            //InterpixelBoundary
            MultiArray<1, P> res(img2.shape());
        
            static const P desired[] = {P(2.5), P(1.5), P(0.5), P(-0.5), P(-0.5), P(-1.5), P(-2.5) };
            boundaryVectorDistance(img2, res, false, InterpixelBoundary);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
        {
            // InnerBoundary
            MultiArray<1, P> res(img2.shape());
        
            static const P desired[] = {P(2), P(1), P(0), P(0), P(0), P(-1), P(-2)};
            boundaryVectorDistance(img2, res, false, InnerBoundary);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
        {
            //OuterBoundary and image border
            MultiArray<1, P> res(img2.shape());
        
            static const P desired[] = {P(-1), P(2), P(1), P(1), P(-1), P(2), P(1) };
            boundaryVectorDistance(img2, res, true, OuterBoundary);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
        {
            //InterpixelBoundary and image border
            MultiArray<1, P> res(img2.shape());
        
            static const P desired[] = {P(-0.5), P(1.5), P(0.5), P(-0.5), P(-0.5), P(1.5), P(0.5)};
            boundaryVectorDistance(img2, res, true, InterpixelBoundary);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
        {
            //InnerBoundary and image border
            MultiArray<1, P> res(img2.shape());
        
            static const P desired[] = {P(0), P(1), P(0), P(0), P(0), P(1), P(0) };
            boundaryVectorDistance(img2, res, true, InnerBoundary);
            shouldEqualSequence(res.begin(), res.end(), desired);
        }
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
        add( testCase( &MultiDistanceTest::testDistanceVolumesAnisotropic));
        add( testCase( &MultiDistanceTest::distanceTransform2DCompare));
        add( testCase( &MultiDistanceTest::distanceTest1D));
        add( testCase( &BoundaryMultiDistanceTest::distanceTest1D));
        add( testCase( &BoundaryMultiDistanceTest::testDistanceVolumes));
        add( testCase( &BoundaryMultiDistanceTest::vectorDistanceTest1D));
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

