#include <iostream>
#include <functional>
#include <cmath>
#include "unittest.hxx"

#include "vigra/watersheds3d.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_convolution.hxx"
#include "list"

#include <stdlib.h>
#include <time.h>


bool sortPredicate(const vigra::TinyVector<int,3>& vec1, const vigra::TinyVector<int,3>& vec2){
    if(vec1[2]!=vec2[2])
        return (vec1[2]<vec2[2]);
    if(vec1[1]!=vec2[1])
        return (vec1[1]<vec2[1]);
    if(vec1[0]!=vec2[0])
        return (vec1[0]<vec2[0]);
    return false;
}

using namespace vigra;
struct Watersheds3dTest
{
    typedef vigra::MultiArray<3,int> IntVolume;
    typedef vigra::MultiArray<3,double> DVolume;
    typedef vigra::TinyVector<int,3> IntVec;

    static const int WIDTH    =   100, // 
                     HEIGHT   =   100, // Volume-Dimensionen
                     DEPTH    =   100;

    DVolume volume;
    IntVolume shouldVol;

    Watersheds3dTest()
    : volume(DVolume::difference_type(WIDTH,HEIGHT,DEPTH)),
      shouldVol(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH))
    {}

    void initTest(){
        /* Seed the random-number generator with current time so that
         * the numbers will be different every time we run.
         */
        srand( (unsigned)time( NULL ) );
    }

    void testDistanceVolumeSix()
    {
        std::list<IntVec> points;

        int count=1+(rand()%10);
        for(int i=0; i<count; ++i){
            bool newpoint=true;
            IntVec temppoint;
            do{
                temppoint = IntVec(rand()%WIDTH,rand()%HEIGHT,rand()%DEPTH);

                for(std::list<IntVec>::iterator iter=points.begin(); iter!=points.end(); ++iter){
                    if( (temppoint-*iter).magnitude() < 2 ){
                        newpoint=false;
                        break;
                    }
                }
            }while(!newpoint);
            points.push_back(temppoint);
        }

        points.sort(sortPredicate);
        
        IntVec temp;
        for(int z=0; z<DEPTH; ++z)
            for(int y=0; y<HEIGHT; ++y)
                for(int x=0; x<WIDTH; ++x){
                    temp = vigra::TinyVector<float,3>(x,y,z);
                    double tempVal=10000000;
                    for(std::list<IntVec>::iterator iter=points.begin(); iter!=points.end(); ++iter){
                        if((double)(temp-*iter).squaredMagnitude()<tempVal){
                            tempVal = (double)(temp-*iter).squaredMagnitude();
                        }
                    }
                    volume(x,y,z)=tempVal;
                }

        //Watersheds3D
        IntVolume labelVolume(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH));

        int max_region_label = vigra::watersheds3DSix( vigra::srcMultiArrayRange(volume),
                                                       vigra::destMultiArray(labelVolume));
        should(max_region_label == count);
    }

    void testDistanceVolumeTwentySix()
    {
        typedef vigra::TinyVector<int,3> IntVec;
        std::list<IntVec> points;

        int count=1+(rand()%10);
        for(int i=0; i<count; ++i){
            bool newpoint=true;
            IntVec temppoint;
            do{
                temppoint = IntVec(rand()%WIDTH,rand()%HEIGHT,rand()%DEPTH);

                for(std::list<IntVec>::iterator iter=points.begin(); iter!=points.end(); ++iter){
                    if( (temppoint-*iter).magnitude() < 2 ){
                        newpoint=false;
                        break;
                    }
                }
            }while(!newpoint);
            points.push_back(temppoint);
        }
        
        points.sort(sortPredicate);
        
        IntVec temp;
        for(int z=0; z<DEPTH; ++z)
            for(int y=0; y<HEIGHT; ++y)
                for(int x=0; x<WIDTH; ++x){
                    temp = vigra::TinyVector<float,3>(x,y,z);
                    double tempVal=10000000;
                    for(std::list<IntVec>::iterator iter=points.begin(); iter!=points.end(); ++iter){
                        if((double)(temp-*iter).squaredMagnitude()<tempVal){
                                tempVal = (double)(temp-*iter).squaredMagnitude();
                        }
                    }
                    volume(x,y,z)=tempVal;
                }

        //Watersheds3D
        IntVolume labelVolume(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH));

        int max_region_label = vigra::watersheds3DTwentySix( vigra::srcMultiArrayRange(volume),
                                                             vigra::destMultiArray(labelVolume));
        should(max_region_label == count);
    }

    void testWatersheds3dSix1()
    {
        static const int data[] = {
            6, 5, 4, 5,    
            7, 6, 5, 6,    
           10, 9, 8, 9,  
            6, 5, 4, 5,

            5, 3, 2, 3,  
            6, 4, 3, 4,    
            9, 7, 6, 7,  
            5, 3, 2, 3,

            4, 2, 1, 2,  
            5, 3, 2, 3,    
            8, 6, 5, 6,  
            4, 2, 1, 2,

            5, 3, 2, 3,  
            6, 4, 3, 4,    
            9, 7, 6, 7,  
            5, 3, 2, 3
        };

        IntVolume vol(IntVolume::difference_type(4,4,4));
        const int *i=data;
        for(IntVolume::iterator iter=vol.begin(); iter!=vol.end(); ++iter, ++i){
            *iter=*i;
        }

        static const int desired[] = {
            1, 1, 1, 1,  
            1, 1, 1, 1,    
            1, 1, 2, 2, 
            2, 2, 2, 2,

            1, 1, 1, 1, 
            1, 1, 1, 1,    
            2, 2, 2, 2, 
            2, 2, 2, 2,

            1, 1, 1, 1, 
            1, 1, 1, 1,    
            2, 2, 2, 2,  
            2, 2, 2, 2,

            1, 1, 1, 1, 
            1, 1, 1, 1,    
            2, 2, 2, 2, 
            2, 2, 2, 2
        };

        IntVolume labelVolume(IntVolume::difference_type(4,4,4));

        int max_region_label = vigra::watersheds3DTwentySix( vigra::srcMultiArrayRange(vol),
                                                             vigra::destMultiArray(labelVolume));

        shouldEqualSequence(labelVolume.begin(), labelVolume.end(), desired);
    }

    void testWatersheds3dSix2()
    {
        static const int data[] = {
            3, 2, 3,    
            7, 6, 7,      
            2, 3, 4,

            2, 1, 2, 
            6, 5, 6,
            3, 2, 3,

            3, 2, 3, 
            7, 6, 7,    
            4, 3, 4
        };

        IntVolume vol(IntVolume::difference_type(3,3,3));
        const int *i=data;
        for(IntVolume::iterator iter=vol.begin(); iter!=vol.end(); ++iter, ++i){
            *iter=*i;
        }

        static const int desired[] = {
            1, 1, 1,
            2, 1, 1,
            2, 2, 2,

            1, 1, 1,
            1, 1, 1,
            2, 3, 3,

            1, 1, 1,
            1, 1, 1,
            2, 3, 3
        };

        IntVolume labelVolume(IntVolume::difference_type(3,3,3));

        int max_region_label = vigra::watersheds3DSix( vigra::srcMultiArrayRange(vol),
                                                       vigra::destMultiArray(labelVolume));

        shouldEqualSequence(labelVolume.begin(), labelVolume.end(), desired);
    }

    void testWatersheds3dGradient1()
    {
        int w=11,h=11,d=11;

        DVolume vol(DVolume::difference_type(w,h,d));

        for(int z=0; z<d; ++z)
            for(int y=0; y<h; ++y)
                for(int x=0; x<w; ++x)
                    vol(x,y,z)=(double)(z<d/2.0);
        
        vigra::MultiArray<3, vigra::TinyVector<float,3> > dest(IntVolume::difference_type(w,h,d));
        vigra::gaussianGradientMultiArray(srcMultiArrayRange(vol),destMultiArray(dest),.5);

        for(int z=0; z<d; ++z)
            for(int y=0; y<h; ++y)
                for(int x=0; x<w; ++x)
                    vol(x,y,z) = norm(dest(x,y,z));


        IntVolume labelVolume(IntVolume::difference_type(w,h,d));
        IntVolume labelVolume26(IntVolume::difference_type(w,h,d));

        int max_region_label = vigra::watersheds3DSix( vigra::srcMultiArrayRange(vol),
                                                       vigra::destMultiArray(labelVolume));
        int max_region_label26 = vigra::watersheds3DTwentySix( vigra::srcMultiArrayRange(vol),
                                                               vigra::destMultiArray(labelVolume26));
        should(max_region_label==max_region_label26);

        for(int z=0; z<d; ++z)
            for(int y=0; y<h; ++y)
                for(int x=0; x<w; ++x){
                    if(z<d/2.0){
                        shouldEqual( labelVolume(x,y,z), 1 );
                        shouldEqual( labelVolume26(x,y,z), 1 );
                    }
                    else{
                        shouldEqual( labelVolume(x,y,z), 2 );
                        shouldEqual( labelVolume26(x,y,z), 2 );
                    }
                }
    }

    void testWatersheds3dGradient2()
    {
        int w=11,h=11,d=11;

        DVolume vol(DVolume::difference_type(w,h,d));
        DVolume src(DVolume::difference_type(w,h,d));

        for(int z=0; z<d; ++z){
            for(int y=0; y<h; ++y){
                for(int x=0; x<w; ++x){
                    vol(x,y,z)=( ((z>=d/2.0)<<2) | ((y>=h/2.0)<<1) | ((x>=w/2.0)<<0) );
                }
            }
        }
        
        vigra::MultiArray<3, vigra::TinyVector<float,3> > dest(IntVolume::difference_type(w,h,d));
        vigra::gaussianGradientMultiArray(srcMultiArrayRange(vol),destMultiArray(dest),.5);

        for(int z=0; z<d; ++z)
            for(int y=0; y<h; ++y)
                for(int x=0; x<w; ++x)
                    src(x,y,z) = norm(dest(x,y,z));


        IntVolume labelVolume(IntVolume::difference_type(w,h,d));
        IntVolume labelVolume26(IntVolume::difference_type(w,h,d));

        int max_region_label = vigra::watersheds3DSix( vigra::srcMultiArrayRange(src),
                                                       vigra::destMultiArray(labelVolume));

        for(int z=0; z<d; ++z){
            for(int y=0; y<h; ++y){
                for(int x=0; x<w; ++x){
                    shouldEqual( labelVolume(x,y,z)-1, (int)vol(x,y,z) );
                }
            }
        }
    }


};


struct SimpleAnalysisTestSuite
: public vigra::test_suite
{
    SimpleAnalysisTestSuite()
    : vigra::test_suite("SimpleAnalysisTestSuite")
    {
        add( testCase( &Watersheds3dTest::initTest));
        add( testCase( &Watersheds3dTest::testDistanceVolumeSix));
        add( testCase( &Watersheds3dTest::testDistanceVolumeSix));
        add( testCase( &Watersheds3dTest::testDistanceVolumeSix));
        add( testCase( &Watersheds3dTest::testDistanceVolumeSix));
        add( testCase( &Watersheds3dTest::testDistanceVolumeTwentySix));
        add( testCase( &Watersheds3dTest::testDistanceVolumeTwentySix));
        add( testCase( &Watersheds3dTest::testDistanceVolumeTwentySix));
        add( testCase( &Watersheds3dTest::testDistanceVolumeTwentySix));
        add( testCase( &Watersheds3dTest::testWatersheds3dSix1));
        add( testCase( &Watersheds3dTest::testWatersheds3dSix2));
        add( testCase( &Watersheds3dTest::testWatersheds3dGradient1));
        add( testCase( &Watersheds3dTest::testWatersheds3dGradient2));

    }
};

int main()
{
    SimpleAnalysisTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;
    return (failed != 0);
}

