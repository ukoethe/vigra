/************************************************************************/
/*                                                                      */
/*       Copyright 2004 by F. Heinrich, B. Seppke, Ullrich Koethe       */
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
#include "vigra/unittest.hxx"

#include "vigra/watersheds3d.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_watersheds.hxx"
#include "list"

#include <stdlib.h>
#include <time.h>


using namespace vigra;

struct Watersheds3dTest
{
    typedef MultiArray<3,int> IntVolume;
    typedef MultiArray<3,double> DVolume;
    typedef TinyVector<int,3> IntVec;

    enum { WIDTH    =   100, // 
           HEIGHT   =   100, // Volume-Dimensionen
           DEPTH    =   100 }; //

    DVolume volume;
    IntVolume shouldVol;
    std::list<std::list<IntVec> > pointslists;

    Watersheds3dTest()
    : volume(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH)),
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
        temp.push_back(IntVec(      0,        3,       4));
        temp.push_back(IntVec(      0,        3,       6));
        temp.push_back(IntVec(      0,        5,       4));
        temp.push_back(IntVec(      0,        5,       6));
        pointslists.push_back(temp);

        temp.clear();
        temp.push_back(IntVec(     80,       98,      77));
        temp.push_back(IntVec(     33,        8,      97));
        temp.push_back(IntVec(     93,       85,      30));
        temp.push_back(IntVec(     73,       62,      43));
        temp.push_back(IntVec(     93,       84,      14));
        temp.push_back(IntVec(     87,       31,      61));
        temp.push_back(IntVec(     44,       12,      20));
        temp.push_back(IntVec(      1,       29,      86));
        temp.push_back(IntVec(      1,       65,      87));
        temp.push_back(IntVec(     26,       40,      71));
        pointslists.push_back(temp);

        temp.clear();
        temp.push_back(IntVec(      0, HEIGHT/2, DEPTH/2));
        temp.push_back(IntVec(WIDTH/2, HEIGHT/2,       0));
        temp.push_back(IntVec(WIDTH/2,        0, DEPTH/2));
        temp.push_back(IntVec(WIDTH-1, HEIGHT/2, DEPTH/2));
        temp.push_back(IntVec(WIDTH/2, HEIGHT/2, DEPTH-1));
        temp.push_back(IntVec(WIDTH/2, HEIGHT-1, DEPTH/2));
        pointslists.push_back(temp);
    }

    void testDistanceVolumesSix()
    {    
        for(std::list<std::list<IntVec> >::iterator list_iter=pointslists.begin(); list_iter!=pointslists.end(); ++list_iter){
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
                        volume(x,y,z)=tempVal;
                    }

            //Watersheds3D
            IntVolume labelVolume(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH));

            int max_region_label = watersheds3DSix( srcMultiArrayRange(volume),
                                                        destMultiArray(labelVolume));
            should(max_region_label == (int)(*list_iter).size());
        }
    }

    void testDistanceVolumesTwentySix()
    {
        for(std::list<std::list<IntVec> >::iterator list_iter=pointslists.begin(); list_iter!=pointslists.end(); ++list_iter){
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
                        volume(x,y,z)=tempVal;
                    }

            //Watersheds3D
            IntVolume labelVolume(IntVolume::difference_type(WIDTH,HEIGHT,DEPTH));

            int max_region_label = watersheds3DTwentySix( srcMultiArrayRange(volume),
                                                                destMultiArray(labelVolume));
            should(max_region_label == (int)(*list_iter).size());
        }
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
            2, 2, 2, 2, 
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

        IntVolume labelVolume(vol.shape()), labelVolume2(vol.shape());

        int max_region_label = watersheds3DSix( srcMultiArrayRange(vol),
                                                       destMultiArray(labelVolume));
        
        shouldEqual(max_region_label, 2);
        
        //int c=1;
        //for(IntVolume::iterator iter=labelVolume.begin(); iter!=labelVolume.end(); ++iter, ++c){
            //std::cerr << *iter << ", ";
            //if(c%4==0) std::cerr << std::endl;
            //if(c%16==0) std::cerr << std::endl;
        //}

        shouldEqualSequence(labelVolume.begin(), labelVolume.end(), desired);

        should(2 == watershedsMultiArray(vol, labelVolume2, DirectNeighborhood, WatershedOptions().unionFind()));
        should(labelVolume == labelVolume2);
    }

    void testWatersheds3dSix2()
    {
        static const double data[] = {
            3.0, 2.0, 3.0,    
            7.0, 6.0, 7.0,      
            2.0, 3.0, 4.0,

            2.0, 1.0, 2.0, 
            6.0, 5.0, 6.0,
            3.0, 2.1, 3.1,

            3.0, 2.0, 3.0, 
            7.0, 6.0, 7.0,    
            4.0, 3.1, 4.0
        };

        DVolume vol(Shape3(3,3,3));
        const double *i=data;
        for(DVolume::iterator iter=vol.begin(); iter!=vol.end(); ++iter, ++i){
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

        IntVolume labelVolume(vol.shape()), labelVolume2(vol.shape());

        int max_region_label = watersheds3DSix( vol, labelVolume);

        shouldEqual(max_region_label, 3);
        shouldEqualSequence(labelVolume.begin(), labelVolume.end(), desired);

        should(3 == watershedsMultiArray(vol, labelVolume2, DirectNeighborhood, WatershedOptions().unionFind()));
        should(labelVolume == labelVolume2);
    }

    void testWatersheds3dGradient1()
    {
        int w=4,h=4,d=8;

        IntVolume vol(IntVolume::difference_type(w,h,d));

        static const int data[] = {
            0, 0, 0, 0, 
            0, 0, 0, 0, 
            0, 0, 0, 0, 
            0, 0, 0, 0,

            0, 0, 0, 0, 
            0, 0, 0, 0, 
            0, 0, 0, 0, 
            0, 0, 0, 0,

            4, 4, 4, 4, 
            4, 4, 4, 4, 
            4, 4, 4, 4, 
            4, 4, 4, 4,

            9, 9, 9, 9,  
            9, 9, 9, 9,    
            9, 9, 9, 9,  
            9, 9, 9, 9,

            9, 9, 9, 9,  
            9, 9, 9, 9,    
            9, 9, 9, 9,  
            9, 9, 9, 9,

            4, 4, 4, 4, 
            4, 4, 4, 4, 
            4, 4, 4, 4, 
            4, 4, 4, 4,

            0, 0, 0, 0, 
            0, 0, 0, 0, 
            0, 0, 0, 0, 
            0, 0, 0, 0,

            0, 0, 0, 0, 
            0, 0, 0, 0, 
            0, 0, 0, 0, 
            0, 0, 0, 0
        };

        const int *i=data;
        for(IntVolume::iterator iter=vol.begin(); iter!=vol.end(); ++iter, ++i){
            *iter=*i;
        }

 
        IntVolume labelVolume(vol.shape()), labelVolume2(vol.shape());

        int max_region_label = watersheds3DSix( srcMultiArrayRange(vol),
                                                       destMultiArray(labelVolume));
        shouldEqual(2, max_region_label);
        using namespace multi_math;
        should(all(labelVolume.subarray(Shape3(0), Shape3(w,h,d/2)) == 1));
        should(all(labelVolume.subarray(Shape3(w,h,d/2), Shape3(w,h,d)) == 2));

        max_region_label = watersheds3DTwentySix( srcMultiArrayRange(vol),
                                                  destMultiArray(labelVolume2));

        shouldEqual(2, max_region_label);        
        should(all(labelVolume == labelVolume2));

        labelVolume2.init(0);
        should(2 == watershedsMultiArray(vol, labelVolume2, DirectNeighborhood, WatershedOptions().unionFind()));
        should(all(labelVolume == labelVolume2));

        labelVolume2.init(0);
        should(2 == watershedsMultiArray(vol, labelVolume2, IndirectNeighborhood, WatershedOptions().unionFind()));
        should(all(labelVolume == labelVolume2));
    }

    void testWatersheds3dGradient2()
    {
        int w=8,h=8,d=8;

        IntVolume vol(IntVolume::difference_type(w,h,d));

        static const int data[] = {

            0, 0, 4, 9, 9, 4, 0, 0,  
            0, 0, 4, 9, 9, 4, 0, 0,  
            4, 4, 4, 9, 9, 4, 4, 4,  
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,  
            4, 4, 4, 9, 9, 4, 4, 4,  
            0, 0, 4, 9, 9, 4, 0, 0,  
            0, 0, 4, 9, 9, 4, 0, 0,

            0, 0, 4, 9, 9, 4, 0, 0,  
            0, 0, 4, 9, 9, 4, 0, 0,  
            4, 4, 4, 9, 9, 4, 4, 4,   
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,   
            4, 4, 4, 9, 9, 4, 4, 4,  
            0, 0, 4, 9, 9, 4, 0, 0,  
            0, 0, 4, 9, 9, 4, 0, 0,

            4, 4, 4, 9, 9, 4, 4, 4,  
            4, 4, 4, 9, 9, 4, 4, 4,  
            4, 4, 4, 9, 9, 4, 4, 4,  
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,   
            4, 4, 4, 9, 9, 4, 4, 4,  
            4, 4, 4, 9, 9, 4, 4, 4,  
            4, 4, 4, 9, 9, 4, 4, 4,

            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,  
           10,10,10,12,12,10,10,10,  
           10,10,10,12,12,10,10,10, 
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,

            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,  
           10,10,10,12,12,10,10,10,  
           10,10,10,12,12,10,10,10, 
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,

            4, 4, 4, 9, 9, 4, 4, 4,  
            4, 4, 4, 9, 9, 4, 4, 4,  
            4, 4, 4, 9, 9, 4, 4, 4,  
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,   
            4, 4, 4, 9, 9, 4, 4, 4,  
            4, 4, 4, 9, 9, 4, 4, 4,  
            4, 4, 4, 9, 9, 4, 4, 4,

            0, 0, 4, 9, 9, 4, 0, 0,  
            0, 0, 4, 9, 9, 4, 0, 0,  
            4, 4, 4, 9, 9, 4, 4, 4,  
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,   
            4, 4, 4, 9, 9, 4, 4, 4,  
            0, 0, 4, 9, 9, 4, 0, 0,  
            0, 0, 4, 9, 9, 4, 0, 0,

            0, 0, 4, 9, 9, 4, 0, 0,  
            0, 0, 4, 9, 9, 4, 0, 0,  
            4, 4, 4, 9, 9, 4, 4, 4,  
            9, 9, 9,10,10, 9, 9, 9,
            9, 9, 9,10,10, 9, 9, 9,   
            4, 4, 4, 9, 9, 4, 4, 4,  
            0, 0, 4, 9, 9, 4, 0, 0,  
            0, 0, 4, 9, 9, 4, 0, 0

        };

        const int *i=data;
        for(IntVolume::iterator iter=vol.begin(); iter!=vol.end(); ++iter, ++i){
            *iter=*i;
        }

        IntVolume labelVolume(Shape3(w,h,d)), labelVolume2(Shape3(w,h,d));

        int max_region_label = watersheds3DSix( srcMultiArrayRange(vol),
                                                       destMultiArray(labelVolume));
        shouldEqual(8, max_region_label);

        // each octand should be one region
        for(int z=0; z<d; ++z){
            for(int y=0; y<h; ++y){
                for(int x=0; x<w; ++x){
                    shouldEqual( labelVolume(x,y,z)-1, ( ((z>=d/2.0)<<2) | ((y>=h/2.0)<<1) |((x>=w/2.0)<<0) ) );
                }
            }
        }

        max_region_label = watersheds3DTwentySix( vol, labelVolume2);
        shouldEqual(8, max_region_label);
        should(labelVolume == labelVolume2);

        labelVolume2.init(0);
        max_region_label = watershedsMultiArray(vol, labelVolume2, DirectNeighborhood, WatershedOptions().unionFind());
        shouldEqual(8, max_region_label);
        should(labelVolume == labelVolume2);

        labelVolume2.init(0);
        max_region_label = watershedsMultiArray(vol, labelVolume2, IndirectNeighborhood, WatershedOptions().unionFind());
        shouldEqual(8, max_region_label);
        should(labelVolume == labelVolume2);
    }
};


struct Watershed3DTestSuite
: public test_suite
{
    Watershed3DTestSuite()
    : test_suite("Watershed3DTestSuite")
    {
        add( testCase( &Watersheds3dTest::testDistanceVolumesSix));
        add( testCase( &Watersheds3dTest::testDistanceVolumesTwentySix));
        add( testCase( &Watersheds3dTest::testWatersheds3dSix1));
        add( testCase( &Watersheds3dTest::testWatersheds3dSix2));
        add( testCase( &Watersheds3dTest::testWatersheds3dGradient1));
        add( testCase( &Watersheds3dTest::testWatersheds3dGradient2));
    }
};

int main(int argc, char ** argv)
{
    Watershed3DTestSuite test;

    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;
    return (failed != 0);
}

