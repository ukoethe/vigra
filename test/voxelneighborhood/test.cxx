/************************************************************************/
/*                                                                      */
/*       Copyright 2004 by F. Heinrich, B. Seppke, Ullrich Koethe       */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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
#include "unittest.hxx"

#include "vigra/voxelneighborhood.hxx"
#include "vigra/multi_array.hxx"
#include "set"

using namespace vigra;
struct NeighborhoodTraverserTest
{
    typedef vigra::NeighborhoodTraverser<vigra::MultiIterator<3,int>, vigra::NeighborCode3DSix> SixTraverser;
    typedef vigra::NeighborhoodTraverser<vigra::MultiIterator<3,int>, vigra::NeighborCode3DTwentySix> TwentySixTraverser;
    typedef vigra::MultiArray<3,int> IntVolume;

    IntVolume vol;
        
    SixTraverser sixTrav;
    TwentySixTraverser twentySixTrav;

    static const int w=5,h=5,d=5;
        

    NeighborhoodTraverserTest()
    : vol(IntVolume::difference_type(w,h,d)),
      sixTrav(vol.traverser_begin() + vigra::Diff3D(1,1,1)),         // set sixTrav to voxel 31=(1,1,1)
      twentySixTrav(vol.traverser_begin() + vigra::Diff3D(1,1,1))    // set twentySixTrav to voxel 31=(1,1,1)
    {
        int i=0;
        for(vigra::MultiArray<3,int>::iterator iter = vol.begin(); iter!= vol.end(); ++iter, ++i){
                *iter=i;
        }
    }

    void testInit()
    {
        should(*sixTrav == 6);
        should(!sixTrav.isDiagonal());
        should(*(sixTrav.center()) == 31);
        should(*twentySixTrav == 0);
        should(twentySixTrav.isDiagonal());
        should(*(twentySixTrav.center()) == 31);
    }

    void testSixTraverserForward()
    {
        sixTrav++;
        shouldEqual(*sixTrav, 26);
        sixTrav++;
        shouldEqual(*sixTrav, 30);
        sixTrav++;
        shouldEqual(*sixTrav, 56);
        sixTrav++;
        shouldEqual(*sixTrav, 36);
        sixTrav++;
        shouldEqual(*sixTrav, 32);
    }

    void testSixTraverserBackward()
    {
        sixTrav--;
        shouldEqual(*sixTrav, 32);
        sixTrav--;
        shouldEqual(*sixTrav, 36);
        sixTrav--;
        shouldEqual(*sixTrav, 56);
        sixTrav--;
        shouldEqual(*sixTrav, 30);
        sixTrav--;
        shouldEqual(*sixTrav, 26);
    }

    void testTwentySixTraverserForward()
    {
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 1);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 2);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 5);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 6);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 7);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 10);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 11);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 12);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 25);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 26);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 27);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 30);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 32);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 35);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 36);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 37);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 50);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 51);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 52);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 55);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 56);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 57);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 60);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 61);
        twentySixTrav++;
        shouldEqual(*twentySixTrav, 62);
    }

    void testTwentySixTraverserBackward()
    {
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 62);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 61);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 60);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 57);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 56);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 55);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 52);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 51);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 50);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 37);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 36);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 35);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 32);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 30);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 27);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 26);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 25);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 12);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 11);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 10);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 7);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 6);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 5);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 2);
        twentySixTrav--;
        shouldEqual(*twentySixTrav, 1);
    }

    void testIsDiagonal()
    {
        should(*twentySixTrav == 0);
        should(*sixTrav == 6);
        for(int i=0; i<27; i++, sixTrav++, twentySixTrav++)
        {
            switch(i){
                case  4 : 
                case 10 :
                case 12 :
                case 13 :
                case 15 :
                case 21 : should(!twentySixTrav.isDiagonal()); break;
                default : should(twentySixTrav.isDiagonal());
            }
            should(!sixTrav.isDiagonal());
        }
    }

    void testEquality()
    {
        TwentySixTraverser twentySixTrav2 = twentySixTrav;
        should(twentySixTrav == twentySixTrav2);
        twentySixTrav2++;
        should(twentySixTrav != twentySixTrav2);
        twentySixTrav2++;
        should(twentySixTrav != twentySixTrav2);
        twentySixTrav++;
        should(twentySixTrav != twentySixTrav2);
        twentySixTrav2--;
        should(twentySixTrav == twentySixTrav2);

        SixTraverser sixTrav2(vol.traverser_begin() + vigra::Diff3D(1,1,1));
        should(sixTrav == sixTrav2);
        sixTrav--;
        should(sixTrav != sixTrav2);

        twentySixTrav2 = twentySixTrav + 3;
        twentySixTrav += 3;
        should(twentySixTrav == twentySixTrav2);

        sixTrav2 = sixTrav + 3;
        sixTrav += 3;
        should(sixTrav == sixTrav2);

        twentySixTrav2 = twentySixTrav - 5;
        twentySixTrav -= 5;
        should(twentySixTrav == twentySixTrav2);

        sixTrav2 = sixTrav - 5;
        sixTrav -= 5;
        should(sixTrav == sixTrav2);
    }

    void testTurning()
    {
        SixTraverser sixTrav2 = sixTrav;
        sixTrav2.turnRound();
        should(sixTrav != sixTrav2);
        sixTrav2.turnRound();
        should(sixTrav == sixTrav2);

        TwentySixTraverser twentySixTrav2 = twentySixTrav;
        twentySixTrav2.turnRound();
        should(twentySixTrav != twentySixTrav2);
        twentySixTrav2.turnRound();
        should(twentySixTrav == twentySixTrav2);

        sixTrav++; // turn to North
        sixTrav2.turnTo(Neighborhood3DSix::North);
        should(sixTrav == sixTrav2);
        sixTrav-=2; // turn to East
        sixTrav2.turnTo(Neighborhood3DSix::East);
        should(sixTrav == sixTrav2);

        twentySixTrav++; // turn to InFrontNorth
        twentySixTrav2.turnTo(Neighborhood3DTwentySix::InFrontNorth);
        should(twentySixTrav == twentySixTrav2);
        twentySixTrav+=12; // turn to East
        twentySixTrav2.turnTo(Neighborhood3DTwentySix::East);
        should(twentySixTrav == twentySixTrav2);
    }

    void testMoving()
    {
        sixTrav.turnTo(Neighborhood3DSix::Behind);
        sixTrav.swapCenterNeighbor();
        should(*(sixTrav.center()) == 56 );
        shouldEqual(*sixTrav, 31); // looking Behind from 6 now
        should(sixTrav.direction() == Neighborhood3DSix::InFront);

        twentySixTrav.turnTo(Neighborhood3DTwentySix::BehindSouthEast);
        twentySixTrav.swapCenterNeighbor();
        should(*(twentySixTrav.center()) == 62 );
        shouldEqual(*twentySixTrav, 31); // looking Behind from 6 now
        should(twentySixTrav.direction() == Neighborhood3DTwentySix::InFrontNorthWest);
    }

    void testMiscellaneous()
    {
        SixTraverser sixTrav2 = sixTrav;
        TwentySixTraverser twentySixTrav2 = twentySixTrav;
        // test operator[]
        for(int i=0; i<6; i++, sixTrav++, twentySixTrav++){
                should(sixTrav2[i] == *sixTrav);
                should(twentySixTrav2[i] == *twentySixTrav);
        }

        twentySixTrav.turnTo(Neighborhood3DTwentySix::InFrontNorthWest);

        // test base()
        should(twentySixTrav.base() == vol.traverser_begin());
        twentySixTrav += Neighborhood3DTwentySix::BehindEast - twentySixTrav.direction();
        twentySixTrav.moveCenterToNeighbor();
        should(*(twentySixTrav.base()) == 83);
        twentySixTrav.turnRound();
        twentySixTrav.swapCenterNeighbor();
        should(*(twentySixTrav.base()) == 57);
    }
};

struct RestrictedNeighborhoodTraverserTest
{
    typedef vigra::RestrictedNeighborhoodTraverser<vigra::MultiIterator<3,int>, vigra::NeighborCode3DSix> SixTraverser;
    typedef vigra::RestrictedNeighborhoodTraverser<vigra::MultiIterator<3,int>, vigra::NeighborCode3DTwentySix> TwentySixTraverser;
    typedef vigra::MultiArray<3,int> IntVolume;

    IntVolume vol;
        
    SixTraverser sixTrav;
    TwentySixTraverser twentySixTrav;

    static const int w=3,h=3,d=3;
        

    RestrictedNeighborhoodTraverserTest()
    : vol(IntVolume::difference_type(w,h,d)),
      sixTrav(vol.traverser_begin(), TopLeftFrontBorder),
      twentySixTrav(vol.traverser_begin(), TopLeftFrontBorder)
    {
        int i=0;
        for(vigra::MultiArray<3,int>::iterator iter = vol.begin(); iter!= vol.end(); ++iter, ++i){
            *iter=i;
        }
    }

    void testInit()
    {
        should(*sixTrav == 9);
        should(!sixTrav.isDiagonal());
        should(*(sixTrav.center()) == 0);
        should(*twentySixTrav == 1);
        should(!twentySixTrav.isDiagonal());
        should(*(twentySixTrav.center()) == 0);
    }

    void testBordersSix()
    {
        std::set<int> directions;
        std::set<int>::iterator dir_Iter;
        int x,y,z,start;
        AtVolumeBorder atBorder;

////////test first plane in volume

        //-------------TopLeftFrontBorder--------------//
        x=0,y=0,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 9);
        should(atBorder == TopLeftFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DSix::InFront); // insert "0"   \
        //directions.insert(Neighborhood3DSix::North);   // insert "1"   |- border directions
        //directions.insert(Neighborhood3DSix::West);    // insert "2"   /
        directions.insert(Neighborhood3DSix::Behind); // insert "3"
        directions.insert(Neighborhood3DSix::South);  // insert "4"
        directions.insert(Neighborhood3DSix::East);   // insert "5"

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------TopFrontBorder--------------//
        x=1,y=0,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 0);
        should(atBorder == TopFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(2);
        directions.insert(3);
        directions.insert(4);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------TopRightFrontBorder--------------//
        x=2,y=0,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 1);
        should(atBorder == TopRightFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(2);
        directions.insert(3);
        directions.insert(4);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------FrontLeftBorder--------------//
        x=0,y=1,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 0);
        should(atBorder == FrontLeftBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(3);
        directions.insert(4);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------FrontBorder--------------//
        x=1,y=1,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 1);
        should(atBorder == FrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(2);
        directions.insert(3);
        directions.insert(4);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------FrontRightBorder--------------//
        x=2,y=1,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 2);
        should(atBorder == FrontRightBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(2);
        directions.insert(3);
        directions.insert(4);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------BottomLeftFrontBorder--------------//
        x=0,y=2,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 3);
        should(atBorder == BottomLeftFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(3);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------BottomFrontBorder--------------//
        x=1,y=2,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 4);
        should(atBorder == BottomFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(2);
        directions.insert(3);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------BottomRightFrontBorder--------------//
        x=2,y=2,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 5);
        should(atBorder == BottomRightFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(2);
        directions.insert(3);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);


////////test second plane in volume

        //-------------TopLeftBorder--------------//
        x=0,y=0,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 0);
        should(atBorder == TopLeftBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(Neighborhood3DSix::InFront); // insert "1"
        directions.insert(Neighborhood3DSix::Behind); // insert "3"
        directions.insert(Neighborhood3DSix::South);  // insert "4"
        directions.insert(Neighborhood3DSix::East);   // insert "5"

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------TopBorder--------------//
        x=1,y=0,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 1);
        should(atBorder == TopBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(0);
        directions.insert(2);
        directions.insert(3);
        directions.insert(4);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------TopRightBorder--------------//
        x=2,y=0,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 2);
        should(atBorder == TopRightBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(0);
        directions.insert(2);
        directions.insert(3);
        directions.insert(4);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------LeftBorder--------------//
        x=0,y=1,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 3);
        should(atBorder == LeftBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(0);
        directions.insert(1);
        directions.insert(3);
        directions.insert(4);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------NotAtBorder--------------//
        x=1,y=1,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 4);
        should(atBorder == NotAtBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(0);
        directions.insert(1);
        directions.insert(2);
        directions.insert(3);
        directions.insert(4);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------RightBorder--------------//
        x=2,y=1,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 5);
        should(atBorder == RightBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(0);
        directions.insert(1);
        directions.insert(2);
        directions.insert(3);
        directions.insert(4);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------BottomLeftBorder--------------//
        x=0,y=2,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 6);
        should(atBorder == BottomLeftBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(0);
        directions.insert(1);
        directions.insert(3);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------BottomBorder--------------//
        x=1,y=2,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 7);
        should(atBorder == BottomBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(0);
        directions.insert(1);
        directions.insert(2);
        directions.insert(3);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------BottomRightBorder--------------//
        x=2,y=2,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 8);
        should(atBorder == BottomRightBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(0);
        directions.insert(1);
        directions.insert(2);
        directions.insert(3);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);


////////test third plane in volume

        //-------------TopLeftRearBorder--------------//
        x=0,y=0,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 9);
        should(atBorder == TopLeftRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(Neighborhood3DSix::InFront);// insert "0"
        directions.insert(Neighborhood3DSix::South);  // insert "4"
        directions.insert(Neighborhood3DSix::East);   // insert "5"

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------TopRearBorder--------------//
        x=1,y=0,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 10);
        should(atBorder == TopRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(2);
        directions.insert(0);
        directions.insert(4);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------TopRightRearBorder--------------//
        x=2,y=0,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 11);
        should(atBorder == TopRightRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(2);
        directions.insert(0);
        directions.insert(4);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------RearLeftBorder--------------//
        x=0,y=1,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 12);
        should(atBorder == RearLeftBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(0);
        directions.insert(4);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------RearBorder--------------//
        x=1,y=1,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 13);
        should(atBorder == RearBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(2);
        directions.insert(0);
        directions.insert(4);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------RearRightBorder--------------//
        x=2,y=1,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 14);
        should(atBorder == RearRightBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(2);
        directions.insert(0);
        directions.insert(4);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------BottomLeftRearBorder--------------//
        x=0,y=2,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 15);
        should(atBorder == BottomLeftRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(0);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------BottomRearBorder--------------//
        x=1,y=2,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 16);
        should(atBorder == BottomRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(2);
        directions.insert(0);
        directions.insert(5);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);

        //-------------BottomRightRearBorder--------------//
        x=2,y=2,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        sixTrav = SixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*sixTrav == 17);
        should(atBorder == BottomRightRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(1);
        directions.insert(2);
        directions.insert(0);

        start = *sixTrav;
        dir_Iter=directions.begin();
        do{
            should(*dir_Iter == sixTrav.direction());
            sixTrav++;
            dir_Iter++;
        }while(*sixTrav!=start);
    }


    void testBordersTwentySix()
    {
        std::set<int> directions;
        std::set<int>::iterator dir_Iter;
        int x,y,z,start;
        AtVolumeBorder atBorder;

////////test first plane in volume

        //-------------TopLeftFrontBorder--------------//
        x=0,y=0,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 1);
        should(atBorder == TopLeftFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        //directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        //directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        //directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------TopFrontBorder--------------//
        x=1,y=0,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 0);
        should(atBorder == TopFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        //directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        //directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------TopRightFrontBorder--------------//
        x=2,y=0,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 1);
        should(atBorder == TopRightFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        //directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        //directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        //directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------FrontLeftBorder--------------//
        x=0,y=1,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 0);
        should(atBorder == FrontLeftBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        //directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        //directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------FrontBorder--------------//
        x=1,y=1,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 0);
        should(atBorder == FrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        //directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------FrontRightBorder--------------//
        x=2,y=1,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 1);
        should(atBorder == FrontRightBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        //directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        //directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------BottomLeftFrontBorder--------------//
        x=0,y=2,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 3);
        should(atBorder == BottomLeftFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        //directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        //directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        //directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------BottomFrontBorder--------------//
        x=1,y=2,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 3);
        should(atBorder == BottomFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        //directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        //directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------BottomRightFrontBorder--------------//
        x=2,y=2,z=0;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 4);
        should(atBorder == BottomRightFrontBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        //directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        //directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        //directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);


////////test middle plane in volume

        //-------------TopLeftBorder--------------//
        x=0,y=0,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 0);
        should(atBorder == TopLeftBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        //directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        //directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------TopBorder--------------//
        x=1,y=0,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 0);
        should(atBorder == TopBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        //directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------TopRightBorder--------------//
        x=2,y=0,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 1);
        should(atBorder == TopRightBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        //directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        //directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------LeftBorder--------------//
        x=0,y=1,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 0);
        should(atBorder == LeftBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        //directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------NotAtBorder--------------//
        x=1,y=1,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 0);
        should(atBorder == NotAtBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------RightBorder--------------//
        x=2,y=1,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 1);
        should(atBorder == RightBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        //directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------BottomLeftBorder--------------//
        x=0,y=2,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 3);
        should(atBorder == BottomLeftBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        //directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        //directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------BottomBorder--------------//
        x=1,y=2,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 3);
        should(atBorder == BottomBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        //directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------BottomRightBorder--------------//
        x=2,y=2,z=1;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 4);
        should(atBorder == BottomRightBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        //directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        //directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        directions.insert(Neighborhood3DTwentySix::BehindWest);
        directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

////////test back plane in volume

        //-------------TopLeftRearBorder--------------//
        x=0,y=0,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 9);
        should(atBorder == TopLeftRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        //directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        //directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        //directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------TopRearBorder--------------//
        x=1,y=0,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 9);
        should(atBorder == TopRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        //directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        //directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------TopRightRearBorder--------------//
        x=2,y=0,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 10);
        should(atBorder == TopRightRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        //directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        //directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        //directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------RearLeftBorder--------------//
        x=0,y=1,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 9);
        should(atBorder == RearLeftBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        //directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        //directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------RearBorder--------------//
        x=1,y=1,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 9);
        should(atBorder == RearBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        //directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------RearRightBorder--------------//
        x=2,y=1,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 10);
        should(atBorder == RearRightBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        //directions.insert(Neighborhood3DTwentySix::East);
        directions.insert(Neighborhood3DTwentySix::SouthWest);
        directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        //directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------BottomLeftRearBorder--------------//
        x=0,y=2,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 12);
        should(atBorder == BottomLeftRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        //directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        //directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        //directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        //directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------BottomRearBorder--------------//
        x=1,y=2,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 12);
        should(atBorder == BottomRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        //directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        //directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

        //-------------BottomRightRearBorder--------------//
        x=2,y=2,z=2;
        atBorder = isAtVolumeBorder(x,y,z,w,h,d);
        twentySixTrav = TwentySixTraverser(vol.traverser_begin() + vigra::Diff3D(x,y,z), atBorder);
        should(*twentySixTrav == 13);
        should(atBorder == BottomRightRearBorder);

        //clear direction list and fill to compare
        directions.clear();
        directions.insert(Neighborhood3DTwentySix::InFrontNorthWest);
        directions.insert(Neighborhood3DTwentySix::InFrontNorth);
        //directions.insert(Neighborhood3DTwentySix::InFrontNorthEast);
        directions.insert(Neighborhood3DTwentySix::InFrontWest);
        directions.insert(Neighborhood3DTwentySix::InFront);
        //directions.insert(Neighborhood3DTwentySix::InFrontEast);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthWest);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouth);
        //directions.insert(Neighborhood3DTwentySix::InFrontSouthEast);

        directions.insert(Neighborhood3DTwentySix::NorthWest);
        directions.insert(Neighborhood3DTwentySix::North);
        //directions.insert(Neighborhood3DTwentySix::NorthEast);
        directions.insert(Neighborhood3DTwentySix::West);
        //directions.insert(Neighborhood3DTwentySix::East);
        //directions.insert(Neighborhood3DTwentySix::SouthWest);
        //directions.insert(Neighborhood3DTwentySix::South);
        //directions.insert(Neighborhood3DTwentySix::SouthEast);

        //directions.insert(Neighborhood3DTwentySix::BehindNorthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindNorth);
        //directions.insert(Neighborhood3DTwentySix::BehindNorthEast);
        //directions.insert(Neighborhood3DTwentySix::BehindWest);
        //directions.insert(Neighborhood3DTwentySix::Behind);
        //directions.insert(Neighborhood3DTwentySix::BehindEast);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthWest);
        //directions.insert(Neighborhood3DTwentySix::BehindSouth);
        //directions.insert(Neighborhood3DTwentySix::BehindSouthEast);


        start = *twentySixTrav;
        dir_Iter=directions.begin();
        do{
            //std::cerr << *dir_Iter << " =?= " << twentySixTrav.direction() << std::endl;
            should(*dir_Iter == twentySixTrav.direction());
            twentySixTrav++;
            dir_Iter++;
        }while(*twentySixTrav!=start);

    }
};

struct SimpleAnalysisTestSuite
: public vigra::test_suite
{
    SimpleAnalysisTestSuite()
    : vigra::test_suite("SimpleAnalysisTestSuite")
    {
        add( testCase( &NeighborhoodTraverserTest::testInit));
        add( testCase( &NeighborhoodTraverserTest::testSixTraverserForward));
        add( testCase( &NeighborhoodTraverserTest::testSixTraverserBackward));
        add( testCase( &NeighborhoodTraverserTest::testTwentySixTraverserForward));
        add( testCase( &NeighborhoodTraverserTest::testTwentySixTraverserBackward));
        add( testCase( &NeighborhoodTraverserTest::testIsDiagonal));
        add( testCase( &NeighborhoodTraverserTest::testEquality));
        add( testCase( &NeighborhoodTraverserTest::testTurning));
        add( testCase( &NeighborhoodTraverserTest::testMoving));
        add( testCase( &NeighborhoodTraverserTest::testMiscellaneous));

        add( testCase( &RestrictedNeighborhoodTraverserTest::testInit));
        add( testCase( &RestrictedNeighborhoodTraverserTest::testBordersSix));
        add( testCase( &RestrictedNeighborhoodTraverserTest::testBordersTwentySix));
    }
};

int main()
{
    SimpleAnalysisTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;
    return (failed != 0);
}

