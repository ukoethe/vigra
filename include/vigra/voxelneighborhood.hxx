/************************************************************************/
/*                                                                      */
/*     Copyright 1998-2006 by F. Heinrich, B. Seppke, Ullrich Koethe    */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_VOXELNEIGHBORHOOD_HXX
#define VIGRA_VOXELNEIGHBORHOOD_HXX

#include "vigra/tinyvector.hxx"

namespace vigra {
    
    //Define Diff3D for vigraext
    typedef vigra::TinyVector<int, 3> Diff3D;
   
/********************************************************/
/*                                                      */
/*                      AtVolumeBorder                  */
/*                                                      */
/********************************************************/

/** \brief Encode whether a voxel is near the volume border.

        //    This enum is used with \ref isAtVolumeBorder() and 
        //    \ref vigra::RestrictedNeighborhoodTraverser.

        //<b>\#include</b> "<a href="pixelneighborhood_8hxx-source.html">vigra/pixelneighborhood.hxx</a>"<br>
        //Namespace: vigra
*/
enum AtVolumeBorder
{
        NotAtBorder       = 0,     ///< &nbsp;
        RightBorder       = 1,     ///< &nbsp;
        LeftBorder        = 2,     ///< &nbsp;
        TopBorder         = 4,     ///< &nbsp;
        BottomBorder      = 8,     ///< &nbsp;
        FrontBorder       = 16,    ///< &nbsp;
        RearBorder        = 32,
        TopRightBorder    = TopBorder    | RightBorder,   //5
        TopLeftBorder     = TopBorder    | LeftBorder,    //6
        TopFrontBorder    = TopBorder    | FrontBorder,   //20
        TopRearBorder     = TopBorder    | RearBorder,    //36
        BottomLeftBorder  = BottomBorder | LeftBorder,    //10
        BottomRightBorder = BottomBorder | RightBorder,   //9
        BottomFrontBorder = BottomBorder | FrontBorder,   //24
        BottomRearBorder  = BottomBorder | RearBorder,    //40
        FrontLeftBorder   = FrontBorder  | LeftBorder,    //18
        FrontRightBorder  = FrontBorder  | RightBorder,   //17
        RearLeftBorder    = RearBorder   | LeftBorder,    //34
        RearRightBorder   = RearBorder   | RightBorder,   //33
        
        TopRightFrontBorder    = TopBorder    | RightBorder | FrontBorder,    //21
        TopLeftFrontBorder     = TopBorder    | LeftBorder  | FrontBorder,    //22
        BottomLeftFrontBorder  = BottomBorder | LeftBorder  | FrontBorder,    //26
        BottomRightFrontBorder = BottomBorder | RightBorder | FrontBorder,    //25
        TopRightRearBorder     = TopBorder    | RightBorder | RearBorder,     //37
        TopLeftRearBorder      = TopBorder    | LeftBorder  | RearBorder,     //38
        BottomLeftRearBorder   = BottomBorder | LeftBorder  | RearBorder,     //42
        BottomRightRearBorder  = BottomBorder | RightBorder | RearBorder      //41
};

/** \brief Find out whether a voxel is at the volume border.

    This function checks if \a x == 0 or \a x == \a width - 1 and 
    \a y == 0 or \a y == \a height - 1 and so on and returns the appropriate value
    of \ref vigraext::AtVolumeBorder, or zero when the voxel is not at te volume border.
    The behavior of the function is undefined if (x,y,z) is not inside the volume.
*/
inline AtVolumeBorder isAtVolumeBorder(int x, int y, int z, int width, int height, int depth)
{
    return static_cast<AtVolumeBorder>((x == 0 
                                         ? LeftBorder
                                         : x == width-1
                                             ? RightBorder
                                             : NotAtBorder) |
                                       (y == 0 
                                         ? TopBorder
                                         : y == height-1
                                             ? BottomBorder
                                             : NotAtBorder) |
                                       (z == 0 
                                         ? FrontBorder
                                         : z == depth-1
                                             ? RearBorder
                                             : NotAtBorder));
};
/** \brief Find out whether a voxel is at a scan-order relevant volume border.
    This function checks if \a x == 0 or \a y == 0 or \a z == \a 0 and returns the 
        appropriate value of \ref vigraext::AtVolumeBorder, or zero when the voxel is
        not at te volume border.
    The behavior of the function is undefined if (x,y,z) is not inside the volume.
*/
inline AtVolumeBorder isAtVolumeBorderCausal(int x, int y, int z, int width, int height, int depth)
{
    return static_cast<AtVolumeBorder>((x == 0 
                                         ? LeftBorder
                                         : NotAtBorder) |
                                       (y == 0 
                                         ? TopBorder
                                         : NotAtBorder) |
                                       (z == 0 
                                         ? FrontBorder
                                         : NotAtBorder));
};
/** TODO: Write new comment \brief Find out whether a voxel is at a scan-order relevant volume border.
    This function checks if \a x == 0 or \a y == 0 or \a z == \a 0 and returns the 
        appropriate value of \ref vigraext::AtVolumeBorder, or zero when the voxel is
        not at te volume border.
    The behavior of the function is undefined if (x,y,z) is not inside the volume.
*/
inline AtVolumeBorder isAtVolumeBorderAntiCausal(int x, int y, int z, int width, int height, int depth)
{
    return static_cast<AtVolumeBorder>((x == width-1 
                                         ? RightBorder
                                         : NotAtBorder) |
                                       (y == height-1 
                                         ? BottomBorder
                                         : NotAtBorder) |
                                       (z == depth-1 
                                         ? RearBorder
                                         : NotAtBorder));
};

/********************************************************/
/*                                                      */
/*                   Neighborhood3DSix                  */
/*                                                      */
/********************************************************/

/** 3D 6-Neighborhood. */
namespace Neighborhood3DSix
{

/** \brief Encapsulation of direction management of neighbors for a 3D 6-neighborhood.
*/
class NeighborCode3D
{
    public:
    /** provides enumeration of all directions.
        DirectionCount may be used for portable loop termination conditions.
    */
    enum Direction {
        Error = -1,    
        InFront= 0,         
        North,          
        West ,      
        Behind,
        South,
        East,         
        DirectionCount,
        CausalFirst = InFront,
        CausalLast  = West,
        AntiCausalFirst = Behind,
        AntiCausalLast  = East,

        OppositeDirPrefix = 1,
        OppositeOffset = 3
    };

    static unsigned int directionBit(Direction d) 
    {
        static unsigned int b[] = { 1 << (InFront + 1),
                                    1 << (North + 1), 
                                    1 << (West + 1), 
                                    1 << (Behind + 1),
                                    1 << (South + 1), 
                                    1 << (East + 1)
                                  };
        return b[d];
    };


    /** The number of valid neighbors if the current center is at the volume border.
    */
    static unsigned int nearBorderDirectionCount(AtVolumeBorder b)
    {
        static unsigned int c[] = { 6, 5, 5, 0, 5, 4, 4, 0, 5, 4, 
                                    4, 0, 0, 0, 0, 0, 5, 4, 4, 0,
                                    4, 3, 3, 0, 4, 3, 3, 0, 0, 0,
                                    0, 0, 5, 4, 4, 0, 4, 3, 3, 0,
                                    4, 3, 3};
        return c[b];
    }
    
    /** The valid direction codes when the center is at the volume border.
        \a index must be in the range <tt>0...nearBorderDirectionCount(b)-1</tt>.
    */
    static Direction nearBorderDirections(AtVolumeBorder b, int index)
    {
        static Direction c[43][6] = {
                { InFront, North, West, Behind, South, East},   // 0 - NotAtBorder
                { InFront, North, West, Behind, South, Error},  // 1 - AtRightBorder
                { InFront, North, Behind, South, East, Error},  // 2 - AtLeftBorder
                { Error, Error, Error, Error, Error, Error},    
                { InFront, West, Behind, South, East, Error},   // 4 - AtTopBorder    
                { InFront, West, Behind, South, Error, Error},  // 5 - AtTopRightBorder
                { InFront, Behind, South, East, Error, Error},  // 6 - AtTopLeftBorder
                { Error, Error, Error, Error, Error, Error},    
                { InFront, North, West, Behind, East, Error},   // 8 - AtBottomBorder
                { InFront, North, West, Behind, Error, Error},  // 9 - AtBottomRightBorder
                { InFront, North, Behind, East, Error, Error},  //10- AtBottomLeftBorder
                { Error, Error, Error, Error, Error, Error},    
                { Error, Error, Error, Error, Error, Error},    
                { Error, Error, Error, Error, Error, Error},    
                { Error, Error, Error, Error, Error, Error},    
                { Error, Error, Error, Error, Error, Error},    
                { North, West, Behind, South, East, Error},     //16 - AtFrontBorder
                { North, West, Behind, South, Error, Error},    //17 - AtFrontRightBorder
                { North, Behind, South, East, Error, Error},    //18 - AtFrontLeftBorder
                { Error, Error, Error, Error, Error, Error},    
                { West, Behind, South, East, Error,  Error},    //20 - AtTopFrontBorder
                { West, Behind, South, Error, Error, Error},    //21 - AtTopRightFrontBorder
                { Behind, South, East, Error, Error,  Error},   //22 - AtTopLeftFrontBorder
                { Error, Error, Error, Error, Error, Error},    
                { North, West, Behind, East, Error, Error},     //24 - AtBottomFrontBorder
                { North, West, Behind, Error, Error, Error},    //25 - AtBottomRightFrontBorder
                { North, Behind, East, Error, Error, Error},    //26 - AtBottomLeftFrontBorder
                { Error, Error, Error, Error, Error, Error},    
                { Error, Error, Error, Error, Error, Error},    
                { Error, Error, Error, Error, Error, Error},    
                { Error, Error, Error, Error, Error, Error},    
                { Error, Error, Error, Error, Error, Error},    
                { InFront, North, West, South, East,Error},     //32 - AtRearBorder
                { InFront, North, West, South, Error, Error},   //33 - AtRearRightBorder
                { InFront, North, South, East, Error, Error},   //34 - AtRearLeftBorder
                { Error, Error, Error, Error, Error, Error},    
                { InFront, West, South, East, Error, Error},    //36 - AtTopRearBorder    
                { InFront, West, South, Error, Error, Error},   //37 - AtTopRightRearBorder    
                { InFront, South, East, Error, Error, Error},   //38 - AtTopLeftRearBorder
                { Error, Error, Error, Error, Error, Error},    
                { InFront, North, West, East, Error, Error},    //40 - AtBottomRearBorder
                { InFront, North, West, Error, Error, Error},   //41 - AtBottomRightRearBorder
                { InFront, North, East, Error, Error, Error}    //42 - AtBottomLeftRearBorder
               };
        return c[b][index];
    }
        
    /** The valid direction three codes in anti causal direction (means: look back in scanline
        direction)when the center is at the volume border.
        Should be used with isAtVolumeBorderCausal to determine the Directions, as this
        avoids using of the nonesense border ids (e.g. 0,1,8,9...) of this table.
        \a index must be in the range <tt>0...nearBorderDirectionCount(b)-1</tt>.
    */
    static Direction nearBorderDirectionsCausal(AtVolumeBorder b, int index)
    {
        static Direction c[43][3] = {
            { InFront, North, West},                    // 0 - NotAtBorder -----> should never be used
            { InFront, North, West},                    // 1 - AtRightBorder -----> should never be used
            { InFront, North, Error},                   // 2 - AtLeftBorder
            { Error, Error, Error},    
            { InFront, West, Error},                    // 4 - AtTopBorder    
            { InFront, West, Error},                    // 5 - AtTopRightBorder
            { InFront, Error,Error},                    // 6 - AtTopLeftBorder
            { Error, Error, Error},    
            { InFront, North, West},                    // 8 - AtBottomBorder -----> should never be used
            { InFront, North, West},                    // 9 - AtBottomRightBorder -----> should never be used
            { InFront, North, Error},                   //10- AtBottomLeftBorder
            { Error, Error, Error},    
            { Error, Error, Error},    
            { Error, Error, Error},    
            { Error, Error, Error},    
            { Error, Error, Error},    
            { North, West, Error},                      //16 - AtFrontBorder
            { North, West, Error},                      //17 - AtFrontRightBorder
            { North, Error, Error},                     //18 - AtFrontLeftBorder
            { Error, Error, Error},    
            { West, Error, Error},                      //20 - AtTopFrontBorder
            { West, Error, Error},                      //21 - AtTopRightFrontBorder
            { Error, Error,  Error},                    //22 - AtTopLeftFrontBorder
            { Error, Error, Error},    
            { North, West, Error},                      //24 - AtBottomFrontBorder
            { North, West, Error},                      //25 - AtBottomRightFrontBorder
            { North, Error, Error},                     //26 - AtBottomLeftFrontBorder
            { Error, Error, Error},    
            { Error, Error, Error},    
            { Error, Error, Error},    
            { Error, Error, Error},    
            { Error, Error, Error},                      
            { InFront, North, West},                    //32 - AtRearBorder -----> should never be used
            { InFront, North, West},                    //33 - AtRearRightBorder -----> should never be used
            { InFront, North, Error},                   //34 - AtRearLeftBorder
            { Error, Error, Error},    
            { InFront, West, Error},                    //36 - AtTopRearBorder    
            { InFront, West, Error},                    //37 - AtTopRightRearBorder    
            { InFront, Error, Error},                   //38 - AtTopLeftRearBorder
            { Error, Error, Error},    
            { InFront, North, West},                    //40 - AtBottomRearBorder -----> should never be used
            { InFront, North, West},                    //41 - AtBottomRightRearBorder -----> should never be used
            { InFront, North, Error}                    //42 - AtBottomLeftRearBorder
        };
        return c[b][index];
    }
    
    /** transform direction code into corresponding Diff3D offset.
        (note: there is no bounds checking on the code you pass.)
    */
    static Diff3D const & diff(Direction code)
    {
        static Diff3D d[] = {
                    Diff3D(  0,  0, -1),  //InFront
                    Diff3D(  0, -1,  0),  //North
                    Diff3D( -1,  0,  0),  //West
                    Diff3D(  0,  0,  1),  //Behind
                    Diff3D(  0,  1,  0),  //South
                    Diff3D(  1,  0,  0)   //East
                };
        return d[code];
    }

    /** Equivalent to <tt>diff(static_cast<Direction>(code))</tt>.
        (note: there is no bounds checking on the code you pass.)
    */
    static Diff3D const & diff(int code) { return diff(static_cast<Direction>(code)); }

    /**  Equivalent to <tt>diff(code)[dim]</tt> */
    static int diff(Direction code, int dim) { return diff(code)[dim]; }
   
    /** Get the relative offset from one neighbor to the other.
        For example, <tt>relativeDiff(East, West) == multi_differencetype(-2,0,0)</tt>.
        (note: there is no bounds checking on the code you pass.)
    */
    static Diff3D const & relativeDiff(Direction fromCode, Direction toCode)
    {
      static Diff3D d[6][6] = 
          {    
                //     InFront      -      North         -           West     -         Behind     -      South        -        East       
                { Diff3D( 0, 0, 0), Diff3D(0, -1, 1), Diff3D(-1, 0, 1), Diff3D( 0, 0, 2), Diff3D( 0, 1, 1),  Diff3D( 1, 0, 1)}, //InFront
                { Diff3D( 0, 1,-1), Diff3D( 0, 0, 0), Diff3D(-1, 1, 0), Diff3D( 0, 1, 1), Diff3D( 0, 2, 0),  Diff3D( 1, 1, 0)}, //North
                { Diff3D( 1, 0,-1), Diff3D( 1,-1, 0), Diff3D( 0, 0, 0), Diff3D( 1, 0, 1), Diff3D( 1, 1, 0),  Diff3D( 2, 0, 0)}, //West
                { Diff3D( 0, 0,-2), Diff3D( 0,-1,-1), Diff3D(-1, 0,-1), Diff3D( 0, 0, 0), Diff3D( 0, 1,-1),  Diff3D( 1, 0,-1)}, //Behind
                { Diff3D( 0,-1,-1), Diff3D( 0,-2, 0), Diff3D(-1,-1, 0), Diff3D( 0,-1, 1), Diff3D( 0, 0, 0),  Diff3D( 1,-1, 0)}, //South
                { Diff3D(-1, 0,-1), Diff3D(-1,-1, 0), Diff3D(-2, 0, 0), Diff3D(-1, 0, 1), Diff3D(-1, 1, 0), Diff3D( 0, 0, 0) }  //East
          };

        return d[fromCode][toCode];
    }

    /** Equivalent to relativeDiff(static_cast<Direction>(fromCode), static_cast<Direction>(toCode)).
        (note: there is no bounds checking on the code you pass.)
    */
    static Diff3D const & relativeDiff(int fromCode, int toCode)
    {
        return relativeDiff(static_cast<Direction>(fromCode), static_cast<Direction>(toCode));
    }

    /**  X-component of diff() */
    static int dX(Direction code) { return diff(code)[0]; }
    /**  Y-component of diff() */
    static int dY(Direction code) { return diff(code)[1]; }
    /**  Z-component of diff() */
    static int dZ(Direction code) { return diff(code)[2]; }
    
    /**  X-component of diff() */
    static int dX(int code) { return diff(code)[0]; }
    /**  Y-component of diff() */
    static int dY(int code) { return diff(code)[1]; }
    /**  Z-component of diff() */
    static int dZ(int code) { return diff(code)[2]; }
    

    /** transform Diff3D offset into corresponding direction code.
        The code <tt>Direction::Error</tt> will be returned if <tt>diff</tt>
        is not in the 3DSix-Neighborhood.
    */
    static Direction code(Diff3D const & diff)
    {
        switch(diff[0]) {
            case  0:
            {
                switch(diff[1]) {
                    case 0:
                        switch(diff[2]) {
                            case 1:
                                return Behind;
                            case -1:
                                return InFront;
                            default:
                                return Error;
                        }
                                
                    case 1:
                        return (diff[2] == 0) ? South : Error;
                    case -1:
                        return (diff[2] == 0) ? North : Error;
                    default:
                        return Error;
                }
        }
        case -1:
            return ((diff[1] == 0) && (diff[2] == 0)) ? West : Error;
        case  1:
            return ((diff[1] == 0) && (diff[2] == 0)) ? East : Error;
        }
        return Error;
    }
  
    /** Check whether a code refers to a diagonal direction.
        Useful if you want to abstract the differences between 6- and 26-neighborhood.
        Always <tt>false</tt> for 6-neighborhood.
    */
    static bool isDiagonal(Direction) { return false; }

    static Diff3D const & right()           { return diff(East); }          /**<  Offset to the right neighbor */
    static Diff3D const & top()             { return diff(North); }         /**<  Offset to the top neighbor */
    static Diff3D const & left()            { return diff(West); }          /**<  Offset to the left neighbor */
    static Diff3D const & bottom()          { return diff(South); }         /**<  Offset to the bottom neighbor */
    static Diff3D const & rear()            { return diff(Behind); }        /**<  Offset to the rear neighbor */
    static Diff3D const & front()           { return diff(InFront); }       /**<  Offset to the neighbor in front */

    static Diff3D const & east()            { return diff(East); }          /**<  Offset to the east neighbor */
    static Diff3D const & north()           { return diff(North); }         /**<  Offset to the north neighbor */
    static Diff3D const & west()            { return diff(West); }          /**<  Offset to the west neighbor */
    static Diff3D const & south()           { return diff(South); }         /**<  Offset to the south neighbor */
    static Diff3D const & behind()          { return diff(Behind); }        /**<  Offset to the rear neighbor */
    static Diff3D const & infront()         { return diff(InFront); }       /**<  Offset to the neighbor in front */

}; // class Neighborhood3DSix


/** Export NeighborCode3D::Direction into the scope of namespace Neighborhood3DSix.
*/
typedef NeighborCode3D::Direction Direction;

static const Direction East           = NeighborCode3D::East;               /**<  Export NeighborCode3D::East to namespace Neighborhood3DSix */
static const Direction North          = NeighborCode3D::North;              /**<  Export NeighborCode3D::North to namespace Neighborhood3DSix */
static const Direction West           = NeighborCode3D::West;               /**<  Export NeighborCode3D::West to namespace Neighborhood3DSix */
static const Direction South          = NeighborCode3D::South;              /**<  Export NeighborCode3D::South to namespace Neighborhood3DSix */
static const Direction Behind         = NeighborCode3D::Behind;             /**<  Export NeighborCode3D::Behind to namespace Neighborhood3DSix */
static const Direction InFront        = NeighborCode3D::InFront;            /**<  Export NeighborCode3D::InFront to namespace Neighborhood3DSix */
static const Direction DirectionCount = NeighborCode3D::DirectionCount;     /**<  Export NeighborCode3D::DirectionCount to namespace Neighborhood3DSix */


}//namespace Neighborhood3DSix
    
/** Export \ref vigraext::Neighborhood3DSix::NeighborCode3D into the scope of namespace vigraext.
*/
typedef Neighborhood3DSix::NeighborCode3D NeighborCode3DSix;

/********************************************************/
/*                                                      */
/*                   Neighborhood3DTwentySix            */
/*                                                      */
/********************************************************/
/** 3D 26-Neighborhood. */
namespace Neighborhood3DTwentySix
{

/** \brief Encapsulation of direction management of neighbors for a 3D 26-neighborhood.
*/
class NeighborCode3D
{
    public:
   /** provides enumeration of all directions.
       DirectionCount may be used for portable loop termination conditions.
     */
    enum Direction {
        Error = -1,
            InFrontNorthWest = 0,
            InFrontNorth,
            InFrontNorthEast,
            InFrontWest,
        InFront,
            InFrontEast,
            InFrontSouthWest,
            InFrontSouth,
            InFrontSouthEast,
        
            NorthWest,
            North,
            NorthEast,
        West,
            East,
            SouthWest,
            South,
            SouthEast,

            BehindNorthWest,
            BehindNorth,
            BehindNorthEast,
            BehindWest,
        Behind,
            BehindEast,
            BehindSouthWest,
            BehindSouth,
            BehindSouthEast,

        DirectionCount,
            CausalFirst = InFrontNorthWest,
            CausalLast  = West,
            AntiCausalFirst = BehindSouthEast,
            AntiCausalLast  = East,

            OppositeDirPrefix = -1,
            OppositeOffset = 25
    };

    static unsigned int directionBit(Direction d) 
    {
        static unsigned int b[] = { 
                1 <<    (InFrontNorthWest+1),
                1 <<    (InFrontNorth+1),
                1 <<  (InFrontNorthEast+1),
                1 <<  (InFrontWest+1),
                1 <<  (InFront+1),
                1 <<  (InFrontEast+1),
                1 <<  (InFrontSouthWest+1),
                1 <<  (InFrontSouth+1),
                1 <<  (InFrontSouthEast+1),

                1 <<  (NorthWest+1),
                1 <<  (North+1),
                1 <<  (NorthEast+1),
                1 <<  (West+1),
                1 <<  (East+1),
                1 <<  (SouthWest+1),
                1 <<  (South+1),
                1 <<  (SouthEast+1),

                1 <<  (BehindNorthWest+1),
                1 <<  (BehindNorth+1),
                1 <<  (BehindNorthEast+1),
                1 <<  (BehindWest+1),
                1 <<  (Behind+1),
                1 <<  (BehindEast+1),
                1 <<  (BehindSouthWest+1),
                1 <<  (BehindSouth+1),
                1 <<  (BehindSouthEast+1)
            };
        return b[d];
    };


    /** The number of valid neighbors if the current center is at the volume border.
    */
    static unsigned int nearBorderDirectionCount(AtVolumeBorder b)
    {
        static unsigned int c[] = { 26, 17, 17,  0, 17, 11, 11,  0, 17, 11, 
                                    11,  0,  0,  0,  0,  0, 17, 11, 11,  0,
                                    11,  7,  7,  0, 11,  7,  7,  0,  0,  0,
                                     0,  0, 17, 11, 11,  0, 11,  7,  7,  0,
                                    11,  7,  7};
        return c[b];
    }
    
    /** The valid direction codes when the center is at the volume border.
        \a index must be in the range <tt>0...nearBorderDirectionCount(b)-1</tt>.
    */
    static Direction nearBorderDirections(AtVolumeBorder b, int index)
    {
        static Direction c[43][26] = {
                   //0 - NotAtBorder
                   {    InFrontNorthWest,   InFrontNorth,   InFrontNorthEast,
                        InFrontWest,        InFront,        InFrontEast,
                        InFrontSouthWest,   InFrontSouth,   InFrontSouthEast,
      
                        NorthWest, North, NorthEast,
                        West,             East,
                        SouthWest, South, SouthEast, 

                        BehindNorthWest, BehindNorth, BehindNorthEast,
                        BehindWest,      Behind,      BehindEast,
                        BehindSouthWest, BehindSouth, BehindSouthEast},    
                    //1 - AtRightBorder
                    {   InFrontNorthWest, InFrontNorth, /*InFrontNorthEast,*/
                        InFrontWest,      InFront,      /*InFrontEast,*/
                        InFrontSouthWest, InFrontSouth, /*InFrontSouthEast,*/
      
                        NorthWest, North, /*NorthEast,*/
                        West,             /*East,*/
                        SouthWest, South, /*SouthEast,*/ 

                        BehindNorthWest, BehindNorth, /*BehindNorthEast,*/
                        BehindWest,      Behind,      /*BehindEast,*/
                        BehindSouthWest, BehindSouth, /*BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //2 - AtLeftBorder
                    {   /*InFrontNorthWest,*/   InFrontNorth, InFrontNorthEast,
                        /*InFrontWest,*/       InFront,      InFrontEast,
                        /*InFrontSouthWest,*/   InFrontSouth, InFrontSouthEast,
                  
                        /*NorthWest,*/ North, NorthEast,
                        /*West,*/            East,
                        /*SouthWest,*/ South, SouthEast, 

                        /*BehindNorthWest,*/ BehindNorth, BehindNorthEast,
                        /*BehindWest,*/     Behind,      BehindEast,
                        /*BehindSouthWest,*/ BehindSouth, BehindSouthEast, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},
                    //3 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error,        Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //4 - AtTopBorder    
                    {   /*InFrontNorthWest, InFrontNorth, InFrontNorthEast,*/
                        InFrontWest,        InFront,      InFrontEast,
                        InFrontSouthWest, InFrontSouth,   InFrontSouthEast,
                  
                        /*NorthWest, North, NorthEast,*/
                        West,             East,
                        SouthWest, South, SouthEast, 

                        /*BehindNorthWest, BehindNorth, BehindNorthEast,*/
                        BehindWest,                 Behind,                 BehindEast,
                        BehindSouthWest, BehindSouth, BehindSouthEast, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //5 - AtTopRightBorder
                    {   /*InFrontNorthWest, InFrontNorth,   InFrontNorthEast,*/
                        InFrontWest,        InFront,        /*InFrontEast,*/
                        InFrontSouthWest, InFrontSouth,     /*InFrontSouthEast,*/
      
                        /*NorthWest, North, NorthEast,*/
                        West,             /*East,*/
                        SouthWest, South, /*SouthEast,*/ 

                        /*BehindNorthWest, BehindNorth, BehindNorthEast,*/
                        BehindWest,        Behind,      /*BehindEast,*/
                        BehindSouthWest, BehindSouth,   /*BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},    
                    //6 - AtTopLeftBorder
                    {   /*InFrontNorthWest,     InFrontNorth,   InFrontNorthEast,*/
                        /*InFrontWest,*/        InFront,        InFrontEast,
                        /*InFrontSouthWest,*/   InFrontSouth,   InFrontSouthEast,
       
                        /*NorthWest,    North,  NorthEast,*/
                        /*West,*/               East,
                        /*SouthWest,*/  South,  SouthEast, 

                        /*BehindNorthWest,      BehindNorth, BehindNorthEast,*/
                        /*BehindWest, */        Behind,      BehindEast,
                        /*BehindSouthWest,*/    BehindSouth, BehindSouthEast, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},
                    //7 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //8 - AtBottomBorder
                    {   InFrontNorthWest, InFrontNorth,    InFrontNorthEast,
                        InFrontWest,      InFront,         InFrontEast,
                        /*InFrontSouthWest, InFrontSouth, InFrontSouthEast,*/
                  
                        NorthWest,      North,  NorthEast,
                        West,                   East,
                        /*SouthWest,    South,  SouthEast,*/

                        BehindNorthWest,    BehindNorth, BehindNorthEast,
                        BehindWest,         Behind,      BehindEast,
                        /*BehindSouthWest,  BehindSouth, BehindSouthEast*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},                    
                    //9 - AtBottomRightBorder
                    {   InFrontNorthWest, InFrontNorth,    /*InFrontNorthEast,*/
                        InFrontWest,      InFront,         /*InFrontEast,*/
                        /*InFrontSouthWest, InFrontSouth,  InFrontSouthEast,*/
      
                        NorthWest, North,   /*NorthEast,*/
                        West,               /*East,*/
                        /*SouthWest, South, SouthEast,*/

                        BehindNorthWest, BehindNorth,   /*BehindNorthEast,*/
                        BehindWest,      Behind,        /*BehindEast,*/
                        /*BehindSouthWest, BehindSouth, BehindSouthEast*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},        
                    //10 - AtBottomLeftBorder
                    {   /*InFrontNorthWest,*/   InFrontNorth,   InFrontNorthEast,
                        /*InFrontWest,*/        InFront,        InFrontEast,
                        /*InFrontSouthWest,     InFrontSouth,   InFrontSouthEast,*/
      
                        /*NorthWest,*/  North,  NorthEast,
                        /*West,*/               East,
                        /*SouthWest,    South,  SouthEast,*/

                        /*BehindNorthWest,*/ BehindNorth,   BehindNorthEast,
                        /*BehindWest,*/      Behind,        BehindEast,
                        /*BehindSouthWest, BehindSouth,     BehindSouthEast*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},
                    //11 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //12 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //13 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //14 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //15 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //16 - AtFrontBorder
                    {   /*InFrontNorthWest, InFrontNorth,    InFrontNorthEast,
                        InFrontWest,        InFront,         InFrontEast,
                        InFrontSouthWest, InFrontSouth, InFrontSouthEast,*/
      
                        NorthWest, North, NorthEast,
                        West,             East,
                        SouthWest, South, SouthEast, 

                        BehindNorthWest, BehindNorth, BehindNorthEast,
                        BehindWest,      Behind,      BehindEast,
                        BehindSouthWest, BehindSouth, BehindSouthEast, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //17 - AtFrontRightBorder
                    {   /*InFrontNorthWest, InFrontNorth,    InFrontNorthEast,
                        InFrontWest,              InFront,                 InFrontEast,
                        InFrontSouthWest, InFrontSouth, InFrontSouthEast,*/
      
                        NorthWest, North, /*NorthEast,*/
                        West,             /*East,*/
                        SouthWest, South, /*SouthEast,*/ 

                        BehindNorthWest, BehindNorth, /*BehindNorthEast,*/
                        BehindWest,      Behind,      /*BehindEast,*/
                        BehindSouthWest, BehindSouth, /*BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},    
                    //18 - AtFrontLeftBorder
                    {   /*InFrontNorthWest, InFrontNorth,   InFrontNorthEast,
                        InFrontWest,        InFront,        InFrontEast,
                        InFrontSouthWest,   InFrontSouth,   InFrontSouthEast,*/
      
                        /*NorthWest,*/ North, NorthEast,
                        /*West,*/             East,
                        /*SouthWest,*/ South, SouthEast, 

                        /*BehindNorthWest,*/ BehindNorth, BehindNorthEast,
                        /*BehindWest,*/      Behind,      BehindEast,
                        /*BehindSouthWest,*/ BehindSouth, BehindSouthEast, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},    
                    //19 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},
                    //20 - AtTopFrontBorder
                    {   /*InFrontNorthWest, InFrontNorth,   InFrontNorthEast,
                        InFrontWest,        InFront,        InFrontEast,
                        InFrontSouthWest,   InFrontSouth,   InFrontSouthEast,*/
      
                        /*NorthWest,    North,  NorthEast,*/
                        West,                   East,
                        SouthWest,      South,  SouthEast, 

                        /*BehindNorthWest,  BehindNorth,    BehindNorthEast,*/
                        BehindWest,         Behind,         BehindEast,
                        BehindSouthWest,    BehindSouth,    BehindSouthEast, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},    
                    //21 - AtTopRightFrontBorder
                    {   /*InFrontNorthWest, InFrontNorth,  InFrontNorthEast,
                        InFrontWest,        InFront,       InFrontEast,
                        InFrontSouthWest,   InFrontSouth,  InFrontSouthEast,*/
      
                        /*NorthWest, North, NorthEast,*/
                        West,               /*East,*/
                        SouthWest,   South, /*SouthEast,*/ 

                        /*BehindNorthWest, BehindNorth, BehindNorthEast,*/
                        BehindWest,        Behind,      /*BehindEast,*/
                        BehindSouthWest, BehindSouth,   /*BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error},    
                    //22 - AtTopLeftFrontBorder
                    {   /*InFrontNorthWest, InFrontNorth,   InFrontNorthEast,
                        InFrontWest,        InFront,        InFrontEast,
                        InFrontSouthWest,   InFrontSouth,   InFrontSouthEast,*/
      
                        /*NorthWest,    North, NorthEast,*/
                        /*West,*/              East,
                        /*SouthWest,*/  South, SouthEast, 

                        /*BehindNorthWest,      BehindNorth, BehindNorthEast,*/
                        /*BehindWest,*/         Behind,      BehindEast,
                        /*BehindSouthWest,*/    BehindSouth, BehindSouthEast, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error},    
                    //23 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},
                    //24 - AtBottomFrontBorder
                    {   /*InFrontNorthWest, InFrontNorth, InFrontNorthEast,
                        InFrontWest,        InFront,      InFrontEast,
                        InFrontSouthWest,   InFrontSouth, InFrontSouthEast,*/
      
                        NorthWest,      North, NorthEast,
                        West,                  East,
                        /*SouthWest,    South, SouthEast,*/

                        BehindNorthWest,    BehindNorth, BehindNorthEast,
                        BehindWest,         Behind,      BehindEast,
                        /*BehindSouthWest,  BehindSouth, BehindSouthEast*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},    
                    //25 - AtBottomRightFrontBorder
                    {   /*InFrontNorthWest, InFrontNorth,    InFrontNorthEast,
                        InFrontWest,        InFront,         InFrontEast,
                        InFrontSouthWest, InFrontSouth, InFrontSouthEast,*/
      
                        NorthWest,      North,  /*NorthEast,*/
                        West,                   /* East,*/
                        /*SouthWest,    South,  SouthEast,*/

                        BehindNorthWest,    BehindNorth, /*BehindNorthEast,*/
                        BehindWest,         Behind,      /*BehindEast,*/
                        /*BehindSouthWest,  BehindSouth, BehindSouthEast*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error},    
                    //26 - AtBottomLeftFrontBorder
                    { /*InFrontNorthWest, InFrontNorth, InFrontNorthEast,
                        InFrontWest,      InFront,      InFrontEast,
                        InFrontSouthWest, InFrontSouth, InFrontSouthEast,*/
      
                        /*NorthWest,*/ North, NorthEast,
                        /*West,*/             East,
                        /*SouthWest,   South, SouthEast,*/

                        /*BehindNorthWest,*/ BehindNorth, BehindNorthEast,
                        /*BehindWest,*/      Behind,      BehindEast,
                        /*BehindSouthWest,   BehindSouth, BehindSouthEast*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error},    
                    //27 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //28 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //29 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //30 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //31 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},
                    //32 - AtRearBorder
                    {   InFrontNorthWest, InFrontNorth, InFrontNorthEast,
                        InFrontWest,      InFront,      InFrontEast,
                        InFrontSouthWest, InFrontSouth, InFrontSouthEast,
      
                        NorthWest, North, NorthEast,
                        West,             East,
                        SouthWest, South, SouthEast, 

                        /*BehindNorthWest, BehindNorth, BehindNorthEast,
                        BehindWest,        Behind,      BehindEast,
                        BehindSouthWest,   BehindSouth, BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},    
                    //33 - AtRearRightBorder
                    {   InFrontNorthWest, InFrontNorth, /*InFrontNorthEast,*/
                        InFrontWest,      InFront,      /*InFrontEast,*/
                        InFrontSouthWest, InFrontSouth, /*InFrontSouthEast,*/
      
                        NorthWest, North, /*NorthEast,*/
                        West,             /*East,*/
                        SouthWest, South, /*SouthEast,*/

                        /*BehindNorthWest, BehindNorth, BehindNorthEast,
                        BehindWest,        Behind,      BehindEast,
                        BehindSouthWest,   BehindSouth, BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},    
                    //34 - AtRearLeftBorder
                    {   /*InFrontNorthWest,*/ InFrontNorth, InFrontNorthEast,
                        /*InFrontWest,*/      InFront,      InFrontEast,
                        /*InFrontSouthWest,*/ InFrontSouth, InFrontSouthEast,
      
                        /*NorthWest,*/ North, NorthEast,
                        /*West,*/             East,
                        /*SouthWest,*/ South, SouthEast, 

                        /*BehindNorthWest, BehindNorth,   BehindNorthEast,
                        BehindWest,        Behind,        BehindEast,
                        BehindSouthWest,   BehindSouth,   BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},        
                    //35 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},
                    //36 - AtTopRearBorder
                    {   /*InFrontNorthWest, InFrontNorth,   InFrontNorthEast,*/
                        InFrontWest,        InFront,        InFrontEast,
                        InFrontSouthWest,   InFrontSouth,   InFrontSouthEast,
      
                        /*NorthWest,    North, NorthEast,*/
                        West,                  East,
                        SouthWest,      South, SouthEast, 

                        /*BehindNorthWest, BehindNorth,   BehindNorthEast,
                        BehindWest,        Behind,        BehindEast,
                        BehindSouthWest,   BehindSouth,   BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},            
                    //37 - AtTopRightRearBorder    
                    {   /*InFrontNorthWest, InFrontNorth,   InFrontNorthEast,*/
                        InFrontWest,        InFront,        /*InFrontEast,*/
                        InFrontSouthWest,   InFrontSouth,   /*InFrontSouthEast,*/
      
                        /*NorthWest, North, NorthEast,*/
                        West,               /*East,*/
                        SouthWest,   South, /*SouthEast,*/ 

                        /*BehindNorthWest,  BehindNorth, BehindNorthEast,
                        BehindWest,         Behind,      BehindEast,
                        BehindSouthWest,    BehindSouth, BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error},    
                    //38 - AtTopLeftRearBorder
                    {   /*InFrontNorthWest,     InFrontNorth,    InFrontNorthEast,*/
                        /*InFrontWest,*/        InFront,         InFrontEast,
                        /*InFrontSouthWest,*/   InFrontSouth,   InFrontSouthEast,
      
                        /*NorthWest,    North,  NorthEast,*/
                        /*West,*/               East,
                        /*SouthWest,*/  South,  SouthEast, 

                        /*BehindNorthWest,  BehindNorth,    BehindNorthEast,
                        BehindWest,         Behind,         BehindEast,
                        BehindSouthWest,    BehindSouth,    BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error},    
                    //39 - Nothin'
                    {   Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error, Error, Error, 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error},
                    //40 - AtBottomRearBorder
                    {   InFrontNorthWest,   InFrontNorth,   InFrontNorthEast,
                        InFrontWest,        InFront,        InFrontEast,
                        /*InFrontSouthWest, InFrontSouth,   InFrontSouthEast,*/
      
                        NorthWest,      North, NorthEast,
                        West,                  East,
                        /*SouthWest,    South, SouthEast,*/ 

                        /*BehindNorthWest,  BehindNorth, BehindNorthEast,
                        BehindWest,         Behind,      BehindEast,
                        BehindSouthWest,    BehindSouth, BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error},    
                    //41 - AtBottomRightRearBorder
                    {   InFrontNorthWest,   InFrontNorth, /*InFrontNorthEast,*/
                        InFrontWest,        InFront,      /*InFrontEast,*/
                        /*InFrontSouthWest, InFrontSouth, InFrontSouthEast,*/
      
                        NorthWest,   North, /*NorthEast,*/
                        West,               /*East,*/
                        /*SouthWest, South, SouthEast,*/ 

                        /*BehindNorthWest, BehindNorth, BehindNorthEast,
                        BehindWest,        Behind,      BehindEast,
                        BehindSouthWest,   BehindSouth, BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error},    
                    //42 - AtBottomLeftRearBorder
                    {   /*InFrontNorthWest,*/   InFrontNorth,   InFrontNorthEast,
                        /*InFrontWest,*/        InFront,        InFrontEast,
                        /*InFrontSouthWest,     InFrontSouth,   InFrontSouthEast,*/
      
                        /*NorthWest,*/  North,  NorthEast,
                        /*West,*/               East,
                        /*SouthWest,    South,  SouthEast,*/ 

                        /*BehindNorthWest,  BehindNorth, BehindNorthEast,
                        BehindWest,         Behind,      BehindEast,
                        BehindSouthWest,    BehindSouth, BehindSouthEast,*/ 
                        Error, Error, Error, Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error, Error, Error,
                        Error, Error, Error, Error}    
               };
        return c[b][index];
    }
        
    /** The valid direction three codes in anti causal direction (means: look back in scanline
        direction)when the center is at the volume border.
            Should be used with isAtVolumeBorderCausal to determine the Directions, as this
            avoids using of the nonesense border ids (e.g. 0,1,8,9...) of this table.
        \a index must be in the range <tt>0...nearBorderDirectionCount(b)-1</tt>.
     */
    static Direction nearBorderDirectionsCausal(AtVolumeBorder b, int index)
    {
        static Direction c[43][13] = {
            //0 - NotAtBorder -----> should never be used
                                { InFrontNorthWest, InFrontNorth,    InFrontNorthEast,
                                  InFrontWest,      InFront,         InFrontEast,
                                  InFrontSouthWest, InFrontSouth,    InFrontSouthEast,
                  
                                  NorthWest,        North,           NorthEast,
                                  West},                    
            //1 - AtRightBorder -----> should never be used
                                { InFrontNorthWest, InFrontNorth,    InFrontNorthEast,
                                  InFrontWest,      InFront,         InFrontEast,
                                  InFrontSouthWest, InFrontSouth,    InFrontSouthEast,
                  
                                  NorthWest,        North,           NorthEast,
                                  West},    
            //2 - AtLeftBorder
                                { /*InFrontNorthWest,*/ InFrontNorth,InFrontNorthEast,
                                  /*InFrontWest,*/  InFront,         InFrontEast,
                                  /*InFrontSouthWest,*/InFrontSouth, InFrontSouthEast,
                  
                                  /*NorthWest,*/    North,           NorthEast,
                                  /*West*/
                                  Error, Error, Error, Error, Error},    
            //3 - Nothin'
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //4 - AtTopBorder
                                { /*InFrontNorthWest,InFrontNorth,   InFrontNorthEast,*/
                                  InFrontWest,       InFront,        InFrontEast,
                                  InFrontSouthWest,  InFrontSouth,   InFrontSouthEast,
                  
                                  /*NorthWest,       North,          NorthEast,*/
                                  West,
                                  Error, Error, Error, Error, Error, Error},    
            //5 - AtTopRightBorder
                                { /*InFrontNorthWest,InFrontNorth,   InFrontNorthEast,*/
                                  InFrontWest,       InFront,        /*InFrontEast,*/
                                  InFrontSouthWest,  InFrontSouth,   /*InFrontSouthEast,*/
                  
                                  /*NorthWest, North, NorthEast,*/
                                  West,
                                  Error, Error, Error, Error, Error, Error, Error, Error},    
            //6 - AtTopLeftBorder
                                { /*InFrontNorthWest,InFrontNorth,    InFrontNorthEast,*/
                                  /*InFrontWest,*/   InFront,         InFrontEast,
                                  /*InFrontSouthWest,*/InFrontSouth,  InFrontSouthEast,
                  
                                  /*NorthWest,       North,           NorthEast,*/
                                  /*West,*/
                                  Error, Error, Error, Error, Error, Error, Error, Error, Error},
            //7 - Nothin'
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //8 - AtBottomBorder -----> should never be used
                                { InFrontNorthWest,  InFrontNorth,    InFrontNorthEast,
                                  InFrontWest,       InFront,         InFrontEast,
                                  InFrontSouthWest,  InFrontSouth,    InFrontSouthEast,
                  
                                  NorthWest,         North,           NorthEast,
                                  West},    
            //9 - AtBottomRightBorder -----> should never be used
                                { InFrontNorthWest, InFrontNorth,    InFrontNorthEast,
                                  InFrontWest,      InFront,         InFrontEast,
                                  InFrontSouthWest, InFrontSouth,    InFrontSouthEast,
                  
                                  NorthWest,        North,           NorthEast,
                                  West},    
            //10 - AtBottomLeftBorder
                                { /*InFrontNorthWest,*/InFrontNorth, InFrontNorthEast,
                                  /*InFrontWest,*/  InFront,         InFrontEast,
                                  /*InFrontSouthWest,*/InFrontSouth, InFrontSouthEast,
                  
                                  /*NorthWest,*/    North,           NorthEast,
                                  /*West*/
                                  Error, Error, Error, Error, Error},    
            //11 - Nothin'
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //12 - Nothin'
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //13 - Nothin'
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //14 - Nothin'
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //15 - Nothin'
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //16 - AtFrontBorder
                                { /*InFrontNorthWest,InFrontNorth,   InFrontNorthEast,
                                  InFrontWest,      InFront,         InFrontEast,
                                  InFrontSouthWest, InFrontSouth,    InFrontSouthEast,*/
                  
                                  NorthWest,        North,           NorthEast,
                                  West,
                                  Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //17 - AtFrontRightBorder
                                { /*InFrontNorthWest,InFrontNorth,   InFrontNorthEast,
                                  InFrontWest,      InFront,         InFrontEast,
                                  InFrontSouthWest, InFrontSouth,    InFrontSouthEast,*/
                  
                                  NorthWest,        North,           /*NorthEast,*/
                                  West,
                                  Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //18 - AtFrontLeftBorder
                                { /*InFrontNorthWest,InFrontNorth,   InFrontNorthEast,
                                  InFrontWest,      InFront,         InFrontEast,
                                  InFrontSouthWest, InFrontSouth,    InFrontSouthEast,*/
                  
                                  /*NorthWest,*/    North,           NorthEast,
                                  /*West,*/
                                  Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //19 - Nothin'
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //20 - AtTopFrontBorder
                                { /*InFrontNorthWest,InFrontNorth,   InFrontNorthEast,
                                  InFrontWest,      InFront,         InFrontEast,
                                  InFrontSouthWest, InFrontSouth,    InFrontSouthEast,*/
                  
                                  /*NorthWest, North, NorthEast,*/
                                  West,
                                  Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //21 - AtTopRightFrontBorder
                                { /*InFrontNorthWest, InFrontNorth,  InFrontNorthEast,
                                  InFrontWest,        InFront,       InFrontEast,
                                  InFrontSouthWest,   InFrontSouth,  InFrontSouthEast,*/
                  
                                  /*NorthWest,        North,         NorthEast,*/
                                  West,
                                  Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //22 - AtTopLeftFrontBorder
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //23 - Nothin
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //24 - AtBottomFrontBorder
                                { /*InFrontNorthWest, InFrontNorth,  InFrontNorthEast,
                                  InFrontWest,        InFront,       InFrontEast,
                                  InFrontSouthWest,   InFrontSouth,  InFrontSouthEast,*/
                  
                                  NorthWest,          North,         NorthEast,
                                  West,
                                  Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //25 - AtBottomRightFrontBorder
                                { /*InFrontNorthWest, InFrontNorth,  InFrontNorthEast,
                                  InFrontWest,        InFront,       InFrontEast,
                                  InFrontSouthWest,   InFrontSouth,  InFrontSouthEast,*/
                  
                                  NorthWest,          North,         /*NorthEast,*/
                                  West,
                                  Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //26 - AtBottomLeftFrontBorder
                                { /*InFrontNorthWest, InFrontNorth,  InFrontNorthEast,
                                  InFrontWest,        InFront,       InFrontEast,
                                  InFrontSouthWest,   InFrontSouth,  InFrontSouthEast,*/
                  
                                  /*NorthWest,*/      North,         NorthEast,
                                  West,
                                  Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //27 - Nothin
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //28 - Nothin
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //29 - Nothin
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //30 - Nothin
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //31 - Nothin
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //32 - AtRearBorder -----> should never be used
                                { InFrontNorthWest,   InFrontNorth,  InFrontNorthEast,
                                  InFrontWest,        InFront,       InFrontEast,
                                  InFrontSouthWest,   InFrontSouth,  InFrontSouthEast,
                  
                                  NorthWest,          North,         NorthEast,
                                  West},    
            //33 - AtRearRightBorder -----> should never be used
                                { InFrontNorthWest,   InFrontNorth,  InFrontNorthEast,
                                  InFrontWest,        InFront,       InFrontEast,
                                  InFrontSouthWest,   InFrontSouth,  InFrontSouthEast,
                  
                                  NorthWest,          North,         NorthEast,
                                  West},    
            //34 - AtRearLeftBorder
                                { /*InFrontNorthWest,*/InFrontNorth, InFrontNorthEast,
                                  /*InFrontWest,*/    InFront,       InFrontEast,
                                  /*InFrontSouthWest,*/InFrontSouth, InFrontSouthEast,
                  
                                  /*NorthWest,*/      North,         NorthEast,
                                  /*West*/
                                  Error, Error, Error, Error, Error},    
            //35 - Nothin
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //36 - AtTopRearBorder    
                                { /*InFrontNorthWest, InFrontNorth,  InFrontNorthEast,*/
                                  InFrontWest,        InFront,       InFrontEast,
                                  InFrontSouthWest,   InFrontSouth,  InFrontSouthEast,
                  
                                  /*NorthWest,        North,         NorthEast,*/
                                  West,
                                  Error, Error, Error, Error, Error, Error},    
            //37 - AtTopRightRearBorder        
                                { /*InFrontNorthWest, InFrontNorth,  InFrontNorthEast,*/
                                  InFrontWest,        InFront,       /*InFrontEast,*/
                                  InFrontSouthWest,   InFrontSouth,  /*InFrontSouthEast,*/
                  
                                  /*NorthWest,        North,         NorthEast,*/
                                  West,
                                  Error, Error, Error, Error, Error, Error, Error, Error},
            //38 - AtTopLeftRearBorder    
                                { /*InFrontNorthWest, InFrontNorth,  InFrontNorthEast,*/
                                  /*InFrontWest,*/    InFront,       InFrontEast,
                                  /*InFrontSouthWest,*/InFrontSouth, InFrontSouthEast,
                  
                                  /*NorthWest, North, NorthEast,*/
                                  /*West,*/
                                  Error, Error, Error, Error, Error, Error, Error, Error, Error},
            //39 - Nothin
                                { Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error, Error},    
            //40 - AtBottomRearBorder -----> should never be used
                                { InFrontNorthWest,   InFrontNorth,  InFrontNorthEast,
                                  InFrontWest,        InFront,       InFrontEast,
                                  InFrontSouthWest,   InFrontSouth,  InFrontSouthEast,
                  
                                  NorthWest,          North,         NorthEast,
                                  West},    
            //41 - AtBottomRightRearBorder -----> should never be used
                                { InFrontNorthWest,  InFrontNorth,   InFrontNorthEast,
                                  InFrontWest,       InFront,        InFrontEast,
                                  InFrontSouthWest,  InFrontSouth,   InFrontSouthEast,
                  
                                  NorthWest,         North,          NorthEast,
                                  West},    
            //42 - AtBottomLeftRearBorder
                                { /*InFrontNorthWest,*/InFrontNorth, InFrontNorthEast,
                                  /*InFrontWest,*/   InFront,        InFrontEast,
                                  /*InFrontSouthWest,InFrontSouth,   InFrontSouthEast,*/
                  
                                  /*NorthWest,*/     North,          NorthEast,
                                  /*West*/
                                  Error, Error, Error, Error, Error, Error, Error}    
        };
        return c[b][index];
    }
    
    /** transform direction code into corresponding Diff3D offset.
        (note: there is no bounds checking on the code you pass.)
    */
    static Diff3D const & diff(Direction code)
    {
        static Diff3D d[] = {   Diff3D( -1, -1, -1),  //InFrontNorthWest
                                Diff3D(  0, -1, -1),  //InFrontNorth
                                Diff3D(  1, -1, -1),  //InFrontNorthEast
                                Diff3D( -1,  0, -1),  //InFrontWest
                                Diff3D(  0,  0, -1),  //InFront
                                Diff3D(  1,  0, -1),  //InFrontEast
                                Diff3D( -1,  1, -1),  //InFrontSouthWest
                                Diff3D(  0,  1, -1),  //InFrontSouth
                                Diff3D(  1,  1, -1),  //InFrontSouthEast

                                Diff3D( -1, -1,  0),  //NorthWest
                                Diff3D(  0, -1,  0),  //North
                                Diff3D(  1, -1,  0),  //NorthEast
                                Diff3D( -1,  0,  0),  //West
                                Diff3D(  1,  0,  0),  //East
                                Diff3D( -1,  1,  0),  //SouthWest
                                Diff3D(  0,  1,  0),  //South
                                Diff3D(  1,  1,  0),  //SouthEast

                                Diff3D( -1, -1,  1),  //BehindNorthWest
                                Diff3D(  0, -1,  1),  //BehindNorth
                                Diff3D(  1, -1,  1),  //BehindNorthEast
                                Diff3D( -1,  0,  1),  //BehindWest
                                Diff3D(  0,  0,  1),  //Behind
                                Diff3D(  1,  0,  1),  //BehindEast
                                Diff3D( -1,  1,  1),  //BehindSouthWest
                                Diff3D(  0,  1,  1),  //BehindSouth
                                Diff3D(  1,  1,  1),  //BehindSouthEast
                            };
        return d[code];
    }

    /** Equivalent to <tt>diff(static_cast<Direction>(code))</tt>.
        (note: there is no bounds checking on the code you pass.)
    */
    static Diff3D const & diff(int code) { return diff(static_cast<Direction>(code)); }

    /**  Equivalent to <tt>diff(code)[dim]</tt> */
    static int diff(Direction code, int dim) { return diff(code)[dim]; }
   
    /** Get the relative offset from one neighbor to the other.
    For example, <tt>relativeDiff(East, West) == multi_differencetype(-2,0,0)</tt>.
    (note: there is no bounds checking on the code you pass.)
    */
    static Diff3D const relativeDiff(Direction fromCode, Direction toCode)
    {
        //Uncomment the following lines may cause the program to crash because of insufficient
        //static allocatable memory on the stack
        /*
            static Diff3D d[][] = {
                //   InFront-NW    ---     InFront-N   ---     InFront-NE  ---    Infront-W    ---     InFront     ---   InFront-E     ---   InFront-SW    ---      InFront-S  ---   InFront-SE    ---    NorthWest   ---       North      ---    NorthEast    ---      West      ---       East       ---    SouthWest    ---    South        ---     SouthEast  ---     Behind-NW    ---     Behind-N    ---      Behind-NE  ---     Behind-W    ---      Behind     ---    Behind-E     ---    Behind-SW    ---      Behind-S   ---    Behind-SE  
                {    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D( 2,  0,  0),    Diff3D( 0,  1,  0),    Diff3D( 1,  1,  0),    Diff3D( 2,  1,  0),    Diff3D( 0,  2,  0),    Diff3D( 1,  2,  0),    Diff3D( 2, 2,  0),    Diff3D( 0,  0,  1),    Diff3D( 1,  0,  1),    Diff3D( 2,  0,  1),    Diff3D( 0,  1,  1),    Diff3D( 2,  1,  1),    Diff3D( 0,  2,  1),    Diff3D( 1,  2,  1),    Diff3D( 2,  2,  1),    Diff3D( 0,  0,  2),    Diff3D( 1,  0,  2),    Diff3D( 2,  0,  2),    Diff3D( 0,  1,  2),    Diff3D( 1,  1,  2), Diff3D( 2,  1,  2),    Diff3D( 0,  2,  2),    Diff3D( 1,  2,  2),    Diff3D( 2,  2,  2)    },    //InFront-NW  
                {    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D(-1,  1,  0),    Diff3D( 0,  1,  0),    Diff3D( 1,  1,  0),    Diff3D(-1,  2,  0),    Diff3D( 0,  2,  0),    Diff3D( 1, 2,  0),    Diff3D(-1,  0,  1),    Diff3D( 0,  0,  1),    Diff3D( 1,  0,  1),    Diff3D(-1,  1,  1),    Diff3D( 1,  1,  1),    Diff3D(-1,  2,  1),    Diff3D( 0,  2,  1),    Diff3D( 1,  2,  1),    Diff3D(-1,  0,  2),    Diff3D( 0,  0,  2),    Diff3D( 1,  0,  2),    Diff3D(-1,  1,  2),    Diff3D( 0,  1,  2), Diff3D( 1,  1,  2),    Diff3D(-1,  2,  2),    Diff3D( 0,  2,  2),    Diff3D( 1,  2,  2)    },    //InFront-N  
                {    Diff3D(-2,  0,  0),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D(-2,  1,  0),    Diff3D(-1,  1,  0),    Diff3D( 0,  1,  0),    Diff3D(-2,  2,  0),    Diff3D(-1,  2,  0),    Diff3D( 0, 2,  0),    Diff3D(-2,  0,  1),    Diff3D(-1,  0,  1),    Diff3D( 0,  0,  1),    Diff3D(-2,  1,  1),    Diff3D( 0,  1,  1),    Diff3D(-2,  2,  1),    Diff3D(-1,  2,  1),    Diff3D( 0,  2,  1),    Diff3D(-2,  0,  2),    Diff3D(-1,  0,  2),    Diff3D( 0,  0,  2),    Diff3D(-2,  1,  2),    Diff3D(-1,  1,  2), Diff3D( 0,  1,  2),    Diff3D(-2,  2,  2),    Diff3D(-1,  2,  2),    Diff3D( 0,  2,  2)    },    //InFront-NE
                {    Diff3D(0,  -1,  0),    Diff3D( 1, -1,  0),    Diff3D( 2, -1,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D( 2,  0,  0),    Diff3D( 0,  1,  0),    Diff3D( 1,  1,  0),    Diff3D( 2, 1,  0),    Diff3D( 0, -1,  1),    Diff3D( 1, -1,  1),    Diff3D( 2, -1,  1),    Diff3D( 0,  0,  1),    Diff3D( 2,  0,  1),    Diff3D( 0,  1,  1),    Diff3D( 1,  1,  1),    Diff3D( 2,  1,  1),    Diff3D( 0, -1,  2),    Diff3D( 1, -1,  2),    Diff3D( 2, -1,  2),    Diff3D( 0,  0,  2),    Diff3D( 1,  0,  2), Diff3D( 2,  0,  2),    Diff3D( 0,  1,  2),    Diff3D( 1,  1,  2),    Diff3D( 2,  1,  2)    },    //Infront-W
                {    Diff3D(-1, -1,  0),    Diff3D( 0, -1,  0),    Diff3D( 1, -1,  0),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D(-1,  1,  0),    Diff3D( 0,  1,  0),    Diff3D( 1, 1,  0),    Diff3D(-1, -1,  1),    Diff3D( 0, -1,  1),    Diff3D( 1, -1,  1),    Diff3D(-1,  0,  1),    Diff3D( 1,  0,  1),    Diff3D(-1,  1,  1),    Diff3D( 0,  1,  1),    Diff3D( 1,  1,  1),    Diff3D(-1, -1,  2),    Diff3D( 0, -1,  2),    Diff3D( 1, -1,  2),    Diff3D(-1,  0,  2),    Diff3D( 0,  0,  2), Diff3D( 1,  0,  2),    Diff3D(-1,  1,  2),    Diff3D( 0,  1,  2),    Diff3D( 1,  1,  2)    },    //InFront
                {    Diff3D(-2, -1,  0),    Diff3D(-1, -1,  0),    Diff3D( 0, -1,  0),    Diff3D(-2,  0,  0),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D(-2,  1,  0),    Diff3D(-1,  1,  0),    Diff3D( 0, 1,  0),    Diff3D(-2, -1,  1),    Diff3D(-1, -1,  1),    Diff3D( 0, -1,  1),    Diff3D(-2,  0,  1),    Diff3D( 0,  0,  1),    Diff3D(-2,  1,  1),    Diff3D(-1,  1,  1),    Diff3D( 0,  1,  1),    Diff3D(-2, -1,  2),    Diff3D(-1, -1,  2),    Diff3D( 0, -1,  2),    Diff3D(-2,  0,  2),    Diff3D(-1,  0,  2), Diff3D( 0,  0,  2),    Diff3D(-2,  1,  2),    Diff3D(-1,  1,  2),    Diff3D( 0,  1,  2)    },    //InFront-E
                {    Diff3D( 0, -2,  0),    Diff3D( 1, -2,  0),    Diff3D( 2, -2,  0),    Diff3D( 0, -1,  0),    Diff3D( 1, -1,  0),    Diff3D( 2, -1,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D( 2, 0,  0),    Diff3D( 0, -2,  1),    Diff3D( 1, -2,  1),    Diff3D( 2, -2,  1),    Diff3D( 0, -1,  1),    Diff3D( 2, -1,  1),    Diff3D( 0,  0,  1),    Diff3D( 1,  0,  1),    Diff3D( 2,  0,  1),    Diff3D( 0, -2,  2),    Diff3D( 1, -2,  2),    Diff3D( 2, -2,  2),    Diff3D( 0, -1,  2),    Diff3D( 1, -1,  2), Diff3D( 2, -1,  2),    Diff3D( 0,  0,  2),    Diff3D( 1,  0,  2),    Diff3D( 2,  0,  2)    },    //InFront-SW
                {    Diff3D(-1, -2,  0),    Diff3D( 0, -2,  0),    Diff3D( 1, -2,  0),    Diff3D(-1, -1,  0),    Diff3D( 0, -1,  0),    Diff3D( 1, -1,  0),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D( 1, 0,  0),    Diff3D(-1, -2,  1),    Diff3D( 0, -2,  1),    Diff3D( 1, -2,  1),    Diff3D(-1, -1,  1),    Diff3D( 1, -1,  1),    Diff3D(-1,  0,  1),    Diff3D( 0,  0,  1),    Diff3D( 1,  0,  1),    Diff3D(-1, -2,  2),    Diff3D( 0, -2,  2),    Diff3D( 1, -2,  2),    Diff3D(-1, -1,  2),    Diff3D( 0, -1,  2), Diff3D( 1, -1,  2),    Diff3D(-1,  0,  2),    Diff3D( 0,  0,  2),    Diff3D( 1,  0,  2)    },    //InFront-S 
                {    Diff3D(-2, -2,  0),    Diff3D(-1, -2,  0),    Diff3D( 0, -2,  0),    Diff3D(-2, -1,  0),    Diff3D(-1, -1,  0),    Diff3D( 0, -1,  0),    Diff3D(-2,  0,  0),    Diff3D(-1,  0,  0),    Diff3D( 0, 0,  0),    Diff3D(-2, -2,  1),    Diff3D(-1, -2,  1),    Diff3D( 0, -2,  1),    Diff3D(-2, -1,  1),    Diff3D( 0, -1,  1),    Diff3D(-2,  0,  1),    Diff3D(-1,  0,  1),    Diff3D( 0,  0,  1),    Diff3D(-2, -2,  2),    Diff3D(-1, -2,  2),    Diff3D( 0, -2,  2),    Diff3D(-2, -1,  2),    Diff3D(-1, -1,  2), Diff3D( 0, -1,  2),    Diff3D(-2,  0,  2),    Diff3D(-1,  0,  2),    Diff3D( 0,  0,  2)    },    //InFront-SE
                {    Diff3D( 0,  0, -1),    Diff3D( 1,  0, -1),    Diff3D( 2,  0, -1),    Diff3D( 0,  1, -1),    Diff3D( 1,  1, -1),    Diff3D( 2,  1, -1),    Diff3D( 0,  2, -1),    Diff3D( 1,  2, -1),    Diff3D( 2, 2, -1),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D( 2,  0,  0),    Diff3D( 0,  1,  0),    Diff3D( 2,  1,  0),    Diff3D( 0,  2,  0),    Diff3D( 1,  2,  0),    Diff3D( 2,  2,  0),    Diff3D( 0,  0,  1),    Diff3D( 1,  0,  1),    Diff3D( 2,  0,  1),    Diff3D( 0,  1,  1),    Diff3D( 1,  1,  1), Diff3D( 2,  1,  1),    Diff3D( 0,  2,  1),    Diff3D( 1,  2,  1),    Diff3D( 2,  2,  1)    },    //NorthWest
                {    Diff3D(-1,  0, -1),    Diff3D( 0,  0, -1),    Diff3D( 1,  0, -1),    Diff3D(-1,  1, -1),    Diff3D( 0,  1, -1),    Diff3D( 1,  1, -1),    Diff3D(-1,  2, -1),    Diff3D( 0,  2, -1),    Diff3D( 1, 2, -1),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D(-1,  1,  0),    Diff3D( 1,  1,  0),    Diff3D(-1,  2,  0),    Diff3D( 0,  2,  0),    Diff3D( 1,  2,  0),    Diff3D(-1,  0,  1),    Diff3D( 0,  0,  1),    Diff3D( 1,  0,  1),    Diff3D(-1,  1,  1),    Diff3D( 0,  1,  1), Diff3D( 1,  1,  1),    Diff3D(-1,  2,  1),    Diff3D( 0,  2,  1),    Diff3D( 1,  2,  1)    },    //North
                {    Diff3D(-2,  0, -1),    Diff3D(-1,  0, -1),    Diff3D( 0,  0, -1),    Diff3D(-2,  1, -1),    Diff3D(-1,  1, -1),    Diff3D( 0,  1, -1),    Diff3D(-2,  2, -1),    Diff3D(-1,  2, -1),    Diff3D( 0, 2, -1),    Diff3D(-2,  0,  0),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D(-2,  1,  0),    Diff3D( 0,  1,  0),    Diff3D(-2,  2,  0),    Diff3D(-1,  2,  0),    Diff3D( 0,  2,  0),    Diff3D(-2,  0,  1),    Diff3D(-1,  0,  1),    Diff3D( 0,  0,  1),    Diff3D(-2,  1,  1),    Diff3D(-1,  1,  1), Diff3D( 0,  1,  1),    Diff3D(-2,  2,  1),    Diff3D(-1,  2,  1),    Diff3D( 0,  2,  1)    },    //NortEast
                {    Diff3D( 0, -1, -1),    Diff3D( 1, -1, -1),    Diff3D( 2, -1, -1),    Diff3D( 0,  0, -1),    Diff3D( 1,  0, -1),    Diff3D( 2,  0, -1),    Diff3D( 0,  1, -1),    Diff3D( 1,  1, -1),    Diff3D( 2, 1, -1),    Diff3D( 0, -1,  0),    Diff3D( 1, -1,  0),    Diff3D( 2, -1,  0),    Diff3D( 0,  0,  0),    Diff3D( 2,  0,  0),    Diff3D( 0,  1,  0),    Diff3D( 1,  1,  0),    Diff3D( 2,  1,  0),    Diff3D( 0, -1,  1),    Diff3D( 1, -1,  1),    Diff3D( 2, -1,  1),    Diff3D( 0,  0,  1),    Diff3D( 1,  0,  1), Diff3D( 2,  0,  1),    Diff3D( 0,  1,  1),    Diff3D( 1,  1,  1),    Diff3D( 2,  1,  1)    },    //West
                {    Diff3D(-2, -1, -1),    Diff3D(-1, -1, -1),    Diff3D( 0, -1, -1),    Diff3D(-2,  0, -1),    Diff3D(-1,  0, -1),    Diff3D( 0,  0, -1),    Diff3D(-2,  1, -1),    Diff3D(-1,  1, -1),    Diff3D( 0, 1, -1),    Diff3D(-2, -1,  0),    Diff3D(-1, -1,  0),    Diff3D( 0, -1,  0),    Diff3D(-2,  0,  0),    Diff3D( 0,  0,  0),    Diff3D(-2,  1,  0),    Diff3D(-1,  1,  0),    Diff3D( 0,  1,  0),    Diff3D(-2, -1,  1),    Diff3D(-1, -1,  1),    Diff3D( 0, -1,  1),    Diff3D(-2,  0,  1),    Diff3D(-1,  0,  1), Diff3D( 0,  0,  1),    Diff3D(-2,  1,  1),    Diff3D(-1,  1,  1),    Diff3D( 0,  1,  1)    },    //East
                {    Diff3D( 0, -2, -1),    Diff3D( 1, -2, -1),    Diff3D( 2, -2, -1),    Diff3D( 0, -1, -1),    Diff3D( 1, -1, -1),    Diff3D( 2, -1, -1),    Diff3D( 0,  0, -1),    Diff3D( 1,  0, -1),    Diff3D( 2, 0, -1),    Diff3D( 0, -2,  0),    Diff3D( 1, -2,  0),    Diff3D( 2, -2,  0),    Diff3D( 0, -1,  0),    Diff3D( 2, -1,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D( 2,  0,  0),    Diff3D( 0, -2,  1),    Diff3D( 1, -2,  1),    Diff3D( 2, -2,  1),    Diff3D( 0, -1,  1),    Diff3D( 1, -1,  1), Diff3D( 2, -1,  1),    Diff3D( 0,  0,  1),    Diff3D( 1,  0,  1),    Diff3D( 2,  0,  1)    },    //SouthWest
                {    Diff3D(-1, -2, -1),    Diff3D( 0, -2, -1),    Diff3D( 1, -2, -1),    Diff3D(-1, -1, -1),    Diff3D( 0, -1, -1),    Diff3D( 1, -1, -1),    Diff3D(-1,  0, -1),    Diff3D( 0,  0, -1),    Diff3D( 1, 0, -1),    Diff3D(-1, -2,  0),    Diff3D( 0, -2,  0),    Diff3D( 1, -2,  0),    Diff3D(-1, -1,  0),    Diff3D( 1, -1,  0),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D(-1, -2,  1),    Diff3D( 0, -2,  1),    Diff3D( 1, -2,  1),    Diff3D(-1, -1,  1),    Diff3D( 0, -1,  1), Diff3D( 1, -1,  1),    Diff3D(-1,  0,  1),    Diff3D( 0,  0,  1),    Diff3D( 1,  0,  1)    },    //South
                {    Diff3D(-2, -2, -1),    Diff3D(-1, -2, -1),    Diff3D( 0, -2, -1),    Diff3D(-2, -1, -1),    Diff3D(-1, -1, -1),    Diff3D( 0, -1, -1),    Diff3D(-2,  0, -1),    Diff3D(-1,  0, -1),    Diff3D( 0, 0, -1),    Diff3D(-2, -2,  0),    Diff3D(-1, -2,  0),    Diff3D( 0, -2,  0),    Diff3D(-2, -1,  0),    Diff3D( 0, -1,  0),    Diff3D(-2,  0,  0),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D(-2, -2,  1),    Diff3D(-1, -2,  1),    Diff3D( 0, -2,  1),    Diff3D(-2, -1,  1),    Diff3D(-1, -1,  1), Diff3D( 0, -1,  1),    Diff3D(-2,  0,  1),    Diff3D(-1,  0,  1),    Diff3D( 0,  0,  1)    },    //SouthEast
                {    Diff3D( 0,  0, -2),    Diff3D( 1,  0, -2),    Diff3D( 2,  0, -2),    Diff3D( 0,  1, -2),    Diff3D( 1,  1, -2),    Diff3D( 2,  1, -2),    Diff3D( 0,  2, -2),    Diff3D( 1,  2, -2),    Diff3D( 2, 2, -2),    Diff3D( 0,  0, -1),    Diff3D( 1,  0, -1),    Diff3D( 2,  0, -1),    Diff3D( 0,  1, -1),    Diff3D( 2,  1, -1),    Diff3D( 0,  2, -1),    Diff3D( 1,  2, -1),    Diff3D( 2,  2, -1),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D( 2,  0,  0),    Diff3D( 0,  1,  0),    Diff3D( 1,  1,  0),    Diff3D( 2,  1,  0),    Diff3D( 0,  2,  0),    Diff3D( 1,  2,  0),    Diff3D( 2,  2,  0)    },    //Behind-NW
                {    Diff3D(-1,  0, -2),    Diff3D( 0,  0, -2),    Diff3D( 1,  0, -2),    Diff3D(-1,  1, -2),    Diff3D( 0,  1, -2),    Diff3D( 1,  1, -2),    Diff3D(-1,  2, -2),    Diff3D( 0,  2, -2),    Diff3D( 1, 2, -2),    Diff3D(-1,  0, -1),    Diff3D( 0,  0, -1),    Diff3D( 1,  0, -1),    Diff3D(-1,  1, -1),    Diff3D( 1,  1, -1),    Diff3D(-1,  2, -1),    Diff3D( 0,  2, -1),    Diff3D( 1,  2, -1),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D(-1,  1,  0),    Diff3D( 0,  1,  0),    Diff3D( 1,  1,  0),    Diff3D(-1,  2,  0),    Diff3D( 0,  2,  0),    Diff3D( 1,  2,  0)    },    //Behind-N
                {    Diff3D(-2,  0, -2),    Diff3D(-1,  0, -2),    Diff3D( 0,  0, -2),    Diff3D(-2,  1, -2),    Diff3D(-1,  1, -2),    Diff3D( 0,  1, -2),    Diff3D(-2,  2, -2),    Diff3D(-1,  2, -2),    Diff3D( 0, 2, -2),    Diff3D(-2,  0, -1),    Diff3D(-1,  0, -1),    Diff3D( 0,  0, -1),    Diff3D(-2,  1, -1),    Diff3D( 0,  1, -1),    Diff3D(-2,  2, -1),    Diff3D(-1,  2, -1),    Diff3D( 0,  2, -1),    Diff3D(-2,  0,  0),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D(-2,  1,  0),    Diff3D(-1,  1,  0),    Diff3D( 0,  1,  0),    Diff3D(-2,  2,  0),    Diff3D(-1,  2,  0),    Diff3D( 0,  2,  0)    },    //Behind-NE
                {    Diff3D( 0, -1, -2),    Diff3D( 1, -1, -2),    Diff3D( 2, -1, -2),    Diff3D( 0,  0, -2),    Diff3D( 1,  0, -2),    Diff3D( 2,  0, -2),    Diff3D( 0,  1, -2),    Diff3D( 1,  1, -2),    Diff3D( 2, 1, -2),    Diff3D( 0, -1, -1),    Diff3D( 1, -1, -1),    Diff3D( 2, -1, -1),    Diff3D( 0,  0, -1),    Diff3D( 2,  0, -1),    Diff3D( 0,  1, -1),    Diff3D( 1,  1, -1),    Diff3D( 2,  1, -1),    Diff3D( 0, -1,  0),    Diff3D( 1, -1,  0),    Diff3D( 2, -1,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0), Diff3D( 2,  0,  0),    Diff3D( 0,  1,  0),    Diff3D( 1,  1,  0),    Diff3D( 2,  1,  0)    },    //Behind-W
                {    Diff3D(-1, -1, -2),    Diff3D( 0, -1, -2),    Diff3D( 1, -1, -2),    Diff3D(-1,  0, -2),    Diff3D( 0,  0, -2),    Diff3D( 1,  0, -2),    Diff3D(-1,  1, -2),    Diff3D( 0,  1, -2),    Diff3D( 1, 1, -2),    Diff3D(-1, -1, -1),    Diff3D( 0, -1, -1),    Diff3D( 1, -1, -1),    Diff3D(-1,  0, -1),    Diff3D( 1,  0, -1),    Diff3D(-1,  1, -1),    Diff3D( 0,  1, -1),    Diff3D( 1,  1, -1),    Diff3D(-1, -1,  0),    Diff3D( 0, -1,  0),    Diff3D( 1, -1,  0),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0), Diff3D( 1,  0,  0),    Diff3D(-1,  1,  0),    Diff3D( 0,  1,  0),    Diff3D( 1,  1,  0)    },    //Behind
                {    Diff3D(-2, -1, -2),    Diff3D(-1, -1, -2),    Diff3D( 0, -1, -2),    Diff3D(-2,  0, -2),    Diff3D(-1,  0, -2),    Diff3D( 0,  0, -2),    Diff3D(-2,  1, -2),    Diff3D(-1,  1, -2),    Diff3D( 0, 1, -2),    Diff3D(-2, -1, -1),    Diff3D(-1, -1, -1),    Diff3D( 0, -1, -1),    Diff3D(-2,  0, -1),    Diff3D( 0,  0, -1),    Diff3D(-2,  1, -1),    Diff3D(-1,  1, -1),    Diff3D( 0,  1, -1),    Diff3D(-2, -1,  0),    Diff3D(-1, -1,  0),    Diff3D( 0, -1,  0),    Diff3D(-2,  0,  0),    Diff3D(-1,  0,  0), Diff3D( 0,  0,  0),    Diff3D(-2,  1,  0),    Diff3D(-1,  1,  0),    Diff3D( 0,  1,  0)    },    //Behind-E
                {    Diff3D( 0, -2, -2),    Diff3D( 1, -2, -2),    Diff3D( 2, -2, -2),    Diff3D( 0, -1, -2),    Diff3D( 1, -1, -2),    Diff3D( 2, -1, -2),    Diff3D( 0,  0, -2),    Diff3D( 1,  0, -2),    Diff3D( 2, 0, -2),    Diff3D( 0, -2, -1),    Diff3D( 1, -2, -1),    Diff3D( 2, -2, -1),    Diff3D( 0, -1, -1),    Diff3D( 2, -1, -1),    Diff3D( 0,  0, -1),    Diff3D( 1,  0, -1),    Diff3D( 2,  0, -1),    Diff3D( 0, -2,  0),    Diff3D( 1, -2,  0),    Diff3D( 2, -2,  0),    Diff3D( 0, -1,  0),    Diff3D( 1, -1,  0), Diff3D( 2, -1,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0),    Diff3D( 2,  0,  0)    },    //Behind-SW
                {    Diff3D(-1, -2, -2),    Diff3D( 0, -2, -2),    Diff3D( 1, -2, -2),    Diff3D(-1, -1, -2),    Diff3D( 0, -1, -2),    Diff3D( 1, -1, -2),    Diff3D(-1,  0, -2),    Diff3D( 0,  0, -2),    Diff3D( 1, 0, -2),    Diff3D(-1, -2, -1),    Diff3D( 0, -2, -1),    Diff3D( 1, -2, -1),    Diff3D(-1, -1, -1),    Diff3D( 1, -1, -1),    Diff3D(-1,  0, -1),    Diff3D( 0,  0, -1),    Diff3D( 1,  0, -1),    Diff3D(-1, -2,  0),    Diff3D( 0, -2,  0),    Diff3D( 1, -2,  0),    Diff3D(-1, -1,  0),    Diff3D( 0, -1,  0), Diff3D( 1, -1,  0),    Diff3D(-1,  0,  0),    Diff3D( 0,  0,  0),    Diff3D( 1,  0,  0)    },    //Behind-S
                {    Diff3D(-2, -2, -2),    Diff3D(-1, -2, -2),    Diff3D( 0, -2, -2),    Diff3D(-2, -1, -2),    Diff3D(-1, -1, -2),    Diff3D( 0, -1, -2),    Diff3D(-2,  0, -2),    Diff3D(-1,  0, -2),    Diff3D( 0, 0, -2),    Diff3D(-2, -2, -1),    Diff3D(-1, -2, -1),    Diff3D( 0, -2, -1),    Diff3D(-2, -1, -1),    Diff3D( 0, -1, -1),    Diff3D(-2,  0, -1),    Diff3D(-1,  0, -1),    Diff3D( 0,  0, -1),    Diff3D(-2, -2,  0),    Diff3D(-1, -2,  0),    Diff3D( 0, -2,  0),    Diff3D(-2, -1,  0),    Diff3D(-1, -1,  0), Diff3D( 0, -1,  0),    Diff3D(-2,  0,  0),    Diff3D(-1,  0,  0), Diff3D( 0,  0,  0)    }        //Behind-SE
            };
            return d[fromCode][toCode];
        */
        return diff(toCode)-diff(fromCode);
    }

   /** Equivalent to relativeDiff(static_cast<Direction>(fromCode), static_cast<Direction>(toCode)).
       (note: there is no bounds checking on the code you pass.)
   */
    static Diff3D const relativeDiff(int fromCode, int toCode)
    {
        return relativeDiff(static_cast<Direction>(fromCode), static_cast<Direction>(toCode));
    }
    
    /**  X-component of diff() */
    static int dX(Direction code) { return diff(code)[0]; }
    /**  Y-component of diff() */
    static int dY(Direction code) { return diff(code)[1]; }
    /**  Z-component of diff() */
    static int dZ(Direction code) { return diff(code)[2]; }
    
    /**  X-component of diff() */
    static int dX(int code) { return diff(code)[0]; }
    /**  Y-component of diff() */
    static int dY(int code) { return diff(code)[1]; }
    /**  Z-component of diff() */
    static int dZ(int code) { return diff(code)[2]; }
    
    /** transform 6-neighborhood code into 26-neighborhood code.
     */
    static Direction code(Neighborhood3DSix::Direction d)
    {
        switch (d){
            case Neighborhood3DSix::InFront :
                    return InFront;
            case Neighborhood3DSix::North :
                    return North;
            case Neighborhood3DSix::West :
                    return West;
            case Neighborhood3DSix::East :
                    return East;
            case Neighborhood3DSix::South :
                    return South;
            case Neighborhood3DSix::Behind :
                    return Behind;
        }
        return Error;
    }

    /** transform Diff3D offset into corresponding direction code.
       The code <tt>Direction::Error</tt> will be returned if <tt>diff</tt>
       is not in the 3DTwentySix-Neighborhood.
    */
    static Direction code(Diff3D const & diff)
    {
        switch(diff[0]){
            case -1:
                switch(diff[1]){
                    case -1:
                        switch(diff[2]){
                            case -1: return InFrontNorthWest; // ( -1, -1, -1)
                            case  0: return NorthWest;        // ( -1, -1,  0)
                            case  1: return BehindNorthWest;  // ( -1, -1,  1)
                        }
                    case  0:
                        switch(diff[2]){
                            case -1: return InFrontWest;      // ( -1,  0, -1)
                            case  0: return West;             // ( -1,  0,  0)
                            case  1: return BehindWest;       // ( -1,  0,  1)
                        }
                    case  1:
                        switch(diff[2]){
                            case -1: return InFrontSouthWest; // ( -1,  1, -1)
                            case  0: return SouthWest;        // ( -1,  1,  0)
                            case  1: return BehindSouthWest;  // ( -1,  1,  1)
                        }
                }
            case  0:
                switch(diff[1]){
                    case -1:
                        switch(diff[2]){
                            case -1: return InFrontNorth;     // (  0,  0, -1)
                            case  0: return North;            // (  0, -1,  0)
                            case  1: return BehindNorth;      // (  0, -1,  1)
                        }
                    case  0:
                        switch(diff[2]){
                            case -1: return InFront;          // (  0,  0, -1)
                            case  1: return Behind;           // (  0,  0,  1)
                        }
                    case  1: 
                        switch(diff[2]){
                            case -1: return InFrontSouth;     // (  0,  1, -1)
                            case  0: return South;            // (  0,  1,  0)
                            case  1: return BehindSouth;      // (  0,  1,  1)
                        }
                }                    
            case  1:
                switch(diff[1]){
                    case -1:
                        switch(diff[2]){
                            case -1: return InFrontNorthEast;  // (  1, -1, -1)
                            case  0: return NorthEast;         // (  1, -1,  0)
                            case  1: return BehindNorthEast;   // (  1, -1,  1)
                        }
                    case  0:
                        switch(diff[2]){
                            case -1: return InFrontEast;       // (  1,  0, -1)
                            case  0: return East;              // (  1,  0,  0)
                            case  1: return BehindEast;        // (  1,  0,  1)
                        }
                    case  1:
                        switch(diff[2]){
                            case -1: return InFrontSouthEast;  // (  1,  1, -1)
                            case  0: return SouthEast;         // (  1,  1,  0)
                            case  1: return BehindSouthEast;   // (  1,  1,  1)
                        }
                }
        }
        return Error; // better safe than sorry
    }
  
    /** Check whether a code refers to a diagonal direction.
        Useful if you want to abstract the differences between 6- and 26-neighborhood.
        Always <tt>false</tt> for 6-neighborhood.
    */
    static bool isDiagonal(Direction dir) {
        Diff3D d = diff(dir);
        if (abs(d[0])+abs(d[1])+abs(d[2])==1)
            return false;
        else
            return true;
    }

    static Diff3D const & frontTopLeft()        { return diff(InFrontNorthWest); }    /**<  Offset to the front-top-left neighbor */
    static Diff3D const & frontTop()            { return diff(InFrontNorth); }        /**<  Offset to the front-top neighbor */
    static Diff3D const & frontTopRight()       { return diff(InFrontNorthEast); }    /**<  Offset to the front-top-right neighbor */
    static Diff3D const & frontLeft()           { return diff(InFrontWest); }         /**<  Offset to the front-left neighbor */
    static Diff3D const & front()               { return diff(InFront); }             /**<  Offset to the front neighbor */
    static Diff3D const & frontRight()          { return diff(InFrontEast); }         /**<  Offset to the front-right neighbor */
    static Diff3D const & frontBottomLeft()     { return diff(InFrontSouthWest); }    /**<  Offset to the front-bottom-left neighbor */
    static Diff3D const & frontBottom()         { return diff(InFrontSouth); }        /**<  Offset to the front-bottom neighbor */
    static Diff3D const & frontBottomRight()    { return diff(InFrontSouthEast); }    /**<  Offset to the front-bottom-right neighbor */
    
    static Diff3D const & topLeft()             { return diff(NorthWest); }           /**<  Offset to the top-left neighbor */
    static Diff3D const & top()                 { return diff(North); }               /**<  Offset to the top neighbor */
    static Diff3D const & topRight()            { return diff(NorthEast); }           /**<  Offset to the top-right neighbor */
    static Diff3D const & left()                { return diff(West); }                /**<  Offset to the left neighbor */
    static Diff3D const & right()               { return diff(East); }                /**<  Offset to the right neighbor */
    static Diff3D const & bottomLeft()          { return diff(SouthWest); }           /**<  Offset to the bottom-left neighbor */
    static Diff3D const & bottom()              { return diff(South); }               /**<  Offset to the bottom neighbor */
    static Diff3D const & bottomRight()         { return diff(SouthEast); }           /**<  Offset to the bottom-right neighbor */

    static Diff3D const & rearTopLeft()         { return diff(BehindNorthWest); }     /**<  Offset to the rear-top-left neighbor */
    static Diff3D const & rearTop()             { return diff(BehindNorth); }         /**<  Offset to the rear-top neighbor */
    static Diff3D const & rearTopRight()        { return diff(BehindNorthEast); }     /**<  Offset to the rear-top-right neighbor */
    static Diff3D const & rearLeft()            { return diff(BehindWest); }          /**<  Offset to the rear-left neighbor */
    static Diff3D const & rear()                { return diff(Behind); }              /**<  Offset to the rear neighbor */
    static Diff3D const & rearRight()           { return diff(BehindEast); }          /**<  Offset to the rear-right neighbor */
    static Diff3D const & rearBottomLeft()      { return diff(BehindSouthWest); }     /**<  Offset to the rear-bottom-left neighbor */
    static Diff3D const & rearBottom()          { return diff(BehindSouth); }         /**<  Offset to the rear-bottom neighbor */
    static Diff3D const & rearBottomRight()     { return diff(BehindSouthEast); }     /**<  Offset to the rear-bottom-right neighbor */

    //----- other namings

    static Diff3D const & infrontNorthWest()    { return diff(InFrontNorthWest); }    /**<  Offset to the infront-north-west neighbor */
    static Diff3D const & infrontNorth()        { return diff(InFrontNorth); }        /**<  Offset to the infront-north neighbor */
    static Diff3D const & infrontNorthEast()    { return diff(InFrontNorthEast); }    /**<  Offset to the infront-north-east neighbor */
    static Diff3D const & infrontWest()         { return diff(InFrontWest); }         /**<  Offset to the infront-west neighbor */
    static Diff3D const & infront()             { return diff(InFront); }             /**<  Offset to the infront neighbor */
    static Diff3D const & infrontEast()         { return diff(InFrontEast); }         /**<  Offset to the infront-east neighbor */
    static Diff3D const & infrontSouthWest()    { return diff(InFrontSouthWest); }    /**<  Offset to the infront-south-west neighbor */
    static Diff3D const & infrontSouth()        { return diff(InFrontSouth); }        /**<  Offset to the infront-south neighbor */
    static Diff3D const & infrontSouthEast()    { return diff(InFrontSouthEast); }    /**<  Offset to the infront-south-east neighbor */
    
    static Diff3D const & northWest()           { return diff(NorthWest); }            /**<  Offset to the north-west neighbor */
    static Diff3D const & north()               { return diff(North); }                /**<  Offset to the north neighbor */
    static Diff3D const & northEast()           { return diff(NorthEast); }            /**<  Offset to the north-east neighbor */
    static Diff3D const & west()                { return diff(West); }                 /**<  Offset to the west neighbor */
    static Diff3D const & east()                { return diff(East); }                 /**<  Offset to the right neighbor */
    static Diff3D const & southWest()           { return diff(SouthWest); }            /**<  Offset to the south-west neighbor */
    static Diff3D const & south()               { return diff(South); }                /**<  Offset to the south neighbor */
    static Diff3D const & southEast()           { return diff(SouthEast); }            /**<  Offset to the south-east neighbor */

    static Diff3D const & behindNorthWest()     { return diff(BehindNorthWest); }      /**<  Offset to the behind-north-west neighbor */
    static Diff3D const & behindNorth()         { return diff(BehindNorth); }          /**<  Offset to the behind-north neighbor */
    static Diff3D const & behindNorthEast()     { return diff(BehindNorthEast); }      /**<  Offset to the behind-north-east neighbor */
    static Diff3D const & behindEast()          { return diff(BehindWest); }           /**<  Offset to the behind-west neighbor */
    static Diff3D const & behind()              { return diff(Behind); }               /**<  Offset to the behind neighbor */
    static Diff3D const & behindWest()          { return diff(BehindEast); }           /**<  Offset to the behind-right neighbor */
    static Diff3D const & behindSouthWest()     { return diff(BehindSouthWest); }      /**<  Offset to the behind-south-west neighbor */
    static Diff3D const & behindSouth()         { return diff(BehindSouth); }          /**<  Offset to the behind-south neighbor */
    static Diff3D const & behindSouthEast()     { return diff(BehindSouthEast); }      /**<  Offset to the behind-south-east neighbor */
}; // class Neighborhood3D


/** Export NeighborCode3D::Direction into the scope of namespace Neighborhood3DSix.
 */
typedef NeighborCode3D::Direction Direction;

static const Direction InFrontNorthWest   = NeighborCode3D::InFrontNorthWest;     /**<  Export NeighborCode3D::InFrontNorthWest to namespace Neighborhood3DTwentySix */
static const Direction InFrontNorth       = NeighborCode3D::InFrontNorth;         /**<  Export NeighborCode3D::InFrontNorth to namespace Neighborhood3DTwentySix */
static const Direction InFrontNorthEast   = NeighborCode3D::InFrontNorthEast;     /**<  Export NeighborCode3D::InFrontNorthEast to namespace Neighborhood3DTwentySix */
static const Direction InFrontWest        = NeighborCode3D::InFrontWest;          /**<  Export NeighborCode3D::InFrontWest to namespace Neighborhood3DTwentySix */
static const Direction InFront            = NeighborCode3D::InFront;              /**<  Export NeighborCode3D::InFront to namespace Neighborhood3DTwentySix */
static const Direction InFrontEast        = NeighborCode3D::InFrontEast;          /**<  Export NeighborCode3D::InFrontEast to namespace Neighborhood3DTwentySix */
static const Direction InFrontSouthWest   = NeighborCode3D::InFrontSouthWest;     /**<  Export NeighborCode3D::InFrontSouthWest to namespace Neighborhood3DTwentySix */
static const Direction InFrontSouth       = NeighborCode3D::InFrontSouth;         /**<  Export NeighborCode3D::InFrontSouth to namespace Neighborhood3DTwentySix */
static const Direction InFrontSouthEast   = NeighborCode3D::InFrontSouthEast;     /**<  Export NeighborCode3D::InFrontSouthEast to namespace Neighborhood3DTwentySix */

static const Direction NorthWest          = NeighborCode3D::NorthWest;            /**<  Export NeighborCode3D::NorthWest to namespace Neighborhood3DTwentySix */
static const Direction North              = NeighborCode3D::North;                /**<  Export NeighborCode3D::North to namespace Neighborhood3DTwentySix */
static const Direction NorthEast          = NeighborCode3D::NorthEast;            /**<  Export NeighborCode3D::NorthEast to namespace Neighborhood3DTwentySix */
static const Direction West               = NeighborCode3D::West;                 /**<  Export NeighborCode3D::West to namespace Neighborhood3DTwentySix */
static const Direction East               = NeighborCode3D::East;                 /**<  Export NeighborCode3D::East to namespace Neighborhood3DTwentySix */
static const Direction SouthWest          = NeighborCode3D::SouthWest;            /**<  Export NeighborCode3D::SouthWest to namespace Neighborhood3DTwentySix */
static const Direction South              = NeighborCode3D::South;                /**<  Export NeighborCode3D::South to namespace Neighborhood3DTwentySix */
static const Direction SouthEast          = NeighborCode3D::SouthEast;            /**<  Export NeighborCode3D::SouthEast to namespace Neighborhood3DTwentySix */

static const Direction BehindNorthWest    = NeighborCode3D::BehindNorthWest;      /**<  Export NeighborCode3D::BehindNorthWest to namespace Neighborhood3DTwentySix */
static const Direction BehindNorth        = NeighborCode3D::BehindNorth;          /**<  Export NeighborCode3D::BehindNorth to namespace Neighborhood3DTwentySix */
static const Direction BehindNorthEast    = NeighborCode3D::BehindNorthEast;      /**<  Export NeighborCode3D::BehindNorthEast to namespace Neighborhood3DTwentySix */
static const Direction BehindWest         = NeighborCode3D::BehindWest;           /**<  Export NeighborCode3D::BehindWest to namespace Neighborhood3DTwentySix */
static const Direction Behind             = NeighborCode3D::Behind;               /**<  Export NeighborCode3D::Behind to namespace Neighborhood3DTwentySix */
static const Direction BehindEast         = NeighborCode3D::BehindEast;           /**<  Export NeighborCode3D::BehindEast to namespace Neighborhood3DTwentySix */
static const Direction BehindSouthWest    = NeighborCode3D::BehindSouthWest;      /**<  Export NeighborCode3D::BehindSouthWest to namespace Neighborhood3DTwentySix */
static const Direction BehindSouth        = NeighborCode3D::BehindSouth;          /**<  Export NeighborCode3D::BehindSouth to namespace Neighborhood3DTwentySix */
static const Direction BehindSouthEast    = NeighborCode3D::BehindSouthEast;      /**<  Export NeighborCode3D::BehindSouthEast to namespace Neighborhood3DTwentySix */

static const Direction DirectionCount     = NeighborCode3D::DirectionCount;       /**<  Export NeighborCode3D::DirectionCount to namespace Neighborhood3DTwentySix */

}//namespace Neighborhood3DTwentySix
    
/** Export \ref vigraext::Neighborhood3DTwentySix::NeighborCode3D into the scope of namespace vigraext.
 */
typedef Neighborhood3DTwentySix::NeighborCode3D NeighborCode3DTwentySix;



/********************************************************/
/*                                                      */
/*              NeighborOffsetTraverser                 */
/*                                                      */
/********************************************************/

/** \brief Traverser that walks around a given location.

    The template parameter defines the kind of neighborhood used, e.g.

    \code
    NeighborOffsetTraverser<NeighborCode3DSix> traverser_3dsix;
    \endcode

    Since this traverser doesn't know about the voxels in any particular volume,
    you usually doesn't use it directly but rather as a base class or helper for
    neighborhood traversers refering to a particular volume (e.g. NeighborhoodTraverser)
*/

template<class NEIGHBORCODE3D>
class NeighborOffsetTraverser
: public NEIGHBORCODE3D
{
public:
    typedef NEIGHBORCODE3D NeighborCode;

    /** return type of direction()
    */
    typedef typename NEIGHBORCODE3D::Direction Direction;

    /** the traverser's value type
    */
    typedef Diff3D value_type;

    /** the traverser's reference type (return type of <TT>*trav</TT>)
      */
    typedef Diff3D const & reference;

    /** the traverser's pointer type (return type of <TT>operator-></TT>)
     */
    typedef Diff3D const * pointer;

    /** the traverser's difference type (argument type of <TT>trav[diff]</TT>)
     */
    typedef int difference_type;

protected:
    Direction direction_;

public:
    /** Create traverser refering to the given direction.
     */
    NeighborOffsetTraverser(Direction dir = NEIGHBORCODE3D::CausalFirst)
        : direction_(dir)
    {
    }

    /** pre-increment */
    NeighborOffsetTraverser & operator++()
    {
        direction_ = static_cast<Direction>((direction_+1) % NEIGHBORCODE3D::DirectionCount);
        return *this;
    }

    /** pre-decrement */
    NeighborOffsetTraverser & operator--()
    {
        direction_ = static_cast<Direction>((direction_ + NEIGHBORCODE3D::DirectionCount-1) % NEIGHBORCODE3D::DirectionCount);
        return *this;
    }

    /** post-increment */
    NeighborOffsetTraverser operator++(int)
    {
        NeighborOffsetTraverser ret(*this);
        operator++();
        return ret;
    }

    /** post-decrement */
    NeighborOffsetTraverser operator--(int)
    {
        NeighborOffsetTraverser ret(*this);
        operator--();
        return ret;
    }

    /** add-assignment */
    NeighborOffsetTraverser & operator+=(difference_type d)
    {
        direction_ = static_cast<Direction>((direction_ + d) % NEIGHBORCODE3D::DirectionCount);
        if(direction_ < 0)
        direction_ = static_cast<Direction>(direction_ + NEIGHBORCODE3D::DirectionCount);
        return *this;
    }

    /** subtract-assignment */
    NeighborOffsetTraverser & operator-=(difference_type d)
    {
        direction_ = static_cast<Direction>((direction_ - d) % NEIGHBORCODE3D::DirectionCount);
        if(direction_ < 0)
        direction_ = static_cast<Direction>(direction_ + NEIGHBORCODE3D::DirectionCount);
        return *this;
    }

    /** addition */
    NeighborOffsetTraverser operator+(difference_type d) const
    {
        return NeighborOffsetTraverser(*this) += d;
    }

    /** subtraction */
    NeighborOffsetTraverser operator-(difference_type d) const
    {
        return NeighborOffsetTraverser(*this) -= d;
    }

    /** Move to the opposite direction of the current direction.
     */
    NeighborOffsetTraverser & turnRound()
    {
                direction_ = opposite();
        return *this;
    }

    /** Move to the given direction.
     */
    NeighborOffsetTraverser & turnTo(Direction d)
    {
        direction_ = d;
        return *this;
    }

    /** equality */
    bool operator==(NeighborOffsetTraverser const & o) const
    {
        return direction_ == o.direction_;
    }

    /** unequality */
    bool operator!=(NeighborOffsetTraverser const & o) const
    {
        return direction_ != o.direction_;
    }

    /** subtraction */
    difference_type operator-(NeighborOffsetTraverser const & o) const
    {
        return direction_ - o.direction_;
    }

    /** dereference */
    reference operator*() const
    {
        return diff();
    }

    /** index */
    reference operator[](difference_type d) const
    {
        return NEIGHBORCODE3D::diff(direction(d));
    }

    /** member access */
    pointer operator->() const
    {
        return &diff();
    }

    /** Get Diff3D offset from center to current neighbor.
     */
    Diff3D const & diff() const
    {
        return NEIGHBORCODE3D::diff(direction_);
    }

    /** Get Diff3D offset to given direction.
     */
    static Diff3D const & diff(Direction dir)
    {
        return NEIGHBORCODE3D::diff(dir);
    }

    /** Get relative distance (Diff3D) from current neighbor to neighbor
        at given offset.
     */
    Diff3D const relativeDiff(difference_type offset) const
    {
        Direction toDir = static_cast<Direction>((direction_ + offset) % NEIGHBORCODE3D::DirectionCount);
        if(toDir < 0)
        toDir = static_cast<Direction>(toDir + NEIGHBORCODE3D::DirectionCount);
        return NEIGHBORCODE3D::relativeDiff(direction_, toDir);
    }

    /** X-component of diff()  */
    int dX() const
    {
        return NEIGHBORCODE3D::dX(direction_);
    }

    /** Y-component of diff() */
    int dY() const
    {
        return NEIGHBORCODE3D::dY(direction_);
    }
        
        /** Z-component of diff() */
    int dZ() const
    {
        return NEIGHBORCODE3D::dZ(direction_);
    }

    /** Check whether current direction is a diagonal one.
     */
    bool isDiagonal() const
    {
        return NEIGHBORCODE3D::isDiagonal(direction_);
    }

    /** Get current direction.
     */
    Direction direction() const
    {
        return direction_;
    }

    /** Get current direction bit.
     */
    unsigned int directionBit() const
    {
        return NEIGHBORCODE3D::directionBit(direction_);
    }

    /** Get opposite of current direction.
     */
    Direction opposite() const
    {
        return static_cast<Direction>((NEIGHBORCODE3D::OppositeDirPrefix*direction_ + NEIGHBORCODE3D::OppositeOffset) % NEIGHBORCODE3D::DirectionCount);
    }

    /** Get opposite bit of current direction.
     */
    unsigned int oppositeDirectionBit() const
    {
        return NEIGHBORCODE3D::directionBit(opposite());
    }

    /** Get direction code at offset of current direction.
     */
    Direction direction(difference_type offset) const
    {
        int result = (direction_ + offset) % NEIGHBORCODE3D::DirectionCount;
        if(result < 0)
        result += NEIGHBORCODE3D::DirectionCount;
        return static_cast<Direction>(result);
    }
};//class NeighborOffsetTraverser

/** Specialization of NeighborOffsetTraverser for Neighborhood3DSix.
*/
typedef NeighborOffsetTraverser<NeighborCode3DSix> NeighborOffsetTraverser3DSix;

/** Specialization of NeighborOffsetTraverser for Neighborhood3DTwentySix.
*/
typedef NeighborOffsetTraverser<NeighborCode3DTwentySix> NeighborOffsetTraverser3DTwentySix;


/********************************************************/
/*                                                      */
/*                NeighborhoodTraverser                 */
/*                                                      */
/********************************************************/

/** \brief Traverser that walks around a given location in a given volume.

    The template parameters define the kind of neighborhood used and the underlying
    volume. The access functions return the value of the current neighbor voxel. 
    Use <tt>center()</tt> to access the center voxel of the neighborhood. 
    The center can be changed by calling <tt>moveCenterToNeighbor()</tt> 
    or <tt>swapCenterNeighbor()</tt>. Note that this traverser cannot
    when the center is at the volume border. You must then use 
    \ref vigra::RestrictedNeighborhoodTraverser
*/
template <class VOLUMEITERATOR, class NEIGHBORCODE3D>
class NeighborhoodTraverser : private VOLUMEITERATOR
{
    typedef NeighborOffsetTraverser<NEIGHBORCODE3D> NEIGHBOROFFSETTRAVERSER;

public:
    /** type of the underlying volume iterator
     */
    typedef VOLUMEITERATOR base_type;

    /** type of the used neighbor code
     */
    typedef NEIGHBORCODE3D NeighborCode;
 
    /** the traverser's value type
     */
    typedef typename VOLUMEITERATOR::value_type value_type;
    
    /** type of the direction code
     */
    typedef typename NEIGHBORCODE3D::Direction Direction;
    
    /** the traverser's reference type (return type of <TT>*trav</TT>)
     */
    typedef typename VOLUMEITERATOR::reference reference;

    /** the traverser's pointer type (return type of <TT>operator-></TT>)
     */
    typedef typename VOLUMEITERATOR::pointer pointer;

    /** the traverser's difference type (argument type of <TT>trav[diff]</TT>)
     */
    typedef typename NEIGHBOROFFSETTRAVERSER::difference_type difference_type;

    /** Construct traverser with given <tt>center</tt> pixel, pointing to the neighbor
         at the given direction <tt>d</tt>.
     */
    NeighborhoodTraverser(VOLUMEITERATOR const & center = VOLUMEITERATOR(),
                          Direction d = NEIGHBOROFFSETTRAVERSER::CausalFirst)
        : VOLUMEITERATOR(center), neighborCode_(d)
    {
        VOLUMEITERATOR::operator+=(neighborCode_.diff());
    }

    /** pre-increment */
    NeighborhoodTraverser & operator++()
    {
        return operator+=(1);
    }

    /** pre-decrement */
    NeighborhoodTraverser operator++(int)
    {
        NeighborhoodTraverser ret(*this);
        operator++();
        return ret;
    }

    /** post-increment */
    NeighborhoodTraverser & operator--()
    {
        return operator+=(-1);
    }

    /** post-decrement */
    NeighborhoodTraverser operator--(int)
    {
        NeighborhoodTraverser ret(*this);
        operator--();
        return ret;
    }

    /** add-assignment */
    NeighborhoodTraverser & operator+=(difference_type d)
    {
        VOLUMEITERATOR::operator+=(neighborCode_.relativeDiff(d));
        neighborCode_+= d;
        return *this;
    }

    /** subtract-assignment */
    NeighborhoodTraverser & operator-=(difference_type d)
    {
        return operator+=(-d);
    }

    /** addition */
    NeighborhoodTraverser operator+(difference_type d) const
    {
        NeighborhoodTraverser result(*this);
        result+= d;
        return result;
    }

    /** subtraction */
    NeighborhoodTraverser operator-(difference_type d) const
    {
        NeighborhoodTraverser result(*this);
        result-= d;
        return result;
    }

    /** Move to the opposite direction of the current direction.
     */
    NeighborhoodTraverser & turnRound()
    {
        Direction oldDirection = neighborCode_.direction();
        neighborCode_.turnRound();
        VOLUMEITERATOR::operator+=(NeighborCode::relativeDiff
                                  (oldDirection, neighborCode_.direction()));
        return *this;
    }

    /** Move to the given direction.
     */
    NeighborhoodTraverser & turnTo(Direction d)
    {
        Direction oldDirection = neighborCode_.direction();
        neighborCode_.turnTo(d);
        VOLUMEITERATOR::operator+=(NeighborCode::relativeDiff
                                  (oldDirection, neighborCode_.direction()));
        return *this;
    }

    /** Move the center in the current direction.
        The current neighbor becomes the new center, the direction does not change.
     */
    NeighborhoodTraverser & moveCenterToNeighbor()
    {
        VOLUMEITERATOR::operator+=(neighborCode_.diff());
        return *this;
    }

    /** Exchange the center with the current neighbor.
         Equivalent to <tt>circ.moveCenterToNeighbor().turnRound()</tt>
         (but shorter and more efficient).
     */
    NeighborhoodTraverser & swapCenterNeighbor()
    {
        neighborCode_.turnRound();
        VOLUMEITERATOR::operator+=(neighborCode_.diff());
        return *this;
    }

    /** equality */
    bool operator==(NeighborhoodTraverser const & rhs) const
    {
        return neighborCode_ == rhs.neighborCode_ &&
        VOLUMEITERATOR::operator==(rhs);
    }

    /** inequality */
    bool operator!=(NeighborhoodTraverser const & rhs) const
    {
        return neighborCode_ != rhs.neighborCode_ ||
        VOLUMEITERATOR::operator!=(rhs);
    }

    /** subtraction */
    difference_type operator-(NeighborhoodTraverser const & rhs) const
    {
        return neighborCode_ - rhs.neighborCode_;
    }

    /** dereference */
    reference operator*() const
    {
        return VOLUMEITERATOR::operator*();
    }

    /** index */
    reference operator[](difference_type d) const
    {
        return VOLUMEITERATOR::operator[](neighborCode_.relativeDiff(d));
    }

    /** member access */
    pointer operator->() const
    {
        return VOLUMEITERATOR::operator->();
    }

    /** Get the base iterator for the current neighbor. */
    base_type const & base() const
    {
        return *this;
    }

    /** Get the base iterator for the center of the traverser. */
    base_type center() const
    {
        return (base_type)*this - neighborCode_.diff();
    }

    /** Get the current direction. */
    Direction direction() const
    {
        return neighborCode_.direction();
    }

    /** Get the current direction bit. */
    unsigned int directionBit() const
    {
        return neighborCode_.directionBit();
    }

    /** Get opposite bit of current direction.
     */
    unsigned int oppositeDirectionBit() const
    {
        return neighborCode_.oppositeDirectionBit();
    }


    /** Get the difference vector (Diff3D) from the center to the current neighbor. */
    Diff3D const & diff() const
    {
        return neighborCode_.diff();
    }

    /** Is the current neighbor a diagonal neighbor? */
    bool isDiagonal() const
    {
        return neighborCode_.isDiagonal();
    }

private:
    NEIGHBOROFFSETTRAVERSER neighborCode_;
};//class NeighborhoodTraverser

/********************************************************/
/*                                                      */
/*            RestrictedNeighborhoodTraverser          */
/*                                                      */
/********************************************************/

/** \brief Traverser that walks around a given location in a given volume,
    unsing a restricted neighborhood.

    This traverser behaves essentially like \ref vigra::NeighborhoodTraverser,
    but can also be used near the volume border, where some of the neighbor voxel
    would be outside the volume und must not be accessed. 
    The template parameters define the kind of neighborhood used (only 3DSix)
    and the underlying volume, whereas the required neighborhood restriction is 
    given by the last constructor argument. This below for typical usage.

    The access functions return the value of the current neighbor pixel. Use <tt>center()</tt> to
    access the center pixel of the neighborhood.
*/
template <class VOLUMEITERATOR, class NEIGHBORCODE3D>
class RestrictedNeighborhoodTraverser 
: private NeighborhoodTraverser<VOLUMEITERATOR, NEIGHBORCODE3D>
{
    typedef NeighborhoodTraverser<VOLUMEITERATOR, NEIGHBORCODE3D> BaseType;

    public:
    /** type of the underlying volume iterator
     */
    typedef VOLUMEITERATOR base_type;

    /** type of the used neighbor code
     */
    typedef NEIGHBORCODE3D NeighborCode;

    /** the traverser's value type
     */
    typedef typename BaseType::value_type value_type;

    /** type of the direction code
     */
    typedef typename BaseType::Direction Direction;

    /** the traverser's reference type (return type of <TT>*trav</TT>)
     */
    typedef typename BaseType::reference reference;

    /** the traverser's pointer type (return type of <TT>operator-></TT>)
     */
    typedef typename BaseType::pointer pointer;

    /** the traverser's difference type (argument type of <TT>trav[diff]</TT>)
     */
    typedef typename BaseType::difference_type difference_type;

    /** Construct traverser with given <tt>center</tt> voxel, using the restricted
         neighborhood given by \a atVolumeBorder.
     */
    RestrictedNeighborhoodTraverser(VOLUMEITERATOR const & center = VOLUMEITERATOR(),
                                     AtVolumeBorder atBorder = NotAtBorder)
        : BaseType(center, NEIGHBORCODE3D::nearBorderDirections(atBorder, 0)),
          whichBorder_(atBorder),
          count_(NEIGHBORCODE3D::nearBorderDirectionCount(atBorder)),
          current_(0)
    {}

    /** pre-increment */
    RestrictedNeighborhoodTraverser & operator++()
    {
        return operator+=(1);
    }

    /** post-increment */
    RestrictedNeighborhoodTraverser operator++(int)
    {
        RestrictedNeighborhoodTraverser ret(*this);
        operator++();
        return ret;
    }

    /** pre-decrement */
    RestrictedNeighborhoodTraverser & operator--()
    {
        return operator+=(-1);
    }

    /** post-decrement */
    RestrictedNeighborhoodTraverser operator--(int)
    {
        RestrictedNeighborhoodTraverser ret(*this);
        operator--();
        return ret;
    }

    /** add-assignment */
    RestrictedNeighborhoodTraverser & operator+=(difference_type d)
    {
        current_ = static_cast<Direction>((current_ + count_ + d) % count_);
        BaseType::turnTo(NEIGHBORCODE3D::nearBorderDirections(whichBorder_, current_));
        return *this;
    }

    /** subtract-assignment */
    RestrictedNeighborhoodTraverser & operator-=(difference_type d)
    {
        return operator+=(-d);
    }

    /** addition */
    RestrictedNeighborhoodTraverser operator+(difference_type d) const
    {
        RestrictedNeighborhoodTraverser result(*this);
        result+= d;
        return result;
    }

    /** subtraction */
    RestrictedNeighborhoodTraverser operator-(difference_type d) const
    {
        RestrictedNeighborhoodTraverser result(*this);
        result-= d;
        return result;
    }

    /** equality */
    bool operator==(RestrictedNeighborhoodTraverser const & rhs) const
    {
        return current_ == rhs.current_;
    }

    /** inequality */
    bool operator!=(RestrictedNeighborhoodTraverser const & rhs) const
    {
        return current_ != rhs.current_;
    }

    /** subtraction */
    difference_type operator-(RestrictedNeighborhoodTraverser const & rhs) const
    {
        return (current_ - rhs.current_) % count_;
    }

    /** dereference */
    reference operator*() const
    {
        return BaseType::operator*();
    }

    /** member access */
    pointer operator->() const
    {
        return BaseType::operator->();
    }

    /** Get the base iterator for the current neighbor. */
    base_type const & base() const
    {
        return BaseType::base();
    }

    /** Get the base iterator for the center of the traverser. */
    base_type center() const
    {
        return BaseType::center();
    }

    /** Get the current direction. */
    Direction direction() const
    {
        return BaseType::direction();
    }

    /** Get the current direction bit. */
    unsigned int directionBit() const
    {
        return BaseType::directionBit();
    }

    /** Get the difference vector (Diff3D) from the center to the current neighbor. */
    Diff3D const & diff() const
    {
        return BaseType::diff();
    }

    /** Is the current neighbor a diagonal neighbor? */
    bool isDiagonal() const
    {
        return BaseType::isDiagonal();
    }

    private:
    AtVolumeBorder whichBorder_;
    signed char count_, current_;
}; //class RestrictedNeighborhoodTraverser

} // namespace vigra

#endif /* VIGRA_VOXELNEIGHBORHOOD_HXX */
