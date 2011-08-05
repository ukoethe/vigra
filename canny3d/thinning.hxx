#ifndef THINNING_HXX
#define THINNING_HXX

#include <cmath>
#include <algorithm>
#include <vector>
#include <queue>

#include "vigra/multi_array.hxx"
#include "vigra/multi_impex.hxx"
#include "vigra/multi_convolution.hxx"
#include "vigra/multi_pointoperators.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/mathutil.hxx"
#include "vigra/transformimage.hxx"

#include "thinningMasks.hxx"

namespace vigra {

  /* --------- structs --------------*/

  namespace detail {
    
    template<typename Volume>
    struct SimplePoint3dTester
    {
      typedef vigra::UInt32 BitMask;
      typedef BitMask BitMasks [2];
      typedef BitMasks MaskSet [totalNumbOfMasks];
      typedef typename Volume::difference_type Diff3d;
      
      const Volume & vol;
      unsigned int numbOfMasks;
      MaskSet masks;

      SimplePoint3dTester(const Volume & vol,
                          bool fvols, bool csurfs, bool ccurvs,
                          bool osurfs, bool ocurvs,
                          bool point) 
        : vol(vol)
      { init(fvols, csurfs, ccurvs, osurfs, ocurvs, point); }

      SimplePoint3dTester(const Volume & vol) : vol(vol)
      { init(true, true, true, true, true, true); }


      void init(bool fvols, bool csurfs, bool ccurvs,
                bool osurfs, bool ocurvs,
                bool point)
        
      {
//         std::cerr << fvols << csurfs << ccurvs
//                   << osurfs << ocurvs << point << std::endl;

        //std::cerr << "totalNumbOfMasks" << totalNumbOfMasks << std::endl;
        unsigned int i = 0;
        unsigned int s = 0;
        
        if(fvols) 
          {
            s = 1; i = 1;
            masks[0][0] = erosionMask6ctF;
            masks[0][1] = erosionMask6ctB;
          }
        if(csurfs) 
          {
            s += sizeof(closedSurfacesMasks)/sizeof(BitMasks);
            for(unsigned int j = 0; i < s; ++i,++j) {
              masks[i][0] = closedSurfacesMasks[j][0];
              masks[i][1] = closedSurfacesMasks[j][1];
            }
          }
        if(ccurvs) 
          {
            s += sizeof(closedCurvesMasks)/sizeof(BitMasks);
            for(unsigned int j = 0; i < s; ++i,++j) {
              masks[i][0] = closedCurvesMasks[j][0];
              masks[i][1] = closedCurvesMasks[j][1];
            }
          }
        if(osurfs) 
          {
            s += sizeof(openSurfacesMasks)/sizeof(BitMasks);
            for(unsigned int j = 0; i < s; ++i,++j) {
              masks[i][0] = openSurfacesMasks[j][0];
              masks[i][1] = openSurfacesMasks[j][1];
            }
          }
        if(ocurvs) 
          {
            s += sizeof(openCurvesMasks)/sizeof(BitMasks);
            for(unsigned int j = 0; i < s; ++i,++j) {
              masks[i][0] = openCurvesMasks[j][0];
              masks[i][1] = openCurvesMasks[j][1];
            }
          }
        if(point) 
          {
            masks[s][0] = singlePointMaskF;
            masks[s][1] = singlePointMaskB;
            s += 1;
          }
        numbOfMasks = s;
        //        std::cerr << "numbOfMasks" << numbOfMasks << std::endl;
      }

      bool operator()(const Diff3d & x) const
      {
        
//         std::cerr << x << std::endl; 
//         std::cerr << int(vol[x].mark) << std::endl; 
        using namespace vigra;
        BitMask roi =
          ((((BitMask)(vol[x+Diff3d(-1,-1,-1)]>0))<<0 )|
           (((BitMask)(vol[x+Diff3d(-1,-1, 0)]>0))<<1 )|
           (((BitMask)(vol[x+Diff3d(-1,-1, 1)]>0))<<2 )|

           (((BitMask)(vol[x+Diff3d(-1, 0,-1)]>0))<<3 )|
           (((BitMask)(vol[x+Diff3d(-1, 0, 0)]>0))<<4 )|
           (((BitMask)(vol[x+Diff3d(-1, 0, 1)]>0))<<5 )|

           (((BitMask)(vol[x+Diff3d(-1, 1,-1)]>0))<<6 )|
           (((BitMask)(vol[x+Diff3d(-1, 1, 0)]>0))<<7 )|
           (((BitMask)(vol[x+Diff3d(-1, 1, 1)]>0))<<8 )|

     
           (((BitMask)(vol[x+Diff3d( 0,-1,-1)]>0))<<9 )|
           (((BitMask)(vol[x+Diff3d( 0,-1, 0)]>0))<<10)|
           (((BitMask)(vol[x+Diff3d( 0,-1, 1)]>0))<<11)|

           (((BitMask)(vol[x+Diff3d( 0, 0,-1)]>0))<<12)|
           (((BitMask)(vol[x+Diff3d( 0, 0, 0)]>0))<<13)|
           (((BitMask)(vol[x+Diff3d( 0, 0, 1)]>0))<<14)|

           (((BitMask)(vol[x+Diff3d( 0, 1,-1)]>0))<<15)|
           (((BitMask)(vol[x+Diff3d( 0, 1, 0)]>0))<<16)|
           (((BitMask)(vol[x+Diff3d( 0, 1, 1)]>0))<<17)|

     
           (((BitMask)(vol[x+Diff3d( 1,-1,-1)]>0))<<18)|
           (((BitMask)(vol[x+Diff3d( 1,-1, 0)]>0))<<19)|
           (((BitMask)(vol[x+Diff3d( 1,-1, 1)]>0))<<20)|

           (((BitMask)(vol[x+Diff3d( 1, 0,-1)]>0))<<21)|
           (((BitMask)(vol[x+Diff3d( 1, 0, 0)]>0))<<22)|
           (((BitMask)(vol[x+Diff3d( 1, 0, 1)]>0))<<23)|

           (((BitMask)(vol[x+Diff3d( 1, 1,-1)]>0))<<24)|
           (((BitMask)(vol[x+Diff3d( 1, 1, 0)]>0))<<25)|
           (((BitMask)(vol[x+Diff3d( 1, 1, 1)]>0))<<26)
           );
        
//         int c = 0;
//         for(int z = -1; z <= 1; z++)
//           for(int y = -1; y <= 1; y++)
//             for(int x_ = -1; x_ <= 1; x_++)
//               c+= (BitMask)vol[x+Diff3d( x_,y,z)].mark;
//         std::cerr << c << std::endl;
        

//         std::cerr << roi << std::endl;
//         std::cerr << masks[0][0] << std::endl;
//         std::cerr << (masks[0][0]& roi) << std::endl;
        
//         BitMask value = roi;
//         const int SHIFT = 26;
//         const unsigned MASK = 1 << SHIFT;
//         for ( int i = 1; i <= SHIFT + 1; i++ ) 
//           {
//             std::cerr << ( value & MASK ? '1' : '0' );
//             value <<= 1;
            
//             if ( i % 9 == 0 )
//               std::cerr << ' ';
//           }
        
//         std::cerr << std::endl;
//         std::cerr << std::endl;

        for(unsigned int mc = 0; mc < numbOfMasks;++mc)
          if(((masks[mc][0]& roi)==masks[mc][0])
             &&
             ((masks[mc][1]&~roi)==masks[mc][1])
             )
            return true;
        return false;
      }
    };
  }

  template<typename POSITIONTYPE, typename VALUETYPE>
  struct Surfel
  {
    typedef POSITIONTYPE PositionType;
    typedef VALUETYPE ValueType;

    PositionType p;
    ValueType strength;
 
    Surfel()
      : p(vigra::NumericTraits<PositionType>::zero()), 
        strength(vigra::NumericTraits<ValueType>::zero())
    {}
  
    Surfel(const PositionType & p, ValueType is)
      : p(p), 
        strength(is)
    {}
  
    Surfel(const PositionType & p)
      : p(p), 
        strength(vigra::NumericTraits<ValueType>::zero())
    {}

    template<typename SURFEL>
    bool operator<(const SURFEL & o) const
    {
      return strength < o.strength; 
    }
    
    template<typename SURFEL>
    bool operator>(const SURFEL  & o) const
    {
      return strength > o.strength; 
    }

    template<typename SURFEL>
    bool operator==(SURFEL const& x) const 
    { return p == x.p && strength == x.strength;}

    template<typename SURFEL>
    bool operator!=(SURFEL const& x) const 
    { return p != x.p || strength != x.strength; }
  };

  template<typename VALUETYPE>
  struct SimpleVoxel
  {
    typedef VALUETYPE ValueType;
    typedef float SquaredNormType;
    ValueType strength;
    unsigned char mark;
 
    SimpleVoxel()
      : strength(vigra::NumericTraits<ValueType>::zero()),
        mark(0)
    {}
  
    SimpleVoxel(ValueType s, unsigned char m)
      : strength(s), 
        mark(m)
    {}

    void operator=(unsigned char m)
    {
      mark = m; 
    }

    bool operator==(unsigned char m) const
    {
      return mark == m; 
    }

    bool operator!=(unsigned char m) const
    {
      return mark != m; 
    }
  
    bool operator<(unsigned char m) const
    {
      return mark < m; 
    }
    
    bool operator>(unsigned char m) const
    {
      return mark > m; 
    }
  };

  template <class Iterable, class BackInsertable, typename ShapeType>
  void surfelThinning(const Iterable & surfels,  BackInsertable & ret, 
                      ShapeType shape, 
                      bool fvols, bool csurfs, bool ccurvs,
                      bool osurfs, bool ocurvs,
                      bool point)
  {
    typedef typename Iterable::value_type Surfel;
    typedef typename Surfel::ValueType ValueType;
    typedef typename Surfel::PositionType PositionType;
    typedef SimpleVoxel<ValueType> Voxel;
    typedef MultiArray<3, Voxel> Volume;
    typedef typename Volume::difference_type Diff3D;
  
    // 18 if skeletonizing surfaces, 26 if skeletonizing curves as well
    const unsigned char neighbourhood = 26; 

    unsigned int scount = surfels.size(); 

//     std::cerr << "size: " << surfels.size() << std::endl;
//     std::cerr << "shape: " << shape << std::endl;

    MultiArray<3, Voxel> surfaceVolume(shape,Voxel(0,0));
    initMultiArray(destMultiArrayRange(surfaceVolume), Voxel(0,0));//FIXME: Remove

    detail::SimplePoint3dTester<MultiArray<3, Voxel> > 
      isSimplePoint(surfaceVolume,
                    fvols, csurfs, ccurvs,
                    osurfs, ocurvs,
                    point);
    
//     MultiArray<3, unsigned char> surfaceVolume(shape);
//     initMultiArray(destMultiArrayRange(surfaceVolume), 0);//FIXME: Remove


    std::priority_queue<Surfel, 
      std::vector<Surfel>, 
      std::greater<Surfel> > 
      pqueue;
    
    for (unsigned int i = 0; i < scount; ++i) {
      //surfaceVolume[surfels[i].p] = Voxel(surfels[i].strength,1);
      surfaceVolume[surfels[i].p] = 1;
      //std::cerr << (surfels[i].p -Diff3D(1,1,1)) << std::endl;
    }

//     for (unsigned int i = 0; i < scount; ++i) {
//       Diff3D x = surfels[i].p;
      
//       int c = 0;
//       for(int z = -1; z <= 1; z++)
//         for(int y = -1; y <= 1; y++)
//           for(int x1 = -1; x1 <= 1; x1++)
//             c+= (BitMask)surfaceVolume[x + Diff3D(x1,y,z)];
//       std::cerr << "neigs: " << c << std::endl;
//     }

    for (unsigned int i = 0; i < scount; ++i) {
      Diff3D x = surfels[i].p;

      if (!isSimplePoint(x)) //FIXME: change
        {
          pqueue.push(surfels[i]);
          surfaceVolume[x] = 2;
        }
    }
  
    static const Diff3D dist[] = {
      // 6 neighbourhood
      Diff3D(-1,0,0), Diff3D(1,0,0), Diff3D(0,-1,0), Diff3D(0,1,0), Diff3D(0,0,-1), Diff3D(0,0,1),
      // 18 neighbourhood
      Diff3D(0,-1,-1), Diff3D(0,-1,1), Diff3D(0,1,-1), Diff3D(0,1,1), Diff3D(-1,0,-1), Diff3D(-1,0,1),
      Diff3D(1,0,-1), Diff3D(1,0,1), Diff3D(-1,-1,0), Diff3D(-1,1,0), Diff3D(1,-1,0), Diff3D(1,1,0),
      // 26 neighbourhood
      Diff3D(-1,-1,-1), Diff3D(-1,-1,1), Diff3D(-1,1,-1), Diff3D(-1,1,1), Diff3D(1,-1,-1), Diff3D(1,-1,1), Diff3D(1,1,-1), Diff3D(1,1,1)
    };
                                
    while(pqueue.size())
     {
        Diff3D x = pqueue.top().p; pqueue.pop(); 

        if(isSimplePoint(x))
          {
            surfaceVolume[x] = 1;
            continue;
          }

        surfaceVolume[x] = 0; // delete simple point
        
         for(unsigned int i = 0; i < neighbourhood; ++i)
           {
             Diff3D neigh = x + dist[i];
             Voxel voxel = surfaceVolume[neigh];
             if(voxel == 1)
              {
                if(!isSimplePoint(neigh))
                  {
                    pqueue.push(Surfel(neigh, voxel.strength));
                    surfaceVolume[neigh] = 2;
                  }
              }
          }
     }
    
    Diff3D z(1),y,x;
    z[0]=z[1]=z[2]=1;
    for(; z[2] < shape[2]-1; ++z[2])
      for(y=z; y[1] < shape[1]-1; ++y[1])
        for(x=y; x[0] < shape[0]-1; ++x[0])
          if(surfaceVolume[x] > 0)
          //if(surfaceVolume[x] == 1)
            ret.push_back(Surfel(x,surfaceVolume[x].strength));
            //ret.push_back(Surfel(x,1.0));
  }
  
  template<typename SrcIterator, typename SrcShape, typename SrcAccessor, 
           typename BackInsertable, typename GradValue>
  void cannySurfelList(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                       BackInsertable & surfels, double sigma, GradValue treshold)
  {
    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote ScalarType;
    typedef vigra::TinyVector<ScalarType,3> GradientValueType;
    typedef typename AccessorTraits<GradientValueType>::default_accessor GradientAccessor;

    MultiArray<3, GradientValueType> grad(shape);
    MultiArray<3, ScalarType> magn(shape);
  
    gaussianGradientMultiArray(siter, shape, src, grad.traverser_begin(), GradientAccessor(), sigma);
    transformMultiArray(srcMultiArrayRange(grad), destMultiArray(magn), 
                        VectorNormFunctor<GradientValueType>());
  
    cannyFindSurfels(grad.traverser_begin(), shape, GradientAccessor(), 
                     magn.traverser_begin(),
                     typename vigra::AccessorTraits<ScalarType>::default_accessor(), 
                     surfels, treshold);
  }

  template<typename SrcIterator, typename SrcShape, typename SrcAccessor, 
           typename BackInsertable, typename GradValue>
  void cannySurfelList(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
                       BackInsertable & surfels, double sigma, GradValue threshold) 
  {
    cannySurfelList(src.first, src.second, src.third, 
                    surfels, sigma, threshold);
  }
} // end namespace vigra
#endif // THINNING_HXX
