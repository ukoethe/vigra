#ifndef BFGS_HXX
#define BFGS_HXX

#include "vigra/tinyvector.hxx"
#include "vigra/matrix.hxx"
#include <vigra/linear_algebra.hxx>
#include "vigra/mathutil.hxx"
#include "tools.hxx"
#include "linesearch.hxx"

namespace vigra {

  template<typename FUNCTOR>
  struct InvertFunctor
  {
    typedef FUNCTOR Functor;
    typedef typename FUNCTOR::difference_type difference_type;
    typedef typename FUNCTOR::value_type value_type;
    
    Functor f; 
    InvertFunctor(const Functor & f) : f(f) {}

    value_type operator()(difference_type const & d) const
    { return -f(d[0], d[1], d[2]); }
    
    difference_type d(difference_type const & d) const
    { return difference_type(-f.dx(d), -f.dy(d), -f.dz(d)); }    

    difference_type dir(difference_type const & v) const
    { return d(v); }    

    bool isInside(difference_type const & d) const
    { return f.isInside(d); }
  };

  template<typename FUNCTOR>
  struct TestFunctor
  {
    typedef FUNCTOR Functor;
    typedef typename FUNCTOR::difference_type difference_type;
    typedef typename FUNCTOR::value_type value_type;
    
    Functor f; 
    Functor gm; 
    TestFunctor(const Functor & f, const Functor & gm) 
      : f(f),
    gm(gm)
    {}
    
    value_type operator()(difference_type const & d) const
    { return -gm(d); }
    
    difference_type d(difference_type const & d) const
    { return -gm.d(d); }    
    
    difference_type dir(difference_type const & v) const
    { 
      difference_type df(f.d(v));
      df /= norm(df);
      double o = dot(df,d(v));
      if(o < 0.0) {
    //std::cerr << df << std::endl;
    return -df;
      }
      else {
    vigra_precondition(dot(-df,d(v))<0.0,"fail");
    return  df;
      }
      //return difference_type(1,0,0);
    }    

    bool isInside(difference_type const & d) const
    { return f.isInside(d); }
  };

  /**
   * BFGS (Broydon-Fletcher-Goldfarb-Shanno) Method from Numerical Recipes
   * (2001, p. 430 et seqq) with modifications to a generic template based style
   * using linear algebra from vigra (concerning types and simple functions)
   * and for small dimensions. 
   * The user must supply two functors, one for the to be omptimized function and 
   * one for its gradient.
   */
  template<typename LINESEARCH>
  class BFGS {
  public:
    typedef LINESEARCH LineSearch;    
    typedef typename LineSearch::Functor Functor;
    typedef typename Functor::difference_type Vector;
    typedef typename Vector::value_type ValueType;
    typedef typename NumericTraits<ValueType>::RealPromote Real;
    typedef linalg::Matrix<Real> Matrix;

    enum {N = Vector::static_size};
    
  private:
    Functor f;
    LineSearch lineSearch;
    int max_iter;
    Real eps;
    Real tolx;
    Real tolg;
    Real stpmax;
    
  public:
    BFGS(const Functor & f, int max_iter) 
      : f(f),
    lineSearch(f),
    max_iter(max_iter),
    tolx(4*std::numeric_limits<Real>::epsilon()),
    tolg(1E-5),
    stpmax(100.0)
    {}

    int operator()(Vector&, int&);
  }; // class BFGS

  namespace detail {
    template<typename MATRIX, typename VECTOR, typename REAL, int N>
    inline  void updateIverseHessianBFGS
    (MATRIX & H, 
     const VECTOR & xi, const VECTOR & dg, const VECTOR & hdg, 
     REAL fac, REAL fad, REAL fae)
    {
      //std::cerr <<  N << std::endl;
   
      for(unsigned int i=0; i < N; i++) {
    for(unsigned int j=i; j < N; j++) {
      H(i,j) += 
        fac*xi[i]*xi[j] - 
        fad*hdg[i]*hdg[j] + 
        fae*dg[i]*dg[j];
      H(j,i) = H(i,j);
    }
      }
    }
    
//     template<typename MATRIX, typename VECTOR, typename REAL>
//     inline  void updateHessianBFGS<MATRIX, VECTOR, REAL, 3>
//     (MATRIX & H, 
//      VECTOR & xi, VECTOR & dg, VECTOR & hdg, 
//      REAL fac, REAL fad, REAL fae)
//     {
//       std::cout << "jo" << std::endl;
//       H(0,0) += 
// 	fac*xi[0]*xi[0] - 
// 	fad*hdg[0]*hdg[0] + 
// 	fae*dg[0]*dg[0];
//       H(0,1) += 
// 	fac*xi[0]*xi[1] - 
// 	fad*hdg[0]*hdg[1] + 
// 	fae*dg[0]*dg[1];
//       H(0,2) += 
// 	fac*xi[0]*xi[2] - 
// 	fad*hdg[0]*hdg[2] + 
// 	fae*dg[0]*dg[2];

//       H(1,1) += 
// 	fac*xi[1]*xi[1] - 
// 	fad*hdg[1]*hdg[1] + 
// 	fae*dg[1]*dg[1];
//       H(1,2) += 
// 	fac*xi[1]*xi[2] - 
// 	fad*hdg[1]*hdg[2] + 
// 	fae*dg[1]*dg[2];

//       H(2,2) += 
// 	fac*xi[2]*xi[2] - 
// 	fad*hdg[2]*hdg[2] + 
// 	fae*dg[2]*dg[2];

//       H(1,0) = H(0,1);
//       H(2,0) = H(0,2);
//       H(2,1) = H(1,2);
//     }
  }

  template<typename LINESEARCH>
  int BFGS<LINESEARCH>::operator()
    (typename BFGS<LINESEARCH>::Vector &p, int &iter) 
  {
    Real den, fret, maxi, fac, fae, fad, sumdg, sumxi,alpha;
    Vector dg, hdg, pnew;

    const Real   zero = NumericTraits<Real>::zero();
    const Real   one  = NumericTraits<Real>::one();
    const Vector ones = NumericTraits<Vector>::one();

    Matrix H(N,N); for(int i = 0; i < N;++i) H(i,i)=one;
    Real fp = f(p);
    Vector gp = f.d(p);
    Vector xi = -gp;
    //Real dp = dot(p,p);

    lineSearch.init(p);

    for(unsigned int i = 1; true; ++i) {
      iter = i;
      //alpha = lineSearch(p, xi, lstpmax);
      alpha = lineSearch(p, xi);
      pnew = p + alpha*xi;
      //fp = lineSearch.lastFunctionValue;
      
      fp = f(pnew);
      xi = pnew - p;
      p  = pnew;

      //test for convergence
      if(i == max_iter) {
    //std::cerr << 1 << std::endl;
    return 0;
      }
      maxi = max(abs(xi)*max(abs(p),ones));
      if( maxi < tolx && maxi > zero) {
    //std::cerr << 2 << std::endl;
    return 0;
      }
      dg = gp;
      gp = f.d(p);
      den = VIGRA_CSTD::max(fret, one);
      maxi = max(abs(gp)*max(abs(p),ones)/den);
      if( maxi < tolg && maxi > zero) {
    //std::cerr << 3 << std::endl;
    return 0;
      }

      dg = gp - dg;
      hdg = H * dg;
      fac = dot(dg,xi);
      fae = dot(dg,hdg);
      fae = dot(dg,hdg);
      sumdg = dot(dg,dg);
      sumxi = dot(xi,xi);
      if(fac > sqrt(eps*sumdg*sumxi)) {
    fac = one/fac;
    fad = one/fae;
    dg  = fac*xi - fad*hdg;
    detail::updateIverseHessianBFGS<Matrix,Vector,Real,N>
      (H,xi,dg,hdg,fac,fad,fae);
      }
      xi = -(H*gp);
//       std::cerr << H  << std::endl;
//       std::cerr << gp  << std::endl;
    }
  }

  template<typename LINESEARCH>
  class CLM {
  public:
    typedef LINESEARCH LineSearch;    
    typedef typename LineSearch::Functor Functor;
    typedef typename Functor::difference_type Vector;
    typedef typename Vector::value_type ValueType;
    typedef typename NumericTraits<ValueType>::RealPromote Real;
    typedef linalg::Matrix<Real> Matrix;

    enum {N = Vector::static_size};
    
  private:
    Functor f;
    LineSearch lineSearch;
    int max_iter;
    Real eps;
    Real tolx;
    Real tolg;
    Real stpmax;
    
  public:
    CLM(const Functor & f, int max_iter) 
      : f(f),
    lineSearch(f),
    max_iter(max_iter),
    tolx(4*std::numeric_limits<Real>::epsilon()),
    tolg(1E-5),
    stpmax(100.0)
    {}

    int operator()(Vector&, int&);
  }; // class  CLM

  template<typename LINESEARCH>
  int CLM<LINESEARCH>::operator()
    (typename CLM<LINESEARCH>::Vector &p, int &iter) 
  {
    Real den, fret, maxi, fac, fae, fad, sumdg, sumxi,alpha;
    Vector dg, hdg, pnew;

    const Real   zero = NumericTraits<Real>::zero();
    const Real   one  = NumericTraits<Real>::one();
    const Vector ones = NumericTraits<Vector>::one();

    Matrix H(N,N); for(int i = 0; i < N;++i) H(i,i)=one;
    Real fp = f(p);
    Vector gp = f.d(p);
    Vector xi = -gp;
    //Real dp = dot(p,p);

    lineSearch.init(p);

    for(unsigned int i = 1; true; ++i) {
      iter = i;
      //alpha = lineSearch(p, xi, lstpmax);
      alpha = lineSearch(p, xi);
      pnew = p + alpha*xi;
      //fp = lineSearch.lastFunctionValue;
      
      fp = f(pnew);
      xi = pnew - p;
      p  = pnew;

      //test for convergence
      if(i == max_iter) {
    //std::cerr << 1 << std::endl;
    return 0;
      }
      maxi = max(abs(xi)*max(abs(p),ones));
      if( maxi < tolx && maxi > zero) {
    //std::cerr << 2 << std::endl;
    return 0;
      }
      dg = gp;
      gp = f.d(p);
      den = VIGRA_CSTD::max(fret, one);
      maxi = max(abs(gp)*max(abs(p),ones)/den);
      if( maxi < tolg && maxi > zero) {
    //std::cerr << 3 << std::endl;
    return 0;
      }

      dg = gp - dg;
      hdg = H * dg;
      fac = dot(dg,xi);
      fae = dot(dg,hdg);
      fae = dot(dg,hdg);
      sumdg = dot(dg,dg);
      sumxi = dot(xi,xi);
      if(fac > sqrt(eps*sumdg*sumxi)) {
    fac = one/fac;
    fad = one/fae;
    dg  = fac*xi - fad*hdg;
    detail::updateIverseHessianBFGS<Matrix,Vector,Real,N>
      (H,xi,dg,hdg,fac,fad,fae);
      }
      xi = -(H*gp);
//       std::cerr << H  << std::endl;
//       std::cerr << gp  << std::endl;
    }
  }

  struct NoCall {
    template<typename TYPE>
    void operator()(TYPE const & p){}
  };

  /**
   */
  template<typename LINESEARCH, typename CALLABLE = NoCall>
  class CG {
  public:
    typedef LINESEARCH LineSearch;    
    typedef typename LineSearch::Functor Functor;
    typedef typename Functor::difference_type Vector;
    typedef typename Vector::value_type ValueType;
    typedef typename NumericTraits<ValueType>::RealPromote Real;
    typedef linalg::Matrix<Real> Matrix;

    enum {N = Vector::static_size};
    
  private:
    Functor f;
    LineSearch lineSearch;
    int max_iter;
    Real eps;
    Real tol;
    CALLABLE * call;

    void record(Vector const & p) {
      if(call != NULL)
    (*call)(p);
    }
    
  public:
    CG(const Functor & f, int max_iter, CALLABLE * call_ = NULL) 
      : f(f),
    lineSearch(f),
    max_iter(max_iter),
    tol(0.000001),
    call(call_)
    {}

    int operator()(Vector&, int&);
  }; // class ConjugateGradient
  
  
  template<typename LINESEARCH, typename CALLABLE>
  int CG<LINESEARCH, CALLABLE>::operator()
    (typename CG<LINESEARCH,CALLABLE>::Vector &p, int &iter) 
  {
    const Real EPS=1.0e-18;
    int j,its;
    Real gg,gam,fp,dgg, alpha, fret;

    fp=f(p);
    
    Vector rp,pnew;
    Vector xi = -f.dir(p);
    Vector g = xi;
    Vector h = xi;

    lineSearch.init(p);
    record(p);
    
    for(its = 1; its <= max_iter; ++its) {
      iter=its;
      //       if(squaredNorm(xi)<ftol) {
      // 	return 0;
      //       }
      
      //dlinmin(p, xi, fret);
      alpha = lineSearch(p, xi);
      p    += alpha*xi;
      record(p);
      fret  = f(p);

      if (2.0*VIGRA_CSTD::abs(fret-fp) <= 
      tol*(VIGRA_CSTD::abs(fret)+VIGRA_CSTD::abs(fp)+EPS)) {
    return 0;
      }
      fp=fret;
      
      xi = f.dir(p);
      
      gg = dot(g,g);
      dgg = dot(xi,xi);
      
      if (gg == 0.0) {
    return 1;
      }
      
      gam=dgg/gg;
      g = -xi;
      xi=h=g + gam*h;
    }
    return -1;
  }

  //critical point search; extended version of Hans' code
  template<typename FUNCTOR>
  class NewtonCriticalPoint
  {
  public:
    typedef FUNCTOR Functor;
    typedef typename Functor::difference_type Vector;
    typedef typename Vector::value_type ValueType;
    typedef typename NumericTraits<ValueType>::RealPromote Real;

    enum PointType { Minimum = -1, Saddle, Maximum, Failed = 999, Seen };

  private:
    const Functor & func;
    const Real stepEpsilon;
  public:

    NewtonCriticalPoint(const Functor & func, Real stepEpsilon)
      : func(func),
    stepEpsilon(stepEpsilon)
    {}
    
    // 	     bool (*func)(double, double, double, void *)=NULL, 
    // 	     void * data = NULL)
    
    
    int operator()(Vector& p,  int& interc)
    {
      Real zero = NumericTraits<Real>::zero();
      
      int width = func.width();
      int height = func.height();
      int depth = func.depth();
      
//     xx = x;
//     yy = y;
//     zz = z;
      Real & xx = p[0];
      Real & yy = p[1];
      Real & zz = p[2];
      Real sxx, syy, szz, stepEpsilon2 = stepEpsilon * stepEpsilon;
      for(int i=0; i<100; ++i) // do at most 100 iterations
    {
      Real dx = func.dx(xx, yy, zz);
      Real dy = func.dy(xx, yy, zz);
      Real dz = func.dz(xx, yy, zz);
      
      Real dxx = func.dxx(xx, yy, zz);
      Real dxy = func.dxy(xx, yy, zz);
      Real dxz = func.dxz(xx, yy, zz);
      Real dyy = func.dyy(xx, yy, zz);
      Real dyz = func.dyz(xx, yy, zz);
      Real dzz = func.dzz(xx, yy, zz);

      Real det = 
        dxx*dyy*dzz + dxy*dyz*dxz + dxz*dxy*dyz -
        dxz*dyy*dxz - dyz*dyz*dxx - dzz*dxy*dxy;
    
      if(det != zero)
        {
          sxx = ((dyz*dyz-dyy*dzz)*dx + (dxy*dzz-dxz*dyz)*dy + (dyy*dxz-dxy*dyz)*dz) / det;
          syy = ((dxy*dzz-dyz*dxz)*dx + (dxz*dxz-dxx*dzz)*dy + (dxx*dyz-dxz*dxy)*dz) / det;
          szz = ((dyy*dxz-dxy*dyz)*dx + (dxx*dyz-dxy*dxz)*dy + (dxy*dxy-dxx*dyy)*dz) / det;
          xx += sxx;
          yy += syy;
          zz += szz;
          if(!func.isInside(xx, yy, zz)) // changed from isValid
        return Failed; // coordinates out of range
        }
      else
        {
          sxx = syy = szz = 0.0;
        }
    
      // 	if(func != NULL) {
      // 	  bool seen = func(xx,yy,zz, data);
      // 	  if(seen)
      // 	    return Seen;
      // 	}
      
      double dist2 = sxx*sxx + syy*syy + szz*szz;
      if(dist2 < stepEpsilon2) // convergence
        {
          // FIXME: don't reuse stepEpsilon for a third purpose (I
          // already refactored the other reuses away)
          if(xx < -stepEpsilon || xx > (double)(width)-1.0+stepEpsilon ||
         yy < -stepEpsilon || yy > (double)(height)-1.0+stepEpsilon ||
         zz < -stepEpsilon || zz > (double)(depth)-1.0+stepEpsilon)
        {
          return Failed; // coordinates out of range
        }
          if(det == zero)
        {
          if(dx == zero && dy == zero && dz == zero) //FIXMENEW: zero precession
            {
              return Saddle;
            }
          return Failed;
          }
          else if(det < zero)
        {
          return Saddle;
          }
          else if(dxx + dyy + dzz > zero)
        {
          return Minimum;
        }
          else
        {
          return Maximum;
          }
        }
    }
      return Failed;
    }
  };
} // namespace vigra
#endif
