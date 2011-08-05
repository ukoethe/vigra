#ifndef VIGRA_ANALYTICVOLUMEVIEW_HXX
#define VIGRA_ANALYTICVOLUMEVIEW_HXX

#include "vigra/mathutil.hxx"
#include "vigra/tinyvector.hxx"
#include "vigra/fixedpoint.hxx"
#include <vigra/multi_pointoperators.hxx>
#include <vigra/multi_impex.hxx>
#include "quaternion.hxx"

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

// #include "boost/units/pow.hpp"
// #include "boost/math/special_functions/bessel.hpp"

#define pow2(x) gsl_pow_2(x) //FIXME:replace with better evaluation
#define pow3(x) gsl_pow_3(x) //FIXME:replace with better evaluation
#define pow4(x) gsl_pow_4(x) //FIXME:replace with better evaluation
#define cyl_bessel_i(n,x) gsl_sf_bessel_In(n,x)

namespace vigra {

  template <typename VALUETYPE>
  class GaussianVolume
  {
  public:
    typedef typename NumericTraits<VALUETYPE>::RealPromote real_type;
    typedef VALUETYPE value_type;
    typedef TinyVector<double, 3> size_type;
    typedef size_type::value_type size_value_type;
    typedef TinyVector<double, 3> difference_type;    
    typedef difference_type::value_type difference_value_type;    
    typedef Quaternion<real_type> quaternion;

    template<typename SrcShape>
    GaussianVolume(real_type s_, SrcShape shape, //real_type w_, 
           difference_type v_, difference_type trans) 
      : s(s_),
    ss(s_*s_),
    w(shape[0]), 
    h(shape[1]), 
    t(shape[2]),
    //c(difference_type((shape-NumericTraits<SrcShape>::one())/2.)),
    c(trans),
    //q(cos(w_/2.), sin(w_/2.)*v_)
    q(cos(v_[0]/2.)*cos(v_[1]/2.)*cos(v_[2]/2.)+sin(v_[0]/2.)*sin(v_[1]/2.)*sin(v_[2]/2.),
      difference_type(sin(v_[0]/2.)*cos(v_[1]/2.)*cos(v_[2]/2.)-cos(v_[0]/2.)*sin(v_[1]/2.)*sin(v_[2]/2.),
              cos(v_[0]/2.)*sin(v_[1]/2.)*cos(v_[2]/2.)+sin(v_[0]/2.)*cos(v_[1]/2.)*sin(v_[2]/2.),
              cos(v_[0]/2.)*cos(v_[1]/2.)*sin(v_[2]/2.)-sin(v_[0]/2.)*sin(v_[1]/2.)*cos(v_[2]/2.)))

    {std::cerr << q.w() << " " << q.x() << " " << q.y() << " " << q.z();;}

    template<typename SrcShape>
    GaussianVolume(real_type s_, SrcShape shape, difference_type trans) 
      : s(s_),
    ss(s_*s_),
    w(shape[0]), 
    h(shape[1]), 
    t(shape[2]),
    //c(difference_type((shape-NumericTraits<SrcShape>::one())/2.)),
    c(trans),
    q(NumericTraits<real_type>::one(), 
      NumericTraits<difference_type>::zero())
    {}

    template<typename SrcShape>
    GaussianVolume(real_type s_, SrcShape shape) 
      : s(s_),
    ss(s_*s_),
    w(shape[0]), 
    h(shape[1]), 
    t(shape[2]),
    //c(difference_type((shape-NumericTraits<SrcShape>::one())/2.)),
    c(difference_type((shape-NumericTraits<SrcShape>::one())/2.)),
    q(NumericTraits<real_type>::one(), 
      NumericTraits<difference_type>::zero())
    {}


    virtual real_type operator()(difference_value_type x, 
                 difference_value_type y, 
                 difference_value_type z) const = 0;

    virtual real_type dx(double x, double y, double z) const = 0;

    virtual real_type dy(double x, double y, double z) const = 0;

    virtual real_type dz(double x, double y, double z) const = 0;

    virtual real_type g2(double x, double y, double z) const = 0;

    real_type operator()(const difference_type & d) const
    { return operator()(d[0],d[1],d[2]); }

    real_type operator[](difference_type const & d) const
    { return operator()(d); }

    difference_type d(difference_type const & d) const
     { return difference_type(dx(d),dy(d),dz(d)); }

    real_type dx(difference_type const & d) const
    { return dx(d[0],d[1],d[2]); }
    
    real_type dy(difference_type const & d) const
    { return dy(d[0],d[1],d[2]); }
    
    real_type dz(difference_type const & d) const
    { return dz(d[0],d[1],d[2]); }

    real_type g2(difference_type const & d) const
    { return g2(d[0], d[1], d[2]); }
        
    size_type shape() const
    { return size_type(width(),height(),depth()); }

    size_value_type width() const
    { return w; }
    
    size_value_type height() const
    { return h; }
    
    size_value_type depth() const
    { return t; }

    bool isInsideX(difference_value_type x) const
    { return x >= 0.0 && x <= width()-1.0; }

    bool isInsideY(difference_value_type y) const
    { return y >= 0.0 && y <= height()-1.0; }
    
    bool isInsideZ(difference_value_type z) const
    { return z >= 0.0 && z <= depth()-1.0; }

    bool isInside(difference_value_type x, 
          difference_value_type y, 
          difference_value_type z) const
    { return isInsideX(x)  && isInsideY(y) && isInsideZ(z); }
    
    bool isInside(difference_type const & p) const
    { return isInside(p[0],p[1],p[2]); }
   
  protected:
    real_type s,ss;
    size_value_type w, h, t;
    difference_type c;
    quaternion q;

    real_type gerr(real_type x) const
    { return 0.5 + 0.5 * gsl_sf_erf(x/(M_SQRT2*s)); }

    real_type gerr(real_type x, real_type serr) const
    { return 0.5 + 0.5 * gsl_sf_erf(x/(M_SQRT2*serr)); }
    
    void trans(difference_value_type &x, 
           difference_value_type &y, 
           difference_value_type &z) const
    {
//       quaternion v(0,difference_type(x-c[0],y-c[1],z-c[2]));
//       v = q*v*conj(q);

//       x = v.x();
//       y = v.y();
//       z = v.z();
       x = x - c[0];
       y = y - c[1];
       z = z - c[2];

    }

    void rotate(difference_value_type &x, 
        difference_value_type &y, 
        difference_value_type &z) const
    {
      quaternion v(0,difference_type(x,y,z));
      v = q*v*conj(q);
      
      x = v.x();
      y = v.y();
      z = v.z();
    }

    difference_type  rotate(difference_type const & d) const
    {
      difference_type ret = d;
      std::cerr << "b " <<ret << std::endl;
      rotate(ret[0],ret[1],ret[2]);
      std::cerr << "a " << ret << std::endl;
      return ret;
    }

  };
  
  template <typename VALUETYPE>
  class GaussianPlane : public GaussianVolume<VALUETYPE>
  {
  public:   
    typedef GaussianVolume<VALUETYPE> base_type;
    typedef GaussianPlane<VALUETYPE>  this_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::difference_value_type difference_value_type;
    
  private:
    typedef Quaternion<real_type> quaternion;

    difference_type dir;
    value_type t;
    value_type p;
    value_type q;
    
    real_type dv(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z,
         difference_value_type v) const
    {
      real_type m = dir[0]*x + dir[1]*y + dir[2]*z;
      return
    ((-p + q)*v)/(exp(pow2(t + m)/(2*base_type::ss))*sqrt(2*M_PI)*base_type::ss);
    }
    
    real_type g2v(difference_value_type x, 
          difference_value_type y, 
          difference_value_type z,
          difference_value_type v) const
    {
      real_type m = dir[0]*x + dir[1]*y + dir[2]*z;
      return
    -((pow2(p - q)*v*squaredNorm(dir)*(t + m))/
      (exp(pow2(t + m)/base_type::ss)*M_PI*base_type::ss*base_type::ss));
    }

  public:
    template<typename SrcShape>
    GaussianPlane(difference_type dir_, value_type trans, value_type sigma, 
          value_type p, value_type q, 
          SrcShape shape) 
      :  base_type(sigma,shape,
           difference_type((shape-NumericTraits<SrcShape>::one())/2.)),
// 	 dir(difference_type(1,0,0)),
     dir(dir_),
     t(trans),
     p(p),
     q(q)
    {

//       dir = base_type::rotate(dir);
//       dir[0] = cos(dir_[1])*cos(dir_[2]);
//       dir[1] = cos(dir_[1])*sin(dir_[2]);
//       dir[2] = -sin(dir_[1]);
//       dir /= norm(dir);
      //std::cerr << dir_;
    }

    real_type distance(difference_type const & v) const {
      difference_type tv(t,0,0);
//       difference_type x = v-base_type::c;
//       quaternion qc(0,x);
//       qc = base_type::q*qc*conj(base_type::q);

//       x[0] = qc.x();
//       x[1] = qc.y();
//       x[2] = qc.z();



      difference_type x = v-base_type::c;      

//       std::cerr << x << std::endl;
//       std::cerr << dot(x,dir)  << std::endl;
//       std::cerr << norm(dir)  << std::endl;

      return dot(x,dir)+t;
    }
        
    real_type operator()(difference_value_type x, 
             difference_value_type y, 
             difference_value_type z) const
    { 
      base_type::trans(x,y,z);
      //x = x+t;
      return p + (q-p)*base_type::gerr(x*dir[0] + y*dir[1] + z*dir[2]+t,
                       base_type::s);
    }

    real_type dx(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    {
      base_type::trans(x,y,z);     
      return dv(x,y,z,dir[0]);
    }
    
    real_type dy(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    {
      base_type::trans(x,y,z);     
      return dv(x,y,z,dir[1]);
    }

    real_type dz(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    {
      base_type::trans(x,y,z);     
      return dv(x,y,z,dir[2]);
    }

    real_type g2(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    {
      base_type::trans(x,y,z);
//       std::cerr << x << " " 
// 		<< y << " " 
// 		<< z  << " "
// 		<< p-q  << " "
// 		<< squaredNorm(d)  << " "
// 		<< base_type::ss 
// 		<< pow2(t + d[0]*x + d[1]*y + d[2]*z) << std::endl;
      return
    (pow2(p - q)*squaredNorm(dir))/
    (2*exp(pow2(t + dir[0]*x + dir[1]*y + dir[2]*z)/base_type::ss)*
     M_PI*base_type::ss);  
    }
    
    real_type g2x(difference_value_type x, 
          difference_value_type y, 
          difference_value_type z) const
      
    {
      base_type::trans(x,y,z);     
      return g2v(x,y,z,dir[0]);
    }

    real_type g2y(difference_value_type x, 
          difference_value_type y, 
          difference_value_type z) const
      
    {
      base_type::trans(x,y,z);     
      return g2v(x,y,z,dir[1]);
    }

    real_type g2z(difference_value_type x, 
          difference_value_type y, 
          difference_value_type z) const
      
    {
      base_type::trans(x,y,z);     
      return g2v(x,y,z,dir[2]);
    }
    
    using base_type::operator();

    //FIXME: shoudl be inherinted?
//     real_type operator()(difference_type const & d) const
//     { return operator()(d[0],d[1],d[2]); }

//     value_type dx(difference_type const & d) const
//     { return dx(d[0], d[1], d[2]); }

//     value_type dy(difference_type const & d) const
//     { return dy(d[0], d[1], d[2]); }

//     value_type dz(difference_type const & d) const
//     { return dz(d[0], d[1], d[2]); }
    
    real_type g2(difference_type const & d) const
    { return g2(d[0], d[1], d[2]); }

  };

  
  template <typename VALUETYPE>
  class GaussianCylinder : public GaussianVolume<VALUETYPE>
  {
  public:   
    typedef GaussianVolume<VALUETYPE> base_type;
    typedef GaussianPlane<VALUETYPE>  this_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::difference_value_type difference_value_type;
    
  private:
    value_type r;
    value_type p;
    value_type q;
    value_type terr, serr;
    real_type a;
    real_type b;
    real_type trans;

        
    real_type gv(real_type x, real_type y) const
    {
      real_type lsq = x*x + y*y;
      real_type l   = sqrt(x*x + y*y);
      double & ss = base_type::ss;
      return 
    (r*x*(r*cyl_bessel_i(0, (r*l)/ss) - 2*l*cyl_bessel_i(1, (r*l)/ss) + 
          r*cyl_bessel_i(2, (r*l)/ss)))/
    (2*exp((r*r + lsq)/(2*ss))*ss*ss*l);
    };
    
  public:
    template<typename SrcShape>
    GaussianCylinder(real_type r, real_type sigma, 
             real_type p, real_type q, 
             SrcShape shape) 
      :  base_type(sigma,shape),
     r(r),
     p(p),
     q(q),
     terr(1.72),
     serr(0.1),
     a(0.0),
     b(0.0)
    {}

    template<typename SrcShape>
    GaussianCylinder(real_type r, real_type sigma, 
             real_type p, real_type q, 
             SrcShape shape, 
             real_type a_, real_type b_,
             real_type trans_) 
      :  base_type(sigma,shape),
     r(r),
     p(p),
     q(q),
     terr(1.72),
     serr(0.1),
     a(a_),
     b(b_),
     trans(trans_)
    {}

    real_type distance(difference_type const & v) const {
      real_type x,y,z;
      x=v[0]+trans;y=v[1]+trans;z=v[2]+trans;
      base_type::trans(x,y,z);

      x = x*cos(b)+z*sin(b);
      y = y*cos(a)-z*cos(b)*sin(a)+x*sin(b)*sin(a);
      z = z*cos(b)*cos(a)-x*cos(a)*sin(b)+y*sin(a);

      return sqrt(x*x+y*y);
    }

    
    real_type operator()(difference_value_type x, 
             difference_value_type y, 
             difference_value_type z) const
    {
      x=x+trans;y=y+trans;z=z+trans;
      base_type::trans(x,y,z);

      x = x*cos(b)+z*sin(b);
      y = y*cos(a)-z*cos(b)*sin(a)+x*sin(b)*sin(a);
      z = z*cos(b)*cos(a)-x*cos(a)*sin(b)+y*sin(a);
      
      real_type s   = base_type::s;
      real_type ss  = base_type::ss;
      real_type lsq = x*x + y*y;

      real_type bg =  base_type::gerr(r/s - terr, serr);
      real_type bs = 1. - bg;
      real_type ds = (2*r*r)/(4*ss+r*r)*exp(-(2*lsq)/(4*s*s+r*r));
      real_type v  = cbrt((r*r)/(2*s*s + lsq));
      real_type u  = (2./3.)*s*((sqrt(ss + lsq))/(2*s*s + lsq));
      real_type dg = base_type::gerr(((v - 1.)/u) + u, 1.);

      return p + (q-p)*(ds*bs + dg*bg);
    }

    real_type dx(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    { 
      throw new std::string("The method or operation is not implemented.");
    }

    real_type dy(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    { 
      throw new std::string("The method or operation is not implemented.");
    }

    real_type dz(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    { 
      throw new std::string("The method or operation is not implemented.");
    }

    real_type g(difference_value_type x, 
        difference_value_type y, 
        difference_value_type z) const
    {
      base_type::trans(x,y,z);
      
      real_type lsq = x*x + y*y;
      real_type l   = sqrt(x*x + y*y);
      return
    (r/base_type::ss)*exp(-(r*r + lsq)/(2*base_type::ss))*
    cyl_bessel_i(1, (l*r)/(base_type::ss));
    }
    
    real_type g2(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    { return pow2(g(x,y,z)); }
    
    
    real_type gx(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    {
      base_type::trans(x,y,z);
      return(x,y);
    }
    
    real_type gy(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    {
      base_type::trans(x,y,z);
      return(y,x);
    }
    
    real_type gz(difference_value_type x, 
         difference_value_type y, 
         difference_value_type z) const
    {
      return 0.0;
    }

    using base_type::operator();
    using base_type::g2;
    using base_type::d;
    using base_type::dx;
    using base_type::dy;
    using base_type::dz;
  };

  template <typename VALUETYPE>
  class GaussianSphere : public GaussianVolume<VALUETYPE>
  {
  public:   
    typedef GaussianVolume<VALUETYPE> base_type;
    typedef GaussianPlane<VALUETYPE>  this_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::difference_value_type difference_value_type;

    
  private:
    value_type r,rr;
    value_type p;
    value_type q;
    value_type fac;
    
    real_type gerr(real_type x) const
    { return 0.5 + 0.5 * erf(x/(M_SQRT2*base_type::s)); }
    
    real_type gauss(real_type x) const
    { return (1./fac) * exp(-x*x/(2*base_type::ss)); }

    real_type dv(double x, double y, double z) const
    {  
      double m = x*x + y*y + z*z;
      double n = sqrt(m);
      
      return
    ((p-q)*x*(-(exp(pow2(r - n)/(2*base_type::ss))* n) + exp(pow2(r + n)/(2*base_type::ss))* (-2*r + n))) /
    (exp((rr + m)/base_type::ss)*fac*r*n);
    }

    real_type g2v(double x, double y, double z) const 
    {      
      double m = x*x + y*y + z*z;
      double n = sqrt(m);
      double ss = base_type::ss;
      return 
    -((pow2(p - q)*x*(exp(pow2(r - n)/(2*ss))*n - exp(pow2(r + n)/(2*ss))*
              (-2*r + n))*(-ss + m + r*n + exp((2*r*n)/ss)*
                       (-2*rr + ss - x*x - y*y - z*z + 3*r*n)))/
      (exp((3*rr + 2*r*n + 3*m)/(2*ss))*M_PI*rr*ss*ss*n));
    }

    real_type gv(double x, double y, double z) const 
    {      
      double m = x*x + y*y + z*z;
      double R = sqrt(m);
      double ss = base_type::ss;

      return 
    (pow2(p-q)*(-(exp(pow2(r-R)/(2*ss))*R)+exp(pow2(r+R)/(2*ss))*
            (-2*r+R))*(r*R+m-ss+exp((2*r*R)/ss)*
                   (-2*pow2(r)+3*r*R-m+ss))*x)/
    (exp((3*pow2(r)+2*r*R+3*m)/(2*ss))*sqrt(2*M_PI)*pow2(r)*R*
     sqrt((pow2(p-q)*pow2(exp(pow2(r+R)/(2*ss))*(2*r-R)+
                  exp(pow2(r-R)/(2*ss))*R))/
          (exp((2*(pow2(r)+m))/ss)*pow2(r)*ss))*ss*ss);  
    }
    

  public:
    template<typename SrcShape>
    GaussianSphere(value_type radius, value_type sigma, 
           value_type p, value_type q, 
           SrcShape shape, difference_type trans) 
      :  base_type(sigma,shape,trans),
     r(radius),
     rr(radius*radius),
     p(p),
     q(q),
     fac(sqrt(2*M_PI)*sigma)
    {//std::cerr << base_type::c << std::endl;
    }

    template<typename SrcShape>
    GaussianSphere(value_type radius, value_type sigma, 
           value_type p, value_type q, 
           SrcShape shape) 
      :  base_type(sigma,shape),
     r(radius),
     rr(radius*radius),
     p(p),
     q(q),
     fac(sqrt(2*M_PI)*sigma)
    {}



    value_type operator()(double x, double y, double z) const
    { 
      base_type::trans(x,y,z);
      double n = sqrt(x*x + y*y + z*z);
      real_type ss = base_type::ss;
      real_type s = base_type::s;

      if(n < 1e-05) 
    return
      p+(-p+q)*(-((sqrt(2/M_PI)*r)/(exp(r*r/(2*ss))*s)) + erf(r/(M_SQRT2*s)));
      
      return p+(-p+q)*(s/(exp(pow2(r + n)/(2*ss))*sqrt(2*M_PI)*n) - 
               (exp((2*r*n)/ss - pow2(r + n)/(2*ss))*s)/(sqrt(2*M_PI)*n) + 
               erf((r - n)/(M_SQRT2*s))/2 + erf((r + n)/(M_SQRT2*s))/2);
    }

    difference_type d(difference_type const & d) const
    { 
      return difference_type(dx(d),dy(d),dz(d));
    }
    
    value_type dx(double x, double y, double z) const
    {  
      base_type::trans(x,y,z);
      return dv(x,y,z);
    }

    value_type dy(double x, double y, double z) const
    {  
      base_type::trans(x,y,z);
      return dv(y,x,z);
    }

    value_type dz(double x, double y, double z) const
    {  
      base_type::trans(x,y,z);
      return dv(z,x,y);
    }

    real_type g(difference_value_type x, 
        difference_value_type y, 
        difference_value_type z) const
    { 
      base_type::trans(x,y,z);      
      double m = x*x + y*y + z*z;
      double R = sqrt(m);
      real_type ss = base_type::ss;
      return
    sqrt((pow2(p-q)*pow2(exp((2*r*R)/ss)*(2*r-R)+R))/
         (exp(pow2(r+R)/ss)*r*r*ss))/sqrt(2*M_PI);

}

    value_type gx(double x, double y, double z) const 
    {
      base_type::trans(x,y,z);     
      return gv(x,y,z);
    }

    value_type gy(double x, double y, double z) const 
    {
      base_type::trans(x,y,z);     
      return gv(y,x,z);
    }

    value_type gz(double x, double y, double z) const 
    {
      base_type::trans(x,y,z);     
      return gv(z,x,y);
    }
    
    value_type g2(double x, double y, double z) const 
    {
      base_type::trans(x,y,z);      
      double m = x*x + y*y + z*z;
      double n = sqrt(m);
      real_type ss = base_type::ss;
//       return
// 	(m - 2*exp((2*r*n)/ss) * (m - 2*r*n) +  exp((4*r*n)/ss)*(m + 4*r*(r - n)))/
// 	(2*exp((rr + m + 2*r*n)/ss)*M_PI*rr*ss);

      return
    (pow2(p - q)*pow2(exp(pow2(r-n)/(2*ss))*n - exp(pow2(r + n)/(2*ss))* (-2*r + n)))/
    (2*exp((2*(rr + m))/ss)*M_PI*rr*ss);
    }

    value_type g2x(double x, double y, double z) const 
    {
      base_type::trans(x,y,z);     
      return g2v(x,y,z);
    }

    value_type g2y(double x, double y, double z) const 
    {
      base_type::trans(x,y,z);     
      return g2v(y,x,z);
    }

    value_type g2z(double x, double y, double z) const 
    {
      base_type::trans(x,y,z);     
      return g2v(z,x,y);
    }

    using base_type::operator();
    using base_type::g2;
    using base_type::dx;
    using base_type::dy;
    using base_type::dz;
  };

  template <typename VALUETYPE>
  class GaussianRing
  {
  public:   
    typedef typename NumericTraits<VALUETYPE>::RealPromote InternalValue;
    typedef VALUETYPE value_type;
    typedef TinyVector<double, 3> size_type;
    typedef TinyVector<double, 3> difference_type;
    
  private:
    unsigned int w_, h_, d_;
    value_type r_;
    value_type s_;
    difference_type c_;

    value_type ssq;
    value_type rsq;
    value_type rds;
    value_type fac;
    
  public:
    template<typename SrcShape>
    GaussianRing(value_type radius, value_type sigma, SrcShape shape) 
      : r_(radius),
    s_(sigma),
    w_(shape[0]), 
    h_(shape[1]), 
    d_(shape[2]),
    c_(difference_type(shape[0]/2.,shape[1]/2.,shape[2]/2.)),
    ssq(sigma*sigma),
    rsq(radius*radius),
    rds(radius/sigma),
    fac(vigra::sqrt(2*M_PI)*rds)
    {}
    
    value_type operator()(difference_type const & d) const
    { 
      double f1 = fac*exp((vigra::squaredNorm(d) + r_*r_) / (-2.*ssq));
      double f2 = bessi0(rds*vigra::sqrt(d[0]*d[0] + d[1]*d[1]));
//       std::cerr << rds << " " << vigra::sqrt(d[0]*d[0] + d[1]*d[1]) << " " << rds*vigra::sqrt(d[0]*d[0] + d[1]*d[1]) << std::endl;
//       std::cerr << f1 << " " << f2 << std::endl;
      return f1*f2;
    }

    value_type operator[](difference_type const & d) const
    { return operator()(d); }
    
//     value_type dx(double x, double y, double z) const
//     { 
//       difference_type v = d - c_;
//       difference_type w = (v / vigra::norm(v)) * r_;
//       value_type n = vigra::norm(v-w);
//       return fac*exp((n*n)/(-2.*ssq));
//     }

    
//     value_type dy(double x, double y, double z) const
//     { return operator()(x, y, z, 0, 1, 0); }
    
//     value_type dz(double x, double y, double z) const
//     { return operator()(x, y, z, 0, 0, 1); }


    unsigned int width() const
    { return w_; }
    
    unsigned int height() const
    { return h_; }
    
    unsigned int depth() const
    { return d_; }

    bool isInsideX(double x) const
    {
      return x >= 0.0 && x <= width()-1.0;
    }
    bool isInsideY(double y) const
    {
      return y >= 0.0 && y <= height()-1.0;
    }
    
    bool isInsideZ(double z) const
    {
      return z >= 0.0 && z <= depth()-1.0;
    }

    bool isInside(double x, double y, double z) const
    {
      return isInsideX(x)  && isInsideY(y) && isInsideZ(z);
    }
    
    bool isInside(difference_type const & p) const
    {
      return isInside(p[0],p[1],p[2]);
    }
  };

  template <typename VALUETYPE>
  class GaussianLine
  {
  public:   
    typedef typename NumericTraits<VALUETYPE>::RealPromote InternalValue;
    typedef VALUETYPE value_type;
    typedef TinyVector<double, 3> size_type;
    typedef TinyVector<double, 3> difference_type;
    
  private:
    unsigned int w_, h_, d_;
    value_type sx_;
    value_type sy_;
    difference_type c_;
    
  public:
    template<typename SrcShape>
    GaussianLine(value_type sigmax, value_type sigmay, SrcShape shape) 
      : sx_(sigmax),
    sy_(sigmay),
    w_(shape[0]), 
    h_(shape[1]), 
    d_(shape[2]),
    c_(difference_type(shape[0]/2.,shape[1]/2.,shape[2]/2.))
    {}

    difference_type shape()
    {  return difference_type(width(),height(),depth());
    }

    value_type operator()(double x, double y, double z) const
    {
      x = x - c_[0];
      y = y - c_[1];
      return exp(-x*x/(2.*sx_*sx_))*exp(-y*y/(2.*sy_*sy_));
    }
    
    value_type operator()(difference_type const & d) const
    { 
      return operator()(d[0],d[1],d[2]);
    }

    value_type operator[](difference_type const & d) const
    { return operator()(d); }

    difference_type d(difference_type const & d) const
    { 
      return difference_type(dx(d),dy(d),dz(d));
    }

    value_type dx(difference_type const & d) const
    { 
      return dx(d[0],d[1],d[2]);
    }

    value_type dy(difference_type const & d) const
    { 
      return dy(d[0],d[1],d[2]);
    }

    value_type dz(difference_type const & d) const
    { 
      return dz(d[0],d[1],d[2]);
    }
    
    value_type dx(double x, double y, double z) const
    { 
      x = x - c_[0];
      y = y - c_[1];

      return -(exp(-x*x/(2*sx_*sx_) - y*y/(2*sy_*sy_))*x)/(sx_*sx_);
    }

    value_type dy(double x, double y, double z) const
    { 
      x = x - c_[0];
      y = y - c_[1];

      return -(exp(-x*x/(2*sx_*sx_) - y*y/(2*sy_*sy_))*y)/(sy_*sy_);
    }

    value_type dz(double x, double y, double z) const
    { 
      return 0.0;
    }

    value_type dxx(double x, double y, double z) const
    { 
      x = x - c_[0];
      y = y - c_[1];

      return 
    -(exp(-x*x/(2*sx_*sx_) - y*y/(2*sy_*sy_))/(sx_*sx_))
    +(exp(-x*x/(2*sx_*sx_) - y*y/(2*sy_*sy_))*x*x)
    /(sx_*sx_*sx_*sx_); 
    }

    value_type dyy(double x, double y, double z) const
    {       
      x = x - c_[0];
      y = y - c_[1];

      return 
    -(exp(-x*x/(2*sx_*sx_) - y*y/(2*sy_*sy_))/(sy_*sy_))
    +(exp(-x*x/(2*sx_*sx_) - y*y/(2*sy_*sy_))*y*y)
    /(sy_*sy_*sy_*sy_); 
    }

    value_type dzz(double x, double y, double z) const
    { return 0.0; }
   
    value_type dxy(double x, double y, double z) const
    {
      x = x - c_[0];
      y = y - c_[1];
 
      return
    (exp(-x*x/(2*sx_*sx_) - y*y/(2*sy_*sy_))*x*y)
    /(sx_*sx_*sy_*sy_);
    }
    
    value_type dxz(double x, double y, double z) const
    { return 0.0; }
    
    
    value_type dyz(double x, double y, double z) const
    { return 0.0; }
   

    unsigned int width() const
    { return w_; }
    
    unsigned int height() const
    { return h_; }
    
    unsigned int depth() const
    { return d_; }

    bool isInsideX(double x) const
    {
      return x >= 0.0 && x <= width()-1.0;
    }
    bool isInsideY(double y) const
    {
      return y >= 0.0 && y <= height()-1.0;
    }
    
    bool isInsideZ(double z) const
    {
      return z >= 0.0 && z <= depth()-1.0;
    }

    bool isInside(double x, double y, double z) const
    {
      return isInsideX(x)  && isInsideY(y) && isInsideZ(z);
    }
    
    bool isInside(difference_type const & p) const
    {
      return isInside(p[0],p[1],p[2]);
    }
  };

  template <typename VALUETYPE>
  class GaussianBranch
  {
  public:   
    typedef typename NumericTraits<VALUETYPE>::RealPromote InternalValue;
    typedef VALUETYPE value_type;
    typedef TinyVector<double, 3> size_type;
    typedef TinyVector<double, 3> difference_type;
    
  private:
    unsigned int w_, h_, d_;
    value_type s_;
    value_type p_;
    value_type q_;
    value_type qx_;
    value_type qy_;
    difference_type c_;
    
  public:
    template<typename SrcShape>
    GaussianBranch(value_type sigma, value_type p, value_type q, value_type qx, value_type qy, SrcShape shape) 
      : s_(sigma),
    p_(p),
    q_(q),
    qx_(qx),
    qy_(qy),
    w_(shape[0]), 
    h_(shape[1]), 
    d_(shape[2]),
    c_(difference_type(shape[0]/2.,shape[1]/2.,shape[2]/2.))
    {}

    difference_type shape()
    {
      return difference_type(width(),height(),depth());
    }

    value_type operator()(double x, double y, double z) const
    {
      x = x - c_[0];
      y = y - c_[1];
      z = z - c_[2];
      if(z < 0) 
    return exp(-(x*x)/(2.*s_*s_) - (y*y)/(2.*s_*s_));
     
      return  p_*(exp(-(y*y + pow2(x + q_*z*z))/(2.*s_*s_))) + 
    (1 - p_)*(exp(-(pow2(y + qy_*z*z) + pow2(x + qx_*z*z))/(2.*s_*s_)));
    }
    
    value_type operator()(difference_type const & d) const
    { 
      return operator()(d[0],d[1],d[2]);
    }

    value_type operator[](difference_type const & d) const
    { return operator()(d); }

    difference_type d(difference_type const & d) const
    { 
      return difference_type(dx(d),dy(d),dz(d));
    }

    value_type dx(difference_type const & d) const
    { 
      return dx(d[0],d[1],d[2]);
    }

    value_type dy(difference_type const & d) const
    { 
      return dy(d[0],d[1],d[2]);
    }

    value_type dz(difference_type const & d) const
    { 
      return dz(d[0],d[1],d[2]);
    }
    
    value_type dx(double x, double y, double z) const
    { 
      x = x - c_[0];
      y = y - c_[1];
      z = z - c_[2];
      
      return
    -((p_*(x + q_*z*z))/(exp((y*y + pow2(x + q_*z*z))/(2*s_*s_))*s_*s_)) + 
    ((-1 + p_)*(x + qx_*z*z))/(exp((pow2(x + qx_*z*z) + pow2(y + qy_*z*z))/(2*s_*s_))*s_*s_);
    }

    value_type dy(double x, double y, double z) const
    { 
      x = x - c_[0];
      y = y - c_[1];
      z = z - c_[2];

      return
    -((p_*y)/(exp((y*y + pow2(x + q_*z*z))/(2*s_*s_))*s_*s_)) + ((-1 + p_)*(y + qy_*z*z))/
    (exp((pow2(x + qx_*z*z) + pow2(y + qy_*z*z))/(2*s_*s_))*s_*s_);

    }

    value_type dz(double x, double y, double z) const
    { 
      x = x - c_[0];
      y = y - c_[1];
      z = z - c_[2];

      return 
    (-2*p_*q_*z*(x + q_*z*z))/(exp((y*y + pow2(x + q_*z*z))/(2*s_*s_))*s_*s_) + 
    (2*(-1 + p_)*z*(qx_*x + qy_*y + (qx_*qx_ + qy_*qy_)*z*z))/(exp((pow2(x + qx_*z*z) + pow2(y + qy_*z*z))/(2*s_*s_))*s_*s_);
    }

    value_type dxx(double x, double y, double z) const
    { 
      x = x - c_[0];
      y = y - c_[1];
    }

    value_type dyy(double x, double y, double z) const
    {       
      x = x - c_[0];
      y = y - c_[1];
    }

    value_type dzz(double x, double y, double z) const
    { return 0.0; }
   
    value_type dxy(double x, double y, double z) const
    {
      x = x - c_[0];
      y = y - c_[1];
    }
    
    value_type dxz(double x, double y, double z) const
    { return 0.0; }
    
    
    value_type dyz(double x, double y, double z) const
    { return 0.0; }
   

    unsigned int width() const
    { return w_; }
    
    unsigned int height() const
    { return h_; }
    
    unsigned int depth() const
    { return d_; }

    bool isInsideX(double x) const
    {
      return x >= 0.0 && x <= width()-1.0;
    }
    bool isInsideY(double y) const
    {
      return y >= 0.0 && y <= height()-1.0;
    }
    
    bool isInsideZ(double z) const
    {
      return z >= 0.0 && z <= depth()-1.0;
    }

    bool isInside(double x, double y, double z) const
    {
      return isInsideX(x)  && isInsideY(y) && isInsideZ(z);
    }
    
    bool isInside(difference_type const & p) const
    {
      return isInside(p[0],p[1],p[2]);
    }
  };

  template <typename VALUETYPE>
  class GaussianShell
  {
  public:   
    typedef GaussianShell<VALUETYPE> THISTYPE;
    typedef typename NumericTraits<VALUETYPE>::RealPromote InternalValue;
    typedef VALUETYPE value_type;
    typedef TinyVector<double, 3> size_type;
    typedef TinyVector<double, 3> difference_type;
    
  private:
    unsigned int w_, h_, d_;
    value_type r_;
    value_type s_;
    difference_type c_;

    value_type ssq;
    value_type fac;
    
  public:
    template<typename SrcShape>
    GaussianShell(value_type radius, value_type sigma, SrcShape shape) 
      : r_(radius),
    s_(sigma),
    w_(shape[0]), 
    h_(shape[1]), 
    d_(shape[2]),
    c_(difference_type(shape[0]/2.,shape[1]/2.,shape[2]/2.)),
    ssq(sigma*sigma),
    fac(1./(vigra::sqrt(2*M_PI)*sigma))
    {}

    difference_type shape()
    {
      return difference_type(width(),height(),depth());
    }

    value_type operator()(double x, double y, double z) const
    {
      return operator()(difference_type(x,y,z));
    }
    
    value_type operator()(difference_type const & d) const
    { 
      difference_type v = d - c_;
      difference_type w = (v / vigra::norm(v)) * r_;
      value_type n = vigra::norm(v-w);
      return fac*exp((n*n)/(-2.*ssq));
    }

    value_type operator[](difference_type const & d) const
    { return operator()(d); }

    difference_type d(difference_type const & d) const
    { 
      return difference_type(dx(d),dy(d),dz(d));
    }

    value_type dx(difference_type const & d) const
    { 
      return dx(d[0],d[1],d[2]);
    }

    value_type dy(difference_type const & d) const
    { 
      return dy(d[0],d[1],d[2]);
    }

    value_type dz(difference_type const & d) const
    { 
      return dz(d[0],d[1],d[2]);
    }
    
    value_type dx(double x, double y, double z) const
    {  

      x = x - c_[0];
      y = y - c_[1];
      z = z - c_[2];

      double xx=x*x;
      double xxxx=xx*xx;
      double yy=y*y;
      double yyyy=yy*yy;
      double zz=z*z;
      double zzzz=zz*zz;
      double m = xx+yy+zz;
      double mm = m*m;
      double mmm = mm*m;
      double n = vigra::sqrt(m);
      double sq2pi = vigra::sqrt(2*M_PI);
      double r_n4 = pow4(r_-n);      
      double ss = s_*s_;
      double sss=ss*s_;

      return
    -((x*pow3(-r_ + n)*(xx*sqrt(pow3(m)) + r_*(yyyy+zzzz-xx*(yy + zz))))
      /(exp(vigra::sqrt(((xxxx + yyyy + zzzz)*r_n4)/mm)/(2*ss))*sq2pi*sss*mmm*
        vigra::sqrt(((xxxx + yyyy + zzzz)*r_n4)/mm)));
    }

    value_type dy(double x, double y, double z) const
    {  

      x = x - c_[0];
      y = y - c_[1];
      z = z - c_[2];

      double xx=x*x;
      double xxxx=xx*xx;
      double yy=y*y;
      double yyyy=yy*yy;
      double zz=z*z;
      double zzzz=zz*zz;
      double m = xx+yy+zz;
      double mm = m*m;
      double mmm = mm*m;
      double n = vigra::sqrt(m);
      double sq2pi = vigra::sqrt(2*M_PI);
      double r_n4 = pow4(r_-n);      
      double ss = s_*s_;
      double sss=ss*s_;

      return
    -((y*pow3(-r_ + n)*(yy*sqrt(pow3(m)) + r_*(xxxx+xx*yy-yy*zz + zzzz)))
      /(exp(vigra::sqrt(((xxxx + yyyy + zzzz)*r_n4)/mm)/(2*ss))*sq2pi*sss*mmm*
        vigra::sqrt(((xxxx + yyyy + zzzz)*r_n4)/mm)));
    }

    value_type dz(double x, double y, double z) const
    {  
      x = x - c_[0];
      y = y - c_[1];
      z = z - c_[2];

      double xx=x*x;
      double xxxx=xx*xx;
      double yy=y*y;
      double yyyy=yy*yy;
      double zz=z*z;
      double zzzz=zz*zz;
      double m = xx+yy+zz;
      double mm = m*m;
      double mmm = mm*m;
      double n = vigra::sqrt(m);
      double sq2pi = vigra::sqrt(2*M_PI);
      double r_n4 = pow4(r_-n);      
      double ss = s_*s_;
      double sss=ss*s_;

      return
    -((z*pow3(-r_ + n)*(zz*sqrt(pow3(m)) + r_*(xxxx+yyyy-xx*zz - yy*zz)))
      /(exp(vigra::sqrt(((xxxx + yyyy + zzzz)*r_n4)/mm)/(2*ss))*sq2pi*sss*mmm*
        vigra::sqrt(((xxxx + yyyy + zzzz)*r_n4)/mm)));
    }

    unsigned int width() const
    { return w_; }
    
    unsigned int height() const
    { return h_; }
    
    unsigned int depth() const
    { return d_; }

    bool isInsideX(double x) const
    {
      return x >= 0.0 && x <= width()-1.0;
    }
    bool isInsideY(double y) const
    {
      return y >= 0.0 && y <= height()-1.0;
    }
    
    bool isInsideZ(double z) const
    {
      return z >= 0.0 && z <= depth()-1.0;
    }

    bool isInside(double x, double y, double z) const
    {
      return isInsideX(x)  && isInsideY(y) && isInsideZ(z);
    }
    
    bool isInside(difference_type const & p) const
    {
      return isInside(p[0],p[1],p[2]);
    }
  };

  template<typename VOLUME>
  struct G2VolumeView 
  {
    typedef typename VOLUME::value_type value_type;
    typedef typename VOLUME::size_type size_type;
    typedef typename VOLUME::difference_type difference_type;
    
    const VOLUME vol;

    explicit G2VolumeView(VOLUME vol_) : vol(vol_) {}

    difference_type shape() const
    { return vol.shape(); }

    value_type operator()(double x, double y, double z) const
    { return vol.g(x,y,z); }
    
    value_type operator()(difference_type const & d) const
    { return operator()(d[0],d[1],d[2]); }

    value_type operator[](difference_type const & d) const
    { return operator()(d); }

    difference_type d(difference_type const & d) const
    {  return difference_type(dx(d),dy(d),dz(d)); }

    value_type dx(difference_type const & d) const
    { return dx(d[0],d[1],d[2]); }

    value_type dy(difference_type const & d) const
    { return dy(d[0],d[1],d[2]); }

    value_type dz(difference_type const & d) const
    { return dz(d[0],d[1],d[2]); }

    value_type dx(double x, double y, double z) const
    { return vol.gx(x,y,z); }

    value_type dy(double x, double y, double z) const
    { return vol.gy(x,y,z); }

    value_type dz(double x, double y, double z) const
    { return vol.gz(x,y,z); }

    unsigned int width() const
    { return vol.width(); }
    
    unsigned int height() const
    { return vol.height(); }
    
    unsigned int depth() const
    { return vol.depth(); }

    bool isInsideX(double x) const
    { return vol.isInsideX(x); }
    
    bool isInsideY(double y) const
    { return vol.isInsideY(y); }

    bool isInsideZ(double z) const
    { return vol.isInsideZ(z); }

    bool isInside(double x, double y, double z) const
    { return vol.isInside(x,y,z); }
    
    bool isInside(difference_type const & p) const
    { return vol.isInside(p); }
  };
}
  

#endif /* VIGRA_ANALYTICVOLUMEVIEW_HXX */
