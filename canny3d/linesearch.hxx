#ifndef LINESEARCH_HXX
#define LINESEARCH_HXX

#include <iostream>
#include <limits>
#include <vigra/mathutil.hxx>

namespace vigra {

  namespace detail {
    template<class T>
    inline void swap(T &a, T &b)
    {T dum=a; a=b; b=dum;}
  }
    
  template<typename FUNCTOR>
  class LineSearch {
    
  public:
    typedef FUNCTOR Functor;
    typedef typename Functor::difference_type difference_type;
    typedef typename difference_type::value_type value_type;
    typedef typename NumericTraits<value_type>::RealPromote real_type;

    LineSearch(const Functor & arg_func)
      : func(arg_func)
    {};    

    virtual void init(const difference_type & p) = 0;
    real_type operator()(const difference_type & arg_p, const difference_type & arg_dir) 
    { p = arg_p; dir = arg_dir; return impl(); }
    
  protected:
    Functor func;    
    difference_type p;
    difference_type dir;

    virtual real_type impl() = 0;

    real_type f1dim(real_type x) const {
      difference_type v = p + x*dir;
      if(!func.isInside(v)) {
    throw std::runtime_error("value not in range");
      }
      return func(v);
    }

    real_type df1dim(real_type x) const {
      difference_type v = p + x*dir;
      if(!func.isInside(v)) {
    throw std::runtime_error("value not in range");
      }
      return dot(func.d(v),dir);
    }            
  };

  template<typename FUNCTOR>
  class LineSearchBrend : public LineSearch<FUNCTOR> 
  {
  public:
    typedef FUNCTOR Functor;
    typedef LineSearch<FUNCTOR> base_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::real_type real_type;

    LineSearchBrend(const Functor & func, real_type mlimit)
      : base_type(func),
    ax_(0.0),
    xx_(0.125),
    tol(2.0e-6),
    mlimit(mlimit),
    itmax(400)
    {};    

    LineSearchBrend(const Functor & func)
      : base_type(func),
    ax_(0.0),
    xx_(0.125),
    tol(2.0e-10),
    mlimit(100.0),
    itmax(400)
    {};    

    void init(const difference_type & p) {};
    
  private:
    real_type impl();
    void mnbrak(real_type &, real_type &, real_type &, real_type &, real_type &, real_type &);
    real_type dbrent(real_type, real_type, real_type, real_type &);
    
    void shft3(real_type &a, real_type &b, real_type &c, 
           const real_type d) const
    { a=b; b=c; c=d; }
    
    void mov3(real_type &a, real_type &b, real_type &c, 
          real_type d, real_type e, real_type f) const
    { a=d; b=e; c=f; }
    
    real_type ax_;
    real_type xx_;
    real_type tol;
    real_type mlimit;
    int itmax;
  };

  template<typename FUNCTOR>
  typename LineSearchBrend<FUNCTOR>::real_type
  LineSearchBrend<FUNCTOR>::impl()
  {
    int j;
    real_type xmin,fx,fb,fa,bx,ax,xx;

    ax=ax_;
    xx=xx_;
    mnbrak(ax,xx,bx,fa,fx,fb);
    dbrent(ax,xx,bx,xmin);
    //dbrent(0.0,0.25,0.5,xmin);
    return xmin;
  }

  /**
   * find inital bracketting of a functions local minimum between ax and
   * bx. Taken from Numerical Recipes (2001)
   */ 
  template<typename FUNCTOR>
  void 
  LineSearchBrend<FUNCTOR>::mnbrak
  (LineSearchBrend<FUNCTOR>::real_type &ax,   
   LineSearchBrend<FUNCTOR>::real_type &bx, 
   LineSearchBrend<FUNCTOR>::real_type &cx, 
   LineSearchBrend<FUNCTOR>::real_type &fa, 
   LineSearchBrend<FUNCTOR>::real_type &fb,   
   LineSearchBrend<FUNCTOR>::real_type &fc)
  {
    const real_type GOLD=1.618034, GLIMIT=mlimit,TINY=1.0e-20;
    real_type ulim,u,r,q,fu;
    
    fa=f1dim(ax);
    fb=f1dim(bx);
    if (fb > fa) {
      detail::swap(ax,bx);
      detail::swap(fb,fa);
    }
    cx=bx+GOLD*(bx-ax);
    fc=f1dim(cx);
    while (fb > fc) {
      r=(bx-ax)*(fb-fc);
      q=(bx-cx)*(fb-fa);
      u=bx-((bx-cx)*q-(bx-ax)*r)/
    (2.0*sign(VIGRA_CSTD::max(fabs(q-r),TINY),q-r));
      ulim=bx+GLIMIT*(cx-bx);
      if ((bx-u)*(u-cx) > 0.0) {
    fu=f1dim(u);
    if (fu < fc) {
      ax=bx;
      bx=u;
      fa=fb;
      fb=fu;
      return;
    } else if (fu > fb) {
      cx=u;
      fc=fu;
      return;
    }
    u=cx+GOLD*(cx-bx);
    fu=f1dim(u);
      } else if ((cx-u)*(u-ulim) > 0.0) {
    fu=f1dim(u);
    if (fu < fc) {
      shft3(bx,cx,u,cx+GOLD*(cx-bx));
      shft3(fb,fc,fu,f1dim(u)); //???
    }
      } else if ((u-ulim)*(ulim-cx) >= 0.0) {
    u=ulim;
    fu=f1dim(u);
      } else {
    u=cx+GOLD*(cx-bx);
    fu=f1dim(u);
      }
      shft3(ax,bx,cx,u);
      shft3(fa,fb,fc,fu);
    }
  }

  template<typename FUNCTOR>
  typename LineSearchBrend<FUNCTOR>::real_type
  LineSearchBrend<FUNCTOR>::dbrent
  (const LineSearchBrend<FUNCTOR>::real_type ax, 
   const LineSearchBrend<FUNCTOR>::real_type bx, 
   const LineSearchBrend<FUNCTOR>::real_type cx, 
   LineSearchBrend<FUNCTOR>::real_type &xmin)
  {
    const real_type ZEPS=std::numeric_limits<real_type>::epsilon()*1.0e-3;
    bool ok1,ok2;
    int iter;
    real_type a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
    real_type fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;
    fw=fv=fx=f1dim(x);
    dw=dv=dx=df1dim(x);
    for (iter=0;iter<itmax;iter++) {
      xm=0.5*(a+b);
      tol1=tol*fabs(x)+ZEPS;
      tol2=2.0*tol1;
      if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
    xmin=x;
    return fx;
      }
      if (fabs(e) > tol1) {
    d1=2.0*(b-a);
    d2=d1;
    if (dw != dx) d1=(w-x)*dx/(dx-dw);
    if (dv != dx) d2=(v-x)*dx/(dx-dv);
    u1=x+d1;
    u2=x+d2;
    ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
    ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
    olde=e;
    e=d;
    if (ok1 || ok2) {
      if (ok1 && ok2)
        d=(fabs(d1) < fabs(d2) ? d1 : d2);
      else if (ok1)
        d=d1;
      else
        d=d2;
      if (fabs(d) <= fabs(0.5*olde)) {
        u=x+d;
        if (u-a < tol2 || b-u < tol2)
          d=sign(tol1,xm-x);
      } else {
        d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
    } else {
      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
    }
      } else {
    d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
      if (fabs(d) >= tol1) {
    u=x+d;
    fu=f1dim(u);
      } else {
    u=x+sign(tol1,d);
    fu=f1dim(u);
    if (fu > fx) {
      xmin=x;
      return fx;
    }
      }
      du=df1dim(u);
      if (fu <= fx) {
    if (u >= x) a=x; else b=x;
    mov3(v,fv,dv,w,fw,dw);
    mov3(w,fw,dw,x,fx,dx);
    mov3(x,fx,dx,u,fu,du);
      } else {
    if (u < x) a=u; else b=u;
    if (fu <= fw || w == x) {
      mov3(v,fv,dv,w,fw,dw);
      mov3(w,fw,dw,u,fu,du);
    } else if (fu < fv || v == x || v == w) {
      mov3(v,fv,dv,u,fu,du);
    }
      }
    }
    return 0.0;
  }

  template<typename FUNCTOR>
  struct LineSearchDNew : public LineSearch<FUNCTOR> {

    typedef FUNCTOR Functor;
    typedef LineSearch<FUNCTOR> base_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::real_type real_type;
    
    real_type tol;
    int itmax;
    real_type stpmax;
    bool check;

    LineSearchDNew(const Functor & func)
      : base_type(func),
    tol(2.0e-8),
    itmax(200)
    {};    

    enum {N = difference_type::static_size};

    void init(const difference_type & p) 
    { stpmax = 100.0 * VIGRA_CSTD::max(norm(p),real_type(N)); };

    real_type impl();
  };

  //p.388
  template<typename FUNCTOR>
  typename LineSearchDNew<FUNCTOR>::real_type
  LineSearchDNew<FUNCTOR>::impl()
  {
    real_type f; //FIX
    difference_type x;//FIX
    difference_type g=base_type::func.d(base_type::p);//FIX
    real_type fold=base_type::func(base_type::p);

    const real_type ALF=1.0e-4, TOLX=std::numeric_limits<real_type>::epsilon();
    int i;
    real_type a,alpha,alpha2=0.0,alphamin,b,disc,f2=0.0;
    real_type rhs1,rhs2,slope,len,temp,test,tmplam;

    check = false;
    len = norm(base_type::dir);
    if (len > stpmax)
      base_type::dir *= stpmax/len;
    slope=dot(g,base_type::dir);//???
    if (slope >= 0.0) {
      vigra_fail("round off problems");
    }

    test=0.0;
    for(i = 0;i < N;++i) {
      temp = abs(base_type::dir[i])/VIGRA_CSTD::max(abs(base_type::p[i]),1.0);
      if (temp > test) test=temp;
    }

    alphamin = TOLX/test;
    alpha = 1.0;
    for(;;) {
      x = base_type::p + alpha*base_type::dir;
      f = base_type::func(x);
      if(alpha < alphamin) {
    x=base_type::p;//fix interface
    check=true;
    return alpha;
      } 
      else if(f <= fold+ALF*alpha*slope) 
    return alpha;
      else {
    if (alpha == 1.0)
      tmplam = -slope/(2.0*(f-fold-slope));
    else {
      rhs1 = f-fold-alpha*slope;
      rhs2 = f2-fold-alpha2*slope;
      a = (rhs1/(alpha*alpha) - rhs2/(alpha2*alpha2))/(alpha-alpha2);
      b =(-alpha2*rhs1/(alpha*alpha) + alpha*rhs2/(alpha2*alpha2))/(alpha-alpha2);
      if (a == 0.0) 
        tmplam = -slope/(2.0*b);
      else {
        disc=b*b-3.0*a*slope;
        if(disc < 0.0) 
          tmplam=0.5*alpha;
        else if(b <= 0.0) 
          tmplam=(-b+VIGRA_CSTD::sqrt(disc))/(3.0*a);
        else 
          tmplam=-slope/(b+VIGRA_CSTD::sqrt(disc));
      }
      if (tmplam > 0.5*alpha)
        tmplam=0.5*alpha;
    }
      }
      alpha2=alpha;
      f2 = f;
      alpha=VIGRA_CSTD::max(tmplam,0.1*alpha);
    }
  }

  template<typename FUNCTOR>
  struct LineSearchDNew2 : public LineSearchDNew<FUNCTOR> {

    typedef FUNCTOR Functor;
    typedef LineSearchDNew<FUNCTOR> base_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::real_type real_type;
    
    LineSearchDNew2(const Functor & func)
      : base_type(func)
    {};

    difference_type p0;
    difference_type g0;

    void init(const difference_type & p) {
      p0 = p;
      g0 = -base_type::func.d(p);
      g0/= norm(g0);
      //      base_type::init(p);
    }
    
    void computeStpmax()
    { 
      difference_type xi = base_type::p + base_type::dir - p0;
      real_type px = dot(xi,g0);
      difference_type pxg = px*g0;
      difference_type pp = xi - pxg;
      base_type::stpmax = norm(xi)/norm(pxg + (pp/norm(pp)));
      std::cerr << base_type::stpmax << std::endl;
      std::cerr << 1.0/norm(base_type::dir) << std::endl;
//       std::cerr << g0 << std::endl;
//       std::cerr << base_type::dir << std::endl;
//       std::cerr << base_type::p << std::endl;
//       std::cerr << p0 << std::endl;

    };
    
    real_type impl()
    {
      computeStpmax();
      return base_type::impl();
    }
  };
}
#endif // LINESEARCH_HXX
