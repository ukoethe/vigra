#include <utility>
#include <algorithm>

#include <vigra/tuple.hxx>
#include <vigra/mathutil.hxx>

#include "linesearch.hxx"

namespace vigra {

  template<typename FUNCTOR>
  class LineSearchMore : public LineSearch<FUNCTOR> 
  {
  public:
    typedef FUNCTOR Functor;
    typedef LineSearch<FUNCTOR> base_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::real_type real_type;
  
    void init(const difference_type & p) {}

    LineSearchMore(const Functor & func)
      : base_type(func),
    a_init  (  1e-6),
    a_min   (  1e-8),
    a_max   (  sqrt(3.)/ 2.),
    a_tol   (  1e-10), 
    phi_min ( -1e3),//not used
    mu      (  0.00001),
    eta     (  0.1)

    {};    
    
  private:
    typedef pair<real_type,real_type> Pair;
    struct Interval { Pair a,phi,phid; };
    struct Alpha { real_type a,phi,phid; };
    typedef triple<real_type,Interval,bool> Triple;
    
    real_type func(difference_type);
    real_type dfunc(difference_type);

    real_type a_init;
    real_type a_min;
    real_type a_max;
    real_type a_tol; 
    real_type phi_min;
    real_type mu;
    real_type eta;

    real_type cubic(real_type a, real_type b, 
            real_type fa, real_type fb, 
            real_type ga, real_type gb)
    {
      real_type beta1,beta2;
      beta1 = ga + gb - 3.0*(fa - fb)/(a - b);
      if(a < b)
        beta2 = sqrt(beta1*beta1 - ga*gb);
      else
        beta2 = -sqrt(beta1*beta1 - ga*gb);

      real_type denom = gb - ga + 2.0*beta2;
      if(denom == 0.0)
    return -1e99;
      else
        return b - (b - a)*(gb + beta2 - beta1)/denom;
    }

    real_type quadratic(real_type a, real_type b, 
            real_type fa, real_type fb, 
            real_type ga)
    {
      real_type denom = fa - fb - (a - b)*ga;
      if(denom == 0.0)
        return 1e99;
      else
        return a + 0.5 * ga / denom;
    }
    
    real_type secant(real_type a, real_type b, 
             real_type ga, real_type gb)
    {
      real_type denom = ga - gb;
      if(denom == 0.0)
        return 1e99;
      else
        return (b*ga - a*gb) / denom;
    }

    real_type cubic_ext(real_type a, real_type b, 
            real_type fa, real_type fb, 
            real_type ga, real_type gb)
    {
      using std::max;

      real_type beta1,beta2,alpha;
      beta1 = ga + gb - 3.0*(fa - fb)/(a - b);
      if(a < b)
        beta2 = sqrt(max(0.0, beta1*beta1 - ga*gb));
      else
        beta2 = -sqrt(max(0.0, beta1*beta1 - ga*gb));

      real_type denom = gb - ga + 2.0*beta2;
      if(denom == 0.0)
        alpha = -1e99;
      else
        alpha = b - (b - a)*(gb + beta2 - beta1)/denom;

      return alpha;
    }

    Triple update(Alpha a, Interval Ik, real_type at, real_type al, real_type au, 
          real_type ft, real_type fl, real_type fu, real_type gt, 
          real_type gl, real_type gu, bool bracketed, Pair Ik_lim, 
          real_type d = 0.66) 
    {
      using std::max;
      using std::min;

      // Trial value selection.
      // Case 1.
      real_type ac,aq,as,at_new;

      if(ft > fl){
    // The minimum is bracketed.
    bracketed = 1;
    // Interpolation.
    ac = cubic(al, at, fl, ft, gl, gt);
    aq = quadratic(al, at, fl, ft, gl);
    // Return at+.
    if(abs(ac - al) < abs(aq - al))
      at_new = ac;
    else
      at_new = 0.5*(aq + ac);
      }// Case 2.
      else if((gt * gl) < 0.0) {
    // The minimum is bracketed.
    bracketed = 1;
    // Interpolation.
    ac = cubic(al, at, fl, ft, gl, gt);
    as = secant(al, at, gl, gt);
    // Return at+.
    if(abs(ac - at) >= abs(as - at))
      at_new = ac;
    else
      at_new = as;
      }
      // Case 3.
      else if(abs(gt) <= abs(gl)){
    // Interpolation.
    ac = cubic_ext(al, at, fl, ft, gl, gt);
    //beta1 = ret.second; //remove?
        //beta2 = ret.third;
    
    if(at > al)
        // Set ac to the upper limit.
      ac = Ik_lim.second;
    else
      // Set ac to the lower limit.
      ac = Ik_lim.first;
    
    as = secant(al, at, gl, gt);
    
    // Test if(bracketed.
    if(bracketed) {
      if(abs(ac - at) < abs(as - at))
        at_new = ac;
      else
        at_new = as;
      // Redefine at+.
      if(at > al)
        at_new = min(at + d*(au - at), at_new);
      else
        at_new = max(at + d*(au - at), at_new);
    }else{
      if(abs(ac - at) > abs(as - at))
        at_new = ac;
      else
        at_new = as;
      
      // Check limits.
      if(at_new < Ik_lim.first)
        at_new = Ik_lim.first;
      if(at_new > Ik_lim.second)
        at_new = Ik_lim.second;
    }
      }
      // Case 4.
      else {
    if(bracketed)
      at_new = cubic(au, at, fu, ft, gu, gt);
    else if(at > al)
      at_new = Ik_lim.second;
    else
      at_new = Ik_lim.first;
      }
  
      // Interval updating algorithm.
      Interval Ik_new = Ik;

      if(ft > fl) {
    Ik_new.a.second = at;
    Ik_new.phi.second = a.phi;
    Ik_new.phid.second = a.phid;
      }
      else if(gt*(al - at) > 0.0){
    Ik_new.a.first = at;
    Ik_new.phi.first = a.phi;
    Ik_new.phid.first = a.phid;
      }
      else{
    Ik_new.a.first = at;
    Ik_new.phi.first = a.phi;
    Ik_new.phid.first = a.phid;
    Ik_new.a.second = al;
    Ik_new.phi.second = Ik.phi.first;
    Ik_new.phid.second = Ik.phid.first;
      }
      return Triple(at_new, Ik_new, bracketed);
    }
    
    real_type impl()
    {
      using std::max;
      using std::min;
      
      //Initialise values.
      int k          = 0;
      bool mod_flag  = 1;
      bool bracketed = 0;
    
      Alpha a0;
      a0.a = 0.0;
      a0.phid = base_type::df1dim(0.0);
      if(a0.phid > 0.0) {
// 	std::cerr << base_type::p << std::endl;
// 	std::cerr <<  << std::endl;
// 	std::cerr << base_type::func.d(base_type::p) << std::endl;
// 	std::cerr << a0.phid << std::endl;
    base_type::dir *= -1;
// 	vigra_precondition(false, "The gradient at point 0 of this line search is positive,"
// 			   "ie p is not a descent direction && the line search will not work.");
      }      
      a0.phid = base_type::df1dim(0.0);
      a0.phi = base_type::f1dim(0.0);
    
      if(a_max == -1)
    a_max = 4.*max(1.0,a_init);
    
      Pair Ik_lim(0.0, 5.0*a_init);
    
      real_type width = a_max - a_min;
      real_type width2 = 2.0*width;

      //Initialise sequence data.
      Alpha a;
      a.a    = a_init;
      a.phi  = f1dim(a.a);
      a.phid = df1dim(a.a);
  
      // Initialise Interval data.
      Interval Ik;
      Ik.a    = Pair(0.0, 0.0);
      Ik.phi  = Pair(a0.phi, a0.phi);
      Ik.phid = Pair(a0.phid, a0.phid);


      // Test for errors.

      vigra_precondition(a.a >= a_min,"Alpha is less than Alpha_min");
      vigra_precondition(a.a <= a_max,"Alpha is greater than Alpha_max");
  
      while(true) {
    real_type curv,suff_dec;
    // Test values.
    curv = mu * a0.phid;
    suff_dec = a0.phi + a.a * curv;

    // Modification flag, 0 - phi, 1 - psi.
    if(mod_flag)
      if(a.phi <= suff_dec && a.phid >= 0.0)
        mod_flag = 0;

    // Test for convergence using the strong Wolfe conditions.
    if(a.phi <= suff_dec && abs(a.phid) <= eta * abs(a0.phid))
      return a.a;

    // Test if(limits have been reached.
    if(a.a == a_min)
      if(a.phi > suff_dec || a.phid >= curv)
        return a.a;
    if(a.a == a_max)
      if(a.phi <= suff_dec && a.phid <= curv)
        return a.a;

    if(bracketed) {
      // Test f|| roundoff err||.
      if(a.a <= Ik_lim.first || a.a >= Ik_lim.second)
        return a.a;
      // Test to see if(a_tol has been reached.
      if(Ik_lim.second - Ik_lim.first <= a_tol * Ik_lim.second)
        return a.a;
    }
    
    // Choose a safeguarded ak in set Ik which is a subset of [a_min, a_max], && update the Interval Ik.
    Alpha a_new;
    Interval Ik_new;
    if(mod_flag && a.phi <= Ik.phi.first && a.phi > suff_dec){
      // Calculate the modified function values && gradients at at, al, && au.
      real_type psi,psi_l,psi_u,psid,psi_ld,psi_ud;
      psi = a.phi - curv * a.a;
      psi_l = Ik.phi.first - curv * Ik.a.first;
      psi_u = Ik.phi.second - curv * Ik.a.second;
      psid = a.phid - curv;
      psi_ld = Ik.phid.first - curv;
      psi_ud = Ik.phid.second - curv;
      
      Triple ret = update(a, Ik, a.a, Ik.a.first, Ik.a.second, 
                  psi, psi_l, psi_u, psid, psi_ld, psi_ud, bracketed, Ik_lim);
      a_new.a    = ret.first;
      Ik_new     = ret.second;
      bracketed  = ret.third;
    }
    else {
      Triple ret = update(a, Ik, a.a, Ik.a.first, Ik.a.second, a.phi, Ik.phi.first, Ik.phi.second, a.phid, Ik.phid.first, Ik.phid.second, bracketed, Ik_lim);
      a_new.a    = ret.first;
      Ik_new     = ret.second;
      bracketed  = ret.third;
    }
    
    // Bisection step.
    if(bracketed){
      real_type size = abs(Ik_new.a.first - Ik_new.a.second);
      if(size >= 0.66 * width2)
        a_new.a = 0.5 * (Ik_new.a.first + Ik_new.a.second);
      width2 = width;
      width = size;
    }
    // Limit.
    if(bracketed) {
      Ik_lim.first = min(Ik_new.a.first, Ik_new.a.second);
      Ik_lim.second = max(Ik_new.a.first, Ik_new.a.second);
    }
    else {
      Ik_lim.first = a_new.a + 1.1 * (a_new.a - Ik_new.a.first);
      Ik_lim.second = a_new.a + 4.0 * (a_new.a - Ik_new.a.first);
    }
    
    if(bracketed)
      if(a_new.a <= Ik_lim.first || a_new.a >= Ik_lim.second || Ik_lim.second - Ik_lim.first <= a_tol * Ik_lim.second)
        a_new.a = Ik.a.first;
    
    // The step must be between a_min && a_max.
    if(a_new.a < a_min)
      a_new.a = a_min;
    if(a_new.a > a_max)
      a_new.a = a_max;
    
    // Calculate new values.
    a_new.phi = f1dim(a_new.a);
    a_new.phid = df1dim(a_new.a);
    
    // Shift data from k+1 to k.
    k = k + 1;
    a = a_new;
    Ik = Ik_new;
      }
    }    
  };
} // namespace vigra
