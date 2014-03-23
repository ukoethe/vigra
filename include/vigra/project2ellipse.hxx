
/************************************************************************/
/*                                                                      */
/*                 Copyright 2012 by Frank Lenzen &                     */
/*                                           Ullrich Koethe             */
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

#ifndef VIGRA_PROJECT2ELLIPSE_HXX
#define VIGRA_PROJECT2ELLIPSE_HXX
#include <iostream>

namespace vigra{

  namespace detail{

//---------------------------------------------------------------------------

inline void projectEllipse2D_FirstQuad (double &vx, double &vy, double a, double b, const double eps, const int iter_max){                
  
  double t=0,tmin,tmax,d1,d2,f,x1,y1,l;                      
  int i;
  tmax=std::max(2*a*vx,2*b*vy);
  tmin=-b*b;
  d1=a*vx/(tmax+a*a);
  d2=b*vy/(tmax+b*b);
  f=d1*d1+d2*d2-1;
    
  for (i=0;i<iter_max;i++){
  
    t=.5*(tmin+tmax);
    d1=a*vx/(t+a*a);
    d2=b*vy/(t+b*b);
    f=d1*d1+d2*d2-1;
    x1=a*vx/(t+a*a);
    y1=b*vy/(t+b*b);
    l=x1*x1+y1*y1-1;
   
    if (fabs(l)<eps)
      break;
    if(f>0)
      tmin=t;
    else if(f<0)
      tmax=t;
    else
      break;
  }
  d1=vx;
  d2=vy;
  vx=a*a*vx/(t+a*a);
  vy=b*b*vy/(t+b*b);
  d1 = vx - d1;
  d2 = vy - d2;
  return;
}


inline void projectEllipse2D(double &vx, double &vy, const double _a,const double _b,const double eps,const int max_iter){
  
  //double err;
  double a=_a,b=_b;
  
  //check if inside ellipse
  if (((vx/a)*(vx/a)+(vy/b)*(vy/b))<=1){
    return;
  }
  
  // special case of a circle
  if (fabs(a-b) < eps){
    double l = sqrt(vx*vx+vy*vy);
    if (l>(a+b)/2.){
      vx=(a+b)/(2*l)*vx;
      vy=(a+b)/(2*l)*vy;
    }
    return;
  }
  
  // reflect vx -> -vx, if necessary
  bool x_reflect;
  if (vx > eps){
    x_reflect = false;
  }
  else if (vx < -eps){
    x_reflect = true;
    vx = -vx;
  }
  else{
    x_reflect = false;
    vx = 0.0;
  }
  // reflect vy -> vy = -V if necessary
  bool y_reflect;
  if (vy > eps){
    y_reflect = false;
  }
  else if (vy < -eps){
    y_reflect = true;
    vy = -vy;
  }
  else{
    y_reflect = false;
    vy = 0.0;
  }
  
  // swap axes if necessary
  bool swapped;
  if (a >= b){
    swapped = false;
  }
  else{
    swapped = true;
    std::swap(a,b);
    std::swap(vx,vy);
  }
  if (vx != 0.0){
    if (vy != 0.0){
      projectEllipse2D_FirstQuad(vx,vy,a,b,eps,max_iter);
    }
    else{
      if (vx < a - b*b/a){
        double vx_temp = vx;
        vx = a*a*vx/(a*a-b*b);
        vy = vy*sqrt(fabs(1.0-vx_temp/a*vx_temp/a));
      }
      else{
        vx = a;
        vy = 0.0;
      }
    }
  }
  else{
    vx = 0.0;
    vy = b;
  }
  if (swapped){
    std::swap(vx,vy);
  }
  if (y_reflect){
    vy = -vy;
  }
  if (x_reflect){
    vx = -vx;
  }
  return;
}
  } //closing namespace detail
} //closing namespace vigra
#endif // VIGRA_PROJECT2ELLIPSE_HXX
