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

#ifndef VIGRA_TV_FILTER_HXX
#define VIGRA_TV_FILTER_HXX

#include <iostream>
#include <cmath>
#include "config.hxx"
#include "impex.hxx"
#include "separableconvolution.hxx"
#include "multi_array.hxx"
#include "multi_math.hxx"
#include "eigensystem.hxx"
#include "convolution.hxx"
#include "fixedpoint.hxx"
#include "project2ellipse.hxx"

#ifndef VIGRA_MIXED_2ND_DERIVATIVES
#define VIGRA_MIXED_2ND_DERIVATIVES 1
#endif

#define setZeroX(A) A.subarray(Shape2(width-1,0),Shape2(width,height))*=0;
#define setZeroY(A) A.subarray(Shape2(0,height-1),Shape2(width,height))*=0;

namespace vigra{



/** \addtogroup NonLinearDiffusion
*/

//@{

/********************************************************/
/*                                                      */
/*           totalVariationFilter                       */
/*                                                      */
/********************************************************/

/** \brief Performs standard Total Variation Regularization

The algorithm minimizes

\f[
       \min_u \int_\Omega \frac{1}{2} (u-f)^2\;dx + \alpha TV(u)\qquad\qquad (1)
\f]
where <em>\f$ f=f(x)\f$</em> are the two dimensional noisy data,
<em> \f$ u=u(x)\f$</em> are the smoothed data,<em>\f$ \alpha \ge 0 \f$</em>
is the filter parameter and <em>\f$ TV(u)\f$ </em> is the total variation semi-norm.

<b> Declarations:</b>

\code
namespace vigra {
      template <class stride1,class stride2>
      void totalVariationFilter(MultiArrayView<2,double,stride1> data,
                                MultiArrayView<2,double,stride2> out,
                                double alpha,
                                int steps,
                                double eps=0);
      void totalVariationFilter(MultiArrayView<2,double,stride1> data,
                                MultiArrayView<2,double,stride2> weight,
                                MultiArrayView<2,double,stride3> out,
                                double alpha,
                                int steps,
                                double eps=0);
}
\endcode

\ref totalVariationFilter() implements a primal-dual algorithm to solve (1).

Input:
     <table>
     <tr><td><em>data</em>:  </td><td> input data to be smoothed. </td></tr>
     <tr><td><em>alpha</em>: </td><td> smoothing parameter.</td></tr>
     <tr><td><em>steps</em>: </td><td> maximal number of iteration steps. </td></tr>
     <tr><td><em>eps</em>:   </td><td> The algorithm stops, if the primal-dual gap is below the threshold <em>eps</em>.
     </table>

     Output:

     <em>out</em> contains the filtered data.

     In addition, a point-wise weight (\f$ \ge 0 \f$) for the data term can be provided (overloaded function).

    <b> Usage:</b>

    <b>\#include</b> \<vigra/tv_filter.hxx\>

    \code
    MultiArray<2,double> data(Shape2(width,height));  // to be initialized
    MultiArray<2,double> out(Shape2(width,height));
    MultiArray<2,double> weight(Shape2(width,height))); //optional argument in overloaded function, to be initialized if used
    int steps;        // to be initialized
    double alpha,eps; // to be initialized

    totalVariationFilter(data,out,alpha,steps,eps);
    \endcode
    or
    \code
    totalVariationFilter(data,weight,out,alpha,steps,eps);
    \endcode

 */
doxygen_overloaded_function(template <...> void totalVariationFilter)

template <class stride1,class stride2>
void totalVariationFilter(MultiArrayView<2,double,stride1> data,MultiArrayView<2,double,stride2> out, double alpha, int steps, double eps=0){

  using namespace multi_math;
  int width=data.shape(0),height=data.shape(1);

  MultiArray<2,double> temp1(data.shape()),temp2(data.shape()),vx(data.shape()),vy(data.shape()),u_bar(data.shape());
  Kernel1D<double> Lx,LTx;
  Lx.initExplicitly(-1,0)=1,-1;                       // = Right sided finite differences for d/dx and d/dy
  Lx.setBorderTreatment(BORDER_TREATMENT_REFLECT);   //   with hom. Neumann boundary conditions
  LTx.initExplicitly(0,1)=-1,1;                     //  = Left sided finite differences for -d/dx  and -d/dy
  LTx.setBorderTreatment(BORDER_TREATMENT_ZEROPAD);  //   with hom. Dirichlet b. c.

  out=data;
  u_bar=data;

  double tau=1.0 / std::max(alpha,1.) / std::sqrt(8.0) * 0.06;
  double sigma=1.0 / std::sqrt(8.0) / 0.06;

  for (int i=0;i<steps;i++){

    separableConvolveX(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroX(temp1);
    vx+=(sigma*temp1);
    separableConvolveY(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroY(temp1);
    vy+=(sigma*temp1);

    //project to constraint set
    for (int y=0;y<data.shape(1);y++){
      for (int x=0;x<data.shape(0);x++){
        double l=hypot(vx(x,y),vy(x,y));
        if (l>1){
          vx(x,y)/=l;
          vy(x,y)/=l;
        }
      }
    }

    separableConvolveX(srcImageRange(vx),destImage(temp1),kernel1d(LTx));
    separableConvolveY(srcImageRange(vy),destImage(temp2),kernel1d(LTx));
    u_bar=out;
    out-=tau*(out-data+alpha*(temp1+temp2));
    u_bar=2*out-u_bar;   //cf. Chambolle/Pock and Popov's algorithm


    //stopping criterion
    if (eps>0){
      separableConvolveX(srcImageRange(out),destImage(temp1),kernel1d(Lx));setZeroX(temp1);
      separableConvolveY(srcImageRange(out),destImage(temp2),kernel1d(Lx));setZeroY(temp2);

      double f_primal=0,f_dual=0;
      for (int y=0;y<data.shape(1);y++){
        for (int x=0;x<data.shape(0);x++){
          f_primal+=.5*(out(x,y)-data(x,y))*(out(x,y)-data(x,y))+alpha*hypot(temp1(x,y),temp2(x,y));
        }
      }
      separableConvolveX(srcImageRange(vx),destImage(temp1),kernel1d(LTx));
      separableConvolveY(srcImageRange(vy),destImage(temp2),kernel1d(LTx));
      for (int y=0;y<data.shape(1);y++){
        for (int x=0;x<data.shape(0);x++){
          double divv=temp1(x,y)+temp2(x,y);
          f_dual+=-.5*alpha*alpha*(divv*divv)+alpha*data(x,y)*divv;
        }
      }
      if (f_primal>0 && (f_primal-f_dual)/f_primal<eps){
        break;
      }
    }
  }
}

template <class stride1,class stride2, class stride3>
void totalVariationFilter(MultiArrayView<2,double,stride1> data,MultiArrayView<2,double,stride2> weight, MultiArrayView<2,double,stride3> out,double alpha, int steps, double eps=0){

  using namespace multi_math;
  int width=data.shape(0),height=data.shape(1);

  MultiArray<2,double> temp1(data.shape()),temp2(data.shape()),vx(data.shape()),vy(data.shape()),u_bar(data.shape());
  Kernel1D<double> Lx,LTx;
  Lx.initExplicitly(-1,0)=1,-1;                       // = Right sided finite differences for d/dx and d/dy
  Lx.setBorderTreatment(BORDER_TREATMENT_REFLECT);   //   with hom. Neumann boundary conditions
  LTx.initExplicitly(0,1)=-1,1;                     //  = Left sided finite differences for -d/dx  and -d/dy
  LTx.setBorderTreatment(BORDER_TREATMENT_ZEROPAD);  //   with hom. Dirichlet b. c.

  out=data;
  u_bar=data;

  double tau=1.0 / std::max(alpha,1.) / std::sqrt(8.0) * 0.06;
  double sigma=1.0 / std::sqrt(8.0) / 0.06;

  for (int i=0;i<steps;i++){
    separableConvolveX(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroX(temp1);
    vx+=(sigma*temp1);
    separableConvolveY(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroY(temp1);
    vy+=(sigma*temp1);

    //project to constraint set
    for (int y=0;y<data.shape(1);y++){
      for (int x=0;x<data.shape(0);x++){
        double l=hypot(vx(x,y),vy(x,y));
        if (l>1){
          vx(x,y)/=l;
          vy(x,y)/=l;
        }
      }
    }

    separableConvolveX(srcImageRange(vx),destImage(temp1),kernel1d(LTx));
    separableConvolveY(srcImageRange(vy),destImage(temp2),kernel1d(LTx));
    u_bar=out;
    out-=tau*(weight*(out-data)+alpha*(temp1+temp2));
    u_bar=2*out-u_bar;


    //stopping criterion
    if (eps>0){
      separableConvolveX(srcImageRange(out),destImage(temp1),kernel1d(Lx));setZeroX(temp1);
      separableConvolveY(srcImageRange(out),destImage(temp2),kernel1d(Lx));setZeroY(temp2);

      double f_primal=0,f_dual=0;
      for (int y=0;y<data.shape(1);y++){
        for (int x=0;x<data.shape(0);x++){
          f_primal+=.5*weight(x,y)*(out(x,y)-data(x,y))*(out(x,y)-data(x,y))+alpha*hypot(temp1(x,y),temp2(x,y));
        }
      }
      separableConvolveX(srcImageRange(vx),destImage(temp1),kernel1d(LTx));
      separableConvolveY(srcImageRange(vy),destImage(temp2),kernel1d(LTx));
      for (int y=0;y<data.shape(1);y++){
        for (int x=0;x<data.shape(0);x++){
          double divv=temp1(x,y)+temp2(x,y);
          f_dual+=-.5*alpha*alpha*(weight(x,y)*divv*divv)+alpha*data(x,y)*divv;
        }
      }
      if (f_primal>0 && (f_primal-f_dual)/f_primal<eps){
        break;
      }
    }
  }
}
//<!--\f$ \alpha(x)=\beta(x)=\beta_{par}\f$ in homogeneous regions without edges,
//and \f$ \alpha(x)=\alpha_{par}\f$ at edges.-->


/********************************************************/
/*                                                      */
/*         getAnisotropy                                */
/*                                                      */
/********************************************************/

/** \brief Sets up directional data for anisotropic regularization

This routine provides a two-dimensional normalized vector field \f$ v \f$, which is normal to edges in the given data,
found as the eigenvector of the structure tensor belonging to the largest eigenvalue.
\f$ v \f$ is encoded by a scalar field \f$ \varphi \f$ of angles, i.e.
\f$ v(x)=(\cos(\varphi(x)),\sin(\varphi(x)))^\top \f$.

In addition, two scalar fields \f$ \alpha \f$ and  \f$ \beta \f$ are generated from
scalar parameters \f$ \alpha_{par}\f$ and \f$ \beta_{par}\f$, such that
<center>
<table>
<tr> <td>\f$ \alpha(x)= \alpha_{par}\f$ at edges,</td></tr>
<tr> <td>\f$ \alpha(x)= \beta_{par}\f$ in homogeneous regions,</td></tr>
<tr> <td>\f$ \beta(x)=\beta_{par}\f$ .</td></tr>
</table>
</center>

<b> Declarations:</b>

\code
namespace vigra {
void getAnisotropy(MultiArrayView<2,double,stride1> data,
                   MultiArrayView<2,double,stride2> phi,
                   MultiArrayView<2,double,stride3> alpha,
                   MultiArrayView<2,double,stride4> beta,
                   double alpha_par,
                   double beta_par,
                   double sigma_par,
                   double rho_par,
                   double K_par);
}
\endcode

Output:
<table>
  <tr><td>Three scalar fields <em>phi</em>, <em>alpha</em> and <em>beta</em>.</td></tr>
</table>

Input:
<table>
<tr><td><em>data</em>:</td><td>two-dimensional scalar field.</td></tr>
<tr><td><em>alpha_par,beta_par</em>:</td><td>two positive values for setting up the scalar fields alpha and beta</tr>
<tr><td><em>sigma_par</em>:</td><td> non-negative parameter for presmoothing the data.</td></tr>
<tr><td> <em>rho_par</em>:</td><td> non-negative parameter for presmoothing the structure tensor.</td></tr>
<tr><td><em>K_par</em>:</td><td> positive edge sensitivity parameter.</td></tr>
 </table>

(see \ref anisotropicTotalVariationFilter() and \ref secondOrderTotalVariationFilter() for usage in an application).
*/
doxygen_overloaded_function(template <...> void getAnisotropy)

template <class stride1,class stride2,class stride3,class stride4>
void getAnisotropy(MultiArrayView<2,double,stride1> data,MultiArrayView<2,double,stride2> phi,
                    MultiArrayView<2,double,stride3> alpha, MultiArrayView<2,double,stride4> beta,
                    double alpha_par, double beta_par, double sigma_par, double rho_par, double K_par){

  using namespace multi_math;
  int width=data.shape(0),height=data.shape(1);

  MultiArray<2,double> smooth(data.shape()),tmp(data.shape());
  vigra::Kernel1D<double> gauss;


  gauss.initGaussian(sigma_par);
  separableConvolveX(srcImageRange(data), destImage(tmp), kernel1d(gauss));
  separableConvolveY(srcImageRange(tmp), destImage(smooth), kernel1d(gauss));

  MultiArray<2,double> stxx(data.shape()),stxy(data.shape()),styy(data.shape());

  // calculate Structure Tensor at inner scale = sigma and outer scale = rho
  vigra::structureTensor(srcImageRange(smooth),destImage(stxx), destImage(stxy), destImage(styy),1.,1.);

  gauss.initGaussian(rho_par);
  separableConvolveX(srcImageRange(stxx), destImage(tmp), kernel1d(gauss));
  separableConvolveY(srcImageRange(tmp), destImage(stxx), kernel1d(gauss));
  separableConvolveX(srcImageRange(stxy), destImage(tmp), kernel1d(gauss));
  separableConvolveY(srcImageRange(tmp), destImage(stxy), kernel1d(gauss));
  separableConvolveX(srcImageRange(styy), destImage(tmp), kernel1d(gauss));
  separableConvolveY(srcImageRange(tmp), destImage(styy), kernel1d(gauss));

  MultiArray<2,double> matrix(Shape2(2,2)),ev(Shape2(2,2)),ew(Shape2(2,1));

   for (int y=0;y<data.shape(1);y++){
    for (int x=0;x<data.shape(0);x++){

      matrix(0,0)=stxx(x,y);
      matrix(1,1)=styy(x,y);
      matrix(0,1)=stxy(x,y);
      matrix(1,0)=stxy(x,y);
      vigra::symmetricEigensystemNoniterative(matrix,ew,ev);

      phi(x,y)=std::atan2(ev(1,0),ev(0,0));
      double coherence=ew(0,0)-ew(1,0);
      double c=std::min(K_par*coherence,1.);
      alpha(x,y)=alpha_par*c+(1-c)*beta_par;
      beta(x,y)=beta_par;
      }
  }
}

/********************************************************/
/*                                                      */
/*   anisotropicTotalVariationFilter                    */
/*                                                      */
/********************************************************/

/** \brief Performs Anisotropic Total Variation Regularization

The algorithm minimizes
\f[
\min_u \int_\Omega \frac{1}{2} (u-f)^2 + \sqrt{\nabla u^\top A \nabla u}\;dx\qquad\qquad(2)
\f]

where <em> \f$ f=f(x)\f$ </em> are the noisy data, <em> \f$ u=u(x)\f$ </em> are the smoothed data,<em>\f$ \nabla u \f$ </em>
is the image gradient in the sense of Total Variation and <em>\f$ A \f$ </em> is a locally varying symmetric, positive definite  2x2 matrix.

Matrix <em>\f$ A \f$ </em> is described by providing  for each data point a normalized eigenvector (via angle \f$ \phi \f$)
and two eigenvalues \f$ \alpha>0 \f$ and \f$ \beta>0 \f$.

\ref getAnisotropy() can be use to set up such data \f$ \phi,\alpha,\beta \f$ by providing a vector field normal to edges.

<b> Declarations:</b>

\code
namespace vigra {
  template <class stride1,class stride2,class stride3,class stride4,class stride5,class stride6>
  void anisotropicTotalVariationFilter(MultiArrayView<2,double,stride1> data,
                                       MultiArrayView<2,double,stride2> weight,
                                       MultiArrayView<2,double,stride3> phi,
                                       MultiArrayView<2,double,stride4> alpha,
                                       MultiArrayView<2,double,stride5> beta,
                                       MultiArrayView<2,double,stride6> out,
                                       int steps);
}
\endcode

\ref anisotropicTotalVariationFilter() implements a primal-dual algorithm to solve (2).

Input:
<table>
<tr><td><em>data</em>:</td><td>input data to be filtered. </td></tr>
<tr><td><em>steps</em>:</td><td>iteration steps.</td></tr>
<tr><td><em>weight</em> :</td><td>a point-wise weight (\f$ \ge 0 \f$ ) for the data term.</td></tr>
<tr><td><em>phi</em>,<em>alpha</em> and <em>beta</em> :</td><td> describe matrix \f$ A \f$, see above.</td></tr>
</table>

Output:
<table>
<tr><td><em>out</em> :</td><td>contains filtered data.</td></tr>
</table>

<b> Usage:</b>

E.g. with a solution-dependent adaptivity cf. [1], by updating the matrix \f$ A=A(u)\f$
in an outer loop:

<b>\#include</b> \<vigra/tv_filter.hxx\>

\code
MultiArray<2,double> data(Shape2(width,height)); //to be initialized
MultiArray<2,double> out (Shape2(width,height));
MultiArray<2,double> weight(Shape2(width,height));  //to be initialized
MultiArray<2,double> phi  (Shape2(width,height));
MultiArray<2,double> alpha(Shape2(width,height));
MultiArray<2,double> beta (Shape2(width,height));
double alpha0,beta0,sigma,rho,K;  //to be initialized
int outer_steps,inner_steps;//to be initialized

out=data; // data serves as initial value

for (int i=0;i<outer_steps;i++){
  getAnisotropy(out,phi,alpha,beta,alpha0,beta0,sigma,rho,K);  // sets phi, alpha, beta
  anisotropicTotalVariationFilter(data,weight,phi,alpha,beta,out,inner_steps);
  }
\endcode

[1] Frank Lenzen, Florian Becker, Jan Lellmann, Stefania Petra and Christoph Schn&ouml;rr, A Class of Quasi-Variational Inequalities for Adaptive Image Denoising and Decomposition, Computational Optimization and Applications, Springer, 2012.
*/
doxygen_overloaded_function(template <...>  void anisotropicTotalVariationFilter)

template <class stride1,class stride2,class stride3,class stride4,class stride5,class stride6>
void anisotropicTotalVariationFilter(MultiArrayView<2,double,stride1> data,MultiArrayView<2,double,stride2> weight,
                    MultiArrayView<2,double,stride3> phi,MultiArrayView<2,double,stride4> alpha,
                    MultiArrayView<2,double,stride5> beta,MultiArrayView<2,double,stride6> out,
                    int steps){

  using namespace multi_math;
  int width=data.shape(0),height=data.shape(1);

  MultiArray<2,double> temp1(data.shape()),temp2(data.shape()),vx(data.shape()),vy(data.shape()),u_bar(data.shape());
  MultiArray<2,double> rx(data.shape()),ry(data.shape());

  Kernel1D<double> Lx,LTx;
  Lx.initExplicitly(-1,0)=1,-1;                       // = Right sided finite differences for d/dx and d/dy
  Lx.setBorderTreatment(BORDER_TREATMENT_REFLECT);   //   with hom. Neumann boundary conditions
  LTx.initExplicitly(0,1)=-1,1;                     //  = Left sided finite differences for -d/dx  and -d/dy
  LTx.setBorderTreatment(BORDER_TREATMENT_ZEROPAD);  //   with hom. Dirichlet b. c.

  u_bar=out;

  double m=0;
  for (int y=0;y<data.shape(1);y++){
    for (int x=0;x<data.shape(0);x++){
      m=std::max(m,alpha(x,y));
      m=std::max(m,beta (x,y));
    }
  }
  m=std::max(m,1.);
  double tau=.9/m/std::sqrt(8.)*0.06;
  double sigma=.9/m/std::sqrt(8.)/0.06;


  for (int i=0;i<steps;i++){
    separableConvolveX(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroX(temp1);
    vx+=(sigma*temp1);
    separableConvolveY(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroY(temp1);
    vy+=(sigma*temp1);

    //project to constraint set
    for (int y=0;y<data.shape(1);y++){
      for (int x=0;x<data.shape(0);x++){
        double e1,e2,skp1,skp2;

        e1=std::cos(phi(x,y));
        e2=std::sin(phi(x,y));
        skp1=vx(x,y)*e1+vy(x,y)*e2;
        skp2=vx(x,y)*(-e2)+vy(x,y)*e1;
        vigra::detail::projectEllipse2D (skp1,skp2,alpha(x,y),beta(x,y),0.001,100);

        vx(x,y)=skp1*e1-skp2*e2;
        vy(x,y)=skp1*e2+skp2*e1;
      }
    }

    separableConvolveX(srcImageRange(vx),destImage(temp1),kernel1d(LTx));
    separableConvolveY(srcImageRange(vy),destImage(temp2),kernel1d(LTx));
    u_bar=out;
    out-=tau*(weight*(out-data)+(temp1+temp2));
    u_bar=2*out-u_bar;   //cf. Chambolle/Pock and Popov's algorithm
  }
}

/********************************************************/
/*                                                      */
/*   secondOrderTotalVariationFilter                    */
/*                                                      */
/********************************************************/

/** \brief Performs Anisotropic Total Variation Regularization

The algorithm minimizes

\f[
\min_u \int_\Omega \frac{1}{2} (u-f)^2 + \sqrt{\nabla u^\top A \nabla u}  + \gamma |Hu|_F\;dx \qquad\qquad (3)
\f]
where <em> \f$ f=f(x)\f$ </em> are the noisy data, <em> \f$ u=u(x)\f$ </em> are the smoothed data,<em>\f$ \nabla u \f$ </em>
is the image gradient in the sense of Total Variation, <em>\f$ A \f$ </em> is a locally varying
symmetric, positive-definite  2x2 matrix and <em>\f$ |Hu|_F \f$</em> is the Frobenius norm of the Hessian of \f$ u \f$.

Matrix <em>\f$ A \f$ </em> is described by providing  for each data point a normalized eigenvector (via angle \f$ \phi \f$)
and two eigenvalues \f$ \alpha>0 \f$ and \f$ \beta>0 \f$.
\ref getAnisotropy() can be use to set up such data \f$ \phi,\alpha, \beta \f$ by providing a vector field normal to edges.

\f$ \gamma>0 \f$ is the locally varying regularization parameter for second order.

<b> Declarations:</b>

\code
namespace vigra {
  template <class stride1,class stride2,...,class stride9>
  void secondOrderTotalVariationFilter(MultiArrayView<2,double,stride1> data,
                                       MultiArrayView<2,double,stride2> weight,
                                       MultiArrayView<2,double,stride3> phi,
                                       MultiArrayView<2,double,stride4> alpha,
                                       MultiArrayView<2,double,stride5> beta,
                                       MultiArrayView<2,double,stride6> gamma,
                                       MultiArrayView<2,double,stride7> xedges,
                                       MultiArrayView<2,double,stride8> yedges,
                                       MultiArrayView<2,double,stride9> out,
                                       int steps);
}
\endcode

\ref secondOrderTotalVariationFilter() implements a primal-dual algorithm to solve (3).

Input:
<table>
<tr><td><em>data</em>: </td><td> the input data to be filtered. </td></tr>
<tr><td><em>steps</em> : </td><td> number of iteration steps.</td></tr>
<tr><td><em>out</em> : </td><td>contains the filtered data.</td></tr>
<tr><td><em>weight</em> : </td><td> point-wise weight (\f$ \ge 0\f$ ) for the data term.</td></tr>
<tr><td><em>phi</em>,<em>alpha</em>,<em>beta</em>: </td><td> describe matrix \f$ A\f$, see above.</td></tr>
<tr><td><em> xedges </em> and <em> yedges </em>: </td><td> binary arrays indicating the
presence of horizontal (between (x,y) and (x+1,y)) and vertical edges (between (x,y) and (x,y+1)).
These data are considered in the calculation of \f$ Hu\f$, such that
finite differences across edges are artificially set to zero to avoid second order smoothing over edges.</td></tr>
</table>

<b> Usage:</b>

E.g. with a solution-dependent adaptivity (cf.[1]), by updating the matrix \f$ A=A(u)\f$
in an outer loop:

<b>\#include</b> \<vigra/tv_filter.hxx\>

\code
MultiArray<2,double> data(Shape2(width,height)); //to be initialized
MultiArray<2,double> out(Shape2(width,height));
MultiArray<2,double> weight(Shape2(width,height))); //to be initialized
MultiArray<2,double> phi(Shape2(width,height);
MultiArray<2,double> alpha(Shape2(width,height);
MultiArray<2,double> beta(Shape2(width,height));
MultiArray<2,double> gamma(Shape2(width,height));  //to be initialized
MultiArray<2,double> xedges(Shape2(width,height));  //to be initialized
MultiArray<2,double> yedges(Shape2(width,height));  //to be initialized


double alpha0,beta0,sigma,rho,K;  //to be initialized
int outer_steps,inner_steps;//to be initialized

out=data; // data serves as initial value

for (int i=0;i<outer_steps;i++){

  getAnisotropy(out,phi,alpha,beta,alpha0,beta0,sigma,rho,K);  // sets phi, alpha, beta
  secondOrderTotalVariationFilter(data,weight,phi,alpha,beta,gamma,xedges,yedges,out,inner_steps);
}
\endcode


[1] Frank Lenzen, Florian Becker, Jan Lellmann, Stefania Petra and Christoph Schn&ouml;rr, A Class of Quasi-Variational Inequalities for Adaptive Image Denoising and Decomposition, Computational Optimization and Applications, Springer, 2012.
*/
doxygen_overloaded_function(template <...> void secondOrderTotalVariationFilter)

template <class stride1,class stride2,class stride3,class stride4,class stride5,class stride6,class stride7,class stride8,class stride9>
void secondOrderTotalVariationFilter(MultiArrayView<2,double,stride1> data,
                            MultiArrayView<2,double,stride2> weight,MultiArrayView<2,double,stride3> phi,
                            MultiArrayView<2,double,stride4> alpha,MultiArrayView<2,double,stride5> beta,
                            MultiArrayView<2,double,stride6> gamma,
                            MultiArrayView<2,double,stride7> xedges,MultiArrayView<2,double,stride8> yedges,
                    MultiArrayView<2,double,stride9> out,
                            int steps){

  using namespace multi_math;
  int width=data.shape(0),height=data.shape(1);

  MultiArray<2,double> temp1(data.shape()),temp2(data.shape()),vx(data.shape()),vy(data.shape()),u_bar(data.shape());
  MultiArray<2,double> rx(data.shape()),ry(data.shape());
  MultiArray<2,double> wx(data.shape()),wy(data.shape()),wz(data.shape());


  Kernel1D<double> Lx,LTx;
  Lx.initExplicitly(-1,0)=1,-1;                       // = Right sided finite differences for d/dx and d/dy
  Lx.setBorderTreatment(BORDER_TREATMENT_REFLECT);   //   with hom. Neumann boundary conditions
  LTx.initExplicitly(0,1)=-1,1;                     //  = Left sided finite differences for -d/dx  and -d/dy
  LTx.setBorderTreatment(BORDER_TREATMENT_ZEROPAD);  //   with hom. Dirichlet b. c.

  u_bar=out;

  double m=0;
  for (int y=0;y<data.shape(1);y++){
    for (int x=0;x<data.shape(0);x++){
      m=std::max(m,alpha(x,y));
      m=std::max(m,beta (x,y));
      m=std::max(m,gamma(x,y));
     }
  }
  m=std::max(m,1.);
  double tau=.1/m;//std::sqrt(8)*0.06;
  double sigma=.1;//m;/std::sqrt(8)/0.06;

  //std::cout<<"tau= "<<tau<<std::endl;

  for (int i=0;i<steps;i++){

    separableConvolveX(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroX(temp1);
    vx+=(sigma*temp1);
    separableConvolveY(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroY(temp1);
    vy+=(sigma*temp1);


    // update wx
    separableConvolveX(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroX(temp1);
    temp1*=xedges;
    separableConvolveX(srcImageRange(temp1),destImage(temp2),kernel1d(LTx));
    wx-=sigma*temp2;//(-Lx'*(xedges.*(Lx*u)));

    //update wy
    separableConvolveY(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroY(temp1);
    temp1*=yedges;
    separableConvolveY(srcImageRange(temp1),destImage(temp2),kernel1d(LTx));
    wy-=sigma*temp2;//(-Ly'*(yedges.*(Ly*u)));


    //update wz
    #if (VIGRA_MIXED_2ND_DERIVATIVES)
    separableConvolveY(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroY(temp1);
    temp1*=yedges;
    separableConvolveX(srcImageRange(temp1),destImage(temp2),kernel1d(LTx));
    wz-=sigma*temp2;//-Lx'*(yedges.*(Ly*u))

    separableConvolveX(srcImageRange(u_bar),destImage(temp1),kernel1d(Lx));setZeroX(temp1);
    temp1*=xedges;
    separableConvolveY(srcImageRange(temp1),destImage(temp2),kernel1d(LTx));
    wz-=sigma*temp2;//-Ly'*(xedges.*(Lx*u)));

    #endif


    //project to constraint sets
    for (int y=0;y<data.shape(1);y++){
      for (int x=0;x<data.shape(0);x++){
        double e1,e2,skp1,skp2;

        //project v
        e1=std::cos(phi(x,y));
        e2=std::sin(phi(x,y));
        skp1=vx(x,y)*e1+vy(x,y)*e2;
        skp2=vx(x,y)*(-e2)+vy(x,y)*e1;
        vigra::detail::projectEllipse2D (skp1,skp2,alpha(x,y),beta(x,y),0.001,100);
        vx(x,y)=skp1*e1-skp2*e2;
        vy(x,y)=skp1*e2+skp2*e1;

        //project w
        double l=sqrt(wx(x,y)*wx(x,y)+wy(x,y)*wy(x,y)+wz(x,y)*wz(x,y));
        if (l>gamma(x,y)){
          wx(x,y)=gamma(x,y)*wx(x,y)/l;
          wy(x,y)=gamma(x,y)*wy(x,y)/l;
          #if (VIGRA_MIXED_2ND_DERIVATIVES)
          wz(x,y)=gamma(x,y)*wz(x,y)/l;
          #endif
        }
      }
    }

    separableConvolveX(srcImageRange(vx),destImage(temp1),kernel1d(LTx));
    separableConvolveY(srcImageRange(vy),destImage(temp2),kernel1d(LTx));

    u_bar=out;
    out-=tau*(weight*(out-data)+temp1+temp2);


    // update wx
    separableConvolveX(srcImageRange(wx),destImage(temp1),kernel1d(Lx));setZeroX(temp1);
    temp1*=xedges;
    separableConvolveX(srcImageRange(temp1),destImage(temp2),kernel1d(LTx));
    out+=tau*temp2; // (-1)^2


    //update wy
    separableConvolveY(srcImageRange(wy),destImage(temp1),kernel1d(Lx));setZeroY(temp1);
    temp1*=yedges;
    separableConvolveY(srcImageRange(temp1),destImage(temp2),kernel1d(LTx));
    out+=tau*temp2;

    //update wz
    #if (VIGRA_MIXED_2ND_DERIVATIVES)

    separableConvolveY(srcImageRange(wz),destImage(temp1),kernel1d(Lx));setZeroY(temp1);
    temp1*=yedges;
    separableConvolveX(srcImageRange(temp1),destImage(temp2),kernel1d(LTx));
    out+=tau*temp2;

    separableConvolveX(srcImageRange(wz),destImage(temp1),kernel1d(Lx));setZeroX(temp1);
    temp1*=xedges;
    separableConvolveY(srcImageRange(temp1),destImage(temp2),kernel1d(LTx));
    out+=tau*temp2;

    #endif

    u_bar=2*out-u_bar;   //cf. Chambolle/Pock and Popov's algorithm

  }
}

//@}
} // closing namespace vigra

#endif // VIGRA_TV_FILTER_HXX