/************************************************************************/
/*                                                                      */
/*        Copyright 2009-2010 by Ullrich Koethe and Janis Fehr          */
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

#ifndef VIGRA_WIGNER_MATRIX_HXX
#define VIGRA_WIGNER_MATRIX_HXX

#include <complex>
#include "config.hxx"
#include "error.hxx"
#include "utilities.hxx"
#include "mathutil.hxx"
#include "array_vector.hxx"
#include "matrix.hxx"
#include "tinyvector.hxx"
#include "quaternion.hxx"
#include "clebsch-gordan.hxx"
#include "multi_fft.hxx"


namespace vigra {


    /*!
     *  \class WignerMatrix 
     *  \brief computation of Wigner D matrix + rotation functions 
     *         in SH,VH and R³
     *
     * All rotations in Euler zyz' convention	
     *
     * WARNING: not thread safe! use a new instance of WignerMatrix
     * for each thread!!!
     */
template <class Real>
class WignerMatrix
{   
  public:
  
    typedef FFTWComplex<Real> Complex;
    typedef ArrayVector<ArrayVector<ArrayVector<Complex> > > NestedArray;
    
         /*! \brief constructor
          * 
          * \param l_max	maximum expansion band (used to pre-compute 
          *			the D matrix)
          *
          */
    WignerMatrix(int l_max);
    
         /*! \brief Compute D with fixed theta = pi/2, phi=0, psi=0.
          *
          * \param band	expansion band
          *
          */
    void compute_D(int band);

         /*! \brief Compute D for arbitrary rotations.
          *
          * \param l     	expansion band
          * \param phi	rotation angle	
          * \param theta	rotation angle
          * \param psi 	rotation angle
          *
          */
    void compute_D(int l, Real phi, Real theta, Real psi);

         /*! \brief  Get the (n,m) entry of D.
          *
          * \param l        expansion band
          * \param n      	
          * \param m    
          */
    Complex get_d(int l, int h, int m) const
    {
        if (l>0)
        {
            std::string message = std::string("WignerMatrix::get_D(): index out of bounds: l=");
            message << l << " l_max=" << D.size() << " m=" << m << " n=" << h << "\n";
                                  
            vigra_precondition(l < (int)D.size() && m+l <= 2*l+1 &&
                               h+l <= 2*l+1 && m+l >= 0 && h+l >= 0,
                               message.c_str());
            return D[l](h+l, m+l);
        }
        else 
        {
            return Complex(Real(1.0));
        }
    }

         /*! \brief Return the rotation matrix D for the lth band.
          *
          * \param l        expansion band
          */
    Matrix<Complex> const & get_D(int l) const
    {
        std::string message = std::string("WignerMatrix::get_D(): index out of bounds: l=");
        message << l << " l_max=" << l_max << "\n";
                              
        vigra_precondition(l > 0 && l <= l_max, message.c_str());
        return D[l];
    }

         /*! \brief Rotate in PH.
          *
          * \param PH       input PH expansion 
          * \param phi      rotation angle
          * \param theta    rotation angle
          * \param psi      rotation angle
          *
          * \retval PHresult PH expansion   
          */
    void rotatePH(std::vector<std::vector<std::vector<vigra::FFTWComplex<Real> > > > const & PH, Real phi, Real theta, Real psi, std::vector<std::vector<std::vector<vigra::FFTWComplex<Real> > > >  & PHresult);

    void rotateSH(std::vector<std::vector<vigra::FFTWComplex<Real> > >  const & SH, Real phi, Real theta, Real psi, std::vector<std::vector<vigra::FFTWComplex<Real> > > & SHresult);



  private:
    int l_max;
    int cg1cnt;
    ArrayVector<double> CGcoeff;
    ArrayVector<Matrix<Complex> > D;
};

template <class Real>
WignerMatrix<Real>::WignerMatrix(int band)
: l_max(0),
  cg1cnt(0)
{
    //precompute clebschGordan coeffs
    for (int l = 2; l <= band+1; l++)
    {
        for(int m = -l; m <= l ; m++)
        {
            for(int n = -l; n <= l ; n++)
            {
                for (int m2 = -1; m2 <= 1; m2++)
                {
                    for (int n2 = -1; n2 <= 1; n2++)
                    {
                        int m1 = m-m2;
                        int n1 = n-n2;
                        if (m1 > -l && m1 < l && n1 > -l && n1 < l)
                        {
                            CGcoeff.push_back((clebschGordan(l-1,m1,1,m2,l,m))*(clebschGordan(l-1,n1,1,n2,l,n)));
                        }
                    }
                }
            }
        }
    }
}

template <class Real>
void
WignerMatrix<Real>::compute_D(int l, Real phi, Real theta, Real psi)
{
    Real s = std::sin(theta);
    Real c = std::cos(theta);
    
    FFTWComplex<Real> i(0.0, 1.0);
    FFTWComplex<Real> eiphi = exp(i*phi);
    FFTWComplex<Real> emiphi = exp(-i*phi);
    FFTWComplex<Real> eipsi = exp(i*psi);
    FFTWComplex<Real> emipsi = exp(-i*psi);
    
    if (D.size() < (std::size_t)(l+1)) 
        D.resize(l+1);
    D[1].reshape(MultiArrayShape<2>::type(3,3));
    
    D[1](0,0) = emipsi * Complex(Real(0.5*(1.0+c))) * emiphi;
    D[1](1,0) = Complex(Real(-s/M_SQRT2)) * emiphi;
    D[1](2,0) = eipsi * Complex(Real(0.5*(1.0-c))) * emiphi;
    D[1](0,1) = emipsi * Complex(Real(s/M_SQRT2));
    D[1](1,1) = Complex(Real(c));
    D[1](2,1) = eipsi * Complex(Real(-s/M_SQRT2));
    D[1](0,2) = emipsi * Complex(Real(0.5*(1.0-c))) * eiphi; 
    D[1](1,2) = Complex(Real(s/M_SQRT2)) * eiphi;
    D[1](2,2) = eipsi * Complex(Real(0.5*(1.0+c))) * eiphi;

    l_max = 1;
    cg1cnt = 0;
    if(l > 1)
        compute_D( l);
}


template <class Real>
void
WignerMatrix<Real>::compute_D(int l)
{
    if (D.size() < (std::size_t)(l+1) ) 
    {
        D.resize(l+1);
        l_max = 0;
    }

    if (l==1)
    {
        //precompute D0 =1 and D1 = (90 degree rot)
        D[1].reshape(MultiArrayShape<2>::type(3,3));
        D[1](0,0) = Real(0.5);
        D[1](0,1) = Real(0.5*M_SQRT2);
        D[1](0,2) = Real(0.5);
        D[1](1,0) = Real(-0.5*M_SQRT2);
        D[1](1,1) = Real(0.0);
        D[1](1,2) = Real(0.5*M_SQRT2);
        D[1](2,0) = Real(0.5);
        D[1](2,1) = Real(-0.5*M_SQRT2);
        D[1](2,2) = Real(0.5);
        l_max = 1;
        cg1cnt = 0;
    }
    else
    {
        //compute D2-Dl_max recursive
        if (l>l_max+1)
        {
            compute_D(l-1);
        }

        D[l].reshape(MultiArrayShape<2>::type(2*l+1,2*l+1));
        D[l].init(Real(0.0));	   
        
        for(int m = -l; m <= l ; m++)
        {
            for(int n = -l; n <= l ; n++)
            {
                for (int m2 = -1; m2 <= 1; m2++)
                {
                    for (int n2 = -1; n2 <= 1; n2++)
                    {
                        int m1 = m-m2;
                        int n1 = n-n2;
                        if ((m1 > -l) && (m1 < l) && (n1 > -l) && (n1 < l))
                        {
                            D[l](m+l,n+l) += D[1](m2+1,n2+1) * D[l-1](m1+l-1,n1+l-1) * Real(CGcoeff[cg1cnt++]);
                        }
                    }
                }
            }
        }

        l_max = l;
    }
}

template <class Real>
void
WignerMatrix<Real>::rotateSH(std::vector<std::vector< FFTWComplex<Real> > > const  & SH, 
                             Real phi, Real theta, Real psi, 
                             std::vector<std::vector< FFTWComplex<Real> > > & SHresult)
{
    int band = SH.size()-1;
    compute_D(band, phi, theta, psi);

    SHresult.resize(SH.size());

    for (int l=0; l<(int)SH.size(); l++)
    {
        SHresult[l].resize(SH[l].size());
        for(int m=-l; m<=l; m++)
        {
            Complex tmp = 0;
            for (int h=-l; h<=l; h++)
            {
                tmp += get_d(l,h,m) * SH[l][h+l];
            }

            SHresult[l][m+l] = tmp;
        }
    }
}


template <class Real>
void 
WignerMatrix<Real>::rotatePH(std::vector<std::vector<std::vector< FFTWComplex<Real> > > > const  & PH, 
                             Real phi, Real theta, Real psi,
                             std::vector<std::vector<std::vector< FFTWComplex<Real> > > > & PHresult)
{
    int band = PH.size()-1;
    compute_D(band, phi, theta, psi);

    PHresult.resize(PH.size());

    for(int n=1; n<(int)PH.size(); n++)
    {
        PHresult[n].resize(PH[n].size());
        for (int l=0; l<(int)PH[n].size(); l++)
        {
            PHresult[n][l].resize(PH[n][l].size());
            for(int m=-l; m<=l; m++)
            {
                Complex tmp = 0;
                for (int h=-l; h<=l; h++)
                {
                    tmp += get_d(l,h,m) * PH[n][l][h+l];
                }

                PHresult[n][l][m+l] = tmp;
            }
        }
    }
}

} // namespace vigra 

#endif // VIGRA_WIGNER_MATRIX_HXX
