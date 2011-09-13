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

#ifndef VIGRA_INVARIANT_FEATURES3D_HXX
#define VIGRA_INVARIANT_FEATURES3D_HXX

#include <complex>
#include "config.hxx"
#include "error.hxx"
#include "utilities.hxx"
#include "mathutil.hxx"
#include "array_vector.hxx"
#include "matrix.hxx"
#include "tinyvector.hxx"
#include "quaternion.hxx"
#include "harmonics.hxx"
#include "clebsch-gordan.hxx"
#include "multi_fft.hxx"
#include "wigner-matrix.hxx"
#include "functortraits.hxx"
#include <vigra/multi_math.hxx>

namespace vigra {

/** \addtogroup InvariantFeatures local invariant features
*/
//@{

/********************************************************/
/*                                                      */
/*           Spherical Harmonic cross correlation       */
/*                                                      */
/********************************************************/

/** \brief Computes the complete cross correlation over all possible 3D
    rotations of two Spherical Harmonic (SH) representations. See harmonics.hxx
    on how to compute SH representations of 3D volume data.

    Returns (normalized) correlation value as well as the roation offset at maximum correlation
    of the two signals.

    Uses the fast Fourier space corrleation algorithm presented in:

    <i>Fehr, J., Reisert, M. & Burkhardt, H. "Fast and Accurate Rotation Estimation on the 2-Sphere without Correspondences", Computer Vision -- ECCV 2008, LNCS 5303.</i>

    More details and theoretic analysis in:

    <i>2009: Fehr, J. "Local invariant Features for 3D Image Analysis", Dissertation, 2009, pages:47-53</i> 

    NOTE: this is actually not an invariant feature! Used for one-on-one euqivariant comparison. 
    Use SHautocorr if a feature is needed. 

    The Euler angle representation is given in ZYZ' notation. 

    \param SH_A	first SH expansion
    \param SH_B	second SH expansion
    \param padsize size of the padding for FFT sinc interpolation. Larger padding increases the
    accuracy of the rotation estimation at increased cost. Theoretic limit of the angular resolution is: Err >= 180°/(2*b_max+padzize), where b_max is the maximum band of the SH expansion. Reasonable
    values for padsize are [0,32,64,...,512] 
    \param normalize bool:compute normalized cross corelation with perfect match at max==1 or unnormalized
    \param max returns maximum correlation value
    \param phi returns rotation offset in first Euler angle [0...2PI]
    \param theta returns rotation offset in second Euler angle [0...PI]
    \param psi returns rotation offset in third Euler angle [0...2PI]
    \param corMat returns the 3D (theta,phi,psi) correlation space
    \WARNING: there are natural ambiguities of the rotation representation in Euler angles, especially for theta=0  

*/
template <typename REAL>
void 
SHxcorr(std::vector<std::vector<FFTWComplex<REAL> >  > const &SH_A, 
        std::vector<std::vector<FFTWComplex<REAL> >  > const &SH_B,
        unsigned int padsize, 
        bool normalize, 
        REAL &max, REAL &phi, REAL &theta, REAL &psi, 
        MultiArray<3,REAL > &corMat)
{
    int band=SH_A.size()-1;
    TinyVector<int,3> baseshape(2*band+padsize,2*band+padsize,2*band+padsize);
    TinyVector<int,3> padshape(fftwBestPaddedShape(baseshape));
    MultiArray<3,FFTWComplex<REAL> > C_fft;
    C_fft.reshape(padshape,0);

    int size = padshape[0];

    WignerMatrix<REAL> W(band);
    W.compute_D(band);
    FFTWComplex<REAL> Im(0, 1);

    for (int l=0;l<=band;l++)
    {
        for (int m1=-l;m1<=l;m1++)
        {
            for (int m2=-l;m2<=l;m2++)
            {
                for (int m3=-l;m3<=l;m3++)
                {
                    C_fft((m2+size+1)%(size),(m1+size+1)%(size),(m3+size+1)%(size)) +=  
                    (W.get_d(l,m1,m2) * W.get_d(l,m2,m3) * (SH_A[l][l+m1]) * conj(SH_B[l][l+m3])) *
                       pow(Im, (m1 + 2*m2 + m3));
                }
            }
        }
    }

    MultiArray<3,FFTWComplex<REAL> > res(C_fft.shape());

    //FFT^-1
    fourierTransformInverse(C_fft,res);
    
    REAL variance = 1;
    if(normalize)
    {
        REAL varianceA=0, varianceB=0;
        for (int l=0;l<=band;l++)
        {
            REAL tmpA=0, tmpB=0;
            for (int m=-l;m<=l;m++)
            {
                tmpA += (REAL)(std::pow((double)real(SH_A[l][m+l]),2.0) + 
                               std::pow((double)imag(SH_A[l][m+l]),2.0));
                tmpB += (REAL)(std::pow((double)real(SH_B[l][m+l]),2.0) + 
                               std::pow((double)imag(SH_B[l][m+l]),2.0));
            }
            varianceA += tmpA;
            varianceB += tmpB;
        }

        variance = std::sqrt(varianceA*varianceB)/(size*size*size);
        if (variance == 0) 
        {
            variance = FLT_MIN;
        }
    }

    corMat.reshape(res.shape(),0);
    
    TinyVector<int, 3> hotpos1;
    max = 0;
    for (int z=0;z<corMat.shape()[0];z++)
    {
        for (int y=0;y<corMat.shape()[1];y++)
        {
            for (int x=0;x<corMat.shape()[2];x++)
            {
                REAL tmp = abs(res(z,y,x)/variance);
                if (tmp>max)
                {
                    max=tmp;
                    hotpos1[0]=z;
                    hotpos1[1]=y;
                    hotpos1[2]=x;
                }
                corMat(z,y,x) = tmp;
            }
        }
    }
    REAL grid = REAL(2.0*M_PI/(size));

    //phi
    if(hotpos1[1]*grid >M_PI)
        phi = REAL(M_PI + (2*M_PI -hotpos1[1]*grid ));
    else
        phi = REAL(M_PI - hotpos1[1]*grid);

    //theta     
    if(hotpos1[0]*grid >M_PI)
        theta = REAL(2*M_PI - (hotpos1[0])*grid);
    else
        theta = REAL(hotpos1[0]*grid);

    //psi
    if(hotpos1[2]*grid >M_PI)
        psi = REAL(M_PI + (2*M_PI -hotpos1[2]*grid ));
    else
        psi = REAL(M_PI - hotpos1[2]*grid);
}

//@{

/********************************************************/
/*                                                      */
/*           Spherical Harmonic auto correlation       */
/*                                                      */
/********************************************************/

/** \brief Computes integrates over  complete auto correlation over all possible 3D
    rotations of a Spherical Harmonic (SH) representation. See harmonics.hxx
    on how to compute SH representations of 3D volume data.

    Returns scalar invariant feature.

    More details and theoretic analysis in:

    <i>2009: Fehr, J. "Local invariant Features for 3D Image Analysis", Dissertation, 2009, pages:73-74</i> 

    \param SH_A  SH expansion
    \param padsize size of the padding for FFT sinc interpolation (see SHxcorr for details). 
    \param f functor with scalar non-linear mapping of the correlation values before integration. 
    Reasonable functors could be f(x)=x^2, sqrt(x), 1/x, ...
    \param res return feature value
*/
template <typename REAL, typename FUNCTOR>
void 
SHautocorr(std::vector<std::vector<FFTWComplex<REAL> >  > const &SH_A, 
           unsigned int padsize, FUNCTOR &f, REAL &res)
{
    REAL cmax,phi,theta,psi=0;
    MultiArray<3,float > corMat;
    SHxcorr(SH_A,SH_A,padsize,true,cmax,phi,theta,psi,corMat);
    res=0;
    for ( MultiArray<3,float >::iterator p=corMat.begin();p!=corMat.end();p++)
        res *= f(*p);
}

/********************************************************/
/*                                                      */
/*           Spherical Harmonic power-specturm          */
/*                                                      */
/********************************************************/

/** \brief Computes the power spectrum of the SH representation. 
    See harmonics.hxx on how to compute SH representations of 3D volume data.
    
    Returns vectorial invariant feature.

    Well known Idea, one of the first papers using the SH power-spectrum as a feature was:

    <i>M. Kazhdan. Rotation invariant spherical harmonic representation of 3d shape descriptors. Symp. on Geom. Process., 2003.</i>
 
    More details and theoretic analysis in:

    <i>2009: Fehr, J. "Local invariant Features for 3D Image Analysis", Dissertation, 2009, pages:68-69</i> 

    \param SH_A  SH expansion
    \param res return vector of feature values
*/
template <typename REAL>
void 
SHabs(std::vector<std::vector<FFTWComplex<REAL> >  > const &SH_A, 
      std::vector<REAL> &res)
{
    res.resize(SH_A.size());
    for (int l=0;l<(int)SH_A.size();l++)
    {
        res[l]=0;
        for (int m=-l;m<=l;m++)
            res[l] += (REAL)(std::pow((double)real(SH_A[l][m+l]),2.0)+ 
                             std::pow((double)imag(SH_A[l][m+l]),2.0));
    }   
}

//this version computes feature at all positions, use Array2SH to get SH_A 
template <typename REAL>
void 
SHabs(std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > >const &SH_A, 
      std::vector<MultiArray<3,REAL> >&res)
{
    using namespace multi_math;
    res.resize(SH_A.size());
    for (int l=0;l<(int)SH_A.size();l++)
    {
        res[l].reshape(SH_A[l][0].shape(),0);
        for (int m=-l;m<=l;m++)
        {
            typename MultiArray<3,REAL >::iterator q=res[l].begin();
            typename MultiArray<3,FFTWComplex<REAL> >::iterator p=SH_A[l][m+l].begin();
            for (;p<SH_A[l][m+l].end();p++,q++)
            {
                *q += REAL(std::pow((double)real(*p),2.0)+ std::pow((double)imag(*p),2.0));
            }
        }
    }
}


/********************************************************/
/*                                                      */
/*           Spherical Harmonic bi-specturm             */
/*                                                      */
/********************************************************/

/** \brief Computes the bi-spectrum of the SH representation. 
    See harmonics.hxx on how to compute SH representations of 3D volume data.
    
    Returns vectorial invariant feature.

    This is a complete feature! Turns out to be less robust for real world data, hence SH expansion
    should be limited to low bands only. Original paper:

    <i>S. Wenndt and S. Shamsunder. Bispectrum features for robust speaker identification. Acoustics, Speech, and Signal Processing, IEEE International Conference on, 2:1095, 1997.</i>
 
    More details and theoretic analysis in:

    <i>2009: Fehr, J. "Local invariant Features for 3D Image Analysis", Dissertation, 2009, pages:75-76</i> 

    \param SH_A  SH expansion
    \param res return vector of feature values
*/

template <typename REAL>
void 
SHbispec(std::vector<std::vector<FFTWComplex<REAL> >  > const &SH_A, 
         std::vector<FFTWComplex<REAL> > &res)
{
    res.clear();
    for (int l=0;l<(int)SH_A.size();l++)
    {
        for (int l1=0;l1<=l;l1++)
        {
            for (int l2=0;l2<=l;l2++)
            {
                // some l,l1,l2 combinations are not allowed (triangular inequallity) - rejected by Clebschgordan
                try
                {
                    FFTWComplex<REAL> f;
                    f = (REAL)0.0;

                    for(int m=-l;m<=l;m++)
                    {
                        FFTWComplex<REAL> tmp_feature;
                        tmp_feature = (REAL)0.0;
                        for(int m1=std::max(-l1,m-l2);m1<=std::min(l1,m+l2);m1++)
                        {
                            FFTWComplex<REAL> s;
                            s = (REAL)0.0;
                            s.real() = (REAL)clebschGordan(l1, m1, l2, m-m1, l, m);
                            tmp_feature += s * conj(SH_A[l1][m1+l1]) * conj(SH_A[l2][m-m1+l2]);
                        }
                        f += SH_A[l][m+l] * tmp_feature;

                    }

                    res.push_back(f);
                }
                catch (ContractViolation &){}
            }
        }
    }
}

//this version computes feature at all positions, use Array2SH to get SH_A 
template <typename REAL>
void 
SHbispec(std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > >const &SH_A, 
         std::vector<MultiArray<3,FFTWComplex<REAL> > > &res)
{
    res.clear();
    for (int l=0;l<(int)SH_A.size();l++)
    {
        for (int l1=0;l1<=l;l1++)
        {
            for (int l2=0;l2<=l;l2++)
            {
                // some l,l1,l2 combinations are not allowed (triangular inequallity) - rejected by Clebschgordan
                try
                {
                    MultiArray<3,FFTWComplex<REAL> > f;
                    f.reshape(SH_A[0][0].shape(),0);

                    for(int m=-l;m<=l;m++)
                    {
                        MultiArray<3,FFTWComplex<REAL> > tmp_feature;
                        tmp_feature.reshape(SH_A[0][0].shape(),0);;
                        for(int m1=std::max(-l1,m-l2);m1<=std::min(l1,m+l2);m1++)
                        {
                            FFTWComplex<REAL> s;
                            s = (REAL)0.0;
                            s.real() = (REAL)clebschGordan(l1, m1, l2, m-m1, l, m);
                            typename MultiArray<3,FFTWComplex<REAL> >::iterator h = 
                                SH_A[l2][m-m1+l2].begin();
                            typename MultiArray<3,FFTWComplex<REAL> >::iterator q =
                                SH_A[l1][m1+l1].begin();
                            typename MultiArray<3,FFTWComplex<REAL> >::iterator p = 
                                tmp_feature.begin();
                            for (;p<tmp_feature.end();p++,q++,h++)
                                *p += s * conj(*q) * conj(*h);
                        }
                        typename MultiArray<3,FFTWComplex<REAL> >::iterator h=SH_A[l][m+l].begin();
                        typename MultiArray<3,FFTWComplex<REAL> >::iterator q=f.begin();
                        typename MultiArray<3,FFTWComplex<REAL> >::iterator p = tmp_feature.begin();
                        for (;p<tmp_feature.end();p++,q++,h++)
                            *q += *h * *p;

                    }

                    res.push_back(f);

                }
                catch (ContractViolation &){}
            }
        }
    }
}

} // namespace vigra 

#endif // VIGRA_INVARIANT_FEATURES3D_HXX
