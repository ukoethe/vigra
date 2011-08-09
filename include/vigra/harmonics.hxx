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

#ifndef VIGRA_HARMONICS_HXX
#define VIGRA_HARMONICS_HXX

#include <complex>
#include "config.hxx"
#include "error.hxx"
#include "utilities.hxx"
#include "mathutil.hxx"
#include "array_vector.hxx"
#include "matrix.hxx"
#include "tinyvector.hxx"
#include "quaternion.hxx"
#include "wigner-matrix.hxx"
#include "clebsch-gordan.hxx"
#include "multi_fft.hxx"
#include "multi_math.hxx"
#include "bessel.hxx"


namespace vigra {

namespace detail  {

// computes the normalization for SH base functions
inline double realSH(double l, double m)
{
    return std::sqrt((2.0*l + 1.0) / (4.0*M_PI*facLM(l,m)));

}

template <typename REAL>
inline REAL fac(REAL in)
{   
    REAL temp = 1;
    for (int i=2;i<=in;i++)
    {   
        temp *= i;
    }     
    return temp;
}           


template <typename REAL, typename T> 
TinyVector<REAL, 3> centerOfBB(const MultiArray<3,T>  &A) 
{
    return TinyVector<REAL, 3>(A.shape()) /= 2.0;                        
}


template <typename REAL>
void 
eulerAngles(REAL sphereRadius_um,
            MultiArray<3,REAL>& phi,
            MultiArray<3,REAL>& theta,
            MultiArray<3,REAL>& psi,
            REAL gaussWidthAtHalfMaximum_um,
            TinyVector<REAL,3> voxelSize=TinyVector<REAL,3>(1.0))
{
    vigra_fail("eulerAngles is untested");
    REAL radiusLev = sphereRadius_um /voxelSize[0] + gaussWidthAtHalfMaximum_um*3;
    REAL radiusRow = sphereRadius_um /voxelSize[1] + gaussWidthAtHalfMaximum_um*3;
    REAL radiusCol = sphereRadius_um /voxelSize[2] + gaussWidthAtHalfMaximum_um*3;

    int intRadiusLev = (int)std::ceil( radiusLev);
    int intRadiusRow = (int)std::ceil( radiusRow);
    int intRadiusCol = (int)std::ceil( radiusCol);

    MultiArrayShape<3>::type NewShape( intRadiusLev*2 + 1, intRadiusRow*2 + 1, intRadiusCol*2 + 1 );

    TinyVector<int,3> M;
    M[0]= NewShape[0] / 2;
    M[1]= NewShape[1] / 2;
    M[2]= NewShape[2] / 2;

    phi.reshape( NewShape,0 );
    theta.reshape( NewShape,0 );
    psi.reshape( NewShape,0 );


    for (int z=0;z<NewShape[0];++z)
    {
        for (int y=0;y<NewShape[1];++y)
        {
            for (int x=0;x<NewShape[2];++x)
            {
                REAL X = (REAL)x *voxelSize[2];
                REAL Y = (REAL)y *voxelSize[1];
                REAL Z = (REAL)z *voxelSize[0];
                int   iZ = z;

                // calculate psi
                phi(z,y,x ) = std::atan2( Y, X );

                // calculate theta
                REAL r = std::sqrt( X*X + Y*Y + Z*Z );
                REAL alpha = (REAL)std::asin( Z / r );

                if( iZ == 0 )
                {
                    theta( z,y,x ) = M_PI * 0.5;
                }
                else
                {
                    theta( z,y,x ) = M_PI * 0.5 - alpha;
                }
            }
        }
    }
}

} // namespace detail

template <typename REAL>
MultiArray<3,REAL>
binarySphereREAL(REAL radius_um, 
                 REAL gaussWidthAtHalfMaximum_um, 
                 TinyVector<REAL,3> voxelSize=TinyVector<REAL,3>(1.0))
{
    vigra_fail("binarySphereREAL is untested");
    REAL kernelRadius_um = radius_um;// + gaussWidthAtHalfMaximum_um*3;
    REAL radiusLev = kernelRadius_um /voxelSize[0] + gaussWidthAtHalfMaximum_um*3 ;
    REAL radiusRow = kernelRadius_um /voxelSize[1] + gaussWidthAtHalfMaximum_um*3;
    REAL radiusCol = kernelRadius_um /voxelSize[2] + gaussWidthAtHalfMaximum_um*3;

    int intRadiusLev = (int)std::ceil( radiusLev);
    int intRadiusRow = (int)std::ceil( radiusRow);
    int intRadiusCol = (int)std::ceil( radiusCol);

    MultiArrayShape<3>::type outshape(intRadiusLev*2 + 1, intRadiusRow*2 + 1, intRadiusCol*2 + 1);
    MultiArray<3,REAL> output( outshape);


    for( int m = 0; m < outshape[0]; ++m)
    {
        REAL z_um = (m - intRadiusLev) *voxelSize[0];
        REAL sqr_z_um = z_um * z_um;
        for( int r = 0; r < outshape[1]; ++r)
        {
            REAL y_um = (r - intRadiusRow) *voxelSize[1];
            REAL sqr_y_um = y_um * y_um;
            for( int c = 0; c < outshape[2]; ++c)
            {
                REAL x_um = (c - intRadiusCol) *voxelSize[2];
                REAL sqr_x_um = x_um * x_um;
                REAL dist_um = sqrt( sqr_z_um + sqr_y_um + sqr_x_um);

                if( fabs(dist_um - radius_um) <voxelSize[2]/2)
                {
                    output(m,r,c) = 1;
                }
                else
                {
                    output(m,r,c) = 0;
                }
            }
        }
    }
    return output;
}

/** ------------------------------------------------------
        sphereSurfHarmonic
----------------------------------------------------------
\brief computes a Spherical harmonic base function

\pqrqm out returns 3D SH base function
\param sphereRadius_um radius of the base function
\param gaussWidthAtHalfMaximum_um gaussian smothing of the spherical surface
\param l expansion band, l =[0,l_max]
\param m expansion sub-band, m=[-l,l]
\param full bool, if true, the volume of the sphere is filled, otherwise only surface function  
\param voxelsize optional parameter used to compute base functions for non-equdistant volume samplings
*/

namespace detail {

template <class REAL>
inline
void avoidNans(REAL & temp, REAL eps = 0.00000001)
{
    if (temp == 1.0)
        temp -= eps;
    if (temp == -1.0) 
        temp += eps;
}

} // namespace detail

template <typename REAL>
void 
sphereSurfHarmonic(MultiArray<3,FFTWComplex<REAL> >& output, 
                   REAL sphereRadius_um, REAL gaussWidthAtHalfMaximum_um, 
                   int l, int m, 
                   bool full, 
                   TinyVector<REAL,3> voxelSize=TinyVector<REAL,3>(1.0))
{
    if (gaussWidthAtHalfMaximum_um <=1)
        gaussWidthAtHalfMaximum_um =1;

    REAL radiusLev = sphereRadius_um /voxelSize[0];
    REAL radiusRow = sphereRadius_um /voxelSize[1];
    REAL radiusCol = sphereRadius_um /voxelSize[2];

    radiusLev += gaussWidthAtHalfMaximum_um*3;
    radiusRow += gaussWidthAtHalfMaximum_um*3;
    radiusCol += gaussWidthAtHalfMaximum_um*3;

    int intRadiusLev = (int)std::ceil( radiusLev);
    int intRadiusRow = (int)std::ceil( radiusRow);
    int intRadiusCol = (int)std::ceil( radiusCol);

    MultiArrayShape<3>::type outshape(intRadiusLev*2 + 1,intRadiusRow*2 + 1, intRadiusCol*2 + 1);
    output.reshape( outshape, 0);

    REAL sigmaFactor = REAL(-2.0*std::log(0.5) / (gaussWidthAtHalfMaximum_um * gaussWidthAtHalfMaximum_um));
    for( int s = 0; s < outshape[0]; ++s)
    {
        REAL z_um = (s - intRadiusLev) *voxelSize[0];
        REAL sqr_z_um = z_um * z_um;
        for( int r = 0; r < outshape[1]; ++r)
        {
            REAL y_um = (r - intRadiusRow) *voxelSize[1];
            REAL sqr_y_um = y_um * y_um;
            for( int c = 0; c < outshape[2]; ++c)
            {
                REAL x_um = (c - intRadiusCol) *voxelSize[2];
                REAL sqr_x_um = x_um * x_um;
                REAL dist_um = sqrt( sqr_z_um + sqr_y_um + sqr_x_um);
                REAL gauss_x = 0;
                if (!full||(dist_um > sphereRadius_um))
                {
                    gauss_x=(dist_um - sphereRadius_um);
                }
                else
                {
                    gauss_x=1;
                }
                if( x_um*x_um+y_um*y_um == 0)
                {
                    y_um += REAL(0.00001); //avoid nans
                }
                REAL theta;
                REAL temp = z_um/sqrt( (REAL) (x_um*x_um+y_um*y_um+z_um*z_um));
                detail::avoidNans(temp);
                theta = std::acos(temp);

                REAL phi;
                if (y_um>=0)
                {
                    REAL temp = x_um/sqrt( (REAL)(x_um*x_um+y_um*y_um));
                    detail::avoidNans(temp);
                    phi = std::acos(temp);
                }
                else
                {
                    REAL temp = x_um/sqrt( (REAL)(x_um*x_um+y_um*y_um));
                    detail::avoidNans(temp);
                    phi = REAL(2.0*M_PI - std::acos(temp));
                }

                FFTWComplex<REAL> SHfactor;
                SHfactor.real()= REAL(detail::realSH(l,m)*legendre(l,m,std::cos(theta)) * std::cos(m * phi));
                SHfactor.imag()= REAL(detail::realSH(l,m)*legendre(l,m,std::cos(theta)) * std::sin(m * phi));
                output(s,r,c) = ((REAL)std::exp( -0.5 * gauss_x * gauss_x * sigmaFactor)) * SHfactor;
            }
        }
    }
}

template <typename REAL>
void 
sphereVecHarmonic(MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >  & res,
                  REAL radius, REAL gauss, 
                  int l, int k, int m)
{
    vigra_fail("sphereVecHarmonic is untested");
    using namespace multi_math;
    
    //std::cerr<<"computing VH: "<<l<<" "<<k<<" "<<m<<"\n";
    //1-m
    MultiArray<3,FFTWComplex<REAL> > tmpSH;

    InvariantViolation err("");
    try
    {
        if (abs(1-m)>l) 
            throw err;
        FFTWComplex<REAL> cg(clebschGordan(l+k, m, l, 1-m, 1, 1));

        sphereSurfHarmonic(tmpSH,radius, gauss, l, 1-m, false);        
        // res.bindElementChannel(0) = cg * tmpSH;
        typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >::iterator p = res.begin();
        typename MultiArray<3,FFTWComplex<REAL> >::iterator q = tmpSH.begin();
        for (;q!=tmpSH.end();++p,++q)
            (*p)[0] = cg * *q;
    }
    catch(InvariantViolation &) //in case clebsh gordan dilivers invalid combination
    {
        // res.bindElementChannel(0).init(0.0);
        
        // //std::cerr<<"no Z\n";
        typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >::iterator p = res.begin();
        for (; p!=res.end();++p)
            (*p)[0] = (REAL)0;
    }
    //-m
    try
    {
        if (abs(-m)>l) 
            throw err;
        FFTWComplex<REAL> cg(clebschGordan(l+k, m, l, -m, 1, 0));
        sphereSurfHarmonic(tmpSH,radius, gauss, l, 1-m, false);
        // res.bindElementChannel(1) = cg * tmpSH;
        typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >::iterator p = res.begin();
        typename MultiArray<3,FFTWComplex<REAL> >::iterator q = tmpSH.begin();
        for (;q!=tmpSH.end();++p,++q)
            (*p)[1] = cg * *q;
    }
    catch(InvariantViolation &) //in case clebsh gordan dilivers invalid combination
    {
        // res.bindElementChannel(1).init(0.0);
        
        //std::cerr<<"no Y\n";
        typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >::iterator p = res.begin();
        for (; p!=res.end();++p)
            (*p)[1] = (REAL)0;

    }
    //-(1+m)
    try
    {
        if (abs(-(1+m))>l) 
            throw err;
        FFTWComplex<REAL> cg(clebschGordan(l+k, m, l, -(1+m), 1, -1));
        sphereSurfHarmonic(tmpSH,radius, gauss, l, -(m+1), false);
        // res.bindElementChannel(2) = cg * tmpSH;
        
        typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >::iterator p = res.begin();
        typename MultiArray<3,FFTWComplex<REAL> >::iterator q = tmpSH.begin();
        for (;q!=tmpSH.end();++p,++q)
            (*p)[0] = cg * *q;
    }
    catch(InvariantViolation &) //in case clebsh gordan dilivers invalid combination
    {
        // res.bindElementChannel(2).init(0.0);
        
        //std::cerr<<"no X\n";
        typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >::iterator p = res.begin();
        for (; p!=res.end();++p)
            (*p)[2] = (REAL)0;
    }
}

template <typename REAL>
inline REAL bessel_zero_Jnu(unsigned int l, unsigned int n)
{
    vigra_invariant(l <= 10 && n <= 10,
        "bessel_zero_Jnu(): max implemented band is 10.");

    //cache of the first 10 bessel zero points (max l=10)
    static REAL cache[110] = {
        0,0,0,0,0,0,0,0,0,0, 
        REAL(2.4048255576957729), REAL(3.8317059702075125), 
        REAL(5.1356223018406828), REAL(6.3801618959239841), 
        REAL(7.5883424345038035), REAL(8.7714838159599537), 
        REAL(9.9361095242176845), REAL(11.086370019245084), 
        REAL(12.225092264004656), REAL(13.354300477435331), 
        REAL(5.5200781102863106), REAL(7.0155866698156188), 
        REAL(8.4172441403998643), REAL(9.7610231299816697), 
        REAL(11.064709488501185), REAL(12.338604197466944), 
        REAL(13.589290170541217), REAL(14.821268727013171), 
        REAL(16.03777419088771), REAL(17.241220382489129), 
        REAL(8.6537279129110125), REAL(10.173468135062722), 
        REAL(11.61984117214906), REAL(13.015200721698434), 
        REAL(14.37253667161759), REAL(15.700174079711671), 
        REAL(17.003819667816014), REAL(18.287582832481728), 
        REAL(19.554536430997054), REAL(20.807047789264107), 
        REAL(11.791534439014281), REAL(13.323691936314223), 
        REAL(14.795951782351262), REAL(16.223466160318768), 
        REAL(17.615966049804832), REAL(18.98013387517992), 
        REAL(20.320789213566506), REAL(21.6415410198484), 
        REAL(22.945173131874618), REAL(24.233885257750551), 
        REAL(14.930917708487787), REAL(16.470630050877634), 
        REAL(17.959819494987826), REAL(19.409415226435012), 
        REAL(20.826932956962388), REAL(22.217799896561267), 
        REAL(23.586084435581391), REAL(24.934927887673023), 
        REAL(26.266814641176644), REAL(27.583748963573008), 
        REAL(18.071063967910924), REAL(19.615858510468243), 
        REAL(21.116997053021844), REAL(22.582729593104443), 
        REAL(24.01901952477111), REAL(25.430341154222702), 
        REAL(26.820151983411403), REAL(28.1911884594832), 
        REAL(29.54565967099855), REAL(30.885378967696674), 
        REAL(21.211636629879258), REAL(22.760084380592772), 
        REAL(24.270112313573105), REAL(25.748166699294977), 
        REAL(27.19908776598125), REAL(28.626618307291139), 
        REAL(30.033722386570467), REAL(31.422794192265581), 
        REAL(32.795800037341465), REAL(34.154377923855094), 
        REAL(24.352471530749302), REAL(25.903672087618382), 
        REAL(27.420573549984557), REAL(28.908350780921758), 
        REAL(30.371007667117247), REAL(31.811716724047763), 
        REAL(33.233041762847122), REAL(34.637089352069324), 
        REAL(36.025615063869573), REAL(37.400099977156586), 
        REAL(27.493479132040257), REAL(29.046828534916855), 
        REAL(30.569204495516395), REAL(32.06485240709771), 
        REAL(33.53713771181922), REAL(34.988781294559296), 
        REAL(36.422019668258457), REAL(37.838717382853609), 
        REAL(39.240447995178137), REAL(40.628553718964525), 
        REAL(30.634606468431976), REAL(32.189679910974405), 
        REAL(33.716519509222699), REAL(35.218670738610115), 
        REAL(36.699001128744648), REAL(38.15986856196713), 
        REAL(39.603239416075404), REAL(41.030773691585537), 
        REAL(42.443887743273557), REAL(43.84380142033735) 
    };
    return cache[n*10+l];
}

template <typename REAL>
void 
sphereFullHarmonic(MultiArray<3,FFTWComplex<REAL> >& output, 
                   REAL sphereRadius_um, 
                   int n, int l, int m, 
                   TinyVector<REAL,3> voxelSize=TinyVector<REAL,3>(1.0))
{
    REAL radiusLev = sphereRadius_um / voxelSize[0] +3;
    REAL radiusRow = sphereRadius_um / voxelSize[1] +3;
    REAL radiusCol = sphereRadius_um / voxelSize[2] +3;
    int intRadiusLev = (int)std::ceil( radiusLev);
    int intRadiusRow = (int)std::ceil( radiusRow);
    int intRadiusCol = (int)std::ceil( radiusCol);


    //precompute SH parts for l and m
    MultiArray<3,FFTWComplex<REAL> > SH;
    sphereSurfHarmonic<REAL>(SH,sphereRadius_um, 1, l, m, true);

    MultiArrayShape<3>::type outshape(intRadiusLev*2 + 1,intRadiusRow*2 + 1, intRadiusCol*2 + 1);
    output.reshape( outshape, 0);

    REAL xnl=sphereRadius_um;
    if(n>0)
    {
        xnl=bessel_zero_Jnu<REAL>(l,n);
    }
    REAL k=xnl/sphereRadius_um;
    REAL J2=(REAL)std::pow(besselJ(l+1,xnl),2.0);
    REAL N=(sphereRadius_um*sphereRadius_um*sphereRadius_um)/2*J2;

    REAL sigmaFactor = REAL(-2.0*std::log(0.5) / 4.0);
    FFTWComplex<REAL> I(0,1);
    for( int z = 0; z < outshape[0]; ++z)
    {
        REAL z_um = (z - intRadiusLev) * voxelSize[0];
        REAL sqr_z_um = z_um * z_um;
        for( int y = 0; y < outshape[1]; ++y)
        {
            REAL y_um = (y - intRadiusRow) * voxelSize[1];
            REAL sqr_y_um = y_um * y_um;
            for( int x = 0; x < outshape[2]; ++x)
            {
                REAL x_um = (x - intRadiusCol) * voxelSize[2];
                REAL sqr_x_um = x_um * x_um;
                REAL r = sqrt( sqr_z_um+sqr_y_um + sqr_x_um);
                REAL gauss_x = (r - sphereRadius_um);

                FFTWComplex<REAL> Phi = SH(z,y,x);
                REAL J1=(REAL)besselJ(l,k*r);
                FFTWComplex<REAL> R = 1/sqrt(N)*J1;
                FFTWComplex<REAL> Psi = R*Phi;

                if(r<=sphereRadius_um)
                    output(z,y,x) = Psi;
                else
                    output(z,y,x) = (REAL)exp( -0.5 * gauss_x * gauss_x * sigmaFactor) * Psi;
            }
        }
    }
}

template <typename REAL>
void 
sphereFullVecHarmonic(MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >& output, 
                      REAL sphereRadius_um, 
                      int n, int l, int k, int m, 
                      TinyVector<REAL,3> voxelSize=TinyVector<REAL,3>(1.0))
{
    vigra_fail("sphereFullVecHarmonic is untested");
    REAL radiusLev = sphereRadius_um / voxelSize[0] +3;
    REAL radiusRow = sphereRadius_um / voxelSize[1] +3;
    REAL radiusCol = sphereRadius_um / voxelSize[2] +3;
    int intRadiusLev = (int)std::ceil( radiusLev);
    int intRadiusRow = (int)std::ceil( radiusRow);
    int intRadiusCol = (int)std::ceil( radiusCol);


    //precompute VH parts for l and m
    MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > > VH;
    sphereVecHarmonic(VH, sphereRadius_um, 1, l, k, m);

    MultiArrayShape<3>::type outshape(intRadiusLev*2 + 1,intRadiusRow*2 + 1, intRadiusCol*2 + 1);
    TinyVector<FFTWComplex<REAL>,3 > zero(0,0,0);
    output.reshape( outshape, zero);

    REAL xnl=sphereRadius_um;
    if(n>0)
    {
            xnl=0;//gsl_sf_bessel_zero_Jnu(l,n);
            //std::cerr<<n<<" "<<m<<" "<<xnm<<"\n";
    }
    REAL K=xnl/sphereRadius_um;
    REAL J2=std::pow(besselJ(l+1,xnl),2.0);
    REAL N=(sphereRadius_um*sphereRadius_um*sphereRadius_um)/2*J2;

    REAL sigmaFactor = -2*log(0.5) / (4);
    FFTWComplex<REAL> I(0,1);
    for( int z = 0; z < outshape[0]; ++z)
    {
        REAL z_um = (z - intRadiusLev) * voxelSize[0];
        REAL sqr_z_um = z_um * z_um;
        for( int y = 0; y < outshape[1]; ++y)
        {
            REAL y_um = (y - intRadiusRow) * voxelSize[1];
            REAL sqr_y_um = y_um * y_um;
            for( int x = 0; x < outshape[2]; ++x)
            {
                REAL x_um = (x - intRadiusCol) * voxelSize[2];
                REAL sqr_x_um = x_um * x_um;
                REAL r = sqrt( sqr_z_um+sqr_y_um + sqr_x_um);
                REAL gauss_x = (r - sphereRadius_um);
                TinyVector<FFTWComplex<REAL>,3> Phi = VH(z,y,x);
                REAL J1=besselJ(l,K*r);
                FFTWComplex<REAL> R = 1/sqrt(N)*J1;
                TinyVector<FFTWComplex<REAL>,3> Psi;
                Psi[0]= R*Phi[0];
                Psi[1]= R*Phi[1];
                Psi[2]= R*Phi[2];

                if(r<=sphereRadius_um)
                {
                    output(z,y,x) = Psi;
                }
                else
                {
                    Psi[0]*=((FFTWComplex<REAL>)exp( -0.5 * gauss_x * gauss_x * sigmaFactor));
                    Psi[1]*=((FFTWComplex<REAL>)exp( -0.5 * gauss_x * gauss_x * sigmaFactor));
                    Psi[2]*=((FFTWComplex<REAL>)exp( -0.5 * gauss_x * gauss_x * sigmaFactor));
                    output(z,y,x) =  Psi;
                }
            }
        }
    }
}

/** ------------------------------------------------------
        computeSHbaseF
----------------------------------------------------------
\brief pre-computes Spherical harmonic base functions

\param radius radius of the spherical expansion
\param gauss smothing of the spherical surface
\param band maximum expansion band
\param SHbaseF holds precomputed SH base functions  
*/
template <typename REAL>
void 
computeSHbaseF(REAL radius, REAL gauss, unsigned int band, 
               std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > &SHbaseF)
{
    SHbaseF.resize(band + 1);

    for (int l = 0; l <= (int)band; l++)
    {
        SHbaseF[l].resize(2*l + 1);
        for (int m = -l; m <= l; m++)
        {
            MultiArray<3, FFTWComplex<REAL> > Coeff;
            sphereSurfHarmonic(Coeff,radius, gauss, l, m, false);
            SHbaseF[l][m+l].reshape(Coeff.shape());
            SHbaseF[l][m+l]=Coeff;
        }
    }
}

template <typename REAL>
void 
computePHbaseF(REAL radius, unsigned int band, 
               std::vector<std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > > &PHbaseF, 
               bool realdata=false)
{
    PHbaseF.resize(band + 1);
    //n=0 undefined
    int mMin;
    int mOff;
    //for real data one only needs positive coeffs (symmetry)
    for (int n = 1; n <= (int)band; n++)
    {
        PHbaseF[n].resize(band + 1);
        for (int l = 0; l <= (int)band; l++)
        {
            if (realdata)
            {
                mMin = 0;
                mOff = 0;
                PHbaseF[n][l].resize(l + 1);
            }
            else
            {
                mMin = -l;
                mOff = l;
                PHbaseF[n][l].resize(2*l + 1);
            }
            for (int m = mMin; m <= l; m++)
            {
                //std::cerr<<n<<" "<<l<<" "<<m+mOff<<"\n";
                MultiArray< 3,FFTWComplex<REAL> > Coeff;
                sphereFullHarmonic(Coeff, radius, n,l,m);
                PHbaseF[n][l][m+mOff].reshape(Coeff.shape());
                PHbaseF[n][l][m+mOff]=Coeff;
            }
        }
    }
}

template <typename REAL>
void 
computeVHbaseF(REAL radius, REAL gauss, unsigned int band, 
               std::vector<std::vector<std::vector<MultiArray<3, TinyVector<FFTWComplex<REAL>,3> > > > > &VHbaseF)
{
    vigra_fail("computeVHbaseF is untested");
    VHbaseF.resize(band + 1);

    for (int l = 0; l <= band; l++)
    {
        VHbaseF[l].resize(3);
        for (int k=-1;k<=1;k++)
        {
            VHbaseF[l][k+1].resize(2*band + 1);
            for (int m = -l; m <= l; m++)
            {
                MultiArray<3, TinyVector< FFTWComplex<REAL> ,3 > > Coeff;

                sphereVecHarmonic(Coeff,radius, gauss, l, k, m);
                VHbaseF[l][k+1][m+l].reshape(Coeff.shape());
                VHbaseF[l][k+1][m+l]=Coeff;
            }
        }
    }
}

template <typename REAL>
void 
computeVPHbaseF(REAL radius, unsigned int band, 
                std::vector<std::vector<std::vector<std::vector<MultiArray<3, TinyVector< FFTWComplex<REAL>,3 > > > > > >&VHbaseF)
{
    vigra_fail("computeVPHbaseF is untested");
    VHbaseF.resize(band + 1);
    for(int n=0; n>= band; n++)
    {
        VHbaseF[n].resize(band + 1);

        for (int l = 0; l <= band; l++)
        {
            VHbaseF[n][l].resize(3);
            for (int k=-1;k<=1;k++)
            {
                VHbaseF[l][k+1].resize(2*band + 1);
                for (int m = -l; m <= l; m++)
                {
                    MultiArray<3, TinyVector< FFTWComplex<REAL> ,3 > > Coeff;

                    sphereFullVecHarmonic(Coeff,radius, n, l, k, m);
                    VHbaseF[n][l][k+1][m+l].reshape(Coeff.shape());
                    VHbaseF[n][l][k+1][m+l]=Coeff;
                }
            }
        }
    }
}


template <typename REAL>
MultiArray<3,REAL> 
sphereSurfGauss(REAL sphereRadius_um, REAL gaussWidthAtHalfMaximum_um, 
                TinyVector<REAL,3> voxelSize=TinyVector<REAL,3>(1.0))
{
    vigra_fail("sphereSurfGauss is untested");
    REAL kernelRadius_um = sphereRadius_um;;
    REAL radiusLev = kernelRadius_um /voxelSize[0] + gaussWidthAtHalfMaximum_um*3;
    REAL radiusRow = kernelRadius_um /voxelSize[1] + gaussWidthAtHalfMaximum_um*3;
    REAL radiusCol = kernelRadius_um /voxelSize[2] + gaussWidthAtHalfMaximum_um*3;

    int intRadiusLev = (int)std::ceil( radiusLev);
    int intRadiusRow = (int)std::ceil( radiusRow);
    int intRadiusCol = (int)std::ceil( radiusCol);

    MultiArrayShape<3>::type outShape(intRadiusLev*2 + 1, intRadiusRow*2 + 1, intRadiusCol*2 + 1);
    MultiArray<3,REAL> output( outShape );

    REAL sigmaFactor = -2*log(0.5) / (gaussWidthAtHalfMaximum_um * gaussWidthAtHalfMaximum_um);

    for( int m = 0; m < outShape[0]; ++m)
    {
        REAL z_um = (m - intRadiusLev) *voxelSize[0];
        REAL sqr_z_um = z_um * z_um;
        for( int r = 0; r < outShape[1]; ++r)
        {
            REAL y_um = (r - intRadiusRow) *voxelSize[1];
            REAL sqr_y_um = y_um * y_um;
            for( int c = 0; c < outShape[2]; ++c)
            {
                REAL x_um = (c - intRadiusCol) *voxelSize[2];
                REAL sqr_x_um = x_um * x_um;
                REAL dist_um = sqrt( sqr_z_um + sqr_y_um + sqr_x_um);
                REAL gauss_x = (dist_um - sphereRadius_um);
                output(m,r,c) = exp( -0.5 * gauss_x * gauss_x * sigmaFactor);
            }
        }
    }

    REAL kernelSum = 0;
    typename MultiArray<3,REAL>::iterator p=output.begin();
    for (;p!=output.end();++p)
        kernelSum+=*p;
    output *= 1.0/kernelSum;

    return output;
}

template <typename REAL>
void 
reconstSH(REAL radius, REAL gauss, unsigned int band, 
          MultiArray<3,REAL >& reconstruct, 
          const std::vector<std::vector<FFTWComplex<REAL> >  >& SH_A, 
          const std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > >& SHbaseF)
{
    using namespace multi_math;
    
    //FIXME reconstSH currently only works for real data
    reconstruct.reshape(SHbaseF[0][0].shape());

    MultiArray<3,FFTWComplex<REAL> > T;
    std::string name;
    for (int l=0;l<=(int)band;l++)
    {
        for (int m=-l;m<=l;m++)
        {
            // reconstruct += real(SHbaseF[l][l+m] * conj(SH_A[l][l+m]));
            typename MultiArray<3,FFTWComplex<REAL> >::iterator p=SHbaseF[l][l+m].begin();
            typename MultiArray<3,REAL >::iterator q=reconstruct.begin();
            for (;q!=reconstruct.end();++p,++q)
                *q += real(*p * conj(SH_A[l][l+m]));
        }
    }
}

template <typename REAL>
void 
reconstPH(REAL radius, unsigned int band, 
          MultiArray<3, REAL >& reconstruct, 
          const std::vector<std::vector<std::vector<FFTWComplex<REAL> >  > >& PH_A, 
          std::vector<std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > > &PHbaseF)
{
    using namespace multi_math;
    
    //FIXME reconstPH currently only works for real data
    reconstruct.reshape(PHbaseF[1][0][0].shape(), 0);

    for (int n=1;n<(int)PH_A.size();n++)
    {
        for (int l=0;l<(int)PH_A[n].size();l++)
        {
            for (int m=0;m<(int)PH_A[n][l].size();m++)
            {
                // reconstruct += real(PHbaseF[n][l][m] * conj(PH_A[n][l][m]));
                typename MultiArray<3,FFTWComplex<REAL> >::iterator p = PHbaseF[n][l][m].begin();
                typename MultiArray<3, REAL >::iterator q=reconstruct.begin();
                for (;q!=reconstruct.end();++q,++p)
                    *q += real(*p * conj(PH_A[n][l][m]));
            }
        }
    }
}

template <typename REAL>
void 
reconstVPH(REAL radius, unsigned int band, 
           MultiArray<3, TinyVector<REAL,3> >& reconstruct, 
           const std::vector<std::vector<std::vector<std::vector<FFTWComplex<REAL> > > > >& VPH_A)
{
    vigra_fail("reconstVPH is untested");
    using namespace multi_math;
    
    MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > > base_tmp;
    sphereVecHarmonic(base_tmp, radius, 1, 0, 0, 0);
    reconstruct.reshape(base_tmp.shape());
    MultiArray<3, TinyVector<FFTWComplex<REAL>,3> > tmp(reconstruct.shape(), 0);
    FFTWComplex<REAL> zero(0,0);
    TinyVector<FFTWComplex<REAL>,3> Zero(zero,zero,zero);

    for (int n=1;n<=(int)band;n++)
    {
        for (int l=0;l<=(int)band;l++)
        {
            for (int k=-1;k<=1;++k)
            {
                for (int m=-(l+k);m<=(l+k);m++)
                {
                    sphereVecHarmonic(base_tmp, radius, n, l, k, m);
                    
                    // tmp += conj(VPH_A[n][l][k+1][(l+k)+m]) * base_tmp;
                    typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3> >::iterator p=tmp.begin();
                    typename MultiArray< 3, TinyVector<FFTWComplex<REAL>,3 > >::iterator q=base_tmp.begin();
                    for (;q!=base_tmp.end();++p,++q)
                    {
                        (*p)[0] += conj(VPH_A[n][l][k+1][(l+k)+m]) * (*q)[0];
                        (*p)[1] += conj(VPH_A[n][l][k+1][(l+k)+m]) * (*q)[1];
                        (*p)[2] += conj(VPH_A[n][l][k+1][(l+k)+m]) * (*q)[2];
                    }
                }
            }
        }
    }
    
    // reconstruct.bindElementChannel(0) = real(tmp.bindElementChannel(1));
    // reconstruct.bindElementChannel(1) = -sqrt(2.0)*real(tmp.bindElementChannel(0));
    // reconstruct.bindElementChannel(1) = sqrt(2.0)*imag(tmp.bindElementChannel(0));
    //reconst real vec directions
    typename MultiArray<3, TinyVector<REAL,3> >::iterator p=reconstruct.begin();
    typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3> >::iterator q=tmp.begin();
    for (;q!=tmp.end();++p,++q)
    {
        (*p)[0] = real((*q)[1]);
        (*p)[1] = -sqrt(2)*real((*q)[0]);
        (*p)[2] = sqrt(2)*imag((*q)[0]);
    }
}

template <typename REAL>
void 
reconstVH(REAL radius, REAL gauss, unsigned int band, 
          MultiArray<3, TinyVector<REAL,3> >& reconstruct,
          const std::vector<std::vector<std::vector<FFTWComplex<REAL> > > >& VH_A)
{
    vigra_fail("reconstVH is untested");
    using namespace multi_math;
    
    MultiArray< 3, TinyVector<FFTWComplex<REAL>,3 > > base_tmp;
    sphereVecHarmonic(base_tmp, radius, gauss, 0, 0, 0);
    reconstruct.reshape(base_tmp.shape());
    FFTWComplex<REAL> zero(0,0);
    MultiArray<3, TinyVector<FFTWComplex<REAL>,3> > tmp(reconstruct.shape(), 0);
    TinyVector<FFTWComplex<REAL>,3> Zero(zero,zero,zero);

    for (int l=0;l<=band;l++)
    {
        for (int k=-1;k<=1;++k)
        {
            for (int m=-(l+k);m<=(l+k);m++)
            {
                sphereVecHarmonic(base_tmp, radius, gauss, l, k, m);
                // tmp += conj(VH_A[l][k+1][(l+k)+m]) * base_tmp;
                
                typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3> >::iterator p=tmp.begin();
                typename MultiArray< 3, TinyVector<FFTWComplex<REAL>,3 > >::iterator q=base_tmp.begin();
                for (;q!=base_tmp.end();++p,++q)
                {
                    (*p)[0] += conj(VH_A[l][k+1][(l+k)+m]) * (*q)[0];
                    (*p)[1] += conj(VH_A[l][k+1][(l+k)+m]) * (*q)[1];
                    (*p)[2] += conj(VH_A[l][k+1][(l+k)+m]) * (*q)[2];
                }
            }
        }
    }

    // reconstruct.bindElementChannel(0) = real(tmp.bindElementChannel(1));
    // reconstruct.bindElementChannel(1) = -sqrt(2.0)*real(tmp.bindElementChannel(0));
    // reconstruct.bindElementChannel(1) = sqrt(2.0)*imag(tmp.bindElementChannel(0));
    //reconst real vec directions
    typename MultiArray<3, TinyVector<REAL,3> >::iterator p=reconstruct.begin();
    typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3> >::iterator q=tmp.begin();
    for (;q!=tmp.end();++p,++q)
    {
        (*p)[0] = real((*q)[1]);
        (*p)[1] = -sqrt(2)*real((*q)[0]);
        (*p)[2] = sqrt(2)*imag((*q)[0]);
    }
}

/** ------------------------------------------------------
    SHpos
----------------------------------------------------------
\brief computes singe local Spherical harmonic expansion at given 3D position

\param SH_A holds the returned SH coefficients  
\param radius radius of the spherical expansion
\param band maximum expansion band
\param A input volume data
\param pos 3D position of the SH expansion
\param SHbaseF precomputed SH base functions (-> see computeSHbaseF) 
*/
template <typename REAL>
void 
SHpos(std::vector<std::vector<FFTWComplex<REAL> >  > &SH_A, 
      REAL radius, unsigned int band,
      const MultiArray<3,REAL> &A, 
      const TinyVector<REAL, 3> &pos, 
      std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > const &SHbaseF)
{
    SH_A.resize(SHbaseF.size());

    // SH coeffs
    for (int l = 0; l <= (int)band; l++)
    {
        SH_A[l].resize(SHbaseF[l].size());
        for (int m = 0; m < (int)SHbaseF[l].size(); m++)
        {
            MultiArrayShape<3>::type coffShape(SHbaseF[l][m].shape());

            int xa = (int) floor(pos[2] - coffShape[2] / 2);
            int xe = xa + coffShape[2] - 1;
            int ya = (int) floor(pos[1] - coffShape[1] / 2);
            int ye = ya + coffShape[1] - 1;
            int za = (int) floor(pos[0] - coffShape[0] / 2);
            int ze = za + coffShape[0] - 1;

            SH_A[l][m] = (REAL)0;
            int sz=0;
            for (int z=ze;z>=za;z--,sz++)
            {
                int sy=0;
                for (int y=ye;y>=ya;y--,sy++)
                {
                    int sx=0;
                    for (int x=xe;x>=xa;x--,sx++)
                    SH_A[l][m]+=A(z,y,x)*(SHbaseF[l][m](sz,sy,sx));
                }
            }
        }
    }
}

/** ------------------------------------------------------
        SHcenter
----------------------------------------------------------
\brief computes singe local Spherical harmonic expansion at the center of te given 3D volume

\param SH_A holds the returned SH coefficients  
\param radius radius of the spherical expansion
\param band maximum expansion band
\param A input volume data
\param SHbaseF precomputed SH base functions (-> see computeSHbaseF) 
*/

template <typename REAL>
void 
SHcenter(std::vector<std::vector<FFTWComplex<REAL> >  > &SH_A, 
         REAL radius, unsigned int band,
         const MultiArray<3,REAL> &A, 
         std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > const &SHbaseF)
{
    SHpos(SH_A, radius, band, A,  detail::centerOfBB<REAL>(A), SHbaseF);
}

/** ------------------------------------------------------
        SH_Iterator
---------------------------------------------------------- 
\brief iterator class to access the cascaded vector representations of SH base functions and SH coefficients
*/
template <typename T>
class SH_Iterator
{
  public:
    typedef MultiArray<3,T >* pointer;
    typedef MultiArray<3,T >& reference;
    typedef MultiArray<3,T > value_type;
    typedef int difference_type;
    typedef std::random_access_iterator_tag iterator_category;

    SH_Iterator(std::vector<std::vector<MultiArray<3,T > > >& data, 
                 int l, int m, std::string name="name")
    :   _data(data),
        _l(l),
        _m(m),
        _name(name)
    {}
    
    void operator++()
    {
        if (_m< (int)_data[_l].size()-1)
        {
            //std::cerr<<_name<<"++ "<<_m<<" "<<_data[_l].size()<<"\n"<<std::flush;
            _m++;
        }
        else
        {
            if (_l< (int)_data.size())
            {
                _m=0;
                _l++;
            }
        }
    }
    
    MultiArray<3,T >& operator*() const
    {
        return _data[_l][_m];
    }

    const MultiArray<3,T >* operator->() const
    {
        //std::cerr<<_name<<" "<<_l<<" "<<_m<<" "<<_data[_l].size()<<" "<<_data[_l][_m].shape()<<"\n";
        return &_data[_l][_m];
    }

    int getL()
    {
        return _l;
    }
    
    int getM()
    {
        return _m;
    }
    
    bool operator==(SH_Iterator& A)
    {
        return ((this->_l==A.getL())&&(this->_m==A.getM()));
    }
    
    bool operator!=(SH_Iterator& A)
    {
        return !((this->_l==A.getL())&&(this->_m==A.getM()));
    }

  private:
    std::vector<std::vector<MultiArray<3,T > > > &_data;
    int _l;
    int _m;
    std::string _name;
};

template <typename REAL>
class PH_Iterator
{
  public:
    typedef MultiArray<3,FFTWComplex<REAL> >* pointer;
    typedef MultiArray<3,FFTWComplex<REAL> >& reference;
    typedef MultiArray<3,FFTWComplex<REAL> > value_type;
    typedef int difference_type;
    typedef std::random_access_iterator_tag iterator_category;

    PH_Iterator(std::vector<std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > >& data, 
                int k, int l, int m)
    : _data(data),
      _k(k),
      _l(l),
      _m(m)
    {}
    
    void operator++()
    {
        if (_m<_data[_k][_l].size()-1)
        {
            _m++;
        }
        else
        {
            if (_l<_data[_k].size()-1)
            {
                _m=0;
                _l++;
            }
            else
            {
                if (_k<_data.size())
                {
                    _l=0;
                    _k++;
                }
            }
        }
    }
    
    MultiArray<3,FFTWComplex<REAL> >& operator*() const
    {
        return _data[_k][_l][_m];
    }

    const MultiArray<3,FFTWComplex<REAL> >* operator->() const
    {
        return &_data[_k][_l][_m];
    }
    
    int getK()
    {
        return _k;
    }
    
    int getL()
    {
        return _l;
    }
    
    int getM()
    {
        return _m;
    }
    
    bool operator==(PH_Iterator& A)
    {
        return ((this->_k==A.getL())&&(this->_l==A.getL())&&(this->_m==A.getM()));
    }
    
    bool operator!=(PH_Iterator& A)
    {
        return !((this->_k==A.getK())&&(this->_l==A.getL())&&(this->_m==A.getM()));
    }

  private:
    std::vector<std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > > &_data;
    int _k;
    int _l;
    int _m;
};

/** ------------------------------------------------------
        Array2SH
----------------------------------------------------------
\brief computes  local Spherical harmonic expansion at all positions of the given 3D volume

\param SH_A holds the returned SH coefficients  
\param radius radius of the spherical expansion
\param band maximum expansion band
\param A input volume data
\param SHbaseF precomputed SH base functions (-> see computeSHbaseF) 
*/
template <typename REAL>
void 
Array2SH(std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> >  > > &SH_A, 
         unsigned int band, REAL radius, 
         std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > &SHbaseF,
         const MultiArray<3,REAL> &A)
//FIXME add const iterator and make SHbaseF const
{

    SH_A.resize(SHbaseF.size());
    for (int l=0;l<(int)SHbaseF.size();l++)
    {
        SH_A[l].resize(SHbaseF[l].size());
        for (int m=0;m<(int)SHbaseF[l].size();m++)
        {
            SH_A[l][m].reshape(A.shape(),0);
        }
    }

    SH_Iterator< FFTWComplex<REAL> > SHbaseF_Iter(SHbaseF,0,0,"BaseF");
    SH_Iterator< FFTWComplex<REAL> > SHbaseF_Iter_end(SHbaseF,SHbaseF.size(),0);
    SH_Iterator< FFTWComplex<REAL> > SH_A_Iter(SH_A,0,0,"SH_A");
    convolveFFTComplexMany(A, SHbaseF_Iter, SHbaseF_Iter_end, SH_A_Iter,false);
}

template <typename REAL>
void 
PHpos(std::vector<std::vector<std::vector<FFTWComplex<REAL> >  > > &PH_A, 
      REAL radius, unsigned int band,
      const MultiArray<3,FFTWComplex<REAL> > &A, 
      std::vector<std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > > &PHbaseF,
      const TinyVector<REAL, 3> &pos)
{
    using namespace multi_math;
    PH_A.resize(PHbaseF.size());

    for (int n=1;n<(int)PHbaseF.size();n++)
    {
        PH_A[n].resize(PHbaseF[n].size());
        for (int l = 0; l < (int)PHbaseF[n].size(); l++)
        {
            PH_A[n][l].resize(PHbaseF[n][l].size());
            for (int m = 0; m < (int)PHbaseF[n][l].size(); m++)
            {
                MultiArrayShape<3>::type coffShape(PHbaseF[n][l][m].shape());
                int xa = (int) floor(pos[2] - coffShape[2] / 2);
                int xe = xa + coffShape[2] - 1;
                int ya = (int) floor(pos[1] - coffShape[1] / 2);
                int ye = ya + coffShape[1] - 1;
                int za = (int) floor(pos[0] - coffShape[0] / 2);
                int ze = za + coffShape[0] - 1;
                
                // PH_A[n][l][m] = sum(A.subarray(Shape3(za,ya,xa),Shape3(ze,ye,xe))*PHbaseF[n][l][m],
                                    // FFTWComplex<double>());


                PH_A[n][l][m] = (REAL)0;
                int sz=0;
                int sy=0;
                int sx=0;
                for (int z=za;z<ze;z++,sz++)
                {
                    sy=0;
                    for (int y=ya;y<ye;y++,sy++)
                    {
                        sx=0;
                        for (int x=xa;x<xe;x++,sx++)
                            PH_A[n][l][m]+=A(z,y,x)*PHbaseF[n][l][m](sz,sy,sx);
                    }
                }
            }
        }
    }
}

template <typename REAL>
void 
PHcenter(std::vector<std::vector<std::vector<FFTWComplex<REAL> >  > > &PH_A, 
         REAL radius, unsigned int band,
         const MultiArray<3,FFTWComplex<REAL> > &A, 
         std::vector<std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > > &PHbaseF)
{
    PHpos(PH_A, radius, band, A, PHbaseF, detail::centerOfBB<REAL>(A));
}

template <typename REAL>
void 
Array2PH(std::vector<std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> >  > > > &PH_A, 
         unsigned int band, REAL radius, 
         bool realData, 
         const MultiArray<3,FFTWComplex<REAL> > &A, 
         fftwf_plan forward_plan, fftwf_plan backward_plan)
{
    vigra_fail("Array2PH is untested");
    std::vector<std::vector<std::vector<MultiArray<3,FFTWComplex<REAL> > > > > PHbaseF;
    computePHbaseF(radius, band, PHbaseF, realData);

    PH_A.resize(band+1);

    for (int n=1;n<=band;n++)
    {
        PH_A[n].resize(band+1);
        for (int l=0;l<=band;l++)
            PH_A[n][l].resize(2*band+1);
    }

    PH_Iterator<REAL> PHbaseF_Iter(PHbaseF,0,0,0);
    PH_Iterator<REAL> PHbaseF_Iter_end(PHbaseF,3,PHbaseF.size(),0);
    PH_Iterator<REAL> PH_A_Iter(PH_A,0,0,0);
    convolveFFTComplexMany(A, PHbaseF_Iter, PHbaseF_Iter_end, PH_A_Iter,false);
} 

/*
template <typename REAL>
void Array2VH(std::vector<std::vector<std::vector<MultiArray<3, FFTWComplex<REAL> > > > > &VH_A, unsigned int band, REAL gauss, REAL radius, const MultiArray<3, TinyVector<REAL,3> > &A)
{
    //transform input to C^(2j+1)
    MultiArray< 3,FFTWComplex<REAL> >  inputZ(A.shape());
    MultiArray< 3,FFTWComplex<REAL> >  inputY(A.shape());
    MultiArray< 3,FFTWComplex<REAL> >  inputX(A.shape());
    MultiArray< 3,FFTWComplex<REAL> >::iterator z=inputZ.begin();
    MultiArray< 3,FFTWComplex<REAL> >::iterator y=inputY.begin();
    MultiArray< 3,FFTWComplex<REAL> >::iterator x=inputX.begin();

    for (MultiArray<3,TinyVector<REAL,3> >::const_iterator p=A.begin();p!=A.end();++p,++z,++y,++x)
    {
        z->real() = -(*p)[1];
        z->imag() = -(*p)[2];
        *z *= 1/sqrt(2);
        y->real() =  (*p)[0];
        y->imag() =  0;
        x->real() = (*p)[1];
        x->imag() = -(*p)[2];
        *x *= 1/sqrt(2);
    }
   VH_A.resize(band + 1);
   REAL norm = M_PI/2*pow((REAL)radius,2.0);


    for(int k=-1;k<=1;++k)
    {
        VH_A[0][1+k].resize(3);
        for (int m=-1;m<=1;++m)
        {
            MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > > vh;
            sphereVecHarmonic(vh, radius, gauss, 0, k, m);
            MultiArray<3, FFTWComplex<REAL> > vhz(vh.shape());
            MultiArray<3, FFTWComplex<REAL> > vhy(vh.shape());
            MultiArray<3, FFTWComplex<REAL> > vhx(vh.shape());
            MultiArray<3, FFTWComplex<REAL> >::iterator z = vhz.begin();
            MultiArray<3, FFTWComplex<REAL> >::iterator y = vhy.begin();
            MultiArray<3, FFTWComplex<REAL> >::iterator x = vhx.begin();
            for(MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >::iterator p= vh.begin();p!=vh.end();p++,z++,y++,x++)
            {
                *z = (*p)[0];
                *y = (*p)[1];
                *x = (*p)[2];
            }
            VH_A[0][1+k][1+m].reshape(A.shape());
            VH_A[0][1+k][1+m] = convolveComplex(inputZ, vhz);
            VH_A[0][1+k][1+m] += convolveComplex(inputY, vhy);
            VH_A[0][1+k][1+m] += convolveComplex(inputX, vhx);

            for (MultiArray<3, FFTWComplex<REAL> >::iterator p=VH_A[0][1+k][1+m].begin();p!=VH_A[0][1+k][1+m].end();++p)
                *p*= (complex<REAL>)pow(-1.0,(REAL) 0)/(3*norm);
        }
    }

    for (int l=1;l<=band;l++)
    {
        VH_A[l].resize(3);
        for(int k=-1;k<=1;++k)
        {
            VH_A[l][1+k].resize(2*(l+k)+1);
            for (int m=-(l+k);m<=(l+k);++m)
            {
                MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > > vh;
                sphereVecHarmonic(vh, radius, gauss, 0, k, m);
                MultiArray<3, FFTWComplex<REAL> > vhz(vh.shape());
                MultiArray<3, FFTWComplex<REAL> > vhy(vh.shape());
                MultiArray<3, FFTWComplex<REAL> > vhx(vh.shape());
                MultiArray<3, FFTWComplex<REAL> >::iterator z = vhz.begin();
                MultiArray<3, FFTWComplex<REAL> >::iterator y = vhy.begin();
                MultiArray<3, FFTWComplex<REAL> >::iterator x = vhx.begin();
                for(MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >::iterator p= vh.begin();p!=vh.end();p++,z++,y++,x++)
                {
                    *z = (*p)[0];
                    *y = (*p)[1];
                    *x = (*p)[2];
                }

                VH_A[l][1+k][1+m].reshape(A.shape());


                VH_A[l][1+k][1+m] = convolveComplex(inputZ, vhz);
                VH_A[l][1+k][1+m] += convolveComplex(inputY, vhy);
                VH_A[l][1+k][1+m] += convolveComplex(inputX, vhx);

                for (MultiArray<3, FFTWComplex<REAL> >::iterator p=VH_A[l][1+k][1+m].begin();p!=VH_A[l][1+k][1+m].end();++p)
                    *p*= (complex<REAL>)pow(-1.0,(REAL) 0)/(3*norm);

            }

        }
    }
}

template <typename REAL>
void Array2VPH(std::vector<std::vector<std::vector<std::vector<MultiArray<3, FFTWComplex<REAL> > > > > >&VH_A, unsigned int band, REAL radius, const MultiArray<3, TinyVector<REAL,3> > &A)
{
    //transform input to C^(2j+1)
    MultiArray< 3,FFTWComplex<REAL> >  inputZ(A.shape());
    MultiArray< 3,FFTWComplex<REAL> >  inputY(A.shape());
    MultiArray< 3,FFTWComplex<REAL> >  inputX(A.shape());
    MultiArray< 3,FFTWComplex<REAL> >::iterator z=inputZ.begin();
    MultiArray< 3,FFTWComplex<REAL> >::iterator y=inputY.begin();
    MultiArray< 3,FFTWComplex<REAL> >::iterator x=inputX.begin();

    for (MultiArray<3,TinyVector<REAL,3> >::const_iterator p=A.begin();p!=A.end();++p,++z,++y,++x)
    {
        real(*z) = -(*p)[1];
        imag(*z) = -(*p)[2];
        *z *= 1/sqrt(2);
        real(*y) =  (*p)[0];
        imag(*y) =  0;
        real(*x) = (*p)[1];
        imag(*x) = -(*p)[2];
        *x *= 1/sqrt(2);
    }

    for(int n=1;n<=band;n++)
    {
        VH_A.resize(band+1);
        //0th coeff
        VH_A[n][0].resize(3);

        REAL norm = M_PI/2*pow((REAL)radius,2.0);

        for(int k=-1;k<=1;++k)
        {
            VH_A[n][0][1+k].resize(3);
            for (int m=-1;m<=1;++m)
            {
                MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > > vh;
                sphereFullVecHarmonic(vh, radius,  n, 0, k, m);
                MultiArray<3, FFTWComplex<REAL> > vhz(vh.shape());
                MultiArray<3, FFTWComplex<REAL> > vhy(vh.shape());
                MultiArray<3, FFTWComplex<REAL> > vhx(vh.shape());
                MultiArray<3, FFTWComplex<REAL> >::iterator z = vhz.begin();
                MultiArray<3, FFTWComplex<REAL> >::iterator y = vhy.begin();
                MultiArray<3, FFTWComplex<REAL> >::iterator x = vhx.begin();
                for(MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >::iterator p= vh.begin();p!=vh.end();p++,z++,y++,x++)
                {
                    *z = (*p)[0];
                    *y = (*p)[1];
                    *x = (*p)[2];
                }
                VH_A[n][0][1+k][1+m].reshape(A.shape());
                VH_A[n][0][1+k][1+m] = convolveComplex(inputZ, vhz);
                VH_A[n][0][1+k][1+m] += convolveComplex(inputY, vhy);
                VH_A[n][0][1+k][1+m] += convolveComplex(inputX, vhx);

                for (MultiArray<3, FFTWComplex<REAL> >::iterator p=VH_A[n][0][1+k][1+m].begin();p!=VH_A[n][0][1+k][1+m].end();++p)
                    *p*= (complex<REAL>)pow(-1.0,(REAL) 0)/(3*norm);

            }
        }

        for (int l=1;l<=band;l++)
        {
            VH_A[n][l].resize(3);
            for(int k=-1;k<=1;++k)
            {
                VH_A[n][l][1+k].resize(2*(l+k)+1);
                for (int m=-(l+k);m<=(l+k);++m)
                {
                    MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > > vh;
                    sphereFullVecHarmonic(vh, radius,  n, 0, k, m);
                    MultiArray<3, FFTWComplex<REAL> > vhz(vh.shape());
                    MultiArray<3, FFTWComplex<REAL> > vhy(vh.shape());
                    MultiArray<3, FFTWComplex<REAL> > vhx(vh.shape());
                    typename MultiArray<3, FFTWComplex<REAL> >::iterator z = vhz.begin();
                    typename MultiArray<3, FFTWComplex<REAL> >::iterator y = vhy.begin();
                    typename MultiArray<3, FFTWComplex<REAL> >::iterator x = vhx.begin();
                    typename MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > >::iterator p= vh.begin();
                    for(;p!=vh.end();p++,z++,y++,x++)
                    {
                        *z = (*p)[0];
                        *y = (*p)[1];
                        *x = (*p)[2];
                    }

                    VH_A[n][l][1+k][1+m].reshape(A.shape());


                    VH_A[n][l][1+k][1+m] = convolveComplex(inputZ, vhz);
                    VH_A[n][l][1+k][1+m] += convolveComplex(inputY, vhy);
                    VH_A[n][l][1+k][1+m] += convolveComplex(inputX, vhx);

                    typename MultiArray<3, FFTWComplex<REAL> >::iterator p=VH_A[n][l][1+k][1+m].begin();
                    for (;p!=VH_A[n][l][1+k][1+m].end();++p)
                        *p*= (complex<REAL>)pow(-1.0,(REAL) 0)/(3*norm);

                }
            }
        }
    }
}
*/

template <typename REAL>
void 
VHpos(std::vector<std::vector<std::vector<FFTWComplex<REAL> > > > &VH_A, 
      unsigned int band, REAL gauss, REAL radius, 
      const MultiArray<3,TinyVector<REAL,3> > &A, 
      const TinyVector<REAL, 3> &pos)
{
    using namespace multi_math;
    
    vigra_fail("VHpos is untested");
    
    //transform input to C^(2j+1)
    MultiArray< 3,FFTWComplex<REAL> >  inputZ(A.shape());
    MultiArray< 3,FFTWComplex<REAL> >  inputY(A.shape());
    MultiArray< 3,FFTWComplex<REAL> >  inputX(A.shape());
    typename MultiArray< 3,FFTWComplex<REAL> >::iterator z=inputZ.begin();
    typename MultiArray< 3,FFTWComplex<REAL> >::iterator y=inputY.begin();
    typename MultiArray< 3,FFTWComplex<REAL> >::iterator x=inputX.begin();

    // double one_over_sqrt_2 = 1.0 / M_SQRT2;
    
    // inputZ.bindElementChannel(0) = A.bindElementChannel(1) * -one_over_sqrt_2;
    // inputZ.bindElementChannel(1) = A.bindElementChannel(2) * -one_over_sqrt_2;
    // inputY.bindElementChannel(0) = A.bindElementChannel(0);
    // inputX.bindElementChannel(0) = A.bindElementChannel(1) * one_over_sqrt_2;
    // inputX.bindElementChannel(1) = A.bindElementChannel(2) * -one_over_sqrt_2;
    
    for (typename MultiArray<3,TinyVector<REAL,3> >::const_iterator p=A.begin();p!=A.end();++p,++z,++y,++x)
    {
        z->real() = -(*p)[1];
        z->imag() = -(*p)[2];
        *z *= 1/sqrt(2);
        y->real() =  (*p)[0];
        y->imag() =  0;
        x->real() = (*p)[1];
        x->imag() = -(*p)[2];
        *x *= 1/sqrt(2);
    }
    VH_A.resize(band + 1);
    REAL norm = M_PI/2*std::pow((REAL)radius,2.0);

    // SH coeffs
    int xa=0,xe=0,ya=0,ye=0,za=0,ze=0;

    VH_A[0].resize(3);
    for(int k=-1;k<=1;++k)
    {
        VH_A[0][1+k].resize(3);
        for (int m=-1;m<=1;++m)
        {
            MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > > vh;
            sphereVecHarmonic(vh, radius, gauss, 0, k, m);

            xa = (int) floor(pos[2] - vh.shape()[2] / 2);
            xe = xa + vh.shape()[2] - 1;
            ya = (int) floor(pos[1] - vh.shape()[1] / 2);
            ye = ya + vh.shape()[1] - 1;
            za = (int) floor(pos[0] - vh.shape()[0] / 2);
            ze = za + vh.shape()[0] - 1;
            
            // MultiArrayView<4, FFTWComplex<REAL>, StridedArrayTag > zsub = 
                // inputZ.subarray(Shape3(za,ya,xa),Shape3(ze,ye,xe)).insertSingletonDimension(0);
            
            // VH_A[0][1+k][1+m] = sum(zsub * vh.expandElements(0) / norm, FFTWComplex<double>());
            
            int zs=0;
            int ys=0;
            int xs=0;
            for (int z=za;z!=ze;z++,zs++)
            {
                for (int y=ya;y<=ye;y++,ys++)
                {
                    for (int x=xa;x<=xe;x++,xs++)
                    {
                        VH_A[0][1+k][1+m] += inputZ(z,y,x) * vh(zs,ys,xs)[0]/norm + inputZ(z,y,x) * vh(zs,ys,xs)[1]/norm + inputZ(z,y,x) * vh(zs,ys,xs)[2]/norm;
                    }
                }
            }
        }
    }
    for (int l = 1; l <= band; l++)
    {
        VH_A[l].resize(3);
        for (int k=-1;k<=1;++k)
        {
            VH_A[l][k+1].resize(2 * (l+k) + 1);
            for (int m = -(l+k); m <= (l+k); m++)
            {
                MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > > vh;
                sphereVecHarmonic(vh, radius, gauss, 0, k, m);

                xa = (int) floor(pos[2] - vh.shape()[2] / 2);
                xe = xa + vh.shape()[2] - 1;
                ya = (int) floor(pos[1] - vh.shape()[1] / 2);
                ye = ya + vh.shape()[1] - 1;
                za = (int) floor(pos[0] - vh.shape()[0] / 2);
                ze = za + vh.shape()[0] - 1;
                
                // MultiArrayView<4, FFTWComplex<REAL>, StridedArrayTag > zsub = 
                    // inputZ.subarray(Shape3(za,ya,xa),Shape3(ze,ye,xe)).insertSingletonDimension(0);
                
                // VH_A[l][1+k][1+m] = sum(zsub * vh.expandElements(0) / norm, FFTWComplex<double>());

                int zs=0;
                for (int z=za;z!=ze;z++,zs++)
                {
                    int ys=0;
                    for (int y=ya;y<=ye;y++,ys++)
                    {
                        int xs=0;
                        for (int x=xa;x<=xe;x++,xs++)
                        {
                            VH_A[l][1+k][1+m] += inputZ(z,y,x) * vh(zs,ys,xs)[0]/norm + inputZ(z,y,x) * vh(zs,ys,xs)[1]/norm + inputZ(z,y,x) * vh(zs,ys,xs)[2]/norm;
                        }
                    }
                }
            }
        }
    }
}

template <typename REAL>
void 
VHcenter(std::vector<std::vector<std::vector< FFTWComplex<REAL> > > >&VH_A, 
         unsigned int band, REAL gauss, REAL radius, 
         const MultiArray<3, TinyVector<REAL,3> > &A)
{
    VHpos(VH_A, band, gauss, radius, A, detail::centerOfBB(A));
}

template <typename REAL>
void 
VPHpos(std::vector<std::vector<std::vector<std::vector<FFTWComplex<REAL> > > > >&VH_A, 
       unsigned int band, REAL radius, 
       const MultiArray<3, TinyVector<REAL,3> > &A, const TinyVector<REAL, 3> pos)
{
    using namespace multi_math;
    
    vigra_fail("VPHpos is untested");
    
    REAL gauss = 1;
    //transform input to C^(2j+1)
    MultiArray< 3,FFTWComplex<REAL> >  inputZ(A.shape());
    MultiArray< 3,FFTWComplex<REAL> >  inputY(A.shape());
    MultiArray< 3,FFTWComplex<REAL> >  inputX(A.shape());
    
    // double one_over_sqrt_2 = 1.0 / M_SQRT2;
    
    // inputZ.bindElementChannel(0) = A.bindElementChannel(1) * -one_over_sqrt_2;
    // inputZ.bindElementChannel(1) = A.bindElementChannel(2) * -one_over_sqrt_2;
    // inputY.bindElementChannel(0) = A.bindElementChannel(0);
    // inputX.bindElementChannel(0) = A.bindElementChannel(1) * one_over_sqrt_2;
    // inputX.bindElementChannel(1) = A.bindElementChannel(2) * -one_over_sqrt_2;
    
    typename MultiArray< 3,FFTWComplex<REAL> >::iterator z=inputZ.begin();
    typename MultiArray< 3,FFTWComplex<REAL> >::iterator y=inputY.begin();
    typename MultiArray< 3,FFTWComplex<REAL> >::iterator x=inputX.begin();

    typename MultiArray<3,TinyVector<REAL,3> >::const_iterator p=A.begin();
    for (;p!=A.end();++p,++z,++y,++x)
    {
        z->real() = -(*p)[1];
        z->imag() = -(*p)[2];
        *z *= 1/sqrt(2);
        y->real() =  (*p)[0];
        y->imag() =  0;
        x->real() = (*p)[1];
        x->imag() = -(*p)[2];
        *x *= 1/sqrt(2);
    }

    VH_A.resize(band + 1);
    for(int n=1;n<=band;n++)
    {
        VH_A[n].resize(band + 1);
        REAL norm = M_PI/2*std::pow((REAL)radius,2.0);

        // SH coeffs
        int xa=0,xe=0,ya=0,ye=0,za=0,ze=0;

        VH_A[n][0].resize(3);
        for(int k=-1;k<=1;++k)
        {
            VH_A[n][0][1+k].resize(3);
            for (int m=-1;m<=1;++m)
            {
                MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > > vh;
                sphereVecHarmonic(vh, radius, gauss, 0, k, m);

                xa = (int) floor(pos[2] - vh.shape()[2] / 2);
                xe = xa + vh.shape()[2] - 1;
                ya = (int) floor(pos[1] - vh.shape()[1] / 2);
                ye = ya + vh.shape()[1] - 1;
                za = (int) floor(pos[0] - vh.shape()[0] / 2);
                ze = za + vh.shape()[0] - 1;

                // MultiArrayView<4, FFTWComplex<REAL>, StridedArrayTag > zsub = 
                    // inputZ.subarray(Shape3(za,ya,xa),Shape3(ze,ye,xe)).insertSingletonDimension(0);

                // VH_A[n][0][1+k][1+m] = sum(zsub * vh.expandElements(0) / norm, 
                                            // FFTWComplex<double>());
                
                int zs=0;
                for (int z=za;z!=ze;z++,zs++)
                {
                    int ys=0;
                    for (int y=ya;y<=ye;y++,ys++)
                    {
                        int xs=0;
                        for (int x=xa;x<=xe;x++,xs++)
                        {
                            VH_A[n][0][1+k][1+m] += inputZ(z,y,x) * vh(zs,ys,xs)[0]/norm + inputZ(z,y,x) * vh(zs,ys,xs)[1]/norm + inputZ(z,y,x) * vh(zs,ys,xs)[2]/norm;
                        }
                    }
                }
            }
        }
        for (int l = 1; l <= band; l++)
        {
            VH_A[n][l].resize(3);
            for (int k=-1;k<=1;++k)
            {
                VH_A[n][l][k+1].resize(2 * (l+k) + 1);
                for (int m = -(l+k); m <= (l+k); m++)
                {
                    MultiArray<3, TinyVector<FFTWComplex<REAL>,3 > > vh;
                    sphereVecHarmonic(vh, radius, gauss, 0, k, m);

                    xa = (int) floor(pos[2] - vh.shape()[2] / 2);
                    xe = xa + vh.shape()[2] - 1;
                    ya = (int) floor(pos[1] - vh.shape()[1] / 2);
                    ye = ya + vh.shape()[1] - 1;
                    za = (int) floor(pos[0] - vh.shape()[0] / 2);
                    ze = za + vh.shape()[0] - 1;

                    // MultiArrayView<4, FFTWComplex<REAL>, StridedArrayTag > zsub = 
                        // inputZ.subarray(Shape3(za,ya,xa),Shape3(ze,ye,xe)).insertSingletonDimension(0);
                    // VH_A[n][l][1+k][1+m] = sum(zsub * vh.expandElements(0) / norm, 
                                                // FFTWComplex<double>());

                    int zs=0;
                    int ys=0;
                    int xs=0;
                    for (int z=za;z!=ze;z++,zs++)
                    {
                        for (int y=ya;y<=ye;y++,ys++)
                        {
                            for (int x=xa;x<=xe;x++,xs++)
                            {
                                VH_A[n][l][1+k][1+m] += inputZ(z,y,x) * vh(zs,ys,xs)[0]/norm + inputZ(z,y,x) * vh(zs,ys,xs)[1]/norm + inputZ(z,y,x) * vh(zs,ys,xs)[2]/norm;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename REAL>
void 
VPHcenter(std::vector<std::vector<std::vector<std::vector< FFTWComplex<REAL> > > > >&VH_A, 
          unsigned int band, REAL radius, 
          const MultiArray<3, TinyVector<REAL,3> > &A)
{
    VPHpos(VH_A, band, radius, A, detail::centerOfBB(A));
}


} // namespace vigra 

#endif // VIGRA_INVARIANT_FEATURES3D_HXX
