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

#include <iostream>
#include <functional>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include "unittest.hxx"
#include <vigra/wigner-matrix.hxx>
#include <vigra/multi_pointoperators.hxx>
#include "wigner-matrix-reference.hxx"
#include <vigra/harmonics.hxx>
#include <vigra/invariant_features3D.hxx>
#include <vigra/hdf5impex.hxx>

using namespace vigra;

template <class Real>
struct ComplexRealFunctor
{
    Real operator()(FFTWComplex<Real> const & c) const
    {
        return c.real();
    }
};

template <class Real>
struct ComplexImagFunctor
{
    Real operator()(FFTWComplex<Real> const & c) const
    {
        return c.imag();
    }
};

struct InvariantFeaturesTest
{

    void SHbaseTest()
    {
        std::cerr<<"SHbaseTest\n";
        std::vector<std::vector<vigra::MultiArray<3,FFTWComplex<float> > > > SHbaseF;
        vigra::computeSHbaseF<float>(10,2,5,SHbaseF);

	//FIXME: this is a qualitative test, find objective criteria 
        //HDF5File file("SHbaseTest.h5",HDF5File::New);
	for (int l=0;l<SHbaseF.size();l++)
	    for (int m=0;m<SHbaseF[l].size();m++)
	    {
		std::string name = "r_"+vigra::asString(l)+"_"+vigra::asString(-l+m);
		vigra::MultiArray<3, float > tmp;
		tmp.reshape(SHbaseF[l][m].shape(),0);
		for (int z=0;z<SHbaseF[l][m].shape()[0];++z)
		    for (int y=0;y<SHbaseF[l][m].shape()[1];++y)
			for (int x=0;x<SHbaseF[l][m].shape()[2];++x)
			    tmp(z,y,x)=real(SHbaseF[l][m](z,y,x));
		//file.write(name, tmp);
		name = "c_"+vigra::asString(l)+"_"+vigra::asString(-l+m);
		for (int z=0;z<SHbaseF[l][m].shape()[0];++z)
		    for (int y=0;y<SHbaseF[l][m].shape()[1];++y)
			for (int x=0;x<SHbaseF[l][m].shape()[2];++x)
			    tmp(z,y,x)=imag(SHbaseF[l][m](z,y,x));
		//file.write(name, tmp);

	    }
    }

    void SHcenterTest()
    {
        std::cerr<<"SHCenterTest\n";
        typedef MultiArray<3, float>::difference_type Shape;
        vigra::MultiArrayShape<3>::type testShape(51,51,51);
        vigra::MultiArray<3,float> testVolume;
        vigra::MultiArray<3,float> testVolumeRot;
        testVolume.reshape(testShape,0);
        testVolumeRot.reshape(testShape,0);
        vigra::MultiArrayView<3,float> testView = testVolume.subarray(Shape(0,0,0), Shape(25, 25, 25));
        vigra::MultiArrayView<3,float> testViewRot = testVolumeRot.subarray(Shape(26,0,0), Shape(51, 25, 25));
        testView+=1;
        testViewRot+=1;

	std::vector<std::vector<vigra::MultiArray<3,FFTWComplex<float> > > > SHbaseF;
        vigra::computeSHbaseF<float>(10,2,5,SHbaseF);

        std::vector<std::vector<FFTWComplex<float> >  >  SH_A;
        std::vector<std::vector<FFTWComplex<float> >  >  SH_A_rot;
        std::vector<std::vector<FFTWComplex<float> >  >  SH_A_rerot;

        vigra::SHcenter<float>(SH_A, 10, 5, testVolume, SHbaseF);

        //check if all coefficients are real for m=0
        float cm0=0;
        for (int l=0;l<SHbaseF.size();l++)
	    cm0+=imag(SH_A[l][l]);
        should(abs(cm0) < 1e-4f);

        vigra::MultiArray<3,float> reconstVolume;
        vigra::MultiArray<3,float> reconstVolumeRot;
        vigra::MultiArray<3,float> reconstVolumeReRot;
        reconstVolume.reshape(testShape,0);
        reconstVolumeRot.reshape(testShape,0);
        reconstVolumeReRot.reshape(testShape,0);

        vigra::WignerMatrix<float> W(5);

        //check if zero rotation changes coeffs
        W.rotateSH(SH_A, (float) 0.0, (float) 0.0, (float) 0.0, SH_A_rot);
        float diff = 0;
	for (int l=0;l<SHbaseF.size();l++)
	    for (int m=0;m<SHbaseF[l].size();m++)
		diff+=abs(SH_A[l][m]-SH_A_rot[l][m]);
        should(abs(diff) < 1e-4f);

        //check if rotation and unrotation  changes coeffs
        W.rotateSH(SH_A, (float) M_PI/4, (float) M_PI/8, (float) M_PI/2, SH_A_rot);
	W.rotateSH(SH_A_rot, (float) -M_PI/2, (float) -M_PI/8, (float) -M_PI/4, SH_A_rerot);
        diff = 0;
	for (int l=0;l<SHbaseF.size();l++)
	    for (int m=0;m<SHbaseF[l].size();m++)
	    {
		//std::cerr<<l<<" "<<m-l<<" "<<SH_A[l][m]<<" "<<SH_A_rerot[l][m]<<"\n ";
		diff+=abs(SH_A[l][m]-SH_A_rerot[l][m]);
	    }
        should(abs(diff) < 1e-2f);

	
        //check for coeff diff after rot in spacial domain and unrot in harmonics
        vigra::SHcenter<float>(SH_A_rot, 10, 5, testVolumeRot, SHbaseF);
        W.rotateSH(SH_A, (float) M_PI/2, (float) 0.0, (float) 0, SH_A_rerot);
        diff = 0;
	for (int l=0;l<SHbaseF.size();l++)
	    for (int m=0;m<SHbaseF[l].size();m++)
	    {
		//std::cerr<<l<<" "<<m-l<<" "<<SH_A_rot[l][m]<<" "<<SH_A_rerot[l][m]<<"\n ";
		diff+=abs(SH_A_rot[l][m]-SH_A_rerot[l][m]);
	    }
	std::cerr<<"diff:"<<abs(diff)<<"\n";
        //should(abs(diff) < 1);//FIXME

        vigra::reconstSH<float>(10,2,5,reconstVolume,SH_A,SHbaseF);
        vigra::reconstSH<float>(10,2,5,reconstVolumeRot,SH_A_rot,SHbaseF);
        vigra::reconstSH<float>(10,2,5,reconstVolumeReRot,SH_A_rerot,SHbaseF);

        //HDF5File file("SHreconstTest.h5",HDF5File::New);
        //file.write("reconst", reconstVolume);
        //file.write("reconstRot", reconstVolumeRot);
        //file.write("reconstReRot", reconstVolumeReRot);
        //file.write("orig",testVolume);
    }


    struct pow2f
    {
	float operator()(float t) const
	{
	    return pow( (double)t,2.0);
	}
    };

    void SHfeatureTest()
    {
        std::cerr<<"SHCenterTest\n";
        typedef MultiArray<3, float>::difference_type Shape;
        vigra::MultiArrayShape<3>::type testShape(51,51,51);
        vigra::MultiArray<3,float> testVolume;
        vigra::MultiArray<3,float> testVolumeRot;
        testVolume.reshape(testShape,0);
        testVolumeRot.reshape(testShape,0);
        vigra::MultiArrayView<3,float> testView = testVolume.subarray(Shape(0,0,0), Shape(25, 25, 25));
        vigra::MultiArrayView<3,float> testViewRot = testVolumeRot.subarray(Shape(26,0,0), Shape(51, 25, 25));
        testView+=1;
        testViewRot+=1;

        std::vector<std::vector<vigra::MultiArray<3,FFTWComplex<float> > > > SHbaseF;
        vigra::computeSHbaseF<float>(10,2,5,SHbaseF);

        std::vector<std::vector<FFTWComplex<float> >  >  SH_A;
	std::vector<std::vector<FFTWComplex<float> >  >  SH_A_rot;
        std::vector<std::vector<FFTWComplex<float> >  >  SH_A_real_rot;

        vigra::SHcenter<float>(SH_A, 10, 5, testVolume, SHbaseF);
	vigra::WignerMatrix<float> W(5);
	W.rotateSH(SH_A, (float) M_PI/4, (float) M_PI/8, (float) M_PI/2, SH_A_rot);
	vigra::SHcenter<float>(SH_A_real_rot, 10, 5, testVolumeRot, SHbaseF);

	//test if SHabs is the same 
	std::vector<float> res_A, res_A_rot, res_A_real_rot;
	vigra::SHabs(SH_A, res_A);
	vigra::SHabs(SH_A_rot, res_A_rot);
	vigra::SHabs(SH_A_real_rot, res_A_real_rot);
	float diff1=0,diff2=0;
	for (int l=0;l<res_A.size();l++)
	{
	    //std::cerr<<res_A[l]<<" "<<res_A_rot[l]<<" "<<res_A_real_rot[l]<<"\n";
	    diff1+=abs(res_A[l]-res_A_rot[l]);
	    diff2+=abs(res_A[l]-res_A_real_rot[l]);
	}
	should(abs(diff1) < 1e-1f);//should be zero, but numerical problems
	should(abs(diff2) < 1);//real rotations not as accurate (due to band limitation)


        //test if SHbispec is the same 
        std::vector<FFTWComplex<float> > resC_A, resC_A_rot, resC_A_real_rot;
	vigra::SHbispec(SH_A, resC_A);
        vigra::SHbispec(SH_A_rot, resC_A_rot);
        vigra::SHbispec(SH_A_real_rot, resC_A_real_rot);
        diff1=0;
	diff2=0;
        for (int l=0;l<resC_A.size();l++)
        {
            //std::cerr<<resC_A[l]<<" "<<resC_A_rot[l]<<" "<<resC_A_real_rot[l]<<"\n";
            diff1+=abs(resC_A[l]-resC_A_rot[l]);
            diff2+=abs(resC_A[l]-resC_A_real_rot[l]);
        }
        //std::cerr<<abs(diff1)/resC_A.size()<<" "<<abs(diff2)/resC_A.size()<<"\n";
	should(abs(diff1)/resC_A.size() < 1);//should be zero, but numerical problems
        should(abs(diff2)/resC_A.size() < 100000);//feature very sensitive: large values for real data 

	//test SHcorr
        //check if normalzed cross correlation has max=1 at zero rotation offset
        float cmax,phi,theta,psi=0;
        vigra::MultiArray<3,float > corMat;
        SHxcorr(SH_A,SH_A,64,true,cmax,phi,theta,psi, corMat);
        should( (abs(phi+psi)<0.1 || abs(phi+psi-2*M_PI) < 0.1) && (theta<0.1) && (abs(cmax-1.0)< 1e-4f) ); //note: at theta=0 abitrary combinations of phi+psi%2PI are possible
        //std::cerr<<"xcorr: "<<cmax<<" "<<phi<<" "<<theta<<" "<<psi<<" : "<<phi+psi<<"\n";

	//check if normalzed cross correlation predicts correct correlation offset
	W.rotateSH(SH_A, (float) M_PI/2, (float) M_PI/2, (float) M_PI/2, SH_A_rot);
	SHxcorr(SH_A,SH_A_rot,64,true,cmax,phi,theta,psi, corMat);
	//std::cerr<<"xcorr: "<<cmax<<" "<<phi<<" "<<theta<<" "<<psi<<"\n";
	should( (abs(phi-M_PI/2)<1e-1f) && (abs(theta-M_PI/2)<1e-1f) && (abs(psi-M_PI/2)<1e-1f) && (abs(cmax-1)<1e-1f) ); 

	//HDF5File corrfile("SHCorTest.h5",HDF5File::New);
        //corrfile.write("cormat", corMat);
	
	//test invariance of SHautocorr
	float res1,res2;
	pow2f myfunctor;
	SHautocorr(SH_A, 32, myfunctor, res1);
	SHautocorr(SH_A_rot, 32, myfunctor, res2);
	should( abs(res1-res2)<1e-3f);
    }

    void SHallTest()
    {
        std::cerr<<"SHallTest\n";
        typedef MultiArray<3, float>::difference_type Shape;
        vigra::MultiArrayShape<3>::type testShape(51,51,51);
        vigra::MultiArray<3,float> testVolume;
        testVolume.reshape(testShape,0);
        vigra::MultiArrayView<3,float> testView = testVolume.subarray(Shape(0,0,0), Shape(25, 25, 25));
        testView+=1;

        std::vector<std::vector<vigra::MultiArray<3,FFTWComplex<float> > > > SHbaseF;
        vigra::computeSHbaseF<float>(10,2,5,SHbaseF);

        std::vector<std::vector<FFTWComplex<float> >  >  SH_A;
        vigra::SHcenter<float>(SH_A, 10, 5, testVolume, SHbaseF);

	std::vector<std::vector<vigra::MultiArray<3,FFTWComplex<float> >  > > SH_All;
	vigra::Array2SH<float>(SH_All,5, 10, SHbaseF, testVolume);
	
	//test if computation in fourier space equals computation in spacial domain
	float diff = 0;
        for (int l=0;l<SHbaseF.size();l++)
            for (int m=0;m<SHbaseF[l].size();m++)
            {
                //std::cerr<<l<<" "<<m-l<<" "<<SH_A[l][m]<<" "<<SH_All[l][m](25,25,25)<<"\n ";
                diff+=abs(SH_A[l][m]-SH_All[l][m](25,25,25));
            }
        should(abs(diff) < 1e-2f);

	//test SHabs all
	std::vector<vigra::MultiArray<3,float> > res_All;
	SHabs(SH_All,res_All);
	std::vector<float> res_A;
	vigra::SHabs(SH_A, res_A);
	diff=0;
	for (int i=0;i<SH_A.size();i++)
	    diff+=abs(res_A[i]-res_All[i](25,25,25));
	should(abs(diff) < 1);

	//test SHbispec all;
	std::vector<vigra::MultiArray<3,FFTWComplex<float> > > res_All2;
	SHbispec(SH_All,res_All2);
	std::vector<FFTWComplex<float> > res_A2;
	SHbispec(SH_A,res_A2);
	diff=0;
	for (int i=0;i<SH_A.size();i++)
	    diff+=abs(res_A[i]-res_All[i](25,25,25));
	should(abs(diff) < 1);


    }

    void PHbaseTest()
    {
	std::cerr<<"PHbaseTest\n";
	std::vector<std::vector<std::vector<vigra::MultiArray<3,FFTWComplex<float> > > > > PHbaseF;
	vigra::computePHbaseF<float>(10,5,PHbaseF,false);
	
	//FIXME: this is a qualitative test, find objective criteria
	//HDF5File file("PHbaseTest.h5",HDF5File::New);
	for (int n=0;n<PHbaseF.size();n++)
	    for (int l=0;l<PHbaseF[n].size();l++)
		for (int m=0;m<PHbaseF[n][l].size();m++)
		{
		    std::string name = "r_"+vigra::asString(n)+"_"+vigra::asString(l)+"_"+vigra::asString(-l+m);
		    //std::cerr<<name<<" "<<PHbaseF[n][l][m].shape()<<" "<<k<<"\n";
		    vigra::MultiArray<3, float > tmp;
		    tmp.reshape(PHbaseF[n][l][m].shape(),0);
		    for (int z=0;z<PHbaseF[n][l][m].shape()[0];++z)
			for (int y=0;y<PHbaseF[n][l][m].shape()[1];++y)
			    for (int x=0;x<PHbaseF[n][l][m].shape()[2];++x)
				tmp(z,y,x)=real(PHbaseF[n][l][m](z,y,x));
		    //file.write(name, tmp);
		    name = "c_"+vigra::asString(n)+"_"+vigra::asString(l)+"_"+vigra::asString(-l+m);
		    for (int z=0;z<PHbaseF[n][l][m].shape()[0];++z)
		        for (int y=0;y<PHbaseF[n][l][m].shape()[1];++y)
		            for (int x=0;x<PHbaseF[n][l][m].shape()[2];++x)
		                tmp(z,y,x)=imag(PHbaseF[n][l][m](z,y,x));
		    //file.write(name, tmp);

		}
    }

    void PHcenterTest()
    {
	std::cerr<<"PHCenterTest\n";
	typedef MultiArray<3, float>::difference_type Shape;
	vigra::MultiArrayShape<3>::type testShape(51,51,51);
	vigra::MultiArray<3,float> testVolume;
	vigra::MultiArray<3,float> testVolumeRot;
	testVolume.reshape(testShape,0);
	testVolumeRot.reshape(testShape,0);
	vigra::MultiArrayView<3,float> testView = testVolume.subarray(Shape(0,0,0), Shape(25, 25, 25));
	vigra::MultiArrayView<3,float> testViewRot = testVolumeRot.subarray(Shape(26,0,0), Shape(51, 25, 25));
	testView+=1;
	testViewRot+=1;

	std::vector<std::vector<std::vector<vigra::MultiArray<3,FFTWComplex<float> > > > > PHbaseF;
	vigra::computePHbaseF<float>(10,5,PHbaseF);

	std::vector<std::vector<std::vector<FFTWComplex<float> >  > > PH_A;
	std::vector<std::vector<std::vector<FFTWComplex<float> >  > > PH_A_rot;
        std::vector<std::vector<std::vector<FFTWComplex<float> >  > > PH_A_rerot;

	vigra::PHcenter<float>(PH_A, 10, 5, testVolume, PHbaseF);
	
	//check if all coefficients are real for m=0
	float cm0=0;
	for (int n=1;n<PHbaseF.size();n++)
            for (int l=0;l<PHbaseF[n].size();l++)
		    cm0+=imag(PH_A[n][l][l]);
	should(abs(cm0) < 1e-4f);

	vigra::MultiArray<3,float> reconstVolume;
	vigra::MultiArray<3,float> reconstVolumeRot;
	vigra::MultiArray<3,float> reconstVolumeReRot;
	reconstVolume.reshape(testShape,0);
	reconstVolumeRot.reshape(testShape,0);
	reconstVolumeReRot.reshape(testShape,0);
	
	vigra::WignerMatrix<float> W(5);
	
	//check if zero rotation changes coeffs
	W.rotatePH(PH_A, (float) 0.0, (float) 0.0, (float) 0.0, PH_A_rot);
        float diff = 0;
        for (int n=1;n<PHbaseF.size();n++)
            for (int l=0;l<PHbaseF[n].size();l++)
                for (int m=0;m<PHbaseF[n][l].size();m++)
                    diff+=abs(PH_A[n][l][m]-PH_A_rot[n][l][m]);
	should(abs(diff) < 1e-4f);

        //check if rotation keeps all coefficients real for l=0, m=0
        W.rotatePH(PH_A, (float) M_PI/4, (float) M_PI/8, (float) M_PI/2, PH_A_rot);
        cm0=0;
        for (int n=1;n<PHbaseF.size();n++)
	{
	    //std::cerr<<n<<" "<<PH_A_rot[n][0][0]<<"\n ";
            cm0+=imag(PH_A_rot[n][0][0]);
	}
        should(abs(cm0) < 1e-4f);


	//check if rotation and unrotation  changes coeffs
	W.rotatePH(PH_A_rot, (float) -M_PI/2, (float) -M_PI/8, (float) -M_PI/4, PH_A_rerot);
        diff = 0;
        for (int n=1;n<PHbaseF.size();n++)
            for (int l=0;l<PHbaseF[n].size();l++)
                for (int m=0;m<PHbaseF[n][l].size();m++)
                {
                    //std::cerr<<n<<" "<<l<<" "<<m-l<<" "<<PH_A_rot[n][l][m]<<" "<<PH_A_rerot[n][l][m]<<"\n ";
                    diff+=abs(PH_A[n][l][m]-PH_A_rerot[n][l][m]);
                }
        should(abs(diff) < 1e-4f);

	//check for coeff diff after rot in spacial domain and unrot un harmonics
	vigra::PHcenter<float>(PH_A_rot, 10, 5, testVolumeRot, PHbaseF);
	W.rotatePH(PH_A_rot, (float) 0.0, (float) M_PI/2, (float) 0.0, PH_A_rerot);
	diff = 0;
        for (int n=1;n<PHbaseF.size();n++)
            for (int l=0;l<PHbaseF[n].size();l++)
                for (int m=0;m<PHbaseF[n][l].size();m++)
                {
                    //std::cerr<<n<<" "<<l<<" "<<m-l<<" "<<PH_A[n][l][m]<<" "<<PH_A_rerot[n][l][m]<<"\n ";
                    diff+=abs(PH_A[n][l][m]-PH_A_rerot[n][l][m]);
                }
        //should(abs(diff) < 1);
	//FIXME

	vigra::reconstPH<float>(10,5,reconstVolume,PH_A,PHbaseF);
	vigra::reconstPH<float>(10,5,reconstVolumeRot,PH_A_rot,PHbaseF);
	vigra::reconstPH<float>(10,5,reconstVolumeReRot,PH_A_rerot,PHbaseF);
	
	//HDF5File file("PHreconstTest.h5",HDF5File::New);
	//file.write("reconst", reconstVolume);
	//file.write("reconstRot", reconstVolumeRot);
	//file.write("reconstReRot", reconstVolumeReRot);
	//file.write("orig",testVolume);
	
	
    }


    void PHallTest()
    {

    }

    void wignerMatrixTest()
    {
        typedef Matrix<float> M;
        typedef MultiArrayShape<2>::type Shape;
        
        int l_max = 15;
        WignerMatrix<float> wigner(l_max);
        
        M ref[] = { M(), 
                     M(3, 3, wignerRef1),
                     M(5, 5, wignerRef2),
                     M(7, 7, wignerRef3),
                     M(9, 9, wignerRef4),
                     M(11, 11, wignerRef5),
                     M(13, 13, wignerRef6),
                     M(15, 15, wignerRef7),
                     M(17, 17, wignerRef8),
                     M(19, 19, wignerRef9),
                     M(21, 21, wignerRef10),
                     M(23, 23, wignerRef11),
                     M(25, 25, wignerRef12),
                     M(27, 27, wignerRef13),
                     M(29, 29, wignerRef14),
                     M(31, 31, wignerRef15) };
        
        for(int l=1; l<=l_max; ++l)
        {
            wigner.compute_D(l);
            
            shouldEqual(wigner.get_D(l).shape(), Shape(2*l+1, 2*l+1));
            
            M diff(2*l+1, 2*l+1);
            FindMinMax<float> minmax;

	    transformMultiArray(srcMultiArrayRange(wigner.get_D(l)), destMultiArray(diff), ComplexImagFunctor<float>());
            inspectMultiArray(srcMultiArrayRange(diff), minmax);
            shouldEqual(minmax.min, 0.0f);
            shouldEqual(minmax.max, 0.0f);

	    transformMultiArray(srcMultiArrayRange(wigner.get_D(l)), destMultiArray(diff), ComplexRealFunctor<float>());
	    diff -= ref[l];
            inspectMultiArray(srcMultiArrayRange(diff), minmax);
            should(minmax.min > -1e-4f);
            should(minmax.max <  1e-4f);
        }
        
        WignerMatrix<float> wigner2(l_max);
        for(int l=1; l<=l_max; ++l)
        {
            wigner2.compute_D(l, 0.0f, float(M_PI / 2.0), 0.0f);
            
            shouldEqual(wigner2.get_D(l).shape(), Shape(2*l+1, 2*l+1));
            
            M diff(2*l+1, 2*l+1);
            FindMinMax<float> minmax;

			transformMultiArray(srcMultiArrayRange(wigner.get_D(l)),
				                destMultiArray(diff), ComplexImagFunctor<float>());
            inspectMultiArray(srcMultiArrayRange(diff), minmax);
            shouldEqual(minmax.min, 0.0f);
            shouldEqual(minmax.max, 0.0f);
            
	    transformMultiArray(srcMultiArrayRange((wigner2.get_D(l))), destMultiArray(diff), ComplexRealFunctor<float>());
	    diff -= ref[l];
	    inspectMultiArray(srcMultiArrayRange(diff), minmax);
	    should(minmax.min > -1e-4f);
            should(minmax.max <  1e-4f);
        }
        
        // FIXME: compute_D() with arbitrary angles, rot()
    }

};

struct FeaturesTestSuite
: public vigra::test_suite
{
    FeaturesTestSuite()
    : vigra::test_suite("FeaturesTestSuite")
    {
        add( testCase( &InvariantFeaturesTest::wignerMatrixTest));
	add( testCase( &InvariantFeaturesTest::SHbaseTest));
	add( testCase( &InvariantFeaturesTest::SHcenterTest));
	add( testCase( &InvariantFeaturesTest::SHallTest));
	add( testCase( &InvariantFeaturesTest::SHfeatureTest));
	add( testCase( &InvariantFeaturesTest::PHbaseTest));
	add( testCase( &InvariantFeaturesTest::PHcenterTest));
    }
};

int main(int argc, char ** argv)
{
    FeaturesTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

