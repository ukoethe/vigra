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

using namespace vigra;

template <class Real>
struct ComplexRealFunctor
{
    Real operator()(std::complex<Real> const & c) const
    {
        return c.real();
    }
};

template <class Real>
struct ComplexImagFunctor
{
    Real operator()(std::complex<Real> const & c) const
    {
        return c.imag();
    }
};

struct InvariantFeaturesTest
{

    InvariantFeaturesTest()
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

            transformMultiArray(srcMultiArrayRange(wigner.get_D(l)),
                                destMultiArray(diff), ComplexImagFunctor<float>());
            inspectMultiArray(srcMultiArrayRange(diff), minmax);
            shouldEqual(minmax.min, 0.0f);
            shouldEqual(minmax.max, 0.0f);

            transformMultiArray(srcMultiArrayRange(wigner.get_D(l)),
                                destMultiArray(diff), ComplexRealFunctor<float>());
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
            
            // FIXME: transpose() shouldn't be necessary below
            transformMultiArray(srcMultiArrayRange(transpose(wigner2.get_D(l))),
                                destMultiArray(diff), ComplexRealFunctor<float>());
            diff -= ref[l];
            inspectMultiArray(srcMultiArrayRange(diff), minmax);
            should(minmax.min > -1e-4f);
            should(minmax.max <  1e-4f);
        }
        
        // FIXME: compute_D() with arbitrary angles, rot(), and rotatePH() are untested!
    }

};

struct FeaturesTestSuite
: public vigra::test_suite
{
    FeaturesTestSuite()
    : vigra::test_suite("FeaturesTestSuite")
    {
        add( testCase( &InvariantFeaturesTest::wignerMatrixTest));
    }
};

int main(int argc, char ** argv)
{
    FeaturesTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

