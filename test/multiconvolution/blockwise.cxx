/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
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
#include "vigra/unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_resize.hxx"

#include "vigra/multi_blocking.hxx"
#include "vigra/multi_blockwise.hxx"

using namespace vigra;


template <class Iterator>
void myFillRandom(Iterator begin, Iterator end, int maximum)
{
    using namespace std;

    for( ; begin != end; ++begin)
        *begin = rand() % maximum;
}

struct BlockwiseConvolutionTest{




    BlockwiseConvolutionTest(){       

    }


    void test1()
    {
        double sigma = 1.0;
        blockwise::BlockwiseConvolutionOptions<2>   opt;
        TinyVector<double, 2> sigmaV(sigma, sigma);

        opt.setStdDev(sigmaV);
        opt.setBlockShape(TinyVector<int, 2>(5,7));


        typedef MultiArray<2, double> Array;
        typedef Array::difference_type Shape;
        
        // random array
        Shape shape(200,200);
        Array data(shape);
        myFillRandom(data.begin(), data.end(), 2000);

        // blockwise
        Array resB(shape);
        blockwise::gaussianSmoothMultiArray(data, resB, opt);

        Array res(shape);
        gaussianSmoothMultiArray(data, res, sigma);

        shouldEqualSequenceTolerance(
            res.begin(), 
            res.end(), 
            resB.begin(), 
            1e-14
        );

    }
};
 
struct BlockwiseConvolutionTestSuite
: public vigra::test_suite
{
    BlockwiseConvolutionTestSuite()
    : vigra::test_suite("BlockwiseConvolutionTestSuite")
    {   
        add( testCase( &BlockwiseConvolutionTest::test1));
    }
};

int main(int argc, char ** argv)
{
    BlockwiseConvolutionTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

