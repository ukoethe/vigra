/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2014 by                                 */
/*               Ullrich Koethe,                                        */
/*               Esteban Pardo                                          */
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

#define VIGRA_CHECK_BOUNDS

#include <limits>
#include <algorithm>
#include <vigra/unittest.hxx>
#include <vigra/permutation.hxx>


namespace vigra
{

struct PermutationTest
{
    void testN1()
    {
        PlainChangesPermutations<1> permutations;
        shouldEqual(permutations.size(), 1);
        shouldEqual(permutations[0][0], 0);
        shouldEqual(permutations[0].sign(), 1);
    }

    void testN2()
    {
        PlainChangesPermutations<2> permutations;
        shouldEqual(permutations.size(), 2);
        shouldEqual(permutations[0][0], 0);
        shouldEqual(permutations[0][1], 1);
        shouldEqual(permutations[0].sign(), 1);
        shouldEqual(permutations[1][0], 1);
        shouldEqual(permutations[1][1], 0);
        shouldEqual(permutations[1].sign(), -1);
    }

    void testN3()
    {
        PlainChangesPermutations<3> permutations;
        shouldEqual(permutations.size(), 6);
        shouldEqual(permutations[0][0], 0);
        shouldEqual(permutations[0][1], 1);
        shouldEqual(permutations[0][2], 2);
        shouldEqual(permutations[0].sign(), 1);
        shouldEqual(permutations[1][0], 0);
        shouldEqual(permutations[1][1], 2);
        shouldEqual(permutations[1][2], 1);
        shouldEqual(permutations[1].sign(), -1);
        shouldEqual(permutations[2][0], 2);
        shouldEqual(permutations[2][1], 0);
        shouldEqual(permutations[2][2], 1);
        shouldEqual(permutations[2].sign(), 1);
        shouldEqual(permutations[3][0], 2);
        shouldEqual(permutations[3][1], 1);
        shouldEqual(permutations[3][2], 0);
        shouldEqual(permutations[3].sign(), -1);
        shouldEqual(permutations[4][0], 1);
        shouldEqual(permutations[4][1], 2);
        shouldEqual(permutations[4][2], 0);
        shouldEqual(permutations[4].sign(), 1);
        shouldEqual(permutations[5][0], 1);
        shouldEqual(permutations[5][1], 0);
        shouldEqual(permutations[5][2], 2);
        shouldEqual(permutations[5].sign(), -1);
    }
};

struct PermutationTestSuite : public vigra::test_suite
{
    PermutationTestSuite() : vigra::test_suite("PermutationTestSuite")
    {
        add(testCase(&PermutationTest::testN1));
        add(testCase(&PermutationTest::testN2));
        add(testCase(&PermutationTest::testN3));
    }
};

} // namespace vigra

int main(int argc, char** argv)
{
    vigra::PermutationTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cerr << test.report() << std::endl;

    return failed != 0;
}

