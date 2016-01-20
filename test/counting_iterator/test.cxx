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
#include <numeric>
#include "vigra/unittest.hxx"
#include "vigra/counting_iterator.hxx"

using namespace vigra;

struct CountingIteratorTest
{
    void testCountingIterator()
    {
        {
            auto iter1 = range(0,7,3);
            auto iter2 = range(0,8,3);
            auto iter3 = range(0,9,3);
            auto iterEnd   = iter1.end();
            shouldEqual( *iterEnd, 9);
            should( iterEnd == iter2.end());
            should( iterEnd == iter3.end());
            shouldEqual( std::distance(iter1,iterEnd), 3);
            shouldEqual( std::distance(iter2,iterEnd), 3);
            shouldEqual( std::distance(iter3,iterEnd), 3);

            for(int i=0; i<3; ++i){
                shouldEqual( iter1[i], 3*i);
                shouldEqual( iter2[i], 3*i);
                shouldEqual( iter3[i], 3*i);
                shouldEqual( iterEnd[-i], 9-3*i);
            }

            for(int i=0; i<6; ++i){

                shouldEqual( std::distance(iter1,iterEnd), 3-i);
                shouldEqual( std::distance(iter2,iterEnd), 3-i);
                shouldEqual( std::distance(iter3,iterEnd), 3-i);

                should(iter1 == iter2);
                should(iter1 == iter3);
                shouldNot(iter1 != iter2);
                shouldNot(iter1 != iter3);

                if(i != 3)
                {
                    should(iter1 != iterEnd);
                    shouldNot(iter1 == iterEnd);
                }
                else
                {
                    shouldNot(iter1 != iterEnd);
                    should(iter1 == iterEnd);
                }
                if(i < 3)
                {
                    should(iter1 < iterEnd);
                    shouldNot(iter1 >= iterEnd);
                    shouldNot(iter1.empty());
                }
                else
                {
                    shouldNot(iter1 < iterEnd);
                    should(iter1 >= iterEnd);
                    should(iter1.empty());
                }
                if(i < 4)
                {
                    should(iter1 <= iterEnd);
                    shouldNot(iter1 > iterEnd);
                }
                else
                {
                    shouldNot(iter1 <= iterEnd);
                    should(iter1 > iterEnd);
                }

                shouldEqual( *iter1, 3*i);
                shouldEqual( *iter2, 3*i);
                shouldEqual( *iter3, 3*i);

                ++iter1;
                iter2++;
                iter3+=1;
            }

            iter1 = iterEnd;
            --iter1;
            shouldEqual( *iter1, 6);
            iter1--;
            shouldEqual( *iter1, 3);
            iter1-=1;
            shouldEqual( *iter1, 0);
        }

        {
            auto iter1 = range(10, 1, -2);
            auto iter2 = range(10, 0, -2);
            auto iterEnd   = iter1.end();
            shouldEqual( *iterEnd, 0);
            should( iterEnd == iter2.end());
            shouldEqual( std::distance(iter1,iterEnd), 5);
            shouldEqual( std::distance(iter2,iterEnd), 5);

            for(int i=0; i<5; ++i){
                shouldEqual( iter1[i], 10-2*i);
                shouldEqual( iter2[i], 10-2*i);
                shouldEqual( iterEnd[-i], 2*i);
            }

            for(int i=0; i<5; ++i){

                shouldEqual( std::distance(iter1,iterEnd), 5-i);
                shouldEqual( std::distance(iter2,iterEnd), 5-i);

                should(iter1!=iterEnd);
                should(iter2!=iterEnd);

                shouldEqual( *iter1, 10-2*i);
                shouldEqual( *iter2, 10-2*i);

                ++iter1;
                iter2++;
            }

            should(iter1.empty());
            should(iter1==iterEnd);
            should(iter2.empty());
            should(iter2==iterEnd);
            --iter1;
            shouldEqual( *iter1, 2);
            iter1--;
            shouldEqual( *iter1, 4);
        }

        int count = 0;
        for(auto i: range(0, 19, 2))
        {
            shouldEqual( i, count);
            count += 2;
        }
        shouldEqual(20, count);
        for(auto i: range(20, 0, -2))
        {
            shouldEqual( i, count);
            count -= 2;
        }
        shouldEqual(0, count);

        {
            auto iter = range(5),
                 end  = iter.end();
            shouldEqual(std::accumulate(iter, end, 0), 10);
        }

        {
            double c = 1.0;
            for(auto i: range(1.0, 1.6, 0.1)) // 1.6 is excluded
            {
                shouldEqualTolerance(c, i, 1e-15);
                c += 0.1;
            }
            shouldEqualTolerance(c, 1.6, 1e-15);

            c = 1.0;
            for(auto i: range(1.0, 1.61, 0.1)) // 1.6 is included
            {
                shouldEqualTolerance(c, i, 1e-15);
                c += 0.1;
            }
            shouldEqualTolerance(c, 1.7, 1e-15);

            auto iter = range(1.0, 1.6, 0.1),
                 end  = iter.end();
            shouldEqual(end-iter, 6);
            c = 1.0;
            for(int i=0; i < 9; ++i, ++iter)
            {
                if(i != 6)
                {
                    should(iter != end);
                    shouldNot(iter == end);
                }
                else
                {
                    should(iter == end);
                    shouldNot(iter != end);
                }
                if(i < 6)
                {
                    should(iter < end);
                    shouldNot(iter >= end);
                    shouldNot(iter.empty());
                }
                else
                {
                    shouldNot(iter < end);
                    should(iter >= end);
                    should(iter.empty());
                }
                if(i < 7)
                {
                    should(iter <= end);
                    shouldNot(iter > end);
                }
                else
                {
                    shouldNot(iter <= end);
                    should(iter > end);
                }
                shouldEqualTolerance(*iter, c, 1e-15);
                c += 0.1;
            }
            shouldEqualTolerance(c, 1.9, 1e-15);

            c = 1.6;
            iter = range(1.6, 1.0, -0.1);
            end  = iter.end();
            for(; iter <= end; ++iter)
            {
                shouldEqual(*iter, c);
                c -= 0.1;
            }
            shouldEqualTolerance(c, 0.9, 1e-15);
        }
    }
};
 
struct CountingIteratorTestSuite
: public vigra::test_suite
{
    CountingIteratorTestSuite()
    : vigra::test_suite("CountingIteratorTestSuite")
    {   
        add( testCase( &CountingIteratorTest::testCountingIterator));
    }
};

int main(int argc, char ** argv)
{
    CountingIteratorTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

