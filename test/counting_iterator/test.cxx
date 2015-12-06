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
#include "vigra/counting_iterator.hxx"

using namespace vigra;





struct CountingIteratorTest{   


    void test1(){
        
        auto iterBegin = vigra::CountingIterator<int>(0);
        auto iterEnd   = vigra::CountingIterator<int>(10);
        shouldEqual( std::distance(iterBegin,iterEnd), 10);
        shouldEqual( *iterEnd, 10);
        for(int i=0; i<10; ++i){
            shouldEqual( iterBegin[i], i);
        }
        for(int i=0; i<10; ++i){

            shouldEqual( std::distance(iterBegin,iterEnd), 10-i);
            should(iterBegin!=iterEnd);
            shouldEqual( *iterBegin, i);
            ++iterBegin;
        }
        should(iterBegin==iterEnd);
        --iterEnd;
        shouldEqual( *iterEnd, 9);
        --iterEnd;
        shouldEqual( *iterEnd, 8);

    }

};




 
struct CountingIteratorTestSuite
: public vigra::test_suite
{
    CountingIteratorTestSuite()
    : vigra::test_suite("CountingIteratorTestSuite")
    {   
        add( testCase( &CountingIteratorTest::test1));
    }
};

int main(int argc, char ** argv)
{
    CountingIteratorTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

