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
#include "vigra/delegate/delegate.hxx"

using namespace vigra;



int ff1(int a){
    return a;
}


int ff2(int a,int b){
    return a+b;
}




struct DelegateTest{   
    typedef  vigra::delegate1<int,int> D1;
    typedef  vigra::delegate2<int,int,int> D2;


    int mf1(int a){
        return a;
    }


    int mf2(int a,int b){
        return a+b;
    }

    int mf1(int a)const{
        return a+1;
    }


    int mf2(int a,int b)const{
        return a+b+1;
    }

    void test1(){
        {
        D1 d1(D1::from_method< DelegateTest,&DelegateTest::mf1 >(this));
        shouldEqual(d1(2),2);
        shouldEqual(d1(3),3);
        }
        {
        D1 d1(D1::from_const_method< DelegateTest,&DelegateTest::mf1 >(this));
        shouldEqual(d1(2),2+1);
        shouldEqual(d1(3),3+1);
        }
        {
        D1 d1(D1::from_function< &ff1 >());
        shouldEqual(d1(2),2);
        shouldEqual(d1(3),3);
        }
    }
    void test2(){
        {
        D2 d2(D2::from_method< DelegateTest,&DelegateTest::mf2 >(this));
        shouldEqual(d2(2,2),4);
        shouldEqual(d2(3,2),5);
        }
        {
        D2 d2(D2::from_const_method< DelegateTest,&DelegateTest::mf2 >(this));
        shouldEqual(d2(2,2),4+1);
        shouldEqual(d2(3,2),5+1);
        }
        {
        D2 d2(D2::from_function< &ff2 >());
        shouldEqual(d2(2,2),4);
        shouldEqual(d2(3,2),5);
        }
    }
};




 
struct DelegatesTestSuite
: public vigra::test_suite
{
    DelegatesTestSuite()
    : vigra::test_suite("DelegatesTestSuite")
    {   
        add( testCase( &DelegateTest::test1));
        add( testCase( &DelegateTest::test2));
    }
};

int main(int argc, char ** argv)
{
    DelegatesTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

