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
#include "vigra/sparse_array.hxx"

using namespace vigra;




struct SparseVectorTest{   

    typedef SparseMapVector<int> IntVec;
    typedef IntVec::value_type value_type;
    typedef IntVec::coordinate_value_pair coordinate_value_pair;

    
    void testProxy0(){
        IntVec vec(100);
        vec(1)=1;
        int a = vec.asConst()[1];
        shouldEqual(a,1);
        a = vec.asConst()[2];
        shouldEqual(a,0);

        vec(10)=10;
        a=vec(10);
        shouldEqual(a,10);

        vec(11)=11;
        shouldEqual(vec(11),11);
        ++vec(11);
        shouldEqual(vec(11),12);
        vec(11)++;
        shouldEqual(vec(11),13);
        --vec(11);
        shouldEqual(vec(11),12);
        vec(11)--;
        shouldEqual(vec(11),11);

        vec[90]++;
        shouldEqual(vec(90),1);

        ++vec[91];
        shouldEqual(vec(91),1);

        --vec[91];
        shouldEqual(vec(91),0);


        vec[50]+=10;
        shouldEqual(vec(50),10);
    }
    
    #if 0

    void testProxy1(){
        IntVec vec(100);

        //vec.write(20,25);
        //shouldEqual(vec.read(20),25);


        vec(10)=10;
        int a;
        a=vec(10);
        shouldEqual(a,10);

        vec(11)=11;
        shouldEqual(vec(11),11);

        vec(12)=11;
        ++vec(12);
        shouldEqual(vec(12),12);

        vec(13)=12;
        vec(13)++;
        shouldEqual(vec(13),13);

        vec(15)=10;
        vec(15)+=5;
        shouldEqual(vec(15),15);

        vec(15)=10;
        vec(5)=5;
        shouldEqual(vec(15),10);
        shouldEqual(vec(5),5);
        vec(15)+=vec(5);
        shouldEqual(vec(15),15);


        vec(15)-=vec(5);
        shouldEqual(vec(15),10);
        vec(15)++;
        vec(15)++;
        vec(15)++;
        vec(15)++;
        vec(15)++;
        shouldEqual(vec(15),15);


        vec(77)*=vec(77);
        shouldEqual(vec(77),0);

        vec(77)=2;
        vec(77)*=vec(77);
        shouldEqual(vec(77),4);
    }

    void testProxy2(){
        IntVec vec(100);

        shouldEqual(vec(10),0);
        vec(10)++;
        shouldEqual(vec(10),1);


        shouldEqual(vec(11),0);
        ++vec(11);
        shouldEqual(vec(11),1);


        shouldEqual(vec(12),0);
        vec(12)+=1;
        shouldEqual(vec(12),1);

        shouldEqual(vec(13),0);
        vec(13)-=1;
        shouldEqual(vec(13),-1);


        shouldEqual(vec(14),0);
        vec(14)*=2;
        shouldEqual(vec(14),0);

        shouldEqual(vec(14),0);
        vec(14)/=2;
        shouldEqual(vec(14),0);

        vec(14)=2;
        vec(14)*=2;
        shouldEqual(vec(14),4);
        vec(14)/=2;
        shouldEqual(vec(14),2);
    }


    void testProxy3(){
        IntVec vec(100);
        int res;


        res = vec(20)+vec(10);
        shouldEqual(res,0);

        vec(20)=2;
        res = vec(20)+vec(10);
        shouldEqual(res,2);

        vec(10)=3;
        res = vec(20)+vec(10);
        shouldEqual(res,5);

        res = vec(20)*vec(10);
        shouldEqual(res,6);
    }

    void testProxy4(){
        IntVec vec(100);

        std::vector<int> dvec(100,0);

        int res;


        res = vec[20]+dvec[20];
        shouldEqual(res,0);

        vec[20]=2;
        res = vec[20]+dvec[10];
        shouldEqual(res,2);

        dvec[10]=3;
        res =  vec[20]+dvec[10];
        shouldEqual(res,5);

        res = vec[20]*dvec[10];
        shouldEqual(res,6);
    }


    void testBatchConstructor1(){
        int indices[5]={10,20,30,40,50};
        int values[5]={1,2,3,4,5};
        IntVec vec(100,indices,indices+5,values);
        shouldEqual(vec.size(),100);
        for(size_t i=0;i<5;++i){
            shouldEqual(vec[indices[i]],values[i]);
        }
    }
    void testBatchConstructor2(){
        coordinate_value_pair cv[5]={
            coordinate_value_pair(10,1),
            coordinate_value_pair(20,2),
            coordinate_value_pair(30,3),
            coordinate_value_pair(40,4),
            coordinate_value_pair(50,5)
        };

        IntVec vec(100,cv,cv+5);
        shouldEqual(vec.size(),100);
        for(size_t i=0;i<5;++i){
            shouldEqual( vec[cv[i].first], cv[i].second );
        }
    }
    #endif
};





 
struct SparseArrayTestSuite
: public vigra::test_suite
{
    SparseArrayTestSuite()
    : vigra::test_suite("SparseArrayTestSuite")
    {   
        add( testCase( &SparseVectorTest::testProxy0));
        //add( testCase( &SparseVectorTest::testProxy1));
        //add( testCase( &SparseVectorTest::testProxy2));
        //add( testCase( &SparseVectorTest::testProxy3));
        //add( testCase( &SparseVectorTest::testProxy4));
        //add( testCase( &SparseVectorTest::testBatchConstructor1));
        //add( testCase( &SparseVectorTest::testBatchConstructor2));
    }
};

int main(int argc, char ** argv)
{
    SparseArrayTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

