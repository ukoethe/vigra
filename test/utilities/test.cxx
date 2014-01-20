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

#include <cstddef>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <queue>
#include <set>

#include "vigra/unittest.hxx"
#include "vigra/accessor.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/copyimage.hxx"
#include "vigra/sized_int.hxx"
#include "vigra/priority_queue.hxx"
#include "vigra/merge_graph/iterable_partition.hxx"

using namespace vigra;

struct ArrayVectorTest
{
    typedef int value_type;
    typedef ArrayVector<value_type> Vector;
    typedef Vector::iterator Iterator;
    typedef Vector::const_iterator ConstIterator;
    typedef StandardValueAccessor<value_type> Accessor;
    
    Vector vector_;

    ArrayVectorTest()
    {}
    
    void testAccessor()
    {
        vector_.resize(3, 0);
        Iterator i = vector_.begin();
        ConstIterator ci = const_cast<Vector const &>(vector_).begin();
        
        StandardValueAccessor<value_type> sva;
        StandardConstValueAccessor<value_type> scva;
        sva.set(3, i);
        shouldEqual(vector_[0], 3);
        shouldEqual(sva(i), 3);
        shouldEqual(scva(i), 3);
        shouldEqual(scva(ci), 3);
        sva.set(4, i+1);
        shouldEqual(vector_[1], 4);
        shouldEqual(sva(i+1), 4);
        shouldEqual(scva(i+1), 4);
        shouldEqual(scva(ci+1), 4);
        sva.set(5, i, 2);
        shouldEqual(vector_[2], 5);
        shouldEqual(sva(i, 2), 5);
        shouldEqual(scva(i, 2), 5);
        shouldEqual(scva(ci, 2), 5);
        
        StandardAccessor<value_type> sa;
        StandardConstAccessor<value_type> sca;
        sa.set(6, i);
        shouldEqual(vector_[0], 6);
        shouldEqual(sa(i), 6);
        shouldEqual(sca(i), 6);
        shouldEqual(sca(ci), 6);
        sa.set(7, i+1);
        shouldEqual(vector_[1], 7);
        shouldEqual(sa(i+1), 7);
        shouldEqual(sca(i+1), 7);
        shouldEqual(sca(ci+1), 7);
        sa.set(8, i, 2);
        shouldEqual(vector_[2], 8);
        shouldEqual(sa(i, 2), 8);
        shouldEqual(sca(i, 2), 8);
        shouldEqual(sca(ci, 2), 8);
        
        Vector varray[] = { vector_, vector_, vector_ };
        Vector  * v = varray;
        Vector const * cv = varray;
        int k;
        
        VectorComponentAccessor<Vector> vca(0);
        for(k = 0; k<2; ++k, vca.setIndex(k))
        {
            shouldEqual(vca(v), 6 + k);
            shouldEqual(vca(v + 1), 6 + k);
            shouldEqual(vca(v, 2), 6 + k);
            shouldEqual(vca(cv), 6 + k);
            shouldEqual(vca(cv + 1), 6 + k);
            shouldEqual(vca(cv, 2), 6 + k);
            vca.set(3 + k, v);
            vca.set(3 + k, v + 1);
            vca.set(3 + k, v, 2);
            shouldEqual(v[0][k], 3 + k);
            shouldEqual(v[1][k], 3 + k);
            shouldEqual(v[2][k], 3 + k);
        }
        
        VectorComponentValueAccessor<Vector> vcva(0);
        for(k = 0; k<2; ++k, vcva.setIndex(k))
        {
            shouldEqual(vcva(v), 3 + k);
            shouldEqual(vcva(v + 1), 3 + k);
            shouldEqual(vcva(v, 2), 3 + k);
            shouldEqual(vcva(cv), 3 + k);
            shouldEqual(vcva(cv + 1), 3 + k);
            shouldEqual(vcva(cv, 2), 3 + k);
            vcva.set(6 + k, v);
            vcva.set(6 + k, v + 1);
            vcva.set(6 + k, v, 2);
            shouldEqual(v[0][k], 6 + k);
            shouldEqual(v[1][k], 6 + k);
            shouldEqual(v[2][k], 6 + k);
        }

        VectorAccessor<Vector> va;
        VectorElementAccessor<VectorAccessor<Vector> > vea(0, va);
        for(k = 0; k<2; ++k, vea.setIndex(k))
        {
            shouldEqual(vea(v), 6 + k);
            shouldEqual(vea(v + 1), 6 + k);
            shouldEqual(vea(v, 2), 6 + k);
            shouldEqual(vea(cv), 6 + k);
            shouldEqual(vea(cv + 1), 6 + k);
            shouldEqual(vea(cv, 2), 6 + k);
            vea.set(3 + k, v);
            vea.set(3 + k, v + 1);
            vea.set(3 + k, v, 2);
            shouldEqual(v[0][k], 3 + k);
            shouldEqual(v[1][k], 3 + k);
            shouldEqual(v[2][k], 3 + k);
        }

        for(k = 0; k<2; ++k)
        {
            shouldEqual(va.getComponent(v, k), 3 + k);
            shouldEqual(va.getComponent(v + 1, k), 3 + k);
            shouldEqual(va.getComponent(v, 2, k), 3 + k);
            shouldEqual(va.getComponent(cv, k), 3 + k);
            shouldEqual(va.getComponent(cv + 1, k), 3 + k);
            shouldEqual(va.getComponent(cv, 2, k), 3 + k);
            va.setComponent(6 + k, v, k);
            va.setComponent(6 + k, v + 1, k);
            va.setComponent(6 + k, v, 2, k);
            shouldEqual(v[0][k], 6 + k);
            shouldEqual(v[1][k], 6 + k);
            shouldEqual(v[2][k], 6 + k);
        }

        SequenceAccessor<Vector> sqa;
        SequenceAccessor<const Vector> sqca;
        shouldEqual(sqa.size(v), 3u);
        shouldEqual(sqa.size(v + 1), 3u);
        shouldEqual(sqa.size(v, 2), 3u);
        shouldEqual(sqca.size(cv), 3u);
        shouldEqual(sqca.size(cv + 1), 3u);
        shouldEqual(sqca.size(cv, 2), 3u);
        should(sqa.end(v) == v[0].end());
        should(sqa.end(v + 1) == v[1].end());
        should(sqa.end(v, 2) == v[2].end());
        should(sqca.end(cv) == cv[0].end());
        should(sqca.end(cv + 1) == cv[1].end());
        should(sqca.end(cv, 2) == cv[2].end());
        for(k = 0; k<2; ++k)
        {
            shouldEqual(sqa.begin(v)[k], 6 + k);
            shouldEqual(sqa.begin(v + 1)[k], 6 + k);
            shouldEqual(sqa.begin(v, 2)[k], 6 + k);
            shouldEqual(sqca.begin(cv)[k], 6 + k);
            shouldEqual(sqca.begin(cv + 1)[k], 6 + k);
            shouldEqual(sqca.begin(cv, 2)[k], 6 + k);
            sqa.begin(v)[k] = 3 + k;
            sqa.begin(v + 1)[k] = 3 + k;
            sqa.begin(v, 2)[k] = 3 + k;
            shouldEqual(v[0][k], 3 + k);
            shouldEqual(v[1][k], 3 + k);
            shouldEqual(v[2][k], 3 + k);
        }
    }

    void testBackInsertion()
    {
        static value_type data[] = { 0, 1, 2, 3, 4 };
        
        shouldEqual(vector_.size(), 0u);
        
        Accessor a;
        copyLine(data, data + 5, a, std::back_inserter(vector_), a);
        
        shouldEqual(vector_.size(), 5u);
        shouldEqualSequence(vector_.begin(), vector_.end(), data);
    }

    void testAmbiguousConstructor()
    {
        ArrayVector<std::ptrdiff_t> a(2, std::ptrdiff_t(1));
        ArrayVector<std::ptrdiff_t> b(a.begin(), a.end());
    }
};


template<class ID_TYPE>
struct IterablePartitonTest
{
    typedef ID_TYPE IdType;
    typedef vigra::merge_graph_detail::IterablePartition<IdType> PartitionType;
    typedef std::set<IdType> SetType;
    typedef std::vector<IdType> VecType;
    IterablePartitonTest()
    {

    }

    void trueReps(const PartitionType ufd,SetType &set){
        set.clear();
        for(IdType i=0;i<ufd.numberOfElements();++i){
            if(ufd.find(i)==i){
                set.insert(i);
            }
        }
    }

    void trueReps(const PartitionType ufd,SetType &set,SetType & c){
        set.clear();
        for(IdType i=0;i<ufd.numberOfElements();++i){
            if(ufd.find(i)==i){
                set.insert(i);
            }
        }
        for(typename  SetType::const_iterator iter=c.begin();iter!=c.end();++iter){
            const IdType toRemove=*iter;
            should(set.find(toRemove)!=set.end());
            set.erase(toRemove);
        }

    }


    void testReps(const PartitionType ufd,VecType & vec){
        vec.clear();
        vec.assign(ufd.begin(),ufd.end());
    }




    void iteratorTest1(){
        PartitionType ufd(6);
        SetType trueRep;
        VecType testRep;

        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

            

        ufd.merge(0,1);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());
        
        ufd.merge(0,2);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());
        
        ufd.merge(0,3);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        ufd.merge(3,3);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        ufd.merge(4,5);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        ufd.merge(3,5);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());
    }

    void iteratorTest2(){
        PartitionType ufd(6);
        SetType trueRep;
        VecType testRep;

        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        SetType erased;
        erased.insert(0);
        ufd.eraseElement(0);

        trueReps(ufd,trueRep,erased);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        ufd.merge(1,2);
        trueReps(ufd,trueRep,erased);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        IdType rep12 = ufd.find(1);
        erased.insert(rep12);
        ufd.eraseElement(rep12);
        trueReps(ufd,trueRep,erased);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());
    }
    
};


struct BucketQueueTest
{
    struct Priority
    {
        int operator()(double v) const
        {
            return (int)v;
        }
    };
    
    ArrayVector<double> data;
    ArrayVector<int> idata;
    
    BucketQueueTest()
    {
        data.push_back(1.1);
        data.push_back(4.4);
        data.push_back(12.2);
        data.push_back(2.2);
        data.push_back(3.6);
        data.push_back(4.5);
        
        idata.resize(data.size());
        std::transform(data.begin(), data.end(), idata.begin(), Priority());
    }
    
    void testDescending()
    {
        std::priority_queue<int> queue;
        BucketQueue<int> bqueue;
        
        for(unsigned int k=0; k<idata.size(); ++k)
        {
            queue.push(idata[k]);
            bqueue.push(idata[k], idata[k]);
        }
        
        shouldEqual(idata.size(), bqueue.size());
        shouldEqual(false, bqueue.empty());
        
        for(unsigned int k=0; k<idata.size(); ++k)
        {
            shouldEqual(queue.top(), bqueue.top());
            queue.pop();
            bqueue.pop();
        }
        
        shouldEqual(0u, bqueue.size());
        shouldEqual(true, bqueue.empty());        
    }
    
    void testAscending()
    {
        std::priority_queue<int, std::vector<int>, std::greater<int> > queue;
        BucketQueue<int, true> bqueue;
        
        for(unsigned int k=0; k<idata.size(); ++k)
        {
            queue.push(idata[k]);
            bqueue.push(idata[k], idata[k]);
        }
        
        shouldEqual(idata.size(), bqueue.size());
        shouldEqual(false, bqueue.empty());
        
        for(unsigned int k=0; k<idata.size(); ++k)
        {
            shouldEqual(queue.top(), bqueue.top());
            queue.pop();
            bqueue.pop();
        }
        
        shouldEqual(0u, bqueue.size());
        shouldEqual(true, bqueue.empty());        
    }
    
    void testDescendingMapped()
    {
        Priority priority;
        std::priority_queue<int> queue;
        MappedBucketQueue<double, Priority> bqueue;
        
        for(unsigned int k=0; k<data.size(); ++k)
        {
            queue.push(idata[k]);
            bqueue.push(data[k]);
        }
        
        shouldEqual(data.size(), bqueue.size());
        shouldEqual(false, bqueue.empty());
        
        for(unsigned int k=0; k<data.size(); ++k)
        {
            shouldEqual(queue.top(), priority(bqueue.top()));
            switch(k)
            {
              case 1:
                  shouldEqual(4.4, bqueue.top());
                  break;
              case 2:
                  shouldEqual(4.5, bqueue.top());
                  break;
            }
            queue.pop();
            bqueue.pop();
        }
        
        shouldEqual(0u, bqueue.size());
        shouldEqual(true, bqueue.empty());        
    }
    
    void testAscendingMapped()
    {
        Priority priority;
        std::priority_queue<int, std::vector<int>, std::greater<int> > queue;
        MappedBucketQueue<double, Priority, true> bqueue;
        
        for(unsigned int k=0; k<data.size(); ++k)
        {
            queue.push(idata[k]);
            bqueue.push(data[k]);
        }
        
        shouldEqual(data.size(), bqueue.size());
        shouldEqual(false, bqueue.empty());
        
        for(unsigned int k=0; k<data.size(); ++k)
        {
            shouldEqual(queue.top(), priority(bqueue.top()));
            switch(k)
            {
              case 3:
                  shouldEqual(4.4, bqueue.top());
                  break;
              case 4:
                  shouldEqual(4.5, bqueue.top());
                  break;
            }
            queue.pop();
            bqueue.pop();
        }
        
        shouldEqual(0u, bqueue.size());
        shouldEqual(true, bqueue.empty());        
    }
};


struct ChangeablePriorityQueueTest
{
    typedef ChangeablePriorityQueue<float, std::less<float>    > MinQueueType;
    typedef ChangeablePriorityQueue<float, std::greater<float> > MaxQueueType;

    ChangeablePriorityQueueTest()
    {

    }


    void testMinQueue(){
        const float tol=0.001;
        {
            MinQueueType q(4);

            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => NONE
            // 2 => NONE
            // 3 => NONE
            should(q.empty());
            shouldEqual(q.size(),0);
            should(!q.contains(0));
            should(!q.contains(1));
            should(!q.contains(2));
            should(!q.contains(3));



            q.push(0,3.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => NONE
            // 2 => NONE
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),1);
            should( q.contains(0));
            should(!q.contains(1));
            should(!q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            //shouldEqualTolerance(q.priority(1),3.0, tol);
            //shouldEqualTolerance(q.priority(2),3.0, tol);
            //shouldEqualTolerance(q.priority(3),3.0, tol);
            shouldEqual(q.top(),0);
            shouldEqualTolerance(q.topPriority(),3.0,tol);


            q.push(2,2.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => NONE
            // 2 => 2.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),2);
            should( q.contains(0));
            should(!q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            //shouldEqualTolerance(q.priority(1),3.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),3.0, tol);
            shouldEqual(q.top(),2);
            shouldEqualTolerance(q.topPriority(),2.0,tol);


            q.push(1,3.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 3.0
            // 2 => 2.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),3);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),3.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),3.0, tol);
            shouldEqual(q.top(),2);
            shouldEqualTolerance(q.topPriority(),2.0,tol);


            q.push(3,0.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 3.0
            // 2 => 2.0
            // 3 => 0.0
            should(!q.empty());
            shouldEqual(q.size(),4);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should( q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),3.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            shouldEqualTolerance(q.priority(3),0.0, tol);
            shouldEqual(q.top(),3);
            shouldEqualTolerance(q.topPriority(),0.0,tol);

            q.pop();
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 3.0
            // 2 => 2.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),3);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),3.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),0.0, 0.01);
            shouldEqual(q.top(),2);
            shouldEqualTolerance(q.topPriority(),2.0,tol);


            q.push(1,1.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 1.0
            // 2 => 2.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),3);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),1.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),0.0, 0.01);
            shouldEqual(q.top(),1);
            shouldEqualTolerance(q.topPriority(),1.0,tol);

            q.push(1,4.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 4.0
            // 2 => 2.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),3);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),4.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),0.0, 0.01);
            shouldEqual(q.top(),2);
            shouldEqualTolerance(q.topPriority(),2.0,tol);


            q.push(0,1.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 1.0
            // 1 => 4.0
            // 2 => 2.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),3);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),1.0, tol);
            shouldEqualTolerance(q.priority(1),4.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),0.0, 0.01);
            shouldEqual(q.top(),0);
            shouldEqualTolerance(q.topPriority(),1.0,tol);


            q.pop();
            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => 4.0
            // 2 => 2.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),2);
            should(!q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            //shouldEqualTolerance(q.priority(0),1.0, tol);
            shouldEqualTolerance(q.priority(1),4.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),0.0, 0.01);
            shouldEqual(q.top(),2);
            shouldEqualTolerance(q.topPriority(),2.0,tol);


            q.pop();
            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => 4.0
            // 2 => NONE
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),1);
            should(!q.contains(0));
            should( q.contains(1));
            should(!q.contains(2));
            should(!q.contains(3));
            //shouldEqualTolerance(q.priority(0),1.0, tol);
            shouldEqualTolerance(q.priority(1),4.0, tol);
            //shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),0.0, 0.01);
            shouldEqual(q.top(),1);
            shouldEqualTolerance(q.topPriority(),4.0,tol);


            q.pop();
            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => NONE
            // 2 => NONE
            // 3 => NONE
            should(q.empty());
            shouldEqual(q.size(),0);
            should(!q.contains(0));
            should(!q.contains(1));
            should(!q.contains(2));
            should(!q.contains(3));
            //shouldEqualTolerance(q.priority(0),1.0, tol);
            //shouldEqualTolerance(q.priority(1),4.0, tol);
            //shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),0.0, 0.01);
            //shouldEqual(q.top(),1);
            //shouldEqualTolerance(q.topPriority(),4.0,tol);

        }

    }

    void testMaxQueue(){
        const float tol=0.001;
        {
            MaxQueueType q(4);

            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => NONE
            // 2 => NONE
            // 3 => NONE
            should(q.empty());
            shouldEqual(q.size(),0);
            should(!q.contains(0));
            should(!q.contains(1));
            should(!q.contains(2));
            should(!q.contains(3));



            q.push(0,3.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => NONE
            // 2 => NONE
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),1);
            should( q.contains(0));
            should(!q.contains(1));
            should(!q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            //shouldEqualTolerance(q.priority(1),3.0, tol);
            //shouldEqualTolerance(q.priority(2),3.0, tol);
            //shouldEqualTolerance(q.priority(3),3.0, tol);
            shouldEqual(q.top(),0);
            shouldEqualTolerance(q.topPriority(),3.0,tol);


            q.push(2,2.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => NONE
            // 2 => 2.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),2);
            should( q.contains(0));
            should(!q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            //shouldEqualTolerance(q.priority(1),3.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),3.0, tol);
            shouldEqual(q.top(),0);
            shouldEqualTolerance(q.topPriority(),3.0,tol);


            q.push(1,4.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 4.0
            // 2 => 2.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),3);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),4.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            //shouldEqualTolerance(q.priority(3),3.0, tol);
            shouldEqual(q.top(),1);
            shouldEqualTolerance(q.topPriority(),4.0,tol);


            q.push(3,5.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 4.0
            // 2 => 2.0
            // 3 => 5.0
            should(!q.empty());
            shouldEqual(q.size(),4);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should( q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),4.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            shouldEqualTolerance(q.priority(3),5.0, tol);
            shouldEqual(q.top(),3);
            shouldEqualTolerance(q.topPriority(),5.0,tol);


            q.push(3,2.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 4.0
            // 2 => 2.0
            // 3 => 2.0
            should(!q.empty());
            shouldEqual(q.size(),4);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should( q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),4.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            shouldEqualTolerance(q.priority(3),2.0, tol);
            shouldEqual(q.top(),1);
            shouldEqualTolerance(q.topPriority(),4.0,tol);


            q.push(1,0.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 0.0
            // 2 => 2.0
            // 3 => 2.0
            should(!q.empty());
            shouldEqual(q.size(),4);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should( q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),0.0, tol);
            shouldEqualTolerance(q.priority(2),2.0, tol);
            shouldEqualTolerance(q.priority(3),2.0, tol);
            shouldEqual(q.top(),0);
            shouldEqualTolerance(q.topPriority(),3.0,tol);

            q.push(2,5.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 0.0
            // 2 => 5.0
            // 3 => 2.0
            should(!q.empty());
            shouldEqual(q.size(),4);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should( q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),0.0, tol);
            shouldEqualTolerance(q.priority(2),5.0, tol);
            shouldEqualTolerance(q.priority(3),2.0, tol);
            shouldEqual(q.top(),2);
            shouldEqualTolerance(q.topPriority(),5.0,tol);

            q.push(3,1.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 0.0
            // 2 => 5.0
            // 3 => 1.0
            should(!q.empty());
            shouldEqual(q.size(),4);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should( q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),0.0, tol);
            shouldEqualTolerance(q.priority(2),5.0, tol);
            shouldEqualTolerance(q.priority(3),1.0, tol);
            shouldEqual(q.top(),2);
            shouldEqualTolerance(q.topPriority(),5.0,tol);


            q.pop();
            // CURRENT VALUES
            //-----------------
            // 0 => 3.0
            // 1 => 0.0
            // 2 => NONE
            // 3 => 1.0
            should(!q.empty());
            shouldEqual(q.size(),3);
            should( q.contains(0));
            should( q.contains(1));
            should(!q.contains(2));
            should( q.contains(3));
            shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),0.0, tol);
            //shouldEqualTolerance(q.priority(2),5.0, tol);
            shouldEqualTolerance(q.priority(3),1.0, tol);
            shouldEqual(q.top(),0);
            shouldEqualTolerance(q.topPriority(),3.0,tol);



            q.pop();
            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => 0.0
            // 2 => NONE
            // 3 => 1.0
            should(!q.empty());
            shouldEqual(q.size(),2);
            should(!q.contains(0));
            should( q.contains(1));
            should(!q.contains(2));
            should( q.contains(3));
            //shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),0.0, tol);
            //shouldEqualTolerance(q.priority(2),5.0, tol);
            shouldEqualTolerance(q.priority(3),1.0, tol);
            shouldEqual(q.top(),3);
            shouldEqualTolerance(q.topPriority(),1.0,tol);


            q.pop();
            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => 0.0
            // 2 => NONE
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),1);
            should(!q.contains(0));
            should( q.contains(1));
            should(!q.contains(2));
            should(!q.contains(3));
            //shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),0.0, tol);
            //shouldEqualTolerance(q.priority(2),5.0, tol);
            //shouldEqualTolerance(q.priority(3),1.0, tol);
            shouldEqual(q.top(),1);
            shouldEqualTolerance(q.topPriority(),0.0,tol);

            q.pop();
            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => NONE
            // 2 => NONE
            // 3 => NONE
            should(q.empty());
            shouldEqual(q.size(),0);
            should(!q.contains(0));
            should(!q.contains(1));
            should(!q.contains(2));
            should(!q.contains(3));
            //shouldEqualTolerance(q.priority(0),3.0, tol);
            //shouldEqualTolerance(q.priority(1),0.0, tol);
            //shouldEqualTolerance(q.priority(2),5.0, tol);
            //shouldEqualTolerance(q.priority(3),1.0, tol);
            //shouldEqual(q.top(),1);
            //shouldEqualTolerance(q.topPriority(),0.0,tol);

            q.push(2,1.0);
            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => NONE
            // 2 => 1.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),1);
            should(!q.contains(0));
            should(!q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            //shouldEqualTolerance(q.priority(0),3.0, tol);
            //shouldEqualTolerance(q.priority(1),0.0, tol);
            shouldEqualTolerance(q.priority(2),1.0, tol);
            //shouldEqualTolerance(q.priority(3),1.0, tol);
            shouldEqual(q.top(),2);
            shouldEqualTolerance(q.topPriority(),1.0,tol);


            q.push(2,3.0);
            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => NONE
            // 2 => 3.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),1);
            should(!q.contains(0));
            should(!q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            //shouldEqualTolerance(q.priority(0),3.0, tol);
            //shouldEqualTolerance(q.priority(1),0.0, tol);
            shouldEqualTolerance(q.priority(2),3.0, tol);
            //shouldEqualTolerance(q.priority(3),1.0, tol);
            shouldEqual(q.top(),2);
            shouldEqualTolerance(q.topPriority(),3.0,tol);


            q.push(1,4.0);
            // CURRENT VALUES
            //-----------------
            // 0 => NONE
            // 1 => 4.0
            // 2 => 3.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),2);
            should(!q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            //shouldEqualTolerance(q.priority(0),3.0, tol);
            shouldEqualTolerance(q.priority(1),4.0, tol);
            shouldEqualTolerance(q.priority(2),3.0, tol);
            //shouldEqualTolerance(q.priority(3),1.0, tol);
            shouldEqual(q.top(),1);
            shouldEqualTolerance(q.topPriority(),4.0,tol);


            q.push(0,2.0);
            // CURRENT VALUES
            //-----------------
            // 0 => 2.0
            // 1 => 4.0
            // 2 => 3.0
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),3);
            should( q.contains(0));
            should( q.contains(1));
            should( q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),2.0, tol);
            shouldEqualTolerance(q.priority(1),4.0, tol);
            shouldEqualTolerance(q.priority(2),3.0, tol);
            //shouldEqualTolerance(q.priority(3),1.0, tol);
            shouldEqual(q.top(),1);
            shouldEqualTolerance(q.topPriority(),4.0,tol);


            q.deleteItem(2);
            // CURRENT VALUES
            //-----------------
            // 0 => 2.0
            // 1 => 4.0
            // 2 => NONE
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),2);
            should( q.contains(0));
            should( q.contains(1));
            should(!q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),2.0, tol);
            shouldEqualTolerance(q.priority(1),4.0, tol);
            //shouldEqualTolerance(q.priority(2),3.0, tol);
            //shouldEqualTolerance(q.priority(3),1.0, tol);
            shouldEqual(q.top(),1);
            shouldEqualTolerance(q.topPriority(),4.0,tol);


            q.deleteItem(1);
            // CURRENT VALUES
            //-----------------
            // 0 => 2.0
            // 1 => 4.0
            // 2 => NONE
            // 3 => NONE
            should(!q.empty());
            shouldEqual(q.size(),1);
            should( q.contains(0));
            should(!q.contains(1));
            should(!q.contains(2));
            should(!q.contains(3));
            shouldEqualTolerance(q.priority(0),2.0, tol);
            //shouldEqualTolerance(q.priority(1),4.0, tol);
            //shouldEqualTolerance(q.priority(2),3.0, tol);
            //shouldEqualTolerance(q.priority(3),1.0, tol);
            shouldEqual(q.top(),0);
            shouldEqualTolerance(q.topPriority(),2.0,tol);

        }

    }
};



struct SizedIntTest
{
    void testSizedInt()
    {
        shouldEqual(sizeof(Int8), 1u);
        shouldEqual(sizeof(Int16), 2u);
        shouldEqual(sizeof(Int32), 4u);
        shouldEqual(sizeof(UInt8), 1u);
        shouldEqual(sizeof(UInt16), 2u);
        shouldEqual(sizeof(UInt32), 4u);
        should(sizeof(IntBiggest) >= 4u);
        should(sizeof(UIntBiggest) >= 4u);
    }
};

struct MetaprogrammingTest
{
    struct TrueResult {};
    struct FalseResult {};
    
    struct Derived 
    : public TrueResult 
    {
        typedef TrueResult result_type;
        typedef TrueResult value_type;
    };
    
    void testInt()
    {
        shouldEqual(MetaInt<1>::value, 1);
        shouldEqual((MetaMax<1,2>::value), 2);
        shouldEqual((MetaMin<1,2>::value), 1);
    }
 
    void testLogic()
    {
        shouldEqual(VigraTrueType::value, true);
        shouldEqual(VigraFalseType::value, false);        
        should(typeid(Not<VigraFalseType>::type) == typeid(VigraTrueType));
        should(typeid(Not<VigraTrueType>::type) == typeid(VigraFalseType));
        should(typeid(And<VigraTrueType, VigraTrueType>::type) == typeid(VigraTrueType));
        should(typeid(And<VigraTrueType, VigraFalseType>::type) == typeid(VigraFalseType));
        should(typeid(And<VigraFalseType, VigraTrueType>::type) == typeid(VigraFalseType));
        should(typeid(And<VigraFalseType, VigraFalseType>::type) == typeid(VigraFalseType));
        should(typeid(Or<VigraTrueType, VigraTrueType>::type) == typeid(VigraTrueType));
        should(typeid(Or<VigraTrueType, VigraFalseType>::type) == typeid(VigraTrueType));
        should(typeid(Or<VigraFalseType, VigraTrueType>::type) == typeid(VigraTrueType));
        should(typeid(Or<VigraFalseType, VigraFalseType>::type) == typeid(VigraFalseType));
        should(typeid(If<VigraTrueType, TrueResult, FalseResult>::type) == typeid(TrueResult));
        should(typeid(If<VigraFalseType, TrueResult, FalseResult>::type) == typeid(FalseResult));
        should(typeid(IfBool<true, TrueResult, FalseResult>::type) == typeid(TrueResult));
        should(typeid(IfBool<false, TrueResult, FalseResult>::type) == typeid(FalseResult));
    }
 
    void testTypeTools()
    {
        should(typeid(IsSameType<TrueResult, TrueResult>::type) == typeid(VigraTrueType));
        should(typeid(IsSameType<TrueResult, FalseResult>::type) == typeid(VigraFalseType));
        should(typeid(IsDifferentType<TrueResult, TrueResult>::type) == typeid(VigraFalseType));
        should(typeid(IsDifferentType<TrueResult, FalseResult>::type) == typeid(VigraTrueType));
        should(typeid(IsConvertibleTo<int, double>::type) == typeid(VigraTrueType));
        should(typeid(IsConvertibleTo<int, FalseResult>::type) == typeid(VigraFalseType));
        should(typeid(IsDerivedFrom<Derived, TrueResult>::type) == typeid(VigraTrueType));
        should(typeid(IsDerivedFrom<Derived, FalseResult>::type) == typeid(VigraFalseType));
        should(typeid(has_result_type<Derived>::type) == typeid(VigraTrueType));
        should(typeid(has_result_type<FalseResult>::type) == typeid(VigraFalseType));
        should(typeid(has_value_type<Derived>::type) == typeid(VigraTrueType));
        should(typeid(has_value_type<FalseResult>::type) == typeid(VigraFalseType));

        should(typeid(IsIterator<std::reverse_iterator<int*> >::type) == typeid(VigraTrueType));
        should(typeid(IsIterator<int*>::type) == typeid(VigraTrueType));
        should(typeid(IsIterator<int const*>::type) == typeid(VigraTrueType));
        should(typeid(IsIterator<FalseResult>::type) == typeid(VigraFalseType));

        should(typeid(UnqualifiedType<int>::type) == typeid(int));
        should(typeid(UnqualifiedType<const int>::type) == typeid(int));
        should(typeid(UnqualifiedType<int*>::type) == typeid(int));
        should(typeid(UnqualifiedType<const int*>::type) == typeid(int));
        should(typeid(UnqualifiedType<int&>::type) == typeid(int));
        should(typeid(UnqualifiedType<const int&>::type) == typeid(int));
        should(typeid(UnqualifiedType<int**>::type) == typeid(int));
        should(typeid(UnqualifiedType<const int**>::type) == typeid(int));
        should(typeid(UnqualifiedType<int*&>::type) == typeid(int));
        should(typeid(UnqualifiedType<const int*&>::type) == typeid(int));
    }
};

void stringTest()
{
    std::string s;
    s << "Hallo " << 1 << " " << 2.0 << " " << false;
    shouldEqual(s, std::string("Hallo 1 2 0"));

    shouldEqual(asString(1), "1");
    shouldEqual(asString(2.0), "2");
    shouldEqual(asString(false), "0");

    shouldEqual(tolower("AluFr iNsta< Z89>"), "alufr insta< z89>");
    shouldEqual(normalizeString("AluFr iNsta< Z89>"), "alufrinsta<z89>");
}

struct UtilitiesTestSuite
: public vigra::test_suite
{
    UtilitiesTestSuite()
    : vigra::test_suite("UtilitiesTestSuite")
    {
        add( testCase( &ArrayVectorTest::testAccessor));
        add( testCase( &ArrayVectorTest::testBackInsertion));
        add( testCase( &ArrayVectorTest::testAmbiguousConstructor));
        add( testCase( &IterablePartitonTest<UInt32>::iteratorTest1));
        add( testCase( &IterablePartitonTest<UInt32>::iteratorTest2));
        add( testCase( &IterablePartitonTest<Int32>::iteratorTest1));
        add( testCase( &IterablePartitonTest<Int32>::iteratorTest2));
        add( testCase( &BucketQueueTest::testDescending));
        add( testCase( &BucketQueueTest::testAscending));
        add( testCase( &BucketQueueTest::testDescendingMapped));
        add( testCase( &BucketQueueTest::testAscendingMapped));
        add( testCase( &ChangeablePriorityQueueTest::testMinQueue));
        add( testCase( &ChangeablePriorityQueueTest::testMaxQueue));
        add( testCase( &SizedIntTest::testSizedInt));
        add( testCase( &MetaprogrammingTest::testInt));
        add( testCase( &MetaprogrammingTest::testLogic));
        add( testCase( &MetaprogrammingTest::testTypeTools));
        add( testCase( &stringTest));
    }
};

int main(int argc, char ** argv)
{
    UtilitiesTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

