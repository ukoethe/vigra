#ifndef VIGRA_MIN_INDEXED_PQ_HXX
#define VIGRA_MIN_INDEXED_PQ_HXX
/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>


/*
    This is a modified and templated version of code
    which I (Thorsten Beier) found here :
    http://www.codechef.com/viewsolution/2572064
    From this side i am guessing that "Kartik Kukreja "
    should get the credit for the inital implementation
*/

namespace vigra{


    



    template<class T>
    class MinIndexedPQ {
        int NMAX, N, *heap, *index,; //*keys;
        T * keys;
        void swapItems(int i, int j) {
            int t = heap[i]; heap[i] = heap[j]; heap[j] = t;
            index[heap[i]] = i; index[heap[j]] = j;
        }
     
        void bubbleUp(int k)    {
            while(k > 1 && keys[heap[k/2]] > keys[heap[k]])   {
                swapItems(k, k/2);
                k = k/2;
            }
        }
     
        void bubbleDown(int k)  {
            int j;
            while(2*k <= N) {
                j = 2*k;
                if(j < N && keys[heap[j]] > keys[heap[j+1]])
                    j++;
                if(keys[heap[k]] <= keys[heap[j]])
                    break;
                swapItems(k, j);
                k = j;
            }
        }
     
    public:
        // Create an empty MinIndexedPQ which can contain atmost NMAX elements
        MinIndexedPQ(int NMAX)  {
            this->NMAX = NMAX;
            N = 0;
            keys = new T[NMAX + 1];
            heap = new int[NMAX + 1];
            index = new int[NMAX + 1];
            for(int i = 0; i <= NMAX; i++)
                index[i] = -1;
        }
           
        ~MinIndexedPQ() {
            delete [] keys;
            delete [] heap;
            delete [] index;
        }
     
        // check if the PQ is empty
        bool isEmpty() const {
            return N == 0;
        }
     
        // check if i is an index on the PQ
        bool contains(int i)    {
            return index[i] != -1;
        }
     
        // return the number of elements in the PQ
        int size()  {
            return N;
        }
     
        // associate key with index i; 0 < i < NMAX
        void insert(int i, T key) {
            N++;
            index[i] = N;
            heap[N] = i;
            keys[i] = key;
            bubbleUp(N);
        }
     
        // returns the index associated with the minimal key
        int minIndex() const {
            return heap[1];
        }
     
        // returns the minimal key
        T minKey()    {
            return keys[heap[1]];
        }
     
        // delete the minimal key and return its associated index
        // Warning: Don't try to read from this index after calling this function
        int deleteMin() {
            int min = heap[1];
            swapItems(1, N--);
            bubbleDown(1);
            index[min] = -1;
            heap[N+1] = -1;
            return min;
        }
     
        // returns the key associated with index i
        T keyOf(int i)    {
            return keys[i];
        }
     
        // change the key associated with index i to the specified value
        void changeKey(int i, T key)  {
            keys[i] = key;
            bubbleUp(index[i]);
            bubbleDown(index[i]);
        }
     
        // decrease the key associated with index i to the specified value
        void decreaseKey(int i, T key)    {
            keys[i] = key;
            bubbleUp(index[i]);
        }
     
        // increase the key associated with index i to the specified value
        void increaseKey(int i, T key)    {
            keys[i] = key;
            bubbleDown(index[i]);
        }
     
        // delete the key associated with index i
        void deleteKey(int i)   {
            int ind = index[i];
            swapItems(ind, N--);
            bubbleUp(ind);
            bubbleDown(ind);
            index[i] = -1;
        }
    };

    /*
    struct FloatPq{

        FloatPq()
        :  pq_(0)
        {

        }
        template<class ITER>
        FloatPq(ITER begin,ITER end,const float scale=1000000)
        :   scale_(scale),
            initSize_(std::distance(begin,end)),
            values_(begin,end),
            iValues_(std::distance(begin,end)),
            activeKey_(std::distance(begin,end),1),
            pq_(std::distance(begin,end))
        {
            for(size_t index=0;index<initSize_;++index){
                const float rawVal = values_[index];

                const float scaled = rawVal*scale_;
                const vigra::UInt64 iVal  = static_cast<vigra::UInt64>(scaled+0.5);
                iValues_[index]=iVal;

                // register iVal as key for index
                pq_.insert(index,rawVal);
            }
        }

        FloatPq(const size_t size,const float scale=1000000)
        :   scale_(scale),
            initSize_(size),
            values_(size),
            iValues_(size),
            activeKey_(size,1),
            pq_(size)
        {

        }

        template<class ITER>
        void setValues(ITER begin, ITER end)
        {
            values_.assign(begin,end);
            for(size_t index=0;index<initSize_;++index){
                const float rawVal = values_[index];


                const float scaled = rawVal*scale_;
                const vigra::UInt64 iVal  = static_cast<vigra::UInt64>(scaled+0.5);
                iValues_[index]=iVal;

                // register iVal as key for index
                pq_.insert(index,rawVal);
            }
        }

        void changeValue(const size_t index,const float newVal){
            
            CGP_ASSERT_OP(activeKey_[index],<=,1);
            CGP_ASSERT_OP(activeKey_[index],>=,1);
            CGP_ASSERT_OP(pq_.contains(index),==,true);
            const float oldVal =values_[index];
            values_[index]=newVal;
            const float scaled = newVal*scale_;
            const vigra::UInt64 iValNew  = static_cast<vigra::UInt64>(scaled+0.5);
            const vigra::UInt64 iValOld  = iValues_[index];
            pq_.changeKey(index,newVal);
            
            //if(iValNew<iValOld){
            //    pq_.decreaseKey(index,iValNew);
            //}
            //else if (iValNew>iValOld){
            //    pq_.increaseKey(index,iValNew);
            //}
            //else{
            //    // do nothing if no changes
            //}
            

        }

        size_t minIndex()const{
            return pq_.minIndex();
        }
        void deleteIndex(const size_t index){
            CGP_ASSERT_OP(pq_.contains(index),==,true);
            CGP_ASSERT_OP(activeKey_[index],==,1);
            pq_.deleteKey(index);
            activeKey_[index]=0;
        }

        bool hasIndex(const size_t index){
            CGP_ASSERT_OP( bool(pq_.contains(index)) , == , bool(activeKey_[index]==1));
            return pq_.contains(index);
        }

        float scale_;
        size_t                initSize_;

        std::vector<float>          values_;
        std::vector<vigra::UInt64>  iValues_;
        std::vector<unsigned char>  activeKey_;
        MinIndexedPQ<float>         pq_;
    };
    */

}

#endif // VIGRA_MIN_INDEXED_PQ_HXX