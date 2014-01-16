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


    public:
        // Create an empty MinIndexedPQ which can contain atmost maxSize_ elements
        MinIndexedPQ(const size_t maxSize)  
        : maxSize_(maxSize),
          currentSize_(0),
          heap(maxSize_+1),
          index(maxSize_+1),
          values(maxSize_+1)
        {
            for(int i = 0; i <= maxSize_; i++)
                index[i] = -1;
        }
     
        // check if the PQ is empty
        bool isEmpty() const {
            return currentSize_ == 0;
        }
     
        // check if i is an index on the PQ
        bool contains(int i)    {
            return index[i] != -1;
        }
     
        // return the number of elements in the PQ
        int size()  {
            return currentSize_;
        }
     
        // associate value with index i; 0 < i < maxSize_
        void insert(int i, T value) {
            currentSize_++;
            index[i] = currentSize_;
            heap[currentSize_] = i;
            values[i] = value;
            bubbleUp(currentSize_);
        }
     
        // returns the index associated with the minimal value
        int minIndex() const {
            return heap[1];
        }
     
        // returns the minimal value
        T minValue() const {
            return values[heap[1]];
        }
     
        // delete the minimal value and return its associated index
        // Warning: Don't try to read from this index after calling this function
        int deleteMin() {
            int min = heap[1];
            swapItems(1, currentSize_--);
            bubbleDown(1);
            index[min] = -1;
            heap[currentSize_+1] = -1;
            return min;
        }
     
        // returns the value associated with index i
        T valueOf(const int i)    {
            return values[i];
        }
     
        // change the value associated with index i to the specified value
        void changeValue(const int i,const T value)  {
            values[i] = value;
            bubbleUp(index[i]);
            bubbleDown(index[i]);
        }
     
        // decrease the value associated with index i to the specified value
        void decreaseValue(const int i,const T value)    {
            values[i] = value;
            bubbleUp(index[i]);
        }
     
        // increase the value associated with index i to the specified value
        void increaseValue(const int i,const T value)    {
            values[i] = value;
            bubbleDown(index[i]);
        }
     
        // delete the value associated with index i
        void deleteValue(const int i)   {
            int ind = index[i];
            swapItems(ind, currentSize_--);
            bubbleUp(ind);
            bubbleDown(ind);
            index[i] = -1;
        }
    
    private:

        void swapItems(const int i,const  int j) {
            std::swap(heap[i],heap[j]);
            index[heap[i]] = i; 
            index[heap[j]] = j;
        }
     
        void bubbleUp(int k)    {
            while(k > 1 && values[heap[k/2]] > values[heap[k]])   {
                swapItems(k, k/2);
                k = k/2;
            }
        }
     
        void bubbleDown(int k)  {
            int j;
            while(2*k <= currentSize_) {
                j = 2*k;
                if(j < currentSize_ && values[heap[j]] > values[heap[j+1]])
                    j++;
                if(values[heap[k]] <= values[heap[j]])
                    break;
                swapItems(k, j);
                k = j;
            }
        }
     

        size_t maxSize_,currentSize_;
        std::vector<int> heap,index;
        std::vector<T>   values;

    };
}

#endif // VIGRA_MIN_INDEXED_PQ_HXX