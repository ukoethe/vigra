/************************************************************************/
/*                                                                      */
/*        Copyright 2014-2015 by Thorsten Beier, Philip Schill          */
/*                               and Ullrich Koethe                     */
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
#ifndef VIGRA_THREADPOOL_HXX
#define VIGRA_THREADPOOL_HXX

#include <functional>
#include <thread>
#include <atomic>
#include <vector>
#include <future>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <stdexcept>
#include <cmath>
#include "mathutil.hxx"



namespace vigra
{



class ThreadPool {
public:

    /**
     * Create a thread pool with n threads. The constructor just launches some workers.
     */
    ThreadPool(size_t n);

    /**
     * The destructor joins all threads.
     */
    ~ThreadPool();

    /**
     * Enqueue a task that will be executed by the thread pool.
     * The task result can be obtained using the get() function of the returned future.
     * If the task throws an exception, it will be raised on the call to get().
     */
    template<class F>
    std::future<typename std::result_of<F(int)>::type>  enqueueReturning(F&& f) ;

    /**
     * Enqueue function for tasks without return value.
     * This is a special case of the enqueueReturning template function, but
     * some compilers fail on std::result_of<F(int)>::type for void(int)functions.
     */
    template<class F>
    std::future<void> enqueue(F&& f) ;

    /**
     * Block until all tasks are finished.
     */
    void waitFinished()
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        finish_condition.wait(lock, [this](){ return tasks.empty() && (busy == 0); });
    }

    /**
     * Return the number of worker threads.
     */
    size_t nThreads() const
    {
        return workers.size();
    }

private:

    // need to keep track of threads so we can join them
    std::vector<std::thread> workers;

    // the task queue
    std::queue<std::function<void(int)> > tasks;
    
    // synchronization
    std::mutex queue_mutex;
    std::condition_variable worker_condition;
    std::condition_variable finish_condition;
    bool stop;
    std::atomic<unsigned int> busy, processed;
};

inline ThreadPool::ThreadPool(size_t threads)
    :   stop(false),
        busy(0),
        processed(0)
{
    vigra_precondition(threads > 0, "ThreadPool::ThreadPool(): n_threads must not be zero.");
    for(size_t ti = 0; ti<threads; ++ti)
    {
        workers.emplace_back(
            [ti,this]
            {
                for(;;)
                {
                    std::function<void(int)> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);

                        // will wait if : stop == false  AND queue is empty
                        // if stop == true AND queue is empty thread function will return later
                        //
                        // so the idea of this wait, is : If where are not in the destructor
                        // (which sets stop to true, we wait here for new jobs)
                        this->worker_condition.wait(lock, [this]{ return this->stop || !this->tasks.empty(); });
                        if(!this->tasks.empty())
                        {
                            ++busy;
                            task = std::move(this->tasks.front());
                            this->tasks.pop();
                            lock.unlock();
                            task(ti);
                            ++processed;
                            --busy;
                            finish_condition.notify_one();
                        }
                        else if(stop)
                        {
                            return;
                        }
                    }
                }
            }
        );
    }
}

inline ThreadPool::~ThreadPool()
{
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    worker_condition.notify_all();
    for(std::thread &worker: workers)
        worker.join();
}

template<class F>
inline std::future<typename std::result_of<F(int)>::type> 
ThreadPool::enqueueReturning(F&& f)
{
    typedef typename std::result_of<F(int)>::type result_type;
    typedef std::packaged_task<result_type(int)> PackageType;

    auto task = std::make_shared<PackageType>(f);
    auto res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace(
            [task](int tid)
            {
                (*task)(tid);
            }
        );
    }
    worker_condition.notify_one();
    return res;
}

template<class F>
inline std::future<void>
ThreadPool::enqueue(F&& f)
{
    typedef std::packaged_task<void(int)> PackageType;

    auto task = std::make_shared<PackageType>(f);
    auto res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace(
            [task](int tid)
            {
                (*task)(tid);
            }
        );
    }
    worker_condition.notify_one();
    return res;
}



// nItems must be either zero or std::distance(iter, end).
template<class ITER, class F>
inline void parallel_foreach_impl(
    ThreadPool & pool,
    const uint64_t nItems,                   
    ITER iter, 
    ITER end, 
    F && f,
    std::random_access_iterator_tag
){
    uint64_t workload = std::distance(iter, end);
    vigra_precondition(workload == nItems || nItems == 0, "parallel_foreach(): Mismatch between num items and begin/end.");
    const float workPerThread = float(workload)/pool.nThreads();
    const uint64_t chunkedWorkPerThread = std::max<uint64_t>(roundi(workPerThread/3.0), 1);

    std::vector<std::future<void> > futures;
    for( ;iter<end; iter+=chunkedWorkPerThread)
    {
        const size_t lc = std::min(workload, chunkedWorkPerThread);
        workload-=lc;
        futures.emplace_back(
            pool.enqueue(
                [&f, iter, lc]
                (int id)
                {
                    for(size_t i=0; i<lc; ++i)
                        f(id, iter[i]);
                }
            )
        );
    }
    for (auto & fut : futures)
        fut.get();
}



// nItems must be either zero or std::distance(iter, end).
template<class ITER, class F>
inline void parallel_foreach_impl(
    ThreadPool & pool,
    const uint64_t nItems,                   
    ITER iter, 
    ITER end, 
    F && f,
    std::forward_iterator_tag
){
    if (nItems == 0)
        nItems = std::distance(iter, end);

    uint64_t workload = nItems;
    const float workPerThread = float(workload)/pool.nThreads();
    const uint64_t chunkedWorkPerThread = std::max<uint64_t>(roundi(workPerThread/3.0), 1);

    std::vector<std::future<void> > futures;
    for(;;)
    {
        const size_t lc = std::min(chunkedWorkPerThread, workload);
        workload -= lc;
        futures.emplace_back(
            pool.enqueue(
                [&f, iter, lc]
                (int id)
                {
                    auto iterCopy = iter;
                    for(size_t i=0; i<lc; ++i){
                        f(id, *iterCopy);
                        ++iterCopy;
                    }
                }
            )
        );
        for (size_t i = 0; i < lc; ++i)
        {
            ++iter;
            if (iter == end)
            {
                vigra_postcondition(workload == 0, "parallel_foreach(): Mismatch between num items and begin/end.");
                break;
            }
        }
        if(workload==0)
            break;
    }
    for (auto & fut : futures)
        fut.get();
}



// nItems must be either zero or std::distance(iter, end).
template<class ITER, class F>
inline void parallel_foreach_impl(
    ThreadPool & pool,
    const uint64_t nItems,
    ITER iter, 
    ITER end, 
    F && f,
    std::input_iterator_tag
){
    size_t num_items = 0;
    std::vector<std::future<void> > futures;
    for (; iter != end; ++iter)
    {
        auto item = *iter;
        futures.emplace_back(
            pool.enqueue(
                [&f, &item](int id){
                    f(id, item);
                }
            )
        );
        ++num_items;
    }
    vigra_postcondition(num_items == nItems || nItems == 0, "parallel_foreach(): Mismatch between num items and begin/end.");
    for (auto & fut : futures)
        fut.get();
}



/**
 * Runs the parallel foreach on a single thread.
 */
template<class ITER, class F>
inline void parallel_foreach_single_thread(
    ITER begin, 
    ITER end, 
    F && f,
    const uint64_t nItems = 0
){
    size_t n = 0;
    for (; begin != end; ++begin)
    {
        f(0, *begin);
        ++n;
    }
    vigra_postcondition(n == nItems || nItems == 0, "parallel_foreach(): Mismatch between num items and begin/end.");
}



/**
 * Just like the other parallel_foreach overload, but use the given threadpool
 * instead of creating a new one.
 */
template<class ITER, class F>
inline void parallel_foreach(
    ThreadPool & pool,
    ITER begin, 
    ITER end, 
    F && f,
    const uint64_t nItems = 0
){
    if(pool.nThreads()>1)
    {
        parallel_foreach_impl(pool,nItems, begin, end, f,
            typename std::iterator_traits<ITER>::iterator_category());
    }
    else
    {
        parallel_foreach_single_thread(begin, end, f, nItems);
    }
}



/**
 * Create a threadpool to apply the functor F to all items in [begin, end) in
 * parallel. F must be callable with two arguments of type size_t and T, where
 * the first argument is the thread index (starting by 0) and T is the value
 * type of the iterators.
 * 
 * If ITER is a forward iterator (std::forward_iterator_tag) and the optional
 * argument nItems is not set, nItems is computed with std::distance(begin, end).
 * 
 * Example:
 * \code
 * #include <iostream>
 * #include <algorithm>
 * #include <vector>
 * #include <vigra/threadpool.hxx>
 * 
 * using namespace std;
 * using namespace vigra;
 * 
 * int main()
 * {
 *     size_t const n_threads = 4;
 *     size_t const n = 2000;
 *     vector<size_t> input(n);
 *     iota(input.begin(), input.end(), 0);
 *     
 *     // Compute the sum of the elements in the input vector.
 *     vector<size_t> results(n_threads, 0);
 *     parallel_foreach(n_threads, input.begin(), input.end(),
 *         [&results](size_t thread_id, size_t x)
 *         {
 *             results[thread_id] += x;
 *         }
 *     );
 *     size_t const sum = accumulate(results.begin(), results.end(), 0);
 *     
 *     cout << "The sum " << sum << " should be equal to " << (n*(n-1))/2 << endl;
 * }
 * \endcode
 */
template<class ITER, class F>
inline void parallel_foreach(
    int64_t nThreads,
    ITER begin, 
    ITER end, 
    F && f,                  
    const uint64_t nItems = 0
){
    nThreads = nThreads==-1 ? std::thread::hardware_concurrency() : nThreads;
    vigra_precondition(nThreads > 0, "parallel_foreach(): nThreads must be > 0 or -1.");
    
    ThreadPool pool(nThreads);
    parallel_foreach(pool, begin, end, f, nItems);
}



} // namespace vigra

#endif // VIGRA_THREADPOOL_HXX
