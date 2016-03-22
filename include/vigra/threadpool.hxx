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

#include <vector>
#include <queue>
#include <stdexcept>
#include <cmath>
#include "mathutil.hxx"
#include "counting_iterator.hxx"
#include "threading.hxx"


namespace vigra
{

/** \addtogroup ParallelProcessing
*/

//@{

    /**\brief Option base class for parallel algorithms.

        <b>\#include</b> \<vigra/threadpool.hxx\><br>
        Namespace: vigra
    */
class ParallelOptions
{
  public:

        /** Constants for special settings.
        */
    enum {
        Auto       = -1, ///< Determine number of threads automatically (from <tt>threading::thread::hardware_concurrency()</tt>)
        Nice       = -2, ///< Use half as many threads as <tt>Auto</tt> would.
        NoThreads  =  0  ///< Switch off multi-threading (i.e. execute tasks sequentially)
    };

    ParallelOptions()
    :   numThreads_(actualNumThreads(Auto))
    {}

        /** \brief Get desired number of threads.

            <b>Note:</b> This function may return 0, which means that multi-threading
            shall be switched off entirely. If an algorithm receives this value,
            it should revert to a sequential implementation. In contrast, if
            <tt>numThread() == 1</tt>, the parallel algorithm version shall be
            executed with a single thread.
        */
    int getNumThreads() const
    {
        return numThreads_;
    }

        /** \brief Get desired number of threads.

            In contrast to <tt>numThread()</tt>, this will always return a value <tt>>=1</tt>.
        */
    int getActualNumThreads() const
    {
        return std::max(1,numThreads_);
    }

        /** \brief Set the number of threads or one of the constants <tt>Auto</tt>,
                   <tt>Nice</tt> and <tt>NoThreads</tt>.

            Default: <tt>ParallelOptions::Auto</tt> (use system default)

            This setting is ignored if the preprocessor flag <tt>VIGRA_SINGLE_THREADED</tt>
            is defined. Then, the number of threads is set to 0 and all tasks revert to
            sequential algorithm implementations. The same can be achieved at runtime
            by passing <tt>n = 0</tt> to this function. In contrast, passing <tt>n = 1</tt>
            causes the parallel algorithm versions to be executed with a single thread.
            Both possibilities are mainly useful for debugging.
        */
    ParallelOptions & numThreads(const int n)
    {
        numThreads_ = actualNumThreads(n);
        return *this;
    }


  private:
        // helper function to compute the actual number of threads
    static size_t actualNumThreads(const int userNThreads)
    {
        #ifdef VIGRA_SINGLE_THREADED
            return 0;
        #else
            return userNThreads >= 0
                       ? userNThreads
                       : userNThreads == Nice
                               ? threading::thread::hardware_concurrency() / 2
                               : threading::thread::hardware_concurrency();
        #endif
    }

    int numThreads_;
};

/********************************************************/
/*                                                      */
/*                      ThreadPool                      */
/*                                                      */
/********************************************************/

    /**\brief Thread pool class to manage a set of parallel workers.

        <b>\#include</b> \<vigra/threadpool.hxx\><br>
        Namespace: vigra
    */
class ThreadPool
{
  public:

    /** Create a thread pool from ParallelOptions. The constructor just launches
        the desired number of workers. If the number of threads is zero,
        no workers are started, and all tasks will be executed in synchronously
        in the present thread.
     */
    ThreadPool(const ParallelOptions & options)
    :   stop(false)
    {
        init(options);
    }

    /** Create a thread pool with n threads. The constructor just launches
        the desired number of workers. If \arg n is <tt>ParallelOptions::Auto</tt>,
        the number of threads is determined by <tt>threading::thread::hardware_concurrency()</tt>.
        <tt>ParallelOptions::Nice</tt> will create half as many threads.
        If <tt>n = 0</tt>, no workers are started, and all tasks will be executed
        synchronously in the present thread. If the preprocessor flag
        <tt>VIGRA_SINGLE_THREADED</tt> is defined, the number of threads is always set
        to zero (i.e. synchronous execution), regardless of the value of \arg n. This
        is useful for debugging.
     */
    ThreadPool(const int n)
    :   stop(false)
    {
        init(ParallelOptions().numThreads(n));
    }

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
    auto enqueueReturning(F&& f) -> threading::future<decltype(f(0))>;

    /**
     * Enqueue function for tasks without return value.
     * This is a special case of the enqueueReturning template function, but
     * some compilers fail on <tt>std::result_of<F(int)>::type</tt> for void(int) functions.
     */
    template<class F>
    threading::future<void> enqueue(F&& f) ;

    /**
     * Block until all tasks are finished.
     */
    void waitFinished()
    {
        threading::unique_lock<threading::mutex> lock(queue_mutex);
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

    // helper function to init the thread pool
    void init(const ParallelOptions & options);

    // need to keep track of threads so we can join them
    std::vector<threading::thread> workers;

    // the task queue
    std::queue<std::function<void(int)> > tasks;

    // synchronization
    threading::mutex queue_mutex;
    threading::condition_variable worker_condition;
    threading::condition_variable finish_condition;
    bool stop;
    threading::atomic_long busy, processed;
};

inline void ThreadPool::init(const ParallelOptions & options)
{
    busy.store(0);
    processed.store(0);

    const size_t actualNThreads = options.getNumThreads();
    for(size_t ti = 0; ti<actualNThreads; ++ti)
    {
        workers.emplace_back(
            [ti,this]
            {
                for(;;)
                {
                    std::function<void(int)> task;
                    {
                        threading::unique_lock<threading::mutex> lock(this->queue_mutex);

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
        threading::unique_lock<threading::mutex> lock(queue_mutex);
        stop = true;
    }
    worker_condition.notify_all();
    for(threading::thread &worker: workers)
        worker.join();
}

template<class F>
inline auto
ThreadPool::enqueueReturning(F&& f) -> threading::future<decltype(f(0))>
{
    typedef decltype(f(0)) result_type;
    typedef threading::packaged_task<result_type(int)> PackageType;

    auto task = std::make_shared<PackageType>(f);
    auto res = task->get_future();

    if(workers.size()>0){
        {
            threading::unique_lock<threading::mutex> lock(queue_mutex);

            // don't allow enqueueing after stopping the pool
            if(stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");

            tasks.emplace(
                [task](int tid)
                {
                    (*task)(std::move(tid));
                }
            );
        }
        worker_condition.notify_one();
    }
    else{
        (*task)(0);
    }

    return res;
}

template<class F>
inline threading::future<void>
ThreadPool::enqueue(F&& f)
{
#if defined(USE_BOOST_THREAD) && \
    !defined(BOOST_THREAD_PROVIDES_VARIADIC_THREAD)
    // Without variadic templates, boost:thread::packaged_task only
    // supports the signature 'R()' (functions with no arguments).
    // We bind the thread_id parameter to 0, so this parameter
    // must NOT be used in function f (fortunately, this is the case
    // for the blockwise versions of convolution, labeling and
    // watersheds).
    typedef threading::packaged_task<void()> PackageType;
    auto task = std::make_shared<PackageType>(std::bind(f, 0));
#else
    typedef threading::packaged_task<void(int)> PackageType;
    auto task = std::make_shared<PackageType>(f);
#endif

    auto res = task->get_future();
    if(workers.size()>0){
        {
            threading::unique_lock<threading::mutex> lock(queue_mutex);

            // don't allow enqueueing after stopping the pool
            if(stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");

            tasks.emplace(
               [task](int tid)
               {
#if defined(USE_BOOST_THREAD) && \
    !defined(BOOST_THREAD_PROVIDES_VARIADIC_THREAD)
                    (*task)();
#else
                    (*task)(std::move(tid));
#endif
               }
            );
        }
        worker_condition.notify_one();
    }
    else{
#if defined(USE_BOOST_THREAD) && \
    !defined(BOOST_THREAD_PROVIDES_VARIADIC_THREAD)
        (*task)();
#else
        (*task)(0);
#endif
    }
    return res;
}

/********************************************************/
/*                                                      */
/*                   parallel_foreach                   */
/*                                                      */
/********************************************************/

// nItems must be either zero or std::distance(iter, end).
// NOTE: the redundancy of nItems and iter,end here is due to the fact that, for forward iterators,
// computing the distance from iterators is costly, and, for input iterators, we might not know in advance
// how many items there are  (e.g., stream iterators).
template<class ITER, class F>
inline void parallel_foreach_impl(
    ThreadPool & pool,
    const std::ptrdiff_t nItems,
    ITER iter,
    ITER end,
    F && f,
    std::random_access_iterator_tag
){
    std::ptrdiff_t workload = std::distance(iter, end);
    vigra_precondition(workload == nItems || nItems == 0, "parallel_foreach(): Mismatch between num items and begin/end.");
    const float workPerThread = float(workload)/pool.nThreads();
    const std::ptrdiff_t chunkedWorkPerThread = std::max<std::ptrdiff_t>(roundi(workPerThread/3.0), 1);

    std::vector<threading::future<void> > futures;
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
    {
        fut.get();
    }
}



// nItems must be either zero or std::distance(iter, end).
template<class ITER, class F>
inline void parallel_foreach_impl(
    ThreadPool & pool,
    const std::ptrdiff_t nItems,
    ITER iter,
    ITER end,
    F && f,
    std::forward_iterator_tag
){
    if (nItems == 0)
        nItems = std::distance(iter, end);

    std::ptrdiff_t workload = nItems;
    const float workPerThread = float(workload)/pool.nThreads();
    const std::ptrdiff_t chunkedWorkPerThread = std::max<std::ptrdiff_t>(roundi(workPerThread/3.0), 1);

    std::vector<threading::future<void> > futures;
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
    const std::ptrdiff_t nItems,
    ITER iter,
    ITER end,
    F && f,
    std::input_iterator_tag
){
    std::ptrdiff_t num_items = 0;
    std::vector<threading::future<void> > futures;
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

// Runs foreach on a single thread.
// Used for API compatibility when the numbe of threads is 0.
template<class ITER, class F>
inline void parallel_foreach_single_thread(
    ITER begin,
    ITER end,
    F && f,
    const std::ptrdiff_t nItems = 0
){
    std::ptrdiff_t n = 0;
    for (; begin != end; ++begin)
    {
        f(0, *begin);
        ++n;
    }
    vigra_postcondition(n == nItems || nItems == 0, "parallel_foreach(): Mismatch between num items and begin/end.");
}

/** \brief Apply a functor to all items in a range in parallel.

    <b> Declarations:</b>

    \code
    namespace vigra {
        // pass the desired number of threads or ParallelOptions::Auto
        // (creates an internal thread pool accordingly)
        template<class ITER, class F>
        void parallel_foreach(int64_t nThreads,
                              ITER begin, ITER end,
                              F && f,
                              const uint64_t nItems = 0);

        // use an existing thread pool
        template<class ITER, class F>
        void parallel_foreach(ThreadPool & pool,
                              ITER begin, ITER end,
                              F && f,
                              const uint64_t nItems = 0);

        // pass the integers from 0 ... (nItems-1) to the functor f,
        // using the given number of threads or ParallelOptions::Auto
        template<class F>
        void parallel_foreach(int64_t nThreads,
                              uint64_t nItems,
                              F && f);

        // likewise with an existing thread pool
        template<class F>
        void parallel_foreach(ThreadPool & threadpool,
                              uint64_t nItems,
                              F && f);
    }
    \endcode

    Create a thread pool (or use an existing one) to apply the functor \arg f
    to all items in the range <tt>[begin, end)</tt> in parallel. \arg f must
    be callable with two arguments of type <tt>size_t</tt> and <tt>T</tt>, where
    the first argument is the thread index (starting at 0) and T is convertible
    from the iterator's <tt>reference_type</tt> (i.e. the result of <tt>*begin</tt>).

    If the iterators are forward iterators (<tt>std::forward_iterator_tag</tt>), you
    can provide the optional argument <tt>nItems</tt> to avoid the a
    <tt>std::distance(begin, end)</tt> call to compute the range's length.

    Parameter <tt>nThreads</tt> controls the number of threads. <tt>parallel_foreach</tt>
    will split the work into about three times as many parallel tasks.
    If <tt>nThreads = ParallelOptions::Auto</tt>, the number of threads is set to
    the machine default (<tt>std::thread::hardware_concurrency()</tt>).

    If <tt>nThreads = 0</tt>, the function will not use threads,
    but will call the functor sequentially. This can also be enforced by setting the
    preprocessor flag <tt>VIGRA_SINGLE_THREADED</tt>, ignoring the value of
    <tt>nThreads</tt> (useful for debugging).

    <b>Usage:</b>

    \code
    #include <iostream>
    #include <algorithm>
    #include <vector>
    #include <vigra/threadpool.hxx>

    using namespace std;
    using namespace vigra;

    int main()
    {
        size_t const n_threads = 4;
        size_t const n = 2000;
        vector<int> input(n);

        auto iter = input.begin(),
             end  = input.end();

        // fill input with 0, 1, 2, ...
        iota(iter, end, 0);

        // compute the sum of the elements in the input vector.
        // (each thread computes the partial sum of the items it sees
        //  and stores the sum at the appropriate index of 'results')
        vector<int> results(n_threads, 0);
        parallel_foreach(n_threads, iter, end,
            // the functor to be executed, defined as a lambda function
            // (first argument: thread ID, second argument: result of *iter)
            [&results](size_t thread_id, int items)
            {
                results[thread_id] += items;
            }
        );

        // collect the partial sums of all threads
        int sum = accumulate(results.begin(), results.end(), 0);

        cout << "The sum " << sum << " should be equal to " << (n*(n-1))/2 << endl;
    }
    \endcode
 */
doxygen_overloaded_function(template <...> void parallel_foreach)

template<class ITER, class F>
inline void parallel_foreach(
    ThreadPool & pool,
    ITER begin,
    ITER end,
    F && f,
    const std::ptrdiff_t nItems = 0)
{
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

template<class ITER, class F>
inline void parallel_foreach(
    int64_t nThreads,
    ITER begin,
    ITER end,
    F && f,
    const std::ptrdiff_t nItems = 0)
{

    ThreadPool pool(nThreads);
    parallel_foreach(pool, begin, end, f, nItems);
}

template<class F>
inline void parallel_foreach(
    int64_t nThreads,
    std::ptrdiff_t nItems,
    F && f)
{
    auto iter = range(nItems);
    parallel_foreach(nThreads, iter, iter.end(), f, nItems);
}


template<class F>
inline void parallel_foreach(
    ThreadPool & threadpool,
    std::ptrdiff_t nItems,
    F && f)
{
    auto iter = range(nItems);
    parallel_foreach(threadpool, iter, iter.end(), f, nItems);
}

//@}

} // namespace vigra

#endif // VIGRA_THREADPOOL_HXX
