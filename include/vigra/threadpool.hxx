// Thread pool by Thorsten Beier.
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


namespace vigra
{



class ThreadPool {
public:

    ThreadPool(size_t);

    template<class F>
    std::future<typename std::result_of<F(int)>::type>  enqueueReturning(F&& f) ;

    template<class F>
    void enqueue(F&& f) ;

    ~ThreadPool();

    void waitFinished()
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        conditionF.wait(lock, [this](){ return tasks.empty() && (busy == 0); });
    }

    size_t nThreads()const{
        return workers.size();
    }
private:
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers;

    // the task queue
    std::queue< std::function<void(int)> > tasks;
    
    // synchronization
    std::mutex queue_mutex;
    std::condition_variable condition;
    std::condition_variable conditionF;

    bool stop;
    std::atomic_uint busy, processed;
};
 
// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads)
    :   stop(false),
        busy(ATOMIC_VAR_INIT(0U)),
        processed(ATOMIC_VAR_INIT(0U))
{
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
                        this->condition.wait(lock,[this]{ return this->stop || !this->tasks.empty(); });
                        if(!this->tasks.empty()){
                            ++busy;
                            task = std::move(this->tasks.front());
                            this->tasks.pop();
                            lock.unlock();
                            task(ti);
                            ++processed;
                            --busy;
                            conditionF.notify_one();
                        }
                        else if(stop){
                            return;
                        }
                    }
                    
                }
            }
        );
    }
}

// add new work item to the pool
template<class F>
std::future<typename std::result_of<F(int)>::type> 
ThreadPool::enqueueReturning(F&& f)
{
    typedef typename std::result_of<F(int)>::type result_type;
    typedef std::future<result_type>  FutureResType;

    auto task = std::make_shared< std::packaged_task<result_type(int)> >(
        f                                                           
    );

    FutureResType res = task->get_future();
    //std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace([task](int tid){ (*task)(tid); });
    }
    condition.notify_one();
    return res;
}

// add new work item to the pool
template<class F>
void ThreadPool::enqueue(F&& f)
{
    //using return_type = typename std::result_of<F()>::type;
    
    //std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace(f);
    }
    condition.notify_one();
    //return res;
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool()
{
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(std::thread &worker: workers)
        worker.join();
}



template<class ITER, class F>
void parallel_foreach_impl(
    ThreadPool & pool,
    const uint64_t nItems,                   
    ITER iter, 
    ITER end, 
    F && f,
    std::random_access_iterator_tag
){
    // typedef typename std::iterator_traits<ITER>::reference ReferenceType;
    uint64_t workload = std::distance(iter, end);
    const float workPerThread = float(workload)/pool.nThreads();
    const uint64_t chunkedWorkPerThread = std::max(uint64_t(std::round(workPerThread/3.0f)),1ul);

    for( ;iter<end; iter+=chunkedWorkPerThread){

        // ITER localEnd = iter+chunkedWorkPerThread < end ? iter + chunkedWorkPerThread : end;
        const size_t lc = std::min(workload, chunkedWorkPerThread);
        workload-=lc;
        pool.enqueue(
            [&f, iter, lc]
            (int id)
            {
                for(size_t i=0; i<lc; ++i)
                    f(id, iter[i]);
            }
        );
    }
    pool.waitFinished();
}



template<class ITER, class F>
void parallel_foreach_impl(
    ThreadPool & pool,
    const uint64_t nItems,                   
    ITER iter, 
    ITER end, 
    F && f,
    std::forward_iterator_tag
){

    // typedef typename std::iterator_traits<ITER>::reference ReferenceType;
    uint64_t workload = nItems;
    const float workPerThread = float(workload)/pool.nThreads();
    const uint64_t chunkedWorkPerThread = std::max(uint64_t(std::round(workPerThread/3.0f)),1ul);

    for(;;){

        const size_t lc = std::min(chunkedWorkPerThread, workload);
        workload -= lc;
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
        );
        if(workload==0)
            break;
        std::advance(iter, lc);
    }
    pool.waitFinished();
}



template<class ITER, class F>
void parallel_foreach_single_thread(
    const uint64_t nItems,
    ITER begin, 
    ITER end, 
    F && f
){
    for(size_t i=0; i<nItems; ++i){
        f(0, *begin);
        ++begin;
    }
}



template<class ITER, class F>
void parallel_foreach(
    ThreadPool & pool,
    const uint64_t nItems,
    ITER begin, 
    ITER end, 
    F && f
){

    if(pool.nThreads()!=1){
        parallel_foreach_impl(pool,nItems, begin, end, f,
            typename std::iterator_traits<ITER>::iterator_category());
    }
    else{
        parallel_foreach_single_thread(nItems, begin, end, f);
    }
}



template<class ITER, class F>
void parallel_foreach(
    const int64_t nThreads,                  
    const uint64_t nItems,
    ITER begin, 
    ITER end, 
    F && f
){
    ThreadPool pool(nThreads==-1 ?  std::thread::hardware_concurrency() : nThreads);

    if(pool.nThreads()!=1){
        parallel_foreach_impl(pool, nItems, begin, end, f,
            typename std::iterator_traits<ITER>::iterator_category());
    }
    else{
        parallel_foreach_single_thread(nItems, begin, end, f);
    }
}



} // namespace vigra

#endif // VIGRA_THREADPOOL_HXX
