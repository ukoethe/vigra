#ifndef VIGRA_THREAD_POOL
#define VIGRA_THREAD_POOL


#include <queue>
#include <boost/bind.hpp>
#include <boost/thread.hpp>

namespace vigra{


class ThreadPool
{
private:
  std::queue< boost::function< void() > > tasks_;
  boost::thread_group threads_;
  std::size_t available_;
  boost::mutex mutex_;
  boost::condition_variable condition_;
  bool running_;
public:

  /// @brief Constructor.
  ThreadPool( std::size_t pool_size )
    : available_( pool_size ),
      running_( true )
  {
    for ( std::size_t i = 0; i < pool_size; ++i )
    {
      threads_.create_thread( boost::bind( &ThreadPool::pool_main, this ) ) ;
    }
  }

  /// @brief Destructor.
  ~ThreadPool()
  {
    // Set running flag to false then notify all threads.
    {
      boost::unique_lock< boost::mutex > lock( mutex_ );
      running_ = false;
      condition_.notify_all();
    }

    try
    {
      threads_.join_all();
    }
    // Suppress all exceptions.
    catch ( ... ) {}
  }

  /// @brief Add task to the thread pool if a thread is currently available.
  template < typename Task >
  void run_task( Task task )
  {
    boost::unique_lock< boost::mutex > lock( mutex_ );

    // If no threads are available, then return.
    if ( 0 == available_ ) return;

    // Decrement count, indicating thread is no longer available.
    --available_;

    // Set task and signal condition variable so that a worker thread will
    // wake up andl use the task.
    tasks_.push( boost::function< void() >( task ) );
    condition_.notify_one();
  }

private:
  /// @brief Entry point for pool threads.
  void pool_main()
  {
    while( running_ )
    {
      // Wait on condition variable while the task is empty and the pool is
      // still running.
      boost::unique_lock< boost::mutex > lock( mutex_ );
      while ( tasks_.empty() && running_ )
      {
        condition_.wait( lock );
      }
      // If pool is no longer running, break out.
      if ( !running_ ) break;

      // Copy task locally and remove from the queue.  This is done within
      // its own scope so that the task object is destructed immediately
      // after running the task.  This is useful in the event that the
      // function contains shared_ptr arguments bound via bind.
      {
        boost::function< void() > task = tasks_.front();
        tasks_.pop();

        lock.unlock();

        // Run the task.
        try
        {
          task();
        }
        // Suppress all exceptions.
        catch ( ... ) {}
      }

      // Task has finished, so increment count of available threads.
      lock.lock();
      ++available_;
    } // while running_
  }
};


} // end namespace vigra


#endif