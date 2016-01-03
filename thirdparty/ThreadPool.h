// https://github.com/progschj/ThreadPool/blob/master/ThreadPool.h
//
// modified so that a thread_id integer in [0..nb_threads] is passed to each task
//
#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

class ThreadPool {
public:
    ThreadPool(size_t);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) 
        -> std::future<typename std::result_of<F(int, Args...)>::type>;
    //~ThreadPool();
    void join();
private:
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers;
    // the task queue
    std::queue< std::function<void(int)> > tasks;
    
    // synchronization
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};
 
// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads)
    :   stop(false)
{
    for(size_t thread_id = 0; thread_id<threads; ++ thread_id)
        workers.emplace_back(
            [this, thread_id]
            {
                for(;;)
                {
                    std::function<void(int)> task;

                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                            [this]{ return this->stop || !this->tasks.empty(); });
                        if(this->stop && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }

                    task(thread_id);
                }
            }
        );
}

// add new work item to the pool
template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args) 
    -> std::future<typename std::result_of<F(int, Args...)>::type>
{
    using return_type = typename std::result_of<F(int, Args...)>::type;

    auto task = std::make_shared< std::packaged_task<return_type(int)> >(
            std::bind(std::forward<F>(f), placeholders::_1, std::forward<Args>(args)...)
        );
        
    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace([task](int thread_id){ (*task)(thread_id); });
    }
    condition.notify_one();
    return res;
}

// the destructor joins all threads
// rayan: slightly modified, now explicit join
void ThreadPool::join()
{
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(std::thread &worker: workers)
        worker.join();
}

#endif
