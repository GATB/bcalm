// just a thread safe queue, the most simple ever

// adapted from:
// https://raw.githubusercontent.com/cameron314/concurrentqueue/master/benchmarks/stdqueue.h
// ©2014 Cameron Desrochers.

#pragma once

#include <queue>


// Simple wrapper around std::queue (not thread safe) - RC: made it thread safe
template<typename T>
class LockStdQueue 
{
    
public:
    template<typename U>
    inline bool enqueue(U&& item)
    {
		std::lock_guard<std::mutex> guard(mutex);
        q.push(std::forward<U>(item));
        return true;
    }
    
    inline bool try_dequeue(T& item)
    {
		std::lock_guard<std::mutex> guard(mutex);
        if (q.empty()) {
            return false;
        }
        
        item = std::move(q.front());
        q.pop();
        return true;
    }
    
    unsigned long size_approx()
    {
        return q.size(); 
    }
	
    unsigned long overhead_per_element()
    {
        return 0; // I don't think anymore that's true. there must be some overhead
    }

private:
    std::queue<T> q;
	mutable std::mutex mutex;
};
