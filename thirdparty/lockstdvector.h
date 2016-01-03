// just a thread safe queue
// actually using std::vector instead of std::queue to limit overheads
// see http://stackoverflow.com/questions/14784551/c-stl-queue-memory-usage-compared-to-vector

// adapted from:
// https://raw.githubusercontent.com/cameron314/concurrentqueue/master/benchmarks/stdqueue.h
// ©2014 Cameron Desrochers.

#pragma once

#include <vector>


// Simple wrapper around std::queue (not thread safe) - RC: made it thread safe
template<typename T>
class LockStdVector
{
    
public:

    LockStdVector()
    {
        nb_items = 0;
    }

    template<typename U>
    inline bool enqueue(U&& item)
    {
		std::lock_guard<std::mutex> guard(mutex);
        unsigned long size = v.size();
        v.push_back(std::forward<U>(item));
        nb_items ++;
        return true;
    }
    
    inline bool try_dequeue(T& item)
    {
		std::lock_guard<std::mutex> guard(mutex);
        if (nb_items == 0) {
            return false;
        }
        
        item = std::move(v.back());
        v.pop_back();
        nb_items--;
        return true;
    }
    
    unsigned long size_approx()
    {
        return v.size(); 
    }
	
    unsigned long overhead_per_element()
    {
        return 0 ; 
    }
private:
    std::vector<T> v;
	mutable std::mutex mutex;
    unsigned long nb_items;
};
