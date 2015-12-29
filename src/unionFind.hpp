#if !defined(__UNIONFIND_H)
#define __UNIONFIND_H

#include <vector>
#include <atomic>
#include <iostream>

#include <mutex>
#include <unordered_map>

/**
 * Lock-free parallel disjoint set data structure (aka UNION-FIND)
 * with path compression and union by rank
 *
 * Supports concurrent find(), same() and unite() calls as described
 * in the paper
 *
 * "Wait-free Parallel Algorithms for the Union-Find Problem"
 * by Richard J. Anderson and Heather Woll
 *
 * In addition, this class supports optimistic locking (try_lock/unlock)
 * of disjoint sets and a combined unite+unlock operation.
 *
 * \author Wenzel Jakob
 *
 * modified by Rayan Chikhi to use an unordered_map
 */
template<typename T>
class unionFind {
public:
    unionFind (uint32_t size) {
    }

    // still mutexes are needed because mData may resize itself
    std::mutex set_mutex;


    uint32_t find(uint32_t id) {
        while (id != parent(id)) {
            set_mutex.lock();
            if ( mData.find(id) == mData.end())
                mData[id] = id;
            uint64_t value = mData[id];
            set_mutex.unlock();
            uint32_t new_parent = parent((uint32_t) value);
            uint64_t new_value =
                (value & 0xFFFFFFFF00000000ULL) | new_parent;
            /* Try to update parent (may fail, that's ok) */
            set_mutex.lock();
            if (value != new_value)
                mData[id].compare_exchange_weak(value, new_value);
            set_mutex.unlock();
            id = new_parent;
        }
        return id;
    }

    bool same(uint32_t id1, uint32_t id2) const {
        for (;;) {
            id1 = find(id1);
            id2 = find(id2);
            if (id1 == id2)
                return true;
            if (parent(id1) == id1)
                return false;
        }
    }

    uint32_t union_(uint32_t id1, uint32_t id2) {
        for (;;) {
            id1 = find(id1);
            id2 = find(id2);

            if (id1 == id2)
                return id1;

            uint32_t r1 = rank(id1), r2 = rank(id2);

            if (r1 > r2 || (r1 == r2 && id1 < id2)) {
                std::swap(r1, r2);
                std::swap(id1, id2);
            }

            uint64_t oldEntry = ((uint64_t) r1 << 32) | id1;
            uint64_t newEntry = ((uint64_t) r1 << 32) | id2;

            set_mutex.lock();
            if (!mData[id1].compare_exchange_strong(oldEntry, newEntry))
            {
            set_mutex.unlock();
                continue;
            }
            set_mutex.unlock();

            if (r1 == r2) {
                oldEntry = ((uint64_t) r2 << 32) | id2;
                newEntry = ((uint64_t) (r2+1) << 32) | id2;
                /* Try to update the rank (may fail, that's ok) */
                set_mutex.lock();
                mData[id2].compare_exchange_weak(oldEntry, newEntry);
                set_mutex.unlock();
            }

            break;
        }
        return id2;
    }
    

    uint32_t size() const { return (uint32_t) mData.size(); }

    uint32_t rank(uint32_t id) {
        set_mutex.lock();
        uint32_t res = ((uint32_t) (mData[id] >> 32LL)) & 0x7FFFFFFFu;
        set_mutex.unlock();
        return res;
    }

    uint32_t parent(uint32_t id) {
        set_mutex.lock();
        if ( mData.find(id) == mData.end())
            mData[id] = id;
        uint32_t res = (uint32_t) mData[id];
        set_mutex.unlock();
        return res;
    }

    // compatibility with original unionFind.cpp
    uint32_t getSet(uint32_t key) { return find(key); }
    void printStats(std::string prefix) 
    {
        std::unordered_map<T, std::set<T>> reverseData;
        uint64_t getNumKeys, getNumSets;
        for ( auto it = mData.begin(); it != mData.end(); ++it )
            reverseData[find(it->second)].insert(it->first); // need to call find, because element isn't necessarily the representant of its equivalence class
        unsigned int mean = 0, max = 0;
        for ( auto it = reverseData.begin(); it != reverseData.end(); ++it )
        {
            max = std::max((unsigned int) it->second.size(), max);
            mean += it->second.size();
        }
        if (reverseData.size() > 0)
            mean /= reverseData.size();
        getNumSets = reverseData.size();
        getNumKeys = mData.size();
        std::cout << prefix + " data structure has " << getNumKeys << " inserted elements, and made " << getNumSets << " partitions." << std::endl;
        std::cout << "mean/max number of elements in partitions: " << mean << "/" << max << std::endl;
        std::cout << "raw space of UF hash data: " << ( 2*getNumKeys * sizeof(T)  ) /1024/1024 << " MB" << std::endl; // 2x because each key of type T is associated to a value of type T

    }


    mutable std::unordered_map<uint32_t, std::atomic<uint64_t>> mData;
};

#endif /* __UNIONFIND_H */



// old code

#if 0
/* seems to work ok in single thread, but not multithread
 *
 * there is still a bug, see e.g. fluctations in the number of partitions:
 *
 * uf kmers, std::string data structure has 2074559 inserted elements, and made 695780 partitions.
 * mean/max number of elements in partitions: 2/4817
 *
 * uf kmers data structure has 2074559 inserted elements, and made 695779 partitions.
 * mean/max number of elements in partitions: 2/4817
 *
 * such tricky parallelism.. I'll use another library for now.
 */

#include <iostream>
#include <cassert>
#include <mutex>
#include <vector>
#include <set>
#include <thread>
#include <string.h>
#include <unordered_map>


//#define SPARSEHASH
#ifdef SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;
#endif

template<class T> 
class unionFind {
public:
    unionFind()
    {
        numLocks=0;
        set_mutex=NULL;
        data.resize(1);
    }
    // for spasehash, one would think we need to call data.set_deleted_key("") 
    // but actually not. it's only if we used data.erase()

	unionFind(unsigned long nlocks)
	{
		numLocks = nlocks;
		if (nlocks > 0)
		{
			set_mutex = new std::mutex[numLocks];
		}
		else
		{
			set_mutex = NULL;
		}
        data.resize(numLocks);
	}

    ~unionFind()
    {
        //cout << "destructor" << endl;
        if (set_mutex != NULL)
        {
            delete[] set_mutex;
        }
        //cout << "destructor end" << endl;
    }
    
      
    T getSet(T set)
    {
		bool found = dataFind(set);
		if (!found) {
            std::cerr << "Invalid reference to unionFind getSet()" << std::endl;
			exit(1);
		}
		return find(set);
    }

    // O(1)
    unsigned long getNumKeys()
    {
        unsigned long res = 0;
        unsigned int bound  = (numLocks == 0) ? 1 : numLocks; 
        for (unsigned int i = 0; i < bound; i++)
            res += data[i].size();
        return res;
    }

    // O(|keys|)
    unsigned long getNumSets()
    {
        std::set<T> finalSets;
        unsigned int bound  = (numLocks == 0) ? 1 : numLocks; 
        for (unsigned int i = 0; i < bound; i++)
        for ( auto it = data[i].begin(); it != data[i].end(); ++it )
        {
		        finalSets.insert(find(it->second)); // need to call find, because element isn't necessarily the representant of its equivalence class
        }

        return finalSets.size();
    }

    // O(|keys|)
    void setStats(unsigned long &mean, unsigned long &max)
    {
        std::unordered_map<T, std::set<T>> reverseData;
        unsigned int bound  = (numLocks == 0) ? 1 : numLocks; 
        for (unsigned int i = 0; i < bound; i++)
        for ( auto it = data[i].begin(); it != data[i].end(); ++it )
        {
		        reverseData[find(it->second)].insert(it->first); // need to call find, because element isn't necessarily the representant of its equivalence class
        }

        mean = 0; max = 0;

        for ( auto it = reverseData.begin(); it != reverseData.end(); ++it )
        {
            max = std::max(it->second.size(), max);
            mean += it->second.size();
        }
        if (reverseData.size() > 0)
            mean /= reverseData.size();
    }

    void printStats(std::string prefix)
    {
        std::cout << prefix + " data structure has " << getNumKeys() << " inserted elements, and made " << getNumSets() << " partitions." << std::endl;
        unsigned long mean, max;
        setStats(mean, max);
        std::cout << "mean/max number of elements in partitions: " << mean << "/" << max << std::endl;
        std::cout << "raw space of UF hash data: " << ( 2*getNumKeys() * sizeof(T)  ) /1024/1024 << " MB" << std::endl; // 2x because each key of type T is associated to a value of type T
        //std::cout << "number of UF hash buckets: " << data.bucket_count() << std::endl;
        //std::cout << "estimate of UF hash space: " << ( 2*data.bucket_count() * sizeof(T)  ) /1024/1024 << " MB" << std::endl; // 2x because each key of type T is associated to a value of type T
        std::cout << "but don't forget there are additional data strutures constructed in setStats(), they're quite expensive" << std::endl; // 2x because each key of type T is associated to a value of type T
        
    }


    T dataRead(T key) // low level function for reading from the UF hash
                    // uses mutexes for maximum safery
    {
        T res;

        unsigned long mutexIndex = 0;
        if (numLocks > 0)
        {
            mutexIndex = getMutex(key);
            set_mutex[mutexIndex].lock(); 
        }

        res = data[mutexIndex][key];
        
        if (numLocks > 0)
            set_mutex[mutexIndex].unlock();
        
        return res;
    }
 
    void dataWrite(T key, T val) // low level function for reading from the UF hash
                    // uses mutexes for maximum safery
    {
        unsigned long mutexIndex = 0;
        if (numLocks > 0)
        {
            mutexIndex = getMutex(key);
            set_mutex[mutexIndex].lock(); 
        }

        data[mutexIndex][key] = val;

        if (numLocks > 0)
            set_mutex[mutexIndex].unlock();
        
    }
    
    bool dataFind(T key) // low level function for reading from the UF hash
                    // uses mutexes for maximum safery
    {
        unsigned long mutexIndex = 0;
        if (numLocks > 0)
        {
            mutexIndex = getMutex(key);
            set_mutex[mutexIndex].lock(); 
        }


        bool res = (data[mutexIndex].find(key) != data[mutexIndex].end());
 
        if (numLocks > 0)
            set_mutex[mutexIndex].unlock();
     
        return res;
    }


    T find(T key){
        if (! dataFind(key)) {
            dataWrite(key,key);
        }

        //vector<T> chain;
        
        T tmpKey = key;
        while (dataRead(tmpKey) != tmpKey) {
            //chain.push_back(tmpKey);
            tmpKey = dataRead(tmpKey);
        }

        dataWrite(key,tmpKey);

        // do "compress" on the fly
        // - disabled because I didn't observe any speedup
        //for (auto it = chain.begin(); it != chain.end(); it++)
        //    dataWrite(*it, tmpKey);


        return tmpKey;
    }

	void union_(T set1, T set2){

        T big,  small;
        T findSet1 = find(set1);
        T findSet2 = find(set2);
        if (findSet1 < findSet2) {
            small = set1;
            big = set2;
        } else {
            small = set2;
            big = set1;
        }

        //    cout << "small " << small << " big " << big << endl;

        dataWrite(big, small);

        assert(getSet(big) == getSet(small));
    }


private:
#ifdef SPARSEHASH
	typedef sparse_hash_map<T, T> ufDS;
#else
	typedef std::unordered_map<T, std::atomic<T>> ufDS;
#endif

	std::vector<ufDS> data; 
	
    unsigned long numLocks;
    std::mutex *set_mutex;

    //The following function was implemented but not used in the array implementation.
    //I (PM) did not port it because it is not used and its not immediately clear how to port it to a hash table implementation.
    /*
       void unionFind::compress (unsigned long start, unsigned long end)
       {
       unsigned long key=0,key1=0;
       for (unsigned long i=start; i<end; i++)
       {
       key=i;
       while (key!=sets[key])
       {
       key=sets[key];
       }
       key1=i;
       while (key!=sets[key1])
       {
       unsigned long tmpKey=sets[key1];
       if ((key1>=start) && (key1<end))
       {
       sets[key1]=key;
       }
       key1=tmpKey;
       }
       }
       }
       */

    //unsigned long getMutex(T key) {
		//TODO this function needs to be written for parallelism to work efficiently.
		//When T was an int type, this would return key%numLocks
		//But, in case T is not an int type, we need another way to get an arbitrary mutex.
	//	return 0;
	//}

    // good idea. I got you paul, here it is for hashed int's
    unsigned long getMutex(uint64_t key /*in bcalm, key is a hashed int, so no need to hash it further*/) {
		return key%numLocks;
	}

    // dummy for string, no parallelization
//    unsigned long getMutex(std::string key ) {
//        return 0;
//    }
};

#endif
