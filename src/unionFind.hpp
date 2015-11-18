#include <iostream>
#include <cassert>
#include <mutex>
#include <thread>
#include <string.h>
#include <unordered_map>




using namespace std;


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
#ifdef SPARSEHASH
		data.set_deleted_key("");
#endif

    }

	unionFind(unsigned long nlocks)
	{
		numLocks = nlocks;
		if (nlocks > 0)
		{
			set_mutex = new mutex[numLocks];
		}
		else
		{
			set_mutex = NULL;
		}
	}

    ~unionFind()
    {
        cout << "destructor" << endl;
        if (set_mutex != NULL)
        {
//            delete set_mutex;
            set_mutex=NULL;
        }
        cout << "destructor end" << endl;
    }
    
    
    T getSet(T set)
    {
		auto it = data.find (set);
		if (it == data.end()) {
			cerr << "Invalid reference to unionFind getSet()" << endl;
			exit(1);
		}
		return it->second;
    }

    // O(1)
    unsigned long getNumKeys()
    {
		return data.size();
    }

    // O(|keys|)
    unsigned long getNumSets()
    {
        std::set<T> finalSets;
        for ( auto it = data.begin(); it != data.end(); ++it )
        {
		        finalSets.insert(it->second);
        }
        return finalSets.size();
    }

    // O(|keys|)
    void setStats(unsigned long &mean, unsigned long &max)
    {
        std::unordered_map<T, std::set<T>> reverseData;
        for ( auto it = data.begin(); it != data.end(); ++it )
        {
		        reverseData[it->second].insert(it->first);
        }

        mean = 0; max = 0;

        for ( auto it = reverseData.begin(); it != reverseData.end(); ++it )
        {
            max = std::max(it->second.size(), max);
            mean += it->second.size();
        }
        mean /= reverseData.size();
    }



    T find(T key){
        if (data.find(key) == data.end()) {
            data[key] = key;
        }

        T tmpKey = key;
        while (data[tmpKey] != tmpKey) {
            tmpKey = data[tmpKey];
        }


        unsigned long mutexIndex = 0;

        if (numLocks > 0)
        {
            //mutexIndex = getMutex(key);
            //set_mutex[mutexIndex].lock(); 
        }
        data[key] = tmpKey;

        if (numLocks > 0)
        {
            //        set_mutex[mutexIndex].unlock();
        }
        return tmpKey;
    }

	void union_(T set1, T set2){
        T big, small;
        if (find(set1) < find(set2)) {
            small = set1;
            big = set2;
        } else {
            small = set2;
            big = set1;
        }

        //    cout << "small " << small << " big " << big << endl;
        unsigned long mutexIndex = 0;
        if (numLocks > 0)
        {
            //mutexIndex = getMutex(big); // FIXME: doesn't look like this variable is used
            //        set_mutex[big].lock();
        }
        data[big] = data[small];
        if (numLocks > 0)
        {
            //        set_mutex[big].unlock();
        }
    }


private:
#ifdef SPARSEHASH
	typedef sparse_hash_map<T, T> ufDS;
#else
	typedef unordered_map<T, T> ufDS;
#endif

	ufDS data; 
	unsigned long numLocks;
	mutex *set_mutex;

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

    unsigned long getMutex(T key) {
		//TODO this function needs to be written for parallelism to work efficiently.
		//When T was an int type, this would return key%numLocks
		//But, in case T is not an int type, we need another way to get an arbitrary mutex.
		return 0;
	}
};


