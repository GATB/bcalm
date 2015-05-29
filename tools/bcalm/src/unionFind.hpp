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
    
    
    unsigned long getNumElements()
    {
		return data.size();
    }
    T getSet(T set)
    {
		auto it  = data.find (set);
		if (it == data.end()) {
			cerr << "Invalid reference to unionFind getSet()" << endl;
			exit(1);
		}
		return *it;
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

	T find(T key);
	void union_(T set1, T set2);
	//void compress (unsigned long start, unsigned long end);
	
	unsigned long getMutex(T key) {
		//TODO this function needs to be written for parallelism to work efficiently.
		//When T was an int type, this would return key%numLocks
		//But, in case T is not an int type, we need another way to get an arbitrary mutex.
		return 0;
	}
};

