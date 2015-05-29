#include <iostream>
#include <cassert>
#include <mutex>
#include <thread>
#include <string.h>




using namespace std;

class unionFind
{
private:
    unsigned long numElements;
    unsigned long *sets;
    unsigned long numLocks;
    mutex *set_mutex;
public:
    unionFind()
    {
        numElements=0;
        sets=NULL;
        numLocks=0;
        set_mutex=NULL;
    }

    ~unionFind()
    {
        cout << "destructor" << endl;
        if (sets != NULL)
        {
            delete sets;
            sets=NULL;
        }
        if (set_mutex != NULL)
        {
//            delete set_mutex;
            set_mutex=NULL;
        }
        cout << "destructor end" << endl;
    }
    
    unionFind(unsigned long n);
    unionFind(unsigned long n, unsigned long numLocks);
    
    unsigned long find(unsigned long);
    void union_(unsigned long set1, unsigned long set2);
    void compress (unsigned long start, unsigned long end);
    unsigned long getNumElements()
    {
        return numElements;
    }
    unsigned long getSet(unsigned long set)
    {
        assert(set < numElements);
        return sets[set];
    }
};

