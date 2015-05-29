#include <unionFind.hpp>
#include <iostream>
#include <cassert>
#include <mutex>
#include <thread>
#include <string.h>

using namespace std;

template<class T> 
T unionFind<T>::find (T key) {
	if (data.find(key) == data.end()) {
		data[key] = key;
	}

    T tmpKey = key;
	while (data[tmpKey] == key) {
		tmpKey = data[tmpKey];
	}


    unsigned long mutexIndex = 0;

    if (numLocks > 0)
    {
		mutexIndex = getMutex(key);
		//set_mutex[mutexIndex].lock();
    }
    data[key] = tmpKey;

    if (numLocks > 0)
    {
//        set_mutex[mutexIndex].unlock();
    }
    return tmpKey;
}

template<class T> 
void unionFind<T>::union_ (T set1, T set2) {
	T big, small;
	if (getSet(set1) < getSet(set2)) {
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
        mutexIndex = getMutex(big);
//        set_mutex[big].lock();
    }
    data[big] = data[small];
    if (numLocks > 0)
    {
//        set_mutex[big].unlock();
    }
}

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
