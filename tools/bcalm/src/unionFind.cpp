#include <unionFind.hpp>
#include <iostream>
#include <cassert>
#include <mutex>
#include <thread>
#include <string.h>

using namespace std;

unionFind::unionFind(unsigned long n)
{
    numElements=n;
    sets=new unsigned long[numElements];
    for (unsigned long i=0; i<numElements; i++)
    {
        sets[i]=i;
    }
    numLocks=0;
    set_mutex=NULL;
}

unionFind::unionFind(unsigned long n, unsigned long nlocks)
{
    numElements=n;
    sets=new unsigned long[numElements];
    for (unsigned long i=0; i<numElements; i++)
    {
        sets[i]=i;
    }
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
unsigned long unionFind::find (unsigned long key)
{
    assert(key<numElements);

    unsigned long tmpKey=key;
    while (sets[tmpKey]!=tmpKey)
    {
//        if (key==997768)
//           printf ("tmpkey %ld parent %ld\n", tmpKey, sets[tmpKey]);
        tmpKey=sets[tmpKey];
    }
    unsigned long mutexIndex = 0;
    if (numLocks > 0)
    {
        mutexIndex = key%numLocks;
//        set_mutex[mutexIndex].lock();
    }
    sets[key]=tmpKey;
    if (numLocks > 0)
    {
//        set_mutex[mutexIndex].unlock();
    }
    return tmpKey;
}
void unionFind::union_ (unsigned long set1, unsigned long set2)
{
    unsigned long big, small;
    if (sets[set1] < sets[set2])
    {
        small=set1;
        big=set2;
    }
    else
    {
        small=set2;
        big=set1;
    }
//    cout << "small " << small << " big " << big << endl;
    unsigned long mutexIndex = 0;
    if (numLocks > 0)
    {
        mutexIndex = big%numLocks;
//        set_mutex[big].lock();
    }
    sets[big]=sets[small];
    if (numLocks > 0)
    {
//        set_mutex[big].unlock();
    }
}
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
