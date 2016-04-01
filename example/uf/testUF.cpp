
#include <random>
#include <set>
#include "../../bglue/unionFind.hpp"
    
int main()
{
    
    unionFind<uint32_t> uf(1000); 
    unionFind<uint32_t> uf2(1000); 

    // generate 1000 random integers 
    int nelem=10;
    u_int64_t *data;
	static std::mt19937_64 rng;
	//rng.seed(std::mt19937_64::default_seed);
	rng.seed(static_cast<unsigned long>(time(NULL)));
	data = (u_int64_t * ) calloc(nelem,sizeof(u_int64_t));
	for (u_int64_t i = 0; i < nelem; i++)
		data[i] = rng() % nelem;
    
    #pragma omp parallel for    
	for (u_int64_t i = 0; i <= nelem-2; i+=2)
    {
        uf.union_(data[i],data[i+1]);
#pragma omp critical
        if (nelem <= 10)
        {
            std::cout<< "uf1, union " << data[i] << " " << data[i+1] << std::endl;
        }
    }
    std::cout<<std::endl;
	for (signed int i = nelem-2; i >= 0; i-=2)
    {
        uf2.union_(data[i+1],data[i]);
        if (nelem <= 10)
            std::cout<< "uf2, union " << data[i+1] << " " << data[i] << std::endl;
    }
    uf.printStats("uf1");
    uf2.printStats("uf2");

    // if the UF is small, display it
    if (nelem <= 10)
    {
        std::set<int> already_displayed;
        for (u_int64_t i = 0; i < nelem; i++)
        {
            if (already_displayed.find(data[i]) == already_displayed.end())
                std::cout << "uf element " << data[i] << " has partition " << uf.getSet(data[i]) << std::endl;
            already_displayed.insert(data[i]);
        }
    }

    // check correctness
    for (u_int64_t i = 0; i <= nelem-2; i+=2)
    {
        if (uf.getSet(data[i]) != uf.getSet(data[i+1]))
        {
            std::cout << "inconsistent partitioning, elements " << data[i] << " and " << data[i+1] << " should've been joined" << std::endl;
            exit(1);
        }
    }
    return 0;
}
