#ifndef _COMMON_H_
#define _COMMON_H_


#ifndef PRINT
#define PRINT
#endif

#include <cstdlib>
#include <stdint.h>

#ifdef PRINT
	#define printLog(var) do { std::cerr << "INFO:	" << __FILE__ << "(" << __LINE__ << ")	" << #var << " = " << (var) << std::endl; } while(0)
#else
	#define printLog(var)
#endif

using namespace std;

#define MAX_FILE_PATH 1024

#ifdef _OPENMP
#define MYTHREAD omp_get_thread_num()
#define THREADS omp_get_num_threads()
#define MAXTHREADS omp_get_max_threads()
#else
#define MYTHREAD 0
#define THREADS 1
#define MAXTHREADS 1
#endif

struct Ctparams
{
	uint32_t kmerSize;			// KmerSize
    uint64_t batchsize;
	Ctparams(): kmerSize(21), batchsize(1000000000) {};
};











#endif
