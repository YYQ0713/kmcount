#ifndef _KMERCOUNT_H
#define _KMERCOUNT_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include <omp.h>
#include <algorithm>


#include "common.hpp"
#include "hyperloglog.hpp"


struct fileinfos {
	char filename[MAX_FILE_PATH];
	size_t filesize;
};

vector<fileinfos>  GetFiles(char *filename);

uint64_t hyperloglog(uint32_t kmersize, vector<fileinfos> files, uint64_t chunk_size, uint32_t ts);


void kmer_counting(uint32_t kmersize, dictionary_t_16bit& count_kmer, vector<fileinfos> files, uint64_t CardinalityEstimate, uint64_t chunk_size);


std::string lexsmaller(const char* seq, int i, int j);


#endif