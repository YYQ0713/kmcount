#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <vector>
#include <string>
#include <stdlib.h>

#include <omp.h>

#include "../include/common.hpp"
#include "../include/cxxopts.hpp"
#include "../include/kmc.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  // -------------------- //
  // Program name message
  // -------------------- //
  cxxopts::Options options("Kmcount", "large sequences set kmer counting");


  // -------------------- //
  //    Param inputs
  // -------------------- //
  options.add_options()
  ("k, kmer", "k-mer size",       cxxopts::value<uint32_t>()->default_value("21"))
  ("f, files", "List of Fasta/q(s) (required)", 	cxxopts::value<std::string>())
  ("o, output", "Output Filename (required)", 	cxxopts::value<std::string>())
  ("b, batchsize", "batchsize of file read once", 	cxxopts::value<uint64_t>()->default_value("1000000000"))
  ("h, help", "Usage")
  ;

  auto result = options.parse(argc, argv);

  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  char *inputfofn = NULL;	
  if(result.count("files")) inputfofn = strdup(result["files"].as<std::string>().c_str());
  else
  {
      std::cout << options.help() << std::endl;
      exit(0);		
  }

  char *OutputFile = NULL;	
  if(result.count("output"))
  {
    char* line1 = strdup(result["output"].as<std::string>().c_str());
    char* line2 = strdup(".out");

    unsigned int len1 = strlen(line1);
    unsigned int len2 = strlen(line2);

    OutputFile = (char*)malloc(len1 + len2 + 1);
    if (!OutputFile) abort();

    memcpy(OutputFile,  line1, len1);
    memcpy(OutputFile + len1,  line2, len2);
    OutputFile[len1 + len2] = '\0';

    delete line1, line2;
    
    remove(OutputFile);
  }
  else
  {
      std::cout << options.help() << std::endl;
      exit(0);		
  }

  Ctparams countpars;
	uint64_t CardinalityEstimate;

	countpars.kmerSize 	= result["kmer"].as<uint32_t>();
  countpars.batchsize = result["batchsize"].as<uint64_t>();

  vector<fileinfos> allfiles_path = GetFiles(inputfofn);

  CardinalityEstimate = hyperloglog(countpars.kmerSize, allfiles_path, countpars.batchsize);
  printf("test2 %ld\n", CardinalityEstimate);

}