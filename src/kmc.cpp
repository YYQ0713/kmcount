#include "../include/common.hpp"
#include "../include/kmc.hpp"
#include "../include/bseq.hpp"

/**
 * @brief GetFiles
 * @param filename
 * @return
 */
vector<fileinfos>  GetFiles(char *filename) {
	int64_t totalsize = 0;
	int numfiles = 0;
	std::vector<fileinfos> filesview;
	
	fileinfos fdata;
	ifstream allfiles(filename);
	if(!allfiles.is_open()) {
		std::string someString(filename);
		std::string ErrorMessage = "Could not open " + someString;
		printLog(ErrorMessage);
		exit(1);
	}
	allfiles.getline(fdata.filename, MAX_FILE_PATH);
	while(!allfiles.eof())
	{
		struct stat st;
		stat(fdata.filename, &st);
		fdata.filesize = st.st_size;
		
		filesview.push_back(fdata);
		std::string InputFile = filesview.back().filename;
		std::string str1 = std::to_string((float)filesview.back().filesize / (1024*1024));
		std::string str2 = " MB";
		std::string InputSize = str1 + str2;

		printLog(InputFile);
		printLog(InputSize);
		allfiles.getline(fdata.filename, MAX_FILE_PATH);
		//totalsize += fdata.filesize;
		//numfiles++;
	}
	return filesview;
}

std::string lexsmaller(const char* seq, int i, int j) {
    std::string origin;
    std::string revComp;

    int length = j - i + 1;
    revComp.reserve(length);
    origin.reserve(length);

    for (int k = i; k <= j; k++) {
        char base = seq[k];
        origin.push_back(base);
    }

    for (int k = j; k >= i; --k) {
        char base = seq[k];
        switch (base) {
            case 'A': revComp.push_back('T'); break;
            case 'T': revComp.push_back('A'); break;
            case 'C': revComp.push_back('G'); break;
            case 'G': revComp.push_back('C'); break;
            default: revComp.push_back(base); // 不是有效的碱基，直接添加
        }
    }


    return revComp < origin ? revComp : origin;
}

uint64_t hyperloglog(uint32_t kmersize, vector<fileinfos> files, uint64_t chunk_size) {

    vector <HyperLogLog> hlls(MAXTHREADS, HyperLogLog(12));
	bseq_file_t *fp = 0;
    seqs_t seqs;
    double hlltime = omp_get_wtime();
    printf("test in hll 1\n");

    for(auto itf = files.begin(); itf != files.end(); itf++) {
        fp = bseq_open(itf->filename);
        printf("test in hll 2\n");

        while(1) {
            seqs.seq = bseq_read(fp, chunk_size, &seqs.n_seq);

            if (seqs.n_seq == 0)
            {
                break;
            }
            
            printf("test in hll 3 and %d\n", seqs.n_seq);
            #pragma omp parallel for
            for (auto i = 0; i < seqs.n_seq; i++)
            {
                int next = 0; // kmer's last base
                int prev = 0; // kmer's first base
                int len = seqs.seq[i].l_seq;
                char *s = seqs.seq[i].seq;
                bool last_valid = false;
                bool lex = false;

                if (len < kmersize)
                {
                    continue;
                }
                
                string kmer;
                while(next < len) {
                    char c = s[next];
                    if (c != 'N' && c != 'n') {
                        if (last_valid) {
                            kmer = lexsmaller(s, prev, next);
                            hlls[MYTHREAD].add(kmer.c_str(), kmersize);
                            prev++;
                            next++;
                        } else {
                            if (prev + kmersize -1 == next) {
                                last_valid = true;
                                //next++;
                            } else {
                                next++;
                            }
                        }
                    } else {
                        // invalid character, restart
                        next++;
                        prev = next;
                        last_valid = false;
                    }



                } // while(next < len)
                
            } // for loop n_seqs
            printf("test in hll 3 and %d 2\n", seqs.n_seq);
            free(seqs.seq->seq);
            free(seqs.seq->name);
            free(seqs.seq);
            printf("test in hll 3 and %d 3\n", seqs.n_seq);
            seqs.n_seq = 0;
        } // while()




    } // all files
    printf("test in hll 4\n");
    // HLL reduction (serial for now) to avoid double iteration
    for (int i = 1; i < MAXTHREADS; i++) 
    {
        std::transform(hlls[0].M.begin(), hlls[0].M.end(), hlls[i].M.begin(), hlls[0].M.begin(), [](uint8_t c1, uint8_t c2) -> uint8_t{ return std::max(c1, c2); });
    }

    uint64_t CardinalityEstimate = hlls[0].estimate();

    return CardinalityEstimate;
    printLog(CardinalityEstimate);
}