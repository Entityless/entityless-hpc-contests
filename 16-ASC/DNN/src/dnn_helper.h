#ifndef DNN_HELPER_H
#define DNN_HELPER_H

#include <sys/time.h>
#include <stdio.h>
#include <math.h>

#define MAXLINE 1024
#define MAXCHUNK 102400

#include "dnn_utility.h"

#define DEBUG_YGY 1 

struct WorkPara
{
	int fea_dim;
	int fea_context;
	int targ_offset;
	
	char train_sent_range[MAXLINE];
	char cv_sent_range[MAXLINE];
	
	int traincache;
	int bunchsize;
	
	int init_randem_seed;
	float init_randem_weight_min;
	float init_randem_weight_max;
	float init_randem_bias_max;
	float init_randem_bias_min;
	
	float* indata;
	int* targ;
};

class Interface
{
public:
		Interface();
		~Interface();
public:
		int Initial(CpuArg& cpuArg);
		void get_pfile_info(CpuArg& cpuArg);
		void get_chunk_info(CpuArg& cpuArg);
		void get_chunk_info_cv(CpuArg& cpuArg);
		int Readchunk(int index, CpuArg& cpuArg, ChunkContainer& oneChunk);
        int ReadchunkPara(int index, CpuArg& cpuArg, ChunkContainer& oneChunk, int s_sample, int e_sample);
		int Readchunk_cv(int index, CpuArg& cpuArg, ChunkContainer& oneChunk);
		void GetRandIndex(int *vec, int len);
private:
		struct WorkPara *para;
		
		unsigned int total_frames;
		unsigned int total_sents;
		unsigned int total_chunks;
		unsigned int total_samples;
		unsigned int cv_total_chunks;
		unsigned int cv_total_samples;
		
		int *framesBeforeSent;
		int *chunk_frame_st;
		int *cv_chunk_frame_st;
		
		FILE *fp_log;
		int numlayers;
		int realbunchsize;
		void get_uint(const char* hdr, const char* argname, unsigned int* val);
		void read_tail(FILE *fp, long int file_offset, unsigned int sentnum, int *out);
		
		void GetRandWeight(float *vec, float min, float max, int len);
		
		FILE *fp_data;
		FILE *fp_targ;
		FILE *fp_init_weight;
		FILE *fp_norm;
		FILE *fp_out;
		
		int data_rand_index[MAXLINE];
		
		float *mean;
		float *dVar;

#if DEBUG_YGY
		char *checkFileName;
		FILE *pCheckFile;
		int checkCount;
#endif		
		int sent_st, sent_en;
		int cv_sent_st, cv_sent_en;
		int cur_chunk_index;
		int frames_read;
};

#endif
