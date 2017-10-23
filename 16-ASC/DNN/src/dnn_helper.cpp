#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <assert.h>
#include <time.h>
#include <iostream>
using namespace std;

#define NOTDEBUG 1
#include "dnn_helper.h"


static void swap32(int *val)
{
#if BIGLITTLESWAP
    unsigned int uval;
    unsigned int res;
    int b0, b1, b2, b3;

    uval = (unsigned int) (*val);
    b0 = uval >> 24;
    b1 = (uval >> 8) & 0x0000ff00;
    b2 = (uval << 8) & 0x00ff0000;
    b3 = uval << 24;

    res = b0|b1|b2|b3;

    *val = (int) res;
#endif
    return;
}

Interface::Interface()
{
#if DEBUG_YGY
    checkFileName = (char*)"check.txt";
    pCheckFile = fopen(checkFileName,"w");
    checkCount = 0;
#endif
    para = new struct WorkPara;
    chunk_frame_st = new int[MAXCHUNK];
    cv_chunk_frame_st = new int[MAXCHUNK];
}

Interface::~Interface()
{
    fp_log = NULL;
    fp_data = NULL;
    fp_targ = NULL;
    fp_out = NULL;
    mean = NULL;
    dVar = NULL;

#if DEBUG_YGY
    fclose(pCheckFile);
    pCheckFile = NULL;
#endif
    para->indata = NULL;
    para->targ = NULL;
    delete para;
    delete []chunk_frame_st;
    delete []cv_chunk_frame_st;
    delete []framesBeforeSent;
}

int Interface::Initial(CpuArg& cpuArg)
{
    int i;
    char buff[MAXLINE];
    
    ////Set varibles default if neccessary
    para->init_randem_weight_min = -0.1;
    para->init_randem_weight_max = 0.1;
    para->init_randem_bias_min = -0.1;
    para->init_randem_bias_max = 0.1;
    for(i =0; i< MAXLINE; i++)
    {
        data_rand_index[i] = i;
    }
    
    if(NULL == cpuArg.pLogFile)
    {
        printf("FAIL to open log file:%s\n",cpuArg.logFileName);
        return -1;
    }
    // assign FILE
    fp_log = cpuArg.pLogFile;
    fp_data = cpuArg.pDataFile;
    fp_targ = cpuArg.pLabelFile;
    fp_out = cpuArg.pWtsFile;

    fprintf(fp_log,"parameters input:\n");
    fprintf(fp_log,"fea_file:             %s\n", cpuArg.dataFileName);
    fprintf(fp_log,"norm_file:            %s\n", cpuArg.normFileName);
    fprintf(fp_log,"targ_file:            %s\n", cpuArg.labelFileName);
    fprintf(fp_log,"outwts_file:          %s\n", cpuArg.outputWtsPath);
    fprintf(fp_log,"log_file:             %s\n", cpuArg.logFileName);
    fprintf(fp_log,"initwts_file:         %s\n", cpuArg.weightFileName);
    fprintf(fp_log,"train_sent_range:     %d-%d\n", cpuArg.train_bp_range_beg, cpuArg.train_bp_range_end);
    fprintf(fp_log,"cv_sent_range:        %d-%d\n", cpuArg.train_cv_range_beg, cpuArg.train_cv_range_end);
    fprintf(fp_log,"fea_dim:          %d\n", cpuArg.featDim);
    para->fea_dim = cpuArg.featDim;
    fprintf(fp_log,"fea_context:          %d\n", cpuArg.featCRange);
    para->fea_context = cpuArg.featCRange;
    fprintf(fp_log,"bunchsize:          %d\n", cpuArg.bunchSize);
    para->bunchsize = cpuArg.bunchSize;

      fprintf(fp_log,"train_cache:          %d\n", cpuArg.chunkSize);
    para->traincache = cpuArg.chunkSize;
      fprintf(fp_log,"init_randem_seed:     %d\n", cpuArg.randomSeed);
    para->targ_offset = (cpuArg.featCRange -1)/2;
      fprintf(fp_log,"targ_offset:          %d\n", para->targ_offset);
  
      fprintf(fp_log,"init_randem_weight_max:          %f\n", para->init_randem_weight_max);
      fprintf(fp_log,"init_randem_weight_min:          %f\n", para->init_randem_weight_min);
      fprintf(fp_log,"init_randem_bias_max:          %f\n", para->init_randem_bias_max);
      fprintf(fp_log,"init_randem_bias_min:          %f\n", para->init_randem_bias_min);
    fprintf(fp_log,"learnrate:                  %f\n", cpuArg.lRate);
    fprintf(fp_log,"layersizes:                  ");
    for(int j =0; j < cpuArg.dnnLayerNum; j++)
        fprintf(fp_log,"%d,", cpuArg.dnnLayerArr[j]);
    fprintf(fp_log,"\n");
        
    /******add for jumpframe and discard  ******/
    fprintf(fp_log,"nmod:                          %d\n", cpuArg.nmod);
    fprintf(fp_log,"remainder:                  %d\n", cpuArg.remainder);
    fprintf(fp_log,"discard_prob:                  %.6f\n", cpuArg.discard_prob);
    fprintf(fp_log,"discardLabs:                  ");
    for(int j =0; j < cpuArg.discard_num; j++)
        fprintf(fp_log,"%d,", cpuArg.discard_labs[j]);
    fprintf(fp_log, "\n");
    /************/

    fprintf(fp_log,"Please check...\n");
    para->init_randem_seed = cpuArg.randomSeed;    
    srand48(para->init_randem_seed); 
    mean = cpuArg.normInfo->meanArr;
    dVar = cpuArg.normInfo->dVarArr;
    //// Alloc data and target memory
    if(cpuArg.featDim * cpuArg.featCRange != cpuArg.dnnLayerArr[0])
    {
        fprintf(fp_log,"feadim times context must be equal to layersizes[0]\n");
        return -1;
    }
    
    fflush(fp_log);    
    return 0;
}

///Get frames number,sentence number and frame number per sent in training data and assert them in data and target pfile is consistent
void Interface::get_pfile_info(CpuArg& cpuArg)
{
    char *header = new char[PFILE_HEADER_SIZE]; // Store header here
    long int offset;
    int sizePerFrame;
    int *tmpframepersent;
    unsigned int tmpsentnum;
    unsigned int tmpframenum;
    
      /////Read data pfile
      fseek(fp_data,0,0);
      if (fread(header, PFILE_HEADER_SIZE, 1, fp_data) != 1)
      {
        fprintf(fp_log, "Failed to read data pfile header.\n");
          exit(0);
    }
    get_uint(header, "-num_sentences", &cpuArg.totalSentNum); 
    ///get sent num
    get_uint(header, "-num_frames", &cpuArg.totalFtrNum);
    ///get frame num
    total_frames = cpuArg.totalFtrNum;
    total_sents = cpuArg.totalSentNum;
    framesBeforeSent = new int[cpuArg.totalSentNum];
    sizePerFrame = sizeof(float)*(2+ cpuArg.featDim);
    offset = total_frames * (long int)sizePerFrame + PFILE_HEADER_SIZE;
    read_tail(fp_data, offset, total_sents, framesBeforeSent);  
    /// get frame per sent
    
#if NOTDEBUG
      ///Read target pfile
      fseek(fp_targ,0,0);
      if (fread(header, PFILE_HEADER_SIZE, 1, fp_targ) != 1)
    {
        fprintf(fp_log, "Failed to read target pfile header.\n");
          exit(0);
    }
      get_uint(header, "-num_sentences", &tmpsentnum); ///get sent num
    get_uint(header, "-num_frames", &tmpframenum);   ///get frame num
    
      tmpframepersent = new int[tmpsentnum];
    sizePerFrame = sizeof(int)*(2+ 1);
    offset = tmpframenum *(long int) sizePerFrame + PFILE_HEADER_SIZE;
    read_tail(fp_targ, offset, tmpsentnum, tmpframepersent);  /// get frame per sent
  
      ///assert consistency
      if(tmpsentnum != total_sents || tmpframenum != total_frames)
      {
          fprintf(fp_log, "frames or sentence num in target pfile and data pfile is not consistent.\n");
          exit(0);
      }
      for(int i=0; i<total_sents;i++ )
      {
          if(tmpframepersent[i] != framesBeforeSent[i] )
          {
              fprintf(fp_log, "tails in target pfile and data pfile is not consistent.\n");
              exit(0);
          }
      }
    delete []tmpframepersent;
#endif

  delete []header;
  fprintf(fp_log, "Get pfile info over: Training data has %u frames, %u sentences.\n", total_frames, total_sents);
}

/// get chunk nums and get frame start id for every chunk
void Interface::get_chunk_info(CpuArg& cpuArg)
{
    char *p, *st, *en;
    int sentid;
    int count_chunk = 1;
    int cur_frames_num = 0;
    int cur_frames_lost =0;
    int cur_frame_id = 0;
    int cur_chunk_frames =0;
    int frames_inc;
    int next_st;

    sent_st = cpuArg.train_bp_range_beg;
    sent_en = cpuArg.train_bp_range_end;
    if(sent_en < sent_st || sent_st < 0 || sent_en >= total_sents){
        fprintf(fp_log,"sent range: %d to %d number error.\n",sent_st,sent_en);
        exit(0);
    }
    
    if(sent_st ==0){
        cur_frame_id = 0;
    }
    else{
        cur_frame_id = framesBeforeSent[sent_st -1];
    }
    chunk_frame_st[0] = cur_frame_id;
    
    for(sentid =sent_st; sentid <= sent_en; sentid++)
    {
        frames_inc = framesBeforeSent[sentid] - cur_frame_id;
        cur_frame_id = framesBeforeSent[sentid];
        if(frames_inc >= cpuArg.featCRange){
            cur_frames_lost = cpuArg.featCRange -1;
        }
        else{
            cur_frames_lost = frames_inc;
        }
        cur_frames_num = frames_inc - cur_frames_lost;
        cur_chunk_frames += cur_frames_num;
        while(cur_chunk_frames >= cpuArg.chunkSize ){
            next_st = cur_frame_id -(cur_chunk_frames - cpuArg.chunkSize);
            if(next_st < total_frames){
                chunk_frame_st[count_chunk] = next_st;
                count_chunk++;
                cur_chunk_frames = (cur_frame_id - next_st > cpuArg.featCRange -1)?(cur_frame_id - next_st - cpuArg.featCRange +1):0;
            }
        }
    }
    total_chunks = count_chunk;
    total_samples = (total_chunks -1) *cpuArg.chunkSize + cur_chunk_frames;
    
    fprintf(fp_log, "Get chunk info over: Training sentences have %d chunks, %d samples.\n", total_chunks, total_samples);
    cpuArg.totalChunks = total_chunks;
    cpuArg.totalSamples = total_samples;
    fflush(fp_log);
}

/// get chunk nums and get frame start id for every chunk
void Interface::get_chunk_info_cv(CpuArg& cpuArg)
{
    char *p, *st, *en;
    int sentid;
    int count_chunk = 1;
    int cur_frames_num = 0;
    int cur_frames_lost =0;
    int cur_frame_id = 0;
    int cur_chunk_frames =0;
    int frames_inc;
    int next_st;

    cv_sent_st = cpuArg.train_cv_range_beg;
    cv_sent_en = cpuArg.train_cv_range_end;
    if(cv_sent_en < cv_sent_st || cv_sent_st < 0 || cv_sent_en >= total_sents){
        fprintf(fp_log,"cv sent range: %d to %d number error.\n",cv_sent_st, cv_sent_en);
        exit(0);
    }
    
    if(cv_sent_st ==0){
        cur_frame_id = 0;
    }
    else{
        cur_frame_id = framesBeforeSent[cv_sent_st -1];
    }
    cv_chunk_frame_st[0] = cur_frame_id;
        
    for(sentid = cv_sent_st; sentid <= cv_sent_en; sentid++)
    {
        frames_inc = framesBeforeSent[sentid] - cur_frame_id;
        cur_frame_id = framesBeforeSent[sentid];
        if(frames_inc >= cpuArg.featCRange){
            cur_frames_lost = cpuArg.featCRange -1;
        }
        else{
            cur_frames_lost = frames_inc;
        }
        cur_frames_num = frames_inc - cur_frames_lost;
        cur_chunk_frames += cur_frames_num;
        while(cur_chunk_frames >= cpuArg.chunkSize ){
            next_st = cur_frame_id -(cur_chunk_frames - cpuArg.chunkSize);
            if(next_st < total_frames){
                cv_chunk_frame_st[count_chunk] = next_st;
                count_chunk++;
                cur_chunk_frames = (cur_frame_id - next_st > cpuArg.featCRange -1)?(cur_frame_id - next_st - cpuArg.featCRange +1):0;
            }
        }
    }
    
    cv_total_chunks = count_chunk;
    cv_total_samples = (cv_total_chunks -1) *cpuArg.chunkSize + cur_chunk_frames;
    
    fprintf(fp_log, "Get cv chunk info over: CV sentences have %d chunks, %d samples.\n", cv_total_chunks, cv_total_samples);
    cpuArg.cvTotalChunks = cv_total_chunks;
    cpuArg.cvTotalSamples = cv_total_samples;
    fflush(fp_log);
}

///// read one chunck frames by index
/******change for jumpframe and discard  ******/
int Interface::Readchunk(int chunk_index, CpuArg& cpuArg, ChunkContainer& oneChunk)  
{
    long int offset;
    int frames_need_read;
    int frames_processed;
    int samples_in_chunk, tmp_samples_in_chunk;
    int size_perframe;
    int cur_sent;
    int cur_frame_of_sent;
    int cur_frame_id;
    int cur_sample;
    int *sample_index;
    float *dataori;
    int *targori;
    int i,j,k;
    
    if(chunk_index == total_chunks -1){
        frames_need_read = framesBeforeSent[sent_en] - chunk_frame_st[chunk_index];
        samples_in_chunk = total_samples - cpuArg.chunkSize *chunk_index;
    }
    else{
        samples_in_chunk = cpuArg.chunkSize;
        frames_need_read = chunk_frame_st[chunk_index +1] - chunk_frame_st[chunk_index];
    }
    
    ///Read data pfile
    size_perframe = (cpuArg.featDim +2) *sizeof(float);
    offset = PFILE_HEADER_SIZE + chunk_frame_st[chunk_index]* (long int)size_perframe;
    if(0 != fseek(fp_data, offset,0)){
        fprintf(fp_log,"data pfile cannot fseek to chunk %d.\n",chunk_index);
        exit(0);
    }
    dataori = new float[frames_need_read *(cpuArg.featDim +2)];
    fread((char *)dataori,size_perframe,frames_need_read,fp_data);
    swap32( (int*)&dataori[0]);
    cur_sent = *((int*) &dataori[0]);
    for(i =0;i < frames_need_read; i++){
        for(j =0; j< cpuArg.featDim;j++){
            swap32((int *) (&dataori[2+j +i*(cpuArg.featDim +2)]));
            dataori[2+j +i*(cpuArg.featDim +2)] -= mean[j];
            dataori[2+j +i*(cpuArg.featDim +2)] *= dVar[j];
        }
    }

    ///Read targ pfile
    size_perframe = (1 +2) *sizeof(int);
    offset = PFILE_HEADER_SIZE + chunk_frame_st[chunk_index]* (long int) size_perframe;
        
    if(0 != fseek(fp_targ, offset,0)){
        fprintf(fp_log,"targ pfile cannot fseek to chunk %d.\n",chunk_index);
        exit(0);
    }
    targori = new int[frames_need_read *(1 +2)];
    fread((char *)targori,size_perframe,frames_need_read,fp_targ);
    swap32( &targori[0]);
    cur_sent = targori[0];
    for(i =0;i < frames_need_read; i++){
        swap32(&targori[2 +i*(1 +2)]);
    }

    frames_processed = 0;
    cur_frame_id = chunk_frame_st[chunk_index];
    cur_sample = 0;
    tmp_samples_in_chunk = 0;
    int *discard_flag = new int [samples_in_chunk];
    memset(discard_flag,0,sizeof(int)*samples_in_chunk);

    while(frames_processed != frames_need_read){
        if(framesBeforeSent[cur_sent] > frames_need_read + chunk_frame_st[chunk_index]){
            cur_frame_of_sent = frames_need_read - frames_processed;
        }
        else{
            cur_frame_of_sent = framesBeforeSent[cur_sent] - cur_frame_id;
        }
    
        for(j =0; j<= cur_frame_of_sent - cpuArg.featCRange;j++){
            if (cur_sample%cpuArg.nmod == cpuArg.remainder)
            {
                int sample_lab = targori[(frames_processed +j + para->targ_offset) *(2 + 1) +2];
                for (i=0; i<cpuArg.discard_num; i++)
                {
                    if (sample_lab == cpuArg.discard_labs[i] && drand48()< cpuArg.discard_prob)
                    {
                        discard_flag[cur_sample] = 1;
                        break;
                    }
                }
                if (discard_flag[cur_sample] == 0)
                {
                    tmp_samples_in_chunk++;
                }
            }
            else
            {
                discard_flag[cur_sample] = 1;
            }
            cur_sample++;
        }

        cur_frame_id = framesBeforeSent[cur_sent];
        cur_sent++;
        frames_processed += cur_frame_of_sent;
    }

    samples_in_chunk = tmp_samples_in_chunk;
    sample_index = new int [samples_in_chunk];
    for(i=0; i< samples_in_chunk; i++){
        sample_index[i] = i;
    }
    GetRandIndex(sample_index, samples_in_chunk);

    frames_processed = 0;
    cur_frame_id = chunk_frame_st[chunk_index];
    cur_sample = 0;
    cur_sent = targori[0];
    int randIndex =0;
    int cur_sample_index;

    while(frames_processed != frames_need_read){
        if(framesBeforeSent[cur_sent] > frames_need_read + chunk_frame_st[chunk_index]){
            cur_frame_of_sent = frames_need_read - frames_processed;
        }
        else{
            cur_frame_of_sent = framesBeforeSent[cur_sent] - cur_frame_id;
        }
    
        for(j =0; j<= cur_frame_of_sent - cpuArg.featCRange;j++){
            if (discard_flag[cur_sample] == 0)
            {
                cur_sample_index = sample_index[randIndex];
                oneChunk.labelArr[cur_sample_index] = targori[(frames_processed +j + para->targ_offset) *(2 + 1) +2];
                for(i =0;i< cpuArg.featCRange;i++){
                    for(k=0;k< cpuArg.featDim;k++){
                        oneChunk.dataArr[cur_sample_index* cpuArg.dnnLayerArr[0] +k +i *cpuArg.featDim] = dataori[(frames_processed +j +i) *(2+cpuArg.featDim) +k+2];
                    }
                }
                randIndex++;
            }
            cur_sample++;
        }
        
        cur_frame_id = framesBeforeSent[cur_sent];
        cur_sent++;
        frames_processed += cur_frame_of_sent;
    }

    delete []sample_index;
    delete []dataori;
    delete []targori;
    delete []discard_flag;

    tmp_samples_in_chunk += oneChunk.samples_remain;
    int pre_samples_remain = oneChunk.samples_remain;
    //oneChunk.samples_remain = tmp_samples_in_chunk%cpuArg.bunchSize;
    //samples_in_chunk = tmp_samples_in_chunk - oneChunk.samples_remain;
    if (tmp_samples_in_chunk < cpuArg.bunchSize)
    {
        oneChunk.samples_remain = 0;
        samples_in_chunk = tmp_samples_in_chunk;
    } else {
        oneChunk.samples_remain = tmp_samples_in_chunk % cpuArg.bunchSize;
        samples_in_chunk = tmp_samples_in_chunk - oneChunk.samples_remain;
    }

    oneChunk.labelArr += (oneChunk.samples_remain - pre_samples_remain); 
    oneChunk.dataArr  += (oneChunk.samples_remain - pre_samples_remain)*cpuArg.dnnLayerArr[0];

    return samples_in_chunk;
}
/************/

/*
Modified version : for parallel read in
Huang Hua(huangh223@mail2.sysu.edu.cn), 2016-04-08, SYSU ASC16 Team
This function does not require chunkSize == bunchSize. When reading a chunk, 
it will read samples whoes index mod numN falls in [s_sample, e_sample), because
the mpi process need the samples in this range.
*/
int Interface::ReadchunkPara(int chunk_index, CpuArg& cpuArg, ChunkContainer& oneChunk, int s_sample, int e_sample)  
{
    long int offset;
    int frames_need_read;
    int frames_processed;
    int samples_in_chunk, tmp_samples_in_chunk;
    int size_perframe;
    int cur_sent;
    int cur_frame_of_sent;
    int cur_frame_id;
    int cur_sample;
    int *sample_index;
    float *dataori;
    int *targori;
    int i,j,k;
    
    if(chunk_index == total_chunks -1){
        frames_need_read = framesBeforeSent[sent_en] - chunk_frame_st[chunk_index];
        samples_in_chunk = total_samples - cpuArg.chunkSize *chunk_index;
    }
    else{
        samples_in_chunk = cpuArg.chunkSize;
        frames_need_read = chunk_frame_st[chunk_index +1] - chunk_frame_st[chunk_index];
    }
    
    ///Read data pfile
    size_perframe = (cpuArg.featDim +2) *sizeof(float);
    offset = PFILE_HEADER_SIZE + chunk_frame_st[chunk_index]* (long int)size_perframe;
    if(0 != fseek(fp_data, offset,0)){
        fprintf(fp_log,"data pfile cannot fseek to chunk %d.\n",chunk_index);
        exit(0);
    }
    int fp_data_start_pos = offset;
    dataori = new float[frames_need_read *(cpuArg.featDim +2)];
    //fread((char *)dataori,size_perframe,frames_need_read,fp_data);
    
    // Read the first element for the next two statements, which actually is useless...
    fread((char *)dataori, sizeof(float), 1, fp_data);
    swap32( (int*)&dataori[0]);
    cur_sent = *((int*) &dataori[0]);
    // and we will do the transform later after the read in
    /*
    for (i = 0; i < frames_need_read; i++)
    {
        for (j = 0; j < cpuArg.featDim; j++)
        {
            swap32((int *) (&dataori[2+j +i*(cpuArg.featDim +2)]));
            dataori[2+j +i*(cpuArg.featDim +2)] -= mean[j];
            dataori[2+j +i*(cpuArg.featDim +2)] *= dVar[j];
        }
    }
    */
    
    ///Read targ pfile
    size_perframe = (1 +2) *sizeof(int);
    offset = PFILE_HEADER_SIZE + chunk_frame_st[chunk_index]* (long int) size_perframe;
        
    if(0 != fseek(fp_targ, offset,0)){
        fprintf(fp_log,"targ pfile cannot fseek to chunk %d.\n",chunk_index);
        exit(0);
    }
    targori = new int[frames_need_read *(1 +2)];
    fread((char *)targori,size_perframe,frames_need_read,fp_targ);
    swap32( &targori[0]);
    cur_sent = targori[0];
    for(i =0;i < frames_need_read; i++){
        swap32(&targori[2 +i*(1 +2)]);
    }

    frames_processed = 0;
    cur_frame_id = chunk_frame_st[chunk_index];
    cur_sample = 0;
    tmp_samples_in_chunk = 0;
    int *discard_flag = new int [samples_in_chunk];
    memset(discard_flag,0,sizeof(int)*samples_in_chunk);

    while(frames_processed != frames_need_read){
        if(framesBeforeSent[cur_sent] > frames_need_read + chunk_frame_st[chunk_index]){
            cur_frame_of_sent = frames_need_read - frames_processed;
        }
        else{
            cur_frame_of_sent = framesBeforeSent[cur_sent] - cur_frame_id;
        }
    
        for(j =0; j<= cur_frame_of_sent - cpuArg.featCRange;j++){
            if (cur_sample%cpuArg.nmod == cpuArg.remainder)
            {
                int sample_lab = targori[(frames_processed +j + para->targ_offset) *(2 + 1) +2];
                for (i=0; i<cpuArg.discard_num; i++)
                {
                    if (sample_lab == cpuArg.discard_labs[i] && drand48()< cpuArg.discard_prob)
                    {
                        discard_flag[cur_sample] = 1;
                        break;
                    }
                }
                if (discard_flag[cur_sample] == 0)
                {
                    tmp_samples_in_chunk++;
                }
            }
            else
            {
                discard_flag[cur_sample] = 1;
            }
            cur_sample++;
        }

        cur_frame_id = framesBeforeSent[cur_sent];
        cur_sent++;
        frames_processed += cur_frame_of_sent;
    }

    samples_in_chunk = tmp_samples_in_chunk;
    sample_index = new int [samples_in_chunk];
    for(i=0; i< samples_in_chunk; i++){
        sample_index[i] = i;
    }
    //GetRandIndex(sample_index, samples_in_chunk);

    frames_processed = 0;
    cur_frame_id = chunk_frame_st[chunk_index];
    cur_sample = 0;
    cur_sent = targori[0];
    int randIndex =0;
    int cur_sample_index;

    while(frames_processed != frames_need_read){
        if(framesBeforeSent[cur_sent] > frames_need_read + chunk_frame_st[chunk_index]){
            cur_frame_of_sent = frames_need_read - frames_processed;
        }
        else{
            cur_frame_of_sent = framesBeforeSent[cur_sent] - cur_frame_id;
        }
    
        for (j = 0; j <= cur_frame_of_sent - cpuArg.featCRange; j++)
        {
            int cur_sample_mod = cur_sample % cpuArg.bunchSize;
            if (cur_sample_mod < s_sample || cur_sample_mod >= e_sample)
            {
                cur_sample++;
                continue;
            }
            if (discard_flag[cur_sample] == 0)
            {
                //cur_sample_index = sample_index[randIndex]; 
                cur_sample_index = cur_sample;
                
                oneChunk.labelArr[cur_sample_index] = targori[(frames_processed +j + para->targ_offset) *(2 + 1) +2];
                
                
                for (i = 0; i < cpuArg.featCRange; i++)
                {
                    int dataori_offset = (frames_processed + j + i) * (2 + cpuArg.featDim) + 2;
                    int dataArr_offset = cur_sample_index * cpuArg.dnnLayerArr[0] + i * cpuArg.featDim;
                    float *dataArr_ptr = &oneChunk.dataArr[dataArr_offset];
                    float *dataori_ptr = &dataori[dataori_offset];
                    
                    // Point to the file position according to the offset of dataori and read in
                    fseek(fp_data, fp_data_start_pos + dataori_offset * sizeof(float), 0);
                    fread((char*) dataArr_ptr, sizeof(float), cpuArg.featDim, fp_data);
                    
                    // Now do the transform
                    for (k = 0; k < cpuArg.featDim; k++)
                    {
                        swap32((int *) (&dataArr_ptr[k]));
                        dataArr_ptr[k] -= mean[k];
                        dataArr_ptr[k] *= dVar[k];
                    }
                    
                    /*
                    //for debug use
                    float max_diff = 0.0f;
                    for (k = 0; k < cpuArg.featDim; k++)
                    {
                        if (fabs(dataArr_ptr[k] - dataori_ptr[k]) > max_diff)
                            max_diff = fabs(dataArr_ptr[k] - dataori_ptr[k]);
                    }
                    if (max_diff > 1e-15)
                    {
                        printf("max abs diff = %e\n", max_diff);
                        assert(max_diff < 1e-15);
                    }
                    */
                    
                    //memcpy(dataArr_ptr, dataori_ptr, cpuArg.featDim * sizeof(float)); // An easier way
                    /*
                    for (k = 0; k < cpuArg.featDim; k++)
                    {
                        // The original data copy, equals to dataArr_ptr[k] = dataori_ptr[k];
                        oneChunk.dataArr[cur_sample_index* cpuArg.dnnLayerArr[0] +k +i *cpuArg.featDim] = dataori[(frames_processed +j +i) *(2+cpuArg.featDim) +k+2];
                    }
                    */
                }
                randIndex++;
            }
            cur_sample++;
        }
        
        cur_frame_id = framesBeforeSent[cur_sent];
        cur_sent++;
        frames_processed += cur_frame_of_sent;
    }

    delete []sample_index;
    delete []dataori;
    delete []targori;
    delete []discard_flag;

    tmp_samples_in_chunk += oneChunk.samples_remain;
    int pre_samples_remain = oneChunk.samples_remain;
    //oneChunk.samples_remain = tmp_samples_in_chunk%cpuArg.bunchSize;
    //samples_in_chunk = tmp_samples_in_chunk - oneChunk.samples_remain;
    if (tmp_samples_in_chunk < cpuArg.bunchSize)
    {
        oneChunk.samples_remain = 0;
        samples_in_chunk = tmp_samples_in_chunk;
    } else {
        oneChunk.samples_remain = tmp_samples_in_chunk % cpuArg.bunchSize;
        samples_in_chunk = tmp_samples_in_chunk - oneChunk.samples_remain;
    }

    oneChunk.labelArr += (oneChunk.samples_remain - pre_samples_remain); 
    oneChunk.dataArr  += (oneChunk.samples_remain - pre_samples_remain)*cpuArg.dnnLayerArr[0];

    return samples_in_chunk;
}

///// read one chunck frames by index
int Interface::Readchunk_cv(int chunk_index,CpuArg& cpuArg,ChunkContainer& oneChunk)  
{
    long int offset;
    int frames_need_read;
    int frames_processed;
    int samples_in_chunk;
    int size_perframe;
    int cur_sent;
    int cur_frame_of_sent;
    int cur_frame_id;
    int cur_sample;
    int *sample_index;
    float *dataori;
    int *targori;
    int i,j,k;
    
    para->indata = oneChunk.dataArr;
    para->targ = oneChunk.labelArr;
    ///Read data pfile
    size_perframe = (para->fea_dim +2) *sizeof(float);
    offset = PFILE_HEADER_SIZE + cv_chunk_frame_st[chunk_index]* (long int)size_perframe;

    if(chunk_index == cv_total_chunks -1){
        frames_need_read = framesBeforeSent[cv_sent_en] - cv_chunk_frame_st[chunk_index];
        samples_in_chunk = cv_total_samples - para->traincache *chunk_index;
    }
    else{
        samples_in_chunk = para->traincache;
        frames_need_read = cv_chunk_frame_st[chunk_index +1] - cv_chunk_frame_st[chunk_index];
    }
    
    if(0 != fseek(fp_data, offset,0)){
        fprintf(fp_log,"data pfile cannot fseek to chunk %d.\n",chunk_index);
        exit(0);
    }
    
    sample_index = new int [samples_in_chunk];
    for(i=0; i< samples_in_chunk; i++){
        sample_index[i] = i;
    }
    
    dataori = new float[frames_need_read *(para->fea_dim +2)];
    fread((char *)dataori,size_perframe,frames_need_read,fp_data);
    swap32( (int*)&dataori[0]);
    cur_sent = *((int*) &dataori[0]);
    for(i =0;i < frames_need_read; i++){
        for(j =0; j< para->fea_dim;j++){
            swap32((int *) (&dataori[2+j +i*(para->fea_dim +2)]));
            dataori[2+j +i*(para->fea_dim +2)] -= mean[j];
            dataori[2+j +i*(para->fea_dim +2)] *= dVar[j];
        }
    }

    frames_processed = 0;
    cur_frame_id = cv_chunk_frame_st[chunk_index];
    cur_sample = 0;

    while(frames_processed != frames_need_read){
        if(framesBeforeSent[cur_sent] > frames_need_read + cv_chunk_frame_st[chunk_index]){
            cur_frame_of_sent = frames_need_read - frames_processed;
        }
        else{
            cur_frame_of_sent = framesBeforeSent[cur_sent] - cur_frame_id;
        }
    
        for(j =0; j<= cur_frame_of_sent - para->fea_context;j++){
#if DEBUG_YGY
            fprintf(pCheckFile,"===%d===%d===\n",cur_sample,chunk_index);
#endif
            for(i =0;i< para->fea_context;i++){
                for(k=0;k< para->fea_dim;k++){
                    para->indata[sample_index[cur_sample]* cpuArg.dnnLayerArr[0] +k +i *para->fea_dim] = dataori[(frames_processed +j +i) *(2+para->fea_dim) +k+2];
#if DEBUG_YGY
                    fprintf(pCheckFile,"%1.6f ",dataori[(frames_processed+j+i)*(2+cpuArg.featDim)+k+2]);
#endif                 
                }
            }
#if DEBUG_YGY
            fprintf(pCheckFile,"\n");
            fflush(pCheckFile);
#endif
            cur_sample++;
        }

        cur_frame_id = framesBeforeSent[cur_sent];
        cur_sent++;
        frames_processed += cur_frame_of_sent;
    }

    ///Read targ pfile
    size_perframe = (1 +2) *sizeof(int);
    offset = PFILE_HEADER_SIZE + cv_chunk_frame_st[chunk_index]* (long int) size_perframe;
        
    if(0 != fseek(fp_targ, offset,0)){
        fprintf(fp_log,"targ pfile cannot fseek to chunk %d.\n",chunk_index);
        exit(0);
    }
    
    targori = new int[frames_need_read *(1 +2)];
    fread((char *)targori,size_perframe,frames_need_read,fp_targ);
    swap32( &targori[0]);
    cur_sent = targori[0];
    for(i =0;i < frames_need_read; i++){
        swap32(&targori[2 +i*(1 +2)]);
    }

    frames_processed = 0;
    cur_frame_id = cv_chunk_frame_st[chunk_index];
    cur_sample = 0;

    while(frames_processed != frames_need_read){
        if(framesBeforeSent[cur_sent] > frames_need_read + cv_chunk_frame_st[chunk_index]){
            cur_frame_of_sent = frames_need_read - frames_processed;
        }
        else{
            cur_frame_of_sent = framesBeforeSent[cur_sent] - cur_frame_id;
        }
    
        for(j =0; j<= cur_frame_of_sent - para->fea_context;j++){
            para->targ[sample_index[cur_sample]] = targori[(frames_processed +j + para->targ_offset) *(2 + 1) +2];
#if DEBUG_YGY
            fprintf(pCheckFile,"%4d ",oneChunk.labelArr[sample_index[cur_sample]]);
#endif 
            cur_sample++;
        }
#if DEBUG_YGY
            fprintf(pCheckFile,"\n");
            fflush(pCheckFile);
#endif
        
        cur_frame_id = framesBeforeSent[cur_sent];
        cur_sent++;
        frames_processed += cur_frame_of_sent;
    }

    delete []sample_index;
    delete []dataori;
    delete []targori;

    return samples_in_chunk;
}

void Interface::GetRandWeight(float *vec, float min, float max, int len)  ///// get randem vector with uniform distribution
{
    for(int i =0;i< len;i++)
    {
        vec[i] = drand48()*(max -min) +min;
    }
}

void Interface::GetRandIndex(int *vec, int len)
{
    int i;
    int idx;
    int tmp;
    for(i =0 ;i< len -1;i++){
        idx = lrand48() % (len-i);
        tmp = vec[idx];
        vec[idx] = vec[len -1 -i];
        vec[len -1 -i] = tmp;
    }
}

void Interface::get_uint(const char* hdr, const char* argname, unsigned int* val)
{
    const char* p;        // Pointer to argument
    int count = 0;        // Number of characters scanned

    // Find argument in header
    p = strstr(hdr, argname);
    if (p==NULL){
            fprintf(fp_log, "pfile header format is Not correct.\n");
          exit(0);
        }
    // Go past argument name
    p += strlen(argname);

    // Get value from stfing
    sscanf(p, " %u%n", val, &count);

    // We expect to pass one space, so need >1 characters for success.
    if (count <= 1){
        fprintf(fp_log, "%s num in pfile header is Not correct.\n",argname);
          exit(0);
    }
}

void Interface::read_tail(FILE *fp, long int file_offset, unsigned int sentnum, int *out)
{
    long int offset = file_offset;
    offset += 4;
    fseek(fp,offset,0);
    if(sentnum !=(fread((char*)out, sizeof(int), sentnum,fp))){
        fprintf(fp_log, "pfile tail is Not correct.\n");
        exit(0);
    }
    
    for(int i= 0; i< sentnum; i++){
        swap32((int*) (&out[i]));
    }
}
