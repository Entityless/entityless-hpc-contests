#include "dnn_helper.h"
#include <stdlib.h>
#include <string.h>
#if DEBUG_YGY
#include <iostream>
#endif
#include <time.h>
#include "mkl.h"

static int GetNormData(CpuArg &cpuArg, NormArgument &normData);
static void UninitNormData(NormArgument &normData);
static void isLrtrim_blank(char* line);

int GetInitFileConfig(const char *initFileName, const char* initNum, CpuArg &cpuArg, int myid)
{
    int ret = 0;

    //read config file
    cpuArg.cvFlag = false;
    FILE *pInitFile = fopen(initFileName, "r");
    if(NULL == pInitFile)
    {
        printf("ERROR: Failed to open config file %s\n", initFileName);
        return -1;
    }
    char cfgline[MAXLINE];
    char initWeightBuf[MAXLINE]={0};
    char *p = NULL;

    /******add for jumpframe and discard  ******/
    cpuArg.nmod = 1;
    cpuArg.remainder = 0;
    
    cpuArg.discard_prob = 0.0 ;
    cpuArg.discard_num  = 0 ;
    cpuArg.discard_labs = NULL ;
    /************/

    while(fgets(cfgline,MAXLINE,pInitFile))
    {
        p = strstr(cfgline,"=");
        if(NULL == p)
            continue;
        *p='\0';
        p += 1;
        isLrtrim_blank(cfgline);
        isLrtrim_blank(p);
        if(0 == strcmp(cfgline,"fea_dim"))
        {
            cpuArg.featDim = atoi(p);
        }
        else if (0 == strcmp(cfgline,"fea_context"))
        {
            cpuArg.featCRange = atoi(p);
        }
        else if (0 == strcmp(cfgline,"init_random_seed"))
        {
            cpuArg.randomSeed = atoi(p);
        }
        else if (0 == strcmp(cfgline,"learn_rate"))
        {
            cpuArg.lRate = atof(p);
        }
        else if (0 == strcmp(cfgline,"bunchsize"))
        {
            cpuArg.bunchSize = atoi(p);
        }
        else if (0 == strcmp(cfgline,"chuncksize"))
        {
            cpuArg.chunkSize = atoi(p);
        }
        else if (0 == strcmp(cfgline,"file_norm"))
        {
            cpuArg.pNormFile = fopen(p,"r");
            if(NULL == cpuArg.pNormFile)
            {
                printf("ERROE: FAILED to open norm file %s\n",p);
                return -1;
            }
            cpuArg.normFileName = new char[strlen(p)+1];
            strcpy(cpuArg.normFileName,p);
        }
        else if (0 == strcmp(cfgline,"pfile_fea"))
        {
            cpuArg.pDataFile = fopen(p,"rb");
            if( NULL == cpuArg.pDataFile)
            {
                printf("ERROR: FAILED to open data file %s\n",p);
                return -1;
            }
            cpuArg.dataFileName = new char[strlen(p)+1];
            strcpy(cpuArg.dataFileName,p);
        }
        else if (0 == strcmp(cfgline,"pfile_lab"))
        {
            cpuArg.pLabelFile = fopen(p,"rb");
            if(NULL== cpuArg.pLabelFile)
            {
                printf("ERROR: FAILED to open label file %s\n",p);
                return -1;
            }
            cpuArg.labelFileName = new char[strlen(p)+1];
            strcpy(cpuArg.labelFileName,p);
        }
        else if (0 == strcmp(cfgline,"init_wts"))
        {
            cpuArg.weightFileName = new char[strlen(p)+1];
            strcpy(cpuArg.weightFileName,p);
        }
        else if (0 == strcmp(cfgline,"dnn_layers"))
        {
            char *tempBuff = p;
            p = strstr(p,"{");
            char *pp = p+1;
            int arrLen = 8;
            if(NULL == p)
            {
                printf("ERROR: FAILED to parse struct of %s\n",tempBuff);
                return -1;
            }
            cpuArg.dnnLayerArr = new int[8];
            cpuArg.dnnLayerArr[0] = atoi(pp);
            cpuArg.dnnLayerNum = 1;
            p = strstr(pp,",");
            while(NULL != p)
            {
                cpuArg.dnnLayerNum++;
                pp = p + 1;
                p = strstr(pp,",");
                if(cpuArg.dnnLayerNum>arrLen)
                {
                    int *newArr = new int[2*arrLen];
                    memcpy(newArr,cpuArg.dnnLayerArr,sizeof(int)*arrLen);
                    delete[] cpuArg.dnnLayerArr;
                    cpuArg.dnnLayerArr = newArr;
                    arrLen *= 2;
                }
                cpuArg.dnnLayerArr[cpuArg.dnnLayerNum-1] = atoi(pp);
            }
            cpuArg.outputLayerNum = cpuArg.dnnLayerArr[cpuArg.dnnLayerNum-1];
            if(cpuArg.dnnLayerNum < 3)
            {
                printf("ERROR: FAILED to parse struct of %s\n",tempBuff);
                return -1;
            }
        }
        else if(0 == strcmp(cfgline,"train_bp_range"))
        {
            cpuArg.train_bp_range_beg = atoi(p);
            p = strstr(p,"-") + 1;
            cpuArg.train_bp_range_end = atoi(p);
            if(cpuArg.train_bp_range_beg > cpuArg.train_bp_range_end
                || !(cpuArg.train_bp_range_beg >= 0 && (cpuArg.train_bp_range_end >= 0)))
            {
                printf("ERROR: error train_bp_range_beg %d or train_bp_range_end %d \n",
                cpuArg.train_bp_range_beg, cpuArg.train_bp_range_end);
                return -1;
            }
        }
        else if(0 == strcmp(cfgline,"train_cv_range"))
        {
            cpuArg.train_cv_range_beg = atoi(p);
            p = strstr(p,"-") + 1;
            cpuArg.train_cv_range_end = atoi(p);
            if(cpuArg.train_cv_range_beg > cpuArg.train_cv_range_end 
                || !(cpuArg.train_cv_range_beg >= 0 && cpuArg.train_cv_range_end >= 0))
            {
                printf("ERROR: error train_cv_range_beg %d or train_cv_range_end %d \n",
                cpuArg.train_cv_range_beg, cpuArg.train_cv_range_end);
                return -1;
            }
        }
        else if (0 == strcmp(cfgline,"output_wts"))
        {
            cpuArg.outputWtsPath = new char[MAXLINE];
            memset(cpuArg.outputWtsPath,0,MAXLINE*sizeof(char));
            sprintf(cpuArg.outputWtsPath,"%s.%s.wts",p,initNum);
            sprintf(initWeightBuf,"%s",p);
            cpuArg.pWtsFile = fopen(cpuArg.outputWtsPath,"w");
            if(NULL == cpuArg.pWtsFile)
            {
                printf("ERROR: FAILED to open weight file %s\n",
                cpuArg.outputWtsPath);
                return -1;
            }
        }
        else if (0 == strcmp(cfgline,"output_log"))
        {
            cpuArg.logFileName = new char[MAXLINE];
            memset(cpuArg.logFileName,0,MAXLINE*sizeof(char));
            sprintf(cpuArg.logFileName,"%s.%s.log",p,initNum);
            cpuArg.pLogFile = fopen(cpuArg.logFileName,"w");
            if(NULL == cpuArg.pLogFile)
            {
                printf("ERROR: FAILED to open log file %s\n",
                cpuArg.logFileName);
                return -1;
            }
        }
        else if(0 == strcmp(cfgline,"device_id"))
        {
            cpuArg.devId = atoi(p);
            if(cpuArg.devId < 0 && cpuArg.devId > 7)
            {
                printf("ERROR: error device id %d\n", cpuArg.devId);
                return -1;
            }
        }
        /****** add for jumpframe and discard  ******/
        else if(0 == strcmp(cfgline,"nmod"))
        {
            cpuArg.nmod = atoi(p);
        }
        else if(0 == strcmp(cfgline,"remainder"))
        {
            cpuArg.remainder = atoi(p);
        }

        else if(0 == strcmp(cfgline,"discard_prob"))
        {
            cpuArg.discard_prob = atof(p);
        }
        else if(0 == strcmp(cfgline,"discardLabs"))
        {    
            char *tempBuff = p;
            p = strstr(p,"{");
            char *pp = p+1;
            int arrLen = 4;
            if(NULL == p)
            {
                printf("ERROR: FAILED to parse struct of %s\n",tempBuff);
                return -1;
            }
            cpuArg.discard_labs = new int[arrLen];
            cpuArg.discard_labs[0] = atoi(pp);
            cpuArg.discard_num = 1;
            p = strstr(pp,",");
            while(NULL != p)
            {
                cpuArg.discard_num++;
                pp = p + 1;
                p = strstr(pp,",");
                if(cpuArg.discard_num >arrLen)
                {
                    int *newArr = new int[2*arrLen];
                    memcpy(newArr,cpuArg.discard_labs,sizeof(int)*arrLen);
                    delete[] cpuArg.discard_labs;
                    cpuArg.discard_labs = newArr;
                    arrLen *= 2;
                }
                cpuArg.discard_labs[cpuArg.discard_num -1] = atoi(pp);
            }
        }
        /************/
    }
    //end
    int i = atoi(initNum);
    if(i>0)
    {
        delete[] cpuArg.weightFileName;
        cpuArg.weightFileName = new char[MAXLINE];
        sprintf(cpuArg.weightFileName,"%s.%d.wts",initWeightBuf,i-1);
    }
    //read norm file;
    cpuArg.normInfo = new NormArgument;
    cpuArg.normInfo->meanArr = NULL;
    cpuArg.normInfo->meanArrLength = cpuArg.featDim;
    cpuArg.normInfo->dVarArr = NULL;
    cpuArg.normInfo->dVarArrLength = cpuArg.featDim;
    ret = GetNormData(cpuArg,*(cpuArg.normInfo));
    if(0 != ret)
        return ret;
    //init interface
    cpuArg.interface = new Interface();
    ret = cpuArg.interface->Initial(cpuArg);
    if(0 != ret)
        return ret;
    cpuArg.interface->get_pfile_info(cpuArg);
    cpuArg.interface->get_chunk_info(cpuArg);
    cpuArg.interface->get_chunk_info_cv(cpuArg);
    cpuArg.chunkIndex = new int[cpuArg.totalChunks];
    for (unsigned int i =0;i<cpuArg.totalChunks;i++)
        cpuArg.chunkIndex[i] = i;
    cpuArg.interface->GetRandIndex(cpuArg.chunkIndex,cpuArg.totalChunks);
    cpuArg.curChunkIndex = 0;
    if((0==cpuArg.totalChunks)||(0==cpuArg.cvTotalChunks))
    {
        fprintf(cpuArg.pLogFile,"no chunk train:%d cv:%d\n",
        cpuArg.totalChunks,cpuArg.cvTotalChunks);
    }
    return 0;
}

int GetNormData(CpuArg &cpuArg, NormArgument &normData)
{
    fprintf(cpuArg.pLogFile,"load norm file ......\n");
    char buff[MAXLINE];
    normData.meanArr = new float[normData.meanArrLength];
    normData.dVarArr = new float[normData.dVarArrLength];
    fgets(buff,MAXLINE,cpuArg.pNormFile);
    for (int i = 0;i <cpuArg.featDim;i++)
    {
        fgets(buff,MAXLINE,cpuArg.pNormFile);
        normData.meanArr[i] = atof(buff);
    }
    normData.meanArrLength = cpuArg.featDim;
    fgets(buff,MAXLINE,cpuArg.pNormFile);
    for (int i = 0;i <cpuArg.featDim;i++)
    {
        fgets(buff,MAXLINE,cpuArg.pNormFile);
        normData.dVarArr[i] = atof(buff);
    }
    normData.dVarArrLength = cpuArg.featDim;
    fprintf(cpuArg.pLogFile,"norm file load\n");
    fflush(cpuArg.pLogFile);
    return 0;
}

int InitNodeConfig(const CpuArg &cpuArg, NodeArg &nodeArg)
{    
    //open weight file
    FILE* pWeightFile = fopen(cpuArg.weightFileName, "rb");
    if(NULL == pWeightFile)
    {
        printf("ERROR: FAILED to open weight file %s\n",cpuArg.weightFileName);
        return -1;
    }
    int stat[10];
    char header[256];
    for (int i=1;i<cpuArg.dnnLayerNum;i++)
    {
        fread(stat,sizeof(int),5,pWeightFile);
        fread(header,sizeof(char),stat[4],pWeightFile);

        if((stat[1] != cpuArg.dnnLayerArr[i])||(stat[2] != cpuArg.dnnLayerArr[i-1]))
        {
            printf("ERROR:init weight node numbers do NOT match\n");
            return -1;
        }
        int N = stat[2]*stat[1];
        float *cpuWeightArr = new float[N];
        if(NULL == cpuWeightArr)
        {
            printf("ERROR:FAILED to malloc(cpu) memory\n");
            return -1;
        }
        fread(cpuWeightArr,sizeof(float),N,pWeightFile);

        /*
        nodeArg.d_W[i-1] = new float[N];
        nodeArg.d_Wdelta[i-1] = new float[N];
        */
        
        int arrSize = N * sizeof(float);
        nodeArg.d_W[i - 1] = (float*) _mm_malloc(arrSize, 256);
        nodeArg.d_Wdelta[i-1] =(float*) _mm_malloc(arrSize, 256);
        
        memset(nodeArg.d_Wdelta[i-1], 0, N*sizeof(float));
        memcpy(nodeArg.d_W[i-1],cpuWeightArr,N*sizeof(float));//TODO
        delete[] cpuWeightArr;
        cpuWeightArr = NULL;
        
        //bias
        fread(stat,sizeof(int),5,pWeightFile);
        fread(header,sizeof(char),stat[4],pWeightFile);

        if((stat[2] != cpuArg.dnnLayerArr[i])||(1 != stat[1]))
        {
            printf("ERROR:init bias node numbers do NOT match\n");
            return -1;
        }
        N = stat[2]*stat[1];
        float *cpuBiasArr = new float[N];
        fread(cpuBiasArr,sizeof(float),N,pWeightFile);

        /*
        nodeArg.d_B[i-1] = new float[N];
        nodeArg.d_Bdelta[i-1] = new float[N];
        nodeArg.d_Y[i-1] = new float[N*cpuArg.bunchSize];
        nodeArg.d_E[i-1]= new float[N*cpuArg.bunchSize];
        */
        int arrSize2 = sizeof(float) * N * cpuArg.bunchSize;
        nodeArg.d_B[i-1] = (float*) _mm_malloc(arrSize, 256);
        nodeArg.d_Bdelta[i-1] = (float*) _mm_malloc(arrSize, 256);
        nodeArg.d_Y[i-1] = (float*) _mm_malloc(arrSize2, 256);
        nodeArg.d_E[i-1] = (float*) _mm_malloc(arrSize2, 256);
        if (!nodeArg.d_B[i-1] || !nodeArg.d_Bdelta[i-1] || !nodeArg.d_Y[i-1] || !nodeArg.d_E[i-1])
        {
            printf("Error : mm_malloc for nodeArg failed!");
        }
        
        memset(nodeArg.d_Bdelta[i-1], 0, N*sizeof(float));
        memcpy(nodeArg.d_B[i-1],cpuBiasArr, N*sizeof(float));
        delete[] cpuBiasArr;
        cpuBiasArr = NULL;
    }
    //close weight file
    fclose(pWeightFile);
    pWeightFile = NULL;
    
    /*
    nodeArg.d_X = new float[cpuArg.featDim*cpuArg.featCRange*cpuArg.bunchSize];
    nodeArg.d_T = new int[cpuArg.bunchSize];
    nodeArg.d_R = new int[cpuArg.bunchSize];
    */
    
    nodeArg.d_X = (float*) _mm_malloc(sizeof(float) * cpuArg.featDim * cpuArg.featCRange 
                                      * cpuArg.bunchSize, 256);
    nodeArg.d_T = (int*) _mm_malloc(sizeof(int) * cpuArg.bunchSize, 256);
    nodeArg.d_R = (int*) _mm_malloc(sizeof(int) * cpuArg.bunchSize, 256);
    if (!nodeArg.d_X || !nodeArg.d_X || !nodeArg.d_X)
    {
        printf("Error : mm_malloc for nodeArg failed!");
    }
    
    
    nodeArg.lRate = cpuArg.lRate;
    nodeArg.numN = cpuArg.bunchSize;
    nodeArg.featDim = cpuArg.featDim;
    nodeArg.featCRange = cpuArg.featCRange;
    nodeArg.dnnLayerNum = cpuArg.dnnLayerNum;
    nodeArg.dnnLayerArr = new int[cpuArg.dnnLayerNum];
    memcpy(nodeArg.dnnLayerArr,cpuArg.dnnLayerArr,cpuArg.dnnLayerNum*sizeof(int));
    return 0;
}

int FetchOneChunk(CpuArg &cpuArg, ChunkContainer &oneChunk)
{
    if(1 != oneChunk.initFlag)
    {
        oneChunk.chunkSize = cpuArg.chunkSize;
        oneChunk.bunchSize = cpuArg.bunchSize;

        oneChunk.featDim = cpuArg.featDim;
        oneChunk.featCRange = cpuArg.featCRange;

        //****** change for jumpframe and discard  ******
        // add bunchsize to avoid bug
        oneChunk.dataArr = new float[(oneChunk.chunkSize + cpuArg.bunchSize)*oneChunk.featDim*oneChunk.featCRange];
        oneChunk.labelArr = new int[oneChunk.chunkSize + cpuArg.bunchSize];
        oneChunk.randArr = new int[oneChunk.chunkSize];
        oneChunk.samples_remain = 0 ;
        //************

        oneChunk.offset = 0;
        oneChunk.initFlag =1;
    }
    oneChunk.dataIndex = 0;
    if(cpuArg.cvFlag)
    {
        if(-1 == cpuArg.curChunkIndex)
            return 0;
        oneChunk.dataSize = 
        cpuArg.interface->Readchunk_cv(cpuArg.curChunkIndex,cpuArg,oneChunk);
        if(cpuArg.curChunkIndex == cpuArg.cvTotalChunks-1)
        {
            cpuArg.curChunkIndex = -1;
            return oneChunk.dataSize;
        }
    }
    else
    {
        if(cpuArg.curChunkIndex>=cpuArg.totalChunks)
        {
            cpuArg.curChunkIndex = 0;
            cpuArg.cvFlag = true;

            //****** add for jumpframe and discard  ******
            //reset oneChunk's labelArr dataArr for cross-vaild after training
            oneChunk.labelArr -= oneChunk.samples_remain ;
            oneChunk.dataArr -=  oneChunk.samples_remain*cpuArg.dnnLayerArr[0] ;
            oneChunk.samples_remain = 0 ;
            //************

            return 0;
        }
        oneChunk.dataSize = 
        cpuArg.interface->Readchunk(cpuArg.chunkIndex[cpuArg.curChunkIndex],cpuArg,oneChunk);
    }
    cpuArg.curChunkIndex++;
    return oneChunk.dataSize;
}

int FetchOneChunkPara(CpuArg &cpuArg, ChunkContainer &oneChunk, int s_sample, int e_sample)
{
    if(1 != oneChunk.initFlag)
    {
        oneChunk.chunkSize = cpuArg.chunkSize;
        oneChunk.bunchSize = cpuArg.bunchSize;

        oneChunk.featDim = cpuArg.featDim;
        oneChunk.featCRange = cpuArg.featCRange;

        //****** change for jumpframe and discard  ******
        // add bunchsize to avoid bug
        oneChunk.dataArr = new float[(oneChunk.chunkSize + cpuArg.bunchSize)*oneChunk.featDim*oneChunk.featCRange];
        oneChunk.labelArr = new int[oneChunk.chunkSize + cpuArg.bunchSize];
        oneChunk.randArr = new int[oneChunk.chunkSize];
        oneChunk.samples_remain = 0 ;
        //************

        oneChunk.offset = 0;
        oneChunk.initFlag =1;
    }
    oneChunk.dataIndex = 0;
    if(cpuArg.cvFlag)
    {
        if(-1 == cpuArg.curChunkIndex)
            return 0;
        oneChunk.dataSize = 
        cpuArg.interface->Readchunk_cv(cpuArg.curChunkIndex,cpuArg,oneChunk);
        if(cpuArg.curChunkIndex == cpuArg.cvTotalChunks-1)
        {
            cpuArg.curChunkIndex = -1;
            return oneChunk.dataSize;
        }
    }
    else
    {
        if(cpuArg.curChunkIndex>=cpuArg.totalChunks)
        {
            cpuArg.curChunkIndex = 0;
            cpuArg.cvFlag = true;

            //****** add for jumpframe and discard  ******
            //reset oneChunk's labelArr dataArr for cross-vaild after training
            oneChunk.labelArr -= oneChunk.samples_remain ;
            oneChunk.dataArr -=  oneChunk.samples_remain*cpuArg.dnnLayerArr[0] ;
            oneChunk.samples_remain = 0 ;
            //************

            return 0;
        }
        oneChunk.dataSize = 
        cpuArg.interface->ReadchunkPara(cpuArg.chunkIndex[cpuArg.curChunkIndex],cpuArg,oneChunk, s_sample, e_sample);
    }
    cpuArg.curChunkIndex++;
    return oneChunk.dataSize;
}

int FetchOneBunch(ChunkContainer &oneChunk, NodeArg &nodeArg)
{
    int ret = 0;
    int readSize = 0;
    int ftrSize = oneChunk.featDim*oneChunk.featCRange;
    if(oneChunk.dataIndex>=oneChunk.dataSize)
            return 0;

    if ((oneChunk.dataIndex + oneChunk.bunchSize)>oneChunk.dataSize)
        readSize = oneChunk.dataSize - oneChunk.dataIndex;
    else
        readSize = oneChunk.bunchSize;
    nodeArg.numN = readSize;//adjust nodeArg.numN according to the actual read size

    memcpy(nodeArg.d_T,oneChunk.labelArr+oneChunk.dataIndex,readSize*sizeof(int));
    memcpy(nodeArg.d_X,oneChunk.dataArr+oneChunk.dataIndex*ftrSize,readSize*sizeof(float)*ftrSize);
    
    oneChunk.dataIndex += readSize; 
    return readSize;
}

void UninitNormData(NormArgument &normData)
{
    if(NULL != normData.meanArr)
    {
        delete[] normData.meanArr;
        normData.meanArr = NULL;
    }
    if(NULL != normData.dVarArr)
    {
        delete[] normData.dVarArr;
        normData.dVarArr = NULL;
    }
    return;
}

void UninitProgramConfig(CpuArg &cpuArg, NodeArg &nodeArg, ChunkContainer &oneChunk)
{
    // uninit cpu config
    if(NULL != cpuArg.pDataFile)
    {
        fclose(cpuArg.pDataFile);
        cpuArg.pDataFile = NULL;
    }
    if(NULL != cpuArg.pLabelFile)
    {
        fclose(cpuArg.pLabelFile);
        cpuArg.pLabelFile = NULL;
    }
    if(NULL != cpuArg.pNormFile)
    {
        fclose(cpuArg.pNormFile);
        cpuArg.pNormFile == NULL;
    }
    if(NULL != cpuArg.pWtsFile)
    {
        fclose(cpuArg.pWtsFile);
        cpuArg.pWtsFile == NULL;
    }
    if(NULL != cpuArg.pLogFile)
    {
        fclose(cpuArg.pLogFile);
        cpuArg.pLogFile == NULL;
    }
    if(NULL != cpuArg.dataFileName)
    {
        delete[] cpuArg.dataFileName;
        cpuArg.dataFileName = NULL;
    }
    if(NULL != cpuArg.labelFileName)
    {
        delete[] cpuArg.labelFileName;
        cpuArg.labelFileName = NULL;
    }
    if(NULL != cpuArg.normFileName)
    {
        delete[] cpuArg.normFileName;
        cpuArg.normFileName = NULL;
    }
    if(NULL != cpuArg.weightFileName)
    {
        delete[] cpuArg.normFileName;
        cpuArg.normFileName = NULL;
    }
    if(NULL != cpuArg.logFileName)
    {
        delete[] cpuArg.logFileName;
        cpuArg.logFileName = NULL;
    }
    if(NULL != cpuArg.outputWtsPath)
    {
        delete[] cpuArg.outputWtsPath;
        cpuArg.outputWtsPath = NULL;
    }
    if(NULL != cpuArg.dnnLayerArr)
    {
        delete[] cpuArg.dnnLayerArr;
        cpuArg.dnnLayerArr = NULL;
    }
    if(NULL != cpuArg.normInfo)
    {
        UninitNormData(*(cpuArg.normInfo));
        delete cpuArg.normInfo;
        cpuArg.normInfo = NULL;
    }
    if(NULL != cpuArg.chunkIndex)
    {
        delete[] cpuArg.chunkIndex;
        cpuArg.chunkIndex = NULL;
    }
    if(NULL != cpuArg.interface)
    {
        delete cpuArg.interface;
        cpuArg.interface = NULL;
    }

    if(NULL != nodeArg.d_X)
    {
        //delete[] nodeArg.d_X;
        _mm_free(nodeArg.d_X);
        nodeArg.d_X = NULL;
    }
    if(NULL != nodeArg.d_T)
    {
        //delete[] nodeArg.d_T;
        _mm_free(nodeArg.d_T);
        nodeArg.d_T = NULL;
    }
    for(int i=0;i<NUM_LAYERS;i++)
    {
        if(NULL != nodeArg.d_W[i])
        {
            //delete[] nodeArg.d_W[i];
            _mm_free(nodeArg.d_W[i]);
            nodeArg.d_W[i] = NULL;
        }
        if(NULL != nodeArg.d_Wdelta[i])
        {
            //delete[] nodeArg.d_Wdelta[i];
            _mm_free(nodeArg.d_Wdelta[i]);
            nodeArg.d_Wdelta[i] = NULL;
        }
        if(NULL != nodeArg.d_B[i])
        {
            //delete[] nodeArg.d_B[i];
            _mm_free(nodeArg.d_B[i]);
            nodeArg.d_B[i] = NULL;
        }
        if(NULL != nodeArg.d_Bdelta[i])
        {
            //delete[] nodeArg.d_Bdelta[i];
            _mm_free(nodeArg.d_Bdelta[i]);
            nodeArg.d_Bdelta[i] = NULL;
        }
        if(NULL != nodeArg.d_Y[i])
        {
            //delete[] nodeArg.d_Y[i];
            _mm_free(nodeArg.d_Y[i]);
            nodeArg.d_Y[i] = NULL;
        }
        if(NULL != nodeArg.d_E[i])
        {
            //delete[] nodeArg.d_E[i];
            _mm_free(nodeArg.d_E[i]);
            nodeArg.d_E[i] = NULL;
        }
    }
    if(NULL != nodeArg.dnnLayerArr)
    {
        delete[] nodeArg.dnnLayerArr;
        nodeArg.dnnLayerArr = NULL;
    }

   //uninit chunck memory
    if(NULL == oneChunk.dataArr)
    {
        delete[] oneChunk.dataArr;
        oneChunk.dataArr = NULL;
    }
    if(NULL == oneChunk.labelArr)
    {
        delete[] oneChunk.labelArr;
        oneChunk.labelArr = NULL;
    }
    if(NULL == oneChunk.randArr)
    {
        delete[] oneChunk.randArr;
        oneChunk.randArr = NULL;
    }
    return;
}

int getnodeArg(char *filename, NodeArg &nodeArg)
{
    return 0;                    
}

int fetchMiniBatch(char *ftrfilename, char* labelfilename, NodeArg &nodeArg, long long size, long long offset)
{
    return 0;
}

void isLrtrim_blank(char* line)
{
    int i,j,k;
    i = 0;
    j = strlen(line)-1;
    while( (line[i] == ' ' || line[i] == '\t' || line[i] == '\r' || line[i] == '\n') && i <= j )
    {
        i++;
    }
    while( (line[j] == ' ' || line[j] == '\t' || line[j] == '\r' || line[j] == '\n') && j > i )
    {
        j--;
    }
    if (i == 0)
    {
        line[j-i+1] = '\0';
        return;
    }
    for(k = 0; i <= j; i++)
    {
        line[k++] = line[i];
    }
    line[k] = '\0';
}

extern "C" void WriteWts(NodeArg &nodeArg, CpuArg &cpuArg)
{
    //write W and  B to file
    float **d_W = nodeArg.d_W;
    float **d_B = nodeArg.d_B;

    int numN = nodeArg.numN; //size of minibatch
    int numL = nodeArg.dnnLayerNum - 1; //layer nums
    int numD = nodeArg.dnnLayerArr[0]; //node nums of input layer
    int*numA = nodeArg.dnnLayerArr; //node nums of hiden layer and output layer
    
    int   stat[10];
    char  head[256];

    FILE *fp = cpuArg.pWtsFile;

    for (int i = 0; i < numL; i++) {
        sprintf(head,"weights%d%d",i+1,i+2);
        stat[0] = 10;
        stat[1] = numA[i+1];//para->layersizes[i];
        stat[2] = numA[i];//para->layersizes[i -1];
        stat[3] = 0;
        stat[4] = strlen(head)+1;
        
        fwrite(stat,sizeof(int),5,fp);
        fwrite(head,sizeof(char),stat[4],fp);
        fwrite(d_W[i],sizeof(float),stat[2] *stat[1],fp);

        sprintf(head,"bias%d",i+2);
        stat[0] = 10;
        stat[1] = 1;
        stat[2] = numA[i+1];
        stat[3] = 0;
        stat[4] = strlen(head)+1;

        fwrite(stat,sizeof(int),5,fp);
        fwrite(head,sizeof(char),stat[4],fp);
        fwrite(d_B[i],sizeof(float),stat[2],fp);
    }
}
