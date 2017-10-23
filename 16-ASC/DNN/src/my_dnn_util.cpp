/*------------------------------
@copyright (c) 2016 SYSU ASC16 Team

Last modified : 2016-01
Author(s) : Huang Hua(huangh223@mail2.sysu.edu.cn)

This file is part of DNNTK optimization codes for ASC16.
You can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the 
Free Software Foundation; either version 2 of the License,
or (at your option) any later version.
------------------------------*/

#include "dnn_utility.h"
#include "my_dnn_util.h"
#include <stdlib.h>
#include <string.h>
#include "mkl.h"

int myFetchOneBunch(ChunkContainer &oneChunk, NodeInArg &nodeArg)
{
    int ret = 0;
    int readSize = 0;
    int ftrSize = oneChunk.featDim * oneChunk.featCRange;
    if (oneChunk.dataIndex>=oneChunk.dataSize) return 0;

    if ((oneChunk.dataIndex + oneChunk.bunchSize)>oneChunk.dataSize)
        readSize = oneChunk.dataSize - oneChunk.dataIndex;
    else
        readSize = oneChunk.bunchSize;
    nodeArg.numN = readSize;  //adjust nodeArg.numN according to the actual read size

    memcpy(nodeArg.d_T, oneChunk.labelArr + oneChunk.dataIndex, readSize * sizeof(int));
    memcpy(nodeArg.dd_X, oneChunk.dataArr + oneChunk.dataIndex*ftrSize, readSize * sizeof(float) * ftrSize);
    nodeArg.ftrSize = ftrSize;
    
    oneChunk.dataIndex += readSize; 
    return readSize;
}

void myInitNodeMem(const CpuArg &cpuArg, NodeInArg &nodeArg)
{
    nodeArg.d_X = (double*) _mm_malloc(sizeof(double) * cpuArg.featDim * cpuArg.featCRange
                                      * cpuArg.bunchSize, 512);
    nodeArg.dd_X = (float*) _mm_malloc(sizeof(float) * cpuArg.featDim * cpuArg.featCRange
                                      * cpuArg.bunchSize, 512);
    nodeArg.d_T = (int*) _mm_malloc(sizeof(int) * cpuArg.bunchSize, 512);
}

void myFreeNodeMem(NodeInArg &nodeArg)
{
    if (nodeArg.d_X != NULL) _mm_free(nodeArg.d_X);
    if (nodeArg.dd_X != NULL) _mm_free(nodeArg.dd_X);
    if (nodeArg.d_T != NULL) _mm_free(nodeArg.d_T);
}
