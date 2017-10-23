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


#ifndef MY_DNN_UTIL_H
#define MY_DNN_UTIL_H

#include "dnn_helper.h"
#include "dnn_utility.h"

struct NodeInArg
{
    double *d_X;  //input
    float *dd_X;
    int *d_T;    //target label of input
	int numN;    //miniBatch size
    int ftrSize; //feature size
}; 

extern "C" int myFetchOneBunch(ChunkContainer &oneChunk, NodeInArg &nodeArg);
extern "C" void myInitNodeMem(const CpuArg &cpuArg, NodeInArg &nodeArg);
extern "C" void myFreeNodeMem(NodeInArg &nodeArg);

#endif
