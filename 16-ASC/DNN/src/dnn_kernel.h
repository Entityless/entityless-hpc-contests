/*------------------------------
@copyright (c) 2016 SYSU ASC16 Team

Last modified : 2016-03
Author(s) : Huang Hua(huangh223@mail2.sysu.edu.cn)

This file is part of DNNTK optimization codes for ASC16.
You can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the 
Free Software Foundation; either version 2 of the License,
or (at your option) any later version.
------------------------------*/

#ifndef DNN_KERNEL_H
#define DNN_KERNEL_H

#include "dnn_utility.h"
#include "my_dnn_util.h"

void dnnTrain(CpuArg& cpuArg, NodeArg &nodeArg, NodeInArg *bunchs, int total_bunchs, int _myid, int _nprocs, ChunkContainer& oneChunk);
void* KernelInitial(void* v_paras);
void GetConfig(int* faster_readfile, int total_bunchs, int _myid, int _nprocs, NodeArg &nodeArg);
void* KernelUninitial(void* no_use_ptr);
void OldMain(int _myid, int _nprocs, int argc, char* argv[], CpuArg& cpuArg, NodeArg &nodeArg, ChunkContainer& oneChunk);

struct KernelPara
{
    CpuArg* cpuArgp;
    NodeArg* nodeArgp;
    NodeInArg* bunchsp;
    int total_bunchs;
    int _myid;
    int _nprocs;
};

#endif
