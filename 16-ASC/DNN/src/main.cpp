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

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include "mpi.h"
#include "pthread.h"

#include "dnn_helper.h"
#include "dnn_kernel.h"

#define BLOCK_LOW(id,p,n) ( (id)*(n) / (p) )
#define BLOCK_HIGH(id,p,n) ( BLOCK_LOW((id)+1,p,n) - 1 )
#define BLOCK_SIZE(id,p,n) ( BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1 )

int main(int argc, char* argv[])
{
    int myid, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    puts(argv[0]);

    struct timeval timerStart, timerStop;
    gettimeofday(&timerStart, NULL);

    if((3 != argc )&&( 2 != argc) ){
        printf("ERROR: usage- <command> <config file> <iterator number>\n");
        exit(1);
    }


    CpuArg cpuArg = {0};
    ChunkContainer oneChunk = {0};
    NodeArg nodeArg = {0};

    OldMain(myid, nprocs, argc, argv, cpuArg, nodeArg, oneChunk);
    
    if (myid == 0)
    {
        fprintf(cpuArg.pLogFile,"training over\n");
        fflush(cpuArg.pLogFile);

        WriteWts(nodeArg, cpuArg);
        
        gettimeofday(&timerStop, NULL);
        float timerElapsed = 1000.0 * (timerStop.tv_sec - timerStart.tv_sec) + (timerStop.tv_usec - timerStart.tv_usec) / 1000.0;
        fprintf(cpuArg.pLogFile,"total time cost: %fms\n", timerElapsed);
        printf("total time cost: %fms\n", timerElapsed);
        fflush(cpuArg.pLogFile);

        UninitProgramConfig(cpuArg,nodeArg,oneChunk);
    }


    printf("node %d finish all uninitial\n", myid);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}
