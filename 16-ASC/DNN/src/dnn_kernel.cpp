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

#include "dnn_kernel.h"
#include "kernel_util.h"
#include "mic_kernel_util.h"
#include "sigmoid_util.h"
#include "dnn_helper.h"
#include "dnn_utility.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mkl.h"
#include "mpi.h"
#include "omp.h"
#include "pthread.h"

#define ALLOC  alloc_if (1) free_if (0)
#define REUSE  alloc_if (0) free_if (0)
#define FREE   alloc_if (0) free_if (1)
#define TMPUSE alloc_if (1) free_if (1)

#pragma offload_attribute(push,target(mic))

float* fhW, *fhB;
double *dY, *hW, *hB, *dE, alpha, ori_alpha;
double *d0X, *d1X, *d2X;
int   *d0T, *d1T, *d2T;
double *d0Wdelta, *d1Wdelta, *d2Wdelta;
double *d0Bdelta, *d1Bdelta, *d2Bdelta;
double *src_Wdelta, *dst_Wdelta, *src_Bdelta, *dst_Bdelta;//not to allocate. assigned by the master
int   numN, numL, numD, myid, nprocs, ori_numN, mic_numN, cpu_numN;
int   *Woffset, *Boffset, *numA;
int bunch_count;
NodeInArg* bunchs_ptr;

double *hE, *hY;
double *hWdelta, *hBdelta;
int *Yoffset;

int *start_rows, *row_counts;
int *mic_Yoffset0, *mic_Yoffset1, *mic_Yoffset2;
int *X_sizes, *T_sizes, *W_sizes, *B_sizes, *Y_sizes;

int signal0 = 0, signal1 = 1, signal2 = 2;

double one  = 1.0f;
double zero = 0.0f;

int* mic_identifer;

//flags
int *get_delta_flag0;
int *get_delta_flag1;
int *get_delta_flag2;
int *get_delta_flag_c;
int *forward_flag0;
int *forward_flag1;
int *forward_flag2;
int *backward_flag0;
int *backward_flag1;
int *backward_flag2;
int *out_flag;
int *inner_reduce_flag;
int *reduce_flag;
int *copy_flag;
int *in_flag;
int *assign_flag;
int *finish_flag;
int *receive_flag;
int *sum_flag;
int *special_barrier_assign_flag;
int *special_barrier_finish_flag;
int *read_flag;
int *use_flag0;
int *use_flag1;
int *use_flag2;
int *use_flag_c;
int dead_flag[3] = {0, 0, 0};


//constants
int worker_count = 3;
int dnn_thread_count;
int dnn_start_thread;
double cpu_mic_div = 0.92;
int max_bunch_count = -666;//but be positive if used'
bool dnn_on_cpu;

//debug flags
int *mic_current_bunch_;
int *mic_current_layer_;
int* comm_i_ptr_;
int* comm_j_ptr_;
int** trans_variances_;
int* master_i_ptr_;
int* master_j_ptr_;
int* transfer_stage_;
int master_stage_;

bool printed_ = false;
bool print_line_msg = true;

//other config
int div_reduce = 4;//1~20, the divide a stage of WBdelta into div_reduce parts
unsigned max_2_pow;//>= nprocs, = when nprocs is the exponent of 2
int sendrecv_loop;//if nprocs = 2, this is 1; =4, this is 2......
double learnrate_changerate;
double lr_destination;

//times
double total_time;
double total_comm_time;

int ConsumeSomeTime()
{
    int a = 1, b = 3, c = 2;
    for (int i = 0; i < 5000; i++)
    {
        a = b + c;
        c = a + b;
        b = a - c;
    }
    return a;
}

int ConsumeLittleTime()
{
    int a = 1, b = 3, c = 2;
    for (int i = 0; i < 10; i++)
    {
        a = b + c;
        c = a + b;
        b = a - c;
    }
    return a;
}

int GetSecondElement(int* arr, int deception_var)
{
    arr[0] = deception_var;
    arr[2] = arr[0] * deception_var;
    int new_arr[3];
    memcpy(new_arr, arr, 3 * sizeof(int));
    new_arr[0] = deception_var * new_arr[1];
    new_arr[2] = new_arr[0] * deception_var;
    return *(new_arr + 1);
}

#pragma offload_attribute(pop)

void isLrtrim_blank(char* line)//a function from dnn_utility.cpp, use it to read myconfig.txt
{
    int i, j, k;
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

void SpecialBarrier(int thread_count, int thread_id)
{
    bool round_finished;
    int loop_counter;

    if (thread_id == 0)
    {
        for (int i = 0; i < thread_count * 3; i++)
            special_barrier_assign_flag[i] = 1;

        special_barrier_finish_flag[1] = 1;

        while (true)
        {
            loop_counter++;
            round_finished = true;
            for (int k = 0; k < thread_count; k++)
            {
                if (GetSecondElement(special_barrier_finish_flag + k * 3, loop_counter) <= 0)
                {
                    round_finished = false;
                    break;
                }
            }
            if (round_finished)
                break;
        }
        for (int i = 0; i < thread_count * 3; i++)
            special_barrier_finish_flag[i] = 0;
    }
    else
    {
        while (true)
        {
            loop_counter++;
            if (GetSecondElement(special_barrier_assign_flag + thread_id * 3, loop_counter) > 0)
                break;
        }

        special_barrier_assign_flag[thread_id * 3 + 1] = 0;
        special_barrier_finish_flag[thread_id * 3 + 1] = 1;
    }
}

void CPUCheck2Arr(double* arr1, double* arr2, const char* c, int layer, int bunch)
{
    if (myid != 0)
        return;
    printf("%s, layer = %d, bunch = %d, %f %f %f %f %f %f\n",
            c, layer, bunch, arr1[0], arr1[1], arr1[2], arr2[0], arr2[1], arr2[2]
            );
}

void CPUCheck2Arrf(float* arr1, float* arr2, const char* c, int layer, int bunch)
{
    if (myid != 0)
        return;
    printf("%s, layer = %d, bunch = %d, %f %f %f %f %f %f\n",
            c, layer, bunch, arr1[0], arr1[1], arr1[2], arr2[0], arr2[1], arr2[2]
            );
}

void CPUCheckDelta(double* delta_arr, int my_micid, const char* c, int layer, int bunch)//do by the master line
{
    if (myid != 0)
        return;
    printf("%s, my_micid = %d, layer = %d, bunch = %d, %f %f %f %f %f %f %f %f %f %f\n",
            c, my_micid, layer, bunch, delta_arr[0], delta_arr[1], delta_arr[2], delta_arr[3], delta_arr[4], delta_arr[5], delta_arr[6], delta_arr[7], delta_arr[8], delta_arr[9]);
}

void GetLCrate(int total_bunchs)//get the reduce rate of learnrate, about 0.99xxx
{
    if (lr_destination < 0.99)
    {
        learnrate_changerate = pow(lr_destination, 1.0 / (double)total_bunchs);
        alpha *= ((1 - learnrate_changerate) * total_bunchs) / (1 - lr_destination);
    }
    else
        learnrate_changerate = 1.0f;
}

int BLOCK_SIZE(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) + ((n % total_blocks > block_id) ? 1 : 0);
}

int BLOCK_LOW(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) * block_id + ((n % total_blocks > block_id) ? block_id : n % total_blocks);
}

void CPUArrayInitial()
{
    row_counts   = (int*) malloc(sizeof(int) * 3);
    start_rows   = (int*) malloc(sizeof(int) * 3);
    X_sizes      = (int*) malloc(sizeof(int) * 9);
    T_sizes      = (int*) malloc(sizeof(int) * 9);
    W_sizes      = (int*) malloc(sizeof(int) * 9);
    B_sizes      = (int*) malloc(sizeof(int) * 9);
    Y_sizes      = (int*) malloc(sizeof(int) * 9);
    mic_Yoffset0 = (int*) malloc(sizeof(int) * (numL + 1));
    mic_Yoffset1 = (int*) malloc(sizeof(int) * (numL + 1));
    mic_Yoffset2 = (int*) malloc(sizeof(int) * (numL + 1));
    Yoffset      = (int*) malloc(sizeof(int) * (numL + 1));

    mic_identifer     = (int*) malloc(sizeof(int) * 3);
    get_delta_flag0   = (int*) malloc(sizeof(int) * 3);
    get_delta_flag1   = (int*) malloc(sizeof(int) * 3);
    get_delta_flag2   = (int*) malloc(sizeof(int) * 3);
    get_delta_flag_c  = (int*) malloc(sizeof(int) * 3);
    forward_flag0     = (int*) malloc(sizeof(int) * 3);
    forward_flag1     = (int*) malloc(sizeof(int) * 3);
    forward_flag2     = (int*) malloc(sizeof(int) * 3);
    backward_flag0    = (int*) malloc(sizeof(int) * 3);
    backward_flag1    = (int*) malloc(sizeof(int) * 3);
    backward_flag2    = (int*) malloc(sizeof(int) * 3);
    out_flag          = (int*) malloc(sizeof(int) * 9);
    inner_reduce_flag = (int*) malloc(sizeof(int) * 3);
    reduce_flag       = (int*) malloc(sizeof(int) * 3);
    copy_flag         = (int*) malloc(sizeof(int) * 3);
    in_flag           = (int*) malloc(sizeof(int) * 9);
    assign_flag       = (int*) malloc(sizeof(int) * (worker_count + 1) * 3);
    finish_flag       = (int*) malloc(sizeof(int) * (worker_count + 1) * 3);
    receive_flag      = (int*) malloc(3 * div_reduce * sendrecv_loop * sizeof(int));
    sum_flag          = (int*) malloc(3 * sendrecv_loop * sizeof(int));
    special_barrier_assign_flag = (int*) malloc(3 * dnn_thread_count *  sizeof(int));
    special_barrier_finish_flag = (int*) malloc(3 * dnn_thread_count *  sizeof(int));
    read_flag = (int*) malloc(3 * max_bunch_count * sizeof(int));
    use_flag0 = (int*) malloc(3 * max_bunch_count * sizeof(int));
    use_flag1 = (int*) malloc(3 * max_bunch_count * sizeof(int));
    use_flag2 = (int*) malloc(3 * max_bunch_count * sizeof(int));
    use_flag_c = (int*) malloc(3 * max_bunch_count * sizeof(int));

    mic_current_bunch_ = (int*) malloc(sizeof(int) * 9);
    mic_current_layer_ = (int*) malloc(sizeof(int) * 9);
    trans_variances_   = (int**) malloc(sizeof(int*) * 12);
    transfer_stage_    = (int*) malloc(sizeof(int) * 3);

    for (int i = 0; i < 3 * max_bunch_count; i++)
        read_flag[i] = use_flag0[i] = use_flag1[i] = use_flag2[i] = use_flag_c[i] = 0;
    for (int i = 0; i < 3; i++)
        get_delta_flag_c[i] = get_delta_flag0[i] = get_delta_flag1[i] = get_delta_flag2[i] = forward_flag0[i] = forward_flag1[i] = forward_flag2[i] = backward_flag0[i] = backward_flag1[i] = backward_flag2[i] = inner_reduce_flag[i] = reduce_flag[i] = copy_flag[i] = 0;
    for (int i = 0; i < 9; i++)
        out_flag[i] = in_flag[i] = 0;
    for (int i = 0; i < (worker_count + 1) * 3; i++)
        assign_flag[i] = finish_flag[i] = 0;
    for (int i = 0; i < 3 * div_reduce * sendrecv_loop; i++)
        receive_flag[i] = 0;
    for (int i = 0; i < 3 * sendrecv_loop; i++)
        sum_flag[i] = 0;
    for (int i = 0; i < 3 * dnn_thread_count; i++)
        special_barrier_assign_flag[i] = special_barrier_finish_flag[i] = 0;

    for (int t_mic = 0; t_mic < 3; t_mic++)
    {
        row_counts[t_mic] = BLOCK_SIZE(t_mic, 3, mic_numN);
        start_rows[t_mic] = BLOCK_LOW(t_mic, 3, mic_numN);

        for (int j = 0; j < 3; j++)
        {
            if (j == t_mic)
            {
                X_sizes[t_mic * 3 + j] = row_counts[t_mic] * numD;
                T_sizes[t_mic * 3 + j] = row_counts[t_mic];
                W_sizes[t_mic * 3 + j] = Woffset[numL];
                B_sizes[t_mic * 3 + j] = Boffset[numL];
            }
            else
            {
                X_sizes[t_mic * 3 + j] = 1;
                T_sizes[t_mic * 3 + j] = 1;
                W_sizes[t_mic * 3 + j] = 1;
                B_sizes[t_mic * 3 + j] = 1;
            }
            Y_sizes[t_mic * 3 + j] = 1;
        }
    }

    mic_Yoffset0[0] = 0;
    mic_Yoffset1[0] = 0;
    mic_Yoffset2[0] = 0;
    Yoffset[0] = 0;
    for (int i = 1; i < numL + 1; i++)
    {
        mic_Yoffset0[i] = mic_Yoffset0[i - 1] + row_counts[0] * numA[i - 1];
        mic_Yoffset1[i] = mic_Yoffset1[i - 1] + row_counts[1] * numA[i - 1];
        mic_Yoffset2[i] = mic_Yoffset2[i - 1] + row_counts[2] * numA[i - 1];
        Yoffset[i] = Yoffset[i - 1] + cpu_numN * numA[i - 1];
    }
    Y_sizes[0] = mic_Yoffset0[numL];
    Y_sizes[4] = mic_Yoffset1[numL];
    Y_sizes[8] = mic_Yoffset2[numL];

    dY       = (double*) _mm_malloc(sizeof(double) * mic_Yoffset0[numL], 512);
    dE       = (double*) _mm_malloc(sizeof(double) * mic_Yoffset0[numL], 512);
    hY       = (double*) _mm_malloc(sizeof(double) * Yoffset[numL], 512);
    hE       = (double*) _mm_malloc(sizeof(double) * Yoffset[numL], 512);
    d0Wdelta = (double*) _mm_malloc(sizeof(double) * Woffset[numL], 512);
    d1Wdelta = (double*) _mm_malloc(sizeof(double) * Woffset[numL], 512);
    d2Wdelta = (double*) _mm_malloc(sizeof(double) * Woffset[numL], 512);
    hWdelta  = (double*) _mm_malloc(sizeof(double) * Woffset[numL], 512);
    d0Bdelta = (double*) _mm_malloc(sizeof(double) * Boffset[numL], 512);
    d1Bdelta = (double*) _mm_malloc(sizeof(double) * Boffset[numL], 512);
    d2Bdelta = (double*) _mm_malloc(sizeof(double) * Boffset[numL], 512);
    hBdelta  = (double*) _mm_malloc(sizeof(double) * Boffset[numL], 512);
    d0X      = (double*) _mm_malloc(sizeof(double) * X_sizes[0], 512);
    d1X      = (double*) _mm_malloc(sizeof(double) * X_sizes[4], 512);
    d2X      = (double*) _mm_malloc(sizeof(double) * X_sizes[8], 512);
    d0T      = (int*)   _mm_malloc(sizeof(double) * T_sizes[0], 512);
    d1T      = (int*)   _mm_malloc(sizeof(double) * T_sizes[4], 512);
    d2T      = (int*)   _mm_malloc(sizeof(double) * T_sizes[8], 512);
    assert(dY && dE && hW && hB && hY && hE);
    assert(d0Wdelta && d1Wdelta && d2Wdelta && hWdelta);
    assert(d0Bdelta && d1Bdelta && d2Bdelta && hBdelta);
    assert(d0X && d1X && d2X && d0T && d1T && d2T);
}

void CPUArrayUninitial()
{
    free(row_counts);
    free(start_rows);
    free(X_sizes);
    free(T_sizes);
    free(W_sizes);
    free(B_sizes);
    free(Y_sizes);
    free(mic_Yoffset0);
    free(mic_Yoffset1);
    free(mic_Yoffset2);
    free(Yoffset);

    free(mic_identifer);

    free(get_delta_flag0);
    free(get_delta_flag1);
    free(get_delta_flag2);
    free(get_delta_flag_c);
    free(forward_flag0);
    free(forward_flag1);
    free(forward_flag2);
    free(backward_flag0);
    free(backward_flag1);
    free(backward_flag2);
    free(out_flag);
    free(inner_reduce_flag);
    free(copy_flag);
    free(reduce_flag);
    free(in_flag);
    free(assign_flag);
    free(finish_flag);
    free(receive_flag);
    free(sum_flag);
    free(special_barrier_assign_flag);
    free(special_barrier_finish_flag);

    free(mic_current_bunch_);
    free(mic_current_layer_);

    free(transfer_stage_);

    _mm_free(dY); _mm_free(dE);
    _mm_free(hY); _mm_free(hE);
    _mm_free(hWdelta);  _mm_free(hBdelta);
    _mm_free(d0Wdelta); _mm_free(d1Wdelta); _mm_free(d2Wdelta);
    _mm_free(d0Bdelta); _mm_free(d1Bdelta); _mm_free(d2Bdelta);
    _mm_free(d0X); _mm_free(d1X); _mm_free(d2X);
    _mm_free(d0T); _mm_free(d1T); _mm_free(d2T);
}

void MICDataInitial()
{
    for (int t_mic = 0; t_mic < 3; t_mic++)
    {
        //printf("node %d: small initial on %d\n", myid, t_mic);
        mic_identifer[1] = t_mic;
        #pragma offload_transfer target(mic:t_mic)\
        in     (mic_identifer      : length(3)    ALLOC) \
        in     (X_sizes        : length(9)    ALLOC) \
        in     (T_sizes        : length(9)    ALLOC) \
        in     (W_sizes        : length(9)    ALLOC) \
        in     (B_sizes        : length(9)    ALLOC) \
        in     (Y_sizes        : length(9)    ALLOC) \
        in     (row_counts     : length(3)    ALLOC) \
        in     (start_rows     : length(3)    ALLOC) \
        in     (Woffset        : length(numL + 1) ALLOC) \
        in     (Boffset        : length(numL + 1) ALLOC) \
        in     (numA           : length(numL)     ALLOC) \
        in     (mic_Yoffset0       : length(numL + 1) ALLOC) \
        in     (mic_Yoffset1       : length(numL + 1) ALLOC) \
        in     (mic_Yoffset2       : length(numL + 1) ALLOC) \
        in     (get_delta_flag0    : length(3)    ALLOC) \
        in     (get_delta_flag1    : length(3)    ALLOC) \
        in     (get_delta_flag2    : length(3)    ALLOC) \
        in     (forward_flag0      : length(3)    ALLOC) \
        in     (forward_flag1      : length(3)    ALLOC) \
        in     (forward_flag2      : length(3)    ALLOC) \
        in     (backward_flag0     : length(3)    ALLOC) \
        in     (backward_flag1     : length(3)    ALLOC) \
        in     (backward_flag2     : length(3)    ALLOC) \
        in     (in_flag        : length(9)    ALLOC) \
        nocopy (mic_current_bunch_ : length(9)    ALLOC) \
        nocopy (mic_current_layer_ : length(9)    ALLOC)
    }

    for (int t_mic = 0; t_mic < 3; t_mic++)
    {
        //int t_mic = omp_get_thread_num();
        int* local_Wsz = W_sizes + 3 * t_mic;
        int* local_Bsz = B_sizes + 3 * t_mic;
        int* local_Tsz = T_sizes + 3 * t_mic;
        int* local_Xsz = X_sizes + 3 * t_mic;
        int* local_Ysz = Y_sizes + 3 * t_mic;

        //printf("mic %d in hW with size %d, hB with size %d\n", t_mic, local_Wsz[t_mic], local_Bsz[t_mic]);

        #pragma offload_transfer target(mic:t_mic)\
        nocopy (d0X      : length(local_Xsz[0])      ALLOC  align(512)) \
        nocopy (d0T      : length(local_Tsz[0])      ALLOC  align(512)) \
        nocopy (d1X      : length(local_Xsz[1])      ALLOC  align(512)) \
        nocopy (d1T      : length(local_Tsz[1])      ALLOC  align(512)) \
        nocopy (d2X      : length(local_Xsz[2])      ALLOC  align(512)) \
        nocopy (d2T      : length(local_Tsz[2])      ALLOC  align(512)) \
        nocopy (d0Wdelta : length(local_Wsz[0])      ALLOC  align(512)) \
        nocopy (d0Bdelta : length(local_Bsz[0])      ALLOC  align(512)) \
        nocopy (d1Wdelta : length(local_Wsz[1])      ALLOC  align(512)) \
        nocopy (d1Bdelta : length(local_Bsz[1])      ALLOC  align(512)) \
        nocopy (d2Wdelta : length(local_Wsz[2])      ALLOC  align(512)) \
        nocopy (d2Bdelta : length(local_Bsz[2])      ALLOC  align(512)) \
        in     (hW       : length(local_Wsz[t_mic])  ALLOC  align(512)) \
        in     (hB       : length(local_Bsz[t_mic])  ALLOC  align(512)) \
        nocopy (dY       : length(local_Ysz[t_mic])  ALLOC  align(512)) \
        nocopy (dE       : length(local_Ysz[t_mic])  ALLOC  align(512))
    }

    printf("mic data initial of node %d over\n", myid);

    //for (int i = 0; i < 3; i++)
    //{
    //	#pragma offload_transfer target(mic:i) \
    //	    out(mic_identifer : length(3) REUSE)
    //    printf("node %d: mic %d, identifer = %d\n", myid, i, GetSecondElement(mic_identifer, i));
    //}
}

void MICDataOut()
{
    #pragma offload_transfer target(mic:0)\
    out    (hW       : length(Woffset[numL])  REUSE) \
    out    (hB       : length(Boffset[numL])  REUSE)
}

void MICDataUninitial()//finish this after finish the function above
{
    for (int t_mic = 0; t_mic < 3; t_mic++)
    {
        mic_identifer[1] = t_mic;
        #pragma offload_transfer target(mic:t_mic)\
        nocopy (mic_identifer      : length(3)    FREE) \
        nocopy (X_sizes            : length(9)    FREE) \
        nocopy (T_sizes            : length(9)    FREE) \
        nocopy (W_sizes            : length(9)    FREE) \
        nocopy (B_sizes            : length(9)    FREE) \
        nocopy (Y_sizes            : length(9)    FREE) \
        nocopy (row_counts         : length(3)    FREE) \
        nocopy (start_rows         : length(3)    FREE) \
        nocopy (Woffset            : length(numL + 1) FREE) \
        nocopy (Boffset            : length(numL + 1) FREE) \
        nocopy (numA               : length(numL)     FREE) \
        nocopy (mic_Yoffset0       : length(numL + 1) FREE) \
        nocopy (mic_Yoffset1       : length(numL + 1) FREE) \
        nocopy (mic_Yoffset2       : length(numL + 1) FREE) \
        nocopy (get_delta_flag0    : length(3)    FREE) \
        nocopy (get_delta_flag1    : length(3)    FREE) \
        nocopy (get_delta_flag2    : length(3)    FREE) \
        nocopy (forward_flag0      : length(3)    FREE) \
        nocopy (forward_flag1      : length(3)    FREE) \
        nocopy (forward_flag2      : length(3)    FREE) \
        nocopy (backward_flag0     : length(3)    FREE) \
        nocopy (backward_flag1     : length(3)    FREE) \
        nocopy (backward_flag2     : length(3)    FREE) \
        nocopy (in_flag            : length(9)    FREE) \
        nocopy (mic_current_bunch_ : length(9)    FREE) \
        nocopy (mic_current_layer_ : length(9)    FREE)
    }

    #pragma omp parallel num_threads(3)
    {
        int t_mic = omp_get_thread_num();
        int* local_Wsz = W_sizes + 3 * t_mic;
        int* local_Bsz = B_sizes + 3 * t_mic;
        int* local_Tsz = T_sizes + 3 * t_mic;
        int* local_Xsz = X_sizes + 3 * t_mic;
        int* local_Ysz = Y_sizes + 3 * t_mic;
        #pragma offload_transfer target(mic:t_mic)\
            nocopy (d0X      : length(local_Xsz[0])      FREE) \
            nocopy (d0T      : length(local_Tsz[0])      FREE) \
            nocopy (d1X      : length(local_Xsz[1])      FREE) \
            nocopy (d1T      : length(local_Tsz[1])      FREE) \
            nocopy (d2X      : length(local_Xsz[2])      FREE) \
            nocopy (d2T      : length(local_Tsz[2])      FREE) \
            nocopy (d0Wdelta : length(local_Wsz[0])      FREE) \
            nocopy (d0Bdelta : length(local_Bsz[0])      FREE) \
            nocopy (d1Wdelta : length(local_Wsz[1])      FREE) \
            nocopy (d1Bdelta : length(local_Bsz[1])      FREE) \
            nocopy (d2Wdelta : length(local_Wsz[2])      FREE) \
            nocopy (d2Bdelta : length(local_Bsz[2])      FREE) \
            nocopy (hW       : length(local_Wsz[t_mic])  FREE) \
            nocopy (hB       : length(local_Bsz[t_mic])  FREE) \
            nocopy (dY       : length(local_Ysz[t_mic])  FREE) \
            nocopy (dE       : length(local_Ysz[t_mic])  FREE)
    }
}

void copyDataToLocal(CpuArg& cpuArg, NodeArg &nodeArg)
{
    // !!! update the sample number to train each time with MPI worker id !!!
    Woffset = (int*) malloc(sizeof(int) * (numL + 1));
    Boffset = (int*) malloc(sizeof(int) * (numL + 1));
    Boffset[0] = Woffset[0] = 0;

    for (int i = 1; i <= numL; i++)
    {
        if (i == 1)
        {
            Woffset[i] = numD * numA[0];
            Boffset[i] = numA[0];
        } else {
            Woffset[i] = Woffset[i - 1] + numA[i - 2] * numA[i - 1];
            Boffset[i] = Boffset[i - 1] + numA[i - 1];
        }
    }

    fhW      = (float*) _mm_malloc(sizeof(float) * Woffset[numL], 512);
    fhB      = (float*) _mm_malloc(sizeof(float)* Boffset[numL], 512);

    // Allocate memory for local arrays
    hW      = (double*) _mm_malloc(sizeof(double) * Woffset[numL], 512);
    hB      = (double*) _mm_malloc(sizeof(double)* Boffset[numL], 512);



    assert(hB != NULL && hW != NULL);

    // Copy W B
    for (int i = 0; i < numL; i++)
    {
        if (i == 0)
        {
            memcpy(fhW, nodeArg.d_W[i], sizeof(float) * numD * numA[i]);
            memcpy(fhB, nodeArg.d_B[i], sizeof(float) * numA[i]);
        } else {
            memcpy(fhW + Woffset[i], nodeArg.d_W[i], sizeof(float) * numA[i - 1] * numA[i]);
            memcpy(fhB + Boffset[i], nodeArg.d_B[i], sizeof(float) * numA[i]);

        }
    }

    CPUCheck2Arrf(fhW, fhB, "initial fWB ", 0, 0);

    for(int i = 0; i < Woffset[numL]; i++)
        hW[i] = (double)fhW[i];
    for(int i = 0; i < Boffset[numL]; i++)
        hB[i] = (double)fhB[i];

    CPUCheck2Arr(hW, hB, "initial WB ", 0, 0);

}

void copyResultBack(NodeArg &nodeArg)
{
    for(int i = 0; i < Woffset[numL]; i++)
        fhW[i] = (float)hW[i];
    for(int i = 0; i < Boffset[numL]; i++)
        fhB[i] = (float)hB[i];

    CPUCheck2Arrf(fhW, fhB, "last fWB ", 0, 0);

    for (int i = 0; i < numL; i++)
    {
        if (i == 0)
        {
            memcpy(nodeArg.d_W[i], fhW, sizeof(float) * numD * numA[i]);
            memcpy(nodeArg.d_B[i], fhB, sizeof(float) * numA[i]);
        } else {
            memcpy(nodeArg.d_W[i], fhW + Woffset[i], sizeof(float) * numA[i - 1] * numA[i]);
            memcpy(nodeArg.d_B[i], fhB + Boffset[i], sizeof(float) * numA[i]);
        }
    }



    // Free local arrays
    if (!Boffset) free(Boffset);
    if (!Woffset) free(Woffset);
    if (!dY)      _mm_free(dY);
    if (!dE)      _mm_free(dE);
    if (!hW)      _mm_free(hW);
    if (!hB)      _mm_free(hB);

    _mm_free(fhW);
    _mm_free(fhB);
}

__attribute__((target(mic)))
void mic_dnnForward(
        int numL, int numN, int numD, int *numA,
        double *X, double *Y, double *W, double *B,
        int *Yoffset, int *Boffset, int *Woffset,
        double one, int my_micid, double* my_Wdelta, double* my_Bdelta, int* my_forward_flag, int* mic_current_bunch
        )
{
    __assume_aligned(W, 512);
    __assume_aligned(B, 512);
    __assume_aligned(X, 512);
    __assume_aligned(Y, 512);
    __assume_aligned(my_Wdelta, 512);
    __assume_aligned(my_Bdelta, 512);

    int flag_target, flag_value, loop_counter = 0;

    for (int i = 0; i < numL; i++)
    {
        mic_current_layer_[1 + 3 * my_micid] = i * 10;

        if (mic_current_bunch[1 + 3 * my_micid] > 0)
        {
            flag_target = 1 + i + mic_current_bunch[1 + 3 * my_micid] * numL;
            while (true)
            {
                loop_counter++;
                flag_value = GetSecondElement(in_flag + 3 * my_micid, loop_counter);
                if (flag_value >= flag_target)
                    break;
                ConsumeLittleTime();
            }
        }

        // Init the nodes with bias
        mic_setmatY(Y + Yoffset[i], B + Boffset[i], numN, numA[i]);

        if (i == 0)
        {
            cblas_dgemm(
                CblasRowMajor, CblasNoTrans, CblasNoTrans, numN, numA[i], numD,
                one, X, numD, W + Woffset[i], numA[i], one, dY + Yoffset[i], numA[i]
           );
        } else {
            cblas_dgemm(
                CblasRowMajor, CblasNoTrans, CblasNoTrans, numN, numA[i], numA[i-1],
                one, Y + Yoffset[i - 1], numA[i - 1], W + Woffset[i], numA[i], one, Y + Yoffset[i], numA[i]
           );
        }

        my_forward_flag[1] = mic_current_bunch_[1 + 3 * my_micid] + 1;

        if (i == numL - 1) { // softmax on output layer
            mic_softmaxZ(Y + Yoffset[i], Y + Yoffset[i], numN, numA[i]);
        } else {         // sigmod on hiden layers
            mic_sigmoidY(Y + Yoffset[i], numN, numA[i]);
        }
    }
}

    __attribute__((target(mic)))
void mic_dnnBackward(
        int numL, int numN, int *numA,
        double *E, double *Y, int *T, double *W,
        int *Eoffset, int *Yoffset, int *Woffset,
        double one, double zero, int my_micid, int* my_backward_flag, int* mic_current_bunch
        )
{
    __assume_aligned(W, 512);
    __assume_aligned(T, 512);
    __assume_aligned(E, 512);
    __assume_aligned(Y, 512);
    mic_errorOutput(E + Eoffset[numL - 1], Y + Yoffset[numL - 1], T, numN, numA[numL - 1]);
    my_backward_flag[1] = mic_current_bunch[1 + 3 * my_micid] + 1;
    for (int i = numL - 2; i >= 0; i--)
    {
        mic_current_layer_[1 + 3 * my_micid] = i * 1000;
        cblas_dgemm(
            CblasRowMajor, CblasNoTrans, CblasTrans, numN, numA[i], numA[i + 1],
            one, E + Eoffset[i + 1], numA[i + 1], W + Woffset[i + 1], numA[i + 1], zero, E + Eoffset[i], numA[i]
       );
        mic_errorTrans(E + Eoffset[i], Y + Yoffset[i], numN, numA[i]);
    }
}

    __attribute__((target(mic)))
void mic_dnnGetDelta(
        int numL, int numD, int numN, int *numA,
        double *X, double *E, double *Y, double *Wdelta, double *Bdelta,
        int *Eoffset, int *Yoffset, int *Woffset,
        double zero, double alpha, int my_micid, int* my_delta_flag, int* mic_current_bunch
        )
{
    __assume_aligned(Wdelta, 512);
    __assume_aligned(E, 512);
    __assume_aligned(X, 512);
    __assume_aligned(Y, 512);
    // Get Wdelta and Bdelta of input layer
    cblas_dgemm(
            CblasRowMajor, CblasTrans, CblasNoTrans, numD, numA[0], numN,
            alpha, X, numD, E + Eoffset[0], numA[0], zero, Wdelta + Woffset[0], numA[0]
           );
    mic_getBdelta(E + Eoffset[0], Bdelta, numN, numA[0], alpha);

    my_delta_flag[1] = 1 + 0 + mic_current_bunch[1 + 3 * my_micid] * numL;
    // Get Wdelta and Bdelta of hidden layer
    for (int i = 1; i < numL; i++)
    {
        mic_current_layer_[1 + 3 * my_micid] = i * 100000;
        cblas_dgemm(
            CblasRowMajor, CblasTrans, CblasNoTrans, numA[i - 1], numA[i], numN,
            alpha, Y + Yoffset[i - 1], numA[i - 1], E + Eoffset[i], numA[i], zero, Wdelta + Woffset[i], numA[i]
       );
        mic_getBdelta(E + Eoffset[i], Bdelta + Boffset[i], numN, numA[i], alpha);

        my_delta_flag[1] = 1 + i + mic_current_bunch[1 + 3 * my_micid] * numL;
    }
}

int TransferOutForwardFlag(int t_mic)
{
    if (t_mic == 0)
    {
        #pragma offload_transfer target(mic:0) \
        out    (forward_flag0      : length(3)           REUSE)
        return GetSecondElement(forward_flag0, t_mic);
    }
    else if (t_mic == 1)
    {
        #pragma offload_transfer target(mic:1) \
        out    (forward_flag1      : length(3)           REUSE)
        return GetSecondElement(forward_flag1, t_mic);
    }
    else
    {
        #pragma offload_transfer target(mic:2) \
        out    (forward_flag2      : length(3)           REUSE)
        return GetSecondElement(forward_flag2, t_mic);
    }
}

int TransferOutBackwardFlag(int t_mic)
{
    if (t_mic == 0)
    {
        #pragma offload_transfer target(mic:0) \
        out    (backward_flag0      : length(3)           REUSE)
        return GetSecondElement(backward_flag0, t_mic);
    }
    else if (t_mic == 1)
    {
        #pragma offload_transfer target(mic:1) \
        out    (backward_flag1      : length(3)           REUSE)
        return GetSecondElement(backward_flag1, t_mic);
    }
    else
    {
        #pragma offload_transfer target(mic:2) \
        out    (backward_flag2      : length(3)           REUSE)
        return GetSecondElement(backward_flag2, t_mic);
    }
}

int TransferOutGetDeltaFlag(int t_mic)
{
    if (t_mic == 0)
    {
        #pragma offload_transfer target(mic:0) \
        out    (get_delta_flag0    : length(3)           REUSE)
        return GetSecondElement(get_delta_flag0, t_mic);
    }
    else if (t_mic == 1)
    {
        #pragma offload_transfer target(mic:1) \
        out    (get_delta_flag1    : length(3)           REUSE)
        return GetSecondElement(get_delta_flag1, t_mic);
    }
    else
    {
        #pragma offload_transfer target(mic:2) \
        out    (get_delta_flag2    : length(3)           REUSE)
        return GetSecondElement(get_delta_flag2, t_mic);
    }
}

void TransferInX(int t_mic, NodeInArg *bunchs, int current_bunch)
{
    int MPI_numN_startpos = BLOCK_LOW(myid, nprocs, ori_numN);
    if (t_mic == 0)
    {
        memcpy(d0X, bunchs[current_bunch].d_X + (MPI_numN_startpos + start_rows[0]) * numD, X_sizes[0] * sizeof(double));
        #pragma offload_transfer target(mic:0) \
        in     (d0X        : length(X_sizes[0])      REUSE  align(512))
    }
    else if (t_mic == 1)
    {
        memcpy(d1X, bunchs[current_bunch].d_X + (MPI_numN_startpos + start_rows[1]) * numD, X_sizes[4] * sizeof(double));
        #pragma offload_transfer target(mic:1) \
        in     (d1X        : length(X_sizes[4])      REUSE  align(512))
    }
    else
    {
        memcpy(d2X, bunchs[current_bunch].d_X + (MPI_numN_startpos + start_rows[2]) * numD, X_sizes[8] * sizeof(double));
        #pragma offload_transfer target(mic:2) \
        in     (d2X        : length(X_sizes[8])      REUSE  align(512))
    }
}

void TransferInT(int t_mic, NodeInArg *bunchs, int current_bunch)
{
    int MPI_numN_startpos = BLOCK_LOW(myid, nprocs, ori_numN);
    if (t_mic == 0)
    {
        memcpy(d0T, bunchs[current_bunch].d_T + MPI_numN_startpos + start_rows[0], T_sizes[0] * sizeof(int));
        #pragma offload_transfer target(mic:0) \
        in     (d0T        : length(T_sizes[0])      REUSE  align(512))
    }
    else if (t_mic == 1)
    {
        memcpy(d1T, bunchs[current_bunch].d_T + MPI_numN_startpos + start_rows[1], T_sizes[4] * sizeof(int));
        #pragma offload_transfer target(mic:1) \
        in     (d1T        : length(T_sizes[4])      REUSE  align(512))
    }
    else
    {
        memcpy(d2T, bunchs[current_bunch].d_T + MPI_numN_startpos + start_rows[2], T_sizes[8] * sizeof(int));
        #pragma offload_transfer target(mic:2) \
        in     (d2T        : length(T_sizes[8])      REUSE  align(512))
    }
}

void TransferInDelta(int current_bunch, int current_layer, int t_mic)
{
#pragma offload_transfer target(mic:t_mic)\
    in  (hW[Woffset[current_layer]:Woffset[current_layer + 1] - Woffset[current_layer]] : REUSE) \
    in  (hB[Boffset[current_layer]:Boffset[current_layer + 1] - Boffset[current_layer]] : REUSE)

    in_flag[1 + 3 * t_mic] = 1 + current_layer + numL * (current_bunch + 1);

    #pragma offload_transfer target(mic:t_mic)\
        in  (in_flag[3 * t_mic : 3]    :  REUSE)
}

void TransferOutDelta(int current_layer, int t_mic, int current_bunch)
{
    if (t_mic == 0)
    {
        #pragma offload_transfer target(mic:0)\
        out  (d0Wdelta[Woffset[current_layer]:Woffset[current_layer + 1] - Woffset[current_layer]] : REUSE) \
        out  (d0Bdelta[Boffset[current_layer]:Boffset[current_layer + 1] - Boffset[current_layer]] : REUSE)
    }
    else if (t_mic == 1)
    {
        #pragma offload_transfer target(mic:1)\
        out  (d1Wdelta[Woffset[current_layer]:Woffset[current_layer + 1] - Woffset[current_layer]] : REUSE) \
        out  (d1Bdelta[Boffset[current_layer]:Boffset[current_layer + 1] - Boffset[current_layer]] : REUSE)
    }
    else
    {
        #pragma offload_transfer target(mic:2)\
        out  (d2Wdelta[Woffset[current_layer]:Woffset[current_layer + 1] - Woffset[current_layer]] : REUSE) \
        out  (d2Bdelta[Boffset[current_layer]:Boffset[current_layer + 1] - Boffset[current_layer]] : REUSE)
    }
}

void CPUCopyMission(double* dst, double* src, int length, int thread_count, int thread_id)
{
    int my_start_point = BLOCK_LOW(thread_id, thread_count, length);
    int my_length = BLOCK_SIZE(thread_id, thread_count, length);
    double* my_dst = dst + my_start_point;
    double* my_src = src + my_start_point;
    memcpy(my_dst, my_src, my_length * sizeof(double));
}

void CPUSumMission(double* dst, double* src, int length, int thread_count, int thread_id)
{
    int my_start_point = BLOCK_LOW(thread_id, thread_count, length);
    int my_length = BLOCK_SIZE(thread_id, thread_count, length);
    double* my_dst = dst + my_start_point;
    double* my_src = src + my_start_point;
    for (int i = 0; i < my_length; i++)
    {
        my_dst[i] += my_src[i];
    }
}
void CPUPrintErrorMessage()
{
    if (printed_)
        return;
    printed_ = true;

    int local_get_delta_flag0 = get_delta_flag0[1];
    int local_get_delta_flag1 = get_delta_flag1[1];
    int local_get_delta_flag2 = get_delta_flag2[1];
    int local_in_flag0 = in_flag[1];
    int local_in_flag1 = in_flag[4];
    int local_in_flag2 = in_flag[7];

    //printf("ready\n");

    for (int i = 0; i < 3; i++)
    {
        #pragma offload_transfer target(mic:i) \
        out(mic_identifer : length(3) REUSE)
        printf("node %d: mic %d, identifer = %d\n", myid, i, GetSecondElement(mic_identifer, i));
    }

    #pragma offload_transfer target(mic:0) \
    out    (get_delta_flag0    : length(3)           REUSE) \
    out    (in_flag[3 * 0 : 3]    :  REUSE) \
    out    (mic_current_bunch_[3 * 0 : 3]    :  REUSE) \
    out    (mic_current_layer_[3 * 0 : 3]    :  REUSE)

    #pragma offload_transfer target(mic:1) \
    out    (get_delta_flag1    : length(3)           REUSE) \
    out    (in_flag[3 * 1 : 3]    :  REUSE) \
    out    (mic_current_bunch_[3 * 1 : 3]    :  REUSE) \
    out    (mic_current_layer_[3 * 1 : 3]    :  REUSE)

    #pragma offload_transfer target(mic:2) \
    out    (get_delta_flag2    : length(3)           REUSE) \
    out    (in_flag[3 * 2 : 3]    :  REUSE) \
    out    (mic_current_bunch_[3 * 2 : 3]    :  REUSE) \
    out    (mic_current_layer_[3 * 2 : 3]    :  REUSE)

    //printf("ready2\n");

    printf("mic:0 bunch = %d, layer = %d, get_delta_flag = %d, in_flag = %d\nmic:1 bunch = %d, layer = %d, get_delta_flag = %d, in_flag = %d\nmic:2 bunch = %d, layer = %d, get_delta_flag = %d, in_flag = %d\nTransfer 0: out_bunch = %d, out_layer = %d, in_bunch = %d, in_layer = %d, copy_flag = %d, get_delta_flag = %d, stage = %d\nTransfer 1: out_bunch = %d, out_layer = %d, in_bunch = %d, in_layer = %d,  get_delta_flag = %d, stage = %d\nTransfer 2: out_bunch = %d, out_layer = %d, in_bunch = %d, in_layer = %d,  get_delta_flag = %d, stage = %d\nmaster: i = %d, j = %d, out_flag = %d, %d, %d, reduce_flag = %d, stage = %d\ncommline: i = %d, j = %d, inner_reduce_flag = %d\n\n",
            mic_current_bunch_[1], mic_current_layer_[1], get_delta_flag0[1], in_flag[1],
            mic_current_bunch_[4], mic_current_layer_[4], get_delta_flag1[1], in_flag[4],
            mic_current_bunch_[7], mic_current_layer_[7], get_delta_flag2[1], in_flag[7],
            trans_variances_[0][0], trans_variances_[1][0], trans_variances_[2][0], trans_variances_[3][0], copy_flag[1], local_get_delta_flag0, transfer_stage_[0],
            trans_variances_[4][0], trans_variances_[5][0], trans_variances_[6][0], trans_variances_[7][0], local_get_delta_flag1, transfer_stage_[1],
            trans_variances_[8][0], trans_variances_[9][0], trans_variances_[10][0], trans_variances_[11][0], local_get_delta_flag2, transfer_stage_[2],
            master_i_ptr_[0], master_j_ptr_[0], out_flag[1], out_flag[4], out_flag[7], reduce_flag[1], master_stage_,
            comm_i_ptr_[0], comm_j_ptr_[0], inner_reduce_flag[1]
          );

    printf("use flags: %d, %d, %d, %d\n", use_flag0[1], use_flag2[1], use_flag2[1], use_flag_c[1]);

    printf("I am dead now. Kill Me!\n");

    assert(0 > 1);
}

void meshC_RM_dgemm_threaded(
    int nthreads, int tid, int mesh_rows,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    int m, int n, int k,
    double alpha,
    double *a, int lda,
    double *b, int ldb,
    double beta,
    double *c, int ldc
)
{
    if (nthreads % mesh_rows)
    {
        fprintf(stderr, "meshC_RM_dgemm_threaded received wrong threads and row number of mesh\n");
        return;
    }

    int mesh_cols = nthreads / mesh_rows;
    int thread_mesh_col = tid % mesh_cols;
    int thread_mesh_row = tid / mesh_cols;
    int mesh_srow = BLOCK_LOW (thread_mesh_row, mesh_rows, m);
    int mesh_nrow = BLOCK_SIZE(thread_mesh_row, mesh_rows, m);
    int mesh_scol = BLOCK_LOW (thread_mesh_col, mesh_cols, n);
    int mesh_ncol = BLOCK_SIZE(thread_mesh_col, mesh_cols, n);

    double *top_left_a, *top_left_b, *top_left_c;
    if (transA == CblasNoTrans) top_left_a = a + mesh_srow * lda;
    else            top_left_a = a + mesh_srow;
    if (transB == CblasNoTrans) top_left_b = b + mesh_scol;
    else            top_left_b = b + mesh_scol * ldb;
    top_left_c = c + mesh_srow * ldc + mesh_scol;

    cblas_dgemm(
        CblasRowMajor, transA, transB,
        mesh_nrow, mesh_ncol, k,
        alpha,
        top_left_a, lda,
        top_left_b, ldb,
        beta,
        top_left_c, ldc
   );
}

void tCPUdnnForward(
    int nthreads, int tid, int mesh_rows,
    int numL, int numN, int numD, int *numA,
    double *X, double *Y, double *W, double *B,
    int *Yoffset, int *Boffset, int *Woffset,
    double one, int current_bunch
)
{
    int srow = BLOCK_LOW (tid, nthreads, numN);
    int nrow = BLOCK_SIZE(tid, nthreads, numN);
    int _2d_divN = mesh_rows;
    int n_mm_cols = nthreads / _2d_divN;
    int mm_srow = BLOCK_LOW (tid / n_mm_cols, _2d_divN, numN);
    int mm_nrow = BLOCK_SIZE(tid / n_mm_cols, _2d_divN, numN);

    int flag_target, flag_value, loop_counter = 0;

    for (int i = 0; i < numL; i++)
    {
        int mm_scol = BLOCK_LOW (tid % n_mm_cols, n_mm_cols, numA[i]);
        int mm_ncol = BLOCK_SIZE(tid % n_mm_cols, n_mm_cols, numA[i]);

        if (current_bunch != 0)
        {
            flag_target = 1 + i + current_bunch * numL;

            while (true)
            {
                loop_counter++;
                flag_value = GetSecondElement(copy_flag, loop_counter) + numL;
                if (flag_value >= flag_target)
                    break;
                ConsumeLittleTime();
            }
        }

        SpecialBarrier(dnn_thread_count, tid);

        CPUsetmatY(Y + Yoffset[i] + srow * numA[i], B + Boffset[i], nrow, numA[i]);

        if (i == 0)
        {
            meshC_RM_dgemm_threaded(
                nthreads, tid, mesh_rows,
                CblasNoTrans, CblasNoTrans,
                numN, numA[i], numD,
                one,
                X, numD,
                W + Woffset[i], numA[i],
                one,
                Y + Yoffset[i], numA[i]
            );
        } else {
            meshC_RM_dgemm_threaded(
                nthreads, tid, mesh_rows,
                CblasNoTrans, CblasNoTrans,
                numN, numA[i], numA[i - 1],
                one,
                Y + Yoffset[i - 1], numA[i - 1],
                W + Woffset[i], numA[i],
                one,
                Y + Yoffset[i], numA[i]
            );
        }

        SpecialBarrier(dnn_thread_count, tid);

        if (i < numL - 1)
        {
            CPUsigmoidY(Y + Yoffset[i] + mm_srow * numA[i] + mm_scol,
                    mm_nrow, mm_ncol, numA[i]);
        }
    }

    SpecialBarrier(dnn_thread_count, tid);

    CPUsoftmaxZ(
        Y + Yoffset[numL - 1] + srow * numA[numL - 1],
        Y + Yoffset[numL - 1] + srow * numA[numL - 1],
        nrow, numA[numL - 1]
    );
}

void tCPUdnnBackward(
    int nthreads, int tid, int mesh_rows,
    int numL, int numN, int *numA,
    double *E, double *Y, int *T, double *W,
    int *Eoffset, int *Yoffset, int *Woffset,
    double one, double zero
)
{
    int _2d_divN = mesh_rows;
    int srow = BLOCK_LOW (tid, nthreads, numN);
    int nrow = BLOCK_SIZE(tid, nthreads, numN);
    CPUerrorOutput(
        E + Eoffset[numL - 1] + srow * numA[numL - 1],
        Y + Yoffset[numL - 1] + srow * numA[numL - 1],
        T + srow, nrow, numA[numL - 1]
    );

    SpecialBarrier(dnn_thread_count, tid);

    int n_mm_cols = nthreads / _2d_divN;
    int mm_srow = BLOCK_LOW (tid / n_mm_cols, _2d_divN, numN);
    int mm_nrow = BLOCK_SIZE(tid / n_mm_cols, _2d_divN, numN);
    for (int i = numL - 2; i >= 0; i--)
    {
        int mm_scol = BLOCK_LOW (tid % n_mm_cols, n_mm_cols, numA[i]);
        int mm_ncol = BLOCK_SIZE(tid % n_mm_cols, n_mm_cols, numA[i]);

        SpecialBarrier(dnn_thread_count, tid);

        meshC_RM_dgemm_threaded(
            nthreads, tid, mesh_rows,
            CblasNoTrans, CblasTrans,
            numN, numA[i], numA[i + 1],
            one,
            E + Eoffset[i + 1], numA[i + 1],
            W + Woffset[i + 1], numA[i + 1],
            zero,
            E + Eoffset[i], numA[i]
        );

        SpecialBarrier(dnn_thread_count, tid);

        CPUerrorTrans(
            E + Eoffset[i] + mm_srow * numA[i] + mm_scol,
            Y + Yoffset[i] + mm_srow * numA[i] + mm_scol,
            mm_nrow, mm_ncol, numA[i]
        );
    }
}

void tCPUdnnGetDelta(
    int nthreads, int tid, int mesh_rows,
    int numL, int numD, int numN, int *numA,
    double *W, double *B, double *X, double *E,
    double *Y, double *Wdelta, double *Bdelta,
    int *Eoffset, int *Yoffset, int *Woffset,
    double zero, double alpha, int current_bunch
)
{
    int i = 0;

    meshC_RM_dgemm_threaded(
        nthreads, tid, mesh_rows,
        CblasTrans, CblasNoTrans,
        numD, numA[0], numN,
        alpha,
        X, numD,
        E + Eoffset[0], numA[0],
        zero,
        Wdelta + Woffset[0], numA[0]
    );

    int scol = BLOCK_LOW (tid, nthreads, numA[i]);
    int ncol = BLOCK_SIZE(tid, nthreads, numA[i]);
    CPUgetBdelta(
        E + Eoffset[i], Bdelta + Boffset[i],
        numN, numA[i], alpha, scol, scol + ncol
    );

    SpecialBarrier(dnn_thread_count, tid);

    if (tid == 0)
        get_delta_flag_c[1] = 1 + i + current_bunch * numL;

    for (i = 1; i < numL; i++)
    {
        meshC_RM_dgemm_threaded(
            nthreads, tid, mesh_rows,
            CblasTrans, CblasNoTrans,
            numA[i - 1], numA[i], numN,
            alpha,
            Y + Yoffset[i - 1], numA[i - 1],
            E + Eoffset[i], numA[i],
            zero,
            Wdelta + Woffset[i], numA[i]
        );

        scol = BLOCK_LOW (tid, nthreads, numA[i]);
        ncol = BLOCK_SIZE(tid, nthreads, numA[i]);
        CPUgetBdelta(
            E + Eoffset[i], Bdelta + Boffset[i],
            numN, numA[i], alpha, scol, scol + ncol
        );

        SpecialBarrier(dnn_thread_count, tid);

        if (tid == 0)
            get_delta_flag_c[1] = 1 + i + current_bunch * numL;

    }
}

void CPUDnnLine()
{
    NodeInArg *bunchs = bunchs_ptr;
    //int thread_id = omp_get_thread_num();//starts from 0
    int thread_id = omp_get_thread_num() - dnn_start_thread;
    if (myid == 0 && print_line_msg)
        printf("CPUDnnline, thread_num = %d, thread_id = %d\n", omp_get_thread_num(), thread_id);
    int mesh_rows = 4;
    int MPI_numN_startpos = BLOCK_LOW(myid, nprocs, ori_numN);

    double *my_hX;
    int *my_hT;

    int bunch_iter, bunch_overlap;
    int loop_counter = 0, flag_target;
    int flag_value;

    for (int i = 0; i < bunch_count; i++)
    {
        bunch_iter = i % max_bunch_count;
        bunch_overlap = i / max_bunch_count;

        flag_target = bunch_overlap;
        while(true)
        {
            loop_counter++;
            flag_value = GetSecondElement(read_flag + 3 * bunch_iter, loop_counter);
            if(flag_value >= flag_target)
                break;
        }

        my_hT = bunchs[bunch_iter].d_T + (MPI_numN_startpos + mic_numN);
        my_hX = bunchs[bunch_iter].d_X + (MPI_numN_startpos + mic_numN) * numD;

        SpecialBarrier(dnn_thread_count, thread_id);

        tCPUdnnForward (
            dnn_thread_count, thread_id, mesh_rows,
            numL, cpu_numN, numD, numA,
            my_hX, hY, hW, hB,
            Yoffset, Boffset, Woffset,
            one, i
       );

        tCPUdnnBackward(
            dnn_thread_count, thread_id, mesh_rows,
            numL, cpu_numN, numA,
            hE, hY, my_hT, hW,
            Yoffset, Yoffset, Woffset,
            one, zero
        );

        tCPUdnnGetDelta(
            dnn_thread_count, thread_id, mesh_rows,
            numL, numD, cpu_numN, numA,
            hW, hB, my_hX, hE,
            hY, hWdelta, hBdelta,
            Yoffset, Yoffset, Woffset,
            zero, alpha, i
        );

        use_flag_c[1 + 3 * bunch_iter] = bunch_overlap + 1;

        if(lr_destination < 0.98f)
            alpha *= learnrate_changerate;
    }

    //printf("thread %d escape from dnn line\n", omp_get_thread_num());
}

void CPUWorkerLine(void* ptr)
{

    int thread_id = *((int*)(ptr));//the master is with thread_id 0

    printf("thread %d enter worker with thread id = %d\n", omp_get_thread_num(), thread_id);

    int flag_target, flag_value, loop_counter = 0;
    int j_Woffset, j_Boffset;
    int j_Wlength, j_Blength;
    int div_Woffset, div_Boffset;
    int div_Wlength, div_Blength;

    for (int i = 0; i < bunch_count; i++)
    {
        for (int j = 0; j < numL; j++)
        {
            j_Woffset = Woffset[j];
            j_Wlength = Woffset[j + 1] - Woffset[j];
            j_Boffset = Boffset[j];
            j_Blength = Boffset[j + 1] - Boffset[j];
            loop_counter = 0;

            while (true)//1
            {
                loop_counter++;
                flag_value = GetSecondElement(assign_flag + thread_id * 3, loop_counter);
                if (flag_value > j + i * numL)
                    break;
            }
            //printf("worker %d has been given mission 1, i = %d, j = %d\n", thread_id, i, j);
            CPUSumMission(dst_Wdelta + j_Woffset, src_Wdelta + j_Woffset, j_Wlength, worker_count + 1, thread_id);
            CPUSumMission(dst_Bdelta + j_Boffset, src_Bdelta + j_Boffset, j_Blength, worker_count + 1, thread_id);
            //printf("worker %d has finish mission 1, i = %d, j = %d\n", thread_id, i, j);
            assign_flag[thread_id * 3 + 1] = j + i * numL;
            finish_flag[thread_id * 3 + 1] = j + i * numL + 1;

            //printf("worker %d is waiting for mission 2, i = %d, j = %d\n", thread_id, i, j);
            while (true)//2
            {
                loop_counter++;
                flag_value = GetSecondElement(assign_flag + thread_id * 3, loop_counter);
                if (flag_value > j + i * numL)
                    break;
            }
            //printf("worker %d has been given mission 2, i = %d, j = %d\n", thread_id, i, j);
            CPUSumMission(dst_Wdelta + j_Woffset, src_Wdelta + j_Woffset, j_Wlength, worker_count + 1, thread_id);
            CPUSumMission(dst_Bdelta + j_Boffset, src_Bdelta + j_Boffset, j_Blength, worker_count + 1, thread_id);
            //printf("worker %d has finish mission 2, i = %d, j = %d\n", thread_id, i, j);
            assign_flag[thread_id * 3 + 1] = j + i * numL;
            finish_flag[thread_id * 3 + 1] = j + i * numL + 1;

            //printf("worker %d is waiting for mission 3, i = %d, j = %d\n", thread_id, i, j);
            if(dnn_on_cpu)
            {
                while (true)//3
                {
                    loop_counter++;
                    flag_value = GetSecondElement(assign_flag + thread_id * 3, loop_counter);
                    if (flag_value > j + i * numL)
                        break;
                }
                //printf("worker %d has been given mission 3, i = %d, j = %d\n", thread_id, i, j);
                CPUSumMission(dst_Wdelta + j_Woffset, src_Wdelta + j_Woffset, j_Wlength, worker_count + 1, thread_id);
                CPUSumMission(dst_Bdelta + j_Boffset, src_Bdelta + j_Boffset, j_Blength, worker_count + 1, thread_id);
                //printf("worker %d has finish mission 3, i = %d, j = %d\n", thread_id, i, j);
                assign_flag[thread_id * 3 + 1] = j + i * numL;
                finish_flag[thread_id * 3 + 1] = j + i * numL + 1;
            }


            if (max_2_pow == nprocs && nprocs > 1)
            {
                for (int k = 0; k < sendrecv_loop; k++)
                {
                    for (int div_no = 0; div_no < div_reduce; ++div_no)
                    {
                        while (true)
                        {
                            loop_counter++;
                            flag_value = GetSecondElement(assign_flag + thread_id * 3, loop_counter);
                            if (flag_value > j + i * numL)
                                break;
                        }

                        div_Wlength = BLOCK_SIZE(div_no, div_reduce, j_Wlength);
                        div_Woffset = BLOCK_LOW(div_no, div_reduce, j_Wlength);
                        div_Blength = BLOCK_SIZE(div_no, div_reduce, j_Blength);
                        div_Boffset = BLOCK_LOW(div_no, div_reduce, j_Blength);

                        CPUSumMission(d0Wdelta + j_Woffset + div_Woffset, d1Wdelta + div_Woffset, div_Wlength, worker_count + 1, thread_id);
                        CPUSumMission(d0Bdelta + j_Boffset + div_Boffset, d1Bdelta + div_Boffset, div_Blength, worker_count + 1, thread_id);

                        assign_flag[thread_id * 3 + 1] = j + i * numL;
                        finish_flag[thread_id * 3 + 1] = j + i * numL + 1;
                    }
                }
            }

            while (true)
            {
                loop_counter++;
                flag_value = GetSecondElement(assign_flag + thread_id * 3, loop_counter);
                if (flag_value > j + i * numL)
                    break;
            }

            CPUSumMission(hW + j_Woffset, d0Wdelta + j_Woffset, j_Wlength, worker_count + 1, thread_id);
            CPUSumMission(hB + j_Boffset, d0Bdelta + j_Boffset, j_Blength, worker_count + 1, thread_id);

            assign_flag[thread_id * 3 + 1] = j + i * numL;
            finish_flag[thread_id * 3 + 1] = j + i * numL + 1;
        }
    }

    printf("thread %d escape from worker line\n", omp_get_thread_num());
}

void CPUCommunicateLine()//thread_num = 0
{
    double rst;
    unsigned bitmask;
    int partner;//sendrecv partner
    int div_Woffset, div_Boffset;
    int div_Wlength, div_Blength;
    int j_Woffset, j_Boffset;
    int j_Wlength, j_Blength;
    int i, j, k, div_no, flag_target, flag_value, loop_counter = 0;
    comm_i_ptr_ = &i;
    comm_j_ptr_ = &j;

    for (i = 0; i < bunch_count; i++)
    {
        for (j = 0; j < numL; j++)
        {
            j_Woffset = Woffset[j];
            j_Wlength = Woffset[j + 1] - Woffset[j];
            j_Boffset = Boffset[j];
            j_Blength = Boffset[j + 1] - Boffset[j];

            flag_target = j + 1 + i * numL;
            while (true)
            {
                loop_counter++;
                flag_value = GetSecondElement(inner_reduce_flag, loop_counter);
                if (flag_value >= flag_target)
                    break;
            }
            rst = omp_get_wtime();
            if (nprocs > 1 && max_2_pow == nprocs)//sendrecv
            {
                bitmask = 1;
                for (k = 0; k < sendrecv_loop; k++)
                {
                    partner = myid ^ bitmask;

                    for (div_no = 0; div_no < div_reduce; ++div_no)
                    {
                        div_Wlength = BLOCK_SIZE(div_no, div_reduce, j_Wlength);
                        div_Woffset = BLOCK_LOW(div_no, div_reduce, j_Wlength);
                        div_Blength = BLOCK_SIZE(div_no, div_reduce, j_Blength);
                        div_Boffset = BLOCK_LOW(div_no, div_reduce, j_Blength);

                        MPI_Sendrecv(d0Wdelta + j_Woffset + div_Woffset, div_Wlength, MPI_DOUBLE, partner, 0,
                                d1Wdelta + div_Woffset, div_Wlength, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD,
                                MPI_STATUS_IGNORE);
                        MPI_Sendrecv(d0Bdelta + j_Boffset + div_Boffset, div_Blength, MPI_DOUBLE, partner, 0,
                                d1Bdelta + div_Boffset, div_Blength, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD,
                                MPI_STATUS_IGNORE);
                        receive_flag[(k * div_reduce + div_no) * 3 + 1] = flag_target;//tell the master line
                    }
                    bitmask <<= 1;

                    while (true)//make sure that master line has finished its job
                    {
                        loop_counter++;
                        if (GetSecondElement(sum_flag + 3 * k, loop_counter) >= flag_target)
                            break;
                    }
                }
            }
            else if (nprocs > 1)
            {
                MPI_Allreduce(MPI_IN_PLACE, d0Bdelta + j_Boffset, Boffset[j + 1] - Boffset[j], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, d0Wdelta + j_Woffset, Woffset[j + 1] - Woffset[j], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            }
            total_comm_time += omp_get_wtime() - rst;

            reduce_flag[1] = j + 1 + i * numL;
        }
    }

    printf("thread %d escape from comm line\n", omp_get_thread_num());
}

void CopyBunchToDouble(NodeInArg * curr_bunch_ptr)
{
    for(int i = 0; i < numN * numD; i++)
        curr_bunch_ptr->d_X[i] = (double)curr_bunch_ptr->dd_X[i];
}

void CPUReadLine(CpuArg& cpuArg, NodeArg &nodeArg, ChunkContainer& oneChunk, NodeInArg *bunchs)
{
    //do nothing now
    printf("thread %d enter readline\n", omp_get_thread_num());

    int bunch_iter = 0, bunch_overlap = 0;
    int loop_counter = 0, flag_target = 0;
    int flag_value0;
    int flag_value1;
    int flag_value2;
    int flag_value3;
    int dead_flag_value;
    int readSize;
    int chunkCnt = 0;

    while ((readSize = FetchOneChunk(cpuArg, oneChunk)) && readSize)
    {
        int bunch_in_chunk = readSize / cpuArg.bunchSize;
        if (myid == 0)
        {
            printf("  -_- Chunk %d, samples = %d, bunchs in chunk = %d -_-  \n", chunkCnt + 1, readSize, bunch_in_chunk);
            fprintf(cpuArg.pLogFile,"--chunk(%d) : containing samples %d\n", chunkCnt++, readSize);
            fflush(cpuArg.pLogFile);
        }
        for (int i = 0; i < bunch_in_chunk; i++)
        {
            flag_target = bunch_overlap;
            while(true)
            {
                loop_counter++;
                flag_value0 = GetSecondElement(use_flag0 + bunch_iter * 3, loop_counter);
                flag_value1 = GetSecondElement(use_flag1 + bunch_iter * 3, loop_counter);
                flag_value2 = GetSecondElement(use_flag2 + bunch_iter * 3, loop_counter);
                flag_value3 = GetSecondElement(use_flag_c + bunch_iter * 3, loop_counter);
                dead_flag_value = GetSecondElement(dead_flag, loop_counter);

                if(GetSecondElement(dead_flag, loop_counter) == -2016)
                    break;

                if(flag_value0 >= flag_target && flag_value1 >= flag_target && flag_value2 >= flag_target && flag_value3 >= flag_target)
                    break;
            }
            if(GetSecondElement(dead_flag, loop_counter) == -2016)
                break;

            myFetchOneBunch(oneChunk, bunchs[bunch_iter]);

            CopyBunchToDouble(bunchs + bunch_iter);

            read_flag[1 + 3 * bunch_iter] = bunch_overlap + 1;

            bunch_iter++;

            if(bunch_iter == max_bunch_count)
            {
                bunch_iter = 0;
                bunch_overlap++;

                flag_target = bunch_overlap;

                //printf("readline has reach the end of one overlap %d, waiting for flag %d\n", bunch_overlap, flag_target);
            }
        }
        if(GetSecondElement(dead_flag, loop_counter) == -2016)
            break;
    }



    for(int i = 0; i < max_bunch_count; i++)
        read_flag[1 + i * 3] = bunch_overlap + 2;

    printf("thread %d escape from read line\n", omp_get_thread_num());
}

void CPUOffloadLine(NodeInArg *bunchs, int t_mic, int *t_signal)
{
    int* local_Wsz = W_sizes + 3 * t_mic;
    int* local_Bsz = B_sizes + 3 * t_mic;
    int* local_Tsz = T_sizes + 3 * t_mic;
    int* local_Xsz = X_sizes + 3 * t_mic;
    int* local_Ysz = Y_sizes + 3 * t_mic;
    int my_signal = t_mic;
    int MPI_numN_startpos = BLOCK_LOW(myid, nprocs, ori_numN);

    int loop_counter = 0, flag_target;
    int flag_value;

    flag_target = 1;
    while(true)
    {
        loop_counter++;
        flag_value = GetSecondElement(read_flag, loop_counter);
        if(flag_value >= flag_target)
            break;
    }

    if (t_mic == 0)
    {
        memcpy(d0T, bunchs[0].d_T + MPI_numN_startpos + start_rows[0], T_sizes[0] * sizeof(int));
        memcpy(d0X, bunchs[0].d_X + (MPI_numN_startpos + start_rows[0]) * numD, X_sizes[0] * sizeof(double));
        #pragma offload_transfer target(mic:t_mic) \
        in     (d0X        : length(local_Xsz[0])      REUSE  align(512)) \
        in     (d0T        : length(local_Tsz[0])      REUSE  align(512))
    }
    else if (t_mic == 1)
    {
        memcpy(d1T, bunchs[0].d_T + MPI_numN_startpos + start_rows[1], T_sizes[4] * sizeof(int));
        memcpy(d1X, bunchs[0].d_X + (MPI_numN_startpos + start_rows[1]) * numD, X_sizes[4] * sizeof(double));
        #pragma offload_transfer target(mic:t_mic) \
        in     (d1X        : length(local_Xsz[1])      REUSE  align(512)) \
        in     (d1T        : length(local_Tsz[1])      REUSE  align(512))
    }
    else
    {
        memcpy(d2T, bunchs[0].d_T + MPI_numN_startpos + start_rows[2], T_sizes[8] * sizeof(int));
        memcpy(d2X, bunchs[0].d_X + (MPI_numN_startpos + start_rows[2]) * numD, X_sizes[8] * sizeof(double));
        #pragma offload_transfer target(mic:t_mic) \
        in     (d2X        : length(local_Xsz[2])      REUSE  align(512)) \
        in     (d2T        : length(local_Tsz[2])      REUSE  align(512))
    }

    if(t_mic == 0)
    {
        use_flag0[1] = 1;
        if(!dnn_on_cpu)
            use_flag_c[1] = 1;
    }
    if(t_mic == 1)
    {
        use_flag1[1] = 1;
    }
    if(t_mic == 2)
    {
        use_flag2[1] = 1;
    }

//?
    #pragma offload target(mic:t_mic) signal(t_signal) \
    nocopy (d0X          : length(local_Xsz[0])      REUSE  align(512)) \
    nocopy (d0T          : length(local_Tsz[0])      REUSE  align(512)) \
    nocopy (d1X          : length(local_Xsz[1])      REUSE  align(512)) \
    nocopy (d1T          : length(local_Tsz[1])      REUSE  align(512)) \
    nocopy (d2X          : length(local_Xsz[2])      REUSE  align(512)) \
    nocopy (d2T          : length(local_Tsz[2])      REUSE  align(512)) \
    nocopy (d0Wdelta     : length(local_Wsz[0])      REUSE  align(512)) \
    nocopy (d0Bdelta     : length(local_Bsz[0])      REUSE  align(512)) \
    nocopy (d1Wdelta     : length(local_Wsz[1])      REUSE  align(512)) \
    nocopy (d1Bdelta     : length(local_Bsz[1])      REUSE  align(512)) \
    nocopy (d2Wdelta     : length(local_Wsz[2])      REUSE  align(512)) \
    nocopy (d2Bdelta     : length(local_Bsz[2])      REUSE  align(512)) \
    nocopy (hW           : length(local_Wsz[t_mic])      REUSE  align(512)) \
    nocopy (hB           : length(local_Bsz[t_mic])      REUSE  align(512)) \
    nocopy (dY           : length(local_Ysz[t_mic])      REUSE  align(512)) \
    nocopy (dE           : length(local_Ysz[t_mic])      REUSE  align(512)) \
    nocopy (mic_identifer    : length(3)    REUSE) \
    nocopy (X_sizes      : length(9)    REUSE) \
    nocopy (T_sizes      : length(9)    REUSE) \
    nocopy (W_sizes      : length(9)    REUSE) \
    nocopy (B_sizes      : length(9)    REUSE) \
    nocopy (Y_sizes      : length(9)    REUSE) \
    nocopy (row_counts       : length(3)    REUSE) \
    nocopy (start_rows       : length(3)    REUSE) \
    nocopy (Woffset      : length(numL + 1) REUSE) \
    nocopy (Boffset      : length(numL + 1) REUSE) \
    nocopy (numA         : length(numL)     REUSE) \
    nocopy (mic_Yoffset0     : length(numL + 1) REUSE) \
    nocopy (mic_Yoffset1     : length(numL + 1) REUSE) \
    nocopy (mic_Yoffset2     : length(numL + 1) REUSE) \
    nocopy (get_delta_flag0  : length(3)    REUSE) \
    nocopy (get_delta_flag1  : length(3)    REUSE) \
    nocopy (get_delta_flag2  : length(3)    REUSE) \
    nocopy (forward_flag0    : length(3)    REUSE) \
    nocopy (forward_flag1    : length(3)    REUSE) \
    nocopy (forward_flag2    : length(3)    REUSE) \
    nocopy (backward_flag0   : length(3)    REUSE) \
    nocopy (backward_flag1   : length(3)    REUSE) \
    nocopy (backward_flag2   : length(3)    REUSE) \
    nocopy (in_flag      : length(9)    REUSE) \
    nocopy (mic_current_bunch_ : length(9)      REUSE) \
    nocopy (mic_current_layer_ : length(9)      REUSE) \
    in(numL) in(numD) in(one) in(zero) in(alpha) in(bunch_count) in(learnrate_changerate) in(lr_destination)
    {
        int my_micid = GetSecondElement(mic_identifer, numL);
        int *my_T, *my_Yoffset, *my_delta_flag, *my_forward_flag, *my_backward_flag;
        double *my_X, *my_Wdelta, *my_Bdelta;

        if (my_micid == 0)
        {
            my_T = d0T;
            my_X = d0X;
            my_Yoffset = mic_Yoffset0;
            my_Wdelta = d0Wdelta;
            my_Bdelta = d0Bdelta;
            my_delta_flag = get_delta_flag0;
            my_forward_flag = forward_flag0;
            my_backward_flag = backward_flag0;
        }
        else if (my_micid == 1)
        {
            my_T = d1T;
            my_X = d1X;
            my_Yoffset = mic_Yoffset1;
            my_Wdelta = d1Wdelta;
            my_Bdelta = d1Bdelta;
            my_delta_flag = get_delta_flag1;
            my_forward_flag = forward_flag1;
            my_backward_flag = backward_flag1;
        }
        else
        {
            my_T = d2T;
            my_X = d2X;
            my_Yoffset = mic_Yoffset2;
            my_Wdelta = d2Wdelta;
            my_Bdelta = d2Bdelta;
            my_delta_flag = get_delta_flag2;
            my_forward_flag = forward_flag2;
            my_backward_flag = backward_flag2;
        }

        for (int i = 0; i < bunch_count; i++)
        {
            mic_current_bunch_[1 + 3 * my_micid] = i;

            mic_dnnForward(
                        numL, row_counts[my_micid], numD, numA,
                        my_X, dY, hW, hB,
                        my_Yoffset, Boffset, Woffset,
                        one, my_micid, my_Wdelta, my_Bdelta, my_forward_flag, mic_current_bunch_
                        );

            mic_dnnBackward(
                numL, row_counts[my_micid], numA,
                dE, dY, my_T, hW,
                my_Yoffset, my_Yoffset, Woffset,
                one, zero, my_micid, my_backward_flag, mic_current_bunch_
            );

            mic_dnnGetDelta(
                numL, numD, row_counts[my_micid], numA,
                my_X, dE, dY, my_Wdelta, my_Bdelta,
                my_Yoffset, my_Yoffset, Woffset,
                zero, alpha, my_micid, my_delta_flag, mic_current_bunch_
            );

            if (lr_destination < 0.98f)
                alpha *= learnrate_changerate;
        }
    }

    printf("Node %d  escape from offload line %d\n", myid, t_mic);
}

void CPUTransferLine(void* t_mic_ptr)//t_mic: the specific mic number
{
    int t_mic = *((int*)t_mic_ptr);
    NodeInArg *bunchs = bunchs_ptr;
    printf("transfer line of mic %d is ready.\n", t_mic);
    int in_layer = 0, in_bunch = 0, out_layer = 0, out_bunch = 0;
    trans_variances_[4 * t_mic] = &out_bunch;
    trans_variances_[4 * t_mic + 1] = &out_layer;
    trans_variances_[4 * t_mic + 2] = &in_bunch;
    trans_variances_[4 * t_mic + 3] = &in_layer;
    bool new_X_transfered = false;
    bool new_T_transfered = false;
    int phase = 666666;
    int flag_target, flag_value;
    int bunch_flag_value, bunch_flag_target;
    int bunch_iter, bunch_overlap;
    int loop_count = 0;
    bool last_time_do_job = true;

    double wait_start_time = omp_get_wtime();//for break out, reset when do any jobs

    while(true)
    {
        transfer_stage_[t_mic] = 1;
        loop_count++;
        if (in_bunch >= bunch_count - 1 && out_bunch >= bunch_count)
            break;
        if (omp_get_wtime() - wait_start_time > 10.0)
            CPUPrintErrorMessage();
        if (phase % 2 == 0)//out
        {
            transfer_stage_[t_mic] = 2;
            if (!last_time_do_job)
                ConsumeSomeTime();
            flag_target = out_layer + 1 + numL * out_bunch;
            phase++;
            transfer_stage_[t_mic] = 3;
            flag_value = TransferOutGetDeltaFlag(t_mic);
            if (flag_value >= flag_target)
            {
                TransferOutDelta(out_layer, t_mic, out_bunch);
                out_flag[1 + 3 * t_mic] = flag_target;

                //printf("node %d: get delta bunch = %d, layer = %d\n", t_mic, out_bunch, out_layer);

                if (out_layer == numL - 1)
                {
                    out_layer = 0;
                    out_bunch++;
                }
                else
                {
                    out_layer++;
                }
                wait_start_time = omp_get_wtime();
            }
        }
        else//in
        {
            transfer_stage_[t_mic] = 4;
            last_time_do_job = false;
            flag_target = in_layer + 1 + numL * in_bunch;
            phase--;
            if (in_bunch == bunch_count - 1)
                continue;
            if (!new_X_transfered)//try to transfer in X T after the first layer forward
            {
                transfer_stage_[t_mic] = 5;
                flag_value = TransferOutForwardFlag(t_mic);
                bunch_iter = (in_bunch + 1) % max_bunch_count;
                bunch_overlap = (in_bunch + 1) / max_bunch_count;
                bunch_flag_target = bunch_overlap + 1;
                bunch_flag_value = GetSecondElement(read_flag + 3 * bunch_iter, loop_count);
                if (flag_value > in_bunch && bunch_flag_value >= bunch_flag_target)
                {
                    TransferInX(t_mic, bunchs, bunch_iter);
                    wait_start_time = omp_get_wtime();
                    last_time_do_job = true;
                    new_X_transfered = true;
                    //printf("node %d transfer new X of bunch %d\n", myid, in_bunch);
                }
            }
            else if (!new_T_transfered)
            {
                transfer_stage_[t_mic] = 6;
                flag_value = TransferOutBackwardFlag(t_mic);
                if (flag_value > in_bunch)
                {
                    TransferInT(t_mic, bunchs, bunch_iter);
                    wait_start_time = omp_get_wtime();
                    last_time_do_job = true;
                    new_T_transfered = true;

                    if(t_mic == 0)
                    {
                        use_flag0[1 + 3 * bunch_iter] = 1 + bunch_overlap;
                        if(!dnn_on_cpu)
                            use_flag_c[1 + 3 * bunch_iter] = 1 + bunch_overlap;
                    }
                    if(t_mic == 1)
                    {
                        use_flag1[1 + 3 * bunch_iter] = 1 + bunch_overlap;
                    }
                    if(t_mic == 2)
                    {
                        use_flag2[1 + 3 * bunch_iter] = 1 + bunch_overlap;
                    }

                    //printf("node %d transfer new T of bunch %d\n", myid, in_bunch);
                }
            }
            else
            {
                transfer_stage_[t_mic] = 7;
                flag_value = GetSecondElement(copy_flag, loop_count);
                if (flag_value >= flag_target)//then add layer
                {
                    TransferInDelta(in_bunch, in_layer, t_mic);
                    //printf("node %d: in new delta bunch = %d, layer = %d\n", t_mic, in_bunch, in_layer);
                    if (in_layer == numL - 1)
                    {
                        in_layer = 0;
                        in_bunch++;
                        new_X_transfered = false;
                        new_T_transfered = false;
                    }
                    else
                    {
                        in_layer++;
                    }
                    last_time_do_job = true;
                    wait_start_time = omp_get_wtime();
                }
            }
        }
    }

    if(t_mic == 0)
        dead_flag[1] = -2016;

    printf("thread %d escape from transfer line\n", omp_get_thread_num());
}

void CPUMasterLine()
{
    //printf("enter master line with bunch_count = %d\n", bunch_count);

    double about_total_time;
    double about_time_left;
    int thread_id = 0;
    int i, j, k, l, div_no;
    master_i_ptr_ = &i;
    master_j_ptr_ = &j;
    int flag_target;
    int flag_value0, flag_value1, flag_value2;
    int loop_counter = 0;
    int calculated_statu;//1:0+1  2:1+2  3:0+2
    bool round_finished;
    int j_Woffset, j_Boffset;
    int j_Wlength, j_Blength;
    int div_Woffset, div_Boffset;
    int div_Wlength, div_Blength;

    int print_stride = 50000 / numN + 1;

    double st, et;

    for (i = 0; i < bunch_count; i++)
    {
        st = omp_get_wtime();
        for (j = 0; j < numL; j++)
        {
            j_Woffset = Woffset[j];
            j_Wlength = Woffset[j + 1] - Woffset[j];
            j_Boffset = Boffset[j];
            j_Blength = Boffset[j + 1] - Boffset[j];

            master_stage_ = 1;
            flag_target = j + i * numL + 1;
            //before writing the inner_reduce_flag, I need to set the sendrecv
            for (k = 0; k < 3 * div_reduce * sendrecv_loop; k++)
            {
                receive_flag[k] = 0;
            }

            while (true)//waiting for out 1
            {
                loop_counter++;
                flag_value0 = GetSecondElement(out_flag, loop_counter);
                flag_value1 = GetSecondElement(out_flag + 3, loop_counter);
                flag_value2 = GetSecondElement(out_flag + 6, loop_counter);

                if (flag_value0 > j + i * numL && flag_value1 > j + i * numL)
                {
                    dst_Wdelta = d0Wdelta;
                    dst_Bdelta = d0Bdelta;
                    src_Wdelta = d1Wdelta;
                    src_Bdelta = d1Bdelta;
                    calculated_statu = 1;
                    break;
                }
                if (flag_value2 > j + i * numL && flag_value1 > j + i * numL)
                {
                    dst_Wdelta = d1Wdelta;
                    dst_Bdelta = d1Bdelta;
                    src_Wdelta = d2Wdelta;
                    src_Bdelta = d2Bdelta;
                    calculated_statu = 2;
                    break;
                }
                if (flag_value0 > j + i * numL && flag_value2 > j + i * numL)
                {
                    dst_Wdelta = d0Wdelta;
                    dst_Bdelta = d0Bdelta;
                    src_Wdelta = d2Wdelta;
                    src_Bdelta = d2Bdelta;
                    calculated_statu = 3;
                    break;
                }
            }//if any 2 delta of 3 mic is outed, begin the first step of inner reduce

            //printf("master, bunch = %d, layer = %d, assign mission 1\n", i, j);

            for (k = 0; k < worker_count + 1; k++)//this should be reseted on waiters
                assign_flag[1 + k * 3] = j + i * numL + 1;
            CPUSumMission(dst_Wdelta + j_Woffset, src_Wdelta + j_Woffset, j_Wlength, worker_count + 1, thread_id);
            CPUSumMission(dst_Bdelta + j_Boffset, src_Bdelta + j_Boffset, j_Blength, worker_count + 1, thread_id);
            finish_flag[1] = j + i * numL + 1;
            //waiting for the round 1 finished
            master_stage_ = 2;

            while (true)
            {
                loop_counter++;
                round_finished = true;
                for (int k = 0; k < worker_count + 1; k++)
                {
                    if (GetSecondElement(finish_flag + k * 3, loop_counter) <= j + i * numL)
                    {
                        round_finished = false;
                        break;
                    }
                }
                if (round_finished)
                    break;
            }
            //printf("master, bunch = %d, layer = %d, report for finishing mission 1\n", i, j);
            for (k = 0; k < (worker_count + 1) * 3; k++)
                finish_flag[k] = j + i * numL;
            //printf("statu = %d, end\n", calculated_statu);
            //CPUCheckDelta(src_Bdelta + j_Boffset, 666, "src_Bdelta", j, i);
            //CPUCheckDelta(dst_Bdelta + j_Boffset, 666, "dst_Bdelta", j, i);

            master_stage_ = 3;
            while (true)//waiting for out 2
            {
                loop_counter++;
                flag_value0 = GetSecondElement(out_flag, loop_counter);
                flag_value1 = GetSecondElement(out_flag + 3, loop_counter);
                flag_value2 = GetSecondElement(out_flag + 6, loop_counter);

                if (flag_value0 > j + i * numL && flag_value1 > j + i * numL && flag_value2 > j + i * numL)
                {
                    if (calculated_statu == 1)
                    {
                        dst_Wdelta = d0Wdelta;
                        dst_Bdelta = d0Bdelta;
                        src_Wdelta = d2Wdelta;
                        src_Bdelta = d2Bdelta;
                    }
                    else
                    {
                        dst_Wdelta = d0Wdelta;
                        dst_Bdelta = d0Bdelta;
                        src_Wdelta = d1Wdelta;
                        src_Bdelta = d1Bdelta;
                    }
                    break;
                }
            }
            //printf("master, bunch = %d, layer = %d, assign mission 2\n", i, j);
            for (k = 0; k < (worker_count + 1) * 3; k++)//this should be reseted on waiters
                assign_flag[k] = j + i * numL + 1;
            CPUSumMission(dst_Wdelta + j_Woffset, src_Wdelta + j_Woffset, j_Wlength, worker_count + 1, thread_id);
            CPUSumMission(dst_Bdelta + j_Boffset, src_Bdelta + j_Boffset, j_Blength, worker_count + 1, thread_id);
            finish_flag[1] = j + i * numL + 1;

            master_stage_ = 4;
            while (true)//the second inner reduce mission
            {
                loop_counter++;
                round_finished = true;
                for (k = 0; k < worker_count + 1; k++)
                {
                    if (GetSecondElement(finish_flag + k * 3, loop_counter) <= j + i * numL)
                    {
                        round_finished = false;
                        break;
                    }
                }
                if (round_finished)
                    break;
            }//after this, we can do the reduce mission
            //printf("master, bunch = %d, layer = %d, report for finishing mission 2\n", i, j);
            for (k = 0; k < (worker_count + 1) * 3; k++)
                finish_flag[k] = j + i * numL;

            if(dnn_on_cpu)
            {
                master_stage_ = 5;
                while (true)//waiting for cpu
                {
                    loop_counter++;
                    flag_value0 = GetSecondElement(get_delta_flag_c, loop_counter);

                    if (flag_value0 > j + i * numL)
                    {
                        dst_Wdelta = d0Wdelta;
                        dst_Bdelta = d0Bdelta;
                        src_Wdelta = hWdelta;
                        src_Bdelta = hBdelta;
                        break;
                    }
                }
                //printf("master, bunch = %d, layer = %d, assign mission 3\n", i, j);
                for (k = 0; k < (worker_count + 1) * 3; k++)
                    assign_flag[k] = j + i * numL + 1;
                CPUSumMission(dst_Wdelta + j_Woffset, src_Wdelta + j_Woffset, j_Wlength, worker_count + 1, thread_id);
                CPUSumMission(dst_Bdelta + j_Boffset, src_Bdelta + j_Boffset, j_Blength, worker_count + 1, thread_id);
                finish_flag[1] = j + i * numL + 1;

                master_stage_ = 6;
                while (true)
                {
                    loop_counter++;
                    round_finished = true;
                    for (k = 0; k < worker_count + 1; k++)
                    {
                        if (GetSecondElement(finish_flag + k * 3, loop_counter) <= j + i * numL)
                        {
                            round_finished = false;
                            break;
                        }
                    }
                    if (round_finished)
                        break;
                }
                //printf("master, bunch = %d, layer = %d, report for finishing mission 3\n", i, j);
                for (k = 0; k < (worker_count + 1) * 3; k++)
                    finish_flag[k] = j + i * numL;
            }

            inner_reduce_flag[1] = j + i * numL + 1;

            master_stage_ = 7;
            if (max_2_pow != nprocs)
                while (true)//wait for MPI_Allreduce
                {
                    flag_value1 = GetSecondElement(reduce_flag, loop_counter);

                    if (flag_value1 > j + i * numL)
                        break;
                }
            else if (nprocs > 1)//sendrecv mission
            {
                for (k = 0; k < sendrecv_loop; k++)
                {
                    for (div_no = 0; div_no < div_reduce; ++div_no)
                    {
                        loop_counter = 0;
                        while (true)
                        {
                            loop_counter++;
                            if (GetSecondElement(receive_flag + (k * div_reduce + div_no) * 3, loop_counter) >= flag_target)//the receive finish flag
                            {
                                break;
                            }
                        }

                        for (l = 0; l < (worker_count + 1) * 3; l++)//this should be reseted on waiters
                            assign_flag[l] = j + i * numL + 1;

                        div_Wlength = BLOCK_SIZE(div_no, div_reduce, j_Wlength);
                        div_Woffset = BLOCK_LOW(div_no, div_reduce, j_Wlength);
                        div_Blength = BLOCK_SIZE(div_no, div_reduce, j_Blength);
                        div_Boffset = BLOCK_LOW(div_no, div_reduce, j_Blength);

                        CPUSumMission(d0Wdelta + j_Woffset + div_Woffset, d1Wdelta + div_Woffset, div_Wlength, worker_count + 1, thread_id);
                        CPUSumMission(d0Bdelta + j_Boffset + div_Boffset, d1Bdelta + div_Boffset, div_Blength, worker_count + 1, thread_id);
                        finish_flag[1] = j + i * numL + 1;

                        while (true)//waiting for the finish of outer reduce
                        {
                            loop_counter++;
                            round_finished = true;
                            for (l = 0; l < worker_count + 1; l++)
                            {
                                if (GetSecondElement(finish_flag + l * 3, loop_counter) <= j + i * numL)
                                {
                                    round_finished = false;
                                    break;
                                }
                            }
                            if (round_finished)
                                break;
                        }

                        for (l = 0; l < (worker_count + 1) * 3; l++)
                            finish_flag[l] = j + i * numL;
                    }

                    sum_flag[1 + 3 * k] = flag_target;
                }
            }

            for (k = 0; k < (worker_count + 1) * 3; k++)
                assign_flag[k] = j + i * numL + 1;
            CPUSumMission(hW + j_Woffset, d0Wdelta + j_Woffset, j_Wlength, worker_count + 1, thread_id);
            CPUSumMission(hB + j_Boffset, d0Bdelta + j_Boffset, j_Blength, worker_count + 1, thread_id);
            finish_flag[1] = j + i * numL + 1;

            master_stage_ = 7;
            while (true)
            {
                loop_counter++;
                round_finished = true;
                for (k = 0; k < worker_count + 1; k++)
                {
                    if (GetSecondElement(finish_flag + k * 3, loop_counter) <= j + i * numL)
                    {
                        round_finished = false;
                        break;
                    }
                }
                if (round_finished)
                    break;
            }
            for (k = 0; k < (worker_count + 1) * 3; k++)
                finish_flag[k] = j + i * numL;

            copy_flag[1] = j + i * numL + 1;
        }

        et = omp_get_wtime();
        total_time += et - st;
        if (myid == 0 && i % print_stride == 0)
        {
            CPUCheck2Arr(hW, hB, "3 hW and 3 hB   ", 0, i);
            about_total_time = bunch_count * (total_time / (i + 1));
            about_time_left = about_total_time - total_time;
            printf("train bunch %d, time = %.4f, total time = %.2f, up to now average = %.4f, about %.2f(s) left\n", i, et - st, total_time, total_time / (i + 1), about_time_left);
        }
    }

    printf("thread %d escape from master line\n", omp_get_thread_num());
}

void PureMICTrain(NodeInArg *bunchs, int total_bunchs, CpuArg& cpuArg, NodeArg &nodeArg, ChunkContainer& oneChunk)
{
    bunch_count = total_bunchs;
    if (myid == 0)
    {
        printf("Data initial on MIC is ready.\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int int_values[24];
    for (int i = 0; i < 24; i++)
        int_values[i] = i;

    int thread_count = 6 + worker_count + dnn_thread_count;

    printf("okokok, thread_count = %d\n", thread_count);

    if(!dnn_on_cpu)
    {
        //thread_count -= dnn_thread_count;//then no threads will enter dnnline
        //worker_count += dnn_thread_count;
    }

#pragma omp parallel num_threads(thread_count)
    {
        int my_tid = omp_get_thread_num();

        if(my_tid == 0)
            CPUOffloadLine(bunchs, 0, &signal0);
        if(my_tid == 1)
            CPUOffloadLine(bunchs, 1, &signal1);
        if(my_tid == 2)
            CPUOffloadLine(bunchs, 2, &signal2);

        if(my_tid == 0)
            CPUCommunicateLine();
        else if(my_tid >= 1 && my_tid <= 3)
            CPUTransferLine((void*)(int_values + my_tid - 1));
        else if(my_tid == 4)
            CPUMasterLine();
        else if(my_tid <= 5 + worker_count - 1)
            CPUWorkerLine((void*)(int_values + my_tid - 4));
        else if(my_tid == 5 + worker_count)
            CPUReadLine(cpuArg, nodeArg, oneChunk, bunchs);
        else if(dnn_on_cpu)
            CPUDnnLine();
        else
            CPUWorkerLine((void*)(int_values + my_tid - 5));
    }

    #pragma offload_wait target(mic:0) wait(&signal0)
    #pragma offload_wait target(mic:1) wait(&signal1)
    #pragma offload_wait target(mic:2) wait(&signal2)
    MICDataOut();

    //mic_updateArr(hW, d0Wdelta, Woffset[numL]);
    //mic_updateArr(hB, d0Bdelta, Boffset[numL]);

    mic_updateArr(hW, d0Wdelta, Woffset[numL]);
    mic_updateArr(hB, d0Bdelta, Boffset[numL]);

    printf("all time = %lf(s), average all time = %lf(s), average Allreduce time = %lf(s)\n", total_time, total_time / total_bunchs, total_comm_time / total_bunchs);
    fflush(stdout);
}

void dnnTrain(CpuArg& cpuArg, NodeArg &nodeArg, NodeInArg *bunchs, int total_bunchs, int _myid, int _nprocs, ChunkContainer &oneChunk)
{
    CPUCheckDelta(hW, 666, "initial hW", 0, 0);
    CPUCheckDelta(hB, 666, "initial hB", 0, 0);
    PureMICTrain(bunchs, total_bunchs, cpuArg, nodeArg, oneChunk);
    CPUCheckDelta(hW, 666, "last hW", 0, 0);
    CPUCheckDelta(hB, 666, "last hB", 0, 0);
    copyResultBack(nodeArg);
}

void GetConfig(int* faster_readfile, int total_bunchs, int _myid, int _nprocs, NodeArg &nodeArg)//read myconfig.txt
{
    bunch_count = total_bunchs;
    myid = _myid;
    nprocs = _nprocs;

    numN = nodeArg.numN;          //size of minibatch
    numL = nodeArg.dnnLayerNum - 1;   //layer nums
    numD = nodeArg.dnnLayerArr[0];    //node nums of input layer
    numA = &(nodeArg.dnnLayerArr[1]); //node nums of hiden layer and output layer
    alpha = -nodeArg.lRate / numN;
    ori_alpha = alpha;
    ori_numN = nodeArg.numN;

    max_bunch_count = 1500000 / numN;

    if (myid == 0)
    {
        FILE *pInitFile = fopen("myconfig.txt", "r");

        if (NULL == pInitFile)
        {
            printf("ERROR: Failed to open config file %s\n", "myconfig.txt");
            //fflush(stdout);
            lr_destination = 1.0;
            div_reduce = 4;
            cpu_mic_div = 0.91;
            worker_count = 3;
        }
        else
        {
            char cfgline[MAXLINE];
            char *p = NULL;

            while(fgets(cfgline,MAXLINE,pInitFile))
            {
                p = strstr(cfgline,"=");
                if(NULL == p)
                    continue;
                *p='\0';
                p += 1;
                isLrtrim_blank(cfgline);
                isLrtrim_blank(p);

                if(0 == strcmp(cfgline,"lr_destination"))
                {
                    float tmpfloat;
                    tmpfloat = atof(p);
                    lr_destination = tmpfloat;
                }
                if(0 == strcmp(cfgline,"div_reduce"))
                    div_reduce = atoi(p);
                if(0 == strcmp(cfgline,"worker_count"))
                    worker_count = atoi(p);
                if(0 == strcmp(cfgline,"cpu_mic_div"))
                    cpu_mic_div = atof(p);
                if(0 == strcmp(cfgline,"fast_readfile"))
                    faster_readfile[0] = atoi(p);
            }
        }

        //in case of obvious wrong or empty config
        if (!(lr_destination > 0.001 && lr_destination < 0.9)) lr_destination = 1.0;
        if (div_reduce < 0 || div_reduce > 10)
            div_reduce = 4;
        if(cpu_mic_div < 0.84)
            cpu_mic_div = 0.91;
        if(cpu_mic_div >= 1)
            cpu_mic_div = 100;
    }

    worker_count = 2;

    MPI_Bcast(&lr_destination, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cpu_mic_div, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&div_reduce, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(faster_readfile, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&worker_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(cpu_mic_div >= 1)
    {
        numN = BLOCK_SIZE(myid, nprocs, numN);
        mic_numN = numN;
        cpu_numN = 1;
        dnn_on_cpu = false;
    }
    else
    {
        numN = BLOCK_SIZE(myid, nprocs, numN);
        mic_numN = numN * cpu_mic_div;
        cpu_numN = numN - mic_numN;
        dnn_on_cpu = true;
    }



    dnn_thread_count = 24 - (worker_count + 6);
    dnn_start_thread = worker_count + 6;

    if(!dnn_on_cpu)
    {
        //thread_count -= dnn_thread_count;//then no threads will enter dnnline
        worker_count += dnn_thread_count;
        dnn_thread_count = 0;
    }

    GetLCrate(bunch_count);
    printf("\nconfig on node %d: lr_destination = %.8f , ori_alpha = %.8f, alpha = %.8f, learnrate_change_rate = %.8f, cpu_mic_div = %.4f, div_reduce = %d, faster_readfile = %d, dnn_thread_count = %d, dnn_start_thread = %d, max_bunch_count = %d, total_bunchs = %d\n\n",
           myid, lr_destination, ori_alpha, alpha, learnrate_changerate, cpu_mic_div, div_reduce, faster_readfile, dnn_thread_count, dnn_start_thread, max_bunch_count, total_bunchs
           );

    printf("MPI %d : %d samples, starts from %d, mic: %d, cpu: %d\n", myid, numN, BLOCK_LOW(myid, nprocs, nodeArg.numN), mic_numN, cpu_numN);


    max_2_pow = 1;
    sendrecv_loop = 0;
    while (max_2_pow < nprocs)
    {
        sendrecv_loop += 1;
        max_2_pow *= 2;
    }//if max_2_pow == nprocs, then do sendrecv
}

void* KernelInitial(void* v_paras)
{
    KernelPara* paras = (KernelPara*)v_paras;

    //CpuArg& cpuArg, NodeArg &nodeArg, NodeInArg *bunchs, int total_bunchs, int _myid, int _nprocs;

    //myid = paras->_myid;
    //nprocs = paras->_nprocs;
    MPI_Status status;
    bunchs_ptr = paras->bunchsp;
    //bunch_count = paras->total_bunchs;

    copyDataToLocal((*paras->cpuArgp), (*paras->nodeArgp));

    CPUArrayInitial();
    MICDataInitial();

    pthread_exit(NULL);
}

void* KernelUninitial(void *no_use_ptr)
{
    MICDataUninitial();
    CPUArrayUninitial();

    pthread_exit(NULL);
}

void OldMain(int _myid, int _nprocs, int argc, char* argv[], CpuArg &cpuArg, NodeArg &nodeArg , ChunkContainer& oneChunk)
{
    int myid = _myid, nprocs = _nprocs;
    int faster_readfile;//if this == 2016, then launch FetchOneChunkPara

    const char* initNum = NULL;
    if(2 == argc)
        initNum = "0";
    else
        initNum = argv[2];

    // Open output file IFF myid == 0
    int ret = GetInitFileConfig(argv[1],initNum,cpuArg, myid);

    if(0 != ret) {
        printf("Error happens: GetInitFileCofig\n");
        exit(ret);
    }

    InitNodeConfig(cpuArg,nodeArg);

    //training process
    fprintf(cpuArg.pLogFile,"start training:\n");
    fflush(cpuArg.pLogFile);
    int  chunkCnt  = 0;
    int  readSize  = 0;
    int  countCnt = 0;

    int s_sample = BLOCK_LOW(myid    , nprocs, nodeArg.numN);
    int e_sample = BLOCK_LOW(myid + 1, nprocs, nodeArg.numN);

    double input_st = MPI_Wtime();

    int total_bunchs = cpuArg.totalSamples / cpuArg.bunchSize;

    NodeInArg* bunchs = (NodeInArg*) malloc(sizeof(NodeInArg) * total_bunchs);


    GetConfig(&faster_readfile, total_bunchs, myid, nprocs, nodeArg);

    pthread_t initial_thread;

    KernelPara kernel_para;
    kernel_para._myid = myid;
    kernel_para._nprocs = nprocs;
    kernel_para.total_bunchs = total_bunchs;
    kernel_para.cpuArgp = &cpuArg;
    kernel_para.nodeArgp = &nodeArg;
    kernel_para.bunchsp = bunchs;

    //total_bunchs = 0;

    //KernelInitial((void*)(&kernel_para));
    pthread_create(&initial_thread, NULL, KernelInitial, (void*)(&kernel_para));

    /*while ((readSize = FetchOneChunk(cpuArg, oneChunk)) && readSize)
    {

        int bunch_in_chunk = readSize / cpuArg.bunchSize;
        if (myid == 0)
        {
            printf("Chunk %d, samples = %d, bunchs in chunk = %d\n", chunkCnt + 1, readSize, bunch_in_chunk);
            fprintf(cpuArg.pLogFile,"--chunk(%d) : containing samples %d\n", chunkCnt++, readSize);
            fflush(cpuArg.pLogFile);
        }
        for (int i = 0; i < bunch_in_chunk; i++)
        {
            myInitNodeMem(cpuArg, bunchs[total_bunchs]);
            myFetchOneBunch(oneChunk, bunchs[total_bunchs++]);
        }
    }*/

    for(int i = 0; i < max_bunch_count; i++)
        myInitNodeMem(cpuArg, bunchs[i]);

    pthread_join(initial_thread, NULL);

    double input_et = MPI_Wtime();
    printf("MPI %d read in time = %lf(s)\n", myid, input_et - input_st);

    if (myid == 0)
    {
        printf("Read workload over : %d chunks, %d samples, %d bunch in total\n",
                cpuArg.totalChunks, cpuArg.totalSamples, total_bunchs + 1);
    }

    dnnTrain(cpuArg, nodeArg, bunchs, total_bunchs, myid, nprocs, oneChunk);

    pthread_t uninitial_thread;

    pthread_create(&uninitial_thread, NULL, KernelUninitial, (void*)(NULL));

    //for (int i = 0; i < max_bunch_count; i++)
    //    myFreeNodeMem(bunchs[i]);
    //free(bunchs);

    pthread_join(uninitial_thread, NULL);
}
