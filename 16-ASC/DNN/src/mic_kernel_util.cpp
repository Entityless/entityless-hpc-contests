/*------------------------------
@copyright (c) 2016 SYSU ASC16 Team

Last modified : 2016-03
Author(s) : Huang Chenghuan(peasharminds@163.com)

This file is part of DNNTK optimization codes for ASC16.
You can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License,
or (at your option) any later version.
------------------------------*/

#include "mic_kernel_util.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

__attribute__((target(mic)))
void mic_setmatY(double *Y, double *B, int row, int col)  //expand vector B to matrix Y
{
#pragma omp parallel for
    for (int i = 0; i < row; i++)
        memcpy(Y + i * col, B, sizeof(double) * col);
}

__attribute__((target(mic)))
void mic_softmaxZ(double* in_vec, double* out_vec, int row, int col) //handle output layer
{
    int base;
    double max, tmp;
    double sumexp = 0.0f;
#pragma omp parallel for private(base,max,tmp,sumexp)
    for (int i = 0; i < row; i++)
    {
        base = i*col;
        max = in_vec[base];
        for (int j = 1; j < col; j++)
            if (in_vec[base + j] >max)
                max = in_vec[base + j];

        sumexp = 0.0f;
        for (int j = 0; j < col; j++)
        {
            in_vec[base + j] = expf(in_vec[base + j] - max);
            sumexp += in_vec[base + j];
        }
        tmp = 1.0f / sumexp;

        for (int j = 0; j < col; j++)
            out_vec[base + j] = in_vec[base + j] * tmp;
    }
}

__attribute__((target(mic)))
void mic_errorTrans(double *E, double *Y, int row, int col)
{
    double tmp;
    int idx;
#pragma omp parallel for private(idx,tmp)
    for (int i = 0; i < row; i++)
    {
        idx = i * col;
        for (int j = 0; j < col; j++)
        {
            tmp = Y[idx + j];
            E[idx + j] = E[idx + j] * tmp * (1 - tmp);
        }
    }
}

__attribute__((target(mic)))
void mic_errorOutput(double *E, double *Z, int *T, int row, int col)
{
    int idx;
#pragma omp parallel for private(idx)
    for(int i = 0; i < row; i++)
    {
        idx = i * col;
        memcpy(E + idx, Z + idx, sizeof(double) * col);
        E[idx + T[i]] -= 1.0f;
    }
}

void mic_getBdelta(double *E, double *Bdelta, int row, int col, double alpha)
{
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();

        int my_col_count = (col / nthreads) + ((tid < col % nthreads) ? 1 : 0);
        int my_start_col = (col / nthreads) * tid + ((tid < col % nthreads) ? tid : col % nthreads);

        int i, j;

        double* my_Bdelta = Bdelta + my_start_col;
        double* my_E = E + my_start_col;

        for(i = 0; i < my_col_count; i++)
            my_Bdelta[i] = 0;

        for(i = 0; i < row; i++)
        {
            for(j = 0; j < my_col_count; j++)
                my_Bdelta[j] += my_E[j];
            my_E = my_E + col;
        }

        for(i = 0; i < my_col_count; i++)
            my_Bdelta[i] *= alpha;
    }
}

__attribute__((target(mic)))
void mic_updateArr(double *W, double *Wdelta, int length)
{
#pragma omp parallel for
    for(int i = 0; i < length; i++)
        W[i] += Wdelta[i];
}
