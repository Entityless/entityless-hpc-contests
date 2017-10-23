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

#include "kernel_util.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"

void CPUsetmatY(double *Y, double *B, int row, int col)
{
    for (int i = 0; i < row; i++)
        memcpy(Y + i * col, B, sizeof(double) * col);
}



void CPUsoftmaxZ(double* in_vec, double* out_vec, int row, int col)
{
    int base;
    double max, tmp;
    double sumexp = 0.0f;
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

void CPUerrorTrans(double *E, double *Y, int row, int col, int ldy)
{
    double tmp;
    int idx;
    for (int i = 0; i < row; i++)
    {
        idx = i * ldy;
        for (int j = 0; j < col; j++)
        {
            tmp = Y[idx + j];
            E[idx + j] = E[idx + j] * tmp * (1 - tmp);
        }
    }
}

void CPUerrorOutput(double *E, double *Z, int *T, int row, int col)
{
    double tmp;
    int idx;
    for(int i = 0; i < row; i++)
    {
        idx = i * col;
        memcpy(E + idx, Z + idx, sizeof(double) * col);
        E[idx + T[i]] -= 1.0f;
    }
}

void CPUgetBdelta(double *E, double *Bdelta, int row,
                  int col, double alpha, int scol, int ecol)
{
    for (int i = scol; i < ecol; i++) Bdelta[i] = 0.0f;
    for (int j = 0; j < row; j++)
    {
        int idx = j * col;
        for (int i = scol; i < ecol; i++) Bdelta[i] += E[idx + i];
    }
    for (int i = scol; i < ecol; i++) Bdelta[i] *= alpha;
}
