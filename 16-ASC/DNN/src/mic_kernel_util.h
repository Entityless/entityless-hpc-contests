/*------------------------------
@copyright (c) 2016 SYSU ASC16 Team

Last modified : 2016-02
Author(s) : Huang Hua(huangh223@mail2.sysu.edu.cn)

This file is part of DNNTK optimization codes for ASC16.
You can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the 
Free Software Foundation; either version 2 of the License,
or (at your option) any later version.
------------------------------*/

#ifndef MIC_KERNEL_UTIL_H
#define MIC_KERNEL_UTIL_H

#include "dnn_utility.h"
#include "my_dnn_util.h"
#include "omp.h"



__attribute__((target(mic))) void mic_setmatY    (double *Y, double *B, int row, int col);
__attribute__((target(mic))) void mic_softmaxZ   (double* in_vec, double* out_vec, int row, int col);
__attribute__((target(mic))) void mic_errorTrans (double *E, double *Y, int row, int col);
__attribute__((target(mic))) void mic_errorOutput(double *E, double *Z, int *T, int row, int col);
__attribute__((target(mic))) void mic_getBdelta  (double *E, double *Bdelta, int row, int col, double alpha);
__attribute__((target(mic))) void mic_updateArr    (double *W, double *Wdelta, int length);


#endif
