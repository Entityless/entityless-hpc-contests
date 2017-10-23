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

#ifndef KERNEL_UTIL_H
#define KERNEL_UTIL_H

#include "dnn_utility.h"
#include "my_dnn_util.h"

void CPUsetmatY    (double *Y, double *B, int row, int col);
void CPUsoftmaxZ   (double *in_vec, double *out_vec, int row, int col);
void CPUerrorTrans (double *E, double *Y, int row, int col, int ldy);
void CPUerrorOutput(double *E, double *Z, int *T, int row, int col);
void CPUgetBdelta  (double *E, double *Bdelta, int row, int col, double alpha, int scol, int ecol);
void updateArr  (double *X, double *Xdelta, int length);

#endif
