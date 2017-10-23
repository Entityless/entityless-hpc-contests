#ifndef SIGMOID_UTIL_H
#define SIGMOID_UTIL_H

#include "dnn_utility.h"
#include "my_dnn_util.h"

#include "omp.h"

void CPUsigmoidY   (double *Y, int row, int col, int ldy);
__attribute__((target(mic))) void mic_sigmoidY   (double *Y, int row, int col);

#endif // SIGMOID_UTIL_H

