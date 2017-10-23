#include "sigmoid_util.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void CPUsigmoidY(double *Y, int row, int col, int ldy)
{
    for (int r = 0; r < row; r++)
    {
        int base = r * ldy;
        for (int c = 0; c < col; c++)
            Y[base + c] = 1.0f / (1.0f + exp(-Y[base + c]));
    }
}

__attribute__((target(mic)))
void mic_sigmoidY(double *Y, int row, int col)
{
#pragma omp parallel for
    for (int i = 0; i < row * col; i++)
        Y[i] = 1.0f / (1.0f + exp(-Y[i]));
}
