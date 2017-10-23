#ifndef MPE_FFT_H
#define MPE_FFT_H

#include "c_public_var.h"

#include "athread.h"

extern void SLAVE_FUN(cpe_athread_daemon)();

void c_fft_iter_();

void c_fft_1d_rev_(my_fc* x1, my_fc* x2);

void c_cffts3(int is, int d1, int d2, int d3, my_fc* x, my_fc* xout, my_fc* y);

void c_cffts2(int is, int d1, int d2, int d3, my_fc* x, my_fc* xout, my_fc* y);

void c_cffts1(int is, int d1, int d2, int d3, my_fc* x, my_fc* xout, my_fc* y);

void c_transpose_x_yz(int l1, int l2, my_fc* xin, my_fc* xout);

void c_transpose2_local(int n1, int n2, my_fc* xin, my_fc* xout);

void c_transpose2_global(my_fc* xin, my_fc* xout);

void c_transpose2_finish(int n1, int n2, my_fc* xin, my_fc* xout);

void c_cfftz(int is, int m, int n, my_fc* x, my_fc* y);

void c_evolve(my_fc* u0, my_fc* u1, my_ft* twiddle, int d1, int d2, int d3);

void c_fftz2(int is, int l, int m, int n, int ny, int ny1,
              dfc* u, my_fc* x, my_fc* y);

void c_checksum(int i, my_fc* u1, int d1, int d2, int d3);

void c_checksum_rotated(int i, my_fc* u1, int d1, int d2, int d3);

void c_checksum_merge(int i);

void checksum_offset_get();

#endif // MPE_FFT_H
