#ifndef MPE_FFT_H
#define MPE_FFT_H

#include "c_public_var.h"

#include "athread.h"

extern void SLAVE_FUN(cpe_athread_daemon)();

void c_fft_iter_();

void c_fft_1d_rev_(double* x1_real, double* x1_imag, double* x2_real, double* x2_imag);

//void c_cffts3(int is, int d1, int d2, int d3, my_fc* x, my_fc* xout, my_fc* y);

//void c_cffts2(int is, int d1, int d2, int d3, my_fc* x, my_fc* xout, my_fc* y);

void c_cffts1(int is, int d1, int d2, int d3, double* x_real, double* x_imag, double* xout_real, double* xout_imag, double* y_real, double* y_imag);

void c_cffts2(int is, int d1, int d2, int d3, double* x_real, double* x_imag, double* xout_real, double* xout_imag, double* y_real, double* y_imag);

void c_transpose_x_yz(int l1, int l2, double* xin_real, double* xin_imag, double* xout_real, double* xout_imag);

void c_transpose2_local(int n1, int n2, double* xin_real, double* xin_imag, double* xout_real, double* xout_imag);

void c_transpose2_global(double* xin_real, double* xin_imag, double* xout_real, double* xout_imag);

void c_transpose2_finish(int n1, int n2, double* xin_real, double* xin_imag, double* xout_real, double* xout_imag);

void c_cfftz(int is, int m, int n, double* x_real, double* x_imag, double* y_real, double* y_imag);

void c_evolve(double* u0_real, double* u0_imag, double* u1_real, double* u1_imag, double* twiddle, int d1, int d2, int d3);

void c_fftz2(int is, int l, int m, int n, int ny, int ny1,
              double* u_real, double* u_imag, double* x_real, double* x_imag, double* y_real, double* y_imag);

void c_checksum(int i, double* u1_real, double* u1_imag, int d1, int d2, int d3);

void c_checksum_d1d2_change(int i, double* u1_real, double* u1_imag, int d1, int d2, int d3);

//void c_checksum_rotated(int i, my_fc* u1, int d1, int d2, int d3);

void c_checksum_merge(int i);

void checksum_offset_get();

#endif // MPE_FFT_H
