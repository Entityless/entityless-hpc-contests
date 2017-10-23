#ifndef ATHREAD_EMULATOR_H
#define ATHREAD_EMULATOR_H

#include "../c_header/slave_kernel.h"

float emu_bilinear_interpolation_qr(float aa, float bb, float cc, float dd, float q, float r);

void  emu_cpe_propagat_kernel(float *cpe_ppg_packed_data, float *cpe_e);

void emp_cpe_implsch_init_param(struct cpe_init_param *cpe_param);

void emu_cpe_implsch_kernel(int ig, int core_id);
#endif
