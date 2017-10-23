#ifndef MPE_WENO_H
#define MPE_WENO_H

#include "c_public_var.h"

#include "athread.h"

extern SLAVE_FUN(cpe_athread_daemon)();

void copy_back_boundary_u();
void copy_from_boundary_u();
void exchange_boundary();
void weno_c_core_();
void weno7_c_halo_();
void weno7_c_inner_();

#endif // MPE_WENO_H
