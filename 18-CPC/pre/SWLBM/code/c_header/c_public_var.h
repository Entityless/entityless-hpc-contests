#ifndef C_PUBLIC_VAR_H
#define C_PUBLIC_VAR_H

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include "mpi.h"

#include "c_public_const.h"
#include "slave_kernel.h"
#include <assert.h>
#include <memory.h>

void wait_slave_flag();
void terminate_athread_daemon();
void athread_handshake();
extern volatile long host_flag[FLAG_SIZE];
extern volatile long slave_flag[FLAG_SIZE];
extern volatile long local_cc[FLAG_SIZE];
extern volatile long flag_to_wait;
extern int my_rank, comm_sz;

////FFF
//ptr

////mpe benchmark
extern unsigned long mpe_cc_cur_[100];
extern unsigned long mpe_cc_total_[100];

//static inline unsigned long rpcc()
//{
//    unsigned long time;
//    asm("rtc %0": "=r" (time) : );
//    return time;
//}

static inline double ccts(unsigned long cc)
{
    return cc * 1.0 / CCPS;
}

#endif
