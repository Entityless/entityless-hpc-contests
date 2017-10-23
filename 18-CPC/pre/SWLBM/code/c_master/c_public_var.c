#include "../c_header/c_public_var.h"
#include "../c_header/c_public_const.h"

volatile long host_flag[FLAG_SIZE];
volatile long slave_flag[FLAG_SIZE];
volatile long local_cc[FLAG_SIZE];
volatile long flag_to_wait;
unsigned long mpe_cc_cur_[100];
unsigned long mpe_cc_total_[100];
int my_rank, comm_sz;


void wait_slave_flag()
{
    while(slave_flag[0] < flag_to_wait);
    flag_to_wait++;
}

void terminate_athread_daemon()
{
    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = EXIT_FLAG;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
    wait_slave_flag();

    athread_join();
}

void athread_handshake()
{
    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = HANDSHAKE_FLAG;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
}
