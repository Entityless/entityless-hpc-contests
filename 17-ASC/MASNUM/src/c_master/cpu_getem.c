#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../c_header/c_public_var.h"
#include "../c_header/cpu_implsch.h"

//extern struct cpe_init_param cpe_param;

void mean1_inner_()
{
    //host_flag[INDEXS_PTR] = mean1_packed_int;
    //if(my_rank == 0)
    //    printf("    @@@@    host mean1 ptr = %x, %x, %x, %ld\n", mean1_packed_float, mean1_out_buffer, mean1_packed_int, host_flag[INDEXS_PTR]);
    int i, ig, buff_index, val_index;
    for(i = 0; i < mean1_group_count; i++)
    {
        host_flag[MPI_RANK] = my_rank;
        host_flag[KERNEL_ACTION] = MEAN1_FLAG;
        host_flag[GROUP_SIZE] = mean1_group_count;
        host_flag[REMAIN_POINT] = mean1_group_remain;
        host_flag[PACKED_PTR] = (long)mean1_packed_float;
        host_flag[INDEXS_PTR] = (long)mean1_packed_int;
        host_flag[OUT_DATA_PTR] = (long)mean1_out_buffer;
        host_flag[IAS_PTR] = (long)mean1_ia_set;
        host_flag[ICS_PTR] = (long)mean1_ic_set;
        host_flag[NSPS_PTR] = (long)mean1_nsp_set;
        host_flag[TOTAL_NSP_PTR] =(long) _nsp;
        asm volatile ("#nop":::"memory");
        host_flag[0] = host_flag[0] + 1;
    }

    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = MEAN1_FLAG;
    host_flag[GROUP_SIZE] = mean1_group_count;
    host_flag[REMAIN_POINT] = mean1_group_remain;
    host_flag[PACKED_PTR] = (long)mean1_packed_float;
    host_flag[INDEXS_PTR] = (long)mean1_packed_int;
    host_flag[OUT_DATA_PTR] = (long)mean1_out_buffer;
    host_flag[IAS_PTR] = (long)mean1_ia_set;
    host_flag[ICS_PTR] = (long)mean1_ic_set;
    host_flag[NSPS_PTR] = (long)mean1_nsp_set;
    host_flag[TOTAL_NSP_PTR] = (long)_nsp;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;

    //cpu_mean1_emulator();
    //return;

}

void wait_mean1_inner_()
{
    long tmpcc = rpcc();
    int i, ig, buff_index, val_index;

    for(ig = 0; ig < mean1_group_count; ig++)
    {
        wait_slave_flag();
    }

    wait_slave_flag();

    for(ig = 0; ig < mean1_group_count; ig++)
    {
        for (i = 64 * ig; i < 64 * (ig + 1); i++)
        {
            if(mean1_nsp_set[i] != 1)
                continue;
            buff_index = i * GETEM_OUT_DATA;
            val_index = _iaic_index(mean1_ia_set[i], mean1_ic_set[i]);
            _h1_3[val_index] = mean1_out_buffer[buff_index + H1_3_OFFSET];
            _tpf[val_index] = mean1_out_buffer[buff_index + TPF_OFFSET];
            _ape[val_index] = mean1_out_buffer[buff_index + APE_OFFSET];
            _aet[val_index] = mean1_out_buffer[buff_index + AET_OFFSET];
        }
    }

    for (i = 64 * ig; i < 64 * ig + mean1_group_remain; i++)
    {
        if(mean1_nsp_set[i] != 1)
            continue;
        buff_index = i * GETEM_OUT_DATA;
        val_index = _iaic_index(mean1_ia_set[i], mean1_ic_set[i]);
        _h1_3[val_index] = mean1_out_buffer[buff_index + H1_3_OFFSET];
        _tpf[val_index] = mean1_out_buffer[buff_index + TPF_OFFSET];
        _ape[val_index] = mean1_out_buffer[buff_index + APE_OFFSET];
        _aet[val_index] = mean1_out_buffer[buff_index + AET_OFFSET];
    }
    mean1_wait_total += rpcc() - tmpcc;
}
