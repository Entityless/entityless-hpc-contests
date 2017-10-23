#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../c_header/c_public_var.h"
#include "../c_header/cpu_propagat.h"

// 缓存起来的坐标和中间变量
void c_propagat_init_(int *idxs_, float *tmp_values_)
{
    _idxs = idxs_;
    _tmp_values = tmp_values_;
}

// 双线性插值四合一系数 by 吴梓铭
void get_bilinear_interpolation_qr_(
    float *u1, float *u2, float *v1, float *v2,
    float *xt, float *yt, float *q,  float *r
)
{
    float dox, doy;

    dox = *u2 - *u1;
    doy = *v2 - *v1;

    if (dox == 0.0) // vertical line as a grid
    {
        (*q) = 0.5;
    } else {
        (*q) = (*xt - *u1) / dox;
    }

    if (doy == 0.0) // horizen line as a grid
    {
        (*r) = 0.5;
    } else {
        (*r) = (*yt - *v1) / doy;
    }
}

// CPU 端计算外圈的 propagat
void c_propagat_(int *ia_start, int *ia_end, int *ic_start, int *ic_end)
{
    if(halo_ppg_collect_flag != 0)
        return;

    int ia, ic;

    for (ic = (*ic_start); ic <= (*ic_end); ic++)
    {
        for (ia = (*ia_start); ia <= (*ia_end); ia++)
        {
            if (_nsp[(ic - _iys) * ix_size + (ia - _ixs)] != 1.0) continue;
            halo_ppg_ia_set[halo_ppg_marine_count] = ia;
            halo_ppg_ic_set[halo_ppg_marine_count] = ic;
            halo_ppg_nsp_set[halo_ppg_marine_count] = 1;
            halo_ppg_marine_count++;
        }
    }
}

//在这个函数结束收集。这个函数放在外圈ppg后面。
void collect_halo_ppg_marine_finalize_()
{
    //if(my_rank == 0)
    //    printf("    @@@@    ppg nsp set ptr = %x\n", halo_ppg_nsp_set);

    int ig, core_id, compute_unit, iaic_index, i, buff_pexx_index, pexx_index;
    float *ppg_packed_data_ptr;
    if(halo_ppg_collect_flag == 0)//整理收集成果
    {
        halo_ppg_group_count = halo_ppg_marine_count / 64;
        halo_ppg_group_remain = halo_ppg_marine_count - (halo_ppg_group_count * 64);

        halo_ppg_collect_flag = 1;
        //return;

        halo_ppg_packed_data = (float*) malloc(sizeof(float) * halo_ppg_marine_count * PPG_SV_DATA);
        assert(halo_ppg_packed_data != NULL);
        halo_ppg_out_buffer = (float*) malloc(sizeof(float) * halo_ppg_marine_count * NUM_PEXX);
        assert(halo_ppg_out_buffer != NULL);


        //拷贝不会改变的传入变量U__和D
        for (i = 0; i < halo_ppg_marine_count; i++)
        {
            iaic_index = _iaic_index(halo_ppg_ia_set[i], halo_ppg_ic_set[i]);
            //if(my_rank == 0)
            //    printf("    @@@@    collect halo ppg begin, ia = %d,  ic = %d, iaic = %d\n", halo_ppg_ia_set[i], halo_ppg_ic_set[i], iaic_index);

            ppg_packed_data_ptr = halo_ppg_packed_data + i * PPG_SV_DATA;
            ppg_packed_data_ptr[UXX_OFFSET] = _uxx[iaic_index];
            ppg_packed_data_ptr[UXY_OFFSET] = _uxy[iaic_index];
            ppg_packed_data_ptr[UYX_OFFSET] = _uyx[iaic_index];
            ppg_packed_data_ptr[UYY_OFFSET] = _uyy[iaic_index];
            ppg_packed_data_ptr[  D_OFFSET] =   _d[iaic_index];
        }
    }

    //放flag
    for(ig = 0; ig < halo_ppg_group_count; ig++)
    {
        for (core_id = 0; core_id < 64; core_id++)
        {
            compute_unit = ig * 64 + core_id;
            iaic_index = _iaic_index(halo_ppg_ia_set[compute_unit], halo_ppg_ic_set[compute_unit]);
            ppg_packed_data_ptr = halo_ppg_packed_data + compute_unit * PPG_SV_DATA;
            ppg_packed_data_ptr[  W_OFFSET] =   _w[iaic_index];
            ppg_packed_data_ptr[ WX_OFFSET] =  _wx[iaic_index];
            ppg_packed_data_ptr[ WY_OFFSET] =  _wy[iaic_index];
        }
/*
        host_flag[MPI_RANK] = my_rank;
        host_flag[KERNEL_ACTION] = HALO_PPG_FLAG;
        host_flag[GROUP_SIZE] = halo_ppg_group_count;
        host_flag[REMAIN_POINT] = halo_ppg_group_remain;
        host_flag[PACKED_PTR] = halo_ppg_packed_data;
        host_flag[OUT_DATA_PTR] = halo_ppg_out_buffer;
        host_flag[IAS_PTR] = halo_ppg_ia_set;
        host_flag[ICS_PTR] = halo_ppg_ic_set;
        host_flag[NSPS_PTR] = halo_ppg_nsp_set;
        asm volatile ("#nop":::"memory");
        host_flag[0] = host_flag[0] + 1;*/
    }

    for (core_id = 0; core_id < halo_ppg_group_remain; core_id++)
    {
        compute_unit = ig * 64 + core_id;
        iaic_index = _iaic_index(halo_ppg_ia_set[compute_unit], halo_ppg_ic_set[compute_unit]);
        ppg_packed_data_ptr = halo_ppg_packed_data + compute_unit * PPG_SV_DATA;
        ppg_packed_data_ptr[  W_OFFSET] =   _w[iaic_index];
        ppg_packed_data_ptr[ WX_OFFSET] =  _wx[iaic_index];
        ppg_packed_data_ptr[ WY_OFFSET] =  _wy[iaic_index];
    }

    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = PROPAGAT_FLAG;
    host_flag[GROUP_SIZE] = halo_ppg_group_count;
    host_flag[REMAIN_POINT] = halo_ppg_group_remain;
    host_flag[PACKED_PTR] = (long)halo_ppg_packed_data;
    host_flag[OUT_DATA_PTR] = (long)halo_ppg_out_buffer;
    host_flag[IAS_PTR] = (long)halo_ppg_ia_set;
    host_flag[ICS_PTR] = (long)halo_ppg_ic_set;
    host_flag[NSPS_PTR] = (long)halo_ppg_nsp_set;
    host_flag[TOTAL_NSP_PTR] = (long)_nsp;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + halo_ppg_group_count + 1;

    //等很多次flag
    for(i = 0; i < halo_ppg_group_count; i++)
    {
        wait_slave_flag();/*
        for (i = 64 * ig; i < 64 * (ig + 1); i++)
        {
            buff_pexx_index = i * NUM_PEXX;
            pexx_index = _iaic_index(halo_ppg_ia_set[i], halo_ppg_ic_set[i]);
            _pein[pexx_index] = halo_ppg_out_buffer[buff_pexx_index + PEIN_OFFSET];
            _pebo[pexx_index] = halo_ppg_out_buffer[buff_pexx_index + PEBO_OFFSET];
            _peds[pexx_index] = halo_ppg_out_buffer[buff_pexx_index + PEDS_OFFSET];
        }*/
    }



    wait_slave_flag();/*
    for (i = 64 * ig; i < 64 * ig + halo_ppg_group_remain; i++)
    {
        buff_pexx_index = i * NUM_PEXX;
        pexx_index = _iaic_index(halo_ppg_ia_set[i], halo_ppg_ic_set[i]);
        _pein[pexx_index] = halo_ppg_out_buffer[buff_pexx_index + PEIN_OFFSET];
        _pebo[pexx_index] = halo_ppg_out_buffer[buff_pexx_index + PEBO_OFFSET];
        _peds[pexx_index] = halo_ppg_out_buffer[buff_pexx_index + PEDS_OFFSET];
    }*/
}


void check_halo_ppg_nsp_(int* identifer)
{
    char str[100];
    int i;
    for(i = 0; i < halo_ppg_marine_count; i++)
    {
        if(halo_ppg_nsp_set[i] != 1)
        {
            sprintf(&str[0], "rank = %d, identifer = %d, ia = %d, ic = %d, halo_nsp = %d\n", my_rank, *identifer, halo_ppg_ia_set[i], halo_ppg_ic_set[i], halo_ppg_nsp_set[i]);
            printf("%s", str);
        }
    }

/*
    for(i = 0; i < ppg_inner_total_cu; i++)
    {
        if(ppg_nsp_set[i] != 1)
        {
            sprintf(&str[0], "rank = %d, identifer = %d, ia = %d, ic = %d, inner_nsp = %d\n", my_rank, *identifer, ppg_ia_set[i], ppg_ic_set[i], ppg_nsp_set[i]);
            printf("%s", str);
        }
    }*/
}
