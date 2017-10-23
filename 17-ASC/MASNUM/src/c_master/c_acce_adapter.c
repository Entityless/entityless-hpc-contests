#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <athread.h>
#include <assert.h>
#include <math.h>
#include "mpi.h"
#include "../c_header/c_public_var.h"
#include "../c_header/cpu_propagat.h"
#include "../c_header/cpu_implsch.h"
#include "../c_header/c_acce_adapter.h"
#include "../c_header/slave_kernel.h"

extern SLAVE_FUN(cpe_athread_daemon)();

struct cpe_init_param cpe_param;


int cpe_athread_launched = 0;
int mpi_rank, ips_round = 0;

// propagat inner 和 implsch inner 总共需要计算的海洋格点数和 CPE 要计算的格点数

//int ips_inner_total_cu = 0, ips_inner_cpe_cu = 0, ips_inner_cpe_groups = 0, ips_inner_mpe_cu = 0;
//int ppg_inner_total_cu = 0, ppg_inner_cpe_cu = 0, ppg_inner_cpe_groups = 0, ppg_inner_mpe_cu = 0;

float *pexx_output_buf = NULL;
float *ppg_packed_data = NULL;

float cpe_packed_consts[NUM_HPD];
int   cpe_packed_ranges[NUM_HPR];

// nsp(ia, ic) == 1.0 是海洋点，需要计算
int collect_marine_point_sets(int ia_start, int ia_end, int ic_start, int ic_end, int *ia_set, int *ic_set, int* nsp_set)
{
	int ia, ic, total_compute_unit_num = 0;
	
	for (ic = ic_start; ic <= ic_end; ic++)
	{
		for (ia = ia_start; ia <= ia_end; ia++)
		{
			if (_nsp[(ic - _iys) * ix_size + (ia - _ixs)] != 1.0) continue;
			ic_set[total_compute_unit_num] = ic;
			ia_set[total_compute_unit_num] = ia;
            nsp_set[total_compute_unit_num] = 1;
			total_compute_unit_num++;
		}
	}
	
	return total_compute_unit_num;
}

int collect_nozero_point_sets(int ia_start, int ia_end, int ic_start, int ic_end, int *ia_set, int *ic_set, int* nsp_set)
{
    int ia, ic, total_compute_unit_num = 0;
    int count1, count2;
    count1 = count2 = 0;
    
    for (ic = ic_start; ic <= ic_end; ic++)
    {
        for (ia = ia_start; ia <= ia_end; ia++)
        {
            if (_nsp[(ic - _iys) * ix_size + (ia - _ixs)] == 0.0) continue;

            if (_nsp[(ic - _iys) * ix_size + (ia - _ixs)] != 1.0)
            {
                nsp_set[total_compute_unit_num] = 2;
                count1++;
            }
            else
            {
                nsp_set[total_compute_unit_num] = 1;
                count2++;
            }

            ic_set[total_compute_unit_num] = ic;
            ia_set[total_compute_unit_num] = ia;
            total_compute_unit_num++;
        }
    }

    //printf("myid = %d, mean1 count1 = %d, count2 = %d\n", my_rank, count1, count2);
    
    return total_compute_unit_num;
}

// 重新设定需要计算的内圈 propagat 数据
void reset_ppg_param(int ppg_ia_start, int ppg_ia_end, int ppg_ic_start, int ppg_ic_end)
{
	// 海洋格点数不一定是 64 的倍数，余数交由 CPU 来计算
    ppg_inner_total_cu   = collect_marine_point_sets(ppg_ia_start, ppg_ia_end, ppg_ic_start, ppg_ic_end, ppg_ia_set, ppg_ic_set, ppg_nsp_set);
	
	// 计算的 e 的结果输出缓冲区，需要被重排然后写回到 _e 中
    pexx_output_buf = (float*) malloc(sizeof(float) * ppg_inner_total_cu * NUM_PEXX);
    assert(pexx_output_buf != NULL);

    ppg_packed_data = (float*) malloc(sizeof(float) * ppg_inner_total_cu * PPG_SV_DATA);
    assert(ppg_packed_data != NULL);

    cpe_param.host_ppg_ia_set = ppg_ia_set;
    cpe_param.host_ppg_ic_set = ppg_ic_set;
    cpe_param.host_ppg_inner_groups = ppg_inner_total_cu / 64;
    cpe_param.host_ppg_idxs = _idxs;
    cpe_param.host_ppg_values = _tmp_values;
    cpe_param.host_ppg_inner_left = ppg_inner_total_cu - cpe_param.host_ppg_inner_groups * 64;
    cpe_param.host_halo_ppg_ia_set = halo_ppg_ia_set;
    cpe_param.host_halo_ppg_ic_set = halo_ppg_ic_set;

    // 打包 uxx, uxy, uyx, uyy, d
    float* ppg_packed_data_ptr;
    int iaic_index, i;
    for (i = 0; i < ppg_inner_total_cu; i++)
    {
        iaic_index = _iaic_index(ppg_ia_set[i], ppg_ic_set[i]);
        ppg_packed_data_ptr = ppg_packed_data + i * PPG_SV_DATA;
        ppg_packed_data_ptr[UXX_OFFSET] = _uxx[iaic_index];
        ppg_packed_data_ptr[UXY_OFFSET] = _uxy[iaic_index];
        ppg_packed_data_ptr[UYX_OFFSET] = _uyx[iaic_index];
        ppg_packed_data_ptr[UYY_OFFSET] = _uyy[iaic_index];
        ppg_packed_data_ptr[  D_OFFSET] =   _d[iaic_index];
    }
}

// 重新设定需要计算的内圈 implsch 数据
void reset_ips_param(int ips_ia_start, int ips_ia_end, int ips_ic_start, int ips_ic_end)
{
    int i, iaic_index, ia, ic, k;

    // 海洋格点数不一定是 64 的倍数，余数交由 CPU 来计算
    //ips_inner_total_cu   = collect_marine_point_sets(ips_ia_start, ips_ia_end, ips_ic_start, ips_ic_end, ips_ia_set, ips_ic_set);
    //ips_inner_cpe_groups = ips_inner_total_cu / 64;
    //ips_inner_cpe_cu     = ips_inner_cpe_groups * 64;
    //ips_inner_mpe_cu     = ips_inner_total_cu - ips_inner_cpe_cu;//mpe这种称呼在暂时保留。但意义已经是remainder

    // pein, pebo, peds 输出会连续存放
    //pexx_output_buf = (float*) malloc(sizeof(float) * ips_inner_total_cu * NUM_PEXX);//20170322, 承欢将cpe改成了total
    //assert(pexx_output_buf != NULL);
	
	// 拷贝主机端的指针，准备传给 CPE
	cpe_param.host_ikp  = _ikp;
	cpe_param.host_ikp1 = _ikp1;
	cpe_param.host_ikm  = _ikm;
	cpe_param.host_ikm1 = _ikm1;
	cpe_param.host_jp1  = _jp1;
	cpe_param.host_jp2  = _jp2;
	cpe_param.host_jm1  = _jm1;
	cpe_param.host_jm2  = _jm2;
	cpe_param.host_wp   = _wp;
	cpe_param.host_wm   = _wm;
	cpe_param.host_wk   = _wk;
	cpe_param.host_dwk  = _dwk;
	cpe_param.host_wkh  = _wkh;
	cpe_param.host_wf   = _wf;
	cpe_param.host_ccg  = _ccg;
	cpe_param.host_e    = _e;
    cpe_param.host_ee   = _ee;
	cpe_param.host_wks17      = _wks17;
	cpe_param.host_grolim     = _grolim;
	cpe_param.host_cos_theta  = cos_theta;
    cpe_param.host_sin_theta  = sin_theta;
	cpe_param.host_pexx_output_buf = pexx_output_buf;
	cpe_param.host_packed_ranges   = &cpe_packed_ranges[0];
    cpe_param.host_packed_consts   = &cpe_packed_consts[0];
    cpe_param.host_ssbos = _ssbos;

    //20170422 Huang Chenghuan
    /*
    cpe_param.host_ssbos = (float*)malloc(_kldp1 * _ixl * _iyl);

    int inner_index;
    float sbo, d0, wk0, dk, ssbo;
    sbo = 0.038 * 2.0 / _g;

    inner_index = 0;
    for(ic = _iys; ic <= _iyl; ic++)
    {
        for(ia = _ixs; ia <= _ixl; ia++)
        {
            if (_nsp[(ic - _iys) * ix_size + (ia - _ixs)] != 1.0) continue;
            d0 = _d[_iaic_index(ia, ic)];
            inner_index = (ic - _iys) * ix_size * _kldp1 + (ia - _iys) * _kldp1;
            for(k = 1; k <= _kl; k++)
            {
                wk0 = _wk[k - 1];
                dk  = d0 * wk0;
                if (dk >= 30.0)
                {
                    ssbo = 0.0;
                } else {
                    ssbo = -_abo * sbo * wk0 / sinh(2.0 * dk);//pre tag
                }

                cpe_param.host_ssbos[inner_index + k - 1] = ssbo;
            }
        }
    }*/
	
	cpe_packed_ranges[IXS_OFFSET] = _ixs;
	cpe_packed_ranges[IXL_OFFSET] = _ixl;
	cpe_packed_ranges[IYS_OFFSET] = _iys;
	cpe_packed_ranges[IYL_OFFSET] = _iyl;
	
	cpe_packed_consts[DELTT5_OFFSET] = _deltt5;
	cpe_packed_consts[DELTT_OFFSET]  = _deltt;
	cpe_packed_consts[CONG_OFFSET]   = _cong;
	cpe_packed_consts[AL11_OFFSET]   = _al11;
	cpe_packed_consts[AL21_OFFSET]   = _al21;
	cpe_packed_consts[AL31_OFFSET]   = _al31;
	cpe_packed_consts[AL13_OFFSET]   = _al13;
	cpe_packed_consts[AL23_OFFSET]   = _al23;
	
	// WARNING : 下面这个输出很奇怪，如果去掉（以及 slave_kernel.c中对应的输出），就会出错！
	if (mpi_rank == 0)
	{

    }

    ips_pack_total = ips_write_total = ips_remain_total = ips_call_total = 0;
    mean1_wait_total = mean1_remain_total = mean1_copy_total = mean1_call_total = 0;
}

void mean1_get_n()
{
    int i, ia, ic;
    int _ix2 = _ixl - 1;
    float a, n;
    a = 24;
    int iaic_index;

    for(i = 0; i < (mean1_marine_count + mean1_halo_marine_count); i++)
    {
        ia = mean1_ia_set[i];
        ic = mean1_ic_set[i];
        iaic_index = _iaic_index(ia, ic);

        if(ia != 1 && ia != _ixl && ic != 1 && ic != _iyl)
        {
            n = a;
            if(_nsp[_iaic_index(ia-1,ic)]>0)n=n+1;
            if(_nsp[_iaic_index(ia+1,ic)]>0)n=n+1;
            if(_nsp[_iaic_index(ia,ic-1)]>0)n=n+1;
            if(_nsp[_iaic_index(ia,ic+1)]>0)n=n+1;
            mean1_packed_int[i * GETEM_IN_INT + LEFT_OFFSET] = -dim3;
            mean1_packed_int[i * GETEM_IN_INT + RIGHT_OFFSET] = dim3;
            mean1_packed_int[i * GETEM_IN_INT + ABOVE_OFFSET] = -dim4;
            mean1_packed_int[i * GETEM_IN_INT + BELOW_OFFSET] = dim4;
        }
        else if(_glbflag == 0 && ia == 1 && ic != 1 && ic != _iyl)
        {
            n = a;
            if(_nsp[_iaic_index(_ix2,ic)]>0)n=n+1;
            if(_nsp[_iaic_index(ia+1,ic)]>0)n=n+1;
            if(_nsp[_iaic_index(ia,ic-1)]>0)n=n+1;
            if(_nsp[_iaic_index(ia,ic+1)]>0)n=n+1;
            mean1_packed_int[i * GETEM_IN_INT + LEFT_OFFSET] = (_ix2 - ia) * dim3;
            mean1_packed_int[i * GETEM_IN_INT + RIGHT_OFFSET] = dim3;
            mean1_packed_int[i * GETEM_IN_INT + ABOVE_OFFSET] = -dim4;
            mean1_packed_int[i * GETEM_IN_INT + BELOW_OFFSET] = dim4;
        }
        else if(_glbflag == 0 && ia == _ixl && ic != 1 && ic != _iyl)
        {
            n = a;
            if(_nsp[_iaic_index(ia-1,ic)]>0)n=n+1;
            if(_nsp[_iaic_index(_ixs,ic)]>0)n=n+1;
            if(_nsp[_iaic_index(ia,ic-1)]>0)n=n+1;
            if(_nsp[_iaic_index(ia,ic+1)]>0)n=n+1;
            mean1_packed_int[i * GETEM_IN_INT + LEFT_OFFSET] = -dim3;
            mean1_packed_int[i * GETEM_IN_INT + RIGHT_OFFSET] = (_ixs - ia) * dim3;
            mean1_packed_int[i * GETEM_IN_INT + ABOVE_OFFSET] = -dim4;
            mean1_packed_int[i * GETEM_IN_INT + BELOW_OFFSET] = dim4;
        }
        else
        {
            n = a + 4;
            mean1_packed_int[i * GETEM_IN_INT + LEFT_OFFSET] = 0;
            mean1_packed_int[i * GETEM_IN_INT + RIGHT_OFFSET] = 0;
            mean1_packed_int[i * GETEM_IN_INT + ABOVE_OFFSET] = 0;
            mean1_packed_int[i * GETEM_IN_INT + BELOW_OFFSET] = 0;
        }
        mean1_n[i] = n;

        mean1_packed_float[i * GETEM_IN_FLOAT + EM_N_OFFSET] = n;
        mean1_packed_float[i * GETEM_IN_FLOAT + EM_D_OFFSET] = _d[iaic_index];
    }
}

void reset_mean1_param(int mean1_ia_start, int mean1_ia_end, int mean1_ic_start, int mean1_ic_end)
{
    //mean1_marine_count   = collect_nozero_point_sets(mean1_ia_start + 1, mean1_ia_end - 1, mean1_ic_start + 1, mean1_ic_end - 1, mean1_ia_set, mean1_ic_set, mean1_nsp_set);
    mean1_marine_count   = collect_nozero_point_sets(mean1_ia_start, mean1_ia_end, mean1_ic_start, mean1_ic_end, mean1_ia_set, mean1_ic_set, mean1_nsp_set);
    //mean1_marine_count   = collect_marine_point_sets(mean1_ia_start, mean1_ia_end, mean1_ic_start, mean1_ic_end, mean1_ia_set, mean1_ic_set);
    mean1_group_count = mean1_marine_count / 64;
    mean1_group_remain = mean1_marine_count - (mean1_group_count * 64);

    mean1_halo_marine_count = 0;

    //mean1_halo_marine_count  = collect_nozero_point_sets(mean1_ia_start, mean1_ia_end, mean1_ic_start, mean1_ic_start, mean1_ia_set + mean1_marine_count, mean1_ic_set + mean1_marine_count, mean1_nsp_set + mean1_marine_count);
    //mean1_halo_marine_count  += collect_nozero_point_sets(mean1_ia_start, mean1_ia_end, mean1_ic_end, mean1_ic_end, mean1_ia_set + (mean1_marine_count + mean1_halo_marine_count), mean1_ic_set + (mean1_marine_count + mean1_halo_marine_count), mean1_nsp_set + (mean1_marine_count + mean1_halo_marine_count));
    //mean1_halo_marine_count  += collect_nozero_point_sets(mean1_ia_start, mean1_ia_start, mean1_ic_start + 1, mean1_ic_end - 1, mean1_ia_set + (mean1_marine_count + mean1_halo_marine_count), mean1_ic_set + (mean1_marine_count + mean1_halo_marine_count), mean1_nsp_set + (mean1_marine_count + mean1_halo_marine_count));
    //mean1_halo_marine_count  += collect_nozero_point_sets(mean1_ia_end, mean1_ia_end, mean1_ic_start + 1, mean1_ic_end - 1, mean1_ia_set + (mean1_marine_count + mean1_halo_marine_count), mean1_ic_set + (mean1_marine_count + mean1_halo_marine_count), mean1_nsp_set + (mean1_marine_count + mean1_halo_marine_count));

    mean1_packed_float = (float*) malloc(sizeof(float) * (mean1_marine_count + mean1_halo_marine_count) * GETEM_IN_FLOAT);
    assert(mean1_packed_float != NULL);
    mean1_packed_int = (int*) malloc(sizeof(int) * (mean1_marine_count + mean1_halo_marine_count) * GETEM_IN_INT);
    assert(mean1_packed_int != NULL);

    mean1_out_buffer = (float*) malloc(sizeof(float) * (mean1_marine_count + mean1_halo_marine_count) * GETEM_OUT_DATA);
    assert(mean1_out_buffer != NULL);

    mean1_n = (float*) malloc(sizeof(float) * (mean1_marine_count + mean1_halo_marine_count));
    mean1_get_n();

    cpe_param.host_mean1_in_data    = mean1_packed_float;
    cpe_param.host_mean1_out_buf    = mean1_out_buffer;
    cpe_param.host_dwf =_dwf;
    cpe_param.host_glbflag = _glbflag;
    cpe_param.host_ix2 = _ix2;
    cpe_param.host_mean1_ia_set = mean1_ia_set;
    cpe_param.host_mean1_ic_set = mean1_ic_set;
}

void init_c_acce_kernel_(
        int *_ppg_ia_start, int *_ppg_ia_end, int *_ppg_ic_start, int *_ppg_ic_end,
        int *_ips_ia_start, int *_ips_ia_end, int *_ips_ic_start, int *_ips_ic_end,
        int *_mean1_ia_start, int *_mean1_ia_end, int *_mean1_ic_start, int *_mean1_ic_end,
        int *ipointp, int *jpointp
        )
{
    //char init_info_str[] = "MPI Worker %d CPE daemon launched : %d propagat; MPE remains %d propagat,.\n";
	if (cpe_athread_launched > 0) return;
	cpe_athread_launched++;

#ifdef ASCEXP12
    ipointp[0] = 225;
    jpointp[0] = 50;
    if(my_rank == 0)
        printf("defined ASCEXP12, ipointp[0] = 225; jpointp[0] = 50;\n");
#endif
#ifdef ASCEXP3
    ipointp[0] = 200;
    jpointp[0] = 300;
    if(my_rank == 0)
        printf("defined ASCEXP3, ipointp[0] = 200; jpointp[0] = 300;\n");
#endif


	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	
	// 设定 propagat 内圈计算参数
	reset_ppg_param(*_ppg_ia_start, *_ppg_ia_end, *_ppg_ic_start, *_ppg_ic_end);  
	
	// 设定 implsch 内圈计算参数
    reset_ips_param(*_ips_ia_start, *_ips_ia_end, *_ips_ic_start, *_ips_ic_end);

    // exciting
    reset_mean1_param(*_mean1_ia_start, *_mean1_ia_end, *_mean1_ic_start, *_mean1_ic_end);
	
	// 设定传给 CPE 的参数和初始化信号数组
    cpe_param.host_flag  = (long*)(&host_flag[0]);
    cpe_param.slave_flag  = (long*)(&slave_flag[0]);
    int i;
    for(i = 0; i < FLAG_SIZE; i++)
    {
        host_flag[i] = 0;
        slave_flag[i] = 0;
    }
    cpe_param.host_mpi_rank = my_rank;
    flag_to_wait = 1;
	
	// 启动从核核组守护进程
	athread_init();

//    if (mpi_rank == 0)
//    {
//        printf("Host range : %d %d %d %d\n", _ixs, _ixl, _iys, _iyl);
//        printf("Host const : %e %e %e %e %e %e %e %e \n", _deltt5, _deltt, _cong, _al11, _al21, _al31, _al13, _al23);
//        printf("Host ptr   : %x %x %x %x %x %x\n", pexx_output_buf, _e, _ee, _wf, _ccg, cpe_param.host_ssbos);
//    }

    if(my_rank == -1)
        printf(" ", cpe_param.host_packed_ranges, cpe_param.host_ikp, cpe_param.host_ikp1, cpe_param.host_ikm, cpe_param.host_ikm1, cpe_param.host_jp1,
               cpe_param.host_jp2, cpe_param.host_jm1, cpe_param.host_jm2, cpe_param.host_packed_consts, cpe_param.host_wp, cpe_param.host_wm, cpe_param.host_wks17, cpe_param.host_wk,
               cpe_param.host_dwk, cpe_param.host_grolim, cpe_param.host_wkh, cpe_param.host_cos_theta, cpe_param.host_sin_theta, cpe_param.host_wf, cpe_param.host_ccg, cpe_param.host_e,
               cpe_param.host_ee, cpe_param.host_pexx_output_buf, cpe_param.host_flag, cpe_param.slave_flag, cpe_param.host_mpi_rank, cpe_param.host_ppg_ia_set, cpe_param.host_ppg_ic_set, cpe_param.host_ppg_inner_groups,
               cpe_param.host_ppg_inner_left, cpe_param.host_ppg_idxs, cpe_param.host_ppg_values, cpe_param.host_halo_ppg_ia_set, cpe_param.host_halo_ppg_ic_set, cpe_param.host_mean1_ia_set, cpe_param.host_mean1_ic_set, cpe_param.host_mean1_in_data,
               cpe_param.host_mean1_out_buf, cpe_param.host_dwf, cpe_param.host_glbflag, cpe_param.host_ix2, cpe_param.host_ssbos);

	athread_spawn(cpe_athread_daemon, &cpe_param);
    //printf(init_info_str, mpi_rank, ppg_inner_cpe_groups, ppg_inner_mpe_cu);
}

void pack_ips_wwxwy_for_cpe()
{

}

// 等待 CPE 的数据计算完和传出来，将数据写回去 _e[]
void write_cpe_ppg_result_back(int commandv)
{
    int ig, i, buff_pexx_index, pexx_index;



    int group, remain, offset, len, torrent;


    len = ppg_inner_total_cu;
    group = len / 64;
    remain = len % 64;
    offset = 0;

    torrent = ppg_inner_total_cu * comm_torrent;


    if(commandv == 1)
    {
        len = ppg_inner_total_cu;
        group = len / 64;
        remain = len % 64;
        offset = 0;
    }
    else if(commandv == 2)
    {
        len = torrent;
        group = len / 64;
        remain = len % 64;
        offset = 0;
    }
    else if(commandv == 3)
    {
        len = ppg_inner_total_cu - torrent;
        group = len / 64;
        remain = len % 64;
        offset = torrent;
    }

    for (ig = 0; ig < group; ig++)
    {
        wait_slave_flag();
        /*for (i = 64 * ig; i < 64 * (ig + 1); i++)
        {
            buff_e_index = i * _kl * _jl;
            e_index = _e_4d_index(1, 1, ppg_ia_set[i], ppg_ic_set[i]);
            memcpy(_e + e_index, e_output_buf + buff_e_index, sizeof(float) * _kl * _jl);
        }*/

        /*for (i = 64 * ig; i < 64 * (ig + 1); i++)
        {
            buff_pexx_index = i * NUM_PEXX;
            pexx_index = _iaic_index(ppg_ia_set[i], ppg_ic_set[i]);
            _pein[pexx_index] = pexx_output_buf[buff_pexx_index + PEIN_OFFSET];
            _pebo[pexx_index] = pexx_output_buf[buff_pexx_index + PEBO_OFFSET];
            _peds[pexx_index] = pexx_output_buf[buff_pexx_index + PEDS_OFFSET];
        }*/
    }
    //if(my_rank == 0)
    //    printf("    @@@@    waiting remainder\n");
    wait_slave_flag();

    for (ig = 0; ig < group; ig++)
    {
        for (i = 64 * ig; i < 64 * (ig + 1); i++)
        {
            buff_pexx_index = i * NUM_PEXX + offset * NUM_PEXX;
            pexx_index = _iaic_index(ppg_ia_set[i + offset], ppg_ic_set[i + offset]);
            _pein[pexx_index] = pexx_output_buf[buff_pexx_index + PEIN_OFFSET];
            _pebo[pexx_index] = pexx_output_buf[buff_pexx_index + PEBO_OFFSET];
            _peds[pexx_index] = pexx_output_buf[buff_pexx_index + PEDS_OFFSET];
        }
    }

    for (i = 64 * ig; i < 64 * ig + remain; i++)
    {
        buff_pexx_index = i * NUM_PEXX + offset * NUM_PEXX;
        pexx_index = _iaic_index(ppg_ia_set[i + offset], ppg_ic_set[i + offset]);
        _pein[pexx_index] = pexx_output_buf[buff_pexx_index + PEIN_OFFSET];
        _pebo[pexx_index] = pexx_output_buf[buff_pexx_index + PEBO_OFFSET];
        _peds[pexx_index] = pexx_output_buf[buff_pexx_index + PEDS_OFFSET];
    }

    /*
    for (i = 64 * ig; i < 64 * ig + ppg_inner_mpe_cu; i++)
    {
        buff_pexx_index = i * NUM_PEXX;
        pexx_index = _iaic_index(ppg_ia_set[i], ppg_ic_set[i]);
        _pein[pexx_index] = pexx_output_buf[buff_pexx_index + PEIN_OFFSET];
        _pebo[pexx_index] = pexx_output_buf[buff_pexx_index + PEBO_OFFSET];
        _peds[pexx_index] = pexx_output_buf[buff_pexx_index + PEDS_OFFSET];
    }*/

    /*for (i = 64 * ig; i < 64 * ig + ppg_inner_mpe_cu; i++)
    {
        buff_e_index = i * _kl * _jl;
        e_index = _e_4d_index(1, 1, ppg_ia_set[i], ppg_ic_set[i]);
        memcpy(_e + e_index, e_output_buf + buff_e_index, sizeof(float) * _kl * _jl);
    }*/
}

// 等待 CPE 的数据计算完和传出来，将数据写回去 _pe{ds, bo, in}[]
void write_cpe_ips_pexx_result_back()
{
    return;
    int index_tmp, ia, ic, kj, i;

    for(i = 0; i < ppg_inner_total_cu; i++)
    {
        ia = ppg_ia_set[i];
        ic = ppg_ic_set[i];

        index_tmp = _e_4d_index(1, 1, ia, ic);


        for(kj = 1; kj <= dim3; kj++)
        {
            _ee[index_tmp] = _e[index_tmp];
            index_tmp++;
        }
    }
}

void c_propagat_inner(int commandv)
{

    //if(my_rank == 0)
    //    printf("    @@@@    ppg nsp set ptr = %x\n", ppg_nsp_set);

    // 打包数据等待 CPE 读取和计算
    int ig, core_id, i, iaic_index, compute_unit;
    float *ppg_packed_data_ptr;

    int group, remain, offset, len;
    int torrent;


    len = ppg_inner_total_cu;
    group = len / 64;
    remain = len % 64;
    offset = 0;

    torrent = ppg_inner_total_cu * comm_torrent;


    if(commandv == 1)
    {
        len = ppg_inner_total_cu;
        group = len / 64;
        remain = len % 64;
        offset = 0;
    }
    else if(commandv == 2)
    {
        len = torrent;
        group = len / 64;
        remain = len % 64;
        offset = 0;
    }
    else if(commandv == 3)
    {
        len = ppg_inner_total_cu - torrent;
        group = len / 64;
        remain = len % 64;
        offset = torrent;
    }

    for (ig = 0; ig < group; ig++)
    {
        for (core_id = 0; core_id < 64; core_id++)
        {
            compute_unit = ig * 64 + core_id;
            iaic_index = _iaic_index(ppg_ia_set[compute_unit + offset], ppg_ic_set[compute_unit + offset]);
            ppg_packed_data_ptr = ppg_packed_data + (compute_unit + offset) * PPG_SV_DATA;
            ppg_packed_data_ptr[  W_OFFSET] =   _w[iaic_index];
            ppg_packed_data_ptr[ WX_OFFSET] =  _wx[iaic_index];
            ppg_packed_data_ptr[ WY_OFFSET] =  _wy[iaic_index];
        }
/*
        host_flag[MPI_RANK] = my_rank;
        host_flag[KERNEL_ACTION] = PROPAGAT_FLAG;
        host_flag[GROUP_SIZE] = ppg_inner_cpe_groups;//还是标准分组大小。+1的操作在从核上面完成。
        host_flag[REMAIN_POINT] = ppg_inner_mpe_cu;
        host_flag[PACKED_PTR] = ppg_packed_data;
        host_flag[OUT_DATA_PTR] = pexx_output_buf;
        host_flag[IAS_PTR] = ppg_ia_set;
        host_flag[ICS_PTR] = ppg_ic_set;
        host_flag[NSPS_PTR] = ppg_nsp_set;
        asm volatile ("#nop":::"memory");
        host_flag[0] = host_flag[0] + 1;*/
    }

    for (core_id = 0; core_id < remain; core_id++)
    {
        compute_unit = ig * 64 + core_id;
        iaic_index = _iaic_index(ppg_ia_set[compute_unit + offset], ppg_ic_set[compute_unit + offset]);
        ppg_packed_data_ptr = ppg_packed_data + (compute_unit + offset) * PPG_SV_DATA;
        ppg_packed_data_ptr[  W_OFFSET] =   _w[iaic_index];
        ppg_packed_data_ptr[ WX_OFFSET] =  _wx[iaic_index];
        ppg_packed_data_ptr[ WY_OFFSET] =  _wy[iaic_index];
    }

    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = PROPAGAT_FLAG;
    host_flag[GROUP_SIZE] = group;//还是标准分组大小。+1的操作在从核上面完成。
    host_flag[REMAIN_POINT] = remain;
    host_flag[PACKED_PTR] = (long)(ppg_packed_data + (offset * PPG_SV_DATA));
    host_flag[OUT_DATA_PTR] = (long)(pexx_output_buf + (offset * NUM_PEXX));
    host_flag[IAS_PTR] = (long)(ppg_ia_set + offset);
    host_flag[ICS_PTR] = (long)(ppg_ic_set + offset);
    host_flag[NSPS_PTR] = (long)(ppg_nsp_set + offset);
    host_flag[TOTAL_NSP_PTR] = (long)(_nsp);
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + group + 1;
}

void c_propagat_inner_stop(int commandv)
{
    // 等待 CPE 的数据计算完和传出来，将数据写回去 _e
    write_cpe_ppg_result_back(commandv);
}

void c_implsch_inner()
{
    // 打包 _w[], _wx[], _wy[] 数据
    long tmpcc;

    tmpcc = rpcc();
    pack_ips_wwxwy_for_cpe();
    ips_pack_total += rpcc() - tmpcc;

    // 等待 CPE 的数据计算完和传出来
    tmpcc = rpcc();
    write_cpe_ips_pexx_result_back();
    ips_write_total += rpcc() - tmpcc;

    ips_call_total = ips_call_total + 1;
}

void c_implsch_inner_()  // wrapper for Fortran interface
{
    c_implsch_inner();
}

//Huang Chenghuan, 20170318

void init_main_(int* myid)
{
    my_rank = *myid;
    if(my_rank == -1)
        printf("my_rank = %d\n", my_rank);
}

void c_time_print_(double* times)
{
    int i;
    double tt = times[17];
    double cpe_ppg_av = cpe_param.slave_flag[1] * 1.0 / CCPS, cpe_halo_ppg_av = cpe_param.slave_flag[2] * 1.0 / CCPS, cpe_mean1_av = cpe_param.slave_flag[3] * 1.0 / CCPS;

    double out_time[FLAG_SIZE];
    for(i = 1; i < 32; i++)
    {
        out_time[i] = cpe_param.slave_flag[i] * 1.0 / CCPS;
    }

    if(my_rank == 0)
    {
        //printf("flag 1 = %.5f, 2 = %.5f\n", cpe_ppg_av, cpe_halo_ppg_av);
    }

    ips_pack_total /= ips_call_total;
    ips_write_total /= ips_call_total;
    ips_remain_total /= ips_call_total;
    mean1_wait_total /= mean1_call_total;
    mean1_remain_total /= mean1_call_total;
    mean1_copy_total /= mean1_call_total;

    double mpe_ips_pack_av = ips_pack_total * 1.0 / CCPS, mpe_ips_remain_av = ips_remain_total * 1.0 / CCPS, mpe_ips_write_av = ips_write_total * 1.0 / CCPS;

    double mpe_mean1_wait_av = mean1_wait_total * 1.0 / CCPS, mpe_mean1_remain_av = mean1_remain_total * 1.0 / CCPS, mpe_mean1_copy_av = mean1_copy_total * 1.0 / CCPS;

    for(i = 0; i < 17; i++)
    {
        times[17] -= times[i];
    }

    asm volatile ("#nop":::"memory");
    printf("%d, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n",
           my_rank, tt, times[0], times[1], times[2], times[3], times[4], times[5], times[6], times[7], times[8], times[9], times[10], times[11], times[12], times[13], times[14], times[15], times[16], times[17],
           out_time[PF_PPG_DATA],
            out_time[PF_PPG_KERNEL],
            out_time[PF_IPS_DATA],
            out_time[PF_MEAN2],
            out_time[PF_GETTH],
            out_time[PF_NOLIN],
            out_time[PF_SINPUT],
            out_time[PF_DISSIP],
            out_time[PF_BOTTOM],
            out_time[PF_SCURRENT],
            out_time[PF_UPDATEEE],
            out_time[PF_IPS_WRITE],
            out_time[PF_MEAN1_DATA],
            out_time[PF_MEAN1_E],
            out_time[PF_EM],
            out_time[PF_MEAN1_1],
            out_time[PF_MEAN1_2],
            out_time[PF_MEAN1_WRITE],
            out_time[IPF_NONLIN_1],
            out_time[IPF_NONLIN_2],
            out_time[IPF_NONLIN_3],
            out_time[IPF_NONLIN_4],
            out_time[PF_NOLIN] - out_time[IPF_NONLIN_1] - out_time[IPF_NONLIN_2] - out_time[IPF_NONLIN_3] - out_time[IPF_NONLIN_4],
            out_time[PF_MEAN1_DATA2]);

    times[18] = out_time[PF_PPG_KERNEL];
    times[19] = out_time[PF_MEAN2] + out_time[PF_GETTH]
            + out_time[PF_NOLIN] + out_time[PF_SINPUT]
            + out_time[PF_DISSIP] + + out_time[PF_BOTTOM]
            + out_time[PF_SCURRENT] + + out_time[PF_UPDATEEE];
    times[20] = out_time[PF_PPG_DATA] + out_time[PF_IPS_DATA]
            + out_time[PF_IPS_WRITE] + out_time[PF_MEAN1_DATA]
            + out_time[PF_MEAN1_WRITE];

}

void c_athread_finalize_()
{
    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = EXIT_FLAG;
    host_flag[GROUP_SIZE] = -1;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
    wait_slave_flag();

    if(my_rank == 0)
    {
        //printf("flag 1 = %ld, 2 = %ld\n", cpe_param.slave_flag[1], cpe_param.slave_flag[2]);
    }

    athread_join();
    //printf("athread_join from %d\n", my_rank);
}

void check_inner_ppg_nsp_(int* identifer)
{
    char str[100];
    int i;

    for(i = 0; i < ppg_inner_total_cu; i++)
    {
        if(ppg_nsp_set[i] != 1)
        {
            sprintf(&str[0], "rank = %d, identifer = %d, ia = %d, ic = %d, inner_nsp = %d\n", my_rank, *identifer, ppg_ia_set[i], ppg_ic_set[i], ppg_nsp_set[i]);
            printf("%s", str);
        }
    }
}
