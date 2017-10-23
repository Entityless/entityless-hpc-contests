#ifndef SLAVE_KERNEL_H
#define SLAVE_KERNEL_H

struct cpe_init_param
{
    int   *host_packed_ranges;
	int   *host_ikp, *host_ikp1, *host_ikm, *host_ikm1;
	int   *host_jp1, *host_jp2, *host_jm1, *host_jm2;
    float *host_packed_consts;
    float *host_wp, *host_wm, *host_wks17, *host_wk;
	float *host_dwk, *host_grolim, *host_wkh, *host_cos_theta, *host_sin_theta;
    float *host_wf, *host_ccg, *host_e, *host_ee, *host_pexx_output_buf;
    long  *host_flag, *slave_flag;
    int   host_mpi_rank;
    int   *host_ppg_ia_set, *host_ppg_ic_set;
    int   host_ppg_inner_groups;
    int   host_ppg_inner_left;//如果这个是0，那么所有从核就一起多写一次flag
    int   *host_ppg_idxs;
    float *host_ppg_values;
    int   *host_halo_ppg_ia_set, *host_halo_ppg_ic_set;
    int   *host_mean1_ia_set, *host_mean1_ic_set;
    float *host_mean1_in_data, *host_mean1_out_buf;
    float *host_dwf;
    int   host_glbflag, host_ix2;
    float   *host_ssbos;
};

void cpe_athread_daemon(void *_param);

#define CPE_TOTAL_SYNC 0x0000FFFF


#define PROPAGAT_FLAG 2233
#define EXIT_FLAG     114514
#define MEAN1_FLAG    666666

// 信号数组每一位的作用
#define FLAG_SIZE        32
#define MPI_RANK         1  // 进程的 MPI 编号，用于打印调试信息
#define KERNEL_ACTION    2  // ppg, ips, exit
#define GROUP_SIZE       3
#define REMAIN_POINT     4
#define PACKED_PTR       5
#define OUT_DATA_PTR     6
#define INDEXS_PTR       7
#define IAS_PTR          8
#define ICS_PTR          9
#define NSPS_PTR         10
#define TOTAL_NSP_PTR    11

//输出信号量存储的数据. from 1 to flag_size - 1
#define PF_PPG_DATA      1
#define PF_PPG_KERNEL    2
#define PF_IPS_DATA      3
#define PF_MEAN2         4
#define PF_GETTH         5
#define PF_NOLIN         6
#define PF_SINPUT        7
#define PF_DISSIP        8
#define PF_BOTTOM        9
#define PF_SCURRENT      10
#define PF_UPDATEEE      11
#define PF_IPS_WRITE     12
#define PF_MEAN1_DATA    13
#define PF_MEAN1_E       14
#define PF_EM            15
#define PF_MEAN1_1       16
#define PF_MEAN1_2       17
#define PF_MEAN1_WRITE   18
#define IPF_NONLIN_1     19
#define IPF_NONLIN_2     20
#define IPF_NONLIN_3     21
#define IPF_NONLIN_4     22
#define PF_MEAN1_DATA2   23

#define PPG_SV_DATA    8  // implsch 中每一次 (ia, ic) 用到的 8 个单值
#define UXX_OFFSET     0
#define UXY_OFFSET     1
#define UYX_OFFSET     2
#define UYY_OFFSET     3
#define D_OFFSET       4
#define W_OFFSET       5
#define WX_OFFSET      6
#define WY_OFFSET      7

#define NUM_PEXX       3  // pein, pebo, peds 共三个
#define PEIN_OFFSET    0
#define PEBO_OFFSET    1
#define PEDS_OFFSET    2

#define NUM_HPD        8  // 8 个常量
#define DELTT5_OFFSET  0
#define DELTT_OFFSET   1
#define CONG_OFFSET    2
#define AL11_OFFSET    3
#define AL21_OFFSET    4
#define AL31_OFFSET    5
#define AL13_OFFSET    6
#define AL23_OFFSET    7

#define NUM_HPR        4  // 4 个范围值
#define IXS_OFFSET     0
#define IYS_OFFSET     1
#define IXL_OFFSET     2
#define IYL_OFFSET     3


#define IDXS_IXX_OFF   0
#define IDXS_JYY_OFF   1
#define IDXS_JTH_OFF   2
#define IDXS_IWK_OFF   3
#define IDXS_IWK1_OFF  4

#define VALUES_Q_4_OFF 0
#define VALUES_R_4_OFF 1
#define VALUES_Q_OFF   2
#define VALUES_R_OFF   3
#define VALUES_FIEN_OFF  4
#define VALUES_FIEN1_OFF 5

/*
#define GETEM_OUT_DATA 8
#define ASI_OFFSET     0
#define AWF_OFFSET     1
#define AETC_OFFSET    2
#define AETS_OFFSET    3
#define AWK_OFFSET     4
#define AE_OFFSET      5
#define APE_OFFSET     6*/

#define GETEM_OUT_DATA 4
#define H1_3_OFFSET    0
#define TPF_OFFSET     1
#define APE_OFFSET     2
#define AET_OFFSET     3

#define GETEM_IN_FLOAT  2
//#define EM_DWF_OFFSET  0
//#define EM_WF_OFFSET   1
#define EM_D_OFFSET    0
#define EM_N_OFFSET    1

#define GETEM_IN_INT   4
#define LEFT_OFFSET    0
#define RIGHT_OFFSET   1
#define ABOVE_OFFSET   2
#define BELOW_OFFSET   3

#endif
