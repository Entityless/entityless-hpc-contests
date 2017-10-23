#ifndef C_PUBLIC_CONST_H
#define C_PUBLIC_CONST_H

#define CCPS 1450000000

#define realt float
#define mpi_realt MPI_FLOAT
#define lbm_flag_type int

//CPF* 从核 profile 嵌套第*层
#define CPF1_WAIT              1
#define CPF1_STD_READ1         2
#define CPF1_STD_READ2         3
#define CPF1_STD_STREAM        4
#define CPF1_STD_COLLIDE       5
#define CPF1_STD_WRITE         6
#define CPF1_INSANE_READ1      7
#define CPF1_INSANE_READ2      8
#define CPF1_INSANE_STREAM     9
#define CPF1_INSANE_COLLIDE    10
#define CPF1_INSANE_WRITE      11

//MPF* 主核 profile 嵌套第*层
//使用hch_timer
#define MPF1_SEND_INIT     0
#define MPF1_COMM          1
#define MPF1_WAIT_MPI      2
#define MPF1_UPDATE        3
#define MPF1_WAIT_INNER    4
#define MPF1_HALO          5

//please make sure that 在这一层计时，没有互相嵌套的
//防止RPCC调用开销带来的系统误差
//#define PF_INITIAL_TRANSFER 1
//#define PF_EVV              2
//#define PF_CFFTS_COPY       3
//#define PF_CFFTZ_COPY       4
//#define PF_FFTZ2            5
//#define PF_TP1              6
//#define PF_TP2              7
//#define PF_CHK              8



#endif
