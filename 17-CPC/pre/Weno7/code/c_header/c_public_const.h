#ifndef C_PUBLIC_CONST_H
#define C_PUBLIC_CONST_H

//this is determined by the code itself
//#define WENO_LAP 4
#define CCPS 1450000000


#define PF_OUT_WAIT_FLAG    1
#define PF_WRITE_DATA       2
#define PF_GET_DATA         3
#define PF_MEMCPY_U         4
#define PF_HJ_COMPUTE       5
#define PF_U_OUT_COMPUTE    6

#define MPF_INNER           1
#define MPF_U_U1_COPY       2
#define MPF_EXCHANGE        3
#define MPF_U_COPY_BACK     4
#define MPF_U_COPY_FROM     5
#define MPF_SAVE_DATA       6

#endif
