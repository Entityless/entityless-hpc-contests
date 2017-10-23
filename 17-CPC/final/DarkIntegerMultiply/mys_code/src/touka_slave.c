#include <stdio.h>
#include "slave.h"
#include "simd.h"
#include <memory.h>

#include "public_const.h"

#define CPE_TOTAL_SYNC 0x0000FFFF

//externs
extern hchjity l_num_a[20000] __attribute__((__aligned__(128)));
extern hchjity l_num_b[20000] __attribute__((__aligned__(128)));
extern int num_a_cnt, num_b_cnt;
extern volatile hchjity* local_compute_rets[CB_GTT];

extern hchjity final_val[40000] __attribute__((__aligned__(128)));

__thread_local int my_core, core_x, core_y;
__thread_local int faker_core_x, faker_core_y, faker_my_core;
__thread_local volatile unsigned long get_reply, put_reply;

//calculate
__thread_local int val;
__thread_local hchjity cpe_a[500] __attribute__((__aligned__(128)));
__thread_local hchjity cpe_b[500] __attribute__((__aligned__(128)));
__thread_local hchjity cpe_mul_buf[4000] __attribute__((__aligned__(128)));
__thread_local hchjity cpe_out_cov[1000] __attribute__((__aligned__(128)));
__thread_local volatile hchjity cpe_register_buf[1000] __attribute__((__aligned__(128)));
__thread_local int cpe_a_cnt, cpe_b_cnt;


__thread_local volatile long out_flag[FLAG_SIZE];

extern volatile long slave_flag[FLAG_SIZE];

// static __inline void __attribute__((__always_inline__)) 

void cpe_push_column_stream(int sz256, volatile intv8* bf256, int dst)
{
    int i;
    for(i = 0; i < sz256; i++)
    {
        __builtin_sw64_putC(bf256[i], dst);
    }
}

// static __inline void __attribute__((__always_inline__)) 

void cpe_read_column_stream(int sz256, volatile intv8* bf256)
{
    int i;
    for(i = 0; i < sz256; i++)
    {
        bf256[i] = __builtin_sw64_getxC(bf256[i]);
    }
}

// static __inline void __attribute__((__always_inline__)) 

void cpe_push_row_stream(int sz256, volatile intv8* bf256, int dst)
{
    int i;
    for(i = 0; i < sz256; i++)
    {
        asm volatile ("nop":::"memory");
        __builtin_sw64_putR(bf256[i], dst);
    }
}

// static __inline void __attribute__((__always_inline__)) 

void cpe_read_row_stream(int sz256, volatile intv8* bf256)
{
    int i;
    for(i = 0; i < sz256; i++)
    {
        asm volatile ("nop":::"memory");
        bf256[i] = __builtin_sw64_getxR(bf256[i]);
    }
}

// static __inline void __attribute__((__always_inline__)) 

void cpe_read_row_stream_and_add(int sz256, volatile int256* bf256, volatile int256* dst256)
{
    int i;
    for(i = 0; i < sz256; i++)
    {
        asm volatile ("nop":::"memory");
        bf256[i] = __builtin_sw64_getxR(bf256[i]);
    }

    for(i = 0; i < sz256; i++)
        dst256[i] = (int256)  __builtin_sw64_vaddl2(bf256[i], dst256[i]);
}

// static __inline void __attribute__((__always_inline__)) 

void cpe_read_col_stream_and_add(int sz256, volatile int256* bf256, volatile int256* dst256)
{
    int i;
    for(i = 0; i < sz256; i++)
    {
        asm volatile ("nop":::"memory");
        bf256[i] = __builtin_sw64_getxC(bf256[i]);
    }

    for(i = 0; i < sz256; i++)
        dst256[i] = (int256)  __builtin_sw64_vaddl2(bf256[i], dst256[i]);

}

void standard_write_flag(int my_core)
{
    out_flag[0] = out_flag[0] + 1;

    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);//this is necessary, thanks to Luo Yangze
    //after watching 10+ times of touhou pump ♂ it

    if(my_core == 0)
    {
        put_reply = 0;
        athread_put(
            PE_MODE, (void*)&out_flag[0], (void*)&slave_flag[0],
            sizeof(long) * FLAG_SIZE, (void*)(&put_reply), 0, 0
        );
        while(put_reply != 1);

    }
    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
}

void cpe_athread_daemon(void* __not_to_use_param)
{
    my_core = athread_get_core(-1);
    core_x = my_core % 8;
    core_y = my_core / 8;

//    faker_core_x = core_x;
//    faker_core_y = core_y;
//    faker_my_core = faker_core_y * 8 + faker_core_x;//wrong

    faker_core_x = (8 + core_x - core_y) % 8;//通过faker_core_x，faker_core_y来确定计算的offset
    faker_core_y = core_y;
    faker_my_core = faker_core_y * 8 + faker_core_x;

    out_flag[0] = 0;

    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);

    //gld
    //hchjity* mpe_ret = (hchjity*)local_compute_rets[my_core];
    hchjity* mpe_ret = (hchjity*)local_compute_rets[my_core];
    //my_core == 9, x = 1, y = 1, faker_x = 2, faker_y = 1, faker_my_core = 10
    cpe_a_cnt = num_a_cnt / 8;
    cpe_b_cnt = num_b_cnt / 8;


    int a_len_cb = cpe_a_cnt, b_len_cb = cpe_b_cnt;

    int a_local_offset = a_len_cb * faker_core_x;
    int b_local_offset = b_len_cb * faker_core_y;
    int local_val_len = a_len_cb + b_len_cb - 1;//即使后面使用了向量化，也只需要取这么多个
    int vec_local_val_len = local_val_len / (32 / sizeof(hchjity));//long!!!!
    if(vec_local_val_len * (32 / sizeof(hchjity)) < local_val_len)
        vec_local_val_len++;
    int vec_a_cnt = cpe_a_cnt / (32 / sizeof(hchjity));//long!!!!
    if(vec_a_cnt * (32 / sizeof(hchjity)) < cpe_a_cnt)
        vec_a_cnt++;
    int vec_b_cnt = cpe_b_cnt / (32 / sizeof(hchjity));//long!!!!
    if(vec_b_cnt * (32 / sizeof(hchjity)) < cpe_b_cnt)
        vec_b_cnt++;

    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);

    if(core_x == core_y)//dui jiao xian
    {
        //read data a b
        get_reply = 0;
        athread_get(
            PE_MODE, &l_num_a[a_local_offset], (void*)&cpe_a[0],
            sizeof(hchjity) * a_len_cb, (void*)(&get_reply), 0, 0, 0
        );
        athread_get(
            PE_MODE, &l_num_b[b_local_offset], (void*)&cpe_b[0],
            sizeof(hchjity) * b_len_cb, (void*)(&get_reply), 0, 0, 0
        );
        asm volatile ("#nop":::"memory");
        while(get_reply != 2);
        asm volatile ("#nop":::"memory");

        //bcast a on col
        cpe_push_column_stream(vec_a_cnt, (void*)&cpe_a[0], 0x0000000F);
        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);

        //bcast b on row
        cpe_push_row_stream(vec_b_cnt, (void*)&cpe_b[0], 0x0000000F);
        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    }
    else//bu shi dui jiao xian de
    {
        cpe_read_column_stream(vec_a_cnt, (void*)&cpe_a[0]);
        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);

        cpe_read_row_stream(vec_b_cnt, (void*)&cpe_b[0]);
        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    }

    int i, j;

    for(i = 0; i < local_val_len; i++)
        cpe_out_cov[i] = 0;

    //
//    for(i = 0; i < b_len_cb; i++)
//    {
//        for(j = 0; j < a_len_cb; j++)
//        {
//            cpe_out_cov[i + j] += cpe_a[j] * cpe_b[i];
//        }
//    }

    int256* out_cov_lv4 = (int256*)&cpe_out_cov[0];
    //at last, + 

    hchjity add_buf1[8] __attribute__((__aligned__(32)));
    hchjity add_buf2[8] __attribute__((__aligned__(32)));
    hchjity add_buf3[8] __attribute__((__aligned__(32)));
    hchjity add_buf4[8] __attribute__((__aligned__(32)));

    volatile int256* add_buf1_lv4p = (int256*)&add_buf1[0];
    volatile int256* add_buf2_lv4p = (int256*)&add_buf2[0];
    volatile int256* add_buf3_lv4p = (int256*)&add_buf3[0];
    volatile int256* add_buf4_lv4p = (int256*)&add_buf4[0];

//4 * 4 kernel
    int vec_id_a, vec_id_b;
    vec_id_a = 0;
    for(j = 0; j < a_len_cb; j += 4)
    {
        hchjity a1, a2, a3, a4;
        a1 = cpe_a[j + 0];
        a2 = cpe_a[j + 1];
        a3 = cpe_a[j + 2];
        a4 = cpe_a[j + 3];

        vec_id_b = 0;
        for(i = 0; i < b_len_cb; i += 4)
        {
            int ipj = i + j;

            hchjity b1, b2, b3, b4;
            b1 = cpe_b[i + 0];
            b2 = cpe_b[i + 1];
            b3 = cpe_b[i + 2];
            b4 = cpe_b[i + 3];

            cpe_out_cov[ipj + 0] += a1 * b1;
            cpe_out_cov[ipj + 1] += a2 * b1;
            cpe_out_cov[ipj + 2] += a3 * b1;
            cpe_out_cov[ipj + 3] += a4 * b1;

            cpe_out_cov[ipj + 1] += a1 * b2;
            cpe_out_cov[ipj + 2] += a2 * b2;
            cpe_out_cov[ipj + 3] += a3 * b2;
            cpe_out_cov[ipj + 4] += a4 * b2;

            cpe_out_cov[ipj + 2] += a1 * b3;
            cpe_out_cov[ipj + 3] += a2 * b3;
            cpe_out_cov[ipj + 4] += a3 * b3;
            cpe_out_cov[ipj + 5] += a4 * b3;

            cpe_out_cov[ipj + 3] += a1 * b4;
            cpe_out_cov[ipj + 4] += a2 * b4;
            cpe_out_cov[ipj + 5] += a3 * b4;
            cpe_out_cov[ipj + 6] += a4 * b4;

            vec_id_b++;
        }
        vec_id_a++;
    }


    //本地取余数
    hchjity md = 1;
    for(i = 0; i < KIU; i++)
        md *= 10;

    for(i = 0; i < local_val_len - 1; i++)
    {
        hchjity jhin = cpe_out_cov[i] / md;
        hchjity remain = cpe_out_cov[i] % md;
        cpe_out_cov[i] = remain;
        cpe_out_cov[i + 1] += jhin;
    }

    int write_out_offset = -1;//-1 for not to write out

    int gezi_field, dst_y;
    int id_in_gezi, gezi_volume;

    //下面的格子和上面的格子有不同的id映射

//    //exchange step
    {
        volatile intv8 bcast_tmp256 = {0, 1, 2, 3, 4, 5, 6, 7};//没有意义desu

        if(core_x >= core_y)
        {
            gezi_field = core_x;
            dst_y = 0;

            id_in_gezi = core_y;
            gezi_volume = gezi_field + 1;
        }
        else
        {
            gezi_field = core_x + 8;
            dst_y = 7;

            id_in_gezi = 7 - core_y;//7 always be 0, 6 -> 1
            gezi_volume = 15 - gezi_field;//8 -> 7, 14 -> 1
        }

        //buffer
        volatile intv8* out_cov_256 = (intv8*)&cpe_out_cov[0];
        volatile intv8* buf_256 = (intv8*)&cpe_register_buf[0];

        int reduce_partner, partner_core_y;//要考虑向下的情况
        int reduce_mask;
        reduce_mask = 1;//neighbor

        //进行蝶状收缩
        int reduce_phase;
        int already_sended = -1;
        for(reduce_phase = 0; reduce_phase < 3; reduce_phase++)
        {
            reduce_partner = id_in_gezi ^ reduce_mask;
            reduce_mask *= 2;

            //计算出真实的partner
            if(core_x >= core_y)//右上方
            {
                partner_core_y = reduce_partner;
            }
            else
            {
                //id = 0的总是在最下面
                partner_core_y = 7 - reduce_partner;
            }
            
            if(already_sended > 0 || reduce_partner >= gezi_volume)
            {

            }
            else if(reduce_partner < id_in_gezi)//sender
            {
                asm volatile ("nop":::"memory");
                bcast_tmp256 = __builtin_sw64_getxC(bcast_tmp256);//防止把dst变成灌汤小笼包

                for(j = 0; j < vec_local_val_len; j++)
                {
                    asm volatile ("nop":::"memory");
                    __builtin_sw64_putC(out_cov_256[j], partner_core_y);
                }

                already_sended = 1;
            }
            else//recv
            {
                asm volatile ("nop":::"memory");
                __builtin_sw64_putC(bcast_tmp256, partner_core_y);

                cpe_read_col_stream_and_add(vec_local_val_len, (void*)buf_256, (void*)&cpe_out_cov[0]);
            }

            athread_syn(COL_SCOPE, CPE_TOTAL_SYNC);
        }
        
        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);

        //phase 1
        if(core_y == 0 && core_x >= 1 && core_x <= 6)
        {
            if(core_x % 2 == 1)//攻，向右，发送高位。
            {
                int partner = core_x + 1;
                cpe_push_row_stream(vec_a_cnt, (void*)&cpe_out_cov[a_len_cb], partner);
            }
            else//拔剑四顾心茫然
            {
                cpe_read_row_stream_and_add(vec_a_cnt, (void*)&cpe_register_buf[0], (void*)&cpe_out_cov[0]);
            }
        }
        else if(my_core == 7)//7酱这个时候向63发送
        {
            cpe_push_column_stream(vec_a_cnt, (void*)&cpe_out_cov[a_len_cb], 7);
        }
        else if(my_core == 63)//63君列接收
        {
            //cpe_read_column_stream(vec_a_cnt, (void*)&cpe_register_buf[0]);
            cpe_read_column_stream(vec_a_cnt, (void*)&cpe_out_cov[a_len_cb]);
        }
        else if(core_y == 7 && core_x >= 0 && core_x <= 5)
        {
            if(core_x % 2 == 0)//攻，向右，发送高位。
            {
                int partner = core_x + 1;
                cpe_push_row_stream(vec_a_cnt, (void*)&cpe_out_cov[a_len_cb], partner);
            }
            else//拔剑四顾心茫然
            {
                cpe_read_row_stream_and_add(vec_a_cnt, (void*)&cpe_register_buf[0], (void*)&cpe_out_cov[0]);
            }
        }

        athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);//先对再快
        //phase 2
        if(core_y == 0)
        {
            if(core_x % 2 == 0)//攻，向右，发送高位。
            {
                int partner = core_x + 1;
                cpe_push_row_stream(vec_a_cnt, (void*)&cpe_out_cov[a_len_cb], partner);
            }
            else//拔剑四顾心茫然
            {
                cpe_read_row_stream(vec_a_cnt, (void*)&cpe_register_buf[0]);
                for(i = 0; i < a_len_cb - 1; i++)//向量化able
                {
                    cpe_out_cov[i] += cpe_register_buf[i];
                }
            }
        }
        else if(core_y == 7)
        {
            if(core_x % 2 == 1)//攻，向右，发送高位。记得取余数
            {
                int partner = (core_x + 1) % 8;
                cpe_push_row_stream(vec_a_cnt, (void*)&cpe_out_cov[a_len_cb], partner);
            }
            else//拔剑四顾心茫然
            {
                cpe_read_row_stream(vec_a_cnt, (void*)&cpe_register_buf[0]);
                for(i = 0; i < a_len_cb - 1; i++)//向量化able
                {
                    cpe_out_cov[i] += cpe_register_buf[i];
                }
            }
        }
    }

    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);

    if((core_y == 0 || core_y == 7) && my_core != 63)
    {
        int len = local_val_len;
        if(my_core != 62)
            len = a_len_cb;

        for(i = 0; i < len - 1; i++)//the last value is not moded
        {
            hchjity jhin = cpe_out_cov[i] / md;
            hchjity remain = cpe_out_cov[i] % md;
            cpe_out_cov[i] = remain;
            cpe_out_cov[i + 1] += jhin;
        }

        put_reply = 0;
        athread_put(
            PE_MODE, (void*)&cpe_out_cov[0], (void*)mpe_ret,
            sizeof(hchjity) * len, (void*)(&put_reply), 0, 0
        );
        asm volatile ("#nop":::"memory");
        while(put_reply != 1);
        asm volatile ("#nop":::"memory");
    }

    standard_write_flag(my_core);
}
