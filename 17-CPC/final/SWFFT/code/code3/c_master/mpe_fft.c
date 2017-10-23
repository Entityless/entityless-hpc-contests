#include "../c_header/mpe_fft.h"

//在过渡版代码中，主核u0已经不是evv的地方了。在完成传递的任务之后，就是从核的输出数据的地方。
void total_dt_transpose_x_y(double* in, double* out)
{
#define ux_idx(a, b, c) ((a - 1) + ((b - 1) + (c - 1) * ny_) * nx_)

    int i, j, k;

    for(i = 1; i <= nx_; i++)
    {
        for(j = 1; j <= ny_; j++)
        {
            for(k = 1; k <= nz_; k++)
            {
                out[ux_idx(j, i, k)] = in[ux_idx(i, j, k)];
            }
        }
    }
#undef ux_idx
}

void total_dt_transpose_x_z(double* in, double* out)
{
#define ux_idx(a, b, c) ((a - 1) + ((b - 1) + (c - 1) * ny_) * nx_)

    int i, j, k;

    for(i = 1; i <= nx_; i++)
    {
        for(j = 1; j <= ny_; j++)
        {
            for(k = 1; k <= nz_; k++)
            {
                out[ux_idx(k, j, i)] = in[ux_idx(i, j, k)];
            }
        }
    }
#undef ux_idx
}

//void local_transpose_x_y(double* in, double* out, int d1, int d2, int d3)
//{
//    host_flag[MPI_RANK] = me_;
//    host_flag[KERNEL_ACTION] = TP_X_Y_FLAG;
//    host_flag[IN_PTR_REAL] = (long)in;
//    host_flag[OUT_PTR_REAL] = (long)out;
//    host_flag[FFT_D1] = d1;
//    host_flag[FFT_D2] = d2;
//    host_flag[FFT_D3] = d3;
//    asm volatile ("#nop":::"memory");
//    host_flag[0] = host_flag[0] + 1;
//    wait_slave_flag();


////#define ux_idx(a, b, c) ((a - 1) + ((b - 1) + (c - 1) * d2) * d1)

////    int i, j, k;

////    for(i = 1; i <= d1; i++)
////    {
////        for(j = 1; j <= d2; j++)
////        {
////            for(k = 1; k <= d3; k++)
////            {
////                out[ux_idx(j, i, k)] = in[ux_idx(i, j, k)];
////            }
////        }
////    }
////#undef ux_idx
//}

void c_fft_iter_()
{
    int iter, i, j;
    int fuyi = -1;
    unsigned long stcc, edcc;

    checksum_offset_get();

    struct fft_init_param param;
    param.my_id = me_;
    param.iter = niter_;
    param.host_flag = (long*)&host_flag[0];
    param.slave_flag = (long*)&slave_flag[0];
    param.nx = nx_;
    param.ny = ny_;
    param.nz = nz_;
    param.npc = (((nx_ * ny_) / np2_) * nz_) / 64;
    param.u_real_left = c_u_real_left_;
    param.u_real_right = c_u_real_right_;
    param.last_ureal = c_u_real_[1];
    param.u_imag_left = c_u_imag_left_;
    param.u_imag_right = c_u_imag_right_;
    param.last_uimag = c_u_imag_[1];
    param.twiddle = c_twiddle_;
    param.buf_chk = buf_chk_;
    param.chk_cnt = chk_core_cnt_;
    param.chk_stp = chk_core_stp_;
    param.chk_xs = chk_xs_;
    param.chk_ys = chk_ys_;
    param.chk_zs = chk_zs_;

    if(param.my_id == 0)
    {
        printf("*** MESSAGE FROM MPE ***\n");
        printf("iter = %d, nx = %d, ny = %d, nz = %d, npc = %d\n", param.iter, param.nx, param.ny, param.nz, param.npc);
        printf("ptrs: %x, %x, %x, %x, %x, %x, %x\n%x %x %x %x %x %x\n",
               param.host_flag, param.slave_flag, param.u_real_left, param.u_real_right, param.u_imag_left, param.u_imag_right, param.twiddle,
               param.chk_cnt, param.chk_stp, param.chk_xs, param.chk_ys, param.chk_zs, param.buf_chk);
        printf("*** MESSAGE FROM MPE ***\n");
    }

    athread_init();

    if(me_ == 0)
    {
        for(i = 0; i < 10; i++)
        {
            printf("MPE c_u_[%d] = %.8f %.8f\n", i, c_u_real_[i], c_u_imag_[i]);
        }
    }

    athread_spawn(cpe_athread_daemon, (void*)&param);

//    niter_ = 2;

    fftblockpad_ = 1;
    fftblock_ = 1;
    fftblockpad_default_ = 1;

    double* cur_using_u1_real= c_u1_real_bf1_, *cur_using_u2_real = c_u2_real_bf1_;
    double* cur_pending_u1_real = c_u1_real_bf2_, *cur_pending_u2_real = c_u2_real_bf2_;

    double* cur_using_u1_imag = c_u1_imag_bf1_, *cur_using_u2_imag = c_u2_imag_bf1_;
    double* cur_pending_u1_imag = c_u1_imag_bf2_, *cur_pending_u2_imag = c_u2_imag_bf2_;

    int l1 = 3, l2 = 2;
    int dim_mul = dims(2, l1) * dims(3, l1);

    int n1, n2;
    n1 = dims(1, l1);
    n2 = dim_mul;

    for(iter = 1; iter <= niter_; iter++)
    {
        //fft_(&fuyi, u1_, u2_);
        //c_fft_1d_rev_(cur_using_u1_real, cur_using_u1_imag, cur_using_u2_real, cur_using_u2_imag);

        int fuyi = -1, san = 3, er = 2, yi  = 1;
        stcc = rpcc();
        host_flag[MPI_RANK] = me_;
    //    host_flag[KERNEL_ACTION] = COUTLE_FFT_FLAG;
        host_flag[KERNEL_ACTION] = CP_EVV_FFT_FLAG;
        host_flag[GROUP_SIZE] = (dims(2, 3) * dims(3, 3)) / 32;
        host_flag[REMAIN_POINT] = (dims(2, 3) * dims(3, 3)) % 32;
    //    host_flag[IN_PTR_REAL] = (long)x1_real;
    //    host_flag[IN_PTR_IMAG] = (long)x1_imag;
        host_flag[IN_PTR_REAL] = (long)c_u0_real_;
        host_flag[IN_PTR_IMAG] = (long)c_u0_imag_;
        host_flag[OUT_PTR_REAL] = (long)cur_using_u1_real;
        host_flag[OUT_PTR_IMAG] = (long)cur_using_u1_imag;
        host_flag[IN_STRIDE] = dims(1, 3);
        host_flag[OUT_STRIDE] = dims(1, 3);
        host_flag[FFT_D1] = dims(1, 3);
        host_flag[FFT_D1] = dims(1, 3);
        host_flag[FFT_D2] = dims(2, 3);
        host_flag[FFT_D3] = dims(3, 3);
        host_flag[SUP_PTR_REAL] = (long)cur_using_u2_real;
        host_flag[SUP_PTR_IMAG] = (long)cur_using_u2_imag;
        asm volatile ("#nop":::"memory");
        host_flag[0] = host_flag[0] + 1;
        wait_slave_flag();
        wait_slave_flag();
        wait_slave_flag();
        edcc = rpcc();
        mpe_cc_cur_[MPE_FT1] = edcc - stcc;


        if(iter > 1)
        {
            stcc = rpcc();
            host_flag[MPI_RANK] = me_;
            host_flag[KERNEL_ACTION] = DB_FFT_FLAG;
            host_flag[GROUP_SIZE] = (dims(2, 2) * dims(3, 2)) / 32;
            host_flag[REMAIN_POINT] = (dims(2, 2) * dims(3, 2)) % 32;
            host_flag[IN_PTR_REAL] = (long)cur_pending_u2_real;
            host_flag[IN_PTR_IMAG] = (long)cur_pending_u2_imag;
            host_flag[OUT_PTR_REAL] = (long)cur_pending_u1_real;//buffer
            host_flag[OUT_PTR_IMAG] = (long)cur_pending_u1_imag;//buffer, 2 -> 2,
            host_flag[IN_STRIDE] = dims(1, 2);
            host_flag[OUT_STRIDE] = dims(1, 2);
            host_flag[FFT_D1] = dims(1, 2);
            host_flag[FFT_D2] = dims(2, 2);
            host_flag[FFT_D3] = dims(3, 2);
            host_flag[FFT_D1_FNS] = n2;
            host_flag[FFT_D2_FNS] = n1 / np2_;
            host_flag[FFT_D3_FNS] = np2_;
            asm volatile ("#nop":::"memory");
            host_flag[0] = host_flag[0] + 1;
            //c_cffts1(fuyi, dims(1, 2), dims(2, 2), dims(3, 2), x2_real, x2_imag, x2_real, x2_imag, scratch_real_, scratch_imag_);
            edcc = rpcc();
            mpe_cc_cur_[MPE_FT2] = edcc - stcc;
        }

        stcc = rpcc();
        //transpose_x_yz_(&san, &er, x1, x2);
        c_transpose_x_yz(san, er, cur_using_u1_real, cur_using_u1_imag, cur_using_u2_real, cur_using_u2_imag);
        edcc = rpcc();
        //mpe_cc_cur_[MPE_TP1] = edcc - stcc;

        if(iter > 1)
        {
            stcc = rpcc();
            wait_slave_flag();
            wait_slave_flag();
            wait_slave_flag();
            wait_slave_flag();
            wait_slave_flag();
            wait_slave_flag();
            edcc = rpcc();
            mpe_cc_cur_[MPE_FT2] = edcc - stcc;

            stcc = rpcc();
            //checksum_(&iter, u2_, &dims(1, 1), &dims(2, 1), &dims(3, 1));
            c_checksum_d1d2_change(iter - 1, cur_pending_u2_real, cur_pending_u2_imag, dims(1, 1), dims(2, 1), dims(3, 1));
            edcc = rpcc();
            mpe_cc_cur_[MPE_CHK] = edcc - stcc;
        }

//        stcc = rpcc();
//        //checksum_(&iter, u2_, &dims(1, 1), &dims(2, 1), &dims(3, 1));
//        c_checksum_d1d2_change(iter, cur_using_u2_real, cur_using_u2_imag, dims(1, 1), dims(2, 1), dims(3, 1));
//        edcc = rpcc();
//        mpe_cc_cur_[MPE_CHK] = edcc - stcc;

        if(me_ == 0)
            printf("iter = %d, evv = %.6f, fft1 = %6f, fft2 = %6f, fft3 = %6f, t2l = %6f, t2g = %6f, t2f = %6f, chk = %.6f\n",
                   iter, ccts(mpe_cc_cur_[MPE_EVV]),
                   ccts(mpe_cc_cur_[MPE_FT1]), ccts(mpe_cc_cur_[MPE_FT2]), ccts(mpe_cc_cur_[MPE_FT3]),
                   ccts(mpe_cc_cur_[MPE_T2L]), ccts(mpe_cc_cur_[MPE_T2G]), ccts(mpe_cc_cur_[MPE_T2F]),
                   ccts(mpe_cc_cur_[MPE_CHK]));

        for(i = 0; i < 50; i++)
            mpe_cc_total_[i] += mpe_cc_cur_[i];

        //swap
        {
            double* tmpdp;
            tmpdp = cur_using_u1_real;
            cur_using_u1_real = cur_pending_u1_real;
            cur_pending_u1_real = tmpdp;
            tmpdp = cur_using_u2_real;
            cur_using_u2_real = cur_pending_u2_real;
            cur_pending_u2_real = tmpdp;

            tmpdp = cur_using_u1_imag;
            cur_using_u1_imag = cur_pending_u1_imag;
            cur_pending_u1_imag = tmpdp;
            tmpdp = cur_using_u2_imag;
            cur_using_u2_imag = cur_pending_u2_imag;
            cur_pending_u2_imag = tmpdp;
        }
    }

    host_flag[MPI_RANK] = me_;
    host_flag[KERNEL_ACTION] = DB_FFT_FLAG;
    host_flag[GROUP_SIZE] = (dims(2, 2) * dims(3, 2)) / 32;
    host_flag[REMAIN_POINT] = (dims(2, 2) * dims(3, 2)) % 32;
    host_flag[IN_PTR_REAL] = (long)cur_pending_u2_real;
    host_flag[IN_PTR_IMAG] = (long)cur_pending_u2_imag;
    host_flag[OUT_PTR_REAL] = (long)cur_pending_u1_real;//buffer
    host_flag[OUT_PTR_IMAG] = (long)cur_pending_u1_imag;//buffer, 2 -> 2,
    host_flag[IN_STRIDE] = dims(1, 2);
    host_flag[OUT_STRIDE] = dims(1, 2);
    host_flag[FFT_D1] = dims(1, 2);
    host_flag[FFT_D2] = dims(2, 2);
    host_flag[FFT_D3] = dims(3, 2);
    host_flag[FFT_D1_FNS] = n2;
    host_flag[FFT_D2_FNS] = n1 / np2_;
    host_flag[FFT_D3_FNS] = np2_;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;

    stcc = rpcc();
    wait_slave_flag();
    wait_slave_flag();
    wait_slave_flag();
    wait_slave_flag();
    wait_slave_flag();
    wait_slave_flag();
    edcc = rpcc();
    mpe_cc_cur_[MPE_FT2] = edcc - stcc;

    stcc = rpcc();
    //checksum_(&iter, u2_, &dims(1, 1), &dims(2, 1), &dims(3, 1));
    c_checksum_d1d2_change(iter - 1, cur_pending_u2_real, cur_pending_u2_imag, dims(1, 1), dims(2, 1), dims(3, 1));
    edcc = rpcc();
    mpe_cc_cur_[MPE_CHK] = edcc - stcc;

    terminate_athread_daemon();
}

void c_evolve(double* u0_real, double* u0_imag, double* u1_real, double* u1_imag, double* twiddle, int d1, int d2, int d3)
{
#define ux_idx(a, b, c) ((a - 1) + ((b - 1) + (c - 1) * d2) * d1)
    //double complex u0(d1,d2,d3)
    int i, j, k;

    for(k = 1; k <= d3; k++)
    {
        for(j = 1; j <= d2; j++)
        {
            for(i = 1; i <= d1; i++)
            {
                //u0[ux_idx(i, j, k)] = multiply_dc_d(u0[ux_idx(i, j, k)], twiddle[ux_idx(i, j, k)]);
                //u1[ux_idx(i, j, k)] = u0[ux_idx(i, j, k)];

                u0_real[ux_idx(i, j, k)] = u0_real[ux_idx(i, j, k)] * twiddle[ux_idx(i, j, k)];
                u0_imag[ux_idx(i, j, k)] = u0_imag[ux_idx(i, j, k)] * twiddle[ux_idx(i, j, k)];

                u1_imag[ux_idx(i, j, k)] = u0_imag[ux_idx(i, j, k)];
                u1_real[ux_idx(i, j, k)] = u0_real[ux_idx(i, j, k)];
            }
        }
    }

    //if(me_ == 0)
    //    printf("evv, %.8f, %.8f\n", u0_real[ux_idx(3, 3, 3)], u0_imag[ux_idx(3, 3, 3)]);

#undef ux_idx
}

void c_fft_1d_rev_(double* x1_real, double* x1_imag, double* x2_real, double* x2_imag)
{
    int fuyi = -1, san = 3, er = 2, yi  = 1;
    int i;
    unsigned long stcc, edcc;

    stcc = rpcc();
    host_flag[MPI_RANK] = me_;
//    host_flag[KERNEL_ACTION] = COUTLE_FFT_FLAG;
    host_flag[KERNEL_ACTION] = CP_EVV_FFT_FLAG;
    host_flag[GROUP_SIZE] = (dims(2, 3) * dims(3, 3)) / 32;
    host_flag[REMAIN_POINT] = (dims(2, 3) * dims(3, 3)) % 32;
//    host_flag[IN_PTR_REAL] = (long)x1_real;
//    host_flag[IN_PTR_IMAG] = (long)x1_imag;
    host_flag[IN_PTR_REAL] = (long)c_u0_real_;
    host_flag[IN_PTR_IMAG] = (long)c_u0_imag_;
    host_flag[OUT_PTR_REAL] = (long)x1_real;
    host_flag[OUT_PTR_IMAG] = (long)x1_imag;
    host_flag[IN_STRIDE] = dims(1, 3);
    host_flag[OUT_STRIDE] = dims(1, 3);
    host_flag[FFT_D1] = dims(1, 3);
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();

//    if(me_ == 0)
//    {
//        printf("MPE u1 check in : %.8f %.8f,   %.8f, %.8f\n", x1_real[510], x1_imag[510], x1_real[10 + host_flag[GROUP_SIZE] * dims(1, 3)], x1_imag[10 + host_flag[GROUP_SIZE] * dims(1, 3)]);
//    }


    //for(i = 0; i < )

//    c_cffts1(fuyi, dims(1, 3), dims(2, 3), dims(3, 3), x1_real, x1_imag, x1_real, x1_imag, scratch_real_, scratch_imag_);

//    if(me_ == 0)
//    {
//        printf("MPE u1 check out: %.8f %.8f,   %.8f, %.8f\n", x1_real[510], x1_imag[510], x1_real[10 + host_flag[GROUP_SIZE] * dims(1, 3)], x1_imag[10 + host_flag[GROUP_SIZE] * dims(1, 3)]);
//    }

    edcc = rpcc();
    mpe_cc_cur_[MPE_FT1] = edcc - stcc;

    stcc = rpcc();
    //transpose_x_yz_(&san, &er, x1, x2);
    c_transpose_x_yz(san, er, x1_real, x1_imag, x2_real, x2_imag);
    edcc = rpcc();
    //mpe_cc_cur_[MPE_TP1] = edcc - stcc;

    stcc = rpcc();
    host_flag[MPI_RANK] = me_;
    host_flag[KERNEL_ACTION] = DB_FFT_FLAG;
    host_flag[GROUP_SIZE] = (dims(2, 2) * dims(3, 2)) / 32;
    host_flag[REMAIN_POINT] = (dims(2, 2) * dims(3, 2)) % 32;
    host_flag[IN_PTR_REAL] = (long)x2_real;
    host_flag[IN_PTR_IMAG] = (long)x2_imag;
    host_flag[OUT_PTR_REAL] = (long)x1_real;//buffer
    host_flag[OUT_PTR_IMAG] = (long)x1_imag;//buffer, 2 -> 2,
    host_flag[IN_STRIDE] = dims(1, 2);
    host_flag[OUT_STRIDE] = dims(1, 2);
    host_flag[FFT_D1] = dims(1, 2);
    host_flag[FFT_D2] = dims(2, 2);
    host_flag[FFT_D3] = dims(3, 2);
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
    wait_slave_flag();
    wait_slave_flag();
    wait_slave_flag();
    //c_cffts1(fuyi, dims(1, 2), dims(2, 2), dims(3, 2), x2_real, x2_imag, x2_real, x2_imag, scratch_real_, scratch_imag_);
    edcc = rpcc();
    mpe_cc_cur_[MPE_FT2] = edcc - stcc;


//    stcc = rpcc();

//    host_flag[MPI_RANK] = me_;
//    host_flag[KERNEL_ACTION] = COUTLE_FFT_FLAG;
//    host_flag[GROUP_SIZE] = (dims(2, 2) * dims(3, 2)) / 32;
//    host_flag[REMAIN_POINT] = (dims(2, 2) * dims(3, 2)) % 32;
//    host_flag[IN_PTR_REAL] = (long)x2_real;
//    host_flag[IN_PTR_IMAG] = (long)x2_imag;
//    host_flag[OUT_PTR_REAL] = (long)x2_real;
//    host_flag[OUT_PTR_IMAG] = (long)x2_imag;
//    host_flag[IN_STRIDE] = dims(1, 2);
//    host_flag[OUT_STRIDE] = dims(1, 2);
//    host_flag[FFT_D1] = dims(1, 2);
//    asm volatile ("#nop":::"memory");
//    host_flag[0] = host_flag[0] + 1;
//    wait_slave_flag();
////    c_cffts1(fuyi, dims(1, 2), dims(2, 2), dims(3, 2), x2_real, x2_imag, x2_real, x2_imag, scratch_real_, scratch_imag_);
//    edcc = rpcc();
//    mpe_cc_cur_[MPE_FT2] = edcc - stcc;

//    local_transpose_x_y(x2_real, x1_real, dims(1, 2), dims(2, 2), dims(3, 2));
//    local_transpose_x_y(x2_imag, x1_imag, dims(1, 2), dims(2, 2), dims(3, 2));


//    stcc = rpcc();
////    c_cffts1(fuyi, dims(1, 1), dims(2, 1), dims(3, 1), x1_real, x1_imag, x2_real, x2_imag, scratch_real_, scratch_imag_);
//    host_flag[MPI_RANK] = me_;
//    host_flag[KERNEL_ACTION] = COUTLE_FFT_FLAG;
//    host_flag[GROUP_SIZE] = (dims(2, 1) * dims(3, 1)) / 32;
//    host_flag[REMAIN_POINT] = (dims(2, 1) * dims(3, 1)) % 32;
//    host_flag[IN_PTR_REAL] = (long)x1_real;
//    host_flag[IN_PTR_IMAG] = (long)x1_imag;
//    host_flag[OUT_PTR_REAL] = (long)x2_real;
//    host_flag[OUT_PTR_IMAG] = (long)x2_imag;
//    host_flag[IN_STRIDE] = dims(1, 1);
//    host_flag[OUT_STRIDE] = dims(1, 1);
//    host_flag[FFT_D1] = dims(1, 1);
//    asm volatile ("#nop":::"memory");
//    host_flag[0] = host_flag[0] + 1;
//    wait_slave_flag();
//    edcc = rpcc();
//    mpe_cc_cur_[MPE_FT3] = edcc - stcc;
}

//int __ofs510;

void c_cffts1(int is, int d1, int d2, int d3, double* x_real, double* x_imag, double* xout_real, double* xout_imag, double* y_real, double* y_imag)
{
#define x_idx(a, b, c) ((a - 1) + ((b - 1) + (c - 1) * d2) * d1)
#define y_idx(a, b, c) ((a - 1) + ((b - 1) + (c - 1) * d1) * fftblockpad_)
    //double complex x(d1,d2,d3)
    //double complex xout(d1,d2,d3)
    //double complex y(fftblockpad, d1, 2)
    int logd1, i, j, k, jj;

    logd1 = ilog2(d1);

    for(k = 1; k <= d3; k++)
    {
        for(jj = 0; jj <= d2 - fftblock_; jj += fftblock_)
        {
            for(j = 1; j <= fftblock_; j++)
            {
                for(i = 1; i <= d1; i++)
                {
                    y_real[y_idx(j,i,1)] = x_real[x_idx(i,j+jj,k)];
                    y_imag[y_idx(j,i,1)] = x_imag[x_idx(i,j+jj,k)];
//                    if(x_idx(i,j+jj,k) == 510)
//                        __ofs510 = y_idx(j,i,1);
                }
            }

            //cfftz_(&is, &logd1, &d1, y, &y[y_idx(1, 1, 2)]);
            c_cfftz(is, logd1, d1, y_real, y_imag, &y_real[y_idx(1, 1, 2)], &y_imag[y_idx(1, 1, 2)]);

            for(j = 1; j <= fftblock_; j++)
            {
                for(i = 1; i <= d1; i++)
                {
                    xout_real[x_idx(i,j+jj,k)] = y_real[y_idx(j,i,1)];
                    xout_imag[x_idx(i,j+jj,k)] = y_imag[y_idx(j,i,1)];
                }
            }
        }
    }
#undef x_idx
#undef y_idx
}

void c_cffts2(int is, int d1, int d2, int d3, double* x_real, double* x_imag, double* xout_real, double* xout_imag, double* y_real, double* y_imag)
{
#define x_idx(a, b, c) ((a - 1) + ((b - 1) + (c - 1) * d2) * d1)
#define y_idx(a, b, c) ((a - 1) + ((b - 1) + (c - 1) * d2) * fftblockpad_)
    //double complex x(d1,d2,d3)
    //double complex xout(d1,d2,d3)
    //double complex y(fftblockpad, d2, 2)

    int logd2, i, j, k, ii;

    logd2 = ilog2(d2);

    for(k = 1; k <= d3; k++)
    {
        for(ii = 0; ii <= d1 - fftblock_; ii += fftblock_)
        {
            for(j = 1; j <= d2; j++)
            {
                for(i = 1; i <= fftblock_; i++)
                {
                    y_real[y_idx(i,j,1)] = x_real[x_idx(i+ii,j,k)];
                    y_imag[y_idx(i,j,1)] = x_imag[x_idx(i+ii,j,k)];
                }
            }

            //cfftz_(&is, &logd2, &d2, y, &y[y_idx(1, 1, 2)]);
            c_cfftz(is, logd2, d2, y_real, y_imag, &y_real[y_idx(1, 1, 2)], &y_imag[y_idx(1, 1, 2)]);

            for(j = 1; j <= d2; j++)
            {
                for(i = 1; i <= fftblock_; i++)
                {
                    xout_real[x_idx(i+ii,j,k)] = y_real[y_idx(i,j,1)];
                    xout_imag[x_idx(i+ii,j,k)] = y_imag[y_idx(i,j,1)];
                }
            }
        }
    }
#undef x_idx
#undef y_idx
}

void c_cfftz(int is, int m, int n, double* x_real, double* x_imag, double* y_real, double* y_imag)
{
#define x_idx(a, b) ((a - 1) + (b - 1) * fftblockpad_)
#define y_idx(a, b) ((a - 1) + (b - 1) * fftblockpad_)
    //dimension x(fftblockpad,n), y(fftblockpad,n)
    int i,j,l,lp1,mx;
    //mx = u_[0];
    //... check

//    if(me_ == 0)
//    {
////        printf("MPE u1 check allin_cfftz : %.8f %.8f\n", x_real[__ofs510], x_imag[__ofs510]);
//        for(i = 0; i < 512; i++)
//        {
//            printf("MPE u1 check allin_cfftz i = %d: %.8f %.8f %.8f %.8f\n", i, x_real[i], x_imag[i], y_real[i], y_imag[i]);
//        }
//    }

    for(l = 1; l <= m; l += 2)
    {
        //fftz2_(&is, &l, &m, &n, &fftblock_, &fftblockpad_, u_, x, y);
        c_fftz2(is, l, m, n, fftblock_, fftblockpad_, c_u_real_, c_u_imag_, x_real, x_imag, y_real, y_imag);
//        if(me_ == 0 && l == 1)
//        {
////            printf("MPE u1 check l = %d y : %.8f %.8f\n", l, y_real[__ofs510], y_imag[__ofs510]);
//            for(i = 0; i < 512; i++)
//            {
//                printf("MPE u1 check l i = %d: %.8f %.8f %.8f %.8f\n", i, x_real[i], x_imag[i], y_real[i], y_imag[i]);
//            }
//        }

        if(l == m)
        {
            for(j = 1; j <= n; j++)
            {
                for(i = 1; i <= fftblock_; i++)
                {
                    x_real[x_idx(i, j)] = y_real[y_idx(i, j)];
                    x_imag[x_idx(i, j)] = y_imag[y_idx(i, j)];
                }
            }
            break;
        }

        lp1 = l + 1;
        //fftz2_dbg_(&is, &lp1, &m, &n, &fftblock_, &fftblockpad_, u_, y, x);
        c_fftz2(is, lp1, m, n, fftblock_, fftblockpad_, c_u_real_, c_u_imag_, y_real, y_imag, x_real, x_imag);
        //exit(999);
//        if(me_ == 0 && l == 1)
//        {
////            printf("MPE u1 check l = %d x : %.8f %.8f\n", l, x_real[__ofs510], x_imag[__ofs510]);
//            for(i = 0; i < 512; i++)
//            {
//                printf("MPE u1 check lp1 i = %d: %.8f %.8f %.8f %.8f\n", i, x_real[i], x_imag[i], y_real[i], y_imag[i]);
//            }
//        }
    }

//    if(me_ == 0)
//        printf("__ofs510 = %d\n", __ofs510);


//    exit(23232);

#undef x_idx
#undef y_idx
}

void c_fftz2(int is, int l, int m, int n, int ny, int ny1,
             double* u_real, double* u_imag, double* x_real, double* x_imag, double* y_real, double* y_imag)
{
#define x_idx(a, b) ((a - 1) + (b - 1) * ny1)
#define y_idx(a, b) ((a - 1) + (b - 1) * ny1)
#define u_idx(a) (a - 1)
    //dimension u(n), x(ny1,n), y(ny1,n)

    double u1_real,x11_real,x21_real;
    double u1_imag,x11_imag,x21_imag;

    int k,n1,li,lj,lk,ku,i,j,i11,i12,i21,i22;

    n1 = n / 2;
    lk = int_pow(2, l - 1);
    li = int_pow(2, m - l);
    lj = 2 * lk;
    ku = li + 1;

    //printf("%d %d %d %d %d %d\n", is, l, m, n, ny, ny1);//-1 2 5 32 16 18
    //printf("%d %d %d %d %d\n", n1, lk, li, lj, ku);

    for(i = 0; i <= li - 1; i++)
    {
        i11 = i * lk + 1;
        i12 = i11 + n1;
        i21 = i * lj + 1;
        i22 = i21 + lk;

        u1_real = u_real[u_idx(ku + i)];
        //u1_imag = u_imag[u_idx(ku + i)];
        u1_imag = -u_imag[u_idx(ku + i)];//you need to do this on mpe

        for(k = 0; k <= lk - 1; k++)
        {
            for(j = 1; j <= ny; j++)
            {
                x11_real = x_real[x_idx(j, i11 + k)];
                x11_imag = x_imag[x_idx(j, i11 + k)];
                x21_real = x_real[x_idx(j, i12 + k)];
                x21_imag = x_imag[x_idx(j, i12 + k)];

                y_real[y_idx(j, i21 + k)] = x11_real + x21_real;
                y_imag[y_idx(j, i21 + k)] = x11_imag + x21_imag;

                //minus
                x11_real -= x21_real;
                x11_imag -= x21_imag;//cpe_minus_dc_dc(x11, x21)

                y_real[y_idx(j, i22 + k)] = u1_real * x11_real - u1_imag * x11_imag;
                y_imag[y_idx(j, i22 + k)] = u1_real * x11_imag + u1_imag * x11_real;

//                if(me_ == 0 && l == 1)
//                {
//                    printf("MPE fftz2 dbg: %d %d %d %d %d, %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
//                           u_idx(ku + i), x_idx(j, i11 + k), x_idx(j, i12 + k), y_idx(j, i21 + k), y_idx(j, i22 + k),
//                           u1_real, u1_imag, x11_real, x11_imag, x21_real, x21_imag,
//                           y_real[y_idx(j, i21 + k)], y_imag[y_idx(j, i21 + k)], y_real[y_idx(j, i22 + k)], y_imag[y_idx(j, i22 + k)]
//                            );
//                }
            }
        }
    }

#undef u_idx
#undef x_idx
#undef y_idx
}

void c_transpose_x_yz(int l1, int l2, double* xin_real, double* xin_imag, double* xout_real, double* xout_imag)
{
    int dim_mul = dims(2, l1) * dims(3, l1);
    unsigned long stcc, edcc;

//    stcc = rpcc();
//    //transpose2_local_(&dims(1, l1), &dim_mul, xin, xout);
//    c_transpose2_local(dims(1, l1), dim_mul, xin_real, xin_imag, xout_real, xout_imag);
//    edcc = rpcc();
//    mpe_cc_cur_[MPE_T2L] = edcc - stcc;

    stcc = rpcc();
    //transpose2_global_(xout, xin);
    c_transpose2_global(xout_real, xout_imag, xin_real, xin_imag);
    edcc = rpcc();
    mpe_cc_cur_[MPE_T2G] = edcc - stcc;

//    stcc = rpcc();
//    //transpose2_finish_(&dims(1, l1), &dim_mul, xin, xout);
//    c_transpose2_finish(dims(1, l1), dim_mul, xin_real, xin_imag, xout_real, xout_imag);
//    edcc = rpcc();
//    mpe_cc_cur_[MPE_T2F] = edcc - stcc;
}

//void c_transpose2_local(int n1, int n2, double* xin_real, double* xin_imag, double* xout_real, double* xout_imag)
//{
//    host_flag[MPI_RANK] = me_;
//    host_flag[KERNEL_ACTION] = TP_LOCAL_FLAG;
//    host_flag[IN_PTR_REAL] = (long)xin_real;
//    host_flag[IN_PTR_IMAG] = (long)xin_imag;
//    host_flag[OUT_PTR_REAL] = (long)xout_real;
//    host_flag[OUT_PTR_IMAG] = (long)xout_imag;
//    host_flag[FFT_D1] = dims(1, 3);
//    host_flag[FFT_D2] = dims(2, 3);
//    host_flag[FFT_D3] = dims(3, 3);
//    asm volatile ("#nop":::"memory");
//    host_flag[0] = host_flag[0] + 1;
//    wait_slave_flag();
//    wait_slave_flag();

////#define xin_idx(a, b) (a - 1 + (b - 1) * n1)
////#define xout_idx(a, b) (a - 1 + (b - 1) * n2)
////    //double complex xin(n1, n2), xout(n2, n1)
////    int i, j, ii, jj;

////    //在原来的fortran代码中，有一些非常基础的矩阵转置优化。在这里就不移植了
////    for(i = 1; i <= n1; i++)
////    {
////        for(j = 1; j <= n2; j++)
////        {
////            xout_real[xout_idx(j, i)] = xin_real[xin_idx(i, j)];
////            xout_imag[xout_idx(j, i)] = xin_imag[xin_idx(i, j)];
////        }
////    }
////#undef xin_idx
////#undef xout_idx
//}

void c_transpose2_global(double* xin_real, double* xin_imag, double* xout_real, double* xout_imag)
{
    //call mpi_alltoall(xin, ntdivnp/np, dc_type,
    //     >                  xout, ntdivnp/np, dc_type,
    //     >                  commslice1, ierr)

    MPI_Alltoall(xin_real, ntdivnp_ / np_, MPI_DOUBLE,
                 xout_real, ntdivnp_ / np_, MPI_DOUBLE,
                 commslice1_);
    MPI_Alltoall(xin_imag, ntdivnp_ / np_, MPI_DOUBLE,
                 xout_imag, ntdivnp_ / np_, MPI_DOUBLE,
                 commslice1_);
}

//void c_transpose2_finish(int n1, int n2, double* xin_real, double* xin_imag, double* xout_real, double* xout_imag)
//{
//    host_flag[MPI_RANK] = me_;
//    host_flag[KERNEL_ACTION] = TP_FINISH_FLAG;
//    host_flag[IN_PTR_REAL] = (long)xin_real;
//    host_flag[IN_PTR_IMAG] = (long)xin_imag;
//    host_flag[OUT_PTR_REAL] = (long)xout_real;
//    host_flag[OUT_PTR_IMAG] = (long)xout_imag;
////    host_flag[FFT_D1] = dims(1, 3);//512
////    host_flag[FFT_D2] = dims(3, 3);//8
////    host_flag[FFT_D3] = np2_;//64
//    host_flag[FFT_D1] = n2;
//    host_flag[FFT_D2] = n1 / np2_;
//    host_flag[FFT_D3] = np2_;
//    host_flag[IN_STRIDE] = (long)dims(1, 1);//part len
//    asm volatile ("#nop":::"memory");
//    host_flag[0] = host_flag[0] + 1;
//    wait_slave_flag();
//    wait_slave_flag();

////    int xin_d1 = n2, xin_d2 = n1 / np2_, xin_d3 = np2_;
////    int xout_d1 = n2 * np2_, xout_d2 = n1 / np2_;
////#define xin_idx(a, b, c) (a - 1 + (b - 1 + c * xin_d2) * xin_d1)
////#define xout_idx(a, b) (a - 1 + (b - 1) * xout_d1)
////    //double complex xin(n2, n1/np2, 0:np2-1), xout(n2*np2, n1/np2)
////    int ioff, i, j, p;

////    for(p = 0; p <= np2_ - 1; p++)
////    {
////        ioff = p * n2;
////        for(j = 1; j <= xin_d2; j++)
////        {
////            for(i = 1; i <= n2; i++)
////            {
////                xout_real[xout_idx(i + ioff, j)] = xin_real[xin_idx(i, j, p)];
////                xout_imag[xout_idx(i + ioff, j)] = xin_imag[xin_idx(i, j, p)];
////            }
////        }
////    }
////#undef xin_idx
////#undef xout_idx
//}



void checksum_offset_get()
{
    int j, q, r, s;
    int core_x, core_y, tmpx, tmpy;
    int cpe_cbw;

    chk_core_cnt_ = (int*)malloc(64 * sizeof(int));
    chk_core_stp_ = (int*)malloc(64 * sizeof(int));
    cpe_cbw = nx_ / 8;

    for(j = 0; j < 64; j++)
        chk_core_cnt_[j] = 0;

    chk_cnt_ = 0;
    for(j = 1; j <= 1024; j++)
    {
        q = j % nx_ + 1;
        //q = mod(j, nx)+1
        r = (3 * j) % ny_ + 1;
        //r = mod(3*j,ny)+1
        s = (5 * j) % nz_ + 1;
        //s = mod(5*j,nz)+1

        if(q >= xstart_[0] && q <= xend_[0] &&
                r >= ystart_[0] && r <= yend_[0] &&
                s >= zstart_[0] && s <= zend_[0])
        {
            chk_cnt_++;

            //z, x, y
            //only x and y is needed
            //core_y <=> y
            //core_x <=> x
            tmpx = q-xstart_[0];
            tmpy = r-ystart_[0];
            tmpx /= cpe_cbw;
            tmpy /= cpe_cbw;

//            core_x = tmpx;//???
//            core_y = tmpy;//???

            //fix for lyz
            core_x = tmpy;//???
            core_y = tmpx;//???

            chk_core_cnt_[8 * core_y + core_x]++;
        }
    }

    chk_xs_ = (int*)malloc(sizeof(int) * chk_cnt_);
    chk_ys_ = (int*)malloc(sizeof(int) * chk_cnt_);
    chk_zs_ = (int*)malloc(sizeof(int) * chk_cnt_);

    chk_core_stp_[0] = 0;
    for(j = 1; j < 64; j++)
    {
        chk_core_stp_[j] = chk_core_stp_[j - 1] + chk_core_cnt_[j - 1];
        //printf("core = %d, stp = %d, cnt = %d\n", j, chk_core_stp_[j], chk_core_cnt_[j]);
    }

    int tmp_stp;
    for(j = 1; j <= 1024; j++)
    {
        q = j % nx_ + 1;
        //q = mod(j, nx)+1
        r = (3 * j) % ny_ + 1;
        //r = mod(3*j,ny)+1
        s = (5 * j) % nz_ + 1;
        //s = mod(5*j,nz)+1

        if(q >= xstart_[0] && q <= xend_[0] &&
                r >= ystart_[0] && r <= yend_[0] &&
                s >= zstart_[0] && s <= zend_[0])
        {
            //(chk_xs_[iter], chk_ys_[iter], chk_zs_[iter]) = (q-xstart_[0]+1,r-ystart_[0]+1,s-zstart_[0]+1);

            tmpx = q-xstart_[0];
            tmpy = r-ystart_[0];
            tmpx /= cpe_cbw;
            tmpy /= cpe_cbw;
//            core_x = tmpx;//???
//            core_y = tmpy;//???

            //fix for lyz
            core_x = tmpy;//???
            core_y = tmpx;//???
            tmp_stp = chk_core_stp_[8 * core_y + core_x];

            chk_xs_[tmp_stp] = q-xstart_[0];
            chk_ys_[tmp_stp] = r-ystart_[0];
            chk_zs_[tmp_stp] = s-zstart_[0];
            //yuan ban de x y z zuo biao

            //printf("j = %d, stp = %d, x = %d, y = %d, z = %d, x = %d, y = %d, z = %d, core = %d\n",
            //       j, tmp_stp, q-xstart_[0], r-ystart_[0], s-zstart_[0],
            //        chk_xs_[tmp_stp], chk_ys_[tmp_stp], chk_zs_[tmp_stp],
            //        8 * core_y + core_x);

            chk_core_stp_[8 * core_y + core_x]++;//
        }
    }

    //重新计算一次，刚才临时修改了
    chk_core_stp_[0] = 0;
    for(j = 1; j < 64; j++)
        chk_core_stp_[j] = chk_core_stp_[j - 1] + chk_core_cnt_[j - 1];


    buf_chk_ = (dfc*)malloc(sizeof(dfc) * 64 * niter_);
    for(j = 0; j < 64 * niter_; j++)
        buf_chk_[j].real = buf_chk_[j].imag = 0.0;
    //每个从核都是连续访问的。不使用间杂访问的方法

    //44
    //for(j = 672; j < 672 + 32; j++)
    //{
    //    printf("44, j = %d, x = %d, y  =%d, z = %d\n", j, chk_xs_[j], chk_ys_[j], chk_zs_[j]);
    //}

    if(me_ == 0)
        printf("chk cnt = %d\n", chk_cnt_);
    //then you need to sort them
    //or do it on cpe
}

void c_checksum(int i, double* u1_real, double* u1_imag, int d1, int d2, int d3)
{
#define ux_idx(a, b, c) ((a - 1) + ((b - 1) + (c - 1) * d2) * d1)
    //double complex u1(d1, d2, d3)
    int j, q, r, s;
    struct dcomplex chk, allchk;

    chk.real = 0.0;
    chk.imag = 0.0;

    for(j = 1; j <= 1024; j++)
    {
        q = j % nx_ + 1;
        //q = mod(j, nx)+1
        r = (3 * j) % ny_ + 1;
        //r = mod(3*j,ny)+1
        s = (5 * j) % nz_ + 1;
        //s = mod(5*j,nz)+1

        if(q >= xstart_[0] && q <= xend_[0] &&
                r >= ystart_[0] && r <= yend_[0] &&
                s >= zstart_[0] && s <= zend_[0])
        {
            //chk = add_dc_dc(chk, u1[ux_idx(q-xstart_[0]+1,r-ystart_[0]+1,s-zstart_[0]+1)]);
            chk.real += u1_real[ux_idx(q-xstart_[0]+1,r-ystart_[0]+1,s-zstart_[0]+1)];
            chk.imag += u1_imag[ux_idx(q-xstart_[0]+1,r-ystart_[0]+1,s-zstart_[0]+1)];
            //chk=chk+u1(q-xstart(1)+1,r-ystart(1)+1,s-zstart(1)+1)
        }
    }

    chk.real /= ntotal_f_;
    chk.imag /= ntotal_f_;
    //call MPI_Reduce(chk, allchk, 1, dc_type, MPI_SUM,
         //>                0, MPI_COMM_WORLD, ierr)

    MPI_Reduce(&chk, &allchk, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    f_sums_[i] = allchk;

    if(me_ == 0)
    {
        printf("T =    %d     Checksum =    %.8f    %.8f\n", i, allchk.real, allchk.imag);
        //mmpmmpmmp_(&i, &allchk);
    }

#undef ux_idx
}

void c_checksum_d1d2_change(int i, double* u1_real, double* u1_imag, int d1, int d2, int d3)
{
#define ux_idx(a, b, c) ((b - 1) + ((a - 1) + (c - 1) * d2) * d1)
    //double complex u1(d1, d2, d3)
    int j, q, r, s;
    struct dcomplex chk, allchk;

    chk.real = 0.0;
    chk.imag = 0.0;

    for(j = 1; j <= 1024; j++)
    {
        q = j % nx_ + 1;
        //q = mod(j, nx)+1
        r = (3 * j) % ny_ + 1;
        //r = mod(3*j,ny)+1
        s = (5 * j) % nz_ + 1;
        //s = mod(5*j,nz)+1

        if(q >= xstart_[0] && q <= xend_[0] &&
                r >= ystart_[0] && r <= yend_[0] &&
                s >= zstart_[0] && s <= zend_[0])
        {
            //chk = add_dc_dc(chk, u1[ux_idx(q-xstart_[0]+1,r-ystart_[0]+1,s-zstart_[0]+1)]);
            chk.real += u1_real[ux_idx(q-xstart_[0]+1,r-ystart_[0]+1,s-zstart_[0]+1)];
            chk.imag += u1_imag[ux_idx(q-xstart_[0]+1,r-ystart_[0]+1,s-zstart_[0]+1)];
            //chk=chk+u1(q-xstart(1)+1,r-ystart(1)+1,s-zstart(1)+1)
        }
    }

    chk.real /= ntotal_f_;
    chk.imag /= ntotal_f_;
    //call MPI_Reduce(chk, allchk, 1, dc_type, MPI_SUM,
         //>                0, MPI_COMM_WORLD, ierr)

    MPI_Reduce(&chk, &allchk, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    f_sums_[i] = allchk;

    if(me_ == 0)
    {
        printf("T =    %d     Checksum =    %.8f    %.8f\n", i, allchk.real, allchk.imag);
        //mmpmmpmmp_(&i, &allchk);
    }

#undef ux_idx
}

void c_checksum_merge(int i)
{
    struct dcomplex chk, allchk;
    int j, q, r, s;

    chk.real = 0.0;
    chk.imag = 0.0;

    for(j = 0; j < 64; j++)
    {
        if(chk_core_cnt_[j] > 0)
        {
            chk.real += buf_chk_[j * niter_ + i - 1].real;
            chk.imag += buf_chk_[j * niter_ + i - 1].imag;
        }
    }

    chk.real /= ntotal_f_;
    chk.imag /= ntotal_f_;

    MPI_Reduce(&chk, &allchk, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);

    f_sums_[i].real = allchk.real;
    f_sums_[i].imag = allchk.imag;
    if(me_ == 0)
        printf("T =    %d     Checksum =    %.8f    %.8f\n", i, allchk.real, allchk.imag);
}
