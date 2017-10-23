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
    param.u0_real = c_u0_real_;
    param.u1_real = c_u1_real_;
    param.u2_real = c_u2_real_;
    param.u_real = c_u_real_;
    param.u0_imag = c_u0_imag_;
    param.u1_imag = c_u1_imag_;
    param.u2_imag = c_u2_imag_;
    param.u_imag = c_u_imag_;
    param.twiddle = c_twiddle_;
    param.buf_chk = buf_chk_;
    param.chk_cnt = chk_core_cnt_;
    param.chk_stp = chk_core_stp_;
    param.chk_xs = chk_xs_;
    param.chk_ys = chk_ys_;
    param.chk_zs = chk_zs_;

    // if(param.my_id == 0)
    // {
    //     printf("*** MESSAGE FROM MPE ***\n");
    //     printf("iter = %d, nx = %d, ny = %d, nz = %d, npc = %d\n", param.iter, param.nx, param.ny, param.nz, param.npc);
    //     printf("ptrs: %x, %x, %x, %x, %x, %x, %x, %x, %x, %x, %x\n%x %x %x %x %x %x\n",
    //            param.host_flag, param.slave_flag, param.u0_real, param.u1_real, param.u2_real, param.u_real, param.u0_imag, param.u1_imag, param.u2_imag, param.u_imag, param.twiddle,
    //            param.chk_cnt, param.chk_stp, param.chk_xs, param.chk_ys, param.chk_zs, param.buf_chk);
    //     printf("*** MESSAGE FROM MPE ***\n");
    // }

    athread_init();

//    if(me_ == 0)
//    {
//        for(i = 0; i < nx_; i++)
//        {
//            printf("MPE c_u_[%d] = %.8f %.8f\n", i, c_u_real_[i], c_u_imag_[i]);
//        }
//    }

    athread_spawn(cpe_athread_daemon, (void*)&param);

    //放flag，让从核进入循环

    host_flag[MPI_RANK] = me_;
    host_flag[KERNEL_ACTION] = ALLIN_FFT_FLAG;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;


    //the 万恶的 iteration
    for(iter = 1; iter <= niter_; iter++)
    {
        //c_evolve(c_u0_, c_u2_, c_twiddle_, dims(1, 1), dims(2, 1), dims(3, 1));
        //if(me_ == 0)
        //    printf("MPE finished evv\n");
        //fflush(stdout);
        //asm volatile ("#nop":::"memory");
        //host_flag[0] = host_flag[0] + 1;

        // // 等待第一次转置flag
        // wait_slave_flag();
        // total_dt_transpose_x_y(c_u1_real_, c_u2_real_);
        // total_dt_transpose_x_y(c_u1_imag_, c_u2_imag_);
        // asm volatile ("#nop":::"memory");
        // host_flag[0] = host_flag[0] + 1;

        // 等待第二次转置flag
//        wait_slave_flag();
//        total_dt_transpose_x_z(c_u1_real_, c_u2_real_);
//        total_dt_transpose_x_z(c_u1_imag_, c_u2_imag_);
//        asm volatile ("#nop":::"memory");
//        host_flag[0] = host_flag[0] + 1;

        

        //等待checksum的flag
        wait_slave_flag();
        //if(me_ == 0)
        //    printf("MPE chk\n");
        //c_checksum_rotated(iter, c_u1_, dims(1, 1), dims(2, 1), dims(3, 1));
        c_checksum_merge(iter);
        //asm volatile ("#nop":::"memory");
        //host_flag[0] = host_flag[0] + 1;
    }

    host_flag[MPI_RANK] = me_;
    host_flag[KERNEL_ACTION] = EXIT_FLAG;
    host_flag[GROUP_SIZE] = -1;
    asm volatile ("#nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
    wait_slave_flag();

    athread_join();

    // printf("%d, %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld\n",
           // me_, slave_flag[1], slave_flag[2], slave_flag[3], slave_flag[4], slave_flag[5], slave_flag[6], slave_flag[7], slave_flag[8]);

}



// void c_transpose_x_yz(int l1, int l2, my_fc* xin, my_fc* xout)
// {
//     int dim_mul = dims(2, l1) * dims(3, l1);
//     unsigned long stcc, edcc;

//     stcc = rpcc();
//     transpose2_local_(&dims(1, l1), &dim_mul, xin, xout);
//     edcc = rpcc();
//     mpe_cc_cur_[MPE_T2L] = edcc - stcc;

//     stcc = rpcc();
//     c_transpose2_global(xout, xin);
//     edcc = rpcc();
//     mpe_cc_cur_[MPE_T2G] = edcc - stcc;

//     stcc = rpcc();
//     transpose2_finish_(&dims(1, l1), &dim_mul, xin, xout);
//     edcc = rpcc();
//     mpe_cc_cur_[MPE_T2F] = edcc - stcc;
// }

// void c_transpose2_local(int n1, int n2, my_fc* xin, my_fc* xout)
// {
// #define xin_idx(a, b) (a - 1 + (b - 1) * n1)
// #define xout_idx(a, b) (a - 1 + (b - 1) * n2)
//     //double complex xin(n1, n2), xout(n2, n1)
//     int i, j, ii, jj;

//     //在原来的fortran代码中，有一些非常基础的矩阵转置优化。在这里就不移植了
//     for(i = 1; i <= n1; i++)
//     {
//         for(j = 1; j <= n2; j++)
//         {
//             xout[xout_idx(j, i)] = xin[xin_idx(i, j)];
//         }
//     }
// #undef xin_idx
// #undef xout_idx
// }

// void c_transpose2_global(my_fc* xin, my_fc* xout)
// {
//     //call mpi_alltoall(xin, ntdivnp/np, dc_type,
//     //     >                  xout, ntdivnp/np, dc_type,
//     //     >                  commslice1, ierr)

//     MPI_Alltoall(xin, ntdivnp_ / np_, my_mpift,
//                  xout, ntdivnp_ / np_, my_mpift,
//                  commslice1_);
// }

// void c_transpose2_finish(int n1, int n2, my_fc* xin, my_fc* xout)
// {
//     int xin_d1 = n2, xin_d2 = n1 / np2_, xin_d3 = np2_;
//     int xout_d1 = n2 * np2_, xout_d2 = n1 / np2_;
// #define xin_idx(a, b, c) (a - 1 + (b - 1 + c * xin_d2) * xin_d1)
// #define xout_idx(a, b) (a - 1 + (b - 1) * xout_d1)
//     //double complex xin(n2, n1/np2, 0:np2-1), xout(n2*np2, n1/np2)
//     int ioff, i, j, p;

//     for(p = 0; p <= np2_ - 1; p++)
//     {
//         ioff = p * n2;
//         for(j = 1; j <= xin_d2; j++)
//         {
//             for(i = 1; i <= n2; i++)
//             {
//                 xout[xout_idx(i + ioff, j)] = xin[xin_idx(i, j, p)];
//             }
//         }
//     }
// #undef xin_idx
// #undef xout_idx
// }



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

    // if(me_ == 0)
        // printf("chk cnt = %d\n", chk_cnt_);
    //then you need to sort them
    //or do it on cpe
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
    // if(me_ == 0)
        // printf("T =    %d     Checksum =    %.8f    %.8f\n", i, allchk.real, allchk.imag);
}
