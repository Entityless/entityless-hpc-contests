#include "../c_header/mpe_weno.h"

void exchange_boundary()
{
    //copy_back_boundary_u();
    double* tmpdp;
    tmpdp = c_last_u;
    c_last_u = c_cur_u;
    c_cur_u = tmpdp;
    //just get the cur_u
    //and exchange its boundary
    //next step, it will be the last_u

    exchange_boundary_x_standard_(c_last_u, param_weno.iperiodic + 0);
    //weno is cross-x computation
    //so CFD needs yz exchange, weno need not.
    //exchange_boundary_y_standard_(c_last_u, param_weno.iperiodic + 1);
    //exchange_boundary_z_standard_(c_last_u, param_weno.iperiodic + 2);
    //copy_from_boundary_u();
}

void weno_c_core_()
{
    //check_u_(c_cur_u);
    //check_u_(c_last_u);
    exchange_boundary();
    //check_u_(c_cur_u);
    //check_u_(c_last_u);

    int k, j, i;

    param_weno.host_flag  = &host_flag[0];
    param_weno.slave_flag  = &slave_flag[0];

    for(i = 0; i < FLAG_SIZE; i++)
    {
        host_flag[i] = 0;
        slave_flag[i] = 0;
        local_cc[i] = 0;
    }

    flag_to_wait = 1;

    athread_init();
    athread_spawn(cpe_athread_daemon, &param_weno);

    volatile long stcc = 0, edcc = 0;

    //防止浮点数eps差异导致的超次迭代
    int iter_max;
    double iter_time_double;
    int certain_time;
    if(param_weno.end_time > 0)//具有具体迭代次数
    {
        certain_time = 1;
        iter_time_double = param_weno.end_time / param_weno.dt;
        iter_max = iter_time_double + 0.01;//eps for integer padding, 0.999 -> 1.009, 1.001 -> 1.011, --> 1
    }
    else
    {
        certain_time = 0;
    }

    iter_step = 0;


    int nx = param_weno.nx;
    int ny = param_weno.ny;
    int nz = param_weno.nz;

    int f_dim1 = nx + 2 * param_weno.lap, f_dim2 = ny + 2 * param_weno.lap, f_dim3 = nz + 2 * param_weno.lap;
    int f_s1 = 1 - param_weno.lap, f_s2 = 1 - param_weno.lap, f_s3 = 1 - param_weno.lap;
    int fx_dim1 = nx, fx_dim2 = ny, fx_dim3 = nz;
    int fx_s1 = 1, fx_s2 = 1, fx_s3 = 1;

    fx_dim3 = fx_dim2 * fx_dim1;
    fx_dim2 = fx_dim1;
    f_dim3 = f_dim2 * f_dim1;
    f_dim2 = f_dim1;

    double wend, wstar, wstar0;

    //wstar0 = MPI_Wtime();
    wtime_util_(&wstar0);
    wstar = wstar0;

    ////////stcc = rpcc();

    while(1)
    {
        iter_step++;

        //check_u_(c_cur_u);
        //check_u_(c_last_u);
        //weno7_c_halo_();
        weno7_c_inner_();
        //check_u_(c_cur_u);

        ////////edcc = rpcc();
        local_cc[MPF_INNER] += edcc - stcc;
        stcc = edcc;

        exchange_boundary();
        ////////edcc = rpcc();
        local_cc[MPF_EXCHANGE] += edcc - stcc;
        stcc = edcc;

        param_weno.istep[0]++;
        param_weno.tt[0] += param_weno.dt;

        wtime_util_(&wend);

        if((*param_weno.istep) % param_weno.k_step_show == 0)
        {
            if(param_weno.my_id == 0)
            {
                //wend = MPI_Wtime();
                printf("-------Istep= %d\n", param_weno.istep[0]);
                printf("CPU (wall) time for this step is %.8f\n", wend - wstar);
                printf("Total CPU (wall) time is %.8f\n", wend - wstar0);

                wstar = wend;
            }
        }

        if((*param_weno.istep) % param_weno.k_step_save == 0)
        {
            //memcpy(param_weno.u, c_u, sizeof(double) * sz_c_u);

            ocfd_save_util_(c_last_u);

            ////////edcc = rpcc();
            local_cc[MPF_SAVE_DATA] += edcc - stcc;
            stcc = edcc;

            if(param_weno.end_time < 0)
                break;
        }

        if(certain_time != 0 && iter_max == iter_step)//防止误差导致的迭代
        {
            break;
        }

        if((*param_weno.tt) < param_weno.end_time || param_weno.end_time < 0)
            continue;

        break;
    }

    //finalize
    //{
        host_flag[MPI_RANK] = param_weno.my_id;
        host_flag[KERNEL_ACTION] = EXIT_FLAG;
        host_flag[GROUP_SIZE] = -1;
        asm volatile ("#nop":::"memory");
        host_flag[0] = host_flag[0] + 1;
        wait_slave_flag();
        wait_slave_flag();

        athread_join();
    //}

    //time print
    //{
        double out_time[FLAG_SIZE];
        double local_time[FLAG_SIZE];
        for(i = 1; i < 32; i++)
        {
            out_time[i] = slave_flag[i] * 1.0 / CCPS;
            local_time[i] = local_cc[i] * 1.0 / CCPS;
        }

        asm volatile ("#nop":::"memory");

        printf("%d, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n",
               param_weno.my_id, wend - wstar0, out_time[1] - local_time[6], out_time[2], out_time[3], out_time[4], out_time[5], out_time[6],
                local_time[1], local_time[2], local_time[3], local_time[4], local_time[5], local_time[6]);
    //}

        //if(param_weno.my_id == 0)
        //    printf("ver = 1.015, debug = 2\n");
}

void weno7_c_halo_()
{

}

void weno7_c_inner_()
{
    int i, j, k;


    host_flag[MPI_RANK] = param_weno.my_id;
    host_flag[KERNEL_ACTION] = WENO_INNER_FLAG;
    host_flag[GROUP_SIZE] = (param_weno.ny * param_weno.nz) / 64;
    host_flag[REMAIN_POINT] = (param_weno.ny * param_weno.nz) % 64;
    host_flag[IN_PTR] = c_last_u;
    host_flag[OUT_PTR] = c_cur_u + param_weno.lap;
    host_flag[IN_STRIDE] = param_weno.nx + 2 * param_weno.lap;
    host_flag[OUT_STRIDE] = param_weno.nx + 2 * param_weno.lap;
    host_flag[REQUIRE_IN] = param_weno.nx + 2 * param_weno.lap;
    host_flag[REQUIRE_OUT] = param_weno.nx;
    asm volatile ("#nop":::"memory");

    //printf("release host flag %ld, in_ptr = %x, out_ptr = %x, in_s = %d, out_s = %d\n"
    //       , host_flag[0] + 1, c_u, c_u1, param_weno.nx + 2 * param_weno.lap, param_weno.nx);
    host_flag[0] = host_flag[0] + 1 + host_flag[GROUP_SIZE];
    int ig;
    //for(ig = 0; ig < host_flag[GROUP_SIZE] + 1; ig++)
    //{
        //if(param_weno.my_id == 0)
        //    printf("wait ig = %d\n", ig);
        wait_slave_flag();
    //}

    //printf("get last slave flag of this step\n");

    return;
}

