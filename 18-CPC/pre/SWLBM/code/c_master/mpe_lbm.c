#include "../c_header/mpe_lbm.h"

extern int right_nbr,
left_nbr,
up_nbr,
down_nbr,
lu_nbr,
ru_nbr,
ld_nbr,
rd_nbr;


//ori data
int X, Y, Z, Xst, Xed, Yst, Yed, x_sec, y_sec;
int STEPS;
realt *****ori_nodes;
int ***ori_flags, ****ori_walls;
Real ***right_recv;
Real ***left_recv;
Real ***down_recv;
Real ***up_recv;
Real ***right_send;
Real ***left_send;
Real ***down_send;
Real ***up_send;
Real **lu_recv;
Real **ld_recv;
Real **ru_recv;
Real **rd_recv;
Real **lu_send;
Real **ld_send;
Real **ru_send;
Real **rd_send;

//new data

int halo_point_cnt;
int inner_point_cnt;

int *halo_point_xs, *inner_point_xs;
int *halo_point_ys, *inner_point_ys;

//insane
int *insane_halo_point_xs, *insane_inner_point_xs;
int *insane_halo_point_ys, *insane_inner_point_ys;

int insane_halo_point_cnt;
int insane_inner_point_cnt;

int *point_occu_flag;//x_sec * y_sec

//float**
Real* new_nodes;
int* new_flags, *new_walls;

struct lbm_init_param host_param;

#define DBG01 0

// analyse


#define nodes_idx(d1, d2, d3, d4, d5) (((((d1) * (x_sec + 2) + (d2)) * (y_sec + 2) + (d3)) * Z + (d4)) * 19 + (d5))
#define flags_idx(d2, d3, d4) (((d2) * (y_sec + 2) + (d3)) * Z + (d4))
#define walls_idx(d2, d3, d4, d5) ((((d2) * y_sec + (d3)) * Z + (d4)) * 19 + (d5))

int cur_main_iter;

void mix_inner_nowait(
                      int current,
                      int other)
{
    int i, j, k, l;
    int inv;

    Real rho, u_x, u_y, u_z;
    Real fi;
    Real feq[19], Qo, omegaNew;

    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = STD_LBM_FLAG;
    host_flag[STD_POINT_CNT] = inner_point_cnt;
    host_flag[STD_XS_PTR] = (long)inner_point_xs;
    host_flag[STD_YS_PTR] = (long)inner_point_ys;
    host_flag[STD_CURRENT_HEAD] = (long)(new_nodes + nodes_idx(current, 0, 0, 0, 0));
    host_flag[STD_OTHER_HEAD] = (long)(new_nodes + nodes_idx(other, 0, 0, 0, 0));
    host_flag[INSANE_POINT_CNT] = insane_inner_point_cnt;
    host_flag[INSANE_XS_PTR] = (long)insane_inner_point_xs;
    host_flag[INSANE_YS_PTR] = (long)insane_inner_point_ys;

    //调用从核计算
    asm volatile ("nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
}

void mix_halo(
                              int current,
                              int other)
{
    int i, j, k, l;
    int inv;

    Real rho, u_x, u_y, u_z;
    Real fi;
    Real feq[19], Qo, omegaNew;

    host_flag[MPI_RANK] = my_rank;
    host_flag[KERNEL_ACTION] = STD_LBM_FLAG;
    host_flag[STD_POINT_CNT] = halo_point_cnt;
    host_flag[STD_XS_PTR] = (long)halo_point_xs;
    host_flag[STD_YS_PTR] = (long)halo_point_ys;
    host_flag[STD_CURRENT_HEAD] = (long)(new_nodes + nodes_idx(current, 0, 0, 0, 0));
    host_flag[STD_OTHER_HEAD] = (long)(new_nodes + nodes_idx(other, 0, 0, 0, 0));
    host_flag[INSANE_POINT_CNT] = insane_halo_point_cnt;
    host_flag[INSANE_XS_PTR] = (long)insane_halo_point_xs;
    host_flag[INSANE_YS_PTR] = (long)insane_halo_point_ys;

    //调用从核计算
    asm volatile ("nop":::"memory");
    host_flag[0] = host_flag[0] + 1;
    wait_slave_flag();
}

//事实上遍历的范围，xe--，ye--
void insane_halo_point_collect(int xs, int xe, int ys, int ye)
{
    xe--;
    ye--;
    int x, y;

    for(x = xs; x <= xe; x += 2)
    {
        for(y = ys; y <= ye; y += 2)
        {
            if(point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 1)] == 0 &&
               point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 0)] == 0 &&
               point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 1)] == 0 &&
               point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 0)] == 0
             )
            {
                insane_halo_point_xs[insane_halo_point_cnt] = x;
                insane_halo_point_ys[insane_halo_point_cnt] = y;
                insane_halo_point_cnt++;
                point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 1)] = 1;
                point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 0)] = 1;
                point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 1)] = 1;
                point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 0)] = 1;
            }
        }
    }
}

void insane_inner_point_collect(int xs, int xe, int ys, int ye)
{
    xe--;
    ye--;
    int x, y;

    for(x = xs; x <= xe; x += 2)
    {
        for(y = ys; y <= ye; y += 2)
        {
            if(point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 1)] == 0 &&
               point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 0)] == 0 &&
               point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 1)] == 0 &&
               point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 0)] == 0
             )
            {
                insane_inner_point_xs[insane_inner_point_cnt] = x;
                insane_inner_point_ys[insane_inner_point_cnt] = y;
                insane_inner_point_cnt++;
                point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 1)] = 1;
                point_occu_flag[(x - 1 + 1) * y_sec + (y - 1 + 0)] = 1;
                point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 1)] = 1;
                point_occu_flag[(x - 1 + 0) * y_sec + (y - 1 + 0)] = 1;
            }
        }
    }
}

void halo_point_collect(int xs, int xe, int ys, int ye)
{
    int x, y;
    for(x = xs; x <= xe; x++)
    {
        for(y = ys; y <= ye; y++)
        {
            if(point_occu_flag[(x - 1) * y_sec + (y - 1)] == 0)
            {
                halo_point_xs[halo_point_cnt] = x;
                halo_point_ys[halo_point_cnt] = y;
                halo_point_cnt++;
                point_occu_flag[(x - 1) * y_sec + (y - 1)] = 1;
            }
        }
    }
}

void inner_point_collect(int xs, int xe, int ys, int ye)
{
    int x, y;
    for(x = xs; x <= xe; x++)
    {
        for(y = ys; y <= ye; y++)
        {
            if(point_occu_flag[(x - 1) * y_sec + (y - 1)] == 0)
            {
                inner_point_xs[inner_point_cnt] = x;
                inner_point_ys[inner_point_cnt] = y;
                inner_point_cnt++;
                point_occu_flag[(x - 1) * y_sec + (y - 1)] = 1;
            }
        }
    }
}

void lbm_data_init(int _X, int _Y, int _Z, int _Xst, int _Xed, int _Yst, int _Yed, int _x_sec, int _y_sec,
                   realt *****_ori_nodes, int ***_ori_flags, int ****_ori_walls, int _STEPS,
                   Real ***_temp_right, Real ***_temp_left, Real ***_temp_down, Real ***_temp_up,
                   Real ***_temp_right_send, Real ***_temp_left_send, Real ***_temp_down_send, Real ***_temp_up_send,
                   Real **_temp_lu, Real **_temp_ld, Real **_temp_ru, Real **_temp_rd,
                   Real **_temp_lu_send, Real **_temp_ld_send, Real **_temp_ru_send, Real **_temp_rd_send)
{
    //非常重要！！！
    //忘了加这一段，这才是各种玄学出错的根本原因

    int i;
    for(i = 0; i < FLAG_SIZE; i++)
    {
        host_flag[i] = 0;
        slave_flag[i] = 0;
    }
    flag_to_wait = 1;

    long stcc, edcc;
    stcc = rpcc();


    X = _X;
    Y = _Y;
    Z = _Z;
    Xst = _Xst;
    Xed = _Xed;
    Yst = _Yst;
    Yed = _Yed;
    x_sec = _x_sec;
    y_sec = _y_sec;

    ori_nodes = _ori_nodes;
    ori_walls = _ori_walls;
    ori_flags = _ori_flags;

    STEPS = _STEPS;

    right_recv = _temp_right;
    left_recv = _temp_left;
    down_recv = _temp_down;
    up_recv = _temp_up;
    right_send = _temp_right_send;
    left_send = _temp_left_send;
    down_send = _temp_down_send;
    up_send = _temp_up_send;
    lu_recv = _temp_lu;
    ld_recv = _temp_ld;
    ru_recv = _temp_ru;
    rd_recv = _temp_rd;
    lu_send = _temp_lu_send;
    ld_send = _temp_ld_send;
    ru_send = _temp_ru_send;
    rd_send = _temp_rd_send;

    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    new_flags = &ori_flags[0][0][0];
    new_walls = &ori_walls[0][0][0][0];
    new_nodes = &ori_nodes[0][0][0][0][0];

    //超范围创建，为了方便
    inner_point_cnt = halo_point_cnt = 0;
    halo_point_xs = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    halo_point_ys = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    inner_point_xs = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    inner_point_ys = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));

    insane_inner_point_cnt = insane_halo_point_cnt = 0;
    insane_halo_point_xs = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    insane_halo_point_ys = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    insane_inner_point_xs = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));
    insane_inner_point_ys = malloc(sizeof(int) * (x_sec + 2) * (y_sec + 2));

    point_occu_flag = malloc(sizeof(int) * (x_sec) * (y_sec));

    assert(halo_point_xs && halo_point_ys && inner_point_xs && inner_point_ys &&
           insane_halo_point_xs && insane_halo_point_ys && insane_inner_point_xs && insane_inner_point_ys && point_occu_flag);

    memset(point_occu_flag, 0, sizeof(int) * (x_sec) * (y_sec));

    // inner_point_collect(1, x_sec, 1, y_sec);

    insane_inner_point_collect(1 + 2, x_sec - 2, 1 + 2, y_sec - 2);
    inner_point_collect(1 + 2, x_sec - 2, 1 + 2, y_sec - 2);

    insane_halo_point_collect(1, x_sec, 1, 2);
    insane_halo_point_collect(1, 2, 1 + 2, y_sec - 2);
    insane_halo_point_collect(x_sec - 1, x_sec, 1 + 2, y_sec - 2);
    insane_halo_point_collect(1, x_sec, y_sec - 1, y_sec);
    halo_point_collect(1, x_sec, 1, 2);
    halo_point_collect(1, 2, 1 + 2, y_sec - 2);
    halo_point_collect(x_sec - 1, x_sec, 1 + 2, y_sec - 2);
    halo_point_collect(1, x_sec, y_sec - 1, y_sec);

    MPI_Barrier(MPI_COMM_WORLD);

    host_param.my_id = my_rank;
    host_param.host_flag = (long*)&host_flag[0];
    host_param.slave_flag = (long*)&slave_flag[0];
    host_param.iter = STEPS;
    host_param.x_sec = x_sec;
    host_param.y_sec = y_sec;
    host_param.Z = Z;
    host_param.nu = nu;
    host_param.omega = omega;
    host_param.CSmago = CSmago;

    host_param.walls = new_walls;
    host_param.flags = new_flags;

    if(my_rank == 0)
    {
        printf("### MESSAGE FROM MPE ###\n");
        printf("iter = %d, x_sec = %d, y_sec = %d, Z = %d\n", host_param.iter, host_param.x_sec, host_param.y_sec, host_param.Z);
        printf("val: nu = %f, omega = %f, CSmago = %f\n", host_param.nu, host_param.omega, host_param.CSmago);
        printf("ptrs: %x, %x, %x, %x\n",
               host_param.host_flag, host_param.slave_flag, host_param.walls, host_param.flags);
        printf("sizeof(int) = %d, sizeof(Real) = %d\n", sizeof(int), sizeof(Real));
        printf("### MESSAGE FROM CPE ###\n");
    }

    //faker
    int x, y;

    //seq debug
    for(x = 0; x < comm_sz; x++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if(x == my_rank)
        {
            printf("rank = %d, halo_insane = %d, halo_std = %d, inner_insane = %d, inner_std = %d\n",
                    my_rank, insane_halo_point_cnt * 4, halo_point_cnt, insane_inner_point_cnt * 4, inner_point_cnt);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    edcc = rpcc();
    double ttt = ccts(edcc - stcc);
    if(my_rank == 0)
    {
        printf("total time of lbm_data_init is %f\n", ttt);
    }
}

void main_iter(int* _s)
{
    long stcc, edcc;
    stcc = rpcc();


    athread_init();
    athread_spawn(cpe_athread_daemon, (void*)&host_param);

    int current = 0, other = 1;
    int wait_cnt;
    int n = 0;//10 percent
    MPI_Status sta[16];
    MPI_Request req[16];

    int i;
    for(_s[0] = 0; _s[0] < STEPS; _s[0]++)
    {
        cur_main_iter = _s[0];

        mix_inner_nowait(current, other);

        if(my_rank == 0)
            printf("for(s = 0; s < STEPS; s++) s = %d\n", *_s);
        
        if(DBG01 && my_rank == 0)
            assert(1 + 1 == 3);

        //send init
        bounce_send_init(X,
                         Y,
                         Z,
                         Xst,
                         Xed,
                         Yst,
                         Yed,
                         x_sec,
                         y_sec,
                         current,
                         other,
                         ori_nodes,
                         left_send,
                         right_send,
                         up_send,
                         down_send,
                         ld_send,
                         lu_send,
                         rd_send,
                         ru_send);

        //comm
        bounce_communicate(MPI_COMM_WORLD,
                           NULL,
                           NULL,
                           x_sec,
                           y_sec,
                           Z,
                           &wait_cnt,
                           sta,
                           req,
                           left_send,
                           right_send,
                           up_send,
                           down_send,
                           left_recv,
                           right_recv,
                           up_recv,
                           down_recv,
                           lu_send,
                           ld_send,
                           ru_send,
                           rd_send,
                           lu_recv,
                           ld_recv,
                           ru_recv,
                           rd_recv );


        //waitall and update
        for(i = 0; i < wait_cnt; i++) {
            MPI_Wait(&req[i], &sta[i]);
        }

        bounce_update(X,
                      Y,
                      Z,
                      Xst,
                      Xed,
                      Yst,
                      Yed,
                      my_rank,
                      x_sec,
                      y_sec,
                      other,
                      ori_nodes,
                      left_recv,
                      right_recv,
                      up_recv,
                      down_recv,
                      ld_recv,
                      lu_recv,
                      rd_recv,
                      ru_recv);

        wait_slave_flag();

        //stream and collide
        mix_halo(current, other);
        // stream_and_collide_inner(Z, current, other);


        other = current;
        current = (current+1)%2;
    }

    terminate_athread_daemon();

    edcc = rpcc();
    double ttt = ccts(edcc - stcc);
    if(my_rank == 0)
    {
        printf("total time of main_iter is %f\n", ttt);
    }
}
