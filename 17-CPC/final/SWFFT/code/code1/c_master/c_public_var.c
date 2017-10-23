#include "../c_header/c_public_var.h"
#include "../c_header/c_public_const.h"


volatile long host_flag[FLAG_SIZE];
volatile long slave_flag[FLAG_SIZE];
volatile long local_cc[FLAG_SIZE];
volatile long flag_to_wait;

////FFF
//ptr
struct dcomplex* f_u0_, *f_u1_, *f_u2_;
int* dims_;
double* f_twiddle_;
struct dcomplex* f_u_;
int* xstart_, *ystart_, *zstart_;
int* xend_, *yend_, *zend_;
struct dcomplex* f_sums_;
//val
int niter_;
int fftblockpad_default_, maxdim_, layout_type_;
int fftblock_, fftblockpad_;
int me_, np1_, np2_, np_, ntdivnp_;
int transblock_, transblockpad_;
int nx_, ny_, nz_;
double ntotal_f_;
//local pre-allocated ptr
my_fc* scratch_;
//comm
MPI_Comm commslice1_, commslice2_;
//c_local
//my_fc* c_u0_, *c_u1_, *c_u2_;
//my_ft* c_twiddle_;
//dfc* c_u_;
dfc* c_sums_;
double* c_u0_real_, *c_u0_imag_;
double* c_u1_real_, *c_u1_imag_;
double* c_u2_real_, *c_u2_imag_;
double* c_u_real_, *c_u_imag_;
double* c_twiddle_;
//support
int chk_cnt_;
int* chk_xs_;
int* chk_ys_;
int* chk_zs_;
dfc* buf_chk_;
int* chk_core_cnt_;
int* chk_core_stp_;


////mpe benchmark
unsigned long mpe_cc_cur_[100];
unsigned long mpe_cc_total_[100];


void c_param_init_(struct dcomplex* _u0, struct dcomplex* _u1, struct dcomplex* _u2, double* _twiddle, int* _dims, int* _niter,
                   int* _fftblockpad_default, int* _fftblock, int* _maxdim, int* _layout_type,
                   struct dcomplex* _u, int* _me, int* _np1, int* _np2, int* _ntdivnp,
                   int* _np, int* _transblockpad, int* _transblock, int* _fftblockpad,
                   int* _nx, int* _ny, int* _nz,
                   int* _xstart,int *_ystart,int *_zstart,
                   int* _xend,int *_yend,int *_zend, struct dcomplex* _sums)
{
    //ptr copy
       f_u0_ = _u0;
       f_u1_ = _u1;
       f_u2_ = _u2;
  f_twiddle_ = _twiddle;
       dims_ = _dims;
        f_u_ = _u;
     xstart_ = _xstart;
     ystart_ = _ystart;
     zstart_ = _zstart;
       xend_ = _xend;
       yend_ = _yend;
       zend_ = _zend;
     f_sums_ = _sums;

    //val copy
                  niter_ = *_niter;
    fftblockpad_default_ = *_fftblockpad_default;
            fftblockpad_ = *_fftblockpad;
                 maxdim_ = *_maxdim;
               fftblock_ = *_fftblock;
            layout_type_ = *_layout_type;
                     me_ = *_me;
                    np1_ = *_np1;
                    np2_ = *_np2;
                     np_ = *_np;
                ntdivnp_ = *_ntdivnp;
          transblockpad_ = *_transblockpad;
             transblock_ = *_transblock;
                     nx_ = *_nx;
                     ny_ = *_ny;
                     nz_ = *_nz;


    //allocate
    scratch_ = (my_fc*)malloc(sizeof(my_fc) * fftblockpad_default_ * maxdim_ * 2);

    //compute
    ntotal_f_ = 1.0 * nx_ * ny_ * nz_;

    // if(me_ == 0)
    // {
    //     printf("<check> c\n");
    //     printf("%d %d %d\n %d %d %d\n %d %d %d\n %d\n %d %d %d\n%d %d %d\n",
    //            niter_, fftblockpad_default_, maxdim_,
    //            fftblock_, layout_type_, me_,
    //            np1_, np2_, np_,
    //            ntdivnp_,
    //            transblockpad_, transblock_, fftblockpad_,
    //            nx_, ny_, nz_);
    //     printf("</check> c\n");
    // }

    //array initial
    int i;
    for(i = 0; i < 100; i++)
    {
        mpe_cc_cur_[i] = mpe_cc_total_[i] = 0;
    }

    for(i = 0; i < FLAG_SIZE; i++)
    {
        host_flag[i] = 0;
        slave_flag[i] = 0;
        local_cc[i] = 0;
    }

    flag_to_wait = 1;

    //communication initial
    //me2 = mod(me, np2)  ! goes from 0...np2-1
    //me1 = me/np2        ! goes from 0...np1-1
    //call MPI_Comm_split(MPI_COMM_WORLD, me1, me2, commslice1, ierr)
    //call MPI_Comm_split(MPI_COMM_WORLD, me2, me1, commslice2, ierr)
    int me1, me2;
    me2 = me_ % np2_;//col no
    me1 = me_ / np2_;//row no
    MPI_Comm_split(MPI_COMM_WORLD, me1, me2, &commslice1_);
    MPI_Comm_split(MPI_COMM_WORLD, me2, me1, &commslice2_);

    c_u0_real_ = (double*)f_u1_;
    c_u0_imag_ = (double*)f_u1_;
    c_u1_real_ = (double*)f_u0_;
    c_u1_imag_ = (double*)f_u0_;
    c_u2_real_ = (double*)f_u2_;
    c_u2_imag_ = (double*)f_u2_;
    c_u_real_ = (double*)malloc(sizeof(double) * nx_);
    c_u_imag_ = (double*)malloc(sizeof(double) * nx_);

    c_u0_imag_ = c_u0_imag_ + ntdivnp_;
    c_u1_imag_ = c_u1_imag_ + ntdivnp_;
    c_u2_imag_ = c_u2_imag_ + ntdivnp_;

    //copy value from f to c
    for(i = 0; i < ntdivnp_; i++)
    {
        c_u0_real_[i] = f_u0_[i].real;
        c_u0_imag_[i] = f_u0_[i].imag;
    }
    for(i = 0; i < nx_; i++)
    {
        c_u_real_[i] = f_u_[i].real;
        c_u_imag_[i] = f_u_[i].imag;
    }

    c_twiddle_ = f_twiddle_;

    c_sums_ = f_sums_;
}

//用这个函数来测试复数参数传递
void fortran_c_dbg_(struct dcomplex* dcp, double* dp)
{
    if(me_ == 0)
    {
        /*
        printf("<c dbg\n");
        printf("c: intital d = %.8f, dc: real = %.8f, imag = %.8f\n", *dp, dcp->real, dcp->imag);

        struct dcomplex v, ret1, ret2, ret0, ret3, ret4;

        v.real = -2.0;
        v.imag = -2.0;

        ret0 = dconjg(*dcp);
        ret1 = multiply_dc_d(*dcp, *dp);
        ret2 = multiply_dc_dc(*dcp, v);
        ret3 = add_dc_dc(*dcp, v);
        ret4 = minus_dc_dc(*dcp, v);
        printf("%.8f, %.8f\n", ret0.real, ret0.imag);
        printf("%.8f, %.8f\n", ret1.real, ret1.imag);
        printf("%.8f, %.8f\n", ret2.real, ret2.imag);
        printf("%.8f, %.8f\n", ret3.real, ret3.imag);
        printf("%.8f, %.8f\n", ret4.real, ret4.imag);
        printf("</c dbg\n");*/
    }
}


void c_finalize_()
{
}

void wait_slave_flag()
{
    //if(my_rank == 0)
    //    printf("wait flag %ld\n", flag_to_wait);
    while(slave_flag[0] < flag_to_wait);
    //if(my_rank == 0)
    //    printf("get flag %ld\n", flag_to_wait);
    flag_to_wait++;
}
