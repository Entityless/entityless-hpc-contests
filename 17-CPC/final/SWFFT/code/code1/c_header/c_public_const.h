#ifndef C_PUBLIC_CONST_H
#define C_PUBLIC_CONST_H

#define CCPS 1450000000

struct dcomplex
{
    double real, imag;
};

struct fcomplex
{
    float real, imag;
};

#define dfc struct dcomplex
#define dft double
#define sfc struct fcomplex
#define sft float

#define my_fc struct dcomplex
#define my_ft double
#define my_mpift MPI_DOUBLE_COMPLEX

//please make sure that 在这一层计时，没有互相嵌套的
//防止RPCC调用开销带来的系统误差
#define PF_INITIAL_TRANSFER 1
#define PF_EVV              2
#define PF_CFFTS_COPY       3
#define PF_CFFTZ_COPY       4
#define PF_FFTZ2            5
#define PF_TP1              6
#define PF_TP2              7
#define PF_CHK              8

//#define my_fc struct fcomplex
//#define my_ft float
//#define my_mpift MPI_DOUBLE


#define dctofc(a, b) {b.real = a.real; b.imag = a.imag;}

//code: Ga1ahad and Scientific Witchery

static inline int int_pow(int a, int b)
{
    int i;
    int ret = 1;
    for(i = 0; i < b; i++)
        ret *= a;
    return ret;
}

static inline my_fc dconjg(my_fc v)
{
    my_fc nv;
    nv.real = v.real;
    nv.imag = -v.imag;

    return nv;
}

static inline my_fc add_dc_dc(my_fc v1, my_fc v2)
{
    my_fc nv;
    nv.imag = v1.imag + v2.imag;
    nv.real = v1.real + v2.real;
    return nv;
}

static inline my_fc minus_dc_dc(my_fc v1, my_fc v2)
{
    my_fc nv;
    nv.imag = v1.imag - v2.imag;
    nv.real = v1.real - v2.real;
    return nv;
}

static inline my_fc multiply_dc_dc(my_fc v1, my_fc v2)
{
    my_fc nv;
    nv.imag = v1.imag * v2.real + v2.imag * v1.real;
    nv.real = v1.real * v2.real - v1.imag * v2.imag;
    return nv;
}

static inline my_fc multiply_dc_d(my_fc v1, my_ft v2)
{
    my_fc nv;
    nv.imag = v1.imag * v2;
    nv.real = v1.real * v2;
    return nv;
}

int ilog2(int n)
{
    int nn, lg;
    if(n == 1)
        return 0;
    lg = 1;
    nn = 2;
    while(nn < n)
    {
        nn *= 2;
        lg += 1;
    }
    return lg;
}

#endif
