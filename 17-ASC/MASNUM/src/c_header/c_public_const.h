#ifndef C_PUBLIC_CONST_H
#define C_PUBLIC_CONST_H

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

#define _kl        25
#define _kld       30
#define _klp1      26 // _kl + 1
#define _kldp1     31 // _kldp1
#define _jl        12 
#define _jlp1      13 // _jl + 1
#define _klmjl     300

#define _small     0.000001
#define _alog10pwk 0.08278537031
#define _gc2       7.54515549
#define _cksp      14.0
#define _cksa      4.5
#define _zpi       6.2831852
#define _pi        3.1415926
#define _acu       1.0
#define _g         9.81
#define _abo       1.0
#define _ads       1.0
#define _d1        0.000132
#define _d2        2.61
#define _beta0     1.0
#define _beta10    0.0003125
#define _tztz      1.099314
#define _inv_al11  0.8
#define _inv_al21  1.3333333333333333333333333333333
#define _inv_al31  1.0666666666666666666666666666667
#define _inv_al12  0.64
#define _inv_al22  1.7777777777777777777777777777778
#define _inv_al13  0.512
#define _inv_al23  2.3703703703703703703703703703704

#define CCPS 1450000000

#define comm_torrent 0.5

/*
下面的宏处理 Fortran 多维数组转 C 的一维下标。如果用函数来定义，可能无法进行
内联优化，导致频繁调用函数带来巨大的开销。
*/

// 用于 _ee[], _e[], 维度 (kl, jl, ixs:ixl, iys:iyl)
#define _e_4d_index(k, j, ia, ic) (((ic) - _iys) * dim4 + ((ia) - _ixs) * dim3 + ((j) - 1) * dim2 + ((k) - 1))

// 用于 _se[], _dse[], 维度 (klp1, jl)
#define _se_index(k, j) (((j) - 1) * _klp1 + ((k) - 1))

// 用于 _sein[], _sebo[], _seds[], 维度 (kl, jl)
#define _sexx_index(k, j) (((j) - 1) * _kl + ((k) - 1))

// 用于 _wp[], _wm[], 维度 (kl, 2, 2)
#define _wpwm_index(k, i, j) (((j) - 1) * (2 * _kl) + ((i) - 1) * (_kl) + ((k) - 1))

// 用于 _jp[], _jm[], 维度 (2,jl)
#define _jpjm_index(mr, j) (((j) - 1) * 2 + ((mr) - 1))

// 用于以上没有列出的，并且不是单独写转换表达式的数组，维度 (ixs:ixl, iys:iyl)
#define _iaic_index(ia, ic) (((ic) - _iys) * ix_size + ((ia) - _ixs))

#define _e_index_iaic(ia, ic) (((ic) - _iys) * dim4 + ((ia) - _ixs) * dim3)

#define _e_index_kj(k, j) (((j) - 1) * dim2 + ((k) - 1))

#endif
