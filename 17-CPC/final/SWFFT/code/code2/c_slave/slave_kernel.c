#include "../c_header/slave_kernel.h"
#include "../c_header/c_public_const.h"

#include "slave.h"
#include "dma.h"
#include "simd.h"
#include <stdio.h>
#include <memory.h>

int BLOCK_SIZE(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) + ((n % total_blocks > block_id) ? 1 : 0);
}

int BLOCK_LOW(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) * block_id + ((n % total_blocks > block_id) ? block_id : n % total_blocks);
}

#define SIMD_COMM_SIZE 8

#define simd_4x4_transpose(r1, r2, r3, r4)        \
{                                                 \
    doublev4 temp1 = simd_vshff(r1, r2, 0x22);    \
    doublev4 temp2 = simd_vshff(r1, r2, 0x77);    \
    doublev4 temp3 = simd_vshff(r3, r4, 0x22);    \
    doublev4 temp4 = simd_vshff(r3, r4, 0x77);    \
    r3 = simd_vshff(temp3, temp1, 0x22);          \
    r4 = simd_vshff(temp4, temp2, 0x22);          \
    r1 = simd_vshff(temp3, temp1, 0x77);          \
    r2 = simd_vshff(temp4, temp2, 0x77);          \
}

#define simd_complex_mul(a_real, a_img, b_real, b_img, real, img) \
real = a_img * b_img;                                             \
real = a_real * b_real - real;                                    \
img = a_real * b_img;                                             \
img = a_img * b_real + img;                                       \

#define fft_kernel(a_real, a_img, b_real, b_img, u_real, u_img)     \
{                                                                   \
    doublev4 res_real = a_real + b_real;                            \
    doublev4 res_img = a_img + b_img;                               \
    a_real = a_real - b_real;                                       \
    a_img = a_img - b_img;                                          \
    simd_complex_mul(a_real, a_img, u_real, u_img, b_real, b_img);  \
    a_real = res_real;                                              \
    a_img = res_img;                                                \
}



#define simd_hoshff(a, b)                   \
{                                           \
    doublev4 temp = simd_vshff(b, a, 0x44); \
    b = simd_vshff(b, a, 0xEE);             \
    a = temp;                               \
}


#define simd_butshff(a, b)                  \
{                                           \
    simd_hoshff(a,b);                       \
    a = simd_vshff(a, a, 0xD8);             \
    b = simd_vshff(b, b, 0xD8);             \
}

//12345678 -> 15263748
#define fft512_phase1(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real, \
                      a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag, \
                      u1_real, u2_real, u3_real, u4_real,                                     \
                      u1_imag, u2_imag, u3_imag, u4_imag)                                     \
{                                                                                             \
    fft_kernel(a1_real, a1_imag, a5_real, a5_imag, u1_real, u1_imag);                         \
    simd_butshff(a1_real, a5_real);                                                           \
    simd_butshff(a1_imag, a5_imag);                                                           \
    fft_kernel(a2_real, a2_imag, a6_real, a6_imag, u2_real, u2_imag);                         \
    simd_butshff(a2_real, a6_real);                                                           \
    simd_butshff(a2_imag, a6_imag);                                                           \
    fft_kernel(a3_real, a3_imag, a7_real, a7_imag, u3_real, u3_imag);                         \
    simd_butshff(a3_real, a7_real);                                                           \
    simd_butshff(a3_imag, a7_imag);                                                           \
    fft_kernel(a4_real, a4_imag, a8_real, a8_imag, u4_real, u4_imag);                         \
    simd_butshff(a4_real, a8_real);                                                           \
    simd_butshff(a4_imag, a8_imag);                                                           \
}


//12345678 -> 15263748
#define fft512_phase2(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real, \
                     a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag, \
                     u1_real, u2_real, u3_real, u4_real,                                     \
                     u1_imag, u2_imag, u3_imag, u4_imag)                                     \
{                                                                                            \
    fft_kernel(a1_real, a1_imag, a5_real, a5_imag, u1_real, u1_imag);                        \
    simd_hoshff(a1_real, a5_real);                                                           \
    simd_hoshff(a1_imag, a5_imag);                                                           \
    fft_kernel(a2_real, a2_imag, a6_real, a6_imag, u2_real, u2_imag);                        \
    simd_hoshff(a2_real, a6_real);                                                           \
    simd_hoshff(a2_imag, a6_imag);                                                           \
    fft_kernel(a3_real, a3_imag, a7_real, a7_imag, u3_real, u3_imag);                        \
    simd_hoshff(a3_real, a7_real);                                                           \
    simd_hoshff(a3_imag, a7_imag);                                                           \
    fft_kernel(a4_real, a4_imag, a8_real, a8_imag, u4_real, u4_imag);                        \
    simd_hoshff(a4_real, a8_real);                                                           \
    simd_hoshff(a4_imag, a8_imag);                                                           \
}

//12345678 -> 15263748
//12345678 -> 12563478
//12345678 -> 12345678
#define fft512_phase3(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real,    \
                     a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag,    \
                     u1_real, u2_real, u3_real, u4_real,                                        \
                     u1_imag, u2_imag, u3_imag, u4_imag)                                        \
{                                                                                               \
    fft_kernel(a1_real, a1_imag, a5_real, a5_imag, u1_real, u1_imag);                           \
    fft_kernel(a2_real, a2_imag, a6_real, a6_imag, u2_real, u2_imag);                           \
    fft_kernel(a3_real, a3_imag, a7_real, a7_imag, u3_real, u3_imag);                           \
    fft_kernel(a4_real, a4_imag, a8_real, a8_imag, u4_real, u4_imag);                           \
}

#define fft512_phase4(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real,    \
                     a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag,    \
                     u1_real, u2_real,                                                          \
                     u1_imag, u2_imag)                                                          \
{                                                                                               \
    fft_kernel(a1_real, a1_imag, a5_real, a5_imag, u1_real, u1_imag);                           \
    fft_kernel(a2_real, a2_imag, a6_real, a6_imag, u1_real, u1_imag);                           \
    fft_kernel(a3_real, a3_imag, a7_real, a7_imag, u2_real, u2_imag);                           \
    fft_kernel(a4_real, a4_imag, a8_real, a8_imag, u2_real, u2_imag);                           \
}

#define fft512_phase5(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real,    \
                     a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag,    \
                     u1_real,                                                                   \
                     u1_imag)                                                                   \
{                                                                                               \
    fft_kernel(a1_real, a1_imag, a5_real, a5_imag, u1_real, u1_imag);                             \
    fft_kernel(a2_real, a2_imag, a6_real, a6_imag, u1_real, u1_imag);                             \
    fft_kernel(a3_real, a3_imag, a7_real, a7_imag, u1_real, u1_imag);                             \
    fft_kernel(a4_real, a4_imag, a8_real, a8_imag, u1_real, u1_imag);                             \
}

void simd_fft_phase1(double* real, double* imag, double* ureal, double* uimag, 
                     double* dest_real, double* dest_imag, 
                     int degree, int len)
{
    int load_stride = len >> 1;

    doublev4 a1_real; simd_load(a1_real, real + 0 * 4              );
    doublev4 a2_real; simd_load(a2_real, real + 1 * 4              );
    doublev4 a3_real; simd_load(a3_real, real + 2 * 4              );
    doublev4 a4_real; simd_load(a4_real, real + 3 * 4              );
    doublev4 a5_real; simd_load(a5_real, real + 0 * 4 + load_stride);
    doublev4 a6_real; simd_load(a6_real, real + 1 * 4 + load_stride);
    doublev4 a7_real; simd_load(a7_real, real + 2 * 4 + load_stride);
    doublev4 a8_real; simd_load(a8_real, real + 3 * 4 + load_stride);

    doublev4 a1_imag; simd_load(a1_imag, imag + 0 * 4              );
    doublev4 a2_imag; simd_load(a2_imag, imag + 1 * 4              );
    doublev4 a3_imag; simd_load(a3_imag, imag + 2 * 4              );
    doublev4 a4_imag; simd_load(a4_imag, imag + 3 * 4              );
    doublev4 a5_imag; simd_load(a5_imag, imag + 0 * 4 + load_stride);
    doublev4 a6_imag; simd_load(a6_imag, imag + 1 * 4 + load_stride);
    doublev4 a7_imag; simd_load(a7_imag, imag + 2 * 4 + load_stride);
    doublev4 a8_imag; simd_load(a8_imag, imag + 3 * 4 + load_stride);

    double* ureal_ptr = ureal + load_stride;
    double* uimag_ptr = uimag + load_stride;

    doublev4 u1_real; simd_load(u1_real, ureal_ptr +  0);
    doublev4 u2_real; simd_load(u2_real, ureal_ptr +  4);
    doublev4 u3_real; simd_load(u3_real, ureal_ptr +  8);
    doublev4 u4_real; simd_load(u4_real, ureal_ptr + 12);
    doublev4 u1_imag; simd_load(u1_imag, uimag_ptr +  0);
    doublev4 u2_imag; simd_load(u2_imag, uimag_ptr +  4);
    doublev4 u3_imag; simd_load(u3_imag, uimag_ptr +  8);
    doublev4 u4_imag; simd_load(u4_imag, uimag_ptr + 12);

    fft512_phase1(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real,
                  a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag,
                  u1_real, u2_real, u3_real, u4_real,
                  u1_imag, u2_imag, u3_imag, u4_imag);

    simd_store(a1_real, dest_real + 0 * 4);
    simd_store(a5_real, dest_real + 1 * 4);
    simd_store(a2_real, dest_real + 2 * 4);
    simd_store(a6_real, dest_real + 3 * 4);
    simd_store(a3_real, dest_real + 4 * 4);
    simd_store(a7_real, dest_real + 5 * 4);
    simd_store(a4_real, dest_real + 6 * 4);
    simd_store(a8_real, dest_real + 7 * 4);

    simd_store(a1_imag, dest_imag + 0 * 4);
    simd_store(a5_imag, dest_imag + 1 * 4);
    simd_store(a2_imag, dest_imag + 2 * 4);
    simd_store(a6_imag, dest_imag + 3 * 4);
    simd_store(a3_imag, dest_imag + 4 * 4);
    simd_store(a7_imag, dest_imag + 5 * 4);
    simd_store(a4_imag, dest_imag + 6 * 4);
    simd_store(a8_imag, dest_imag + 7 * 4);

}

void simd_fft_phase2(double* real, double* imag, double* ureal, double* uimag, 
                     double* dest_real, double* dest_imag, 
                     int degree, int len)
{
    int load_stride = len >> 1;
    doublev4 a1_real; simd_load(a1_real, real + 0 * 4              );
    doublev4 a2_real; simd_load(a2_real, real + 1 * 4              );
    doublev4 a3_real; simd_load(a3_real, real + 2 * 4              );
    doublev4 a4_real; simd_load(a4_real, real + 3 * 4              );
    doublev4 a5_real; simd_load(a5_real, real + 0 * 4 + load_stride);
    doublev4 a6_real; simd_load(a6_real, real + 1 * 4 + load_stride);
    doublev4 a7_real; simd_load(a7_real, real + 2 * 4 + load_stride);
    doublev4 a8_real; simd_load(a8_real, real + 3 * 4 + load_stride);
    
    doublev4 a1_imag; simd_load(a1_imag, imag + 0 * 4              );
    doublev4 a2_imag; simd_load(a2_imag, imag + 1 * 4              );
    doublev4 a3_imag; simd_load(a3_imag, imag + 2 * 4              );
    doublev4 a4_imag; simd_load(a4_imag, imag + 3 * 4              );
    doublev4 a5_imag; simd_load(a5_imag, imag + 0 * 4 + load_stride);
    doublev4 a6_imag; simd_load(a6_imag, imag + 1 * 4 + load_stride);
    doublev4 a7_imag; simd_load(a7_imag, imag + 2 * 4 + load_stride);
    doublev4 a8_imag; simd_load(a8_imag, imag + 3 * 4 + load_stride);

    double* ureal_ptr = ureal + (load_stride >> 1);
    double* uimag_ptr = uimag + (load_stride >> 1);

    doublev4 u1_real = simd_set_doublev4(ureal_ptr[0], ureal_ptr[0], ureal_ptr[1], ureal_ptr[1]);
    doublev4 u2_real = simd_set_doublev4(ureal_ptr[2], ureal_ptr[2], ureal_ptr[3], ureal_ptr[3]);
    doublev4 u3_real = simd_set_doublev4(ureal_ptr[4], ureal_ptr[4], ureal_ptr[5], ureal_ptr[5]);
    doublev4 u4_real = simd_set_doublev4(ureal_ptr[6], ureal_ptr[6], ureal_ptr[7], ureal_ptr[7]);

    doublev4 u1_imag = simd_set_doublev4(uimag_ptr[0], uimag_ptr[0], uimag_ptr[1], uimag_ptr[1]);
    doublev4 u2_imag = simd_set_doublev4(uimag_ptr[2], uimag_ptr[2], uimag_ptr[3], uimag_ptr[3]);
    doublev4 u3_imag = simd_set_doublev4(uimag_ptr[4], uimag_ptr[4], uimag_ptr[5], uimag_ptr[5]);
    doublev4 u4_imag = simd_set_doublev4(uimag_ptr[6], uimag_ptr[6], uimag_ptr[7], uimag_ptr[7]);

    fft512_phase2(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real,
                  a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag,
                  u1_real, u2_real, u3_real, u4_real,
                  u1_imag, u2_imag, u3_imag, u4_imag);

    simd_store(a1_real, dest_real + 0 * 4);
    simd_store(a5_real, dest_real + 1 * 4);
    simd_store(a2_real, dest_real + 2 * 4);
    simd_store(a6_real, dest_real + 3 * 4);
    simd_store(a3_real, dest_real + 4 * 4);
    simd_store(a7_real, dest_real + 5 * 4);
    simd_store(a4_real, dest_real + 6 * 4);
    simd_store(a8_real, dest_real + 7 * 4);

    simd_store(a1_imag, dest_imag + 0 * 4);
    simd_store(a5_imag, dest_imag + 1 * 4);
    simd_store(a2_imag, dest_imag + 2 * 4);
    simd_store(a6_imag, dest_imag + 3 * 4);
    simd_store(a3_imag, dest_imag + 4 * 4);
    simd_store(a7_imag, dest_imag + 5 * 4);
    simd_store(a4_imag, dest_imag + 6 * 4);
    simd_store(a8_imag, dest_imag + 7 * 4);

}


void simd_fft_general_3step(double* real, double* imag, double* ureal, double* uimag, 
                            double* dest_real, double* dest_imag, 
                            int load_stride, int store_stride,
                            int u_offset, int u_len)
{
    doublev4 a1_real; simd_load(a1_real, real + 0 * load_stride);
    doublev4 a2_real; simd_load(a2_real, real + 1 * load_stride);
    doublev4 a3_real; simd_load(a3_real, real + 2 * load_stride);
    doublev4 a4_real; simd_load(a4_real, real + 3 * load_stride);
    doublev4 a5_real; simd_load(a5_real, real + 4 * load_stride);
    doublev4 a6_real; simd_load(a6_real, real + 5 * load_stride);
    doublev4 a7_real; simd_load(a7_real, real + 6 * load_stride);
    doublev4 a8_real; simd_load(a8_real, real + 7 * load_stride);
    
    doublev4 a1_imag; simd_load(a1_imag, imag + 0 * load_stride);
    doublev4 a2_imag; simd_load(a2_imag, imag + 1 * load_stride);
    doublev4 a3_imag; simd_load(a3_imag, imag + 2 * load_stride);
    doublev4 a4_imag; simd_load(a4_imag, imag + 3 * load_stride);
    doublev4 a5_imag; simd_load(a5_imag, imag + 4 * load_stride);
    doublev4 a6_imag; simd_load(a6_imag, imag + 5 * load_stride);
    doublev4 a7_imag; simd_load(a7_imag, imag + 6 * load_stride);
    doublev4 a8_imag; simd_load(a8_imag, imag + 7 * load_stride);

    double* ureal_ptr = ureal + u_offset;
    double* uimag_ptr = uimag + u_offset;
    int u_stride = u_len >> 2;

    doublev4 u1_real = ureal_ptr[           0];
    doublev4 u2_real = ureal_ptr[    u_stride];
    doublev4 u3_real = ureal_ptr[2 * u_stride];
    doublev4 u4_real = ureal_ptr[3 * u_stride];

    doublev4 u1_imag = uimag_ptr[           0];
    doublev4 u2_imag = uimag_ptr[    u_stride];
    doublev4 u3_imag = uimag_ptr[2 * u_stride];
    doublev4 u4_imag = uimag_ptr[3 * u_stride];

    fft512_phase3(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real,
                  a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag,
                  u1_real, u2_real, u3_real, u4_real,
                  u1_imag, u2_imag, u3_imag, u4_imag);

    u_offset -= (u_len >> 1);
    ureal_ptr = ureal + u_offset;
    uimag_ptr = uimag + u_offset;

    u1_real = ureal_ptr[       0];
    u2_real = ureal_ptr[u_stride];
    u1_imag = uimag_ptr[       0];
    u2_imag = uimag_ptr[u_stride];

    fft512_phase4(a1_real, a5_real, a2_real, a6_real, a3_real, a7_real, a4_real, a8_real,
                  a1_imag, a5_imag, a2_imag, a6_imag, a3_imag, a7_imag, a4_imag, a8_imag,
                  u1_real, u2_real,
                  u1_imag, u2_imag);

    u_offset -= (u_len >> 2);
    ureal_ptr = ureal + u_offset;
    uimag_ptr = uimag + u_offset;
    u1_real = ureal_ptr[0];
    u1_imag = uimag_ptr[0];

    fft512_phase5(a1_real, a5_real, a3_real, a7_real, a2_real, a6_real, a4_real, a8_real,
                  a1_imag, a5_imag, a3_imag, a7_imag, a2_imag, a6_imag, a4_imag, a8_imag,
                  u1_real,
                  u1_imag);

    simd_store(a1_real, dest_real + 0 * store_stride);
    simd_store(a5_real, dest_real + 1 * store_stride);
    simd_store(a3_real, dest_real + 2 * store_stride);
    simd_store(a7_real, dest_real + 3 * store_stride);
    simd_store(a2_real, dest_real + 4 * store_stride);
    simd_store(a6_real, dest_real + 5 * store_stride);
    simd_store(a4_real, dest_real + 6 * store_stride);
    simd_store(a8_real, dest_real + 7 * store_stride);

    simd_store(a1_imag, dest_imag + 0 * store_stride);
    simd_store(a5_imag, dest_imag + 1 * store_stride);
    simd_store(a3_imag, dest_imag + 2 * store_stride);
    simd_store(a7_imag, dest_imag + 3 * store_stride);
    simd_store(a2_imag, dest_imag + 4 * store_stride);
    simd_store(a6_imag, dest_imag + 5 * store_stride);
    simd_store(a4_imag, dest_imag + 6 * store_stride);
    simd_store(a8_imag, dest_imag + 7 * store_stride);

}


void simd_fft512_32_phase1to5(double* real, double* imag, double* ureal, double* uimag,
                              double* dest_real, double* dest_imag)
{
    int i;
    for(i = 0; i < 16; ++i)
    {
        simd_fft_phase1(real + i * 16, imag + i * 16, ureal + i * 16, uimag + i * 16, 
                     dest_real + i * 32, dest_imag + i * 32, 
                     9, 512);
    }

    for(i = 0; i < 16; ++i)
    {
        simd_fft_phase2(dest_real + i * 16, dest_imag + i * 16, ureal + i * 8, uimag + i * 8, 
                        real + i * 32,  imag + i * 32, 
                        9, 512);
    }

    for(i = 0; i < 16; ++i)
    {
        simd_fft_general_3step(real + i * 4, imag + i * 4, ureal, uimag,
                               dest_real + i * 32, dest_imag + i * 32,
                               64, 4, 
                               64 + i, 64);
    }


}

// void simd_fft512_32_phase1to5(double* real, double* imag, double* ureal, double* uimag,
//                               double* dest_real, double* dest_imag)
// {
//     double* ureal_ptr;
//     double* uimag_ptr;

//     doublev4 a1_real = simd_set_doublev4(real[0 * 16],  real[1 * 16],  real[2 * 16],  real[3 * 16]);
//     doublev4 a2_real = simd_set_doublev4(real[4 * 16],  real[5 * 16],  real[6 * 16],  real[7 * 16]);
//     doublev4 a3_real = simd_set_doublev4(real[8 * 16],  real[9 * 16],  real[10 * 16], real[11 * 16]);
//     doublev4 a4_real = simd_set_doublev4(real[12 * 16], real[13 * 16], real[14 * 16], real[15 * 16]);
//     doublev4 a5_real = simd_set_doublev4(real[16 * 16], real[17 * 16], real[18 * 16], real[19 * 16]);
//     doublev4 a6_real = simd_set_doublev4(real[20 * 16], real[21 * 16], real[22 * 16], real[23 * 16]);
//     doublev4 a7_real = simd_set_doublev4(real[24 * 16], real[25 * 16], real[26 * 16], real[27 * 16]);
//     doublev4 a8_real = simd_set_doublev4(real[28 * 16], real[29 * 16], real[30 * 16], real[31 * 16]);
    
//     doublev4 a1_imag = simd_set_doublev4(imag[0 * 16],  imag[1 * 16],  imag[2 * 16],  imag[3 * 16]);
//     doublev4 a2_imag = simd_set_doublev4(imag[4 * 16],  imag[5 * 16],  imag[6 * 16],  imag[7 * 16]);
//     doublev4 a3_imag = simd_set_doublev4(imag[8 * 16],  imag[9 * 16],  imag[10 * 16], imag[11 * 16]);
//     doublev4 a4_imag = simd_set_doublev4(imag[12 * 16], imag[13 * 16], imag[14 * 16], imag[15 * 16]);
//     doublev4 a5_imag = simd_set_doublev4(imag[16 * 16], imag[17 * 16], imag[18 * 16], imag[19 * 16]);
//     doublev4 a6_imag = simd_set_doublev4(imag[20 * 16], imag[21 * 16], imag[22 * 16], imag[23 * 16]);
//     doublev4 a7_imag = simd_set_doublev4(imag[24 * 16], imag[25 * 16], imag[26 * 16], imag[27 * 16]);
//     doublev4 a8_imag = simd_set_doublev4(imag[28 * 16], imag[29 * 16], imag[30 * 16], imag[31 * 16]);

//     ureal_ptr = ureal + 256;
//     uimag_ptr = uimag + 256;
//     doublev4 u1_real = simd_set_doublev4(ureal_ptr[0 * 16],  ureal_ptr[1 * 16],  ureal_ptr[2 * 16],  ureal_ptr[3 * 16]);
//     doublev4 u2_real = simd_set_doublev4(ureal_ptr[4 * 16],  ureal_ptr[5 * 16],  ureal_ptr[6 * 16],  ureal_ptr[7 * 16]);
//     doublev4 u3_real = simd_set_doublev4(ureal_ptr[8 * 16],  ureal_ptr[9 * 16],  ureal_ptr[10 * 16], ureal_ptr[11 * 16]);
//     doublev4 u4_real = simd_set_doublev4(ureal_ptr[12 * 16], ureal_ptr[13 * 16], ureal_ptr[14 * 16], ureal_ptr[15 * 16]);

//     doublev4 u1_imag = simd_set_doublev4(uimag_ptr[0 * 16],  uimag_ptr[1 * 16],  uimag_ptr[2 * 16],  uimag_ptr[3 * 16]);
//     doublev4 u2_imag = simd_set_doublev4(uimag_ptr[4 * 16],  uimag_ptr[5 * 16],  uimag_ptr[6 * 16],  uimag_ptr[7 * 16]);
//     doublev4 u3_imag = simd_set_doublev4(uimag_ptr[8 * 16],  uimag_ptr[9 * 16],  uimag_ptr[10 * 16], uimag_ptr[11 * 16]);
//     doublev4 u4_imag = simd_set_doublev4(uimag_ptr[12 * 16], uimag_ptr[13 * 16], uimag_ptr[14 * 16], uimag_ptr[15 * 16]);

//     fft512_phase1(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real,
//                   a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag,
//                   u1_real, u2_real, u3_real, u4_real,
//                   u1_imag, u2_imag, u3_imag, u4_imag);

//     ureal_ptr = ureal + 128;
//     uimag_ptr = uimag + 128;

//     u1_real = simd_set_doublev4(ureal_ptr[0* 16],  ureal_ptr[0* 16],  ureal_ptr[1* 16],  ureal_ptr[1* 16]);
//     u2_real = simd_set_doublev4(ureal_ptr[2* 16],  ureal_ptr[2* 16],  ureal_ptr[3* 16],  ureal_ptr[3* 16]);
//     u3_real = simd_set_doublev4(ureal_ptr[4* 16],  ureal_ptr[4* 16],  ureal_ptr[5 * 16], ureal_ptr[5 * 16]);
//     u4_real = simd_set_doublev4(ureal_ptr[6 * 16], ureal_ptr[6 * 16], ureal_ptr[7 * 16], ureal_ptr[7 * 16]);

//     u1_imag = simd_set_doublev4(uimag_ptr[0* 16],  uimag_ptr[0* 16],  uimag_ptr[1* 16],  uimag_ptr[1* 16]);
//     u2_imag = simd_set_doublev4(uimag_ptr[2* 16],  uimag_ptr[2* 16],  uimag_ptr[3* 16],  uimag_ptr[3* 16]);
//     u3_imag = simd_set_doublev4(uimag_ptr[4* 16],  uimag_ptr[4* 16],  uimag_ptr[5 * 16], uimag_ptr[5 * 16]);
//     u4_imag = simd_set_doublev4(uimag_ptr[6 * 16], uimag_ptr[6 * 16], uimag_ptr[7 * 16], uimag_ptr[7 * 16]);

//     fft512_phase2(a1_real, a5_real, a2_real, a6_real, a3_real, a7_real, a4_real, a8_real, 
//                   a1_imag, a5_imag, a2_imag, a6_imag, a3_imag, a7_imag, a4_imag, a8_imag,
//                   u1_real, u2_real, u3_real, u4_real,
//                   u1_imag, u2_imag, u3_imag, u4_imag);

//     ureal_ptr = ureal + 64;
//     uimag_ptr = uimag + 64;

//     u1_real = ureal_ptr[0 * 16];
//     u2_real = ureal_ptr[1 * 16];
//     u3_real = ureal_ptr[2 * 16];
//     u4_real = ureal_ptr[3 * 16];

//     u1_imag = uimag_ptr[0 * 16];
//     u2_imag = uimag_ptr[1 * 16];
//     u3_imag = uimag_ptr[2 * 16];
//     u4_imag = uimag_ptr[3 * 16];

//     fft512_phase3(a1_real, a3_real, a5_real, a7_real, a2_real, a4_real, a6_real, a8_real,
//                   a1_imag, a3_imag, a5_imag, a7_imag, a2_imag, a4_imag, a6_imag, a8_imag,
//                   u1_real, u2_real, u3_real, u4_real,
//                   u1_imag, u2_imag, u3_imag, u4_imag);

//     ureal_ptr = ureal + 32;
//     uimag_ptr = uimag + 32;

//     u1_real = ureal_ptr[0 * 16];
//     u2_real = ureal_ptr[1 * 16];

//     u1_imag = uimag_ptr[0 * 16];
//     u2_imag = uimag_ptr[1 * 16];

//     fft512_phase4(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real,
//                   a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag,
//                   u1_real, u2_real, 
//                   u1_imag, u2_imag);


//     ureal_ptr = ureal + 16;
//     uimag_ptr = uimag + 16;

//     u1_real = ureal_ptr[0 * 16];

//     u1_imag = uimag_ptr[0 * 16];

//     fft512_phase5(a1_real, a2_real, a5_real, a6_real, a3_real, a4_real, a7_real, a8_real,
//                  a1_imag, a2_imag, a5_imag, a6_imag, a3_imag, a4_imag, a7_imag, a8_imag,
//                  u1_real, u1_imag);


//     simd_store(a1_real, dest_real + 0 * 4);
//     simd_store(a2_real, dest_real + 1 * 4);
//     simd_store(a5_real, dest_real + 2 * 4);
//     simd_store(a6_real, dest_real + 3 * 4);
//     simd_store(a3_real, dest_real + 4 * 4);
//     simd_store(a4_real, dest_real + 5 * 4);
//     simd_store(a7_real, dest_real + 6 * 4);
//     simd_store(a8_real, dest_real + 7 * 4);

//     simd_store(a1_imag, dest_imag + 0 * 4);
//     simd_store(a2_imag, dest_imag + 1 * 4);
//     simd_store(a5_imag, dest_imag + 2 * 4);
//     simd_store(a6_imag, dest_imag + 3 * 4);
//     simd_store(a3_imag, dest_imag + 4 * 4);
//     simd_store(a4_imag, dest_imag + 5 * 4);
//     simd_store(a7_imag, dest_imag + 6 * 4);
//     simd_store(a8_imag, dest_imag + 7 * 4);

// }

void simd_fft512_32_phase6to8(double* real, double* imag, double* ureal, double* uimag,
                              double* dest_real, double* dest_imag)
{
    doublev4 a1_real; simd_load(a1_real, real + 0 * 64);
    doublev4 a2_real; simd_load(a2_real, real + 1 * 64);
    doublev4 a3_real; simd_load(a3_real, real + 2 * 64);
    doublev4 a4_real; simd_load(a4_real, real + 3 * 64);
    doublev4 a5_real; simd_load(a5_real, real + 4 * 64);
    doublev4 a6_real; simd_load(a6_real, real + 5 * 64);
    doublev4 a7_real; simd_load(a7_real, real + 6 * 64);
    doublev4 a8_real; simd_load(a8_real, real + 7 * 64);
    
    doublev4 a1_imag; simd_load(a1_imag, imag + 0 * 64);
    doublev4 a2_imag; simd_load(a2_imag, imag + 1 * 64);
    doublev4 a3_imag; simd_load(a3_imag, imag + 2 * 64);
    doublev4 a4_imag; simd_load(a4_imag, imag + 3 * 64);
    doublev4 a5_imag; simd_load(a5_imag, imag + 4 * 64);
    doublev4 a6_imag; simd_load(a6_imag, imag + 5 * 64);
    doublev4 a7_imag; simd_load(a7_imag, imag + 6 * 64);
    doublev4 a8_imag; simd_load(a8_imag, imag + 7 * 64);

    double* ureal_ptr; 
    double* uimag_ptr;

    doublev4 u1_real;
    doublev4 u2_real;
    doublev4 u3_real;
    doublev4 u4_real;
    doublev4 u1_imag;
    doublev4 u2_imag;
    doublev4 u3_imag;
    doublev4 u4_imag;

    ureal_ptr = ureal + 8;
    uimag_ptr = uimag + 8;

    u1_real = ureal_ptr[0 * 2];
    u2_real = ureal_ptr[1 * 2];
    u3_real = ureal_ptr[2 * 2];
    u4_real = ureal_ptr[3 * 2];

    u1_imag = uimag_ptr[0 * 2];
    u2_imag = uimag_ptr[1 * 2];
    u3_imag = uimag_ptr[2 * 2];
    u4_imag = uimag_ptr[3 * 2];

    fft512_phase3(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real,
                  a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag,
                  u1_real, u2_real, u3_real, u4_real,
                  u1_imag, u2_imag, u3_imag, u4_imag);

    ureal_ptr = ureal + 4;
    uimag_ptr = uimag + 4;

    u1_real = ureal_ptr[0 * 2];
    u2_real = ureal_ptr[1 * 2];

    u1_imag = uimag_ptr[0 * 2];
    u2_imag = uimag_ptr[1 * 2];

    fft512_phase4(a1_real, a5_real, a2_real, a6_real, a3_real, a7_real, a4_real, a8_real,
                  a1_imag, a5_imag, a2_imag, a6_imag, a3_imag, a7_imag, a4_imag, a8_imag,
                  u1_real, u2_real, 
                  u1_imag, u2_imag);

    ureal_ptr = ureal + 2;
    uimag_ptr = uimag + 2;

    u1_real = ureal_ptr[0 * 2];

    u1_imag = uimag_ptr[0 * 2];

    fft512_phase5(a1_real, a5_real, a3_real, a7_real, a2_real, a6_real, a4_real, a8_real,
                  a1_imag, a5_imag, a3_imag, a7_imag, a2_imag, a6_imag, a4_imag, a8_imag,
                  u1_real, u1_imag);

    simd_store(a1_real, dest_real + 0 * 32);
    simd_store(a5_real, dest_real + 1 * 32);
    simd_store(a3_real, dest_real + 2 * 32);
    simd_store(a7_real, dest_real + 3 * 32);
    simd_store(a2_real, dest_real + 4 * 32);
    simd_store(a6_real, dest_real + 5 * 32);
    simd_store(a4_real, dest_real + 6 * 32);
    simd_store(a8_real, dest_real + 7 * 32);

    simd_store(a1_imag, dest_imag + 0 * 32);
    simd_store(a5_imag, dest_imag + 1 * 32);
    simd_store(a3_imag, dest_imag + 2 * 32);
    simd_store(a7_imag, dest_imag + 3 * 32);
    simd_store(a2_imag, dest_imag + 4 * 32);
    simd_store(a6_imag, dest_imag + 5 * 32);
    simd_store(a4_imag, dest_imag + 6 * 32);
    simd_store(a8_imag, dest_imag + 7 * 32);
}

void simd_fft512_32_phase9(double* real, double* imag, double* ureal, double* uimag,
                              double* dest_real, double* dest_imag)
{

    doublev4 a1_real; simd_load(a1_real, real + 0 * 4      );
    doublev4 a2_real; simd_load(a2_real, real + 1 * 4      );
    doublev4 a3_real; simd_load(a3_real, real + 2 * 4      );
    doublev4 a4_real; simd_load(a4_real, real + 3 * 4      );
    doublev4 a5_real; simd_load(a5_real, real + 0 * 4 + 256);
    doublev4 a6_real; simd_load(a6_real, real + 1 * 4 + 256);
    doublev4 a7_real; simd_load(a7_real, real + 2 * 4 + 256);
    doublev4 a8_real; simd_load(a8_real, real + 3 * 4 + 256);
    
    doublev4 a1_imag; simd_load(a1_imag, imag + 0 * 4      );
    doublev4 a2_imag; simd_load(a2_imag, imag + 1 * 4      );
    doublev4 a3_imag; simd_load(a3_imag, imag + 2 * 4      );
    doublev4 a4_imag; simd_load(a4_imag, imag + 3 * 4      );
    doublev4 a5_imag; simd_load(a5_imag, imag + 0 * 4 + 256);
    doublev4 a6_imag; simd_load(a6_imag, imag + 1 * 4 + 256);
    doublev4 a7_imag; simd_load(a7_imag, imag + 2 * 4 + 256);
    doublev4 a8_imag; simd_load(a8_imag, imag + 3 * 4 + 256);


    double* ureal_ptr; 
    double* uimag_ptr;

    doublev4 u1_real;
    doublev4 u1_imag;

    ureal_ptr = ureal + 1;
    uimag_ptr = uimag + 1;

    u1_real = ureal_ptr[0 * 1];
    u1_imag = uimag_ptr[0 * 1];


    fft512_phase5(a1_real, a2_real, a3_real, a4_real, a5_real, a6_real, a7_real, a8_real,
                  a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, a6_imag, a7_imag, a8_imag,
                  u1_real, u1_imag);

    simd_store(a1_real, dest_real + 0 * 4      );
    simd_store(a2_real, dest_real + 1 * 4      );
    simd_store(a3_real, dest_real + 2 * 4      );
    simd_store(a4_real, dest_real + 3 * 4      );
    simd_store(a5_real, dest_real + 0 * 4 + 256);
    simd_store(a6_real, dest_real + 1 * 4 + 256);
    simd_store(a7_real, dest_real + 2 * 4 + 256);
    simd_store(a8_real, dest_real + 3 * 4 + 256);

    simd_store(a1_imag, dest_imag + 0 * 4      );
    simd_store(a2_imag, dest_imag + 1 * 4      );
    simd_store(a3_imag, dest_imag + 2 * 4      );
    simd_store(a4_imag, dest_imag + 3 * 4      );
    simd_store(a5_imag, dest_imag + 0 * 4 + 256);
    simd_store(a6_imag, dest_imag + 1 * 4 + 256);
    simd_store(a7_imag, dest_imag + 2 * 4 + 256);
    simd_store(a8_imag, dest_imag + 3 * 4 + 256);

}

void simd_fft512(double* real, double* imag, double* ureal, double* uimag, double* dest_real, double* dest_imag)
{
    // 512*1 -> 16*32
    int i;
    // for(i = 0; i < 16; ++i)
    // {
    //     simd_fft512_32_phase1to5(real + i, imag + i, ureal + i, uimag + i, 
    //                             dest_real + i * 32, dest_imag + i * 32);
    // }
        simd_fft512_32_phase1to5(real, imag, ureal, uimag, 
                             dest_real, dest_imag);

    // 16*32 -> 2 * 256
    for(i = 0; i < 8; ++i)
    {
        simd_fft512_32_phase6to8(dest_real + 4 * i, dest_imag + 4 * i, 
                                 ureal, uimag, 
                                 real + i * 4, imag + i * 4);
    }

    for(i = 0; i < 8; ++i)
    {
        simd_fft512_32_phase6to8(dest_real + 4 * i + 32, dest_imag + 4 * i + 32, 
                                 ureal + 1, uimag + 1, 
                                 real + i * 4 + 256, imag + i * 4 + 256);
    }


    for(i = 0; i < 16; ++i)
    {
        simd_fft512_32_phase9(real + i * 16, imag + i * 16, 
                              ureal, uimag, 
                              dest_real + i * 16, dest_imag + i * 16);
    }
    for(i = 0; i < 512; ++i)
    {
        real[i] = dest_real[i];
        imag[i] = dest_imag[i];
    }
}

// //stride must be 4 de bei shu
// void simd_alltoall_doublev4r(int my_core, double* data, int block_sz, int data_count, int stride, 
//                     double* dest)
// {
//     int i, j, phase;
//     int my_rank = my_core % SIMD_COMM_SIZE;

//     //self
//     for(i = 0; i < data_count; ++i)
//     {
//         for(j = 0; j < block_sz; ++j)
//         {
//             doublev4 temp;
//             double* src_ptr  = data + my_rank * (block_sz * 4);
//             double* dest_ptr = dest + my_rank * (block_sz * 4);
//             simd_load (temp, src_ptr + stride * i + 4 * j);
//             simd_store(temp, dest_ptr + stride * i + 4 * j);
//         }
//     }

//     //other
//     for(phase = 1; phase < SIMD_COMM_SIZE; ++phase)
//     {
//         int partner = my_rank ^ phase;

//         for(i = 0; i < data_count; ++i)
//         {
//             for(j = 0; j < block_sz; ++j)
//             {
//                 volatile doublev4 temp;
//                 volatile doublev4 recv;
//                 // doublev4 temp;
//                 // doublev4 recv;
//                 double* src_ptr  = data + partner * (block_sz * 4);
//                 double* dest_ptr = dest + partner * (block_sz * 4);
//                 simd_load (temp, src_ptr + stride * i + 4 * j);
//                 // simd_putr(temp, partner);
//                 // recv = simd_getr(recv);
//                 asm volatile ("nop":::"memory");
//                 __builtin_sw64_putR(temp, partner);
//                 asm volatile ("nop":::"memory");
//                 recv = __builtin_sw64_getxR(recv);                
//                 simd_store(recv, dest_ptr + stride * i + 4 * j);
//             }
//         }
//         athread_syn(ROW_SCOPE, CPE_TOTAL_SYNC);

//     }
// }

// void simd_transpose_32_4x8r(int my_core, double* src, double* dest)
// {
//     int stride = 32;
//     int block_sz = 1;
//     int data_count = 4;
//     int it;

//     simd_alltoall_doublev4r(my_core, src, block_sz, data_count, stride, dest);

//     for(it = 0; it < SIMD_COMM_SIZE; ++it)
//     {
//         doublev4 r1; simd_load(r1, dest + it * 4 * block_sz             );
//         doublev4 r2; simd_load(r2, dest + it * 4 * block_sz + 1 * stride);
//         doublev4 r3; simd_load(r3, dest + it * 4 * block_sz + 2 * stride);
//         doublev4 r4; simd_load(r4, dest + it * 4 * block_sz + 3 * stride);

//         simd_4x4_transpose(r1, r2, r3, r4);

//         simd_store(r1, dest + it * 4 * block_sz             );
//         simd_store(r2, dest + it * 4 * block_sz + 1 * stride);
//         simd_store(r3, dest + it * 4 * block_sz + 2 * stride);
//         simd_store(r4, dest + it * 4 * block_sz + 3 * stride);
//     }
// }

// void simd_transpose_zyx2zxy_32(int my_core, double* src, double* dest)
// {
//     int stride = 32;
//     int block_sz = 1;
//     int data_count = 4;
//     int chunk_offset = stride * data_count;
//     int i;
//     for(i = 0; i < 4; ++i)
//     {
//         simd_transpose_32_4x8r(my_core, src + i * chunk_offset, dest + i * chunk_offset);
//     }
// }


// //stride must be 4 de bei shu
// void simd_alltoall_doublev4c(int my_core, int block_sz, int data_count, 
//                             double* src,  int src_stride, 
//                             double* dest, int dest_stride)
// {
//     int i, j, phase;
//     int my_rank = my_core / SIMD_COMM_SIZE;

//     //self
//     for(i = 0; i < data_count; ++i)
//     {
//         for(j = 0; j < block_sz; ++j)
//         {
//             doublev4 temp;
//             double* src_ptr  = src + my_rank * (block_sz * 4);
//             double* dest_ptr = dest + my_rank * (block_sz * 4);
//             simd_load (temp,  src_ptr +  src_stride * i + 4 * j);
//             simd_store(temp, dest_ptr + dest_stride * i + 4 * j);
//         }
//     }

//     //other
//     for(phase = 1; phase < SIMD_COMM_SIZE; ++phase)
//     {
//         int partner = my_rank ^ phase;

//         for(i = 0; i < data_count; ++i)
//         {
//             for(j = 0; j < block_sz; ++j)
//             {
//                 // doublev4 temp;
//                 // doublev4 recv;
//                 volatile doublev4 temp;
//                 volatile doublev4 recv;
//                 double* src_ptr  = src + partner * (block_sz * 4);
//                 double* dest_ptr = dest + partner * (block_sz * 4);
//                 simd_load (temp, src_ptr + src_stride * i + 4 * j);
//                 // simd_putc(temp, partner);
//                 // recv = simd_getc(recv);
//                 asm volatile ("nop":::"memory");
//                 __builtin_sw64_putC(temp, partner);
//                 asm volatile ("nop":::"memory");
//                 recv = __builtin_sw64_getxC(recv);
//                 simd_store(recv, dest_ptr + dest_stride * i + 4 * j);
//             }
//         }
//         athread_syn(COL_SCOPE, CPE_TOTAL_SYNC);

//     }
// }

// void simd_transpose_32_4x8c(int my_core, double* src, double* dest)
// {
//     int src_stride = 32 * 4;
//     int dest_stride = 32;
//     int block_sz = 1;
//     int data_count = 4;
//     int it;

//     simd_alltoall_doublev4c(my_core, block_sz, data_count, 
//                             src,  src_stride, 
//                             dest, dest_stride);

//     // simd_alltoall_doublev4c(my_core, src, block_sz, data_count, stride, dest);

//     for(it = 0; it < SIMD_COMM_SIZE; ++it)
//     {
//         doublev4 r1; simd_load(r1, dest + it * 4 * block_sz                  );
//         doublev4 r2; simd_load(r2, dest + it * 4 * block_sz + 1 * dest_stride);
//         doublev4 r3; simd_load(r3, dest + it * 4 * block_sz + 2 * dest_stride);
//         doublev4 r4; simd_load(r4, dest + it * 4 * block_sz + 3 * dest_stride);

//         simd_4x4_transpose(r1, r2, r3, r4);

//         simd_store(r1, dest + it * 4 * block_sz                  );
//         simd_store(r2, dest + it * 4 * block_sz + 1 * dest_stride);
//         simd_store(r3, dest + it * 4 * block_sz + 2 * dest_stride);
//         simd_store(r4, dest + it * 4 * block_sz + 3 * dest_stride);
//     }
// }

// void simd_transpose_zyx2yxz_32(int my_core, double* src, double* dest)
// {
//     // int stride = 32;
//     int block_sz = 1;
//     int data_count = 4;
//     int src_stride = 32;
//     int dest_stride = 32 * 4;
//     int i;
//     for(i = 0; i < 4; ++i)
//     {
//         simd_transpose_32_4x8c(my_core, src + i * src_stride, dest + i * dest_stride);
//     }
// }

void* nw_ldm_malloc(int sz)//will not see FBI WARNING during compile
{
    long ret = ldm_malloc(sz);
    return (void*)ret;
}

inline int cpe_int_pow(int a, int b)
{
    int i;
    int ret = 1;
    for(i = 0; i < b; i++)
        ret *= a;
    return ret;
}

int cpe_ilog2(int n)
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


//param
__thread_local struct fft_init_param param;

//control
__thread_local int my_core, core_x, core_y;
__thread_local int faker_my_core, faker_core_x, faker_core_y;
__thread_local volatile long local_flag[FLAG_SIZE];
__thread_local volatile long out_flag[FLAG_SIZE];
__thread_local long slave_flag_to_wait;
__thread_local volatile intv8 bcast_tmp256;
__thread_local volatile unsigned long get_reply, put_reply;
__thread_local volatile unsigned long get_reply_god, put_reply_god;
__thread_local volatile unsigned long get_reply_target, put_reply_target;
__thread_local volatile long stcc, edcc;
__thread_local volatile long stcc2, edcc2;
__thread_local volatile long stcc3, edcc3;

////FFF
__thread_local double* cpe_u_real_, *cpe_u_imag_;//total
//__thread_local double* cpe_u0_real_, *cpe_u0_imag_;//part
__thread_local double* cpe_u1_real_, *cpe_u1_imag_;//part
__thread_local double* cpe_u2_real_, *cpe_u2_imag_;//part
__thread_local double* cpe_twiddle_;//part
//__thread_local int cpe_cbw_, cpe_cbw_vol_;//checkerboard partition, a * a, of x * x
//__thread_local int arr_len_, part_len_;
__thread_local int x_stp_, x_edp_;//作用是标示当前CPE的负责的3D部分, of this process
__thread_local int y_stp_, y_edp_;
__thread_local int z_stp_, z_edp_;
__thread_local int x_stride_, y_stride_, z_stride_;//这个数值理论上是全局相同的。经过转置之后，就会进行改变，顺便还有标示谁是最高维度的作用。
//但是在之后的代码中，可能会使用不规则叠放，不一定按三维进行叠放。有可能先进行4*4 -> 1收缩，再按三维叠放
__thread_local int my_chk_cnt_;
__thread_local int* my_chk_x_;
__thread_local int* my_chk_y_;
__thread_local int* my_chk_z_;
__thread_local int cpe_ig_;//debug

#define SWAP_DB_PTR(ptra, ptrb) {double* tmp; tmp = ptra; ptra = ptrb; ptrb = tmp;}

//conghe profile module (￣△￣；)
__thread_local long pf_stcc_, pf_edcc_;
inline void stpf()
{
    RPCC(pf_stcc_);
}

inline void edpf(int n)
{
    RPCC(pf_edcc_);
    out_flag[n] += pf_edcc_ - pf_stcc_;
    pf_stcc_ = pf_edcc_;
}

//align
void* _nice_ptr_head(char* ptr, int alignment, int right_shft)
{
    long ptrl = (long)ptr;
    int rem = ptrl % alignment;
    if(rem > 0)
    {
        return ptr + (alignment - rem + right_shft);
    }
    else
    {
        return ptr + right_shft;
    }
}

inline static void wait_all_async_write_data()
{
    while(put_reply_god != put_reply_target);
}

inline static void wait_all_async_get_data()
{
    while(get_reply_god != get_reply_target);
}

void cpe0_bcast_256(int my_core, volatile intv8* bcast_data)
{
    bcast_tmp256 = *bcast_data;

    if (my_core / 8 == 0)
    {
        __builtin_sw64_putC(bcast_tmp256, 0x0000000F);
    } else {
        bcast_tmp256 = __builtin_sw64_getxC(bcast_tmp256);
    }

    if (my_core % 8 == 0)
    {
        asm volatile ("nop":::"memory");
        __builtin_sw64_putR(bcast_tmp256, 0x0000000F);
        *bcast_data = bcast_tmp256;
    } else {
        asm volatile ("nop":::"memory");
        *bcast_data = __builtin_sw64_getxR(bcast_tmp256);
    }
}

void standard_wait_flag(int my_core)
{
    if(my_core == 0)
    {
        //printf("SPE waiting for %ld\n", slave_flag_to_wait);

        while(1)
        {
            get_reply = 0;
            athread_get(
                PE_MODE, param.host_flag, (void*)&local_flag[0],
                sizeof(long) * FLAG_SIZE, (void*)(&get_reply), 0, 0, 0
            );
            while(get_reply != 1);
            asm volatile ("#nop":::"memory");
            if(local_flag[0] >= slave_flag_to_wait)
                break;
        }

        //printf("received flag %ld, action = %ld\n", local_flag[0], local_flag[KERNEL_ACTION]);

        slave_flag_to_wait++;
    }

//    cpe0_bcast_256(my_core, (intv8*)&local_flag[0]);
//    cpe0_bcast_256(my_core, (intv8*)&local_flag[4]);
//    cpe0_bcast_256(my_core, (intv8*)&local_flag[8]);
//    cpe0_bcast_256(my_core, (intv8*)&local_flag[12]);
    int i;
    for(i = 0; i < FLAG_SIZE / 4; i++)
    {
        cpe0_bcast_256(my_core, (intv8*)&local_flag[i * 4]);
    }
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
            PE_MODE, (void*)&out_flag[0], param.slave_flag,
            sizeof(long) * FLAG_SIZE, (void*)(&put_reply), 0, 0
        );
        while(put_reply != 1);

    }
    athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
}

void seq_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int my_core, int write_flag)
{
    int wcid;
    for (wcid = 63; wcid >= 0; wcid--)
    {
        if (my_core == wcid && write_flag != 0)
        {
            put_reply = 0;
            athread_put(
                PE_MODE, cpe_data_ptr, mpe_data_ptr,
                size, (void*)(&put_reply), 0, 0
            );
            while (put_reply != 1);
        }
        //athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    }
}

void async_write_data(void* mpe_data_ptr, void* cpe_data_ptr, int size, int my_core, int write_flag)
{
    int wcid;
    //for (wcid = 63; wcid >= 0; wcid--)
    //{
    //    if (my_core == wcid && write_flag != 0)
    //    {
            //put_reply = 0;
            athread_put(
                PE_MODE, cpe_data_ptr, mpe_data_ptr,
                size, (void*)(&put_reply_god), 0, 0
            );
            put_reply_target++;
            //while (put_reply != 1);
    //    }
    //    //athread_syn(ARRAY_SCOPE, CPE_TOTAL_SYNC);
    //}
}

void allin_allocate()
{
     cpe_u_real_ = (double*)nw_ldm_malloc(sizeof(double) * (param.nx));
     //cpe_u0_real_ = (double*)nw_ldm_malloc(sizeof(double) * (param.nx));
     cpe_u1_real_ = (double*)nw_ldm_malloc(sizeof(double) * (param.nx));
     cpe_u2_real_ = (double*)nw_ldm_malloc(sizeof(double) * (param.nx));
     cpe_u_imag_ = (double*)nw_ldm_malloc(sizeof(double) * (param.nx));
     //cpe_u0_imag_ = (double*)nw_ldm_malloc(sizeof(double) * (param.nx));
     cpe_u1_imag_ = (double*)nw_ldm_malloc(sizeof(double) * (param.nx));
     cpe_u2_imag_ = (double*)nw_ldm_malloc(sizeof(double) * (param.nx));
     cpe_twiddle_ = (double*)nw_ldm_malloc(sizeof(double) * (param.nx));
}

void allin_initial_transfer()
{
    int i;

    get_reply = 0;

    athread_get(
        PE_MODE, (void*)param.u_real, (void*)cpe_u_real_,
        sizeof(double) * param.nx, (void*)(&get_reply), 0, 0, 0
    );
    athread_get(
        PE_MODE, (void*)param.u_imag, (void*)cpe_u_imag_,
        sizeof(double) * param.nx, (void*)(&get_reply), 0, 0, 0
    );

    //while(get_reply != cpe_cbw_ * 3 + 2);
    while(get_reply != 2);

    //对于逆傅里叶变换，预计算共轭u
    for(i = 0; i < param.nx; i++)
    {
        cpe_u_imag_[i] = -cpe_u_imag_[i];
//        if(my_core == 0 && param.my_id == 0)
//            printf("CPE cpe_u_[%d] = %.8f %.8f\n", i, cpe_u_real_[i], cpe_u_imag_[i]);
    }
/*
    //distribute checksum
    // (o゜▽゜)o☆
    int tmp_chk_stp;
    get_reply = 0;
    athread_get(
        PE_MODE, (void*)&param.chk_stp[faker_my_core], (void*)&tmp_chk_stp,
        sizeof(int) * 1, (void*)(&get_reply), 0, 0, 0
    );
    athread_get(
        PE_MODE, (void*)&param.chk_cnt[faker_my_core], (void*)&my_chk_cnt_,
        sizeof(int) * 1, (void*)(&get_reply), 0, 0, 0
    );
    while(get_reply != 2);

    my_chk_x_ = (int*)nw_ldm_malloc(sizeof(int) * (my_chk_cnt_ + 1));
    my_chk_y_ = (int*)nw_ldm_malloc(sizeof(int) * (my_chk_cnt_ + 1));
    my_chk_z_ = (int*)nw_ldm_malloc(sizeof(int) * (my_chk_cnt_ + 1));

    if(my_chk_cnt_ > 0)
    {
        for(i = 0; i <= my_chk_cnt_; i++)
        {
            my_chk_z_[i] = 232000 + i;
            my_chk_y_[i] = 999000 + i;
            my_chk_x_[i] = 666000 + i;
        }

        get_reply = 0;
        athread_get(
            PE_MODE, (void*)&param.chk_xs[tmp_chk_stp], (void*)my_chk_x_,
            sizeof(int) * my_chk_cnt_, (void*)(&get_reply), 0, 0, 0
        );
        while(get_reply != 1);

        get_reply = 0;
        athread_get(
            PE_MODE, (void*)&param.chk_ys[tmp_chk_stp], (void*)my_chk_y_,
            sizeof(int) * my_chk_cnt_, (void*)(&get_reply), 0, 0, 0
        );
        while(get_reply != 1);

        get_reply = 0;
        athread_get(
            PE_MODE, (void*)&param.chk_zs[tmp_chk_stp], (void*)my_chk_z_,
            sizeof(int) * my_chk_cnt_, (void*)(&get_reply), 0, 0, 0
        );
        while(get_reply != 1);


        my_chk_z_[my_chk_cnt_] = 23232;
        my_chk_y_[my_chk_cnt_] = 2333;
        my_chk_x_[my_chk_cnt_] = 6666;

        //if(faker_my_core == 0)
        //{
        //    for(i = 0; i < my_chk_cnt_ + 1; i++)
        //    {
        //        printf("initial i = %d, x = %d, y = %d, z = %d\n", i, my_chk_x_[i], my_chk_y_[i], my_chk_z_[i]);
        //    }
        //}
    }*/
}


void sg_evolve(int d1)
{
    int i;
    for(i = 0; i < d1; i++)
    {
        cpe_u1_real_[i] = cpe_u1_real_[i] * cpe_twiddle_[i];
        cpe_u1_imag_[i] = cpe_u1_imag_[i] * cpe_twiddle_[i];
    }
    // int i;
    // int stride = 16;
    // int block_num = d1 / stride;
    // for(i = 0; i < d1; i += stride)
    // {
    //     doublev4 u0_real_1, u0_real_2, u0_real_3, u0_real_4;
    //     doublev4 u0_imag_1, u0_imag_2, u0_imag_3, u0_imag_4;

    //     doublev4 u1_real_1, u1_real_2, u1_real_3, u1_real_4;
    //     doublev4 u1_imag_1, u1_imag_2, u1_imag_3, u1_imag_4;

    //     doublev4 twiddle_1, twiddle_2, twiddle_3, twiddle_4;

    //     simd_load(u0_real_1, cpe_u1_real_ + i +  0);  simd_load(u0_imag_1, cpe_u1_imag_ + i +  0);
    //     simd_load(u0_real_2, cpe_u1_real_ + i +  4);  simd_load(u0_imag_2, cpe_u1_imag_ + i +  4);
    //     simd_load(u0_real_3, cpe_u1_real_ + i +  8);  simd_load(u0_imag_3, cpe_u1_imag_ + i +  8);
    //     simd_load(u0_real_4, cpe_u1_real_ + i + 12);  simd_load(u0_imag_4, cpe_u1_imag_ + i + 12);

    //     simd_load(twiddle_1, cpe_twiddle_ + i +  0);
    //     simd_load(twiddle_2, cpe_twiddle_ + i +  4);
    //     simd_load(twiddle_3, cpe_twiddle_ + i +  8);
    //     simd_load(twiddle_4, cpe_twiddle_ + i + 12);

    //     u0_real_1 = u0_real_1 * twiddle_1;
    //     u0_imag_1 = u0_imag_1 * twiddle_1;
    //     u0_real_2 = u0_real_2 * twiddle_2;
    //     u0_imag_2 = u0_imag_2 * twiddle_2;
    //     u0_real_3 = u0_real_3 * twiddle_3;
    //     u0_imag_3 = u0_imag_3 * twiddle_3;        
    //     u0_real_4 = u0_real_4 * twiddle_4;
    //     u0_imag_4 = u0_imag_4 * twiddle_4;

    //     u1_real_1 = u0_real_1;
    //     u1_imag_1 = u0_imag_1;
    //     u1_real_2 = u0_real_2;
    //     u1_imag_2 = u0_imag_2;
    //     u1_real_3 = u0_real_3;
    //     u1_imag_3 = u0_imag_3;
    //     u1_real_4 = u0_real_4;
    //     u1_imag_4 = u0_imag_4;

    //     // simd_store(u0_real_1, cpe_u0_real_ + i +  0);  simd_store(u0_imag_1, cpe_u0_imag_ + i +  0);
    //     // simd_store(u0_real_2, cpe_u0_real_ + i +  4);  simd_store(u0_imag_2, cpe_u0_imag_ + i +  4);
    //     // simd_store(u0_real_3, cpe_u0_real_ + i +  8);  simd_store(u0_imag_3, cpe_u0_imag_ + i +  8);
    //     // simd_store(u0_real_4, cpe_u0_real_ + i + 12);  simd_store(u0_imag_4, cpe_u0_imag_ + i + 12);

    //     simd_store(u1_real_1, cpe_u1_real_ + i +  0);  simd_store(u1_imag_1, cpe_u1_imag_ + i +  0);
    //     simd_store(u1_real_2, cpe_u1_real_ + i +  4);  simd_store(u1_imag_2, cpe_u1_imag_ + i +  4);
    //     simd_store(u1_real_3, cpe_u1_real_ + i +  8);  simd_store(u1_imag_3, cpe_u1_imag_ + i +  8);
    //     simd_store(u1_real_4, cpe_u1_real_ + i + 12);  simd_store(u1_imag_4, cpe_u1_imag_ + i + 12);

    // }
}

void total_fftz2(int l, int m, int n, int ny, double* x_real, double* x_imag, double* y_real, double* y_imag)
{
    //4 * 4 = 16?
#define x_idx(a) (a - 1)
#define y_idx(a) (a - 1)
#define u_idx(a) (a - 1)
    //dimension u(n), x(ny1,n), y(ny1,n)

    double u1_real,x11_real,x21_real;
    double u1_imag,x11_imag,x21_imag;


    int k,n1,li,lj,lk,ku,i,i11,i12,i21,i22;//ny = cpe_cbw_vol_

    n1 = n / 2;
    lk = cpe_int_pow(2, l - 1);
    li = cpe_int_pow(2, m - l);
    lj = 2 * lk;
    ku = li + 1;

    for(i = 0; i <= li - 1; i++)//16 -> 1
    {
        i11 = i * lk + 1;
        i12 = i11 + n1;
        i21 = i * lj + 1;
        i22 = i21 + lk;

        u1_real = cpe_u_real_[u_idx(ku + i)];
        u1_imag = cpe_u_imag_[u_idx(ku + i)];

        for(k = 0; k <= lk - 1; k++)//1 -> 16
        {
            x11_real = x_real[i11 + k - 1];
            x11_imag = x_imag[i11 + k - 1];
            x21_real = x_real[i12 + k - 1];
            x21_imag = x_imag[i12 + k - 1];

            //y[y_idx(j, i21 + k)] = cpe_add_dc_dc(x11, x21);
            //y[y_idx(j, i22 + k)] = cpe_multiply_dc_dc(u1, cpe_minus_dc_dc(x11, x21));

            y_real[i21 + k - 1] = x11_real + x21_real;
            y_imag[i21 + k - 1] = x11_imag + x21_imag;

            //minus
            x11_real -= x21_real;
            x11_imag -= x21_imag;//cpe_minus_dc_dc(x11, x21)

            y_real[i22 + k - 1] = u1_real * x11_real - u1_imag * x11_imag;
            y_imag[i22 + k - 1] = u1_real * x11_imag + u1_imag * x11_real;

//            if(param.my_id == 0 && my_core == 0 && cpe_ig_ == 0 && l == 1)
//            {
//                printf("CPE fftz2 dbg: %d %d %d %d %d, %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
//                       u_idx(ku + i), i11 + k - 1, i12 + k - 1, i21 + k - 1, i22 + k - 1,
//                       u1_real, u1_imag, x11_real, x11_imag, x21_real, x21_imag,
//                       y_real[i21 + k - 1], y_imag[i21 + k - 1], y_real[i22 + k - 1], y_imag[i22 + k - 1]
//                       );
//            }
        }
    }

#undef u_idx
#undef x_idx
#undef y_idx
}

//void total_fftz2(int l, int m, int n, int ny, double* x_real, double* x_imag, double* y_real, double* y_imag)
//{
//    //4 * 4 = 16?
//#define x_idx(a, b) ((a - 1) + (b - 1) * ny)
//#define y_idx(a, b) ((a - 1) + (b - 1) * ny)
//#define u_idx(a) (a - 1)
//    //dimension u(n), x(ny1,n), y(ny1,n)

//    double u1_real,x11_real,x21_real;
//    double u1_imag,x11_imag,x21_imag;

//    int k,n1,li,lj,lk,ku,i,j,i11,i12,i21,i22;//ny = cpe_cbw_vol_

//    n1 = n / 2;
//    lk = cpe_int_pow(2, l - 1);
//    li = cpe_int_pow(2, m - l);
//    lj = 2 * lk;
//    ku = li + 1;

//    for(i = 0; i <= li - 1; i++)//16 -> 1
//    {
//        i11 = i * lk + 1;
//        i12 = i11 + n1;
//        i21 = i * lj + 1;
//        i22 = i21 + lk;

//        u1_real = cpe_u_real_[u_idx(ku + i)];
//        u1_imag = cpe_u_imag_[u_idx(ku + i)];

//        for(k = 0; k <= lk - 1; k++)//1 -> 16
//        {
//            for(j = 1; j <= ny; j++)//cpe_cbw_vol_
//            {
//                x11_real = x_real[x_idx(j, i11 + k)];
//                x11_imag = x_imag[x_idx(j, i11 + k)];
//                x21_real = x_real[x_idx(j, i12 + k)];
//                x21_imag = x_imag[x_idx(j, i12 + k)];

//                //y[y_idx(j, i21 + k)] = cpe_add_dc_dc(x11, x21);
//                //y[y_idx(j, i22 + k)] = cpe_multiply_dc_dc(u1, cpe_minus_dc_dc(x11, x21));

//                y_real[y_idx(j, i21 + k)] = x11_real + x21_real;
//                y_imag[y_idx(j, i21 + k)] = x11_imag + x21_imag;

//                //minus
//                x11_real -= x21_real;
//                x11_imag -= x21_imag;//cpe_minus_dc_dc(x11, x21)

//                y_real[y_idx(j, i22 + k)] = u1_real * x11_real - u1_imag * x11_imag;
//                y_imag[y_idx(j, i22 + k)] = u1_real * x11_imag + u1_imag * x11_real;
//            }
//        }
//    }

//#undef u_idx
//#undef x_idx
//#undef y_idx
//}

void allin_cfftz(int m, int n, double* x_real, double* x_imag, double* y_real, double* y_imag)
{
    int i,l;

//    if(param.my_id == 0 && my_core == 0 && cpe_ig_ == 0)
//    {
//        for(i = 0; i < 512; i++)
//        {
//            printf("CPE u1 check allin_cfftz i = %d: %.8f %.8f %.8f %.8f\n", i, x_real[i], x_imag[i], y_real[i], y_imag[i]);
//        }
//    }

    for(l = 1; l <= m; l += 2)
    {
        stpf();
        total_fftz2(l, m, n, 1, x_real, x_imag, y_real, y_imag);
        edpf(PF_FFTZ2);

//        if(param.my_id == 0 && my_core == 0 && cpe_ig_ == 0 && l == 1)
//        {
////            printf("CPE u1 check l = %d y : %.8f %.8f\n", l, y_real[510], y_imag[510]);
//            for(i = 0; i < 512; i++)
//            {
//                printf("CPE u1 check l i = %d: %.8f %.8f %.8f %.8f\n", i, x_real[i], x_imag[i], y_real[i], y_imag[i]);
//            }
//        }

        if(l == m)//early break, copy to output to x
        {
            for(i = 0; i < param.nx; i++)
            {
                x_real[i] = y_real[i];
                x_imag[i] = y_imag[i];
            }
            edpf(PF_CFFTZ_COPY);
            break;
        }

        total_fftz2(l + 1, m, n, 1, y_real, y_imag, x_real, x_imag);
        edpf(PF_FFTZ2);

//        if(param.my_id == 0 && my_core == 0 && cpe_ig_ == 0 && l == 1)
//        {
////            printf("CPE u1 check l = %d x : %.8f %.8f\n", l, x_real[510], x_imag[510]);
//            for(i = 0; i < 512; i++)
//            {
//                printf("CPE u1 check lp1 i = %d: %.8f %.8f %.8f %.8f\n", i, x_real[i], x_imag[i], y_real[i], y_imag[i]);
//            }
//        }
    }
}

void local_cffts1(double* data_real, double* data_imag, double* buffer_real, double* buffer_imag, int d1)//这个函数的最终输出是x，y只是起到了缓存的作用
{
    int logd1;
    logd1 = cpe_ilog2(d1);

    //在这里排序成可以向量化的格式
    stpf();
    int j;
    for(j = 0; j < d1; j++)
    {
        buffer_real[j] = data_real[j];
        buffer_imag[j] = data_imag[j];
    }
    edpf(PF_CFFTS_COPY);

    allin_cfftz(logd1, d1, buffer_real, buffer_imag, data_real, data_imag);

    stpf();
    for(j = 0; j < d1; j++)
    {
        data_real[j] = buffer_real[j];
        data_imag[j] = buffer_imag[j];
    }
    edpf(PF_CFFTS_COPY);
}

void write_out_my_u(double* ux_real, double* ux_imag, double* mpe_ux_real, double* mpe_ux_imag, int mpe_offset)
{
    put_reply = 0;
    athread_put(
                PE_MODE, (void*)&ux_real[0], (void*)&mpe_ux_real[mpe_offset],
            sizeof(double) * param.nx, (void*)(&put_reply), 0, 0
            );
    athread_put(
                PE_MODE, (void*)&ux_imag[0], (void*)&mpe_ux_imag[mpe_offset],
            sizeof(double) * param.nx, (void*)(&put_reply), 0, 0
            );
    while(put_reply != 2);
}

void read_in_my_u(double* ux_real, double* ux_imag, double* mpe_ux_real, double* mpe_ux_imag, int mpe_offset)
{
    get_reply = 0;
    athread_get(
                PE_MODE, (void*)&mpe_ux_real[mpe_offset], (void*)&ux_real[0],
            sizeof(double) * param.nx, (void*)(&get_reply), 0, 0, 0
            );
    athread_get(
                PE_MODE, (void*)&mpe_ux_imag[mpe_offset], (void*)&ux_imag[0],
            sizeof(double) * param.nx, (void*)(&get_reply), 0, 0, 0
            );
    while(get_reply != 2);
}

void read_in_twiddle(int mpe_offset)
{
    get_reply = 0;
    athread_get(
                PE_MODE, (void*)&param.twiddle[mpe_offset], (void*)&cpe_twiddle_[0],
            sizeof(double) * param.nx, (void*)(&get_reply), 0, 0, 0
            );
    while(get_reply != 1);
}

void db_fft_1()
{
    int ig, loop_length, jkid, in_stride, out_stride, my_start, my_len, fft_d1;
    double* ptr_in_real, *ptr_out_real;
    double* ptr_in_imag, *ptr_out_imag;

    int remain_length = local_flag[REMAIN_POINT];
    loop_length = local_flag[GROUP_SIZE];

    in_stride = local_flag[IN_STRIDE];
    out_stride = local_flag[OUT_STRIDE];

    ptr_in_real = (double*)(local_flag[IN_PTR_REAL]);
    ptr_out_real = (double*)(local_flag[IN_PTR_REAL]);
    ptr_in_imag = (double*)(local_flag[IN_PTR_IMAG]);
    ptr_out_imag = (double*)(local_flag[IN_PTR_IMAG]);

    fft_d1 = local_flag[FFT_D1];

//    if(param.my_id == 0 && my_core == 0)
//    {
//        printf("enter sg_fft\n");
//        printf("ftw = %ld, loop_length = %d, in_stride = %d, out_stride = %d, d1 = %d, rp = %d\nptr = %x %x %x %x\n"
//              , slave_flag_to_wait, loop_length, in_stride, out_stride, fft_d1, remain_length,
//               ptr_in_real, ptr_in_imag, ptr_out_real, ptr_out_imag);
//    }

    my_start = BLOCK_LOW(my_core, 64, loop_length * 64 + remain_length);
    my_len = BLOCK_SIZE(my_core, 64, loop_length * 64 + remain_length);


    for(ig = 0; ig < loop_length + 1; ig++)
    {
        cpe_ig_ = ig;
        if(my_len == ig)
        {
            break;
        }

        jkid = my_start + ig;

        read_in_my_u(cpe_u1_real_, cpe_u1_imag_, ptr_in_real, ptr_in_imag, jkid * in_stride);

//        if(param.my_id == 0 && my_core == 0 && ig == 0)
//        {
//            printf("CPE u1 check in : %.8f %.8f\n", cpe_u1_real_[510], cpe_u1_imag_[510]);
//        }

        // local_cffts1(cpe_u1_real_, cpe_u1_imag_, cpe_u2_real_, cpe_u2_imag_, fft_d1);
        simd_fft512(cpe_u1_real_, cpe_u1_imag_, 
                    cpe_u_real_, cpe_u_imag_, 
                    cpe_u2_real_, cpe_u2_imag_);

//        if(param.my_id == 0 && my_core == 0 && ig == 0)
//        {
//            printf("CPE u1 check out: %.8f %.8f\n", cpe_u1_real_[510], cpe_u1_imag_[510]);
//        }

        write_out_my_u(cpe_u1_real_, cpe_u1_imag_, ptr_out_real, ptr_out_imag, jkid * out_stride);
    }
    standard_write_flag(my_core);
}

void db_fft_2()
{
    int ig, loop_length, jkid, in_stride, out_stride, my_start, my_len, fft_d1;
    double* ptr_in_real, *ptr_out_real;
    double* ptr_in_imag, *ptr_out_imag;

    int remain_length = local_flag[REMAIN_POINT];
    loop_length = local_flag[GROUP_SIZE];

    in_stride = local_flag[IN_STRIDE];
    out_stride = local_flag[OUT_STRIDE];

    ptr_in_real = (double*)(local_flag[OUT_PTR_REAL]);
    ptr_out_real = (double*)(local_flag[IN_PTR_REAL]);
    ptr_in_imag = (double*)(local_flag[OUT_PTR_IMAG]);
    ptr_out_imag = (double*)(local_flag[IN_PTR_IMAG]);

    fft_d1 = local_flag[FFT_D1];

//    if(param.my_id == 0 && my_core == 0)
//    {
//        printf("enter sg_fft\n");
//        printf("ftw = %ld, loop_length = %d, in_stride = %d, out_stride = %d, d1 = %d, rp = %d\nptr = %x %x %x %x\n"
//              , slave_flag_to_wait, loop_length, in_stride, out_stride, fft_d1, remain_length,
//               ptr_in_real, ptr_in_imag, ptr_out_real, ptr_out_imag);
//    }

    my_start = BLOCK_LOW(my_core, 64, loop_length * 64 + remain_length);
    my_len = BLOCK_SIZE(my_core, 64, loop_length * 64 + remain_length);


    for(ig = 0; ig < loop_length + 1; ig++)
    {
        cpe_ig_ = ig;
        if(my_len == ig)
        {
            break;
        }

        jkid = my_start + ig;

        read_in_my_u(cpe_u1_real_, cpe_u1_imag_, ptr_in_real, ptr_in_imag, jkid * in_stride);

//        if(param.my_id == 0 && my_core == 0 && ig == 0)
//        {
//            printf("CPE u1 check in : %.8f %.8f\n", cpe_u1_real_[510], cpe_u1_imag_[510]);
//        }

        // local_cffts1(cpe_u1_real_, cpe_u1_imag_, cpe_u2_real_, cpe_u2_imag_, fft_d1);
        simd_fft512(cpe_u1_real_, cpe_u1_imag_,
                    cpe_u_real_, cpe_u_imag_,
                    cpe_u2_real_, cpe_u2_imag_);

//        if(param.my_id == 0 && my_core == 0 && ig == 0)
//        {
//            printf("CPE u1 check out: %.8f %.8f\n", cpe_u1_real_[510], cpe_u1_imag_[510]);
//        }

        write_out_my_u(cpe_u1_real_, cpe_u1_imag_, ptr_out_real, ptr_out_imag, jkid * out_stride);
    }
    standard_write_flag(my_core);
}

void evv_fft()
{
    //read in twiddle
    int ig, loop_length, jkid, in_stride, out_stride, my_start, my_len, fft_d1;
    double* ptr_in_real, *ptr_out_real;
    double* ptr_in_imag, *ptr_out_imag;

    int remain_length = local_flag[REMAIN_POINT];
    loop_length = local_flag[GROUP_SIZE];

    in_stride = local_flag[IN_STRIDE];
    out_stride = local_flag[OUT_STRIDE];

    ptr_in_real = (double*)(local_flag[IN_PTR_REAL]);
    ptr_out_real = (double*)(local_flag[OUT_PTR_REAL]);
    ptr_in_imag = (double*)(local_flag[IN_PTR_IMAG]);
    ptr_out_imag = (double*)(local_flag[OUT_PTR_IMAG]);

    fft_d1 = local_flag[FFT_D1];

//    if(param.my_id == 0 && my_core == 0)
//    {
//        printf("enter sg_fft\n");
//        printf("ftw = %ld, loop_length = %d, in_stride = %d, out_stride = %d, d1 = %d, rp = %d\nptr = %x %x %x %x\n"
//              , slave_flag_to_wait, loop_length, in_stride, out_stride, fft_d1, remain_length,
//               ptr_in_real, ptr_in_imag, ptr_out_real, ptr_out_imag);
//    }

    my_start = BLOCK_LOW(my_core, 64, loop_length * 64 + remain_length);
    my_len = BLOCK_SIZE(my_core, 64, loop_length * 64 + remain_length);


    for(ig = 0; ig < loop_length + 1; ig++)
    {
        cpe_ig_ = ig;
        if(my_len == ig)
        {
            break;
        }

        jkid = my_start + ig;

        read_in_my_u(cpe_u1_real_, cpe_u1_imag_, ptr_in_real, ptr_in_imag, jkid * in_stride);
        read_in_twiddle(jkid * in_stride);

        sg_evolve(fft_d1);

        write_out_my_u(cpe_u1_real_, cpe_u1_imag_, ptr_in_real, ptr_in_imag, jkid * out_stride);

//        if(param.my_id == 0 && my_core == 0 && ig == 0)
//        {
//            printf("CPE u1 check in : %.8f %.8f\n", cpe_u1_real_[510], cpe_u1_imag_[510]);
//        }

        // local_cffts1(cpe_u1_real_, cpe_u1_imag_, cpe_u2_real_, cpe_u2_imag_, fft_d1);
        simd_fft512(cpe_u1_real_, cpe_u1_imag_,
                    cpe_u_real_, cpe_u_imag_,
                    cpe_u2_real_, cpe_u2_imag_);

//        if(param.my_id == 0 && my_core == 0 && ig == 0)
//        {
//            printf("CPE u1 check out: %.8f %.8f\n", cpe_u1_real_[510], cpe_u1_imag_[510]);
//        }

        write_out_my_u(cpe_u1_real_, cpe_u1_imag_, ptr_out_real, ptr_out_imag, jkid * out_stride);
    }
    standard_write_flag(my_core);
}

void simd_16x16_transpose(double* buffer)
{
    doublev4 left1, left2, left3, left4;
    doublev4 right1, right2, right3, right4;

    int i, j;
    for(i = 0; i < 4; ++i)
    {
        simd_load(left1, buffer + 0 * 16 + i * (4 + 4 * 16));
        simd_load(left2, buffer + 1 * 16 + i * (4 + 4 * 16));
        simd_load(left3, buffer + 2 * 16 + i * (4 + 4 * 16));
        simd_load(left4, buffer + 3 * 16 + i * (4 + 4 * 16));
        simd_4x4_transpose(left1, left2, left3, left4);
        simd_store(left1, buffer + 0 * 16 + i * (4 + 4 * 16));
        simd_store(left2, buffer + 1 * 16 + i * (4 + 4 * 16));
        simd_store(left3, buffer + 2 * 16 + i * (4 + 4 * 16));
        simd_store(left4, buffer + 3 * 16 + i * (4 + 4 * 16));
    }
    for(i = 0; i < 3; ++i)
    {
        for(j = i + 1; j < 4; ++j)
        {
            simd_load(left1, buffer + 0 * 16 + j * (4 * 16) + i * 4);
            simd_load(left2, buffer + 1 * 16 + j * (4 * 16) + i * 4);
            simd_load(left3, buffer + 2 * 16 + j * (4 * 16) + i * 4);
            simd_load(left4, buffer + 3 * 16 + j * (4 * 16) + i * 4);
            simd_load(right1, buffer + 0 * 16 + i * (4 * 16) + j * 4);
            simd_load(right2, buffer + 1 * 16 + i * (4 * 16) + j * 4);
            simd_load(right3, buffer + 2 * 16 + i * (4 * 16) + j * 4);
            simd_load(right4, buffer + 3 * 16 + i * (4 * 16) + j * 4);

            simd_4x4_transpose(left1, left2, left3, left4);
            simd_4x4_transpose(right1, right2, right3, right4);

            simd_store(left1, buffer + 0 * 16 + i * (4 * 16) + j * 4);
            simd_store(left2, buffer + 1 * 16 + i * (4 * 16) + j * 4);
            simd_store(left3, buffer + 2 * 16 + i * (4 * 16) + j * 4);
            simd_store(left4, buffer + 3 * 16 + i * (4 * 16) + j * 4);
            simd_store(right1, buffer + 0 * 16 + j * (4 * 16) + i * 4);
            simd_store(right2, buffer + 1 * 16 + j * (4 * 16) + i * 4);
            simd_store(right3, buffer + 2 * 16 + j * (4 * 16) + i * 4);
            simd_store(right4, buffer + 3 * 16 + j * (4 * 16) + i * 4);
        }
    }



}



// buffer size must greater or eqaul to 256
void transpose_kernel_16x16(double* buffer, double* in_hostptr, double* out_hostptr, int load_offset, int load_stride,
                            int store_offset, int store_stride)
{
    int i;
    volatile int rpl = 0;
    for(i = 0; i < 16; ++i)
    {
        // athread_get(
        //     PE_MODE, (void*)(param.mat + load_offset + i * load_stride), (void*)(buffer + i * 16),
        //     sizeof(double) * 16, (void*)(&rpl), 0, 0, 0
        // );
        athread_get(
            PE_MODE, (void*)(in_hostptr + load_offset + i * load_stride), (void*)(buffer + i * 16),
            sizeof(double) * 16, (void*)(&rpl), 0, 0, 0
        );
    }
    while(rpl != 16);

    simd_16x16_transpose(buffer);

    rpl = 0;
    for(i = 0; i < 16; ++i)
    {
        // athread_put(
        //     PE_MODE,  (void*)(buffer + i * 16), (void*)(param.mat_out + store_offset + i * store_stride),
        //     sizeof(double) * 16, (void*)(&rpl), 0, 0
        // );
        athread_put(
            PE_MODE,  (void*)(buffer + i * 16), (void*)(out_hostptr + store_offset + i * store_stride),
            sizeof(double) * 16, (void*)(&rpl), 0, 0
        );
    }
    while(rpl != 16);


}

// buffer size must greater or eqaul to 256
// d1 = d2 !!!
void slave_tanspose_zyx2zxy(int my_core, double* buffer, double* in_hostptr, double* out_hostptr, int d1, int d2, int d3)
{
    int row_id = my_core / SIMD_COMM_SIZE;
    int col_id = my_core % SIMD_COMM_SIZE;
    int local_d = (d1 >> 3);

    int z_offset = d1 * d2;
    int load_offset = row_id * local_d * d1 + col_id * local_d;
    int load_stride = d1;

    int store_offset = col_id * local_d * d1 + row_id * local_d;
    int store_stride = d1;

    int block_sz = (local_d >> 4);
    int i, row, col;

    for(i = 0; i < d3; ++i)
    {

        for(row = 0; row < block_sz; ++row)
        {
            for(col = 0; col < block_sz; ++col)
            {
                transpose_kernel_16x16(buffer, in_hostptr, out_hostptr,
                                        load_offset + row * 16 * d1 + col * 16, load_stride,
                                        store_offset + col * 16 * d1 + row * 16, store_stride);
            }
        }

        load_offset += z_offset;
        store_offset += z_offset;

    }
    standard_write_flag(my_core);
}

//void allin_fft_rev()
//{
//    //too difficult!
//    //x, y, z
//    allin_cffts1(cpe_u1_real_, cpe_u1_imag_, cpe_u2_real_, cpe_u2_imag_);

//    //transpose
//    // allin_cute_transpose(cpe_u1_real_, cpe_u1_imag_, cpe_u2_real_, cpe_u2_imag_);
//    stpf();
//    simd_transpose_zyx2zxy_32(my_core, cpe_u1_real_, cpe_u2_real_);
//    simd_transpose_zyx2zxy_32(my_core, cpe_u1_imag_, cpe_u2_imag_);
//    edpf(PF_TP1);

//    //y in main now
//    y_stp_ = 0;
//    y_stp_ = 31;
//    x_stp_ = core_x * 4;
//    x_edp_ = x_stp_ + 4 - 1;
//    z_stp_ = core_y * 4;
//    z_edp_ = z_stp_ + 4 - 1;
//    y_stride_ = 1;
//    x_stride_ = param.nx;
//    z_stride_ = param.nx * param.nx;

//    //y x z
//    allin_cffts1(cpe_u2_real_, cpe_u2_imag_, cpe_u1_real_, cpe_u1_imag_);

//    //transpose
////    allin_cute_transpose(cpe_u2_real_, cpe_u2_imag_, cpe_u1_real_, cpe_u1_imag_);
//    stpf();
//    simd_transpose_zyx2yxz_32(my_core, cpe_u2_real_, cpe_u1_real_);
//    simd_transpose_zyx2yxz_32(my_core, cpe_u2_imag_, cpe_u1_imag_);
//    edpf(PF_TP2);
////    swap_transpose_cpe(sizeof(double) * param.npc, cpe_u1_real_, cpe_u2_real_, my_core);
////    swap_transpose_cpe(sizeof(double) * param.npc, cpe_u1_imag_, cpe_u2_imag_, my_core);


//    //z in main now
//    z_stp_ = 0;
//    z_edp_ = 31;
//    x_stp_ = faker_core_y * cpe_cbw_;
//    x_edp_ = x_stp_ + cpe_cbw_ - 1;
//    y_stp_ = faker_core_x * cpe_cbw_;
//    y_edp_ = y_stp_ + cpe_cbw_ - 1;
//    z_stride_ = 1;
//    // x_stride_ = param.nx;
//    // y_stride_ = param.nx * param.nx;
//    y_stride_ = param.nx;
//    x_stride_ = param.nx * param.nx;

//    //not z x y, but z y x
//    allin_cffts1(cpe_u1_real_, cpe_u1_imag_, cpe_u2_real_, cpe_u2_imag_);

//}

//void allin_checksum(int it)
//{
//    if(my_chk_cnt_ > 0)
//    {
//        //进行本地的checksum
//        dfc chk;
//        int i, offset_local;
//        chk.real = 0.0;
//        chk.imag = 0.0;

//        for(i = 0; i < my_chk_cnt_; i++)
//        {
//            /*
//            offset_local = (my_chk_x_[i] - x_stp_) * x_stride_
//                    +      (my_chk_y_[i] - y_stp_) * y_stride_
//                    +      (my_chk_z_[i] - z_stp_) * z_stride_;*/
//            //the offset above is global offset.
////            offset_local = (my_chk_x_[i] - x_stp_) * param.nx
////                    +      (my_chk_y_[i] - y_stp_) * param.nx * cpe_cbw_
////                    +      (my_chk_z_[i] - z_stp_) * 1;
//            offset_local = (my_chk_x_[i] - x_stp_) * param.nx  * cpe_cbw_
//                    +      (my_chk_y_[i] - y_stp_) * param.nx
//                    +      (my_chk_z_[i] - z_stp_) * 1;

//            chk.real += cpe_u1_real_[offset_local];
//            chk.imag += cpe_u1_imag_[offset_local];
//        }

//        //printf("local it %d, off = %d\n", it, my_core * param.iter + it - 1);

//        put_reply = 0;
//        athread_put(
//            PE_MODE, (void*)&chk, (void*)&param.buf_chk[faker_my_core * param.iter + it - 1],
//            sizeof(dfc) * 1, (void*)(&put_reply), 0, 0
//        );

//        while(put_reply != 1);
//    }
//    standard_write_flag(my_core);
//    //standard_wait_flag(my_core);//等待checksum完成
//    //不用等 (。・・)ノ


//    /*
//    //全都按顺序写到主核上面
//    //write
//    write_out_my_ux(cpe_u1_);

//    //if(my_core == 0 && param.my_id == 0)
//    //    printf("CPE write data out, wait mpe checksum\n");

//    standard_write_flag(my_core);//开始checksum

//    standard_wait_flag(my_core);//等待checksum完成

//    //if(my_core == 0 && param.my_id == 0)
//    //    printf("CPE find that mpe finished checksum\n");*/
//}

void allin_fft_iter()
{
    x_stp_ = 0;
    x_stp_ = 31;//yz ckb p
    y_stp_ = core_x * 4;
    y_edp_ = y_stp_ + 4 - 1;
    z_stp_ = core_y * 4;
    z_edp_ = z_stp_ + 4 - 1;
    x_stride_ = 1;
    y_stride_ = param.nx;
    z_stride_ = param.nx * param.nx;
    //this changes after transpose<(￣ˇ￣)/

    stpf();
    allin_allocate();
    allin_initial_transfer();
    edpf(PF_INITIAL_TRANSFER);

    //int iter;

//    if(my_core == 0 && param.my_id == 0)
//        printf("begin loop\n");
/*
    for(iter = 1; iter <= param.iter; iter++)
    {
        stpf();
        allin_evolve();
        edpf(PF_EVV);
        //standard_wait_flag(my_core);//wait MPE evv
        //if(my_core == 0 && param.my_id == -1)
        //    printf("CPE iter = %d, found mpe finished evv\n", iter);
        //read_in_my_u2();//start

        allin_fft_rev();//reverse 3D fft in a single CG

        stpf();
        allin_checksum(iter);//暂时交给主核
        edpf(PF_CHK);
        //standard_wait_flag(my_core);//wait MPE for checksum
    }*/
}

void slave_sym_tanspose_asyn(int my_core, double* buffer,
                        double* in_hostptr, double* out_hostptr,
                        int m,
                        int load_offset, int load_stride,
                        int store_offset, int store_stride)
{
    int row_id = my_core / SIMD_COMM_SIZE;
    int col_id = my_core % SIMD_COMM_SIZE;
    int local_d = (m >> 3);

    load_offset  += row_id * local_d *  load_stride + col_id * local_d;

    store_offset += col_id * local_d * store_stride + row_id * local_d;

    int block_sz = (local_d >> 4);
    int i, row, col;



    for(row = 0; row < block_sz; ++row)
    {
        for(col = 0; col < block_sz; ++col)
        {
            transpose_kernel_16x16(buffer, in_hostptr, out_hostptr,
                                     load_offset + row * 16 *  load_stride + col * 16, load_stride,
                                    store_offset + col * 16 * store_stride + row * 16, store_stride);
        }
    }

}


//buffer size must greater or equal to 256
void c_transpose2_local(int my_core, double* buffer/*buffer size must greater or equal to 256
*/, double* in_hostptr, double* out_hostptr, int d1, int d2, int d3)
{
    int load_zoffset = d1 * d2;
    int store_zoffset = d2;

    int load_offset = 0;
    int store_offset = 0;

    int load_stride = d1;
    int store_stride = d2 * d3;

    int i;
    for(i = 0; i < d3; ++i)
    {
        slave_sym_tanspose_asyn(my_core, buffer,
                        in_hostptr, out_hostptr,
                        d1,
                        load_offset, load_stride,
                        store_offset, store_stride);

        load_offset += load_zoffset;
        store_offset += store_zoffset;
    }
    standard_write_flag(my_core);

}

void tp_x_y()
{
    int d1, d2, d3;
    d1 = local_flag[FFT_D1];
    d2 = local_flag[FFT_D2];
    d3 = local_flag[FFT_D3];
    double* out_ptr, *in_ptr;
    in_ptr = (double*)(local_flag[IN_PTR_REAL]);
    out_ptr = (double*)(local_flag[OUT_PTR_REAL]);
    slave_tanspose_zyx2zxy(my_core, cpe_u2_real_, in_ptr, out_ptr, d1, d2, d3);

    in_ptr = (double*)(local_flag[IN_PTR_IMAG]);
    out_ptr = (double*)(local_flag[OUT_PTR_IMAG]);
    slave_tanspose_zyx2zxy(my_core, cpe_u2_real_, in_ptr, out_ptr, d1, d2, d3);
}

void tp_local()
{
    int d1, d2, d3;
    d1 = local_flag[FFT_D1];
    d2 = local_flag[FFT_D2];
    d3 = local_flag[FFT_D3];
    double* out_ptr, *in_ptr;
    in_ptr = (double*)(local_flag[OUT_PTR_REAL]);
    out_ptr = (double*)(local_flag[SUP_PTR_REAL]);
    c_transpose2_local(my_core, cpe_u2_real_, in_ptr, out_ptr, d1, d2, d3);

    in_ptr = (double*)(local_flag[OUT_PTR_IMAG]);
    out_ptr = (double*)(local_flag[SUP_PTR_IMAG]);
    c_transpose2_local(my_core, cpe_u2_real_, in_ptr, out_ptr, d1, d2, d3);
}

//buffer size must greater or equal to d1
void slave_transpose2_finish(int my_core, double* buffer,
                        double* in_hostptr, double* out_hostptr,
                        int d1, int d2, int d3)
{
    //double complex xin(n2, n1/np2, 0:np2-1), xout(n2*np2, n1/np2)
    //d1, d2, d3

    int total_yz_p = d2 * d3;
    //div them between cores
    int my_start, my_len, ig, inner_part, part_total, part_len;
    int yz_id, y_id, z_id;

    my_start = BLOCK_LOW(my_core, 64, total_yz_p);
    my_len = BLOCK_SIZE(my_core, 64, total_yz_p);
    part_len = local_flag[IN_STRIDE];
    part_total = d1 / part_len;

    int off_inner, off_outer_in, off_outer_out;

    int xout_d1 = d1 * d3;

    for(ig = 0; ig < my_len; ig++)
    {
        yz_id = my_start + ig;
        y_id = yz_id % d2;
        z_id = yz_id / d2;

//        off_outer_in = (y_id + z_id * d2) * d1;//z is the ** dim
//        off_outer_out = (y_id * d2 + z_id) * d1;//y is the ** dim

        off_outer_in = (y_id + z_id * d2) * d1;
        off_outer_out = y_id * xout_d1 + z_id * d1;

        for(inner_part = 0; inner_part < part_total; inner_part++)
        {
            off_inner = inner_part * part_len;

            get_reply = 0;
            athread_get(
                PE_MODE, (void*)(in_hostptr + (off_inner + off_outer_in)), (void*)buffer,
                sizeof(double) * part_len, (void*)(&get_reply), 0, 0, 0
            );
            while(get_reply != 1);

            put_reply = 0;
            athread_put(
                PE_MODE, (void*)buffer, (void*)(out_hostptr + (off_inner + off_outer_out)),
                sizeof(double) * part_len, (void*)(&put_reply), 0, 0
            );
            while(put_reply != 1);
        }
    }

    standard_write_flag(my_core);
}

void tp_finish()
{
    int d1, d2, d3;
    d1 = local_flag[FFT_D1_FNS];
    d2 = local_flag[FFT_D2_FNS];
    d3 = local_flag[FFT_D3_FNS];
    double* out_ptr, *in_ptr;
//    in_ptr = (double*)(local_flag[IN_PTR_REAL]);
//    out_ptr = (double*)(local_flag[OUT_PTR_REAL]);
    in_ptr = (double*)(local_flag[OUT_PTR_REAL]);
    out_ptr = (double*)(local_flag[IN_PTR_REAL]);
    slave_transpose2_finish(my_core, cpe_u2_real_, in_ptr, out_ptr, d1, d2, d3);

//    in_ptr = (double*)(local_flag[IN_PTR_IMAG]);
//    out_ptr = (double*)(local_flag[OUT_PTR_IMAG]);
    in_ptr = (double*)(local_flag[OUT_PTR_IMAG]);
    out_ptr = (double*)(local_flag[IN_PTR_IMAG]);
    slave_transpose2_finish(my_core, cpe_u2_real_, in_ptr, out_ptr, d1, d2, d3);
}

void cpe_athread_daemon(void *_param)
{
    my_core = athread_get_id(-1);

    int i, j, k;
    stcc = 0;
    edcc = 0;

    param = *((struct fft_init_param*) _param);

    core_x = my_core % 8;
    core_y = my_core / 8;

    faker_core_x = core_y;
    faker_core_y = core_x;
    faker_my_core = faker_core_y * 8 + faker_core_x;
//    faker_core_x = core_x;
//    faker_core_y = core_y;
//    faker_my_core = my_core;

//    cpe_cbw_ = param.nx / 8;//8 * 8 = 64
//    cpe_cbw_vol_ = cpe_cbw_ * cpe_cbw_;

    if(my_core == 0 && param.my_id == 0)
    {
        printf("@@@ MESSAGE FROM CPE @@@\n");
        printf("iter = %d, nx = %d, ny = %d, nz = %d, npc = %d\n", param.iter, param.nx, param.ny, param.nz, param.npc);
        printf("ptrs: %x %x %x %x %x\n%x %x %x %x %x %x\n",
               param.host_flag, param.slave_flag, param.u_real, param.u_imag, param.twiddle,
               param.chk_cnt, param.chk_stp, param.chk_xs, param.chk_ys, param.chk_zs, param.buf_chk);
//        printf("arr_len_ = %d, part_len_ = %d, cpe_cbw_ = %d, cpe_cbw_vol_ = %d\n",
//               arr_len_, part_len_, cpe_cbw_, cpe_cbw_vol_);
        printf("@@@ MESSAGE FROM CPE @@@\n");
    }

    //must do this
    //if(my_core == 0 && param.my_id == 0)
    //    printf("nx = %d, ny = %d, nz = %d, id = %d, host_flag = %x, slave_flag = %x\n",
    //       param.nx, param.ny, param.nz, param.my_id, param.host_flag, param.slave_flag);

    for (i = 0; i < FLAG_SIZE; i++){ local_flag[i] = 0; out_flag[i] = 0;}
    slave_flag_to_wait = 1;

    get_reply_god = 0, put_reply_god = 0;//the Jesus do not like waiting for someone
    get_reply_target = 0, put_reply_target = 0;

    stpf();
    allin_allocate();
    allin_initial_transfer();
    edpf(PF_INITIAL_TRANSFER);

    if(my_core == 0 && param.my_id == 0)
    {
        printf("CPE ENTER LOOP\n");
    }

    ////////RPCC(stcc);

    while(1)
    {
        standard_wait_flag(my_core);
        ////////RPCC(edcc);
        //out_flag[PF_OUT_WAIT_FLAG] += edcc - stcc;
        stcc = edcc;

        if (local_flag[KERNEL_ACTION] == EXIT_FLAG)
        {
            standard_write_flag(my_core);
            standard_write_flag(my_core);
            break;
        }

        if(local_flag[KERNEL_ACTION] == EVV_FFT_FLAG)
        {
            evv_fft();
            tp_local();
        }

        if(local_flag[KERNEL_ACTION] == DB_FFT_FLAG)
        {
            tp_finish();
            db_fft_1();
            tp_x_y();
            db_fft_2();
        }

//        if(local_flag[KERNEL_ACTION] == TP_LOCAL_FLAG)
//        {
//            tp_local();
//        }

//        if(local_flag[KERNEL_ACTION] == TP_FINISH_FLAG)
//        {
//            tp_finish();
//        }

//        if(local_flag[KERNEL_ACTION] == TP_X_Y_FLAG)
//        {
//            tp_x_y();
//        }
    }
}
