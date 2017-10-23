#ifndef PUBLIC_CONST_H
#define PUBLIC_CONST_H

#define CCPS 1450000000

#define KIU 8
#define CB_GS 8
#define CB_GTT (CB_GS * CB_GS)
//checkerboard div group size

//#define hchjity long
//must be long, for some reason

//谁说的
#define hchjity long

#define FLAG_SIZE 32

static inline int BLOCK_SIZE(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) + ((n % total_blocks > block_id) ? 1 : 0);
}

static inline int BLOCK_LOW(int block_id, int total_blocks, int n)
{
    return (n / total_blocks) * block_id + ((n % total_blocks > block_id) ? block_id : n % total_blocks);
}

#endif
