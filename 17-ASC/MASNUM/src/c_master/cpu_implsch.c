#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../c_header/c_public_var.h"
#include "../c_header/cpu_implsch.h"

void c_implsch_(int *_ia_start, int *_ia_end, int *_ic_start, int *_ic_end)
{
    //if(my_rank == 0)
    //    printf("exit(1); //c_implsch_ \n");
    //exit(1);
    /*
    if(halo_ips_collect_flag != 0)
        return;

    int ic, ia, ic_start, ic_end, ia_start, ia_end, iaic_index;

    ic_start = (*_ic_start);
    ic_end   = (*_ic_end);
    ia_start = (*_ia_start);
    ia_end   = (*_ia_end);

    for (ic = ic_start; ic <= ic_end; ic++)
    {
        iaic_index = _iaic_index(ia_start - 1, ic);
        for (ia = ia_start; ia <= ia_end; ia++)
        {
            iaic_index++;
            if (_nsp[iaic_index] != 1.0) continue;

            //if(halo_ips_collect_flag == 0)
            //{
                halo_ips_ia_set[halo_ips_marine_count] = ia;
                halo_ips_ic_set[halo_ips_marine_count] = ic;
                halo_ips_marine_count++;
            //}
        }
    }*/
}

void collect_halo_ips_marine_finalize_()
{
    return;
    /*if(halo_ips_collect_flag == 0)//整理收集成果
    {
        halo_ips_collect_flag = 1;
    }*/

    int index_tmp, ia, ic, kj, i;

    for(i = 0; i < halo_ppg_marine_count; i++)
    {
        ia = halo_ppg_ia_set[i];
        ic = halo_ppg_ic_set[i];

        index_tmp = _e_4d_index(1, 1, ia, ic);


        for(kj = 1; kj <= dim3; kj++)
        {
            _ee[index_tmp] = _e[index_tmp];
            index_tmp++;
        }
    }
}
