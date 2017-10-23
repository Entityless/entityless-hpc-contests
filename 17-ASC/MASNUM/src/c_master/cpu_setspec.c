#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../c_header/c_public_var.h"
#include "../c_header/cpu_implsch.h"

void c_setspec2_(int *_ia_start, int *_ia_end, int *_ic_start, int *_ic_end)
{
    int ic_start, ic_end, ia_start, ia_end;
    ic_start = (*_ic_start);
    ic_end   = (*_ic_end);
    ia_start = (*_ia_start);
    ia_end   = (*_ia_end);
    int index_tmp;

    int ia, ic, iaic_index;

    const float gama=3.3, sq3=1.7320508075688772935274463415059;
    float vx, vy, ww, xj, xj0, arlfa, wsj, wkj, theta0, sinth, costh, wk0, wf0, ws0, wl, sigma, alpha;
    int j, k;


    for (ic = ic_start; ic <= ic_end; ic++)
    {
        iaic_index = _iaic_index(ia_start - 1, ic);
        for (ia = ia_start; ia <= ia_end; ia++)
        {
            iaic_index++;
            if (_nsp[iaic_index] < 2.0) continue;

            vx = _wx[iaic_index];
            vy = _wy[iaic_index];

            _w[iaic_index] = sqrt(vx * vx + vy * vy);

            if(_w[iaic_index] <= 0.0f)
                _w[iaic_index] = 0.9f;

            xj0 = 200.0 * 1000.0;
            xj = _g * xj0 / (_w[iaic_index] * _w[iaic_index]);
            arlfa = (0.076 * (powf(xj, -0.4f))) / _pi;
            wsj = 22. * (powf(xj, -0.33f)) * _g / _w[iaic_index];
            wkj = (wsj * wsj) / _g;

            for(j = 1; j <= _jl; j++)
            {
                costh = cos_theta[j - 1];
                sinth = sin_theta[j - 1];
                for(k = 1; k <= _kl; k++)
                {
                    wk0 = _wk[k - 1];
                    index_tmp = (ic - _iys) * (_kldp1 * ix_size) + (ia - _ixs) * _kldp1;
                    wf0 = _wf[index_tmp + k - 1];
                    ws0 = _zpi * wf0;
                    wl = vx * costh + vy * sinth;

                    if(wl > 0.)
                    {
                        if(ws0 <= wsj)
                        {
                            sigma = 0.07;
                        }
                        else
                        {
                            sigma = 0.09;
                        }
                        alpha = arlfa / powf(wk0, 4) * expf(-1.25*(wkj/wk0)*(wkj/wk0))
                                * powf(gama, exp(-0.5 * powf(((1.0 - ws0 / wsj) / sigma), 2)))
                                * powf((wl / _w[iaic_index]), 2);
                    }
                    else
                        alpha = 0.0;

                    index_tmp = _e_4d_index(k, j, ia, ic);

                    //if(my_rank == 0 && k == 10 && j == 10)
                    //{
                    //    printf("c:\n  k, j, ia, ic: %d, %d, %d, %d\n", k, j, ia, ic);

                    //}

                    _e[index_tmp] = max(alpha, _small);
                }
            }
        }
    }
}