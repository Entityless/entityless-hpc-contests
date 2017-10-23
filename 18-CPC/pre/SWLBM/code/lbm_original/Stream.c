#include "Argument.h"

/*------------------------------------
 *      Main Computing Part
 *-----------------------------------*/
void stream(Real *****nodes,
		    int ****walls,
			int ***flags,
			int Xst,
			int Xed,
			int Yst,
			int Yed,
			int nz,
		    int current,
			int other)
{
	int i, j, k, l;
	int inv;

	for(i = Xst; i < Xed; i++) {
		for(j = Yst; j < Yed; j++) {
			for(k = 0; k < nz; k++) {
				if(flags[i - Xst + 1][j - Yst + 1][k] == FLUID) {
					for(l = 0; l < 19; l++) {
						inv = dfInv[l];
						nodes[current][i - Xst + 1][j - Yst + 1][k][l] = nodes[other][i - Xst + 1 + e_x[inv]][j - Yst + 1 + e_y[inv]][k + e_z[inv]][l];
					}
				}	
				if(flags[i - Xst + 1][j - Yst + 1][k] == BOUNCE) {
					for(l = 0; l < 19; l++) {
						inv = dfInv[l];
						if(walls[i - Xst][j - Yst][k][l]) {
							nodes[current][i - Xst + 1][j - Yst + 1][k][l] = nodes[other][i - Xst + 1][j - Yst + 1][k][inv];
						} else {
							nodes[current][i - Xst + 1][j - Yst + 1][k][l] = nodes[other][i - Xst + 1 + e_x[inv]][j - Yst + 1 + e_y[inv]][k + e_z[inv]][l];
						}
					}
				}	
			}
		}
	}
}
