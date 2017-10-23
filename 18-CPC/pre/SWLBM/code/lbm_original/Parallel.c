#include "Argument.h"

/*------------------------------------
 *      MPI communication Part
 *-----------------------------------*/

void SetMPI(MPI_Comm mycomm,
		    int *dims, 
			int *coords)
{
    tmp_coords[0] = coords[0] - 1; 
	tmp_coords[1] = coords[1];
    if(tmp_coords[0] >= 0) {
        MPI_Cart_rank(mycomm, tmp_coords, &left_nbr);
	}

    tmp_coords[0] = coords[0] + 1;
	tmp_coords[1] = coords[1];
    if(tmp_coords[0] < dims[0]) {
        MPI_Cart_rank(mycomm, tmp_coords, &right_nbr);
	}

    tmp_coords[0] = coords[0];
	tmp_coords[1] = coords[1] + 1;
    if(tmp_coords[1] < dims[1]) {
        MPI_Cart_rank(mycomm, tmp_coords, &up_nbr);
	}

    tmp_coords[0] = coords[0];
	tmp_coords[1] = coords[1] - 1;
    if(tmp_coords[1] >= 0) {
        MPI_Cart_rank(mycomm, tmp_coords, &down_nbr);
	}

    tmp_coords[0] = coords[0] - 1;
	tmp_coords[1] = coords[1] - 1;
    if(tmp_coords[0] >= 0 && tmp_coords[1] >= 0) {
        MPI_Cart_rank(mycomm, tmp_coords, &ld_nbr);
	}

    tmp_coords[0] = coords[0] - 1;
	tmp_coords[1] = coords[1] + 1;
    if(tmp_coords[0] >= 0 && tmp_coords[1] < dims[1]) {
        MPI_Cart_rank(mycomm, tmp_coords, &lu_nbr);
	}

    tmp_coords[0] = coords[0] + 1;
	tmp_coords[1] = coords[1] - 1;
    if(tmp_coords[0] < dims[0] && tmp_coords[1] >= 0) {
        MPI_Cart_rank(mycomm, tmp_coords, &rd_nbr);
	}

    tmp_coords[0] = coords[0] + 1;
	tmp_coords[1] = coords[1] + 1;
    if(tmp_coords[0] < dims[0] && tmp_coords[1] <  dims[1]) {
        MPI_Cart_rank(mycomm, tmp_coords, &ru_nbr);
	}
}

void bounce_send_init(int X,
					  int Y,
					  int Z,
					  int Xst, 
					  int Xed, 
					  int Yst, 
					  int Yed, 
					  int x_sec,
					  int y_sec,
					  int corrent, 
					  int other, 
					  Real *****nodes,
					  Real ***temp_left_send,
					  Real ***temp_right_send, 
					  Real ***temp_up_send, 
					  Real ***temp_down_send,
					  Real **temp_ld_send, 
					  Real **temp_lu_send, 
					  Real **temp_rd_send, 
					  Real **temp_ru_send)
{
    int i, j, k, l;

    if(Xst != 0) {
		for (i = Yst; i < Yed; i++) {
			for (j = 0; j < Z; j++) {
				for(l = 0; l < 19; l++) {
					temp_left_send[i - Yst][j][l] = nodes[other][1][i - Yst + 1][j][l];
				}
			}
		}
	}

    if(Xed != X) { 
		for (i = Yst; i < Yed; i++) {
			for (j = 0; j < Z; j++) {
				for(l = 0; l < 19; l++) {
					   temp_right_send[i - Yst][j][l] = nodes[other][x_sec][i - Yst + 1][j][l];
				}
			}
		}
	}

    if(Yst != 0) {
		for (i = Xst; i < Xed; i++) {
			for (j = 0; j < Z; j++) {
				for(l = 0; l < 19; l++) {
					 temp_down_send[i - Xst][j][l] = nodes[other][i - Xst + 1][1][j][l];
				}
			}
		}
	}

    if(Yed != Y) { 
		for (i = Xst; i < Xed; i++ ) {
			for (j = 0; j < Z; j++) {
				for(l = 0; l < 19; l++) {
					temp_up_send[i - Xst][j][l] = nodes[other][i - Xst + 1][y_sec][j][l];
				}
			}
		}
	}

    if (Xst != 0 && Yst != 0) {
        for (j = 0; j < Z; j++) {
            for(l = 0; l < 19; l++) {
                temp_ld_send[j][l] = nodes[other][1][1][j][l];
			}
		}
	}

    if (Xst != 0 && Yed != Y) {
        for (j = 0; j < Z; j++) {
            for(l = 0; l < 19; l++) {
                temp_lu_send[j][l] = nodes[other][1][y_sec][j][l];
			}
		}
	}

    if (Xed != X && Yst != 0) {
        for (j = 0; j < Z; j++) {
            for(l = 0; l < 19; l++) {
                temp_rd_send[j][l] = nodes[other][x_sec][1][j][l];
			}
		}
	}

    if (Xed != X && Yed != Y) {
        for (j = 0; j < Z; j++) {
            for(l = 0; l < 19; l++) {
                temp_ru_send[j][l] = nodes[other][x_sec][y_sec][j][l];
			}
		}
	}
}

void bounce_communicate(MPI_Comm mycomm, 
		                int *dims, 
					    int *coords, 
						int x_sec, 
						int y_sec,  
						int Z,
						int *l_count,
						MPI_Status *sta,
						MPI_Request *req,
						Real ***temp_left_send, 
						Real ***temp_right_send, 
						Real ***temp_up_send, 
						Real ***temp_down_send, 
						Real ***temp_left, 
						Real ***temp_right, 
						Real ***temp_up, 
						Real ***temp_down,
						Real **temp_lu_send, 
						Real **temp_ld_send, 
						Real **temp_ru_send, 
						Real **temp_rd_send, 
						Real **temp_lu, 
						Real **temp_ld, 
						Real **temp_ru, 
						Real **temp_rd)
{
	int i, count = 0;

     if(left_nbr != -1) {
         MPI_Irecv(&temp_left[0][0][0], 
				   y_sec * Z * 19, 
				   MPI_FLOAT,
				   left_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
         MPI_Isend(&temp_left_send[0][0][0], 
				   y_sec * Z * 19, 
				   MPI_FLOAT,
				   left_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
     }

     if(right_nbr != -1) {
         MPI_Irecv(&temp_right[0][0][0], 
				   y_sec * Z * 19, 
				   MPI_FLOAT,
				   right_nbr, 
				   0, 
				   MPI_COMM_WORLD,
				   &req[count++]);
         MPI_Isend(&temp_right_send[0][0][0], 
				   y_sec * Z * 19, 
				   MPI_FLOAT,
				   right_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
     }

     if(up_nbr != -1) {
         MPI_Irecv(&temp_up[0][0][0], 
				   x_sec * Z * 19, 
				   MPI_FLOAT,
				   up_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
         MPI_Isend(&temp_up_send[0][0][0], 
				   x_sec * Z * 19, 
				   MPI_FLOAT,
				   up_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
     }

     if(down_nbr != -1) {
         MPI_Irecv(&temp_down[0][0][0], 
				   x_sec * Z * 19, 
				   MPI_FLOAT,
				   down_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
         MPI_Isend(&temp_down_send[0][0][0], 
				   x_sec * Z * 19, 
				   MPI_FLOAT,
				   down_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
     }

     if(lu_nbr != -1) {
         MPI_Irecv(&temp_lu[0][0], 
				   Z * 19, 
				   MPI_FLOAT,
				   lu_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
         MPI_Isend(&temp_lu_send[0][0], 
				   Z * 19, 
				   MPI_FLOAT,
				   lu_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
     }

     if(ld_nbr != -1) {
         MPI_Irecv(&temp_ld[0][0], 
				   Z * 19, 
				   MPI_FLOAT,
				   ld_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
         MPI_Isend(&temp_ld_send[0][0], 
				   Z * 19, 
				   MPI_FLOAT,
				   ld_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]); 
     }

     if(ru_nbr != -1) {
         MPI_Irecv(&temp_ru[0][0], 
				   Z * 19, 
				   MPI_FLOAT,
				   ru_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
         MPI_Isend(&temp_ru_send[0][0], 
				   Z * 19, 
				   MPI_FLOAT,
				   ru_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
     }

     if(rd_nbr != -1) {
         MPI_Irecv(&temp_rd[0][0], 
				   Z * 19, 
				   MPI_FLOAT,
				   rd_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
         MPI_Isend(&temp_rd_send[0][0], 
				   Z * 19, 
				   MPI_FLOAT,
				   rd_nbr, 
				   0, 
				   MPI_COMM_WORLD, 
				   &req[count++]);
     }
	 (*l_count) = count;
}

void bounce_update(int X,
					int Y,
					int Z,
					int Xst, 
		            int Xed, 
					int Yst, 
					int Yed, 
					int myrank,
					int x_sec, 
					int y_sec,
					int other,
					Real *****nodes,
				    Real ***temp_left, 
					Real ***temp_right, 
					Real ***temp_up, 
					Real ***temp_down,
				    Real **temp_ld, 
					Real **temp_lu, 
					Real **temp_rd, 
					Real **temp_ru)
{
    int i, j, k, l;

    if(Xst != 0) { 
		for (i = Yst; i < Yed; i++) {
			for (j = 0; j < Z; j++) {
				for(l = 0; l < 19; l++) {
					nodes[other][0][i - Yst + 1][j][l] = temp_left[i - Yst][j][l];
				}
			}
		}
	}

    if(Xed != X) { 
		for (i = Yst; i < Yed; i++) {
			for (j = 0; j < Z; j++) {
				for(l = 0; l < 19; l++) {
					nodes[other][x_sec + 1][i - Yst + 1][j][l] = temp_right[i - Yst][j][l];
				}
			}
		}
	}

    if(Yst != 0) { 
		for (i = Xst; i < Xed; i++) {
			for (j = 0; j < Z; j++) {
				for(l = 0; l < 19; l++) {
					nodes[other][i - Xst + 1][0][j][l] = temp_down[i - Xst][j][l];
				}
			}
		}
	}

    if(Yed != Y) {
		for (i = Xst; i < Xed; i++) {
			for (j = 0; j < Z; j++) {
				for(l = 0; l < 19; l++) {
					nodes[other][i - Xst + 1][y_sec + 1][j][l] = temp_up[i - Xst][j][l];
				}
			}
		}
	}

    if (Xst != 0 && Yst != 0) {
        for (j = 0; j < Z; j++) {
            for(l = 0; l < 19; l++) {
                nodes[other][0][0][j][l] = temp_ld[j][l];
            }
		}
	}

    if (Xst != 0 && Yed != Y) {
        for (j = 0; j < Z; j++) {
            for(l = 0; l < 19; l++) {
                nodes[other][0][y_sec + 1][j][l] = temp_lu[j][l];
            }
		}
	}

    if (Xed != X && Yst != 0) {
        for (j = 0; j < Z; j++) {
            for(l = 0; l < 19; l++) {
                nodes[other][x_sec + 1][0][j][l] = temp_rd[j][l];
            }
		}
	}

    if (Xed != X && Yed != Y) {
        for (j = 0; j < Z; j++) {
            for(l = 0; l < 19; l++) {
                nodes[other][x_sec + 1][y_sec + 1][j][l] = temp_ru[j][l];
            }
		}
	}
}
