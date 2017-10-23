#ifndef __Argument_h__
#define __Argument_h__
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <string.h>

#ifdef NDEBUG
#define DEBUG(M, ...)
#else
#define DEBUG(M, ...) fprintf(stderr, "[%s  %s] [DEBUG] (%d): " M "", __DATE__, __TIME__, __LINE__, ##__VA_ARGS__)
#endif

#ifdef NMLOG

#define MLOG(M, ...)

#else

#define MLOG(M, ...) {if(myrank == 0) { \
	                  time_t t; \
					  struct tm *tmif; \
	                  t = time(NULL); \
	                  tmif = localtime(&t); \
	                  fprintf(stderr, "[%d-%2d-%2d] [%2d:%2d:%2d] [MLOG] " M "",  \
							  tmif->tm_year + 1900, \
							  tmif->tm_mon + 1, \
							  tmif->tm_mday, \
							  tmif->tm_hour, \
							  tmif->tm_min, \
							  tmif->tm_sec, ##__VA_ARGS__); }} 
#endif

#define TLOG(M, ...)  {if(myrank == 0) { \
	                  fprintf(stderr, "[%s  %s] [TLOG ] (%d): " M "", __DATE__, __TIME__, __LINE__, ##__VA_ARGS__) }}

#define OLOG(RANK, M, ...)  \
	                  fprintf(stderr, "[%s  %s] [OLOG ] (RANK%d  %d): " M "", __DATE__, __TIME__, RANK, __LINE__, ##__VA_ARGS__)



#define NUM_DIMS 2
#define NUM_THREADS 4
#define OUTNAME "result/outlbm"

/*-----------------------------*
 *  different types of cells
 *-----------------------------*/
const int FLUID  = 0,
	      NOSLIP = 1,
		  VELOCITY = 2,
		  BOUNCE = 3,
		  PRESSURE = 4;

typedef float Real ;
int X, Y, Z, STEPS;

double t1, t2;
unsigned long r1, r2;

/*------------------------------------------*
 * Declare some constants and globals ...
 * -----------------------------------------*/
int v_log;
Real CSmago;
Real LDC_VELOCITY[3], V0[3] = {0.0};
Real RE;
Real R;
Real L_cylinder;
Real Rmax_suboff;
Real m_per_ft;
Real ox_suboff,oy_suboff,oz_suboff;
Real rho0;
Real nu;
Real dx, dt;
Real pressureScale, refP;
Real omega;     // viscosity of the fluid, 0..2
Real exp_vel;

int tmp_coords[2];
int right_nbr = -1,
    left_nbr  = -1,
    up_nbr    = -1,
    down_nbr  = -1,
    lu_nbr    = -1,
    ru_nbr    = -1,
    ld_nbr    = -1,
    rd_nbr    = -1;

/*----------------------------------------------------*
 *   the weight for the equilibrium distribution
 *   -------------------------------------------------*/
const Real w[19] = {
		(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),(1./18.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),
		(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./36.),(1./3.) };

/*-----------------------------------------------------*
 *     convenience variables for handling
 *     the distribution functions in a loop
 *     the 19 lattice vectors
 *-----------------------------------------------------*/
const int e_x[19] = {0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1,0};
const int e_y[19] = {1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0,0};
const int e_z[19] = {0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1,0};

/*-----------------------------------------------------*
 * the index of the inverse for each lattice vector
 * ----------------------------------------------------*/
const int dfInv[19] = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 18};

typedef struct ImageInfo
{
	int x;
	int y;
	Real u_x;
	Real u_y;
	Real u_z;
	Real MagV;
	Real rho;
} ImageInfo;

static inline unsigned long rpcc()
{
     unsigned long time;
     asm("rtc %0": "=r" (time) : );
     return time;
}

void TIME_ST();

void TIME_ED();

void setParameter();

Real **array2DF(int row, int col);

void arrayFree2DF(Real **array);

int **array2DI(int row, int col);

void arrayFree2DI(int **array);

Real ***array3DF(int row, int col, int channel);

void arrayFree3DF(Real ***array);

int ***array3DI(int row, int col, int channel);

void arrayFree3DI(int ***array);

int ****array4DI(int row1, int row, int col, int channel);

void arrayFree4DI(int ****array);

Real *****array5DF(int row2, int row1, int row, int col, int channel);

void arrayFree5DF(Real *****array);

void write_result(int X, int Y, int Z, int s, double times, Real **image);

void init_Flag_local(int X, int Y, int Z, int Xst, int Xed, int Yst, int Yed, int ***flags);

void collide(Real *****nodes, int ***flags, int Xst, int Xed, int Yst, int Yed, int nz, int current);

void writeImage_tecplot(int X, int Y, int Z, int Xst, int Xed, int Yst, int Yed, int other, int myrank, Real **local_image, Real *****nodes);

void init_Pointer(int ***flags, Real *****nodes, int ****walls, int Xst, int Xed, int Yst, int Yed, int nz,int x_sec, int y_sec, Real rho,Real* u);

void stream(Real *****nodes, int ****walls, int ***flags, int Xst, int Xed, int Yst, int Yed, int nz, int current, int other);

void SetMPI(MPI_Comm mycomm, int *dims, int *coords);

void INITINPUT(int X, int Y, int Z, int Xst, int Xed, int Yst, int Yed, int x_sec, int y_sec, 
		       int myrank, int size, char *user, int *local_rankinfo, int **rankinfo, int ***flags);

int OUTPUT(int X, int Y, int Z, int Xst, int Xed, int Yst, int Yed, int s, 
		   int myrank, int size, int other, int x_sec, int y_sec, char *user, 
		   Real **local_image, Real **image, int **rankinfo, Real *****nodes);

void bounce_send_init_test (int X, int Y, int Z,
							int Xst, int Xed, int Yst, int Yed, 
							int x_sec, int y_sec, int corrent, int other, 
							Real *****nodes,
						    Real ***temp_left_send,
							Real ***temp_right_send, 
							Real ***temp_up_send, 
							Real ***temp_down_send,
						    Real **temp_ld_send, 
							Real **temp_lu_send, 
							Real **temp_rd_send, 
							Real **temp_ru_send);

void bounce_communicate_Itest(MPI_Comm mycomm, int *dims, int *coords, 
							  int x_sec, int y_sec, int Z, int *l_count,
							  MPI_Status *sta, MPI_Request *req,
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
							  Real **temp_rd);

void bounce_update_test(int X, int Y, int Z,
						int Xst, int Xed, int Yst, int Yed, 
						int myrank, int x_sec, int y_sec,
						int other,
						Real *****nodes,
					    Real ***temp_left, 
						Real ***temp_right, 
						Real ***temp_up, 
						Real ***temp_down,
					    Real **temp_ld, 
						Real **temp_lu, 
						Real **temp_rd, 
						Real **temp_ru);

#endif
