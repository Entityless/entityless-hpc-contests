#F77      = mpiifort
F77      = mpif90
CC       = mpicc
SCC      = sw5cc

EXE = weno7.out

CFLAGS   = -O3 -lm -fno-math-errno -msimd -DTEST -DSWAP -DCHECK -OPT:ieee_arith=1
SCFLAGS  = -O3 -slave -msimd -lm_slave -fno-math-errno -OPT:ieee_arith=1
LF77OPTS = -O3 -lm_slave -lm -fno-math-errno -msimd -OPT:ieee_arith=1
FFLAGS   = -O3 -fno-math-errno -msimd -cpp -DTEST -DSWAP -DCHECK -OPT:ieee_arith=1

##intel, debug
#CFLAGS   = -O2 -lm -fno-math-errno -g -cpp -DTEST -DSWAP -DCHECK -OPT:ieee_arith=1
#LF77OPTS = -O3 -g
#FFLAGS   = -O2 -fno-math-errno -g -DTEST -DSWAP -DCHECK -OPT:ieee_arith=1

SRC =  parameters.f90  main.f90 \
Weno.f90       \
Interfaces.f90


CINC = \
slave_kernel.h c_public_var.h c_public_const.h

COBJS = \
c_public_var.o mpe_weno.o


SOBJS = slave_kernel.o

VPATH = ./f90_source/:./c_header/:./c_master/:./c_slave/

OBJS  = ${SRC:.f90=.o}

.SUFFIXES: .f90 .o

.f90.o:
	$(F77) $(FFLAGS) -c $<

$(EXE)  :$(OBJS) $(COBJS) $(SOBJS)
	$(F77) $(LF77OPTS) -o $@ $(OBJS) $ $(QOBJS) $(COBJS) $(SOBJS)

$(SOBJS) : %.o : %.c
	$(SCC) $(SCFLAGS) -c -o $@ $<

$(COBJS) : %.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXE) *.o work.* *.pc *.mod *.rmod
