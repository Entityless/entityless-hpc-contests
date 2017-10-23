#F77      = mpiifort
F77      = mpif90
CC       = mpicc
SCC      = sw5cc

EXE = fft.exe

CFLAGS   = -O3 -lm -fno-math-errno -msimd -OPT:Ofast
#CFLAGS   = 
SCFLAGS  = -O3 -slave -msimd -lm_slave -fno-math-errno -OPT:Ofast
LF77OPTS = -O3 -lm_slave -lm -fno-math-errno -msimd -OPT:Ofast
#LF77OPTS = 
#FFLAGS   = -O3 -fno-math-errno -msimd -cpp -DTEST -DSWAP -DCHECK -OPT:ieee_arith=1
FFLAGS   = -O3 -msimd -fno-math-errno

##intel, debug
#CFLAGS   = -O2 -lm -fno-math-errno -g -cpp -DTEST -DSWAP -DCHECK -OPT:ieee_arith=1
#LF77OPTS = -O3 -g
#FFLAGS   = -O2 -fno-math-errno -g -DTEST -DSWAP -DCHECK -OPT:ieee_arith=1

#SRC =  global.f ft.f  \
#print_results.f randi8.f       \
#timers.f


CINC = \
slave_kernel.h c_public_var.h c_public_const.h

COBJS = \
c_public_const.o c_public_var.o mpe_fft.o 

FOBJS = \
global.o \
ft.o \
randi8.o \
print_results.o \
timers.o \

SOBJS = slave_kernel.o

VPATH = ./f_source/:./c_header/:./c_master/:./c_slave/

$(EXE)  :$(COBJS) $(SOBJS) #$(OBJS) 
	mpif90 -c $(FFLAGS) -freeform ./f_source/global.f
	mpif90 -c $(FFLAGS) ./f_source/ft.f
	mpif90 -c $(FFLAGS) ./f_source/randi8.f
	mpif90 -c $(FFLAGS) ./f_source/print_results.f
	mpif90 -c $(FFLAGS) ./f_source/timers.f
	$(F77) $(LF77OPTS) -o $@ $(FOBJS) $(COBJS) $(SOBJS)

$(SOBJS) : %.o : %.c
	$(SCC) $(SCFLAGS) -c -o $@ $<

$(COBJS) : %.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXE) *.o work.* *.pc *.mod *.rmod
