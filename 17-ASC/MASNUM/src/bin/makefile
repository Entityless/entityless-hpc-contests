F77      = mpif90
CC       = mpicc
SCC      = sw5cc
CFLAGS   = -O2 -lm  -fno-math-errno -msimd -D ASCEXP12
SCFLAGS  = -O2 -slave -msimd -lm_slave  -fno-math-errno
LF77OPTS = -O2 -lm_slave -lm  -fno-math-errno -msimd
FFLAGS   = -O2 -fno-math-errno -msimd

EXE = masnum.wam.mpi

NETCDF_PATH = /usr/sw-mpp/apps/lib/netcdf

LIBS = -L$(NETCDF_PATH)/lib -lnetcdf -L/usr/local/lib 

INC_PATH = -I$(NETCDF_PATH)/include

SRC = \
netcdf_mod.f90 time_mod.f90  wamvar_mod.f90 wamcpl_mod.f90 \
wamfio_mod.f90 wamcor_mod.f90  ympi_mod.f90 wammpi_mod.f90 \
masnum_wam_mpi.f90

INC = \
get_wind.inc  mixture.inc                   wamfio_mod_init.inc   \
mixture_limit.inc             mixture_wit.inc                     \
implsch.inc   nlweight.inc    propagat.inc                        \
intact.inc    outmix.inc                                          \
inter.inc     output.inc      setspec.inc                         \
io_rest.inc   outwav.inc      settopog.inc                        \
mean1.inc     out_wavmix.inc  set_uv.inc                          \
mean2.inc     outwav_t.inc    setwave.inc  outmix_wit.inc

QSRC=

CINC = \
para_marco.h

COBJS = \
c_threading.o c_acce_adapter.o cpu_implsch.o cpu_propagat.o cpu_mean1.o cpu_getem.o c_public_var.o cpu_setspec.o

SOBJS = \
slave_kernel.o

VPATH= ../wave_cor/:../ympi_mod/:../c_header/:../c_master/:../c_slave/

OBJS = ${SRC:.f90=.o}
QOBJS= ${QSRC:.f90=.o}

.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	$(F77) $(INC_PATH) $(FFLAGS) -c $<


$(EXE)  :$(INC) $(OBJS)  $(QOBJS) $(COBJS) $(SOBJS)
	$(F77) $(LF77OPTS) -o $@ $(OBJS) $ $(QOBJS) $(COBJS) $(SOBJS) $(LIBS) $(F90MIXEDCXX)

$(OBJS) $(QOBJS) : $(INC)

$(SOBJS) : %.o : %.c
	$(SCC) $(SCFLAGS) -c -o $@ $<

$(COBJS) : %.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.h:

clean:
	rm -f  $(EXE) *.o *.mod 

cleano:
	rm -f  *.o *.mod

