#!/bin/bash
set -x
PCG_HOME=`pwd`
NETCDF_CPUINC=${PCG_HOME}/../local/include
NETCDF_CPULIB=${PCG_HOME}/../local/lib
cd ./src
rm -f *.o *.mod

FCFLAG="-align array256byte -ipo -mcmodel=medium -qopenmp -I${NETCDF_CPUINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz -free -c"
CCFLAG="-align array256byte -ipo -mcmodel=medium -qopenmp -I${NETCDF_CPUINC} -I. -O3 -xCORE-AVX2 -fp-model source -ftz -c"
LDFLAG="-align array256byte -ipo -O3 -mcmodel=medium -qopenmp -L${NETCDF_CPULIB} -lnetcdf"

mpiifort $FCFLAG netcdf_mod.f90
mpiifort $FCFLAG kinds_mod.f90
mpiifort $FCFLAG domain_size.f90
mpiifort $FCFLAG constants.f90
mpiifort $FCFLAG communicate.f90
mpiifort $FCFLAG exit_mod.f90
mpiifort $FCFLAG blocks.f90
mpiifort $FCFLAG broadcast.f90
mpiifort $FCFLAG distribution.f90
mpiifort $FCFLAG io_types.f90
mpiifort $FCFLAG boundary.f90
mpiifort $FCFLAG domain.f90
mpiifort $FCFLAG gather_scatter.f90
mpiifort $FCFLAG io_netcdf.f90
mpiifort $FCFLAG io_binary.f90
mpiifort $FCFLAG io.f90
mpiifort $FCFLAG global_reductions.f90
mpiifort $FCFLAG grid.f90
mpiifort $FCFLAG prognostic.f90
mpiifort $FCFLAG time_management.f90
mpiifort $FCFLAG solver_pcg_mod.f90
mpiifort $FCFLAG timers.f90
mpiifort $FCFLAG initial.f90
mpiifort $FCFLAG pcg.f90
mpiicc   $CCFLAG hch_timer.c
mpiicc   $CCFLAG doge.c

# pop_9-off

mpiifort -o ../bin/pcg.x doge.o hch_timer.o blocks.o netcdf_mod.o constants.o distribution.o domain.o exit_mod.o grid.o initial.o io_binary.o io.o io_netcdf.o io_types.o kinds_mod.o pcg.o prognostic.o solver_pcg_mod.o domain_size.o boundary.o  broadcast.o communicate.o gather_scatter.o global_reductions.o time_management.o timers.o $LDFLAG 
