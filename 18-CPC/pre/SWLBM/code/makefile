TARGET = LbmCavity3D
USER = $(shell whoami)
CC = sw5cc
LD = mpicc 

CFLAGS =  -O3 -host -I/usr/sw-mpp/mpi2/include/ -lm -OPT:ieee_arith=3 -keep
#CFLAGS   = -O3 -lm -fno-math-errno -msimd -DTEST -DSWAP -DCHECK -OPT:ieee_arith=3
#-ver 5.421-sw-gy -LNO:simd=1:simd_report=2
SCFLAGS  =  -O3 -slave -msimd -lm_slave -fno-math-errno -OPT:ieee_arith=3  -keep 
LDFLAGS  = -O3 -msimd -lm_slave -fno-math-errno -OPT:ieee_arith=3 -keep

COBJ = LbmCavity3D.o Collide.o Parallel.o Stream.o c_public_var.o hch_timer.o mpe_lbm.o
SOBJ = slave_kernel.o

LIB = lib/liblbm.a

VPATH = ./c_header/:./c_master/:./c_slave/:./lbm_original/

$(TARGET): $(COBJ) $(SOBJ)
	$(LD) $(LDFLAGS) $(COBJ) $(SOBJ) $(LIB) -o $(TARGET) 
	# rm $(COBJ) $(SOBJ)

$(SOBJ) : %.o : %.c
	$(SCC) $(SCFLAGS) -c -o $@ $<

$(COBJ) : %.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

run:
	rm -rf doge.log
	bsub -I -b -q q_sw_expr -cgsp 64 -n 16 -np 4 -share_size 6500 -host_stack 500 -J test -o doge.log ./LbmCavity3D $(USER)

#-------------------------------------*
.PHONY : clean clear
clean:
	-rm -rf $(TARGET) $(COBJ) $(SOBJ)
	
CC       = mpicc
SCC      = sw5cc
