EXE=ocean
FC=ifort
FFLAGS= -fpp -align array64byte -qopt-assume-safe-padding -static-intel -c -qopt-report -xHost -fp-model fast=2 -fno-inline -qopenmp -O3 -ipo -qopt-streaming-cache-evict=0 -qopt-streaming-stores always -DSTREAM_ARRAY_SIZE=8000000
LDFLAGS= -fpp -rdynamic -align array64byte -qopt-assume-safe-padding -qopt-report -static-intel -ltbbmalloc_proxy -xHost -fp-model fast=2 -fno-inline -qopenmp -O3 -ipo -qopt-streaming-cache-evict=0 -qopt-streaming-stores always -DSTREAM_ARRAY_SIZE=8000000
SRC=mod_data.f90 mpdata_adiff.f90 ocean.f90 

OBJ=$(patsubst %.f90,%.o,$(SRC))

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $^ libverify.a -o $@

%.o: %.f90
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm *.o *.mod

