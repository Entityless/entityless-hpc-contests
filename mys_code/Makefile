EXE	= yooooooo
CC	= mpicc
SCC	= sw5cc
HFLAGS	= -host  -std=gnu99 -O3 -msimd -OPT:Ofast
SFLAGS	= -slave -std=gnu99 -O3 -msimd -OPT:Ofast -lmfast_slave
LDFLAGS	= -std=gnu99 -O3 -msimd  -OPT:Ofast

$(EXE):kirishima_host.o touka_slave.o
	$(CC) $(LDFLAGS) *.o -o $@

kirishima_host.o: ./src/kirishima_host.c
	$(CC) -c $(HFLAGS) ./src/kirishima_host.c

touka_slave.o: ./src/touka_slave.c
	$(SCC) -c $(SFLAGS) ./src/touka_slave.c

clean:
	rm *.o
	rm $(EXE)
