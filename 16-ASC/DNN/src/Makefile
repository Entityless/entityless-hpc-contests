#TK_INC_PATH	?= /opt/intel/mkl/include
#TK_LIB_PATH	?= /opt/intel/mkl/lib/intel64
#TK_LIB2_PATH?= /opt/intel/lib/intel64

CC	= mpic++

#flags
CCFLAGS	:= -Werror -mkl -xHost -openmp -O2 -fPIC 
CCFLAGSS	:= -Werror -mkl -xHost -openmp -O1 -fPIC 
#-offload-attribute-target=mic

INCLUDES:= -I$(TK_INC_PATH) -I. 

all: build

build: dnntk

dnn_helper.o: dnn_helper.cpp
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ -c $<

dnn_kernel.o: dnn_kernel.cpp
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ -c $<

dnn_utility.o: dnn_utility.cpp
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ -c $<

kernel_util.o: kernel_util.cpp
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ -c $<

mic_kernel_util.o: mic_kernel_util.cpp
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ -c $<

sigmoid_util.o: sigmoid_util.cpp
	$(CC) $(CCFLAGSS) $(INCLUDES)  -o $@ -c $<

my_dnn_util.o: my_dnn_util.cpp
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ -c $< 
    
main.o: main.cpp
	$(CC) $(CCFLAGS) $(INCLUDES)  -o $@ -c $<
    
dnntk: dnn_helper.o dnn_kernel.o dnn_utility.o main.o my_dnn_util.o kernel_util.o mic_kernel_util.o sigmoid_util.o
	$(CC) $(CCFLAGS)  -o $@ $+ $(LDFLAGS)
	cp -f dnntk ../expsakamoto/       
	cp -f dnntk ../expsakamoto_2/       
	cp -f dnntk ../expsakamoto_3/       
	mv -f dnntk ../expsakamoto_4/       
clean:
	rm -f ../expsakamoto/dnntk ../expsakamoto_2/dnntk ../expsakamoto_3/dnntk ../expsakamoto_4/dnntk  *.o *~ | clear
