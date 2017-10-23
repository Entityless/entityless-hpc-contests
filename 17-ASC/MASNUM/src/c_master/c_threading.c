#ifndef __C_THREAD_FOR_FORTRAN__
#define __C_THREAD_FOR_FORTRAN__

#include <pthread.h>
#include <stdio.h>
#include "../c_header/c_acce_adapter.h"

#define num_of_thread 1

void inner_propagat_begin_(int *command) {
    int v = *command;
    c_propagat_inner(v);
}

void inner_propagat_stop_(int* command) {
	//int i;
	//for (i = 0; i < num_of_thread; ++ i)
	//	pthread_join(thread[i], NULL);

   int v = *command;
    c_propagat_inner_stop(v);
}

#undef num_of_thread

#endif
