//#include "mpi.h"
//#include "stdio.h"
//#include <time.h>
//#include <math.h>
//#include <omp.h>
//
//int main(int argc, char* argv[])
//{
//	int rank, ranksize, i;
//	MPI_Init(&argc, &argv);// 
//	//Определяем свой номер в группе:
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //
//	//Определяем размер группы:
//	MPI_Comm_size(MPI_COMM_WORLD, &ranksize);//
//	double start = MPI_Wtime();
//	for (int i = 0; i < 500000; ++i) {}
//	printf("Hello world!\tprocess %d of %d\t\t", rank, ranksize);
//	double end_time = MPI_Wtime() - start;
//	printf("Time of working process %d = %f.\n", rank, end_time);
//	MPI_Finalize();//
//	return 0;
//}

//#include <math.h>
//#include <omp.h>
//#include <time.h>
//#include <stdlib.h>
//#include <locale.h>
//#include <stdio.h>
//
//int main(int argc, char* argv[])
//{
//	omp_set_num_threads(8);
//	int nTheads, theadNum;
//#pragma omp parallel  private(nTheads, theadNum)
//	{
//		nTheads = omp_get_num_threads();
//		theadNum = omp_get_thread_num();
//		double start = omp_get_wtime();
//		for (int i = 0; i < 500000; ++i) {}
//		printf("OpenMP thread %d from %d threads \t\t", theadNum, nTheads);
//		double end = omp_get_wtime() - start;
//		printf("Time of working process %d = %f \n", theadNum, end);
//	}
//	return 0;
//}

#include <iostream>
#include <time.h>

int main(int argc, char* argv[]) {
	clock_t tStart = clock();
	for (int i = 0; i < 500000; ++i) {}
	printf("Hello world\n");
	for (int i = 0; i < 500000; ++i) {}
	printf("Hello world\n");
	for (int i = 0; i < 500000; ++i) {}
	printf("Hello world\n");
	for (int i = 0; i < 500000; ++i) {}
	printf("Hello world\n");
	for (int i = 0; i < 500000; ++i) {}
	printf("Hello world\n");
	for (int i = 0; i < 500000; ++i) {}
	printf("Hello world\n");
	for (int i = 0; i < 500000; ++i) {}
	printf("Hello world\n");
	for (int i = 0; i < 500000; ++i) {}
	printf("Hello world\n");
	printf("Time taken: %f\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
}
