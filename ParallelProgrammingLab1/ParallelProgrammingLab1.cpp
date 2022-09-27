// LR1-----------------------------------------------------------------------------------------------------------------------------

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

//#include <iostream>
//#include <time.h>
//
//int main(int argc, char* argv[]) {
//	clock_t tStart = clock();
//	for (int i = 0; i < 500000; ++i) {}
//	printf("Hello world\n");
//	for (int i = 0; i < 500000; ++i) {}
//	printf("Hello world\n");
//	for (int i = 0; i < 500000; ++i) {}
//	printf("Hello world\n");
//	for (int i = 0; i < 500000; ++i) {}
//	printf("Hello world\n");
//	for (int i = 0; i < 500000; ++i) {}
//	printf("Hello world\n");
//	for (int i = 0; i < 500000; ++i) {}
//	printf("Hello world\n");
//	for (int i = 0; i < 500000; ++i) {}
//	printf("Hello world\n");
//	for (int i = 0; i < 500000; ++i) {}
//	printf("Hello world\n");
//	printf("Time taken: %f\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
//}


//LR2---------------------------------------------------------------------------------------------------------------------------------




// OpenMP

 
//#include <omp.h>
//#include "stdio.h"
//#include <iostream>
//#include <ctime>
//#define CHUNK 100
//#define NMAX 10000000
//
//double reduction(double* a, double sum, int i) {
//#pragma omp parallel for shared(a) private(i) reduction(+: sum) 
//    for (i = 0; i < NMAX; i++) {
//        sum = sum + a[i];
//    }
//    return sum;
//}
//
//double critical(double* a, double sum, int i) {
//#pragma omp parallel for 
//    for (i = 0; i < NMAX; i++) {
//#pragma omp critical
//        {
//            sum = sum + a[i];
//        }
//    }
//    return sum;
//}
//
//double atomic(double* a, double sum, int i) {
//#pragma omp parallel for shared(a) 
//    for (i = 0; i < NMAX; i++) {
//#pragma omp atomic
//        sum = sum + a[i];
//    }
//    return sum;
//}
//
//
//int main(int argc, char* argv[]) {
//    omp_set_num_threads(1);
//    int i;
//    double  sum;
//    srand(time(0));
//    double* a = (double*)malloc(sizeof(double) * NMAX);
//    for (i = 0; i < NMAX; i++)
//    {
//        a[i] = 1;
//    }
//    double st_time, end_time;
//    st_time = omp_get_wtime();
//    sum = 0;
//
//    sum = reduction(a, sum, i);
//    //sum = critical(a, sum, i);
//    //sum = atomic(a, sum, i);
//
//    end_time = omp_get_wtime();
//    end_time = end_time - st_time;
//    printf("\nTotal Sum = %10.2f", sum);
//    printf("\nTIME OF WORK IS %f \n", end_time);
//    free(a);
//    return 0;
//}






// MPI

#include <stdlib.h>
#include "mpi.h"
#include <stdio.h>
#include <math.h>

#define NMAX 1000000
int main(int argc, char* argv[])
{
	double x[NMAX], TotalSum, ProcSum = 0.0;
	int ProcRank, ProcNum, N = 1000000, i;
	MPI_Status Status;

	double st_time, end_time;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (ProcRank == 0)
	{
		for (i = 0; i < NMAX; i++)
			    {
			        x[i] = 1;
			    }
	}
	st_time = MPI_Wtime();

	int k = N / ProcNum;
	int i1 = k * ProcRank;
	int i2 = k * (ProcRank + 1);

	if (ProcRank == ProcNum - 1) i2 = N;
	for (i = i1; i < i2; i++)
		ProcSum = ProcSum + x[i];

	if (ProcRank == 0)
	{
		TotalSum = ProcSum;
		for (i = 1; i < ProcNum; i++)
		{
			MPI_Recv(&ProcSum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
			TotalSum = TotalSum + ProcSum;
		}
	}
	else
		MPI_Send(&ProcSum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	end_time = MPI_Wtime();
	end_time = end_time - st_time;

	if (ProcRank == 0)
	{
		printf("\nTotal Sum = %10.2f", TotalSum);
		printf("\nTIME OF WORK IS %f ", end_time);
	}

	MPI_Finalize();
	return 0;
}