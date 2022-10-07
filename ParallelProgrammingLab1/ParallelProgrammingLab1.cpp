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

 
#include <omp.h>
#include <cstdio>
#include <iostream>


//#define Critical
//#define Atomic
#define Reduction


#define NMAX 6000000
#define N 3
#define Q 18


double reduction(const float* a, double sum, int i) {
#pragma omp parallel  for shared(a)  private(i) reduction(+: sum)
    for (i = 0; i < NMAX; ++i) { sum = sum + a[i]; }
    return sum;
}
double critical(const float* a, double sum, int i) {
#pragma omp parallel  for
    for (i = 0; i < NMAX; ++i) {
#pragma omp critical
        sum = sum + a[i];
    }
    return sum;
}
double atomic(const float* a, double sum, int i) {
#pragma omp  parallel  for  shared(a)
    for (i = 0; i < NMAX; ++i) {
#pragma omp atomic
        sum += a[i];
    }
    return sum;
}
void OpenMP(int, char* []);
void OpenMpSequientel(int, char* []);

extern const int NORMAL_SPREAD = 12;

int main(int argc, char* argv[]) {
    //OpenMpSequientel(argc, argv);
    OpenMP(argc, argv);
    return 0;
}
void OpenMpSequientel(int argc, char* argv[]) {
    double st_time = omp_get_wtime(), end_time;
    for (int i = 0; i < NORMAL_SPREAD; ++i) {
        float  sum = 0.0;
        auto* a = new float[NMAX];
        std::fill(a, a + NMAX, 1.0);
        for (int j = 0; j < Q; ++j) {
            for (int k = 0; k < NMAX; ++k) sum += a[k];
        }
        sum /= Q;
        printf("\nSum = %10.2f", sum);
    }
    end_time = omp_get_wtime();
    end_time = (end_time - st_time) / NORMAL_SPREAD;
    printf("\nTime = %f \n", end_time);
}
void OpenMP(int argc, char* argv[]) {
    int num = N;
    omp_set_num_threads(num);
    printf("\nProcess num = %d", num);
    int i = 0;
    float  sum = 0.0;
    auto* a = new float[NMAX];
    std::fill(a, a + NMAX, 1);
    double st_time, end_time;
    st_time = omp_get_wtime();

#ifdef Critical
    printf("\nType = Critical");
    for (int j = 0; j < Q; ++j) { sum = critical(a, sum, i); }
    sum /= Q;
#endif
#ifdef Atomic
    printf("\nType = Atomic");
    for (int j = 0; j < Q; ++j) { sum = atomic(a, sum, i); }
    sum /= Q;
#endif
#ifdef Reduction
    printf("\nType = Reduction");
    for (int j = 0; j < Q; ++j) { sum = reduction(a, sum, i); }
    sum /= Q;
#endif

    printf("\nSum = %10.2f", sum);
    end_time = omp_get_wtime();
    end_time = end_time - st_time;
    printf("\nTime = %f \n", end_time);
    delete[] a;
}




// MPI

//#include "mpi.h"
//#include <cstdio>
//#include <algorithm>
//#include <iostream>
//#include <fstream>
//
////#define Point_Point
//#define Collect
//#define Q 18
//#define NMAX 6000000
//
//
//int main(int argc, char* argv[])
//{
//    double TotalSum, ProcSum = 0.0;
//    int ProcRank, ProcNum;
//    float* a = nullptr;
//    MPI_Status Status;
//    double st_time, end_time;
//
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
//    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//    if (ProcRank == 0) {
//        printf("\nPorcess number = %d", ProcNum);
//        a = new float[NMAX];
//        std::fill(a, a + NMAX, 1.0);
//    }
//    st_time = MPI_Wtime();
//
//    int k = NMAX / ProcNum;
//    auto* localArray = new float[k];
//    MPI_Scatter(a, k, MPI_FLOAT, localArray, k, MPI_FLOAT, 0, MPI_COMM_WORLD);
//    if (ProcRank == 0) delete[] a;
//
//    for (int q = 0; q < Q; ++q) {
//        for (int i = 0; i < k; ++i)
//            ProcSum += localArray[i];
//    }
//    delete[] localArray;
//#ifdef Point_Point
//    if (ProcRank == 0) {
//        printf("\nType = Point-Point");
//        TotalSum = ProcSum;
//        for (int i = 1; i < ProcNum; ++i) {
//            MPI_Recv(&ProcSum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
//            TotalSum = TotalSum + ProcSum;
//        }
//        TotalSum /= Q;
//    }
//    else
//        MPI_Send(&ProcSum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//#endif
//#ifdef Collect
//    if (ProcRank == 0) printf("\nType = Collect");
//    MPI_Reduce(&ProcSum, &TotalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    TotalSum /= Q;
//#endif
//    MPI_Barrier(MPI_COMM_WORLD);
//    end_time = MPI_Wtime();
//    end_time = end_time - st_time;
//
//    if (ProcRank == 0) {
//        printf("\nSum = %10.2f", TotalSum);
//        printf("\nTime = %f ", end_time);
//    }
//    MPI_Finalize();
//}
