#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

uint64_t runTime_Nsec;
double runTime_Sec;
struct timespec startTime;
struct timespec endTime;

double GFLOPS;
double numberOfOperations;


// ijk - Blocked version algorithm
void ijkBlocked (double *A, double *B, double *C1, int N, int BB)
{
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k, i1, j1, k1;
    for (i = 0; i < N; i += BB)
        for (j = 0; j < N; j += BB)
            for (k = 0; k < N; k += BB)
                for (i1 = i; i1 < i+BB; i1++)
                    for (j1 = j; j1 < j+BB; j1++)
                    {
                        register double r = C1[i1*N+j1];
                        for (k1 = k; k1 < k+BB; k1++)
                            r += A[i1*N+k1] * B[k1*N+j1];
                        C1[i1*N+j1] = r;
                    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

// jik - Blocked version algorithm
void jikBlocked (double *A, double *B, double *C2, int N, int BB)
{
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k, i1, j1, k1;
    for (j = 0; j < N; j += BB)
        for (i = 0; i < N; i += BB)
            for (k = 0; k < N; k += BB)
                for (j1 = j; j1 < j+BB; j1++)
                    for (i1 = i; i1 < i+BB; i1++)
                    {
                        register double r = C2[i1*N+j1];
                        for (k1 = k; k1 < k+BB; k1++)
                            r += A[i1*N+k1] * B[k1*N+j1];
                        C2[i1*N+j1] = r;
                    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

// ikj - Blocked version algorithm
void ikjBlocked (double *A, double *B, double *C3, int N, int BB)
{
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k, i1, j1, k1;
    for (i = 0; i < N; i += BB)
        for (k = 0; k < N; k += BB)
            for (j = 0; j < N; j += BB)
                for (i1 = i; i1 < i+BB; i1++)
                    for (k1 = k; k1 < k+BB; k1++)
                    {
                        register double r = A[i1*N+k1];
                        for (j1 = j; j1 < j+BB; j1++)
                            C3[i1*N+j1] += r * B[k1*N+j1];
                    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

// kij - Blocked version algorithm
void kijBlocked (double *A, double *B, double *C4, int N, int BB)
{
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k, i1, j1, k1;
    for (k = 0; k < N; k += BB)
        for (i = 0; i < N; i += BB)
            for (j = 0; j < N; j += BB)
                for (k1 = k; k1 < k+BB; k1++)
                    for (i1 = i; i1 < i+BB; i1++)
                    {
                        register double r = A[i1*N+k1];
                        for (j1 = j; j1 < j+BB; j1++)
                            C4[i1*N+j1] += r * B[k1*N+j1];
                    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

// jki - Blocked version algorithm
void jkiBlocked (double *A, double *B, double *C5, int N, int BB)
{
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k, i1, j1, k1;
    for (j = 0; j < N; j += BB)
        for (k = 0; k < N; k += BB)
            for (i = 0; i < N; i += BB)
                for (j1 = j; j1 < j+BB; j1++)
                    for (k1 = k; k1 < k+BB; k1++)
                    {
                        register double r = B[k1*N+j1];
                        for (i1 = i; i1 < i+BB; i1++)
                            C5[i1*N+j1] += A[i1*N+k1] * r;
                    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

// kji - Blocked version algorithm
void kjiBlocked (double *A, double *B, double *C6, int N, int BB)
{
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k, i1, j1, k1;
    for (k = 0; k < N; k += BB)
        for (j = 0; j < N; j += BB)
            for (i = 0; i < N; i += BB)
                for (k1 = k; k1 < k+BB; k1++)
                    for (j1 = j; j1 < j+BB; j1++)
                    {
                        register double r = B[k1*N+j1];
                        for (i1 = i; i1 < i+BB; i1++)
                            C6[i1*N+j1] += A[i1*N+k1] * r;
                    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

void verifyResults(double *A, double *B, double *Matrix1, double *Matrix2, char* examType, int N){

    double Diff = fabs(Matrix1[0]-Matrix2[0]);
    double tempDiff;

    double maxA = fabs(A[0]);
    double maxB = fabs(B[0]);

    long int i;
    for (i=0; i<N*N; i++){
        tempDiff = fabs(Matrix1[i]-Matrix2[i]);
        if (tempDiff>Diff){
            Diff = tempDiff;
        }
        if (fabs(A[i])>maxA){
            maxA = fabs(A[i]);
        }
        if (fabs(B[i])>maxB){
            maxB = fabs(B[i]);
        }
    }

    double maxDiff = Diff / (maxA * maxB);
    printf("Maximum Difference 2 for Experiment %s : %f\n",examType ,maxDiff);
}


void calculateElapsedTimeAndGFLOPs(int N){
    runTime_Nsec = 1000000000L * (endTime.tv_sec-startTime.tv_sec) + endTime.tv_nsec - startTime.tv_nsec;
    runTime_Sec = runTime_Nsec * pow(10,(-9));
    numberOfOperations = 2 * pow(N,3);
    GFLOPS = numberOfOperations / runTime_Nsec;
}


void printRunTime(char* testName, int N){
    calculateElapsedTimeAndGFLOPs(N);
    printf("TEST %s: runTime = %f s \t %llu ns\tnumber of FloatingPoint Operations = %.0f\tGFLOPS = %f\n"
            ,testName
            ,runTime_Sec
            ,(long long unsigned int)runTime_Nsec
            ,numberOfOperations
            ,GFLOPS);
}


void runAndVerify(double *A, double *B, double *C1, double *C2, double *C3, double *C4, double *C5, double *C6, int N, int BB){


    printf("BlockSize: %d\n", BB);
    printf("MatrixSize: %d\n", N);

    printf("\n******************************************\n");

    ijkBlocked(A, B, C1, N, BB);
    printRunTime("ijkBlocked", N);

    jikBlocked(A, B, C2, N, BB);
    printRunTime("jikBlocked", N);

    ikjBlocked(A, B, C3, N, BB);
    printRunTime("ikjBlocked", N);

    kijBlocked(A, B, C4, N, BB);
    printRunTime("kijBlocked", N);

    jkiBlocked(A, B, C5, N, BB);
    printRunTime("jkiBlocked", N);

    kjiBlocked(A, B, C6, N, BB);
    printRunTime("kjiBlocked", N);

    printf("\n******************************************\n");

    verifyResults(A, B, C1, C2, "jikBlocked", N);
    verifyResults(A, B, C1, C3, "ikjBlocked", N);
    verifyResults(A, B, C1, C4, "kijBlocked", N);
    verifyResults(A, B, C1, C5, "jkiBlocked", N);
    verifyResults(A, B, C1, C6, "kjiBlocked", N);

    printf("\n******************************************\n");


}

void initializeMatrices(double *A, double *B, double *C1, double *C2, double *C3, double *C4, double *C5, double *C6, int N) {

    srand(time(NULL));

    int i, j;

    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
                A[i*N+j] = (double)rand()/(double) RAND_MAX;
                B[i*N+j] = (double)rand()/(double) RAND_MAX;
                C1[i*N+j] = 0;
                C2[i*N+j] = 0;
                C3[i*N+j] = 0;
                C4[i*N+j] = 0;
                C5[i*N+j] = 0;
                C6[i*N+j] = 0;
            }
        }
    }

int main (int argc, char* argv[])
{
    int N = 64;
    int BB = 2;

    long int totalElements = N * N;

    double *A, *B, *C1, *C2, *C3, *C4, *C5, *C6;

    uint64_t runTime_Nsec;
    double runTime_Sec;
    struct timespec startTime;
    struct timespec endTime;
    double numberOfOperations;
    double GFLOPS;

    N = atoi(argv[1]); // SET N
    BB = atoi(argv[2]);

    totalElements = N * N;

    A = (double*)calloc(totalElements,sizeof(double));
    B = (double*)calloc(totalElements,sizeof(double));
    C1 = (double*)calloc(totalElements,sizeof(double));
    C2 = (double*)calloc(totalElements,sizeof(double));
    C3 = (double*)calloc(totalElements,sizeof(double));
    C4 = (double*)calloc(totalElements,sizeof(double));
    C5 = (double*)calloc(totalElements,sizeof(double));
    C6 = (double*)calloc(totalElements,sizeof(double));

    initializeMatrices(A, B, C1, C2, C3, C4, C5, C6, N);

    runAndVerify(A, B, C1, C2, C3, C4, C5, C6, N, BB);

    free (A);
    free (B);
    free (C1);
    free (C2);
    free (C3);
    free (C4);
    free (C5);
    free (C6);

    return 0;
}

