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

// ijk Simple triple loop algorithm with simple single register reuse
void ijkSimple (double *A, double *B, double *C1, int N){
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k;
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            register double r = C1[i*N+j];
            for (k=0; k<N; k++){
                r += A[i*N+k] * B[k*N+j];
            }
            C1[i*N+j] = r;
        }
    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

// jik Simple triple loop algorithm with simple single register reuse
void jikSimple (double *A, double *B, double *C2, int N){
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k;
    for (j=0; j<N; j++){
        for (i=0; i<N; i++){
            register double r = C2[i*N+j];
            for (k=0; k<N; k++){
                r += A[i*N+k] * B[k*N+j];
            }
            C2[i*N+j] = r;
        }
    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}



// ikj Simple triple loop algorithm with simple single register reuse
void ikjSimple (double *A, double *B, double *C3, int N){
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k;
    for (i=0; i<N; i++){
        for (k=0; k<N; k++){
            register double r = A[i*N+k];
            for (j=0; j<N; j++){
                C3[i*N+j] += r * B[k*N+j];
            }
        }
    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

// kij Simple triple loop algorithm with simple single register reuse
void kijSimple (double *A, double *B, double *C4, int N){
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k;
    for (k=0; k<N; k++){
        for (i=0; i<N; i++){
            register double r = A[i*N+k];
            for (j=0; j<N; j++){
                C4[i*N+j] += r * B[k*N+j];
            }
        }
    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

// jki Simple triple loop algorithm with simple single register reuse
void jkiSimple (double *A, double *B, double *C5, int N){
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k;
    for (j=0; j<N; j++){
        for (k=0; k<N; k++){
            register double r = B[k*N+j];
            for (i=0; i<N; i++){
                C5[i*N+j] += A[i*N+k] * r;
            }
        }
    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

// kji Simple triple loop algorithm with simple single register reuse
void kjiSimple (double *A, double *B, double *C6, int N){
    clock_gettime(CLOCK_REALTIME, &startTime);
    int i, j, k;
    for (k=0; k<N; k++){
        for (j=0; j<N; j++){
            register double r = B[k*N+j];
            for (i=0; i<N; i++){
                C6[i*N+j] += A[i*N+k] * r;
            }
        }
    }
    clock_gettime(CLOCK_REALTIME,&endTime);
}

void verifyResults(double *A, double *B, double *Matrix1, double *Matrix2, char* examType, int N){

    double maxDiff = 0;
    double tempDiff = 0;

    int i, j;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            tempDiff = Matrix1[i*N+j] - Matrix2[i*N+j];
            if (tempDiff>maxDiff){
                maxDiff = tempDiff;
            }
        }
    }

    printf("Maximum Difference for Experiment %s : %f\n",examType,maxDiff);
}

void verifyResults2(double *A, double *B, double *Matrix1, double *Matrix2, char* examType, int N){

    double Diff = fabs(Matrix1[0]-Matrix2[0]);
    double tempDiff;

    double maxA = fabs(A[0]);
    double maxB = fabs(B[0]);

    int i;
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


void runAndVerify(double *A, double *B, double *C1, double *C2, double *C3, double *C4, double *C5, double *C6, int N){

    printf("MatrixSize: %d\n", N);

    ijkSimple(A, B, C1, N);
    printRunTime("ijkSimple", N);

    jikSimple(A, B, C2, N);
    printRunTime("jikSimple", N);

    ikjSimple(A, B, C3, N);
    printRunTime("ikjSimple", N);

    kijSimple(A, B, C4, N);
    printRunTime("kijSimple", N);

    jkiSimple(A, B, C5, N);
    printRunTime("jkiSimple", N);

    kjiSimple(A, B, C6, N);
    printRunTime("kjiSimple", N);


    printf("\n******************************************\n");

    verifyResults(A, B, C1, C2, "jikSimple", N);
    verifyResults(A, B, C1, C3, "ikjSimple", N);
    verifyResults(A, B, C1, C4, "kijSimple", N);
    verifyResults(A, B, C1, C5, "jkiSimple", N);
    verifyResults(A, B, C1, C6, "kjiSimple", N);

    printf("\n******************************************\n");

    verifyResults2(A, B, C1, C2, "jikSimple", N);
    verifyResults2(A, B, C1, C3, "ikjSimple", N);
    verifyResults2(A, B, C1, C4, "kijSimple", N);
    verifyResults2(A, B, C1, C5, "jkiSimple", N);
    verifyResults2(A, B, C1, C6, "kjiSimple", N);

    printf("\n******************************************\n");


}

void printMatrix (double *Matrix, int N){

    int i,j;
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            printf("%f\t",Matrix[i*N+j]);
        }
        printf("\n");
    }

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
    long int totalElements = N * N;

    double *A, *B, *C1, *C2, *C3, *C4, *C5, *C6;

    uint64_t runTime_Nsec;
    double runTime_Sec;
    struct timespec startTime;
    struct timespec endTime;
    double numberOfOperations;
    double GFLOPS;

    N = atoi(argv[1]); // SET N

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

    runAndVerify(A, B, C1, C2, C3, C4, C5, C6, N);

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

