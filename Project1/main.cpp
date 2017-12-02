#include <iostream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>



using namespace std;

class drsvrTimer {

    private:

        uint64_t runTime_Nsec;
        double runTime_Sec;

        struct timespec startTime;
        struct timespec endTime;


        bool timerActive;
        int N;

        double numberOfOperations;
        double GFLOPS;

        void calculateElapsedTime(){
            runTime_Nsec = 1000000000L * (endTime.tv_sec-startTime.tv_sec) + endTime.tv_nsec - startTime.tv_nsec;
            runTime_Sec = runTime_Nsec * pow(10,(-9));
        }

        void calculateGFLOPs(string testName){
            if (testName.compare("dgemm0")){
                numberOfOperations = 2 * pow(N,3);
            }
            else if (testName.compare("dgemm1")){
                numberOfOperations = 2 * pow(N,3);
            }
            else if (testName.compare("dgemm2")){
                numberOfOperations = 16 * pow(N/2,3);
            }
            else if (testName.compare("dgemm3")){
                numberOfOperations = 54 * pow(N/3,3);
            }
            else {
                numberOfOperations = 2 * pow(N,3);
                printf("Test Name is not defined!\n");
            }

            GFLOPS = numberOfOperations / runTime_Nsec;
        }

    public:

        drsvrTimer(int n){
            runTime_Nsec = 0;
            runTime_Sec = 0;
            GFLOPS = 0;
            timerActive = false;
            N = n;
        }

        void printRunTime(string testName){
            this->calculateElapsedTime();
            this->calculateGFLOPs(testName);
            cout<<testName<<": ";
            printf("runTime = %f s \t %lluns\tnumber of FloatingPoint Operations = %.0f\tGFLOPS = %f\n"
                    ,runTime_Sec
                    ,(long long unsigned int)runTime_Nsec
                    ,numberOfOperations
                    ,GFLOPS);
        }

        void startTimer(){
            clock_gettime(CLOCK_REALTIME,&startTime);
            timerActive = true;
        }

        void stopTimer(){
            if (timerActive){
                clock_gettime(CLOCK_REALTIME,&endTime);
                timerActive = false;
            }
            else {
                printf("TIMER IS NOT RUNNING, I CANNOT STOP IT!");
            }
        }

};


class matrixTestbench {

    private:
        double *A, *B, *C1, *C2, *C3, *C4;   // Pointer to Matrix
        int N;
        int extendedN;
        int totalElements;   // Matrix Dimension
        int boundExtension;  // Matrix Extension filling with zeros

        void initializeMatrices() {

            srand(time(NULL));

            int tempN = N+boundExtension;
            extendedN = tempN;

            for (int i=0; i<tempN; i++){
                for (int j=0; j<tempN; j++){

                    if (i>=N || j>=N){

                        A[i*tempN+j] = 0;
                        B[i*tempN+j] = 0;
                        //printf("if\ti: %d\tj: %d \tA[%d,%d] = %f\n",i,j,i,j,A[i*tempN+j]);
                    }
                    else{
                        A[i*tempN+j] = (double)rand()/(double) RAND_MAX;
                        B[i*tempN+j] = (double)rand()/(double) RAND_MAX;
                        //printf("else\ti: %d\tj: %d \tA[%d,%d] = %f\n",i,j,i,j,A[i*tempN+j]);
                    }

                    C1[i*tempN+j] = 0;
                    C2[i*tempN+j] = 0;
                    C3[i*tempN+j] = 0;
                    C4[i*tempN+j] = 0;
                }
            }
        }

        void printMatrix (double *Matrix){

            for (int i=0; i<N; i++){
                for (int j=0; j<N; j++){
                    printf("%f\t",Matrix[i*N+j]);
                }
                printf("\n");
            }

            printf("\n******************************************\n");
        }

        void printMatrix2 (double *Matrix){
            printf("extendedN:%d: \n", extendedN);
            for (int i=0; i<extendedN; i++){
                for (int j=0; j<extendedN; j++){
                    printf("%f\t",Matrix[i*extendedN+j]);
                }
                printf("\n");
            }

            printf("\n******************************************\n");
        }

        // dgemm0
        // simple ijk version triple loop algorithm
        void dgemm0 (){
            for (int i=0; i<N; i++){
                for (int j=0; j<N; j++){
                    for (int k=0; k<N; k++){
                        C1[i*N+j] += A[i*N+k] * B[k*N+j];
                        //printf("A[%d,%d] = %f\tB[%d,%d] = %f\tC[%d,%d] =%f\n",i,k,A[i*N+k],k,j,B[k*N+j],i,j,C1[i*N+j]);
                    }
                }
            }

        }

        void dgemm1 (){
            for (int i=0; i<N; i++){
                for (int j=0; j<N; j++){
                    register double r = C2[i*N+j];
                    for (int k=0; k<N; k++){
                        r += A[i*N+k] * B[k*N+j];
                    }
                    C2[i*N+j] = r;
                }
            }

        }

        void dgemm2 (){

            for (int i = 0; i < N; i += 2) {
                for (int j = 0; j < N; j += 2) {

                    register int t = i*N+j;
                    register int tt = t+N;

                    register double rc00 = C3[t];
                    register double rc01 = C3[t+1];
                    register double rc10 = C3[tt];
                    register double rc11 = C3[tt+1];

                    for (int k = 0; k < N; k += 2) {

                        register int ta = i*N+k;
                        register int tta = ta+N;
                        register int tb = k*N+j;
                        register int ttb = tb+N;

                        register double ra00 = A[ta];
                        register double ra01 = A[ta+1];
                        register double ra10 = A[tta];
                        register double ra11 = A[tta+1];

                        register double rb00 = B[tb];
                        register double rb01 = B[tb+1];
                        register double rb10 = B[ttb];
                        register double rb11 = B[ttb+1];

                        rc00 += (ra00 * rb00) + (ra01 * rb10);
                        rc01 += (ra00 * rb01) + (ra01 * rb11);
                        rc10 += (ra10 * rb00) + (ra11 * rb10);
                        rc11 += (ra10 * rb01) + (ra11 * rb11);
                    }

                    C3[t] = rc00;
                    C3[t+1] = rc01;
                    C3[tt] = rc10;
                    C3[tt+1] = rc11;

                }
            }
        }

        void dgemm3 (){

            for (int i = 0; i < N; i += 3) {
                for (int j = 0; j < N; j += 3) {

                    register int t = i*N+j; // COLUMN 0 ROW 0
                    register int tt = t+N;  // COLUMN 0 ROW 1
                    register int ttt = tt+N; // COLUMN 0 ROW 2

                    register double rc00 = C4[t];
                    register double rc01 = C4[t+1];
                    register double rc02 = C4[t+2];

                    register double rc10 = C4[tt];
                    register double rc11 = C4[tt+1];
                    register double rc12 = C4[tt+2];

                    register double rc20 = C4[ttt];
                    register double rc21 = C4[ttt+1];
                    register double rc22 = C4[ttt+2];

                    for (int k = 0; k < N; k += 3) {

                        register int ta = i*N+k;
                        register int tta = ta+N;
                        register int ttta = tta+N;

                        register int tb = k*N+j;
                        register int ttb = tb+N;
                        register int tttb = ttb+N;


                        /*
                         *  C00 = A00*B00 + A01*B10 + A02*B20
                         *  C01 = A00*B01 + A01*B11 + A02*B21
                         *  C02 = A00*B02 + A01*B12 + A02*B22
                         *
                         *  C10 = A10*B00 + A11*B10 + A12*B20;
                         *  C11 = A10*B01 + A11*B11 + A12*B21;
                         *  C12 = A10*B02 + A11*B12 + A12*B22;
                         *
                         *  C20 = A20*B00 + A21*B10 + A22*B20;
                         *  C21 = A20*B01 + A21*B11 + A22*B21;
                         *  C22 = A20*B02 + A21*B12 + A22*B22;
                         */



                        register double R1 = A[ta]; // ra00
                        register double R2 = A[tta]; // r10
                        register double R3 = A[ttta]; // 20

                        register double R4 = B[tb]; // rb00
                        register double R5 = B[tb+1]; // rb01
                        register double R6 = B[tb+2]; // rb02

                        rc00 += R1 * R4;
                        rc01 += R1 * R5;
                        rc02 += R1 * R6;

                        rc10 += R2 * R4;
                        rc11 += R2 * R5;
                        rc12 += R2 * R6;

                        rc20 += R3 * R4;
                        rc21 += R3 * R5;
                        rc22 += R3 * R6;

                        R1 = A[ta+1];
                        R2 = A[tta+1];
                        R3 = A[ttta+1];

                        R4 = B[ttb];
                        R5 = B[ttb+1];
                        R6 = B[ttb+2];

                        rc00 += R1 * R4;
                        rc01 += R1 * R5;
                        rc02 += R1 * R6;

                        rc10 += R2 * R4;
                        rc11 += R2 * R5;
                        rc12 += R2 * R6;

                        rc20 += R3 * R4;
                        rc21 += R3 * R5;
                        rc22 += R3 * R6;

                        R1 = A[ta+2];
                        R2 = A[tta+2];
                        R3 = A[ttta+2];

                        R4 = B[tttb];
                        R5 = B[tttb+1];
                        R6 = B[tttb+2];

                        rc00 += R1 * R4;
                        rc01 += R1 * R5;
                        rc02 += R1 * R6;

                        rc10 += R2 * R4;
                        rc11 += R2 * R5;
                        rc12 += R2 * R6;

                        rc20 += R3 * R4;
                        rc21 += R3 * R5;
                        rc22 += R3 * R6;

                    }

                    C4[t] = rc00;
                    C4[t+1] = rc01;
                    C4[t+2] = rc02;

                    C4[tt] = rc10;
                    C4[tt+1] = rc11;
                    C4[tt+2] = rc12;

                    C4[ttt] = rc20;
                    C4[ttt+1] = rc21;
                    C4[ttt+2] = rc22;

                }
            }
        }

        void verifyResults(double *Matrix1, double *Matrix2, string examType){

            double maxDiff = 0;
            double tempDiff = 0;


            for (int i=0; i<N; i++) {
                for (int j=0; j<N; j++) {
                    tempDiff = Matrix1[i*N+j] - Matrix2[i*N+j];
                    if (tempDiff>maxDiff){
                        maxDiff = tempDiff;
                    }
                }
            }

            printf("Maximum Difference for Experiment %s : %f\n",examType.c_str(),maxDiff);
        }

        void verifyResults2(double *Matrix1, double *Matrix2, string examType){

            double Diff = fabs(Matrix1[0]-Matrix2[0]);
            double tempDiff;

            double maxA = fabs(A[0]);
            double maxB = fabs(B[0]);

            for (int i=0; i<N*N; i++){
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
            printf("Maximum Difference 2 for Experiment %s : %f\n",examType.c_str() ,maxDiff);
        }

        void runAndVerify(){
            drsvrTimer myTimer = drsvrTimer(N);

            myTimer.startTimer();
            this->dgemm0();
            myTimer.stopTimer();
            myTimer.printRunTime("dgmm0");


            myTimer.startTimer();
            this->dgemm1();
            myTimer.stopTimer();
            myTimer.printRunTime("dgmm1");

            myTimer.startTimer();
            this->dgemm2();
            myTimer.stopTimer();
            myTimer.printRunTime("dgmm2");

            myTimer.startTimer();
            this->dgemm3();
            myTimer.stopTimer();
            myTimer.printRunTime("dgmm3");

            this->verifyResults(C1, C2,"PART1");
            this->verifyResults(C1, C3,"PART2");
            this->verifyResults(C1, C4,"PART3");

            this->verifyResults2(C1, C2,"PART1");
            this->verifyResults2(C1, C3,"PART2");
            this->verifyResults2(C1, C3,"PART3");
        }

        void setTotalElement(){
            boundExtension = 0;
            int tempN = N;
            while (tempN%3!=0){
                tempN++;
                boundExtension++;
            }
            totalElements = tempN*tempN;
            //printf("BoundExtension: %d\ttotalElements: %d\tNew Dimension: %d * %d\n",boundExtension,totalElements,tempN,tempN);
        }

    public:
        matrixTestbench(int n){

            printf("\n***************************************\n");
            printf("TestBench with n = %d has been started!\n",n);

            N = n;

            //totalElements = N*N; // Number of element inside a Matrix

            this->setTotalElement();

            // Allocate Memory for each Matrix
            A = (double*)calloc(totalElements,sizeof(double));
            B = (double*)calloc(totalElements,sizeof(double));
            C1 = (double*)calloc(totalElements,sizeof(double));
            C2 = (double*)calloc(totalElements,sizeof(double));
            C3 = (double*)calloc(totalElements,sizeof(double));
            C4 = (double*)calloc(totalElements,sizeof(double));

            this->initializeMatrices();

            this->runAndVerify();

            printf("***************************************\n");
        }
};


int main() {

    matrixTestbench HW1_64 = matrixTestbench(64);
    matrixTestbench HW1_128 = matrixTestbench(128);
    matrixTestbench HW1_256 = matrixTestbench(256);
    matrixTestbench HW1_512 = matrixTestbench(512);
    matrixTestbench HW1_1024 = matrixTestbench(1024);
    matrixTestbench HW1_2048 = matrixTestbench(2048);

    return 0;
}
