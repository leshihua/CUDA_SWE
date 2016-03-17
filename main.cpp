#include <math.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include "kernel.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <time.h>

int main() {
    const float x_lower = -5.f;
    const float x_upper = 5.f;
    const float tfinal = 1.f;
    //const float x_lower = -500.f;
    //const float x_upper = 500.f;
    //const float tfinal = 100.f;
    const float dx = (x_upper-x_lower)/MX;
    float* x = (float*) malloc( (MX+2*MBC)*sizeof(float));
    generateMesh(x,x_lower,x_upper);
    //const float cfl = 0.9;

    //I.C.
    const float dt = tfinal/NSTEPS;
    const int caseNo = 4;
    int m = MX+2*MBC;

    //CPU implementation
    clock_t start0 = clock();
    Shallow2* q0 = (Shallow2*) malloc(m*sizeof(Shallow2));
    qinit(q0,m,caseNo);
    godunov_serial(q0, x, dt, NSTEPS, dx);
    free(q0);
    clock_t stop0 = clock();
    double timeInMs0 = ((double) (stop0 - start0)/CLOCKS_PER_SEC*1000.f);
    std::cout<<"Time elapsed for CPU implementation (ms): "<<timeInMs0<<std::endl; 

    //cudaSetDevice(0);
    //query device properties
    //cudaDeviceProp prop;
    //cudaGetDeviceProperties(&prop, 0);
    //std::cout<<"Running on device: "<<prop.name<<std::endl;
    //std::cout<<"Clock frequency in kilohertz: "<<prop.clockRate<<std::endl;
    //std::cout<<"Size of global memmory available on device in MB: "<<prop.totalGlobalMem/1024.0/1024.0<<std::endl;

    //GPU implementation with global memory
    cudaEvent_t start1, stop1;
    cudaEventCreate(&start1);
    cudaEventCreate(&stop1);
    cudaEventRecord(start1);
    Shallow2* q1 = (Shallow2*) malloc(m*sizeof(Shallow2));
    qinit(q1,m,caseNo);
    godunov_parallel_global_memory(q1, x, dt, NSTEPS, dx);
    free(q1);
    cudaEventRecord(stop1);
    cudaEventSynchronize(stop1);
    float timeInMs1 = 0;
    cudaEventElapsedTime(&timeInMs1, start1, stop1);
    std::cout<<"Time elapsed for GPU implementation with global memory (ms): "<<timeInMs1<<std::endl;

    //GPU implementation with shared memory
    cudaEvent_t start2, stop2;
    cudaEventCreate(&start2);
    cudaEventCreate(&stop2);
    cudaEventRecord(start2);
    Shallow2* q2 = (Shallow2*) malloc(m*sizeof(Shallow2));
    qinit(q2,m,caseNo);
    godunov_parallel_shared_memory(q2, x, dt, NSTEPS, dx);
    cudaEventRecord(stop2);
    cudaEventSynchronize(stop2);
    float timeInMs2 = 0;
    cudaEventElapsedTime(&timeInMs2, start2, stop2);
    std::cout<<"Time elapsed for GPU implementation with shared memory (ms): "<<timeInMs2<<std::endl;
    free(q2);

    ////GPU implementation with shared memory
    //Shallow2* q2 = (Shallow2*) malloc(m*sizeof(Shallow2));
    //qinit(q2,m,caseNo);
    //godunov_parallel_shared_memory(q2, x, dt, NSTEPS, dx);
    //free(q2);

    free(x);
}

