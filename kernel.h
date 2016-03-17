#ifndef KERNEL_H
#define KERNEL_H

//For standard testing
//#define MX 20000// number of cells
//#define NSTEPS 40000
//#define OUTPUT_FREQUENCY 40000// write output every OUTPUT_FREQUENCY steps 

#define MX 100// number of cells
#define NSTEPS 400//number of time steps
#define OUTPUT_FREQUENCY 40// write output every OUTPUT_FREQUENCY steps 

#define MBC 1 //number of ghost cells at one side
#define GRAVITY 9.81
#define EFIX true

struct Shallow2 // a 2-d vector for quantities of 1-d shallow water problem 
{
    float h;
    float hu;
};

void qinit(Shallow2* q0, int m, int caseNo);  //should include ghost cells

void rp1(Shallow2* q, Shallow2* amdq, Shallow2* apdq, Shallow2* wave1, Shallow2* wave2, float* s1, float* s2, float dt, float dx, bool efix); // q should include ghost cells
void generateMesh(float* x, float x_lower, float x_upper);
void godunov_serial(Shallow2 *q, const float* const x, const float dt, const int nsteps, const float dx);
void godunov_parallel_global_memory(Shallow2 *q, const float* const x, const float dt, const int nsteps, const float dx);
void godunov_parallel_shared_memory(Shallow2 *q, const float* const x, const float dt, const int nsteps, const float dx);
//void apply_BC(Shallow2* q);



#endif
