#ifndef KERNEL_H
#define KERNEL_H

#define MX 200// number of cells
#define MBC 1 //number of ghost cells at one side
//#define NUM_OUTPUT_TIMES 5 //number of output times
#define OUTPUT_FREQUENCY 40// write output every OUTPUT_FREQUENCY steps 
#define TPB 64
#define RAD 1 // radius of the stencil
#define GRAVITY 9.81
#define EFIX true

struct Shallow2 // a 2-d vector for quantities of 1-d shallow water problem 
{
    float h;
    float hu;
};

Shallow2* qinit(int m, int caseNo);  //should include ghost cells

void rp1(Shallow2* q, Shallow2* amdq, Shallow2* apdq, Shallow2* wave1, Shallow2* wave2, float* s1, float* s2, float dt, float dx, bool efix); // q should include ghost cells
float* generateMesh(float x_lower, float x_upper);
void godunov_serial(Shallow2 *q, const float* const x, const float dt, const int nsteps, const float dx);
void apply_BC(Shallow2* q);


#endif
