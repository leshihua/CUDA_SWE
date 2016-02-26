#include <math.h>
#include <stdio.h>
#include <string>
#include "kernel.h"

int main() {
    const float x_lower = -3.f;
    const float x_upper = 3.f;
    //mesh x does not include ghost cells
    const float dx = (x_upper-x_lower)/MX;
    const float* const x = generateMesh(x_lower,x_upper);
    const float tfinal = 0.5;
    const int nsteps = 400; //number of time steps
    //const float cfl = 0.9;
    //I.C.
    //Shallow2* q = qinit(MX+2*MBC,0);
    Shallow2* q = qinit(MX+2*MBC,3);
    const float dt = tfinal/nsteps;
    godunov_serial(q, x, dt, nsteps, dx);
}
