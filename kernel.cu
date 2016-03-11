#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include "kernel.h"
//#include <cuda_runtime_api.h>
//#include <cuda.h>
#define TPB 256 
//#define TPB 4 



void qinit(Shallow2* q0, int m, int caseNo) //should include ghost cells
{
    float hl  = 0.f; 
    float hr  = 0.f; 
    float uhl = 0.f;
    float uhr = 0.f;
    if (0 == caseNo)//1-shock,2-shock
    {
        hl = 1.f;
        hr = 1.f;
        uhl = 3.f;
        uhr = -3.f;
    }
    else if (1 == caseNo)//1-rarefaction, 2-shocek
    {//no transonic wave
        hl = 3.f;
        hr = 1.f;
        uhl = 2.f;
        uhr = 1.f;
    }
    else if (2 == caseNo)//1-shock, 2-rarefaction
    {//transonic wave in 2-wave
        hl = 0.5;
        hr = 1.5;
        uhl = -1.5;
        uhr = -4.5;
    }
    else if (3 == caseNo)//1-rarefaction, 2-shock
    {//transonic wave in 1-wave
        hl = 7.f;
        hr = 1.f;
        uhl = 2.f;
        uhr = 1.f;
    }
    else if (4 == caseNo)//to test if q is initialized correctly
    {
        for (int i = 0; i < m; i++){
            q0[i].h = i;
            q0[i].hu = 10*i;
        }
        return;
    }
    else
    {
        std::cout<<"Invalid case number in qinit()."<<std::endl;
        exit (EXIT_FAILURE);
    }

    for (int i = 0; i < m; i++) {
        q0[i].h = hr;
        q0[i].hu = uhr;
    }

    for (int i = 0; i < (m/2); i++) {
        q0[i].h = hl;
        q0[i].hu = uhl;
    }
    return;
}
void generateMesh(float* x, float x_lower, float x_upper)
{
    float dx = (x_upper - x_lower)/MX;
    for (int i = MBC; i < (MX+MBC); i++)
    {
        x[i] = x_lower + (i-MBC)*dx + dx/2;
    }
    for (int i = 0; i < MBC; i++)
    {
        x[i] = x_lower - (MBC-i)*dx + dx/2;
        x[MX+MBC+i] = x_upper + (i+1)*dx + dx/2;
    }
    return; 
}


void rp1(Shallow2* q, Shallow2* amdq, Shallow2* apdq, Shallow2* wave1, Shallow2* wave2, float* s1, float* s2, float dt, float dx, bool efix) // q should include ghost cells
{
    /*
        numbering of cell and interface
        |--0--|--1--|--2--|........|--MX--|--MX+1--|
              0     1     2  ...  MX-1    MX  
        Riemann problems are solve at each interface 
        i ranges from 0 to MX; 
        index of q ranges from 0 to MX+1
    */
    //apply B.C.
    q[0].h = q[1].h;
    q[MX+MBC].h = q[MX].h;
    q[0].hu = q[1].hu;
    q[MX+MBC].hu = q[MX].hu;
    for (int i = 0; i < (MX+1); i++) 
    {
        //compute  Roe-averaged quantities:
        float ubar = (q[i].hu/sqrt(q[i].h) + q[i+1].hu/sqrt(q[i+1].h)) / (sqrt(q[i].h) + sqrt(q[i+1].h));  
        float cbar = sqrt(GRAVITY*0.5*(q[i].h + q[i+1].h));
        float dq1 = q[i+1].h - q[i].h;
        float dq2 = q[i+1].hu - q[i].hu;
        float a1 = 0.5*((ubar+cbar)*dq1-dq2)/cbar;
        float a2 = 0.5*(-(ubar-cbar)*dq1+dq2)/cbar;

        //compute the waves
        wave1[i].h = a1;
        wave1[i].hu = a1*(ubar - cbar);
        s1[i] = ubar - cbar;

        wave2[i].h = a2;
        wave2[i].hu = a2*(ubar + cbar);
        s2[i] = ubar + cbar;
    }
    if (efix)// entropy fix
    {
    //compute flux differences amdq and apdq.
    //First compute amdq as sum of s*wave for left going waves.
    //Incorporate entropy fix by adding a modified fraction of wave
    //if s should change sign.
    //todo: print CFL number
        for (int i = 0; i < (MX+1); i++) 
        {

            float s1_l = q[i].hu/q[i].h - sqrt(GRAVITY*q[i].h);//slope of 1-character, u-c, to the left of 1-wave
            //1.
            if ((s1_l >= 0.f) && (s1[i] > 0.f))//Fully supersonic case. Everything is right-going
            {
                amdq[i].h = 0.f;
                amdq[i].hu = 0.f;
                apdq[i].h = s1[i]*wave1[i].h + s2[i]*wave2[i].h;
                apdq[i].hu = s1[i]*wave1[i].hu + s2[i]*wave2[i].hu;
                continue;
            }

            float hr1 = q[i].h + wave1[i].h;
            float uhr1 = q[i].hu + wave1[i].hu;
            float s1_r = uhr1/hr1 - sqrt(GRAVITY*hr1);//slope of 1-character, u-c, to the right of 1-wave
            //2.
            if ((s1_l < 0.f) && (s1_r > 0.f))//transonic rarefaction in the 1-wave
            {
                float beta = (s1_r-s1[i])/(s1_r-s1_l); //eq (15.49) in the FVMHP book
                amdq[i].h = beta*s1_l*wave1[i].h;
                amdq[i].hu = beta*s1_l*wave1[i].hu;
                apdq[i].h = (1-beta)*s1_r*wave1[i].h + s2[i]*wave2[i].h; 
                apdq[i].hu = (1-beta)*s1_r*wave1[i].hu + s2[i]*wave2[i].hu; 
                continue;
            }
            else if (s1[i] < 0.f)//1-wave is left-going
            {
                //then check 2-wave
                float s2_r = q[i+1].hu/q[i+1].h + sqrt(GRAVITY*q[i+1].h);//slope of 2-character, u+c, to the right of 2-wave
                float hl2 = q[i+1].h - wave2[i].h;
                float uhl2 = q[i+1].hu - wave2[i].hu;
                float s2_l = uhl2/hl2 + sqrt(GRAVITY*hl2);//slope of 2-character, u+c, to the left of 2-wave
            //3.
                if ((s2_l < 0.f) && (s2_r > 0.f))//transonic rarefaction in the 2-wave
                {
                    //std::cout<<"transonic rarefaction in the 2-wave"<<std::endl;
                    float beta = (s2_r-s2[i])/(s2_r-s2_l);
                    amdq[i].h = s1[i]*wave1[i].h + beta*s2_l*wave2[i].h;
                    amdq[i].hu = s1[i]*wave1[i].hu + beta*s2_l*wave2[i].hu;
                    apdq[i].h = (1-beta)*s2_r*wave2[i].h;
                    apdq[i].hu = (1-beta)*s2_r*wave2[i].hu;
                    continue;
                }
            //4.
                else if (s2[i] < 0.f)//2-wave is left-going
                {
                    amdq[i].h = s1[i]*wave1[i].h + s2[i]*wave2[i].h;
                    amdq[i].hu = s1[i]*wave1[i].hu + s2[i]*wave2[i].hu;
                    apdq[i].h = 0.f;
                    apdq[i].hu = 0.f;
                    continue;
                }
            //5.
                else if (s2[i] > 0.f)//2-wave is right going
                {
                    amdq[i].h = s1[i]*wave1[i].h;
                    amdq[i].hu = s1[i]*wave1[i].hu;
                    apdq[i].h = s2[i]*wave2[i].h;
                    apdq[i].hu = s2[i]*wave2[i].hu;
                    continue;
                }
                else//this should not happend
                {
                    printf ("Conflict encountered when checking 2-wave");
                    exit (EXIT_FAILURE);
                }
            }
            //1-wave is going-right. This should not happen.
            else
            {
                printf ("Conflict encountered when checking 1-wave\n");
                exit (EXIT_FAILURE);
            }
        }

    }
    else
    {
        printf ("Case without entropy fix is not implemented.");
        exit (EXIT_FAILURE);
    }
    //update q
    for (int i = MBC; i < (MX+MBC); i++) 
    {
        q[i].h = q[i].h - dt/dx*(apdq[i-1].h + amdq[i].h); 
        q[i].hu = q[i].hu - dt/dx*(apdq[i-1].hu + amdq[i].hu); 
    }
}

void godunov_serial(Shallow2 *q, const float* const x, const float dt, const int nsteps, const float dx) 
{
    float t = 0.f;
    Shallow2* wave1 = (Shallow2*) malloc((MX+MBC)*sizeof(Shallow2));//1-wave
    Shallow2* wave2 = (Shallow2*) malloc((MX+MBC)*sizeof(Shallow2));//2-wave
    Shallow2* amdq = (Shallow2*) malloc((MX+MBC)*sizeof(Shallow2)); //A^- \delta Q defined at cell edge
    Shallow2* apdq = (Shallow2*) malloc((MX+MBC)*sizeof(Shallow2)); //A^+ \delta Q defined at cell edge
    float* s1 = (float*) malloc((MX+MBC)*sizeof(float));//speed of 1-wave
    float* s2 = (float*) malloc((MX+MBC)*sizeof(float));//speed of 2-wave
    std::cout<<"Solving SWE in serial."<<std::endl;
    //std::cout<<"dt = "<<std::to_string(dt)<<std::endl;
    //std::cout<<"dx = "<<std::to_string(dx)<<std::endl;
    std::cout<<"Writing data at step 0 to file."<<std::endl;
    std::ofstream outfile;
    outfile.open("output_serial_0");
    for (int j = MBC; j < MX+MBC; j++)
    {
        outfile<<std::to_string(x[j])<<","<<std::to_string(q[j].h)<<","<<std::to_string(q[j].hu)<<std::endl;
    }
    outfile.close();
    for (int n = 0; n < nsteps; n ++)
    {
        t = t + dt;
        //std::cout<<"Step: "<<std::to_string(n+1)<<std::endl;
        //std::cout<<"t = "<<std::to_string(t)<<std::endl;
        rp1(q, amdq, apdq, wave1, wave2, s1, s2, dt, dx, EFIX); // q should include ghost cells
        if ((n+1) % OUTPUT_FREQUENCY == 0 )//write output
        {
            std::cout<<"Writing data at t = "<<std::to_string(t)<<" to file."<<std::endl;
            std::ofstream outfile;
            outfile.open("output_serial_"+std::to_string(t));
            for (int j = MBC; j < MX+MBC; j++)
            {
                outfile<<std::to_string(x[j])<<","<<std::to_string(q[j].h)<<","<<std::to_string(q[j].hu)<<std::endl;
            }
            outfile.close();
        }
    }
    free(wave1);
    free(wave2);
    free(amdq);
    free(apdq);
    free(s1);
    free(s2);
}
        

//***************************************************************************//
//Below are parallel part using global memory
//***************************************************************************//

__global__
void rp1Kernel_global_memory(Shallow2* q, Shallow2* amdq, Shallow2* apdq, Shallow2* wave1, Shallow2* wave2, float* s1, float* s2, float dt, float dx, bool efix) // q should include ghost cells
{
    const int i = threadIdx.x + blockDim.x*blockIdx.x;
    if (i > (MX+MBC-1)) return;

    //update B.C.
    if (i < MBC)
    {
        q[i].h = q[i+1].h;
        q[i].hu = q[i+1].hu;
    }
    if (i > MX-1)
    {
        q[i+1].h = q[i].h;
        q[i+1].hu = q[i].hu;
    }
    __syncthreads();

    /*
       numbering of cell and interface
       |--0--|--1--|--2--|........|--MX--|--MX+1--|
             0     1     2  ...  MX-1    MX  
       Riemann problems are solve at each interface 
    */
    //compute  Roe-averaged quantities:
    float ubar = (q[i].hu/sqrt(q[i].h) + q[i+1].hu/sqrt(q[i+1].h)) / (sqrt(q[i].h) + sqrt(q[i+1].h));  
    float cbar = sqrt(GRAVITY*0.5*(q[i].h + q[i+1].h));
    float dq1 = q[i+1].h - q[i].h;
    float dq2 = q[i+1].hu - q[i].hu;
    float a1 = 0.5*((ubar+cbar)*dq1-dq2)/cbar;
    float a2 = 0.5*(-(ubar-cbar)*dq1+dq2)/cbar;

    //compute the waves
    wave1[i].h = a1;
    wave1[i].hu = a1*(ubar - cbar);
    s1[i] = ubar - cbar;

    wave2[i].h = a2;
    wave2[i].hu = a2*(ubar + cbar);
    s2[i] = ubar + cbar;
    if (efix)// entropy fix
    {
    //compute flux differences amdq and apdq.
    //First compute amdq as sum of s*wave for left going waves.
    //Incorporate entropy fix by adding a modified fraction of wave
    //if s should change sign.
    //todo: print CFL number
        float s1_l = q[i].hu/q[i].h - sqrt(GRAVITY*q[i].h);//slope of 1-character, u-c, to the left of 1-wave
        //1.
        if ((s1_l >= 0.f) && (s1[i] > 0.f))//Fully supersonic case. Everything is right-going
        {
            amdq[i].h = 0.f;
            amdq[i].hu = 0.f;
            apdq[i].h = s1[i]*wave1[i].h + s2[i]*wave2[i].h;
            apdq[i].hu = s1[i]*wave1[i].hu + s2[i]*wave2[i].hu;
        }
        else//check if 1-wave is transonic wave
        {
            float hr1 = q[i].h + wave1[i].h;
            float uhr1 = q[i].hu + wave1[i].hu;
            float s1_r = uhr1/hr1 - sqrt(GRAVITY*hr1);//slope of 1-character, u-c, to the right of 1-wave
        //2.
            if ((s1_l < 0.f) && (s1_r > 0.f))//transonic rarefaction in the 1-wave
            {
                float beta = (s1_r-s1[i])/(s1_r-s1_l); //eq (15.49) in the FVMHP book
                amdq[i].h = beta*s1_l*wave1[i].h;
                amdq[i].hu = beta*s1_l*wave1[i].hu;
                apdq[i].h = (1-beta)*s1_r*wave1[i].h + s2[i]*wave2[i].h; 
                apdq[i].hu = (1-beta)*s1_r*wave1[i].hu + s2[i]*wave2[i].hu; 
            }
            else if (s1[i] < 0.f)//1-wave is left-going. Then check 2-wave
            {
                
                float s2_r = q[i+1].hu/q[i+1].h + sqrt(GRAVITY*q[i+1].h);//slope of 2-character, u+c, to the right of 2-wave
                float hl2 = q[i+1].h - wave2[i].h;
                float uhl2 = q[i+1].hu - wave2[i].hu;
                float s2_l = uhl2/hl2 + sqrt(GRAVITY*hl2);//slope of 2-character, u+c, to the left of 2-wave
        //3.
                if ((s2_l < 0.f) && (s2_r > 0.f))//transonic rarefaction in the 2-wave
                {
                    float beta = (s2_r-s2[i])/(s2_r-s2_l);
                    amdq[i].h = s1[i]*wave1[i].h + beta*s2_l*wave2[i].h;
                    amdq[i].hu = s1[i]*wave1[i].hu + beta*s2_l*wave2[i].hu;
                    apdq[i].h = (1-beta)*s2_r*wave2[i].h;
                    apdq[i].hu = (1-beta)*s2_r*wave2[i].hu;
                }
        //4.
                else if (s2[i] < 0.f)//2-wave is left-going
                {
                    amdq[i].h = s1[i]*wave1[i].h + s2[i]*wave2[i].h;
                    amdq[i].hu = s1[i]*wave1[i].hu + s2[i]*wave2[i].hu;
                    apdq[i].h = 0.f;
                    apdq[i].hu = 0.f;
                }
        //5.
                else if (s2[i] > 0.f)//2-wave is right going
                {
                    amdq[i].h = s1[i]*wave1[i].h;
                    amdq[i].hu = s1[i]*wave1[i].hu;
                    apdq[i].h = s2[i]*wave2[i].h;
                    apdq[i].hu = s2[i]*wave2[i].hu;
                }
                else//this should not happend
                {
                    printf ("Conflict encountered when checking 2-wave\n");
                    //exit (EXIT_FAILURE);//not allowed in kernel
                    asm("trap;");  
                }
            }
            else//this should not happend
            {
                printf ("Conflict encountered when checking 1-wave\n");
                //exit (EXIT_FAILURE);//not allowed in kernel
                asm("trap;");  
            }
        }
    }
    __syncthreads();

        //update q
        //do not update q[0], which is ghost cell
    if (i < MX) //i+1 changes from 1 to MX
    {
        q[i+1].h = q[i+1].h - dt/dx*(apdq[i].h + amdq[i+1].h); 
        q[i+1].hu = q[i+1].hu - dt/dx*(apdq[i].hu + amdq[i+1].hu); 
    }
}

void godunov_parallel_global_memory(Shallow2 *q, const float* const x, const float dt, const int nsteps, const float dx) 
{
    int m = MX+2*MBC;
    float t = 0.f;
    Shallow2 *d_q = 0;
    cudaMalloc(&d_q, m*sizeof(Shallow2));
    cudaMemcpy(d_q,q,m*sizeof(Shallow2),cudaMemcpyHostToDevice);

    Shallow2 *d_wave1 = 0;
    Shallow2 *d_wave2 = 0;
    Shallow2 *d_amdq = 0;
    Shallow2 *d_apdq = 0;
    float *d_s1 = 0;
    float *d_s2 = 0;
    cudaMalloc(&d_wave1, (MX+MBC)*sizeof(Shallow2));
    cudaMalloc(&d_wave2, (MX+MBC)*sizeof(Shallow2));
    cudaMalloc(&d_amdq, (MX+MBC)*sizeof(Shallow2));
    cudaMalloc(&d_apdq, (MX+MBC)*sizeof(Shallow2));
    cudaMalloc(&d_s1, (MX+MBC)*sizeof(float));
    cudaMalloc(&d_s2, (MX+MBC)*sizeof(float));

    std::cout<<"Solving SWE in parallel with global memory."<<std::endl;
    //std::cout<<"Write data at step 0 to file."<<std::endl;
    std::ofstream outfile;
    outfile.open("output_global_0");
    for (int j = MBC; j < MX+MBC; j++)//write initial value of q
    {
        outfile<<std::to_string(x[j])<<","<<std::to_string(q[j].h)<<","<<std::to_string(q[j].hu)<<std::endl;
    }
    outfile.close();

    for (int n = 0; n < nsteps; n ++)//update q
    {
        t = t + dt;
        //std::cout<<"Step: "<<std::to_string(n+1)<<std::endl;
        //std::cout<<"t = "<<std::to_string(t)<<std::endl;
        rp1Kernel_global_memory<<<(m+TPB-1)/TPB,TPB>>>(d_q, d_amdq, d_apdq, d_wave1, d_wave2, d_s1, d_s2, dt, dx, EFIX);// q should include ghost cells
        if ((n+1) % OUTPUT_FREQUENCY == 0 )//write output at specified time
        {
            std::cout<<"Writing data at t = "<<std::to_string(t)<<" to file."<<std::endl;
            std::ofstream outfile;
            outfile.open("output_global_"+std::to_string(t));
            cudaMemcpy(q,d_q,m*sizeof(Shallow2),cudaMemcpyDeviceToHost);
            for (int j = MBC; j < MX+MBC; j++)
            {
                outfile<<std::to_string(x[j])<<","<<std::to_string(q[j].h)<<","<<std::to_string(q[j].hu)<<std::endl;
            }
            outfile.close();
        }
    }
    cudaFree(d_q);
    cudaFree(d_wave1);
    cudaFree(d_wave2);
    cudaFree(d_amdq);
    cudaFree(d_apdq);
    cudaFree(d_s1);
    cudaFree(d_s2);
}

//***************************************************************************//
//Below are parallel part using shared memory for input array
//***************************************************************************//

//shared memory for only input array
__global__
void rp1Kernel_shared_memory(Shallow2* q, Shallow2* amdq, Shallow2* apdq, Shallow2* wave1, Shallow2* wave2, float* s1, float* s2, float dt, float dx, bool efix) // q should include ghost cells
{
    const int i = threadIdx.x + blockDim.x*blockIdx.x;
    /*
    numbering of cell and interface

    /-----------shared memory in block 0------------\
    |--0--|--1--|--2--|.......|----TPB----|-(TPB+1)-|
          0     1     2    ........      TPB    

                              /------------shared memory in block 1--------------\
                              |-----0-----|----1----| ...|-----TPB-----|-(TPB+1)-|
                                          0         1  .........      TPB    

                                                         /----------shared memory in block 2 .... 
                                                         |------0------|----1----| ......
                                                                       0         1 ......

          /------------------------------------global memory------------------------------------------\
          |--0--|--1--|.......|--(TPB-1)--|---TPB---|....|--(2*TPB-1)--|--2*TPB--|....|--MX--|--MX+1--|
                |     |       |           |         |                  |         |           | 
                0     1 .... TPB-2     TPB-1       TPB .........   (2*TPB-1)   2*TPB ...... MX  

        Riemann problems are solve at each interface: i ranges from 0 to MX; 
        In each block (e.g. in 1st block):
        i ranges from 0 to TPB-1 => s_i ranges from 1 to TBP
        TODO: we don't really need info in s_q[0]
    */
    if (i > (MX+MBC-1)) return; 
    //update B.C.
    if (i < MBC)//i = 0
    {
        q[i].h = q[i+1].h;
        q[i].hu = q[i+1].hu;
    }
    if (i > MX-1) //i = MX
    {
        q[i+1].h = q[i].h;
        q[i+1].hu = q[i].hu;
    }
    __syncthreads();//To make the loading from global to shared concurrently

    const int s_i = threadIdx.x+1;//shared index
    extern __shared__ Shallow2 s_q[];//q on shared memory
    //load regular cells
    s_q[s_i] = q[i];
    //load halo cells
    //able to handle cases where MX is not multiple of TPB
    if (i == MX) {//last thread in last block
        s_q[0] = q[i-threadIdx.x-1];
        s_q[2+threadIdx.x] = q[i+1];
    }
    else if ( (i+1)%blockDim.x == 0 ) { //use the last thread of in a block to load halo cell, except for last block
        if (i > blockDim.x) {//in 1st block, don't load q[-1] to s_q[0]
            s_q[0] = q[i-blockDim.x];
        }
        s_q[1+blockDim.x] = q[i+1];//load right halo cell
    }
    __syncthreads();//finish reading all data into shared memory

    //compute  Roe-averaged quantities:
    float ubar = (s_q[s_i].hu/sqrt(s_q[s_i].h) + s_q[s_i+1].hu/sqrt(s_q[s_i+1].h)) / (sqrt(s_q[s_i].h) + sqrt(s_q[s_i+1].h));  
    float cbar = sqrt(GRAVITY*0.5*(s_q[s_i].h + s_q[s_i+1].h));
    float dq1 = s_q[s_i+1].h - s_q[s_i].h;
    float dq2 = s_q[s_i+1].hu - s_q[s_i].hu;
    float a1 = 0.5*((ubar+cbar)*dq1-dq2)/cbar;
    float a2 = 0.5*(-(ubar-cbar)*dq1+dq2)/cbar;

    //compute the waves
    wave1[i].h = a1;
    wave1[i].hu = a1*(ubar - cbar);
    s1[i] = ubar - cbar;

    wave2[i].h = a2;
    wave2[i].hu = a2*(ubar + cbar);
    s2[i] = ubar + cbar;
    if (efix)// entropy fix
    {
    //compute flux differences amdq and apdq.
    //First compute amdq as sum of s*wave for left going waves.
    //Incorporate entropy fix by adding a modified fraction of wave
    //if s should change sign.
    //todo: print CFL number
        float s1_l = s_q[s_i].hu/s_q[s_i].h - sqrt(GRAVITY*s_q[s_i].h);//slope of 1-character, u-c, to the left of 1-wave
        //1.
        if ((s1_l >= 0.f) && (s1[i] > 0.f))//Fully supersonic case. Everything is right-going
        {
            amdq[i].h = 0.f;
            amdq[i].hu = 0.f;
            apdq[i].h = s1[i]*wave1[i].h + s2[i]*wave2[i].h;
            apdq[i].hu = s1[i]*wave1[i].hu + s2[i]*wave2[i].hu;
        }
        else//check if 1-wave is transonic wave
        {
            float hr1 = s_q[s_i].h + wave1[i].h;
            float uhr1 = s_q[s_i].hu + wave1[i].hu;
            float s1_r = uhr1/hr1 - sqrt(GRAVITY*hr1);//slope of 1-character, u-c, to the right of 1-wave
        //2.
            if ((s1_l < 0.f) && (s1_r > 0.f))//transonic rarefaction in the 1-wave
            {
                float beta = (s1_r-s1[i])/(s1_r-s1_l); //eq (15.49) in the FVMHP book
                amdq[i].h = beta*s1_l*wave1[i].h;
                amdq[i].hu = beta*s1_l*wave1[i].hu;
                apdq[i].h = (1-beta)*s1_r*wave1[i].h + s2[i]*wave2[i].h; 
                apdq[i].hu = (1-beta)*s1_r*wave1[i].hu + s2[i]*wave2[i].hu; 
            }
            else if (s1[i] < 0.f)//1-wave is left-going. Then check 2-wave
            {
                
                float s2_r = s_q[s_i+1].hu/s_q[s_i+1].h + sqrt(GRAVITY*s_q[s_i+1].h);//slope of 2-character, u+c, to the right of 2-wave
                float hl2 = q[i+1].h - wave2[i].h;
                float uhl2 = q[i+1].hu - wave2[i].hu;
                float s2_l = uhl2/hl2 + sqrt(GRAVITY*hl2);//slope of 2-character, u+c, to the left of 2-wave
        //3.
                if ((s2_l < 0.f) && (s2_r > 0.f))//transonic rarefaction in the 2-wave
                {
                    float beta = (s2_r-s2[i])/(s2_r-s2_l);
                    amdq[i].h = s1[i]*wave1[i].h + beta*s2_l*wave2[i].h;
                    amdq[i].hu = s1[i]*wave1[i].hu + beta*s2_l*wave2[i].hu;
                    apdq[i].h = (1-beta)*s2_r*wave2[i].h;
                    apdq[i].hu = (1-beta)*s2_r*wave2[i].hu;
                }
        //4.
                else if (s2[i] < 0.f)//2-wave is left-going
                {
                    amdq[i].h = s1[i]*wave1[i].h + s2[i]*wave2[i].h;
                    amdq[i].hu = s1[i]*wave1[i].hu + s2[i]*wave2[i].hu;
                    apdq[i].h = 0.f;
                    apdq[i].hu = 0.f;
                }
        //5.
                else if (s2[i] > 0.f)//2-wave is right going
                {
                    amdq[i].h = s1[i]*wave1[i].h;
                    amdq[i].hu = s1[i]*wave1[i].hu;
                    apdq[i].h = s2[i]*wave2[i].h;
                    apdq[i].hu = s2[i]*wave2[i].hu;
                }
                else//this should not happend
                {
                    printf ("Conflict encountered when checking 2-wave\n");
                    //exit (EXIT_FAILURE);//not allowed in kernel
                    asm("trap;");  
                }
            }
            else//this should not happend
            {
                printf ("Conflict encountered when checking 1-wave\n");
                //exit (EXIT_FAILURE);//not allowed in kernel
                asm("trap;");  
            }
        }
    }
    __syncthreads();

        //update q
        //do not update q[0], which is ghost cell
    if (i < MX) //i+1 changes from 1 to MX
    {
        q[i+1].h = s_q[s_i+1].h - dt/dx*(apdq[i].h + amdq[i+1].h); 
        q[i+1].hu = s_q[s_i+1].hu - dt/dx*(apdq[i].hu + amdq[i+1].hu); 
    }
    __syncthreads();
}

//todo: finish rp1Kernel_shared_memory2
//shared memory for both input and output array

void godunov_parallel_shared_memory(Shallow2 *q, const float* const x, const float dt, const int nsteps, const float dx) 
{
    int m = MX+2*MBC;
    float t = 0.f;
    Shallow2 *d_q = 0;
    //void *d_q = NULL;
    cudaMalloc(&d_q, m*sizeof(Shallow2));
    cudaMemcpy(d_q,q,m*sizeof(Shallow2),cudaMemcpyHostToDevice);

    Shallow2 *d_wave1 = 0;
    Shallow2 *d_wave2 = 0;
    Shallow2 *d_amdq = 0;
    Shallow2 *d_apdq = 0;
    float *d_s1 = 0;
    float *d_s2 = 0;
    cudaMalloc(&d_wave1, (MX+MBC)*sizeof(Shallow2));
    cudaMalloc(&d_wave2, (MX+MBC)*sizeof(Shallow2));
    cudaMalloc(&d_amdq, (MX+MBC)*sizeof(Shallow2));
    cudaMalloc(&d_apdq, (MX+MBC)*sizeof(Shallow2));
    cudaMalloc(&d_s1, (MX+MBC)*sizeof(float));
    cudaMalloc(&d_s2, (MX+MBC)*sizeof(float));

    const size_t smemSize = (TPB + MBC)*sizeof(Shallow2);

    std::cout<<"Solving SWE in parallel with shared memory."<<std::endl;
    //std::cout<<"Write data at step 0 to file."<<std::endl;
    std::ofstream outfile;
    outfile.open("output_shared_0");
    for (int j = MBC; j < MX+MBC; j++)//write initial value of q
    {
        outfile<<std::to_string(x[j])<<","<<std::to_string(q[j].h)<<","<<std::to_string(q[j].hu)<<std::endl;
    }
    outfile.close();

    for (int n = 0; n < nsteps; n ++)//update q
    {
        t = t + dt;
        //std::cout<<"Step: "<<std::to_string(n+1)<<std::endl;
        //std::cout<<"t = "<<std::to_string(t)<<std::endl;
        rp1Kernel_shared_memory<<<((m-1)+TPB-1)/TPB,TPB,smemSize>>>(d_q, d_amdq, d_apdq, d_wave1, d_wave2, d_s1, d_s2, dt, dx, EFIX);// q should include ghost cells
        if ((n+1) % OUTPUT_FREQUENCY == 0 )//write output at specified time
        {
            std::cout<<"Writing data at t = "<<std::to_string(t)<<" to file."<<std::endl;
            std::ofstream outfile;
            outfile.open("output_shared_"+std::to_string(t));
            cudaMemcpy(q,d_q,m*sizeof(Shallow2),cudaMemcpyDeviceToHost);
            for (int j = MBC; j < MX+MBC; j++)
            {
                outfile<<std::to_string(x[j])<<","<<std::to_string(q[j].h)<<","<<std::to_string(q[j].hu)<<std::endl;
            }
            outfile.close();
        }
    }
    cudaFree(d_q);
    cudaFree(d_wave1);
    cudaFree(d_wave2);
    cudaFree(d_amdq);
    cudaFree(d_apdq);
    cudaFree(d_s1);
    cudaFree(d_s2);
}