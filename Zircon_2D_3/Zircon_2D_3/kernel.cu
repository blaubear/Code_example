//#include "Utilities.cuh"
#include "stdio.h"
#include "math.h"

#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }


inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

/***********************************/
/* ITERATION FUNCTION - GPU */
/***********************************/
__global__ 
void Calculator_GPU(float* T_old, float* T_new, const int NX, const int NY, double D, double dx, double dt, double C_sat)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < NX-1 && i > 0 && j < NY-1 && j > 0)
        T_new[i + NY * j] = T_old[i + NY * j] + (dt/(dx*dx))* D * (T_old[i - 1 + NY * j] + T_old[i + 1 + NY * j] + T_old[i + (j - 1) * NY] + T_old[i + (j + 1) * NY] - 4.0 * T_old[i + NY * j]);
    if (i == 0)
        //T_new[i + NY * j] = T_new[i + 1 + NY * j];
        T_new[i + NY * j] = C_sat;
    if (i == (NX - 1) )
        //T_new[i + NY * j] = T_new[i - 1 + NY * j];
        T_new[i + NY * j] = C_sat;
    if (j == 0)
        //T_new[i + NY * j] = T_new[i + NY * (j + 1)];
        T_new[i + NY * j] = C_sat;
    if (j == (NY-1))
        //T_new[i + NY * j] = T_new[i + NY * (j - 1)];
        T_new[i + NY * j] = C_sat;
}


/***********************************/
/* ITERATION FUNCTION - GPU */
/***********************************/
__global__
void Calculator_GPU_X(float* T_old, float* T_new, const int NX, const int NY, double D, double dx, double dt, double C_sat)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    D = 0.1;
    if (i < NX - 1 && i > 0 && j < NY - 1 && j > 0)
        T_new[i + NY * j] = T_old[i + NY * j] + (dt / (dx * dx)) * D * (T_old[i - 1 + NY * j] + T_old[i + 1 + NY * j] + T_old[i + (j - 1) * NY] + T_old[i + (j + 1) * NY] - 4.0 * T_old[i + NY * j]);
    if (i == 0)
        //T_new[i + NY * j] = T_new[i + 1 + NY * j];
        T_new[i + NY * j] = C_sat;
    if (i == (NX - 1))
        T_new[i + NY * j] = T_new[i - 1 + NY * j];
    if (j == 0)
        //T_new[i + NY * j] = T_new[i + NY * (j + 1)];
        T_new[i + NY * j] = C_sat;
    if (j == (NY - 1))
        T_new[i + NY * j] = T_new[i + NY * (j - 1)];
}

/***********************************/
/* ITERATION FUNCTION - GPU */
/***********************************/
__global__
void Calculator_GPU_Y(float* T_old, float* T_new, const int NX, const int NY, float D, float dx, float dt)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    D = 0.1;
    if (i < NX - 1 && i > 0 && j < NY - 1 && j > 0)
        T_new[i + NY * j] = T_old[i + NY * j] + (dt / (dx * dx)) * D * (T_old[i - 1 + NY * j] + T_old[i + 1 + NY * j] + T_old[i + (j - 1) * NY] + T_old[i + (j + 1) * NY] - 4.0 * T_old[i + NY * j]);
    if (i == 0)
        //T_new[i + NY * j] = T_new[i + 1 + NY * j];
        T_new[i + NY * j] = 0;
    if (i == (NX - 1))
        T_new[i + NY * j] = T_new[i - 1 + NY * j];
    if (j == 0)
        //T_new[i + NY * j] = T_new[i + NY * (j + 1)];
        T_new[i + NY * j] = 0;
    if (j == (NY - 1))
        T_new[i + NY * j] = T_new[i + NY * (j - 1)];
}

/***********************************/
/*        Progonka - GPU      */
/***********************************/
__global__
void pronochka(float* T_old, float* T_new, const int NX, const int NY, float D, float dx, float dt)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    D = 0.1;
    if (i < NX - 1 && i > 0 && j < NY - 1 && j > 0)
        T_new[i + NY * j] = T_old[i + NY * j] + (dt / (dx * dx)) * D * (T_old[i - 1 + NY * j] + T_old[i + 1 + NY * j] + T_old[i + (j - 1) * NY] + T_old[i + (j + 1) * NY] - 4.0 * T_old[i + NY * j]);
    if (i == 0)
        //T_new[i + NY * j] = T_new[i + 1 + NY * j];
        T_new[i + NY * j] = 0;
    if (i == (NX - 1))
        T_new[i + NY * j] = T_new[i - 1 + NY * j];
    if (j == 0)
        //T_new[i + NY * j] = T_new[i + NY * (j + 1)];
        T_new[i + NY * j] = 0;
    if (j == (NY - 1))
        T_new[i + NY * j] = T_new[i + NY * (j - 1)];
}


/******************************/
/* TEMPERATURE INITIALIZATION */
/******************************/

void Initialize(float* h_T, const int NX, const int NY)
{
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            //h_T[i * NY +j ] = i/128.0;
            h_T[i * NY + j] = 300;
        }
    }
}

/******************************/
/*      Write in file 1       */
/******************************/
void WriteInFile_1(float* h_T_GPU_result, int NX, int NY) {
    // --- Write in file
    FILE* pointer = fopen("1.txt", "w");
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            fprintf(pointer, "%f ", h_T_GPU_result[j * NX + i]);
        }
        putc('\n', pointer);
    }
    fclose(pointer);
}

/******************************/
/*      Write in file 2       */
/******************************/
void WriteInFile_2(float* h_T_GPU_result, int NX, int NY) {
    // --- Write in file
    FILE* pointer = fopen("2.txt", "w");
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            fprintf(pointer, "%f ", h_T_GPU_result[j * NX + i]);
        }
        putc('\n', pointer);
    }
    fclose(pointer);
}

/********/
/* MAIN */
/********/
int main(){
    const int NX = 512;         // --- Number of discretization points along the x axis
    const int NY = NX;         // --- Number of discretization points along the y axis
    int Nt = 200000, MAX_ITER = 1;
    double Lx, Ly, D = 0.1;
    Lx = 0.1;
    Ly = Lx;
    double dx = Lx / NX;
    double T_start = 950 + 273; // начальная температура в кельвинах;
    double T_end = 750 + 273;  //конечная температура в кельвинах;
    double time_coef = (T_start - T_end);
    double t_end = 10000 * 365 * 24 * 60 * 60; // окончание по времени;
    //double t_end = 1;
    double X_H20 = 2;
    double dt = 1 / (double)Nt;
    double T, M, D_nd;
    double time = 0;
    double C_sat = 100, C_cryst = 490000;
    int Counter = 0;
    cudaEvent_t start, stop;
    float elapsedTime;

    float* GPU_D ;
    cudaMalloc((void**)&GPU_D, sizeof(float));
    cudaEventCreate(&start);
    cudaEventRecord(start, 0);

    // --- GPU temperature distribution
    float* h_T = (float*)calloc(NX * NY, sizeof(float));
    float* h_T_old = (float*)calloc(NX * NY, sizeof(float));
    Initialize(h_T, NX, NY);
    Initialize(h_T_old, NX, NY);
    float* h_T_GPU_result = (float*)malloc(NX * NY * sizeof(float));

    WriteInFile_1(h_T, NX, NY);

    // --- GPU temperature distribution
    float* d_T;     cudaMalloc((void**)&d_T, NX * NY * sizeof(float));
    float* d_T_old; cudaMalloc((void**)&d_T_old, NX * NY * sizeof(float));

    cudaMemcpy(d_T, h_T, NX * NY * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_T_old, d_T, NX * NY * sizeof(float), cudaMemcpyDeviceToDevice);
    
    dim3 threadsPerBlock(32, 32);
    int nbx = (NX / threadsPerBlock.x) + (((NX % threadsPerBlock.x) == 0) ? 0 : 1);
    int nby = (NY / threadsPerBlock.y) + (((NY % threadsPerBlock.y) == 0) ? 0 : 1);
    dim3 numBlocks(nbx, nby);

    // --- Iterations on the device
    for (int i = 0; i < Nt; i++) {
        Counter = Counter + 1;
        time = time + dt;
        T = T_start - time_coef * time;
        M = 4.8 * pow(10, -6) * pow(T, 2) - 8.4 * pow(10, -3) * T + 4.84;
        C_sat = C_cryst / (exp(10108 / T + 1.16 * (M - 1) - 1.48));

        //a = fzero(@(x_nd)pi ^ (1 / 2) * x_nd * exp(x_nd ^ 2) * erfc(x_nd) - (C_bound - C_sat) / (C_cryst - C_sat), 0);% вычисление промежуточной величины для аналитического решения;
        D = (exp(-(11.4 * X_H20 + 3.13) / (0.84 * X_H20 + 1) - ((21.4 * X_H20 + 47) / (1.06 * X_H20 + 1)) * (1000) / T));
        D_nd = D * t_end / pow(Lx, 2);
        Calculator_GPU << <numBlocks, threadsPerBlock >> > (d_T, d_T_old, NX, NY, D_nd, dx, dt, C_sat);   // --- Update d_T_old     starting from data stored in d_T

        //cudaMemcpy(d_T_old, d_T, NX * NY * sizeof(float), cudaMemcpyDeviceToDevice);
        d_T_old = d_T;
        if (fmod(Counter, Nt/10) == 0) {
        printf("%f \n", time);
        }
    }
    gpuErrchk(cudaPeekAtLastError());
    // --- Copy result from device to host
    cudaMemcpy(h_T_GPU_result, d_T, NX * NY * sizeof(float), cudaMemcpyDeviceToHost);

    // --- Write in file 2
    WriteInFile_2(h_T_GPU_result, NX, NY);

    // --- Release device memory
    cudaFree(d_T);
    cudaFree(d_T_old);

    free(h_T);
    free(h_T_GPU_result);

    cudaEventCreate(&stop);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);

    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("Elapsed time : %f ms\n", elapsedTime);


    return 0;
}
