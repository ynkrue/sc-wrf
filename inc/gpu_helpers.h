// Shared CUDA utilities for WRF GPU modules (LW, SW, microphysics, …)
#ifndef GPU_HELPERS_H
#define GPU_HELPERS_H

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#define CUDA_CHECK(call)                                                    \
    do {                                                                    \
        cudaError_t err__ = (call);                                         \
        if (err__ != cudaSuccess) {                                         \
            fprintf(stderr, "CUDA error %s:%d: %s\n",                       \
                    __FILE__, __LINE__, cudaGetErrorString(err__));         \
            exit(1);                                                        \
        }                                                                   \
    } while (0)

int gpu_device_id();

extern "C" void wrf_gpu_init();

// ---- small buffer helpers ------------------------------------------------
static inline float* to_dev(const float *h, size_t n) {
    float *d;
    CUDA_CHECK(cudaMalloc((void**)&d, n * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d, h, n * sizeof(float), cudaMemcpyHostToDevice));
    return d;
}
static inline float* dev_zero(size_t n) {
    float *d;
    CUDA_CHECK(cudaMalloc((void**)&d, n * sizeof(float)));
    CUDA_CHECK(cudaMemset(d, 0, n * sizeof(float)));
    return d;
}
static inline int* dev_zero_int(size_t n) {
    int *d;
    CUDA_CHECK(cudaMalloc((void**)&d, n * sizeof(int)));
    CUDA_CHECK(cudaMemset(d, 0, n * sizeof(int)));
    return d;
}

#endif
