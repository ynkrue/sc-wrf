#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

namespace gpu {

__global__ void wrf_gpu_lw_kernel(float *in, float *out, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n) {
        out[i] = in[i] * 2.0f;
    }
}

}

#define CUDA_CHECK(call)                                                     \
    do {                                                                     \
        cudaError_t err__ = (call);                                          \
        if (err__ != cudaSuccess) {                                          \
            fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__,    \
                    cudaGetErrorString(err__));                              \
            abort();                                                         \
        }                                                                    \
    } while (0)

// Pick this rank's GPU: cyclic over whatever devices are visible.
static int pick_device() {
    int ndev = 0;
    if (cudaGetDeviceCount(&ndev) != cudaSuccess || ndev <= 0) {
        fprintf(stderr, "wrf_gpu: no CUDA devices visible\n");
        abort();
    }
    const char *s = getenv("SLURM_LOCALID");
    if (!s) s = getenv("OMPI_COMM_WORLD_LOCAL_RANK");
    int local = s ? atoi(s) : 0;
    int dev = local % ndev;

    char bus[16] = "?";
    cudaDeviceGetPCIBusId(bus, sizeof(bus), dev);
    fprintf(stderr, "wrf_gpu: local_rank %d -> device %d/%d (%s)\n",
            local, dev, ndev, bus);
    return dev;
}

static int device_id() {
    static const int d = pick_device(); // c++17 static init is thread-safe
    return d;
}

extern "C" void wrf_gpu_lw(const float *x, const int *n) {
    CUDA_CHECK(cudaSetDevice(device_id()));

    int nn = *n;
    float *d_x, *d_y;
    CUDA_CHECK(cudaMalloc((void**)&d_x, nn * sizeof(float)));
    CUDA_CHECK(cudaMalloc((void**)&d_y, nn * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_x, x, nn * sizeof(float), cudaMemcpyHostToDevice));

    int blockSize = 256;
    int numBlocks = (nn + blockSize - 1) / blockSize;
    gpu::wrf_gpu_lw_kernel<<<numBlocks, blockSize>>>(d_x, d_y, nn);
    CUDA_CHECK(cudaGetLastError());

    // CUDA_CHECK(cudaMemcpy(y, d_y, nn * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_y));
}
