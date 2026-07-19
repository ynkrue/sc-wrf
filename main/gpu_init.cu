#include "../inc/gpu_helpers.h"

int gpu_device_id() {
    static const int d = []() -> int {
        int ndev = 0;
        if (cudaGetDeviceCount(&ndev) != cudaSuccess || ndev <= 0) {
            fprintf(stderr, "wrf_gpu: no CUDA devices visible\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
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
    }();
    return d;
}

extern "C" void wrf_gpu_init() {
    CUDA_CHECK(cudaSetDevice(gpu_device_id()));
}
