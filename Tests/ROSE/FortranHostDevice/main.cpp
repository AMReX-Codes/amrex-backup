#include <AMReX.H>
#include <AMReX_Gpu.H>

using namespace amrex;

extern "C" {
    AMREX_GPU_HOST_DEVICE void f (int* x);
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        amrex::launch_global<<<1,1>>>(
        [] AMREX_GPU_DEVICE {
            int x = 1;
            f(&x);
            printf("From device, x = %d (should be 2)\n", x);
        });

        amrex::Gpu::Device::synchronize();

        // But I cannot call f from the host.  It seems that 
        // the function call becomes a no-op.
        int y = -1;
        f(&y);
        printf("From host, y = %d (should be -2)\n", y);
    }

    amrex::Finalize();
}
