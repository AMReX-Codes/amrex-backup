#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>

#include <thread>
#include <future>

using namespace amrex;

void main_main ();

int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();

}

void main_main ()
{
    BL_PROFILE("main");

    int n_cell = 512;
    int max_grid_size = 64;
    int n_files = 256;
    int nwork = 50;
    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);
        pp.query("noutfiles", n_files);
        pp.query("nwork", nwork);
    }

    VisMF::SetNOutFiles(n_files);
//    amrex::ResetRandomSeed(33344455666);
    BoxArray ba(Box(IntVect(0),IntVect(n_cell-1)));
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);
    MultiFab mf(ba, dm, 1, 0);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const auto& arr = mf.array(mfi);
        amrex::ParallelFor (bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = amrex::Random();
        });
        Gpu::streamSynchronize(); // because of random nubmer generator
    }

    { // touch pinned memory so that the one-time cost is removed from timers.
        MultiFab mfcpu(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect(),
                       MFInfo().SetArena(The_Pinned_Arena()));
        amrex::dtoh_memcpy(mfcpu, mf);
    }

    amrex::Print() << " Printing random box with n_cell = " << n_cell
                                    << ", max_grid_size = " << max_grid_size
                                    << ", and noutfiles = " << n_files << std::endl;


    for (int ip = 0; ip < ParallelDescriptor::NProcs(); ++ip) {
        if (ip == ParallelDescriptor::MyProc()) {
            amrex::AllPrint() << "Proc. " << ip << " number of boxes = " << mf.local_size() << std::endl;
        }
        ParallelDescriptor::Barrier();
    }

    amrex::UtilCreateDirectoryDestructive("vismfdata");

// ***************************************************************

    amrex::Print() << " No Async " << std::endl;
    {
        BL_PROFILE_REGION("vismf-orig");
        VisMF::Write(mf, "vismfdata/mf1");
        {
            BL_PROFILE_VAR("vismf-orig-work", blp2);
            for (int i = 0; i < nwork; ++i) {
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
            }
        }
    }
    ParallelDescriptor::Barrier();


// ***************************************************************

    amrex::Print() << " Async-file " << std::endl; 
    WriteAsyncStatus status_file;
    {
        BL_PROFILE_REGION("vismf-async-file-overlap");
        auto wrt_future = VisMF::WriteAsync(mf, "vismfdata/mf2");
        {
            BL_PROFILE_VAR("vismf-async-file-work", blp2);
            for (int i = 0; i < nwork; ++i) {
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-file-wait", blp3);
            wrt_future.wait();
            status_file = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

#ifdef AMREX_MPI_MULTIPLE

    amrex::Print() << " Async-MPI " << std::endl; 
    WriteAsyncStatus status_mpi_basic;
    {
        BL_PROFILE_REGION("vismf-async-mpi-basic-overlap");
        auto wrt_future = VisMF::WriteAsyncMPI(mf, "vismfdata/mf3");
        {
            BL_PROFILE_VAR("vismf-async-mpi-basic-work", blp2);
            for (int i = 0; i < nwork; ++i) {
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-basic-wait", blp3);
            wrt_future.wait();
            status_mpi_basic = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Comm " << std::endl; 
    WriteAsyncStatus status_mpi_comm;
    {
        BL_PROFILE_REGION("vismf-async-mpi-comm-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIComm(mf, "vismfdata/mf4");
        {
            BL_PROFILE_VAR("vismf-async-mpi-comm-work", blp2);
            for (int i = 0; i < nwork; ++i) {
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-comm-wait", blp3);
            wrt_future.wait();
            status_mpi_comm = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Wait " << std::endl; 
    WriteAsyncStatus status_mpi_wait;
    {
        BL_PROFILE_REGION("vismf-async-mpi-wait-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIWait(mf, "vismfdata/mf5");
        {
            BL_PROFILE_VAR("vismf-async-mpi-wait-work", blp2);
            for (int i = 0; i < nwork; ++i) {
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-wait-wait", blp3);
            wrt_future.wait();
            status_mpi_wait = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Fence " << std::endl; 
    WriteAsyncStatus status_mpi_fence;
    {
        BL_PROFILE_REGION("vismf-async-mpi-fence-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIOneSidedFence(mf, "vismfdata/mf6");
        {
            BL_PROFILE_VAR("vismf-async-mpi-fence-work", blp2);
            for (int i = 0; i < nwork; ++i) {
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-fence-wait", blp3);
            wrt_future.wait();
            status_mpi_fence = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Post " << std::endl; 
    WriteAsyncStatus status_mpi_post;
    {
        BL_PROFILE_REGION("vismf-async-mpi-post-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIOneSidedPost(mf, "vismfdata/mf7");
        {
            BL_PROFILE_VAR("vismf-async-mpi-post-work", blp2);
            for (int i = 0; i < nwork; ++i) {
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-post-wait", blp3);
            wrt_future.wait();
            status_mpi_post = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

#endif

    for (int ip = 0; ip < ParallelDescriptor::NProcs(); ++ip) {
        if (ip == ParallelDescriptor::MyProc()) {
            amrex::AllPrint() << "Proc. " << ip << std::endl;
            amrex::AllPrint() << "File: " << status_file << std::endl;
#ifdef AMREX_MPI_MULTIPLE
            amrex::AllPrint() << "MPI-Basic: " << status_mpi_basic << std::endl;
            amrex::AllPrint() << "MPI-Comm: "  << status_mpi_comm  << std::endl;
            amrex::AllPrint() << "MPI-Wait: "  << status_mpi_wait  << std::endl;
            amrex::AllPrint() << "MPI-Fence: " << status_mpi_fence << std::endl;
            amrex::AllPrint() << "MPI-Post: "  << status_mpi_post  << std::endl;
#endif
        }
        ParallelDescriptor::Barrier();
    }
}
