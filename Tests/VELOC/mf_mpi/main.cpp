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
    std::string check_file("chk");
    int restart_step = -1;
    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);
        pp.query("check_file", check_file);
        pp.query("restart_step", restart_step);
        pp.query("noutfiles", n_files);
    }

    VisMF::SetNOutFiles(n_files);

    if (restart_step < 0)
    {
        BoxArray ba(Box(IntVect(0),IntVect(n_cell-1)));
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);
        MultiFab mf(ba, dm, 1, 0);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const auto& arr = mf.array(mfi);
//            amrex::CheckSeedArraySizeAndResize(bx.numPts());
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

        amrex::UtilCreateDirectoryDestructive("vismfdata");

        const int nwork = 50;

        amrex::Print() << " Standard " << std::endl;

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

#ifdef AMREX_MPI_MULTIPLE
        amrex::Print() << " Async-MPI " << std::endl; 

        WriteAsyncStatus status_mpi;
        {
            BL_PROFILE_REGION("vismf-async-mpi-overlap");
            auto wrt_future = VisMF::WriteAsyncMPI(mf, "vismfdata/mf3");
            {
                BL_PROFILE_VAR("vismf-async-mpi-work", blp2);
                for (int i = 0; i < nwork; ++i) {
                    amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                    amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                }
            }
            {
                BL_PROFILE_VAR("vismf-async-mpi-wait", blp3);
                wrt_future.wait();
                status_mpi = wrt_future.get();
            }
        }
#endif

        for (int ip = 0; ip < ParallelDescriptor::NProcs(); ++ip) {
            if (ip == ParallelDescriptor::MyProc()) {
                amrex::AllPrint() << "Proc. " << ip << std::endl;
                amrex::AllPrint() << "File: " << status_file << std::endl;
#ifdef AMREX_MPI_MULTIPLE
                amrex::AllPrint() << "MPI: " << status_mpi << std::endl;
#endif
            }
            ParallelDescriptor::Barrier();
        }
    }
    else
    {
    }
}
