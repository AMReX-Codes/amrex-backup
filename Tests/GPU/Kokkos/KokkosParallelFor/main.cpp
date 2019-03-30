#include <cuda_runtime.h>

#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

#include <Kokkos_Core.hpp>

// Current timers:
// Ignore resetting of value and output. All else included.


using namespace amrex;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

int main (int argc, char* argv[])
{

    amrex::Initialize(argc, argv); {
    {
        Kokkos::InitArguments args;
        args.device_id = amrex::Gpu::Device::deviceId();

        Kokkos::initialize(args);
    }
    {

        amrex::Print() << "amrex::Initialize complete." << "\n";

        // ===================================
        // Simple cuda action to make sure all tests have cuda.
        // Allows nvprof to return data.
        int devices = 0;
#ifdef AMREX_USE_CUDA
        cudaGetDeviceCount(&devices);
#endif
        amrex::Print() << "Hello world from AMReX version " << amrex::Version() << ". GPU devices: " << devices << "\n";
        amrex::Print() << "**********************************\n"; 
        // ===================================

        // What time is it now?  We'll use this to compute total run time.
        Real strt_time = amrex::second();

        // AMREX_SPACEDIM: number of dimensions
        int n_cell, max_grid_size;
        Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

        // inputs parameters
        {
            // ParmParse is way of reading inputs from the inputs file
            ParmParse pp;

            // We need to get n_cell from the inputs file - this is the number of cells on each side of 
            //   a square (or cubic) domain.
            pp.get("n_cell",n_cell);

            // The domain is broken into boxes of size max_grid_size
            pp.get("max_grid_size",max_grid_size);
        }

        // make BoxArray and Geometry
        BoxArray ba;
        Geometry geom;
        {
            IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
            IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
            Box domain(dom_lo, dom_hi);

            // Initialize the boxarray "ba" from the single box "bx"
            ba.define(domain);
            // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
            ba.maxSize(max_grid_size);

            // This defines the physical box, [-1,1] in each direction.
            RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                             {AMREX_D_DECL( 1.0, 1.0, 1.0)});

            // This defines a Geometry object
            geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
        }

        // Nghost = number of ghost cells for each array 
        int Nghost = 1;
    
        // Ncomp = number of components for each array
        int Ncomp  = 1;
  
        // How Boxes are distrubuted among MPI processes
        DistributionMapping dm(ba);

        // Malloc value for setval testing.
        Real* val;
        cudaMallocHost(&val, sizeof(Real));
        *val = 0.0;

        // Create the MultiFab and touch the data.
        // Ensures the data in on the GPU for all further testing.
        MultiFab x(ba, dm, Ncomp, Nghost);
        MultiFab y(ba, dm, Ncomp, Nghost);

        x.setVal(*val);
        y.setVal(*val);

        amrex::Print() << "MultiFab initialized with: " << n_cell << "^3, max_box_length = " << max_grid_size
                       << std::endl << std::endl;


        Real points = ba.numPts();

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Initial launch to remove any unknown costs in HtoD setup. 

        {
            BL_PROFILE("Kokkos One");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({lo.x,lo.y,lo.z},{hi.x+1,hi.y+1,hi.z+1}),
                KOKKOS_LAMBDA (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });

                // ..................

            }

//            amrex::Print() << "Kokkos one sum = " << result << ". Expected = " << points*((*val)+1) << std::endl;
        }

        {
            BL_PROFILE("Amrex One");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });
            }

        }

        {
            BL_PROFILE("Kokkos Two");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({lo.x,lo.y,lo.z},{hi.x+1,hi.y+1,hi.z+1}),
                KOKKOS_LAMBDA (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });

                // ..................

            }

        }

        {
            BL_PROFILE("Amrex Two");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
 
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });
            }

        }

        {
            BL_PROFILE("Kokkos Three");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({lo.x,lo.y,lo.z},{hi.x+1,hi.y+1,hi.z+1}),
                KOKKOS_LAMBDA (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });

                // ..................

            }

        }

        {
            BL_PROFILE("Amrex Three");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
 
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });
            }

        }

        int numStreams = 16;
        Kokkos::Cuda instances[numStreams];

        for(int i=0; i<numStreams; ++i) {
           amrex::Cuda::Device::setStreamIndex(i);
           instances[i] = Kokkos::Cuda( amrex::Cuda::Device::cudaStream() );
        }


        {
            BL_PROFILE("Kokkos Streams One");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);
                int localIndex = mfi.LocalIndex();

                Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>(instances[localIndex%numStreams],
                                                                            {lo.x,lo.y,lo.z},{hi.x+1,hi.y+1,hi.z+1}),
                KOKKOS_LAMBDA (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });

                // ..................

            }

        }

        {
            BL_PROFILE("Amrex Four");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
 
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });
            }

        }

        {
            BL_PROFILE("Kokkos Streams Two");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);
                int localIndex = mfi.LocalIndex();

                Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>(instances[localIndex%numStreams],
                                                                            {lo.x,lo.y,lo.z},{hi.x+1,hi.y+1,hi.z+1}),
                KOKKOS_LAMBDA (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });

                // ..................

            }

        }

        {
            BL_PROFILE("Amrex Five");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
 
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });
            }

        }

        {
            BL_PROFILE("Kokkos Streams Three");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);
                int localIndex = mfi.LocalIndex();

                Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>(instances[localIndex%numStreams],
                                                                            {lo.x,lo.y,lo.z},{hi.x+1,hi.y+1,hi.z+1}),
                KOKKOS_LAMBDA (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });

                // ..................

            }

        }

        {
            BL_PROFILE("Amrex Six");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
 
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);
                Array4<Real> b = y.array(mfi);

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    a(i,j,k) += b(i,j,k);
                });
            }

        }



// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    amrex::Print() << "Test Completed." << std::endl;

    } Kokkos::finalize();
    } amrex::Finalize();

}
