#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

// To do:
// ======
// resize going down the wrong path? (Not clearing, just replacing?) -- Should be fine. Shouldn't clear unless new size > old size.
// check NumCallBacks isn't messing anything up.
// Vector of Boxes and indexing scheme to get rid of synchronize. -- Do I need box? Should be properly captured by callback and launches.
// Need Array4 on the device. (Managed Array4 for first attempt?)


// Resize callback stuff
// ======================

struct FabResize{
    FArrayBox* fab;
    const Box bx;
    Array4<Real>* array;
};

void CUDART_CB amrex_farraybox_resize (void* p)
{
    FabResize* fr = (FabResize*) p;

    (fr->fab)->resize((fr->bx), AMREX_SPACEDIM);
    *(fr->array) = (fr->fab)->array();

    delete fr;
}

// ======================

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}


void main_main ()
{

    int do_baseline = 1;
    int do_elixir = 1;
    int do_reuse = 0;

    int box_size = 128;
    int max_grid_size = 64;
    int num_streams = amrex::Gpu::Device::numGpuStreams();
    {
        ParmParse pp;

        pp.query("baseline", do_baseline);
        pp.query("elixir", do_elixir);
        pp.query("reuse", do_reuse);
        pp.query("size", box_size);
        pp.query("max_grid_size", max_grid_size);
        pp.query("num_streams", num_streams);
    }

    BoxArray ba;
    {
        Box domain_box(IntVect(0), IntVect(box_size-1));
        ba.define(domain_box);
        ba.maxSize(max_grid_size);
    }

    MultiFab mf(ba,DistributionMapping{ba},1,0);
    mf.setVal(1.5);

    MultiFab mf_elix(ba,DistributionMapping{ba},1,0);
    mf_elix.setVal(1.5);

    MultiFab mf_reuse(ba,DistributionMapping{ba},1,0);
    mf_reuse.setVal(1.5);

    long FabBytesBegin = amrex::TotalBytesAllocatedInFabs();

    if (do_baseline)
    {
        amrex::Print() << "Baseline -- start." << std::endl;
        BL_PROFILE("1-Baseline");
        for (MFIter mfi(mf_elix, MFItInfo().EnableTiling(FabArrayBase::mfiter_tile_size).SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                fab(i,j,k) = fab(i,j,k)*fab(i,j,k) + 3.14159;
            });
        }
        amrex::Print() << "Baseline -- memory added: " << amrex::TotalBytesAllocatedInFabsHWM() - FabBytesBegin << std::endl << std::endl;
    }

    if (do_elixir)
    {
        amrex::Print() << "Elixir -- start." << std::endl;
        amrex::ResetTotalBytesAllocatedInFabsHWM();

        BL_PROFILE("2-Elixir");

        FArrayBox temp_fab;

        for (MFIter mfi(mf_elix, MFItInfo().EnableTiling(FabArrayBase::mfiter_tile_size).SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf_elix.array(mfi);

            temp_fab.resize(bx, AMREX_SPACEDIM);
            Elixir temp_elix = temp_fab.elixir();
            const auto temp = temp_fab.array();

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                temp(i,j,k) = fab(i,j,k)*fab(i,j,k);
            });

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                fab(i,j,k) = temp(i,j,k) + 3.14159;
            });
        }
        amrex::Print() << "Elixir -- memory added: " << amrex::TotalBytesAllocatedInFabsHWM() - FabBytesBegin << std::endl << std::endl;

    }

    if (do_reuse)
    {
        amrex::Print() << "Reuse -- start" << std::endl;
        amrex::ResetTotalBytesAllocatedInFabsHWM();

        BL_PROFILE("3-Reuse");

        // Version 1: Work!!!

        // Version 2: Search for largest box on each stream, initialize FArrayBox with that box.
        //            resize will only ever update the box for the new Array4.
        //            FArrayBox data can be entirely left on the device.

        // Version 3: Make it without Managed Memory.

        // Version 4: API for users.

        Vector<FArrayBox> temp_fab(num_streams);
        Gpu::ManagedVector<Array4<Real>> temp_arr(num_streams);

        for (MFIter mfi(mf_elix, MFItInfo().EnableTiling(FabArrayBase::mfiter_tile_size).SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf_reuse.array(mfi);

            // Find FArrayBox associated with this stream.
            // Equivalent of setStreamIndex MFIter value in MFIter operator++.
            //    (minus the CPU or no stream variations, for simplicity).
            int stream_id = mfi.tileIndex()%num_streams;
            FArrayBox& my_temp = temp_fab[stream_id];
            Array4<Real>* temp = &(temp_arr[stream_id]);

            // Make these a callback function. No elixir needed. That should be it. :)
            //    ****** my_temp.resize(bx, AMREX_SPACEDIM);
            //    ****** const auto temp = my_temp.array();

            FabResize* fr = new FabResize{&my_temp, bx, &temp_arr[stream_id]};
            cudaLaunchHostFunc(Gpu::gpuStream(), amrex_farraybox_resize, fr);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {   
                (*temp)(i,j,k) = fab(i,j,k)*fab(i,j,k);
            });

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {   
                fab(i,j,k) = (*temp)(i,j,k) + 3.14159;
            });
        }
        amrex::Print() << "Reuse -- memory added: " << amrex::TotalBytesAllocatedInFabsHWM() - FabBytesBegin << std::endl << std::endl;
    }

    {
        BL_PROFILE("99999-Error Check");

        if (do_elixir)
        {
            mf_elix.minus(mf, 0, 1, 0);
            amrex::Print() << " **Elixir method error: " << mf_elix.sum() << std::endl;
        }

        if (do_reuse)
        {
            mf_reuse.minus(mf, 0, 1, 0);
            amrex::Print() << " **Reuse method error: " << mf_reuse.sum() << std::endl;
        }
    }

}
