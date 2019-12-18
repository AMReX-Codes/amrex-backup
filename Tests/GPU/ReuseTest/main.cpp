#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

// To do:
// ======
// resize going down the wrong path? (Not clearing, just replacing?)
// check NumCallBacks isn't messing anything up.
// Vector of Boxes and indexing scheme to get rid of synchronize.


// Resize callback stuff
// ======================

struct FabResize{
    FArrayBox* fab;
    const Box* bx;
    Array4<Real>* array;
};

void CUDART_CB amrex_farraybox_resize (void* p)
{
    FabResize* fr = (FabResize*) p;

    (fr->fab)->resize(*(fr->bx), AMREX_SPACEDIM);
    *(fr->array) = (fr->fab)->array();
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

    MultiFab mf(ba,DistributionMapping{ba},1,0,MFInfo().SetArena(The_Managed_Arena()));
    mf.setVal(1.5);

    MultiFab mf_elix(ba,DistributionMapping{ba},1,0,MFInfo().SetArena(The_Managed_Arena()));
    mf_elix.setVal(1.5);

    MultiFab mf_reuse(ba,DistributionMapping{ba},1,0,MFInfo().SetArena(The_Managed_Arena()));
    mf_reuse.setVal(1.5);

    if (do_baseline)
    {
        amrex::ResetTotalBytesAllocatedInFabsHWM(); 

        BL_PROFILE("1-Baseline");

        amrex::Print() << "Baseline without temps -- Begin" << std::endl;
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
        amrex::Print() << "Baseline without temps -- End" << std::endl;
        amrex::Print() << "Baseline memory HWM: " << amrex::TotalBytesAllocatedInFabsHWM() << std::endl << std::endl;
    }

    if (do_elixir)
    {
        amrex::ResetTotalBytesAllocatedInFabsHWM(); 

        FArrayBox temp_fab(The_Managed_Arena());

        BL_PROFILE("2-Elixir");

        amrex::Print() << "Elixir temps -- Begin" << std::endl;
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
        amrex::Print() << "Elixir temps -- End" << std::endl;
        amrex::Print() << "Elixir memory HWM: " << amrex::TotalBytesAllocatedInFabsHWM() << std::endl << std::endl;

    }

    if (do_reuse)
    {
        amrex::ResetTotalBytesAllocatedInFabsHWM(); 

        Vector<Box> temp_box( mf_elix.local_size() );
        Vector<Array4<Real>> temp_arr( num_streams ); 
        Vector<FArrayBox*> temp_fab( num_streams );
        for (int i = 0; i<num_streams; ++i)
        {
            temp_fab[i] = new FArrayBox(The_Managed_Arena());
        }

        BL_PROFILE("3-Reuse");
        amrex::Print() << "Reuse temps -- Begin" << std::endl;
        for (MFIter mfi(mf_elix, MFItInfo().EnableTiling(FabArrayBase::mfiter_tile_size).SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf_reuse.array(mfi);
            Array4<Real> temp;

            // Find FArrayBox associated with this stream.
            // Equivalent of setStreamIndex MFIter value in MFIter operator++.
            //    (minus the CPU or no stream variations, for simplicity).
            int stream_id = mfi.tileIndex()%num_streams;
            FArrayBox& my_temp = *(temp_fab[stream_id]);

            // Make these a callback function. No elixir needed. That should be it. :)
            //    ****** my_temp.resize(bx, AMREX_SPACEDIM);
            //    ****** const auto temp = my_temp.array();

            FabResize fr{&my_temp, &bx, &temp};
            cudaLaunchHostFunc(Gpu::gpuStream(), amrex_farraybox_resize, &fr);

            // Currently fails without this, as box and array change outside the callback.
            //    An Array4 can be stored with the matching FArrayBox in a "Temporaries" object and cycled in the callback.
            //    Proper pointers to boxes? std::vector<Box*> for each stream with indexing scheme?
            //    Need to keep FabResize outside?
            amrex::Gpu::Device::synchronize();

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
        amrex::Print() << "Reuse temps -- End" << std::endl;
        amrex::Print() << "Reuse memory HWM: " << amrex::TotalBytesAllocatedInFabsHWM() << std::endl << std::endl;
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
