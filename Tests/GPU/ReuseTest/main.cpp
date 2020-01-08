#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

// Try prepping all Array4's for all launches first,
// moving to the device, and then passing the correct index?
// (e.g. fused launch style).

using namespace amrex;

template <class T>
using TrueDeviceVector = PODVector<T, DeviceArenaAllocator<T> >;

// Resize callback stuff
// ======================

struct FabResize{
    const Box bx;
    int ncomp;
    FArrayBox* fab;
    Array4<Real>* array;
};

void CUDART_CB amrex_farraybox_resize (void* p)
{
    FabResize* fr = (FabResize*) p;

    (fr->fab)->resize(fr->bx, fr->ncomp);
    *(fr->array) = (fr->fab)->array();

    delete fr;
}

struct FabResizeP{
    const Box bx;
    int ncomp;
    FArrayBox* fab;
};

void CUDART_CB amrex_farraybox_resize_prepared (void* p)
{
    FabResizeP* fr = (FabResizeP*) p;

    (fr->fab)->resize(fr->bx, fr->ncomp);

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
    int do_reuse_max = 1;
    int do_arrs = 1;
    int do_device = 1;
    int do_error_check = 1;

    int box_size = 128;
    int max_grid_size = 64;
    int num_streams = amrex::Gpu::Device::numGpuStreams();
    std::string file="";
    {
        ParmParse pp;

        pp.query("baseline", do_baseline);
        pp.query("elixir", do_elixir);
        pp.query("reuse", do_reuse);
        pp.query("reuse_max", do_reuse_max);
        pp.query("array4", do_arrs);
        pp.query("reuse_device", do_device);
        pp.query("error", do_error_check);

        pp.query("size", box_size);
        pp.query("max_grid_size", max_grid_size);
        pp.query("num_streams", num_streams);
        pp.query("file", file);
    }

    amrex::Print() << std::endl;
    amrex::Print() << " ******************************************** " << std::endl;
    amrex::Print() << "Building MultiFabs" << std::endl;

    BoxArray ba;
    
    if(file.size() == 0)
    {
        Box domain_box(IntVect(0), IntVect(box_size-1));
        ba.define(domain_box);
        ba.maxSize(max_grid_size);
    }
    else
    {
        Box domain_box;

        std::ifstream ifs;
        ifs.open(file.c_str(), std::ios::in);
        ifs >> domain_box;
        ifs.ignore(1000, '\n');
        ba.readFrom(ifs);
    }

    MultiFab mf(ba,DistributionMapping{ba},1,0);
    mf.setVal(1.5);

    MultiFab mf_elix(ba,DistributionMapping{ba},1,0);
    mf_elix.setVal(1.5);

    MultiFab mf_reuse(ba,DistributionMapping{ba},1,0);
    mf_reuse.setVal(1.5);

    MultiFab mf_max(ba,DistributionMapping{ba},1,0);
    mf_max.setVal(1.5);

    MultiFab mf_arrs(ba,DistributionMapping{ba},1,0);
    mf_arrs.setVal(1.5);

    MultiFab mf_device(ba,DistributionMapping{ba},1,0);
    mf_device.setVal(1.5);

    long FabBytesBegin = amrex::TotalBytesAllocatedInFabs();

// ==================================================================================================

    if (do_baseline)
    {
        amrex::Print() << " ******************************************** " << std::endl;
        amrex::Print() << "Baseline -- start." << std::endl;
        BL_PROFILE("1-Baseline");
        for (MFIter mfi(mf_elix, MFItInfo().SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
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

// ==================================================================================================

    if (do_elixir)
    {
        amrex::Print() << " ******************************************** " << std::endl;
        amrex::Print() << "Elixir -- start." << std::endl;
        amrex::ResetTotalBytesAllocatedInFabsHWM();

        BL_PROFILE("2-Elixir");

        FArrayBox temp_fab;

        for (MFIter mfi(mf_elix, MFItInfo().SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf_elix.array(mfi);

            temp_fab.resize(bx, 1);
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
    Gpu::synchronize();
    amrex::Print() << std::endl << " Synched outside elixir scope. " << std::endl;

// ==================================================================================================

    if (do_reuse)
    {
        amrex::Print() << " ******************************************** " << std::endl;
        amrex::Print() << "Reuse -- start" << std::endl;
        amrex::ResetTotalBytesAllocatedInFabsHWM();

        BL_PROFILE("3-Reuse");

        // Version 1: Work!!! (DONE!) 

        // Version 2: Search for largest box on each stream, initialize FArrayBox with that box.
        //            resize will only ever update the box for the new Array4.
        //            FArrayBox data can be entirely left on the device.

        // Version 3: Make it without Managed Memory.

        // Version 4: API for users.

        Vector<FArrayBox> temp_fab(num_streams);
        Gpu::PinnedVector<Array4<Real>> temp_arr(num_streams);

        for (MFIter mfi(mf_reuse, MFItInfo().SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
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

            FabResize* fr = new FabResize{bx, 1, &my_temp, &temp_arr[stream_id]};
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
    Gpu::synchronize();
    amrex::Print() << std::endl << " Synched outside reuse scope. " << std::endl;
    amrex::Print() << " ******************************************** " << std::endl;

// ==================================================================================================

    if (do_reuse_max)
    {
        amrex::Print() << " ******************************************** " << std::endl;
        amrex::Print() << "Reuse Max -- start" << std::endl;
        amrex::ResetTotalBytesAllocatedInFabsHWM();

        BL_PROFILE("4-Reuse-Max");

        Vector<FArrayBox> temp_fab(num_streams);
        {
            Vector<Box> max_box(num_streams);

            // Find and store maximum box.
            for (MFIter mfi(mf_max, MFItInfo().SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                int sid = mfi.tileIndex()%num_streams;

                if (bx.numPts() > max_box[sid].numPts())
                {
                    max_box[sid] = bx;
                }
            }

            // Initialize temporary fabs to maximum size.
            for (int i=0; i<num_streams; ++i)
            {
                temp_fab[i].resize(max_box[i], 1);
            }
        }

        Gpu::PinnedVector<Array4<Real>> temp_arr(num_streams);

        for (MFIter mfi(mf_max, MFItInfo().SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf_max.array(mfi);

            int stream_id = mfi.tileIndex()%num_streams;
            FArrayBox& my_temp = temp_fab[stream_id];
            Array4<Real>* temp = &(temp_arr[stream_id]);

            FabResize* fr = new FabResize{bx, 1, &my_temp, &temp_arr[stream_id]};
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
        amrex::Print() << "Reuse Max -- memory added: " << amrex::TotalBytesAllocatedInFabsHWM() - FabBytesBegin << std::endl << std::endl;
    }
    Gpu::synchronize();
    amrex::Print() << std::endl << " Synched outside reuse max scope. " << std::endl;
    amrex::Print() << " ******************************************** " << std::endl;

// ==================================================================================================

    if (do_arrs)
    {
        amrex::Print() << " ******************************************** " << std::endl;
        amrex::Print() << "Prepare Array4s -- start" << std::endl;
        amrex::ResetTotalBytesAllocatedInFabsHWM();

        BL_PROFILE("5-Prepare-Array4s");

        int num_fabs = mf_arrs.local_size();

        Vector<FArrayBox> temp_fab(num_streams);
        Gpu::HostVector<Array4<Real>> h_temp_arr(num_fabs);
        Gpu::DeviceVector<Array4<Real>> d_temp_arr(num_fabs);
        {
            Vector<Box> max_box(num_streams);

            // Find and store maximum box.
            for (MFIter mfi(mf_arrs, MFItInfo().SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                int sid = mfi.tileIndex()%num_streams;
                int index = mfi.tileIndex();

                h_temp_arr[index] = amrex::makeArray4(temp_fab[index].dataPtr(), bx, mf_arrs.nComp());

                if (bx.numPts() > max_box[sid].numPts())
                {
                    max_box[sid] = bx;
                }
            }

            amrex::Gpu::htod_memcpy((void*) (d_temp_arr.data()), (void*) (h_temp_arr.data()), sizeof(Array4<Real>)*num_fabs);

            // Initialize temporary fabs to maximum size.
            for (int i=0; i<num_streams; ++i)
            {
                temp_fab[i].resize(max_box[i], 1);
            }
        }

        for (MFIter mfi(mf_arrs, MFItInfo().SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf_arrs.array(mfi);
            Array4<Real>* arrs = d_temp_arr.data();

            int index = mfi.tileIndex();
            int stream_id = index%num_streams;
            FArrayBox& my_temp = temp_fab[stream_id];

            FabResizeP* fr = new FabResizeP{bx, 1, &my_temp};
            cudaLaunchHostFunc(Gpu::gpuStream(), amrex_farraybox_resize_prepared, fr);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {  
                (*(arrs+index))(i,j,k) = fab(i,j,k)*fab(i,j,k);
            });

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {   
                fab(i,j,k) = (*(arrs+index))(i,j,k) + 3.14159;
            });
        }
        amrex::Print() << "Prepare Array4s -- memory added: " << amrex::TotalBytesAllocatedInFabsHWM() - FabBytesBegin << std::endl << std::endl;
    }
    Gpu::synchronize();
    amrex::Print() << std::endl << " Synched outside Prepare Array4s scope. " << std::endl;
    amrex::Print() << " ******************************************** " << std::endl;

// ==================================================================================================

    if (do_device)
    {
        amrex::Print() << " ******************************************** " << std::endl;
        amrex::Print() << "Device -- start" << std::endl;
        amrex::ResetTotalBytesAllocatedInFabsHWM();

        BL_PROFILE("6-Device");

        Vector<FArrayBox> temp_fab(num_streams);
        // Prep with maximum size box.
        {
            Vector<Box> max_box(num_streams);

            for (MFIter mfi(mf_device, MFItInfo().SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                int sid = mfi.tileIndex()%num_streams;

                if (bx.numPts() > max_box[sid].numPts())
                {
                    max_box[sid] = bx;
                }
            }
            for (int i=0; i<num_streams; ++i)
            {
                temp_fab[i].resize(max_box[i], 1);
            }
        }

        Vector<Array4<Real>> h_temp_arr(num_streams);
        TrueDeviceVector<Array4<Real>> d_temp_arr(num_streams);
        for (MFIter mfi(mf_device, MFItInfo().SetNumStreams(num_streams)); mfi.isValid(); ++mfi)
        {
            int stream_id = mfi.tileIndex()%num_streams;
            
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf_device.array(mfi);

            FArrayBox& my_temp = temp_fab[stream_id];
            Array4<Real>* d_temp = &(d_temp_arr[stream_id]);

            FabResize* fr = new FabResize{bx, 1, &my_temp, &(h_temp_arr[stream_id])};
            cudaLaunchHostFunc(Gpu::gpuStream(), amrex_farraybox_resize, fr);
            // Works if this is not asynchronous. FIX!
            amrex::Gpu::htod_memcpy(&(d_temp_arr[stream_id]), &(h_temp_arr[stream_id]), sizeof(Array4<Real>));

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {   
                (*d_temp)(i,j,k) = fab(i,j,k)*fab(i,j,k);
            });

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {   
                fab(i,j,k) = (*d_temp)(i,j,k) + 3.14159;
            });
        }
        amrex::Print() << "Device Reuse -- memory added: " << amrex::TotalBytesAllocatedInFabsHWM() - FabBytesBegin << std::endl << std::endl;
    }

    Gpu::synchronize();
    amrex::Print() << std::endl << " Synched outside device reuse scope. " << std::endl;
    amrex::Print() << " ******************************************** " << std::endl;

// ==================================================================================================

    if (do_error_check)
    {
        BL_PROFILE("99999-Error Check");

        if (do_elixir)
        {
            amrex::Print() << " ******************************************** " << std::endl;
            amrex::Print() << " Elixir error check -- " << std::endl;

            mf_elix.minus(mf, 0, 1, 0);
            amrex::Print() << " **Elixir method error: " << mf_elix.sum() << std::endl;
        }

        if (do_reuse)
        {
            amrex::Print() << " ******************************************** " << std::endl;
            amrex::Print() << " Reuse error check -- " << std::endl;

            mf_reuse.minus(mf, 0, 1, 0);
            amrex::Print() << " **Reuse method error: " << mf_reuse.sum() << std::endl;
        }

        if (do_reuse_max)
        {
            amrex::Print() << " ******************************************** " << std::endl;
            amrex::Print() << " Reuse-max error check -- " << std::endl;

            mf_max.minus(mf, 0, 1, 0);
            amrex::Print() << " **Reuse-max method error: " << mf_max.sum() << std::endl;
        }

        if (do_arrs)
        {
            amrex::Print() << " ******************************************** " << std::endl;
            amrex::Print() << " Reuse-max error check -- " << std::endl;

            mf_arrs.minus(mf, 0, 1, 0);
            amrex::Print() << " **Reuse-max method error: " << mf_arrs.sum() << std::endl;
        }

        if (do_device)
        {
            amrex::Print() << " ******************************************** " << std::endl;
            amrex::Print() << " Reuse-device error check -- " << std::endl;

            mf_device.minus(mf, 0, 1, 0);
            amrex::Print() << " **Reuse-device method error: " << mf_device.sum() << std::endl;
        }
    }

    amrex::Print() << " ******************************************** " << std::endl;
    amrex::Print() << " " << std::endl;
}
