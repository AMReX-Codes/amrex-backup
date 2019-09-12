#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

struct SetValMem
{
    void* ptr = nullptr;
    int ncomp = 0;
    double val = 0.0;
    Dim3 ptr_begin = {0,0,0};
    Dim3 ptr_end = {0,0,0};

    template <class T>
    AMREX_GPU_HOST_DEVICE
    Array4<T> getArray () { return Array4<T>(static_cast<T *>(ptr), ptr_begin, ptr_end, ncomp); }
};

template <typename T>
SetValMem
makeSetValMem (Array4<T> const& dst, Real val)
{
#if __cplusplus < 201402L
    SetValMem mem;
    mem.ptr = (void*)(dst.p);
    mem.ptr_begin = dst.begin;
    mem.ptr_end = dst.end;
    mem.ncomp = ncomp;
    mem.val = val; 
    return mem;

#else

    return SetValMem{ (void*)(dst.p), dst.ncomp, val, dst.begin, dst.end};
#endif
}

// MFIterLoop
// Written as a seperate function for easy changes/testing.
void MFIterLoopFunc(const Box &bx, double val, Array4<Real> &a)
{
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        a(i,j,k) = val;
    });
}

// MFIterLoop
// Written as a seperate function for easy changes/testing.
void MFIterGraphFunc(const Box &bx, SetValMem* svm)
{
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        auto const a = svm->getArray<Real>();
        double val = svm->val;
        a(i,j,k) = val;
    });
}


// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    amrex::Gpu::GraphSafeGuard gpu_gsg(true);
    {
        amrex::Print() << "amrex::Initialize complete." << "\n";

        // AMREX_SPACEDIM: number of dimensions
        int n_cell, max_grid_size, Nghost, Ncomp;
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

            // The domain is broken into boxes of size max_grid_size
            Nghost = 0;
            pp.query("nghost", Nghost);

            Ncomp = 1;
            pp.query("ncomp", Ncomp);
        }

        // make BoxArray and Geometry
        BoxArray ba;
        {
            IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
            IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
            Box domain(dom_lo, dom_hi);

            // Initialize the boxarray "ba" from the single box "bx"
            ba.define(domain);
            // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
            ba.maxSize(max_grid_size);
        }
  
        // How Boxes are distrubuted among MPI processes
        DistributionMapping dm(ba);

        // Malloc value for setval testing.
        Real val = 0.0;

        // Create the MultiFab and touch the data.
        // Ensures the data in on the GPU for all further testing.
        MultiFab x(ba, dm, Ncomp, Nghost);
        x.setVal(val);

        amrex::Print() << "Number of boxes = " << x.size() << std::endl << std::endl;

        Real points = ba.numPts();

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Initial launch to remove any unknown costs in HtoD setup. 

        {
            val = 0.42;
            {
                BL_PROFILE("Initial Cleaning -- IGNORE");
                for (MFIter mfi(x); mfi.isValid(); ++mfi)
                {
                    const Box bx = mfi.validbox();
                    Array4<Real> a = x.array(mfi);

                    MFIterLoopFunc(bx, val, a);
                }
                amrex::Print() << "Test sum = " << x.sum() << ". Expected = " << points*(val) << std::endl;
            }
        }



// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Launch without graphs

        {
            val = 1.0;
            {
                BL_PROFILE("Standard");
                for (MFIter mfi(x); mfi.isValid(); ++mfi)
                {
                    const Box bx = mfi.validbox();
                    Array4<Real> a = x.array(mfi);

                    MFIterLoopFunc(bx, val, a);
                }
            }
            amrex::Print() << "No Graph sum = " << x.sum() << ". Expected = " << points*(val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create one graph per MFIter iteration and execute them.

        {
            BL_PROFILE("IterPerGraph");
            val = 2.0;

            BL_PROFILE_VAR("CREATE: IterPerGraph", cgc2);
            cudaGraph_t     graph[x.local_size()];
            cudaGraphExec_t graphExec[x.local_size()];

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(), cudaStreamCaptureModeGlobal));

                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);

                AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[mfi.LocalIndex()])));
            }

            BL_PROFILE_VAR_STOP(cgc2);
            BL_PROFILE_VAR("INSTANTIATE: IterPerGraph", cgi2);

            for (int i = 0; i<x.local_size(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec[i], graph[i], NULL, NULL, 0));
            }

            BL_PROFILE_VAR_STOP(cgi2);
            BL_PROFILE_VAR("LAUNCH: IterPerGraph", cgl2);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec[mfi.LocalIndex()], amrex::Gpu::Device::cudaStream())); 
            }

            amrex::Gpu::Device::synchronize();
            BL_PROFILE_VAR_STOP(cgl2);

            BL_PROFILE_VAR("DESTROY: IterPerGraph", cgd2);
            for (int i = 0; i<x.local_size(); ++i)
            {
                cudaGraphDestroy(graph[i]);
                cudaGraphExecDestroy(graphExec[i]);
            }
            BL_PROFILE_VAR_STOP(cgd2);

            amrex::Print() << "Graph-per-iter sum = " << x.sum() << ". Expected = " << points*(val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create one graph per CUDA stream and execute them.

        {
            BL_PROFILE("StreamPerGraph");
            val = 3.0;

            BL_PROFILE_VAR("CREATE: StreamPerGraph", cgc3);

            cudaGraph_t     graph[amrex::Gpu::numGpuStreams()];
            cudaGraphExec_t graphExec[amrex::Gpu::numGpuStreams()];

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                if (mfi.LocalIndex() == 0)
                {
                    for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i);
                        AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(), cudaStreamCaptureModeGlobal));
                    }
                    amrex::Gpu::Device::setStreamIndex(mfi.tileIndex());
                } 

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                if (mfi.LocalIndex() == (x.local_size() - 1) )
                { 
                    for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i); 
                        AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[i])));
                    }  
                }
            }

            BL_PROFILE_VAR_STOP(cgc3);
            BL_PROFILE_VAR("INSTANTIATE: StreamPerGraph", cgi3);

            for (int i = 0; i<amrex::Gpu::numGpuStreams(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec[i], graph[i], NULL, NULL, 0));
            }

            BL_PROFILE_VAR_STOP(cgi3);
            BL_PROFILE_VAR("LAUNCH: StreamPerGraph", cgl3);

            for (int i = 0; i<amrex::Gpu::numGpuStreams(); ++i)
            {
                amrex::Gpu::Device::setStreamIndex(i); 
                AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec[i], amrex::Gpu::Device::cudaStream())); 
            }
            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();
            BL_PROFILE_VAR_STOP(cgl3);
            BL_PROFILE_VAR("DESTROY: StreamPerGraph", cgd3); 

            for (int i = 0; i<amrex::Gpu::numGpuStreams(); ++i)
            {
                cudaGraphDestroy(graph[i]);
                cudaGraphExecDestroy(graphExec[i]);
            }
            BL_PROFILE_VAR_STOP(cgd3);

            amrex::Print() << "Graph-per-stream sum = " << x.sum() << ". Expected = " << points*(val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create a single graph for the MFIter loop:
//        an empty node at the start linked to each individually captured stream graph.

        {
            BL_PROFILE("IterGraph");
            val = 4.0;

            BL_PROFILE_VAR("CREATE: IterGraph", cgc4);

            cudaGraph_t     graph[x.local_size()];
            cudaGraphExec_t graphExec;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(), cudaStreamCaptureModeGlobal));

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[mfi.LocalIndex()])));
            }

            cudaGraph_t     graphFull;
            cudaGraphNode_t emptyNode, placeholder;

            AMREX_GPU_SAFE_CALL(cudaGraphCreate(&graphFull, 0));
            AMREX_GPU_SAFE_CALL(cudaGraphAddEmptyNode(&emptyNode, graphFull, &placeholder, 0));
            for (int i=0; i<x.local_size(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphAddChildGraphNode(&placeholder, graphFull, &emptyNode, 1, graph[i]));
            }

            BL_PROFILE_VAR_STOP(cgc4);
            BL_PROFILE_VAR("INSTANTIATE: IterGraph", cgi4);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graphFull, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi4);
            BL_PROFILE_VAR("LAUNCH: IterGraph", cgl4);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 
            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl4);
            BL_PROFILE_VAR("DESTROY: IterGraph", cgd4);

            for (int i = 0; i<x.local_size(); ++i)
            {
                cudaGraphDestroy(graph[i]);
            }
            cudaGraphDestroy(graphFull);
            cudaGraphExecDestroy(graphExec);

            BL_PROFILE_VAR_STOP(cgd4);

            amrex::Print() << "Full-graph-iter sum = " << x.sum() << ". Expected = " << points*(val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create a single graph for the MFIter loop:
//        an empty node at the start linked to each individually captured stream graph.

        {
            BL_PROFILE("StreamGraph");
            val = 5.0;

            BL_PROFILE_VAR("CREATE: StreamGraph", cgc5);

            cudaGraph_t     graph[amrex::Gpu::numGpuStreams()];
            cudaGraphExec_t graphExec;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                if (mfi.LocalIndex() == 0)
                {
                    for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i);
                        AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(), cudaStreamCaptureModeGlobal));
                    }
                    amrex::Gpu::Device::setStreamIndex(mfi.tileIndex());
                } 

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                if (mfi.LocalIndex() == (x.local_size() - 1) )
                { 
                    for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i); 
                        AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[i])));
                    }
                }
            }

            cudaGraph_t     graphFull;
            cudaGraphNode_t emptyNode, placeholder;

            AMREX_GPU_SAFE_CALL(cudaGraphCreate(&graphFull, 0));
            AMREX_GPU_SAFE_CALL(cudaGraphAddEmptyNode(&emptyNode, graphFull, &placeholder, 0));
            for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphAddChildGraphNode(&placeholder, graphFull, &emptyNode, 1, graph[i]));
            }

            BL_PROFILE_VAR_STOP(cgc5);
            BL_PROFILE_VAR("INSTANTIATE: StreamGraph", cgi5);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graphFull, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi5);
            BL_PROFILE_VAR("LAUNCH: StreamGraph", cgl5);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 

            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl5);
            BL_PROFILE_VAR("DESTROY: StreamGraph", cgd5); 

            for (int i=0; i<amrex::Gpu::numGpuStreams(); ++i)
            {
                cudaGraphDestroy(graph[i]);
            }
            cudaGraphDestroy(graphFull);
            cudaGraphExecDestroy(graphExec);

            BL_PROFILE_VAR_STOP(cgd5);

            amrex::Print() << "Linked-graph-stream sum = " << x.sum() << ". Expected = " << points*(val) << std::endl;

        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        {
            BL_PROFILE("Indep - Streamless");
            val = 6.0;

            BL_PROFILE_VAR("CREATE: Indep - Streamless", cgc6);

            cudaGraph_t     graph[x.local_size()];
            cudaGraphExec_t graphExec;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                amrex::Gpu::Device::setStreamIndex(0);
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(), cudaStreamCaptureModeGlobal));

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[mfi.LocalIndex()])));
            }

            cudaGraph_t     graphFull;
            cudaGraphNode_t emptyNode, placeholder;

            AMREX_GPU_SAFE_CALL(cudaGraphCreate(&graphFull, 0));
            AMREX_GPU_SAFE_CALL(cudaGraphAddEmptyNode(&emptyNode, graphFull, &placeholder, 0));
            for (int i=0; i<x.local_size(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphAddChildGraphNode(&placeholder, graphFull, &emptyNode, 1, graph[i]));
            }

            BL_PROFILE_VAR_STOP(cgc6);
            BL_PROFILE_VAR("INSTANTIATE: Indep - Streamless", cgi6);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graphFull, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi6);
            BL_PROFILE_VAR("LAUNCH: Indep - Streamless", cgl6);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 
            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl6);
            BL_PROFILE_VAR("DESTROY: Indep - Streamless", cgd6);

            for (int i=0; i<x.local_size(); ++i)
            {
                cudaGraphDestroy(graph[i]);
            }
            cudaGraphDestroy(graphFull);
            cudaGraphExecDestroy(graphExec);

            BL_PROFILE_VAR_STOP(cgd6);

            amrex::Print() << "Full-graph-independent-streamless sum = " << x.sum() << ". Expected = " << points*(val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        {
            BL_PROFILE("Depend - Streamless");
            val = 7.0;

            BL_PROFILE_VAR("CREATE: Depend - Streamless", cgc7);

            cudaGraph_t     graph;
            cudaGraphExec_t graphExec;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                amrex::Gpu::Device::setStreamIndex(0);
                if (mfi.LocalIndex() == 0)
                {
                    AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(), cudaStreamCaptureModeGlobal));
                } 

                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                if (mfi.LocalIndex() == (x.local_size() - 1) )
                { 
                    amrex::Gpu::Device::setStreamIndex(0); 
                    AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &graph));
                }
            }

            BL_PROFILE_VAR_STOP(cgc7);
            BL_PROFILE_VAR("INSTANTIATE: Depend - Streamless", cgi7);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graph, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi7);
            BL_PROFILE_VAR("LAUNCH: Depend - Streamless", cgl7);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 

            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl7);

            BL_PROFILE_VAR("DESTROY: Depend - Streamless", cgd7);

            cudaGraphDestroy(graph);
            cudaGraphExecDestroy(graphExec);

            BL_PROFILE_VAR_STOP(cgd7);

            amrex::Print() << "Linked-graph-dependent-streamless sum = " << x.sum() << ". Expected = " << points*(val) << std::endl;

        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Stream Graph using cudaEvents

        {
            BL_PROFILE("cudaEvents");
            val = 8.0;

            BL_PROFILE_VAR("CREATE: cudaEvents", cgc8);

            amrex::Gpu::Device::setStreamIndex(0);
            cudaStream_t graph_stream = amrex::Gpu::gpuStream();

            cudaGraph_t     graph;
            cudaGraphExec_t graphExec;
            cudaEvent_t memcpy_event = {0};
            cudaEvent_t rejoin_event = {0}; 

            cudaEventCreate(&memcpy_event, cudaEventDisableTiming);
            cudaEventCreate(&rejoin_event, cudaEventDisableTiming);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                if (mfi.LocalIndex() == 0)
                {
                    cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(), cudaStreamCaptureModeGlobal);
                    cudaEventRecord(memcpy_event, graph_stream);
                    for (int i=1; i<amrex::Gpu::numGpuStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i);
                        cudaStreamWaitEvent(amrex::Gpu::gpuStream(), memcpy_event, 0);
                    }
                    amrex::Gpu::Device::setStreamIndex(0);
                }



                // ..................
                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);
                // ..................

                if (mfi.LocalIndex() == x.local_size()-1)
                {
                    for (int i=1; i<amrex::Gpu::numGpuStreams(); ++i)
                    {
                        amrex::Gpu::Device::setStreamIndex(i);
                        cudaEventRecord(rejoin_event, amrex::Gpu::gpuStream());
                        cudaStreamWaitEvent(graph_stream, rejoin_event, 0);
                    }
                    amrex::Gpu::Device::setStreamIndex(0);
                    cudaStreamEndCapture(graph_stream, &graph);
                }
            }

            BL_PROFILE_VAR_STOP(cgc8);
            BL_PROFILE_VAR("INSTANTIATE: cudaEvents", cgi8);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graph, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi8);
            BL_PROFILE_VAR("LAUNCH: cudaEvents", cgl8);

            amrex::Gpu::Device::setStreamIndex(0); 
            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 

            amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl8);
            BL_PROFILE_VAR("DESTROY: cudaEvents", cgd8); 

            cudaGraphDestroy(graph);
            cudaGraphExecDestroy(graphExec);
            cudaEventDestroy(memcpy_event);
            cudaEventDestroy(rejoin_event);

            BL_PROFILE_VAR_STOP(cgd8);

            amrex::Print() << "cudaGraphs sum = " << x.sum() << ". Expected = " << points*(val) << std::endl;

        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        {
            BL_PROFILE("noTree");
            val = 9.0;

            BL_PROFILE_VAR("CREATE: noTree", cgc9);
            cudaGraph_t     graph[x.local_size()];
            cudaGraphExec_t graphExec;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                AMREX_GPU_SAFE_CALL(cudaStreamBeginCapture(amrex::Gpu::Device::cudaStream(), cudaStreamCaptureModeGlobal));

                const Box bx = mfi.validbox();
                Array4<Real> a = x.array(mfi);

                MFIterLoopFunc(bx, val, a);

                AMREX_GPU_SAFE_CALL(cudaStreamEndCapture(amrex::Gpu::Device::cudaStream(), &(graph[mfi.LocalIndex()])));
            }

            cudaGraph_t     graphFull;
            cudaGraphNode_t emptyNode, placeholder;

            AMREX_GPU_SAFE_CALL(cudaGraphCreate(&graphFull, 0));
            AMREX_GPU_SAFE_CALL(cudaGraphAddEmptyNode(&emptyNode, graphFull, &placeholder, 0));
            for (int i=0; i<x.local_size(); ++i)
            {
                AMREX_GPU_SAFE_CALL(cudaGraphAddChildGraphNode(&placeholder, graphFull, NULL, 0, graph[i]));
            }

            BL_PROFILE_VAR_STOP(cgc9);
            BL_PROFILE_VAR("INSTANTIATE: noTree", cgi9);

            AMREX_GPU_SAFE_CALL(cudaGraphInstantiate(&graphExec, graphFull, NULL, NULL, 0));

            BL_PROFILE_VAR_STOP(cgi9);
            BL_PROFILE_VAR("LAUNCH: noTree", cgl9);

            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, amrex::Gpu::Device::cudaStream())); 
            amrex::Gpu::Device::synchronize();

            BL_PROFILE_VAR_STOP(cgl9);

            BL_PROFILE_VAR("DESTROY: noTree", cgd9);
            for (int i = 0; i<x.local_size(); ++i)
            {
                cudaGraphDestroy(graph[i]);
            }
            cudaGraphDestroy(graphFull);
            cudaGraphExecDestroy(graphExec);
            BL_PROFILE_VAR_STOP(cgd9);

            amrex::Print() << "noTree = " << x.sum() << ". Expected = " << points*(val) << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        {
            BL_PROFILE("amrexObject");
            val = 10.0;

            BL_PROFILE_VAR("CREATE: amrexObject", cgc10);

            amrex::Gpu::Device::setStreamIndex(0);
            cudaStream_t graph_stream = amrex::Gpu::gpuStream();

            cudaGraphExec_t graphExec;
            cudaEvent_t memcpy_event = {0};
            cudaEvent_t rejoin_event = {0}; 
            cudaEventCreate(&memcpy_event, cudaEventDisableTiming);
            cudaEventCreate(&rejoin_event, cudaEventDisableTiming);

            amrex::CudaGraph<SetValMem> cgraph(x.local_size());
 
            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                int idx = mfi.LocalIndex();
                Gpu::Device::startGraphRecording( idx == 0,
                                                  cgraph.getHostPtr(0),
                                                  cgraph.getDevicePtr(0),
                                                  std::size_t(sizeof(SetValMem)*x.local_size()) );

                // ..................
                SetValMem* svm = cgraph.getDevicePtr(idx);
                const Box bx = mfi.validbox();

                MFIterGraphFunc(bx, svm);
                // ..................

                bool last_iter = (idx == x.local_size()-1);
                graphExec = Gpu::Device::stopGraphRecording(last_iter);
                if (last_iter) { cgraph.setGraph(graphExec); }
            }

            BL_PROFILE_VAR_STOP(cgc10);
            BL_PROFILE_VAR("SETUP: amrexObject", cgs10);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                int idx = mfi.LocalIndex();
                cgraph.setParams(idx, makeSetValMem(x[mfi].array(), val));
            }              

            BL_PROFILE_VAR_STOP(cgs10);
            BL_PROFILE_VAR("LAUNCH: amrexObject", cgl10);

            AMREX_GPU_SAFE_CALL(cudaGraphLaunch(graphExec, graph_stream)); 

	    amrex::Gpu::Device::synchronize();
            amrex::Gpu::Device::resetStreamIndex();

            BL_PROFILE_VAR_STOP(cgl10);

            BL_PROFILE_VAR("CudaGraph::~CudaGraph()", cgd10); 

            cudaEventDestroy(memcpy_event);
            cudaEventDestroy(rejoin_event);

            BL_PROFILE_VAR_STOP(cgd10);

            amrex::Print() << "amrexObject sum = " << x.sum() << ". Expected = " << points*(val) << std::endl;

        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        amrex::Print() << "Test Completed." << std::endl;
    }

    amrex::Finalize();
}
