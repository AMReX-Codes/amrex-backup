#include <PlaceContainer.H>

using namespace amrex;

void PlaceContainer::initRandom ()
{
    static_assert(AMREX_SPACEDIM == 2, "This method assumes 2D \n");

    BL_PROFILE("PlaceContainer::initRandom()");
    
    ParmParse pp;

    long n;
    pp.get("nplaces", n);

    int capacity;
    pp.get("place_capacity", capacity);

    amrex::Print() << "Adding " << n << " places. \n";

    int ibegin, iend;
    {
        int myproc = ParallelDescriptor::MyProc();
        int nprocs = ParallelDescriptor::NProcs();
        int navg = n/nprocs;
        int nleft = n - navg * nprocs;
        if (myproc < nleft) {
            ibegin = myproc*(navg+1);
            iend = ibegin + navg+1;
        } else {
            ibegin = myproc*navg + nleft;
            iend = ibegin + navg;
        }
    }

    auto& ptile = DefineAndReturnParticleTile(0, 0, 0); 
    for (int i = ibegin; i < iend; ++i)
    {
        ParticleType p;
        p.id() = ParticleType::NextID();
        p.cpu() = ParallelDescriptor::MyProc();
        
        p.pos(0) = amrex::Random();
        p.pos(1) = amrex::Random();

        p.idata(PlaceIntID::global_id) = i;
        p.idata(PlaceIntID::capacity) = capacity;
        
        ptile.push_back(p);
    }

    Redistribute();
}

long PlaceContainer::numVacancies () const
{
    BL_PROFILE("PlaceContainer::numVacancies()");

    using PType = typename PlaceContainer::SuperParticleType;
    long num_vacancies = amrex::ReduceSum(*this, 
                                          [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int 
                                          { 
                                              return p.idata(PlaceIntID::capacity); 
                                          });
    amrex::ParallelDescriptor::ReduceLongSum(num_vacancies);        
    return num_vacancies;
}
