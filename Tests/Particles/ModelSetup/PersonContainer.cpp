#include <PersonContainer.H>

using namespace amrex;

void PersonContainer::initRandom ()
{
    static_assert(AMREX_SPACEDIM == 2, "This method assumes 2D \n");

    BL_PROFILE("PersonContainer::initRandom()");
    
    ParmParse pp;

    long n;
    pp.get("npeople", n);

    amrex::Print() << "Adding " << n << " people. \n";

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

        p.rdata(PersonRealID::xorig) = p.pos(0);
        p.rdata(PersonRealID::yorig) = p.pos(1);
        
        p.idata(PersonIntID::place) = -1;

        ptile.push_back(p);
    }

    Redistribute();
}

void PersonContainer::jumpToRandomCell ()
{
    static_assert(AMREX_SPACEDIM == 2, "This method assumes 2D \n");

    BL_PROFILE("PersonContainer::jumpToRandomCell");

    const int lev = 0;

    const Geometry& geom = Geom(lev);
    const auto dx = Geom(lev).CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel
#endif        
    for (PersonIterator pti(*this, lev); pti.isValid(); ++pti)
    {
        auto people = pti.data();
        const auto np = pti.numPeople();

        amrex::ParallelFor( np,
        AMREX_GPU_HOST_DEVICE [=] (const int i)
        {
            auto& p = people[i];
            if (p.idata(PersonIntID::place) >= 0) return;
            p.pos(0) = amrex::Random();
            p.pos(1) = amrex::Random();
        });
    }
}

void PersonContainer::restoreOriginal ()
{
    static_assert(AMREX_SPACEDIM == 2, "This method assumes 2D \n");

    BL_PROFILE("PersonContainer::restoreOriginal");

    const int lev = 0;

    const Geometry& geom = Geom(lev);
    const auto dx = Geom(lev).CellSizeArray();
      
#ifdef _OPENMP
#pragma omp parallel
#endif  
    for (PersonIterator pti(*this, lev); pti.isValid(); ++pti)
    {
        auto people = pti.data();
        const auto np = pti.numPeople();

        amrex::ParallelFor( np,
        AMREX_GPU_HOST_DEVICE [=] (const int i)
        {
            auto& p = people[i];
            p.pos(0) = p.rdata(PersonRealID::xorig);
            p.pos(1) = p.rdata(PersonRealID::yorig);
        });
    }
}

long PersonContainer::numUnassigned () const
{
    BL_PROFILE("PersonContainer::numUnassigned()");

    using PType = typename PersonContainer::SuperParticleType;
    long num_unassigned = amrex::ReduceSum(*this, 
                                           [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int 
                                           { 
                                               return (p.idata(PersonIntID::place) < 0); 
                                           });
    amrex::ParallelDescriptor::ReduceLongSum(num_unassigned);        
    return num_unassigned;
}
