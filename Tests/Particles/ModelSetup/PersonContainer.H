#ifndef PERSON_CONTAINER_H_
#define PERSON_CONTAINER_H_

#include <AMReX.H>
#include <AMReX_Particles.H>

struct PersonRealID
{
    enum 
    {
        xorig = 0,
        yorig, 
        ncomps
    };
};

struct PersonIntID
{
    enum 
    {
        place = 0,
        ncomps
    };
};

class PersonIterator
    : public amrex::ParIter<PersonRealID::ncomps, PersonIntID::ncomps>
{

public:

    using amrex::ParIter<PersonRealID::ncomps, PersonIntID::ncomps>::ParIter;    
    using PersonType = amrex::ParIter<PersonRealID::ncomps, PersonIntID::ncomps>::ParticleType; 

    int numPeople () const { return numParticles(); }

    PersonType* data () { return GetArrayOfStructs()().data(); }
    const PersonType* data () const { return GetArrayOfStructs()().data(); }    
};

class PersonContainer
    : public amrex::ParticleContainer<PersonRealID::ncomps, PersonIntID::ncomps>
{

public:

    PersonContainer (const amrex::Geometry& geom,
                     const amrex::DistributionMapping& dm,
                     const amrex::BoxArray& ba)
        : amrex::ParticleContainer<PersonRealID::ncomps, PersonIntID::ncomps>(geom, dm, ba)
        {
            static_assert(AMREX_SPACEDIM == 2, "This class assumes 2D \n");
            do_tiling = true;
            tile_size = amrex::IntVect(1, 1);
        }
    
    void initRandom ();
    
    void jumpToRandomCell ();

    void restoreOriginal ();

    long numUnassigned () const;
};

#endif
