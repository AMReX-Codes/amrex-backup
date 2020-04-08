#ifndef PLACE_CONTAINER_H_
#define PLACE_CONTAINER_H_

#include <AMReX.H>
#include <AMReX_Particles.H>

struct PlaceRealID
{
    enum {
        ncomps = 0
    };
};

struct PlaceIntID
{
    enum {
        global_id = 0,
        capacity,
        ncomps
    };
};

class PlaceContainer
    : public amrex::ParticleContainer<PlaceRealID::ncomps, PlaceIntID::ncomps>
{

public:

    PlaceContainer (const amrex::Geometry& geom,
                    const amrex::DistributionMapping& dm,
                    const amrex::BoxArray& ba)
        : amrex::ParticleContainer<PlaceRealID::ncomps, PlaceIntID::ncomps>(geom, dm, ba)
        {
            static_assert(AMREX_SPACEDIM == 2, "This class assumes 2D \n");
            do_tiling = true;
            tile_size = amrex::IntVect(1, 1);
        }
    
    void initRandom ();

    long numVacancies () const;
};

#endif
