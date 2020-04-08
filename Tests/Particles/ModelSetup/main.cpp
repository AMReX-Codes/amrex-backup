#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <PersonContainer.H>
#include <PlaceContainer.H>
#include <AssignPeopleToPlaces.H>

using namespace amrex;

void testModelSetup ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    testModelSetup();
    
    amrex::Finalize();
}

void testModelSetup ()
{
    static_assert(AMREX_SPACEDIM == 2, "This function assumes 2D \n");

    BL_PROFILE("testModelSetup");    

    const int nx = 40;
    const int ny = 32;
    const IntVect max_grid_size(20, 16);

    int is_per[AMREX_SPACEDIM] = {1, 1};    
    RealBox real_box(0.0, 0.0, 1.0, 1.0);
    IntVect lo(0, 0);
    IntVect hi(nx, ny);
    const Box domain(lo, hi);

    Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);

    DistributionMapping dm(ba);

    PersonContainer people(geom, dm, ba);
    people.initRandom();
    
    PlaceContainer places(geom, dm, ba);
    places.initRandom();

    assignPeopleToPlaces(people, places);
}
