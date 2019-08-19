
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include "f_F.H"

using namespace amrex;

void test_5d ()
{
    Box bx(IntVect(4),IntVect(7),Dimension::five);
    int nc = 1;
    BaseFab<Real> fab(bx,nc);
    amrex::Print() << "Dimension: " << bx.dimension() << ", " << bx << std::endl;

    auto a = fab.arraynd();
    amrex::ParallelFor(bx,
    [=] (IntVect const& iv) noexcept
    {
        a(iv) = 5.;
    });

    f5(fab.dataPtr(), fab.loVect(), fab.hiVect(), fab.nComp());
}

void test_4d ()
{
    Box bx(IntVect(AMREX_D6_DECL(0,0,0,0,-5,-6)),
           IntVect(AMREX_D6_DECL(7,15,7,15,-1,-1)),
           Dimension::four);
    int nc = 2;
    BaseFab<Real> fab(bx,nc);
    amrex::Print() << "Dimension: " << bx.dimension() << ", " << bx << std::endl;

    auto a = fab.arraynd();
    amrex::ParallelFor(bx, nc,
    [=] (IntVect const& iv, int n) noexcept
    {
        a(iv,n) = 4.;
    });

    f4(fab.dataPtr(), fab.loVect(), fab.hiVect(), fab.nComp());
}

void test_3d ()
{
    Box bx(IntVect(4),IntVect(7),Dimension::three);
    int nc = 1;
    BaseFab<Real> fab(bx,nc);
    amrex::Print() << "Dimension: " << bx.dimension() << ", " << bx << std::endl;

    f3(fab.dataPtr(), fab.loVect(), fab.hiVect());

    auto a = fab.arraynd();
    long atot = 0.0;
    amrex::For(bx,
    [=,&atot] (IntVect const& iv) noexcept
    {
        atot += static_cast<long>(a(iv));
    });

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(atot==bx.numPts() && atot==64,
                                     "test_3d failed");
}

void main_main ()
{
    test_5d();
    test_4d();
    test_3d();
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";
    main_main();
    amrex::Finalize();
}

