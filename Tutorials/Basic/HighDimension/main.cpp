
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include "f_F.H"

using namespace amrex;

void test_5d_fab ()
{
    Box bx(IntVect(4),IntVect(7),Dimension::five);
    int nc = 1;
    BaseFab<Real> fab(bx,nc);
    amrex::Print() << "Dimension: " << bx.dimension() << ", " << bx << std::endl;

    auto a = fab.arraynd();
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (IntVect const& iv) noexcept
    {
        a(iv) = 5.;
    });

    Gpu::synchronize(); // ParalleFor is async to host

    f5(fab.dataPtr(), fab.loVect(), fab.hiVect(), fab.nComp());
}

void test_4d_fab ()
{
    Box bx(IntVect(AMREX_D6_DECL(0,0,0,0,-5,-6)),
           IntVect(AMREX_D6_DECL(7,15,7,15,-1,-1)),
           Dimension::four);
    int nc = 2;
    BaseFab<Real> fab(bx,nc);
    amrex::Print() << "Dimension: " << bx.dimension() << ", " << bx << std::endl;

    auto a = fab.arraynd();
    amrex::ParallelFor(bx, nc,
    [=] AMREX_GPU_DEVICE (IntVect const& iv, int n) noexcept
    {
        a(iv,n) = 4.;
    });

    Gpu::synchronize(); // ParalleFor is async to host

    f4(fab.dataPtr(), fab.loVect(), fab.hiVect(), fab.nComp());
}

void test_3d_fab ()
{
    Box bx(IntVect(4),IntVect(7),Dimension::three);
    int nc = 1;
    BaseFab<Real> fab(bx,nc);
    amrex::Print() << "Dimension: " << bx.dimension() << ", " << bx << std::endl;

    f3(fab.dataPtr(), fab.loVect(), fab.hiVect());

    auto a = fab.arraynd();
    Gpu::DeviceScalar<long> stot(0);
    long* ptot = stot.dataPtr();
    // For GPU, this is not the most efficient way.
    // It's here for testing only.
    amrex::For(bx,
    [=] AMREX_GPU_DEVICE (IntVect const& iv) noexcept
    {
        Gpu::Atomic::Add(ptot, static_cast<long>(a(iv)));
    });
    long atot = stot.dataValue();  // For is async to host, dataValue() has synced.

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(atot==bx.numPts() && atot==64,
                                     "test_3d failed");
}

void main_main ()
{
#if (AMREX_SPACEDIM >= 5)
    test_5d_fab();
#endif
#if (AMREX_SPACEDIM >= 4)
    test_4d_fab();
#endif
#if (AMREX_SPACEDIM >= 3)
    test_3d_fab();
#endif
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";
    main_main();
    amrex::Finalize();
}

