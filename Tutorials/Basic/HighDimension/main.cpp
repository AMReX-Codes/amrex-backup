
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include "fab_test.H"
#include "multifab_test.H"

using namespace amrex;

void main_main ()
{
#if (AMREX_SPACEDIM >= 5)
    test_5d_fab();
#endif
#if (AMREX_SPACEDIM >= 4)
    test_4d_fab();
    test_4d_multifab();
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

