#include "multifab_test.H"
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void test_4d_multifab ()
{
    Box bx(IntVect(0), IntVect(15),  Dimension::four);
    BoxArray ba(bx);
    ba.maxSize(8);
    DistributionMapping dm(ba);

    MultiFab mf(ba,dm,1,1); // 1 component, 1 ghost cells
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        const Box& vbx = mfi.validbox();
        const Box& gbx = mfi.fabbox();
        FArrayBox& fab = mf[mfi];
        const Box& fbx = fab.box();
        amrex::Print() << "  Index and Boxes: " << mfi.index() << ", "
                       << tbx << ", " << vbx << ", " << gbx << ", " << fbx << std::endl;
    }
}
