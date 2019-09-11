#include "multifab_test.H"
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

using namespace amrex;

void test_4d_multifab ()
{
    int ncells = 16;
    Box bx(IntVect(0), IntVect(ncells-1),  Dimension::four);
    BoxArray ba(bx);
    ba.maxSize(8);
    DistributionMapping dm(ba);

    Real deltax = 2.0/ncells;
    Real problo = -1.0;

    MultiFab mf(ba,dm,1,1); // 1 component, 1 ghost cells
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box& b = mfi.fabbox();
        ArrayN<Real> const& fab = mf.arraynd(mfi);
        amrex::ParallelFor(b,
        [=] AMREX_GPU_DEVICE (DimN const& iv) noexcept
        {
            // 3D sphere in 4D space. radius is a function of the 4th coordinates
            Real x = problo + (iv.i0+0.5)*deltax;
            Real y = problo + (iv.i1+0.5)*deltax;
            Real z = problo + (iv.i2+0.5)*deltax;
            Real u = problo + (iv.i3+0.5)*deltax;
            Real R = std::abs(1.0 - std::abs(u));
            if ((x*x+y*y+z*z) < R*R) {
                fab(iv) = 1.0;
            } else {
                fab(iv) = 0.0;
            }
        });
    }

    for (int islice = 0; islice < ncells; ++islice) {
        // Generate 3d slices.
        // If this kind operation is useful, we will move it to Src/Base/.
        BoxList bl(Dimension::three);
        Vector<int> pmap;
        Vector<int> index_map; // Given local index in 3d MF, get global index in 4d MF.
        int myproc = ParallelDescriptor::MyProc();
        for (int ibox = 0; ibox < ba.size(); ++ibox) {
            Box b = ba[ibox];
            if (b.smallEnd(3) <= islice and b.bigEnd(3) >= islice) {
                b.removeLastDimension();
                bl.push_back(b);
                pmap.push_back(dm[ibox]);
                if (dm[ibox] == myproc) {
                    index_map.push_back(ibox);
                }
            }
        }
        BoxArray ba3d(std::move(bl));
        DistributionMapping dm3d(std::move(pmap));
        MultiFab mf3d(ba3d,dm3d,1,0);
        for (MFIter mfi(mf3d); mfi.isValid(); ++mfi) {
            const Box& b = mfi.fabbox();
            int idx_4d = index_map[mfi.LocalIndex()];
            ArrayN<Real const> const& fab4d = mf.const_arraynd(idx_4d);
            ArrayN<Real> const& fab3d = mf3d.arraynd(mfi);
            amrex::ParallelFor(b,
            [=] AMREX_GPU_DEVICE (DimN const& iv3d) noexcept
            {
                DimN iv4d = iv3d;
                iv4d.i3 = islice;
                fab3d(iv3d) = fab4d(iv4d);
            });
        }

        VisMF::Write(mf3d, "mf3d-slice-"+std::to_string(islice));
    }
}
