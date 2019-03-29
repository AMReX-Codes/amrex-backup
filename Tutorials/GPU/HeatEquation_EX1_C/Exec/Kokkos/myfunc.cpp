
#include "myfunc.H"
#include "mykernel.H"

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
	      Array<MultiFab, AMREX_SPACEDIM>& flux,
	      Real dt,
              Geometry const& geom)
{

    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    phi_old.FillBoundary(geom.periodicity());

    //
    // Note that this simple example is not optimized.
    // The following two MFIter loops could be merged
    // and we do not have to use flux MultiFab.
    // 
    // =======================================================

    const Real dxinv = geom.InvCellSize(0);
    const Real dyinv = geom.InvCellSize(1);
    const Real dzinv = geom.InvCellSize(2);

    // Compute fluxes one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& ybx = mfi.nodaltilebox(1);
        const Box& zbx = mfi.nodaltilebox(2);
        auto const& fluxx = flux[0].array(mfi);
        auto const& fluxy = flux[1].array(mfi);
        auto const& fluxz = flux[2].array(mfi);
        auto const& phi = phi_old.array(mfi);

        int xcells = xbx.numPts();
        const auto xlo = lbound(xbx);
        const auto xhi = ubound(xbx);

        int ycells = ybx.numPts();
        const auto ylo = lbound(ybx);
        const auto yhi = ubound(ybx);

        int zcells = zbx.numPts();
        const auto zlo = lbound(zbx);
        const auto zhi = ubound(zbx);

        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({xlo.x,xlo.y,xlo.z},{xhi.x+1,xhi.y+1,xhi.z+1}, {32,4,4}),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            compute_flux_x(i,j,k,fluxx,phi,dxinv);
        });

        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({ylo.x,ylo.y,ylo.z},{yhi.x+1,yhi.y+1,yhi.z+1}, {32,4,4}),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            compute_flux_y(i,j,k,fluxy,phi,dyinv);
        });

        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({zlo.x,zlo.y,zlo.z},{zhi.x+1,zhi.y+1,zhi.z+1}, {32,4,4}),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            compute_flux_z(i,j,k,fluxz,phi,dzinv);
        });
    }

    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
        auto const& fluxx = flux[0].array(mfi);
        auto const& fluxy = flux[1].array(mfi);
        auto const& fluxz = flux[2].array(mfi);
        auto const& phiOld = phi_old.array(mfi);
        auto const& phiNew = phi_new.array(mfi);

        int cells = vbx.numPts();
        const auto lo = lbound(vbx);
        const auto hi = ubound(vbx);

        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({lo.x,lo.y,lo.z},{hi.x+1,hi.y+1,hi.z+1}, {32,4,4}),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            update_phi(i,j,k,phiOld,phiNew,fluxx,fluxy,fluxz,dt,dxinv,dyinv,dzinv);
        });
    }
}

void init_phi(amrex::MultiFab& phi_new, amrex::Geometry const& geom)
{

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
    GpuArray<Real,AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();
    // =======================================
    // Initialize phi_new by calling a Fortran routine.
    // MFIter = MultiFab Iterator
    for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phiNew = phi_new.array(mfi);

        int ncells = vbx.numPts();
        const auto lo = lbound(vbx);
        const auto hi = ubound(vbx);

        Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({lo.x,lo.y,lo.z},{hi.x+1,hi.y+1,hi.z+1}, {32,4,4}),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            init_phi(i,j,k,phiNew,dx,prob_lo); 
        });
   
    }

}

