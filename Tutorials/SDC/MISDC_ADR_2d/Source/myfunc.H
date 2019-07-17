#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include "AMReX_SDCstruct.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

// #include <AMReX_LO_BCTYPES.H>

using namespace amrex;

void main_main ();

void SDC_advance (MultiFab& phi_old,
		  MultiFab& phi_new,
		  std::array<MultiFab, AMREX_SPACEDIM>& flux,
		  Real dt,
		  const Geometry& geom, 
		  const Vector<BCRec>& bc,
        	  MLMG&  mlmg,
		  MLABecLaplacian& mlabec,
		  SDCstruct &sdcmats,Real a,Real d,Real r,
          std::array<MultiFab,AMREX_SPACEDIM>& face_bcoef,
          std::array<MultiFab,AMREX_SPACEDIM>& prod_stor,
                  Real time,
                  Real epsilon,
                  Real k_freq,
                  Real kappa
);
void SDC_feval (std::array<MultiFab, AMREX_SPACEDIM>& flux,
		const Geometry& geom, 
		const Vector<BCRec>& bc,
		SDCstruct &sdcmats,
		Real a,Real d,Real r,
        std::array<MultiFab,AMREX_SPACEDIM>& face_bcoef,
        std::array<MultiFab,AMREX_SPACEDIM>& prod_stor,
		int sdc_m,int npiece, Real time,
                Real epsilon,
                Real k_freq,
                Real kappa
);
void SDC_fcomp(MultiFab& rhs,
           std::array<MultiFab, AMREX_SPACEDIM>& flux,
	       const Geometry& geom, 
	       const Vector<BCRec>& bc,
	       SDCstruct &sdcmats,
	       MLMG &mlmg,
	       MLABecLaplacian& mlabec,
	       Real dt,Real a, Real d,Real r,
           std::array<MultiFab,AMREX_SPACEDIM>& face_bcoef,
           std::array<MultiFab,AMREX_SPACEDIM>& prod_stor,
	       int sdc_m,int npiece);

#endif
