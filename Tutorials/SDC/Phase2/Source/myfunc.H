#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include "AMReX_SDCstruct.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>


class Kerrek
  : public MLABecLaplacian 

{ public:
  Kerrek():MLABecLaplacian () {}
  Kerrek(const Vector<Geometry>& a_geom,
	  const Vector<BoxArray>& a_grids,
	  const Vector<DistributionMapping>& a_dmap,
	  const LPInfo& a_info = LPInfo(),
	  const Vector<FabFactory<FArrayBox> const*>& a_factory = {}):MLABecLaplacian (a_geom,a_grids,a_dmap,a_info,a_factory) {}
  virtual ~Kerrek(){}
  
  void fourthOrderBCFill (MultiFab& in, MultiFab& bdry_values)
  {

    BL_PROFILE("MLCellLinOp::fourthOrderBCFill()");
    MLCellLinOp::setLevelBC(0, &bdry_values);
    MLCellLinOp::applyBC (0,0, in, BCMode::Inhomogeneous, StateMode::Solution,
			  m_bndry_sol[0].get(), 0);
    
  }
};
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
		  Kerrek& mlabec,
		  SDCstruct &sdcmats,Real a,Real d,Real r,
          std::array<MultiFab,AMREX_SPACEDIM>& face_bcoef,
          std::array<MultiFab,AMREX_SPACEDIM>& prod_stor,
                  Real time,
                  Real epsilon,
                  Real k_freq,
                  Real kappa,
          MultiFab& bdry_values,
                  int Nprob, 
                  int Lord,
                  int &totV
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
                Real kappa,
                int Nprob,
                int Lord
);
void SDC_fcomp(MultiFab& rhs,
           std::array<MultiFab, AMREX_SPACEDIM>& flux,
	       const Geometry& geom, 
	       const Vector<BCRec>& bc,
	       SDCstruct &sdcmats,
	       MLMG &mlmg,
	       Kerrek& mlabec,
	       Real dt,Real a, Real d,Real r,
           std::array<MultiFab,AMREX_SPACEDIM>& face_bcoef,
           std::array<MultiFab,AMREX_SPACEDIM>& prod_stor,
	       int sdc_m,int npiece, Real time, Real epsilon,
               Real k_freq, Real kappa, int &totV, int Nprob,int Lord);

#endif
