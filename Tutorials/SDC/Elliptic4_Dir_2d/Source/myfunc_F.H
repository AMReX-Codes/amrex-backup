#ifndef MYFUNC_F_H__ 
#define MYFUNC_F_H__ 

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void init_phi(const int* lo, const int* hi,
		amrex_real* data, const int* dlo, const int* dhi,
		const amrex_real* dx, const amrex_real* prob_lo,
		const amrex_real* prob_hi, amrex_real* k_freq);
void analytic_diffusion(const int* lo, const int* hi,
              amrex_real* data, const int* dlo, const int* dhi,
              const amrex_real* dx, const amrex_real* prob_lo,
              const amrex_real* prob_hi, amrex_real* k_freq, amrex_real* epsilon);
void init_beta(const int* lo, const int* hi,
          amrex_real* data, const int* dlo, const int* dhi,
          const amrex_real* dx, const amrex_real* prob_lo,
          const amrex_real* prob_hi, amrex_real* epsilon, amrex_real* k_freq);
  
  void err_phi(const int* lo, const int* hi,
	       amrex_real* data, const int* dlo, const int* dhi,
	       const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi,
	       amrex_real* a,amrex_real* d,amrex_real* r,amrex_real* time,
           amrex_real* epsilon,amrex_real* k_freq, amrex_real* kappa );
    void print_multifab(amrex_real* data, const int* dlo, const int* dhi);

  void SDC_feval_F(const int* lo, const int* hi,
		   const int* domlo, const int* domhi,
		   amrex_real* phi, const int* philo, const int* phihi,
		   amrex_real* fluxx, const int* fxlo, const int* fxhi,
		   amrex_real* fluxy, const int* fylo, const int* fyhi,
#if (BL_SPACEDIM == 3)   
		   amrex_real* fluxz, const int* fzlo, const int* fzhi,
#endif
		   amrex_real* f, const int* flo, const int* fhi,
		   const amrex_real* dx,const amrex_real* a,
		   const amrex_real* d,const amrex_real* r,
           amrex_real* facex, const int* facexlo, const int* facexhi,
           amrex_real* facey, const int* faceylo, const int* faceyhi,
           amrex_real* prodx, const int* prodxlo, const int* prodxhi,
           amrex_real* prody, const int* prodylo, const int* prodyhi,
           int *n, const amrex_real* time, amrex_real* epsilon,amrex_real* k_freq, amrex_real* kappa);
  void SDC_fcomp_reaction_F(const int* lo, const int* hi,
                      const int* domlo, const int* domhi,
                      amrex_real* phi, const int* philo, const int* phihi,
                      amrex_real* rhs, const int* rhslo, const int* rhshi,			    
                      amrex_real* f, const int* flo, const int* fhi,
   		    const amrex_real* dtq, int *n);
  void SDC_Lresid_F(const int* lo, const int* hi,
		    const int* domlo, const int* domhi,
		    amrex_real* phi, const int* philo, const int* phihi,
		    amrex_real* rhs, const int* rhslo, const int* rhshi,
		    amrex_real* res, const int* reslo, const int* reshi,
		    amrex_real* corr, const int* corrlo, const int* corrhi,
		    const amrex_real* dtq,const amrex_real* dx );
void cc_to_face_loc(const int* lo, const int* hi,
                amrex_real* data, const int* dlo, const int* dhi,
                amrex_real* face_data, const int* face_dlo, const int* face_dhi, const int* dir);
void fill_bdry_values(const int* lo, const int* hi,
                      amrex_real* data, const int* dlo, const int* dhi,
                      const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi,
                      amrex_real* time,
                      amrex_real* epsilon,amrex_real* k_freq, amrex_real* kappa );
    
#ifdef __cplusplus
}
#endif

#endif
