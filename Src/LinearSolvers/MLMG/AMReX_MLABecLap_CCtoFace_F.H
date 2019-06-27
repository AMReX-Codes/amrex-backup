#ifndef AMReX_MLABecLap_CCtoFace_F_H_
#define AMReX_MLABecLap_CCtoFace_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void cc_to_face(const int* lo, const int* hi,
		amrex_real* data, const int* dlo, const int* dhi, 
                amrex_real* face_data, const int* face_dlo, const int* face_dhi, int dir);
  
    
#ifdef __cplusplus
}
#endif

#endif
