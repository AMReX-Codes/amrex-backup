#ifndef MY_F_F_H_
#define MY_F_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void f5 (amrex::Real* p, const int* lo, const int* hi, int n);

    void f4 (amrex::Real* p, const int* lo, const int* hi, int n);

    void f3 (amrex::Real* p, const int* lo, const int* hi);

#ifdef __cplusplus
}
#endif

#endif
