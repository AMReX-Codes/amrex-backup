#ifndef AMREX_EB_MULTIFAB_UTIL_3D_C_H_
#define AMREX_EB_MULTIFAB_UTIL_3D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
EB_interp_2D_cell_centroid_to_face_centroid(Real x_in,Real y_in,Real phi0,Real phi1,Real phi2,Real phi3,
                                            Real x0_in,Real y0_in, Real x1_in,Real y1_in,
                                            Real x2_in,Real y2_in, Real x3_in,Real y3_in)
{
          Real x  = x_in  - x0_in;
          Real x1 = x1_in - x0_in;
          Real x2 = x2_in - x0_in;
          Real x3 = x3_in - x0_in;
          Real x4 = x2 - x1 - x3;

          Real y  = y_in  - y0_in;
          Real y1 = y1_in - y0_in;
          Real y2 = y2_in - y0_in;
          Real y3 = y3_in - y0_in;
          Real y4 = y2 - y1 - y3;

          std::cout << " IN: x  y  " << x_in  << " " << y_in << std::endl;
          std::cout << " IN: x0 y0 " << x0_in << " " << y0_in << std::endl;
          std::cout << " IN: x1 y1 " << x1_in << " " << y1_in << std::endl;
          std::cout << " IN: x2 y2 " << x2_in << " " << y2_in << std::endl;
          std::cout << " IN: x3 y3 " << x3_in << " " << y3_in << std::endl;

          std::cout << " CALC x y  " << x  << " " << y << std::endl;
          std::cout << " IN: x1 y1 " << x1 << " " << y1 << std::endl;
          std::cout << " IN: x2 y2 " << x2 << " " << y2 << std::endl;
          std::cout << " IN: x3 y3 " << x3 << " " << y3 << std::endl;
          std::cout << " IN: x4 y4 " << x4 << " " << y4 << std::endl;

          Real xi_p   = x*y4 - x1*y3 + x3*y1 - x4*y + sqrt( ( x*y4 - x1*y3 + x3*y1 - x4*y)*( x*y4 - x1*y3 + x3*y1 - x4*y)
                                                           -4.*(x4*y1 - x1*y4) * (x*y3 - x3*y) );
               xi_p   = xi_p / 2. / (x1*y4 - x4*y1);

          Real xi_m   = x*y4 - x1*y3 + x3*y1 - x4*y - sqrt( ( x*y4 - x1*y3 + x3*y1 - x4*y)*( x*y4 - x1*y3 + x3*y1 - x4*y) 
                                                           -4.*(x4*y1 - x1*y4) * (x*y3 - x3*y) );
               xi_m   = xi_m / 2. / (x1*y4 - x4*y1);

          Real eta_p  = x*y4 + x1*y3 - x3*y1 - x4*y + sqrt( (-x*y4 - x1*y3 + x3*y1 + x4*y)*(-x*y4 - x1*y3 + x3*y1 + x4*y) 
                                                           -4.*(x3*y4 - x4*y3) * (x1*y - x*y1) );
               eta_p  = eta_p / 2. / (x3*y4 - x4*y3);

          Real eta_m  = x*y4 + x1*y3 - x3*y1 - x4*y - sqrt( (-x*y4 - x1*y3 + x3*y1 + x4*y)*(-x*y4 - x1*y3 + x3*y1 + x4*y) 
                                                           -4.*(x3*y4 - x4*y3) * (x1*y - x*y1) );
               eta_m  = eta_m / 2. / (x3*y4 - x4*y3);

          Real xi, eta;

          if (xi_m >= 0. and xi_m <= 1.0) 
             xi = xi_m;
          else if (xi_p >= 0. and xi_p <= 1.0) 
             xi = xi_p;
          else 
            amrex::Abort("Bad xi !");

          if (eta_m >= 0. and eta_m <= 1.0) 
             eta = eta_m;
          else if (eta_p >= 0. and eta_p <= 1.0) 
             eta = eta_p;
          else 
            amrex::Abort("Bad eta !");

          Real phi_out = (1. - xi - eta) * phi0 + xi*(1.0 - eta) * phi1 + (xi*eta) * phi2 + eta*(1.0 - xi) *phi3;
          return phi_out;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_set_covered_nodes (int i, int j, int k, int n, int icomp, Array4<Real> const& d,
                           Array4<EBCellFlag const> const& f, Real v)
{
    if (f(i-1,j-1,k-1).isCovered() and f(i  ,j-1,k-1).isCovered() and
        f(i-1,j  ,k-1).isCovered() and f(i  ,j  ,k-1).isCovered() and
        f(i-1,j-1,k  ).isCovered() and f(i  ,j-1,k  ).isCovered() and
        f(i-1,j  ,k  ).isCovered() and f(i  ,j  ,k  ).isCovered())
    {
        d(i,j,k,n+icomp) = v;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_set_covered_nodes (int i, int j, int k, int n, int icomp, Array4<Real> const& d,
                           Array4<EBCellFlag const> const& f, Real const * AMREX_RESTRICT v)
{
    if (f(i-1,j-1,k-1).isCovered() and f(i  ,j-1,k-1).isCovered() and
        f(i-1,j  ,k-1).isCovered() and f(i  ,j  ,k-1).isCovered() and
        f(i-1,j-1,k  ).isCovered() and f(i  ,j-1,k  ).isCovered() and
        f(i-1,j  ,k  ).isCovered() and f(i  ,j  ,k  ).isCovered())
    {
        d(i,j,k,n+icomp) = v[n];
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_with_vol (int i, int j, int k,
                          Array4<Real const> const& fine, int fcomp,
                          Array4<Real> const& crse, int ccomp,
                          Array4<Real const> const& fv, Array4<Real const> const& vfrc,
                          Dim3 const& ratio, int ncomp)
{
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real cv = 0.0;
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = fv(ii,jj,kk)*vfrc(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}}
        if (cv > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/cv;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,j*ratio.y,k*ratio.z,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown (int i, int j, int k,
                 Array4<Real const> const& fine, int fcomp,
                 Array4<Real> const& crse, int ccomp,
                 Array4<Real const> const& vfrc,
                 Dim3 const& ratio, int ncomp)
{
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real cv = 0.0;
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = vfrc(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}}
        if (cv > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/cv;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,j*ratio.y,k*ratio.z,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_face_x (int i, int j, int k,
                        Array4<Real const> const& fine, int fcomp,
                        Array4<Real> const& crse, int ccomp,
                        Array4<Real const> const& area,
                        Dim3 const& ratio, int ncomp)
{
    int ii = i*ratio.x;
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real a = 0.0;
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }}
        if (a > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/a;
        } else {
            crse(i,j,k,n+ccomp) = fine(ii,j*ratio.y,k*ratio.z,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_face_y (int i, int j, int k,
                        Array4<Real const> const& fine, int fcomp,
                        Array4<Real> const& crse, int ccomp,
                        Array4<Real const> const& area,
                        Dim3 const& ratio, int ncomp)
{
    int jj = j*ratio.y;
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real a = 0.0;
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }}
        if (a > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/a;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,jj,k*ratio.z,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_face_z (int i, int j, int k,
                        Array4<Real const> const& fine, int fcomp,
                        Array4<Real> const& crse, int ccomp,
                        Array4<Real const> const& area,
                        Dim3 const& ratio, int ncomp)
{
    int kk = k*ratio.z;
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real a = 0.0;
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = area(ii,jj,kk);
            c += tmp*fine(ii,jj,kk,n+fcomp);
            a += tmp;
        }}
        if (a > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/a;
        } else {
            crse(i,j,k,n+ccomp) = fine(i*ratio.x,j*ratio.y,kk,n+fcomp);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avgdown_boundaries (int i, int j, int k,
                            Array4<Real const> const& fine, int fcomp,
                            Array4<Real> const& crse, int ccomp,
                            Array4<Real const> const& ba,
                            Dim3 const& ratio, int ncomp)
{
    for (int n = 0; n < ncomp; ++n) {
        Real c = 0.0;
        Real cv = 0.0;
        for (int kk = k*ratio.z; kk < (k+1)*ratio.z; ++kk) {
        for (int jj = j*ratio.y; jj < (j+1)*ratio.y; ++jj) {
        for (int ii = i*ratio.x; ii < (i+1)*ratio.x; ++ii) {
            Real tmp = ba(ii,jj,kk);
            c += fine(ii,jj,kk,n+fcomp)*tmp;
            cv += tmp;
        }}}
        if (cv > 1.e-30) {
            crse(i,j,k,n+ccomp) = c/cv;
        } else {
            crse(i,j,k,n+ccomp) = 0.0;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_compute_divergence (int i, int j, int k, int n, Array4<Real> const& divu,
                            Array4<Real const> const& u, Array4<Real const> const& v,
                            Array4<Real const> const& w, Array4<int const> const& ccm,
                            Array4<EBCellFlag const> const& flag, Array4<Real const> const& vfrc,
                            Array4<Real const> const& apx, Array4<Real const> const& apy,
                            Array4<Real const> const& apz, Array4<Real const> const& fcx,
                            Array4<Real const> const& fcy, Array4<Real const> const& fcz,
                            GpuArray<Real,3> const& dxinv, bool already_on_centroids)
{
    if (flag(i,j,k).isCovered())
    {
        divu(i,j,k,n) = 0.0;
    }
    else if (flag(i,j,k).isRegular())
    {
        divu(i,j,k,n) = dxinv[0] * (u(i+1,j,k,n)-u(i,j,k,n))
            +           dxinv[1] * (v(i,j+1,k,n)-v(i,j,k,n))
            +           dxinv[2] * (w(i,j,k+1,n)-w(i,j,k,n));
    }
    else if (already_on_centroids)
    {
        divu(i,j,k,n) = (1.0/vfrc(i,j,k)) * (
                        dxinv[0] * (apx(i+1,j,k)*u(i+1,j,k,n)-apx(i,j,k)*u(i,j,k,n))
            +           dxinv[1] * (apy(i,j+1,k)*v(i,j+1,k,n)-apy(i,j,k)*v(i,j,k,n))
            +           dxinv[2] * (apz(i,j,k+1)*w(i,j,k+1,n)-apz(i,j,k)*w(i,j,k,n)) );
    }
    else
    {
        Real fxm = u(i,j,k,n);
        if (apx(i,j,k) != 0.0 and apx(i,j,k) != 1.0) {
            int jj = j + static_cast<int>(amrex::Math::copysign(1.0, fcx(i,j,k,0)));
            int kk = k + static_cast<int>(amrex::Math::copysign(1.0, fcx(i,j,k,1)));
            Real fracy = (ccm(i-1,jj,k) || ccm(i,jj,k)) ? amrex::Math::abs(fcx(i,j,k,0)) : 0.0;
            Real fracz = (ccm(i-1,j,kk) || ccm(i,j,kk)) ? amrex::Math::abs(fcx(i,j,k,1)) : 0.0;
            fxm = (1.0-fracy)*(1.0-fracz)*fxm
                +      fracy *(1.0-fracz)*u(i,jj,k ,n)
                +      fracz *(1.0-fracy)*u(i,j ,kk,n)
                +      fracy *     fracz *u(i,jj,kk,n);
        }

        Real fxp = u(i+1,j,k,n);
        if (apx(i+1,j,k) != 0.0 and apx(i+1,j,k) != 1.0) {
            int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcx(i+1,j,k,0)));
            int kk = k + static_cast<int>(amrex::Math::copysign(1.0,fcx(i+1,j,k,1)));
            Real fracy = (ccm(i,jj,k) || ccm(i+1,jj,k)) ? amrex::Math::abs(fcx(i+1,j,k,0)) : 0.0;
            Real fracz = (ccm(i,j,kk) || ccm(i+1,j,kk)) ? amrex::Math::abs(fcx(i+1,j,k,1)) : 0.0;
            fxp = (1.0-fracy)*(1.0-fracz)*fxp
                +      fracy *(1.0-fracz)*u(i+1,jj,k ,n)
                +      fracz *(1.0-fracy)*u(i+1,j ,kk,n)
                +      fracy *     fracz *u(i+1,jj,kk,n);
        }

        Real fym = v(i,j,k,n);
        if (apy(i,j,k) != 0.0 and apy(i,j,k) != 1.0) {
            int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j,k,0)));
            int kk = k + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j,k,1)));
            Real fracx = (ccm(ii,j-1,k) || ccm(ii,j,k)) ? amrex::Math::abs(fcy(i,j,k,0)) : 0.0;
            Real fracz = (ccm(i,j-1,kk) || ccm(i,j,kk)) ? amrex::Math::abs(fcy(i,j,k,1)) : 0.0;
            fym = (1.0-fracx)*(1.0-fracz)*fym
                +      fracx *(1.0-fracz)*v(ii,j,k ,n)
                +      fracz *(1.0-fracx)*v(i ,j,kk,n)
                +      fracx *     fracz *v(ii,j,kk,n);
        }

        Real fyp = v(i,j+1,k,n);
        if (apy(i,j+1,k) != 0.0 and apy(i,j+1,k) != 1.0) {
            int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j+1,k,0)));
            int kk = k + static_cast<int>(amrex::Math::copysign(1.0,fcy(i,j+1,k,1)));
            Real fracx = (ccm(ii,j,k) || ccm(ii,j+1,k)) ? amrex::Math::abs(fcy(i,j+1,k,0)) : 0.0;
            Real fracz = (ccm(i,j,kk) || ccm(i,j+1,kk)) ? amrex::Math::abs(fcy(i,j+1,k,1)) : 0.0;
            fyp = (1.0-fracx)*(1.0-fracz)*fyp
                +      fracx *(1.0-fracz)*v(ii,j+1,k ,n)
                +      fracz *(1.0-fracx)*v(i ,j+1,kk,n)
                +      fracx *     fracz *v(ii,j+1,kk,n);
        }

        Real fzm = w(i,j,k,n);
        if (apz(i,j,k) != 0.0 and apz(i,j,k) != 1.0) {
            int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcz(i,j,k,0)));
            int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcz(i,j,k,1)));
            Real fracx = (ccm(ii,j,k-1) || ccm(ii,j,k)) ? amrex::Math::abs(fcz(i,j,k,0)) : 0.0;
            Real fracy = (ccm(i,jj,k-1) || ccm(i,jj,k)) ? amrex::Math::abs(fcz(i,j,k,1)) : 0.0;
            fzm = (1.0-fracx)*(1.0-fracy)*fzm
                +      fracx *(1.0-fracy)*w(ii,j ,k,n)
                +      fracy *(1.0-fracx)*w(i ,jj,k,n)
                +      fracx *     fracy *w(ii,jj,k,n);
        }

        Real fzp = w(i,j,k+1,n);
        if (apz(i,j,k+1) != 0.0 and apz(i,j,k+1) != 1.0) {
            int ii = i + static_cast<int>(amrex::Math::copysign(1.0,fcz(i,j,k+1,0)));
            int jj = j + static_cast<int>(amrex::Math::copysign(1.0,fcz(i,j,k+1,1)));
            Real fracx = (ccm(ii,j,k) || ccm(ii,j,k+1)) ? amrex::Math::abs(fcz(i,j,k+1,0)) : 0.0;
            Real fracy = (ccm(i,jj,k) || ccm(i,jj,k+1)) ? amrex::Math::abs(fcz(i,j,k+1,1)) : 0.0;
            fzp = (1.0-fracx)*(1.0-fracy)*fzp
                +      fracx *(1.0-fracy)*w(ii,j ,k+1,n)
                +      fracy *(1.0-fracx)*w(i ,jj,k+1,n)
                +      fracx *     fracy *w(ii,jj,k+1,n);
        }

        divu(i,j,k,n) = (1.0/vfrc(i,j,k)) *
            ( dxinv[0] * (apx(i+1,j,k)*fxp-apx(i,j,k)*fxm)
            + dxinv[1] * (apy(i,j+1,k)*fyp-apy(i,j,k)*fym)
            + dxinv[2] * (apz(i,j,k+1)*fzp-apz(i,j,k)*fzm) );
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_avg_fc_to_cc (int i, int j, int k, int n, Array4<Real> const& cc,
                      Array4<Real const> const& fx, Array4<Real const> const& fy,
                      Array4<Real const> const& fz, Array4<Real const> const& ax,
                      Array4<Real const> const& ay, Array4<Real const> const& az,
                      Array4<EBCellFlag const> const& flag)
{
    if (flag(i,j,k).isCovered()) {
        cc(i,j,k,n+0) = 0.0;
        cc(i,j,k,n+1) = 0.0;
        cc(i,j,k,n+2) = 0.0;
    } else {
        if (ax(i,j,k) == 0.0) {
            cc(i,j,k,n+0) = fx(i+1,j,k);
        } else if (ax(i+1,j,k) == 0.0) {
            cc(i,j,k,n+0) = fx(i,j,k);
        } else {
            cc(i,j,k,n+0) = 0.5 * (fx(i,j,k) + fx(i+1,j,k));
        }

        if (ay(i,j,k) == 0.0) {
            cc(i,j,k,n+1) = fy(i,j+1,k);
        } else if (ay(i,j+1,k) == 0.0) {
            cc(i,j,k,n+1) = fy(i,j,k);
        } else {
            cc(i,j,k,n+1) = 0.5 * (fy(i,j,k) + fy(i,j+1,k));
        }

        if (az(i,j,k) == 0.0) {
            cc(i,j,k,n+2) = fz(i,j,k+1);
        } else if (az(i,j,k+1) == 0.0) {
            cc(i,j,k,n+2) = fz(i,j,k);
        } else {
            cc(i,j,k,n+2) = 0.5 * (fz(i,j,k) + fz(i,j,k+1));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2cent (Box const& box,
                        const Array4<Real>& phicent,
                        Array4<Real const > const& phicc,
                        Array4<EBCellFlag const> const& flag,
                        Array4<Real const> const& cent,
                        int ncomp) noexcept
{
  amrex::Loop(box, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (flag(i,j,k).isCovered())
    {
      phicent(i,j,k,n) = phicc(i,j,k,n); //do nothing
    }
    else
    {
      if (flag(i,j,k).isRegular())
      {
        phicent(i,j,k,n) = phicc(i,j,k,n);
      }
      else
      {                      
        Real gx = cent(i,j,k,0);
        Real gy = cent(i,j,k,1);
        Real gz = cent(i,j,k,2);
        int ii = (gx < 0.0) ? i - 1 : i + 1;
        int jj = (gy < 0.0) ? j - 1 : j + 1;
        int kk = (gz < 0.0) ? k - 1 : k + 1;
        gx = amrex::Math::abs(gx);
        gy = amrex::Math::abs(gy);
        gz = amrex::Math::abs(gz);
        Real gxy = gx*gy;
        Real gxz = gx*gz;
        Real gyz = gy*gz;
        Real gxyz = gx*gy*gz;
        phicent(i,j,k,n)
          = ( 1.0 - gx - gy - gz + gxy + gxz + gyz - gxyz) * phicc(i ,j ,k ,n)
          + (                 gz       - gxz - gyz + gxyz) * phicc(i ,j ,kk,n)
          + (            gy      - gxy       - gyz + gxyz) * phicc(i ,jj,k ,n)
          + (                                  gyz - gxyz) * phicc(i ,jj,kk,n)
          + (       gx           - gxy - gxz       + gxyz) * phicc(ii,j ,k ,n)
          + (                            gxz       - gxyz) * phicc(ii,j ,kk,n)
          + (                      gxy             - gxyz) * phicc(ii,jj,k ,n)
          + (                                        gxyz) * phicc(ii,jj,kk,n);
      }
    }
  });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2facecent (Box const& box,
                            Array4<Real> const& phi,
                            Array4<EBCellFlag const> const& flag,
                            Array4<Real const> const& apx,
                            Array4<Real const> const& apy,
                            Array4<Real const> const& apz,
                            Array4<Real const> const& fcx,
                            Array4<Real const> const& fcy,
                            Array4<Real const> const& fcz,
                            Array4<Real> const& edg_x,
                            Array4<Real> const& edg_y,
                            Array4<Real> const& edg_z,
                            int ncomp,
                            const Box& domain,
                            const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);

  const Box ubx = amrex::surroundingNodes(box,0);
  const Box vbx = amrex::surroundingNodes(box,1);
  const Box wbx = amrex::surroundingNodes(box,2);

  //
  // ===================== X =====================
  //
  amrex::Loop(ubx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (apx(i,j,k) == 0)
    {
      edg_x(i,j,k,n) = 1e40;
    }
    else
    {
      if (apx(i,j,k) == 1)
      {
        if ( (i == domlo.x) and (bc[n].lo(0) == BCType::ext_dir) )
        {
          edg_x(i,j,k,n) = phi(domlo.x-1,j,k,n);
        }
        else if ( (i == domhi.x+1) and (bc[n].hi(0) == BCType::ext_dir) )
        {
          edg_x(i,j,k,n) = phi(domhi.x+1,j,k,n);
        }
        else
        {
          edg_x(i,j,k,n) = 0.5 * ( phi(i,j,k,n) + phi(i-1,j,k,n) );
        }          
      }
      else
      {
        Real gx = 0.5; 
        Real gy = fcx(i,j,k,0);
        Real gz = fcx(i,j,k,1);
        int ii = i - 1;
        int jj = (gy < 0.0) ? j - 1 : j + 1;
        int kk = (gz < 0.0) ? k - 1 : k + 1;
        gy = amrex::Math::abs(gy);
        gz = amrex::Math::abs(gz);
        Real gxy = gx*gy;
        Real gxz = gx*gz;
        Real gyz = gy*gz;
        Real gxyz = gx*gy*gz;
        edg_x(i,j,k,n) =
          (   1.0 - gx - gy - gz + gxy + gxz + gyz - gxyz) * phi(i ,j ,k ,n)
          + (                 gz       - gxz - gyz + gxyz) * phi(i ,j ,kk,n)
          + (            gy      - gxy       - gyz + gxyz) * phi(i ,jj,k ,n)
          + (                                  gyz - gxyz) * phi(i ,jj,kk,n)
          + (       gx           - gxy - gxz       + gxyz) * phi(ii,j ,k ,n)
          + (                            gxz       - gxyz) * phi(ii,j ,kk,n)
          + (                      gxy             - gxyz) * phi(ii,jj,k ,n)
          + (                                        gxyz) * phi(ii,jj,kk,n);
      } 
    }
  });

  //
  // ===================== Y =====================
  //
  amrex::Loop(vbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (apy(i,j,k) == 0)
    {
      edg_y(i,j,k,n) = 1e40;
    }
    else
    {
      if (apy(i,j,k) == 1)
      {
        if ( (j == domlo.y) and (bc[n].lo(1) == BCType::ext_dir) )
        {
          edg_y(i,j,k,n) = phi(i,domlo.y-1,k,n);
        }
        else if ( (j == domhi.y+1) and (bc[n].hi(1) == BCType::ext_dir) )
        {
          edg_y(i,j,k,n) = phi(i,domhi.y+1,k,n);
        }
        else
        {
          edg_y(i,j,k,n) = 0.5 * (phi(i,j  ,k,n) + phi(i,j-1,k,n));
        }
      }
      else
      {    
        Real gx = fcy(i,j,k,0);  
        Real gy = 0.5;
        Real gz = fcy(i,j,k,1);
        int ii = (gx < 0.0) ? i - 1 : i + 1;
        int jj = j - 1;
        int kk = (gz < 0.0) ? k - 1 : k + 1;
        gx = amrex::Math::abs(gx);
        gz = amrex::Math::abs(gz);
        Real gxy = gx*gy;
        Real gxz = gx*gz;
        Real gyz = gy*gz;
        Real gxyz = gx*gy*gz;
        edg_y(i,j,k,n) =
          (   1.0 - gx - gy - gz + gxy + gxz + gyz - gxyz) * phi(i ,j ,k ,n)
          + (                 gz       - gxz - gyz + gxyz) * phi(i ,j ,kk,n)
          + (            gy      - gxy       - gyz + gxyz) * phi(i ,jj,k ,n)
          + (                                  gyz - gxyz) * phi(i ,jj,kk,n)
          + (       gx           - gxy - gxz       + gxyz) * phi(ii,j ,k ,n)
          + (                            gxz       - gxyz) * phi(ii,j ,kk,n)
          + (                      gxy             - gxyz) * phi(ii,jj,k ,n)
          + (                                        gxyz) * phi(ii,jj,kk,n);
      } 
    }
  });

  //
  // ===================== Z =====================
  //
  amrex::Loop(wbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if (apz(i,j,k) == 0)
    {
      edg_z(i,j,k,n) = 1e40;
    }
    else
    {
      if (apz(i,j,k) == 1)
      {
        if ( (k == domlo.z) and (bc[n].lo(2) == BCType::ext_dir) )
        {
          edg_z(i,j,k,n) = phi(i,j,domlo.z-1,n);
        }
        else if ( (k == domhi.z+1) and (bc[n].hi(2) == BCType::ext_dir) )
        {
          edg_z(i,j,k,n) = phi(i,j,domhi.z+1,n);
        }
        else
        {
          edg_z(i,j,k,n) = 0.5 * (phi(i,j  ,k,n) + phi(i,j,k-1,n));
        }
      }
      else
      {
        Real gx = fcz(i,j,k,0);  
        Real gy = fcz(i,j,k,1);
        Real gz = 0.5;
        int ii = (gx < 0.0) ? i - 1 : i + 1;
        int jj = (gy < 0.0) ? j - 1 : j + 1;
        int kk = k -1;
        gx = amrex::Math::abs(gx);
        gy = amrex::Math::abs(gy);
        Real gxy = gx*gy;
        Real gxz = gx*gz;
        Real gyz = gy*gz;
        Real gxyz = gx*gy*gz;
        edg_z(i,j,k,n) =
          (   1.0 - gx - gy - gz + gxy + gxz + gyz - gxyz) * phi(i ,j ,k ,n)
          + (                 gz       - gxz - gyz + gxyz) * phi(i ,j ,kk,n)
          + (            gy      - gxy       - gyz + gxyz) * phi(i ,jj,k ,n)
          + (                                  gyz - gxyz) * phi(i ,jj,kk,n)
          + (       gx           - gxy - gxz       + gxyz) * phi(ii,j ,k ,n)
          + (                            gxz       - gxyz) * phi(ii,j ,kk,n)
          + (                      gxy             - gxyz) * phi(ii,jj,k ,n)
          + (                                        gxyz) * phi(ii,jj,kk,n);
      } 
    }
  });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_centroid2facecent (Box const& bx,
                                  Array4<Real> const& phi,
                                  Array4<EBCellFlag const> const& flag,
                                  Array4<Real const> const& apx,
                                  Array4<Real const> const& apy,
                                  Array4<Real const> const& apz,
                                  Array4<Real const> const& cvol,
                                  Array4<Real const> const& ccent,
                                  Array4<Real const> const& fcx,
                                  Array4<Real const> const& fcy,
                                  Array4<Real const> const& fcz,
                                  Array4<Real> const& phi_x,
                                  Array4<Real> const& phi_y,
                                  Array4<Real> const& phi_z,
                                  int ncomp,
                                  const Box& domain,
                                  const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);

  const Box ubx = amrex::surroundingNodes(bx,0);
  const Box vbx = amrex::surroundingNodes(bx,1);
  const Box wbx = amrex::surroundingNodes(bx,2);

  // Note that ccent holds (x,y,z) of the cell centroids  as components (0/1/2)
  //           fcx   holds (  y,z) of the x-face centroid as components ( /0/1)
  //           fcy   holds (x,  z) of the y-face centroid as components (0/ /1)
  //           fcz   holds (x,y  ) of the z-face centroid as components (0/1  )

  //
  // ===================== X =====================
  //
  amrex::Loop(ubx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
      if (apx(i,j,k) == 0)
      {
          phi_x(i,j,k,n) = 1e40;
      }
      else if ( (i == domlo.x) and (bc[n].lo(0) == BCType::ext_dir) )
      {
          phi_x(i,j,k,n) = phi(domlo.x-1,j,k,n);
      }
      else if ( (i == domhi.x+1) and (bc[n].hi(0) == BCType::ext_dir) )
      {
          phi_x(i,j,k,n) = phi(domhi.x+1,j,k,n);
      }
      else if (apx(i,j,k) == 1 && cvol(i,j,k) == 1 && cvol(i-1,j,k) == 1)
      {
          phi_x(i,j,k,n) = 0.5 * ( phi(i,j,k,n) + phi(i-1,j,k,n) );
      } 
      else 
      {
          int ii = i - 1;
          int jj = (fcx(i,j,k,0) < 0.0) ? j - 1 : j + 1;
          int kk = (fcx(i,j,k,1) < 0.0) ? k - 1 : k + 1;

          std::cout << "IN CELL (i,j,k) " << i << " " << j << " " << k << std::endl;
          std::cout << "APX             " << apx(i,j,k) <<  std::endl;
          std::cout << "CVOL            " << cvol(i,j,k) << " " << cvol(i-1,j,k) << std::endl;
          std::cout << "CCENT           " << ccent(i,j,k,0) << " " << ccent(i,j,k,1) << " " << ccent(i,j,k,2) << std::endl; 
          std::cout << "fCENT           " << fcx(i,j,k,0) << " " << fcx(i,j,k,1) << std::endl; 

          Real d0       = 0.5 + ccent( i, j, k,0);                     // distance in x-dir from centroid(i  , j, k) to face(i,j,k)
          Real d1       = 0.5 - ccent(ii, j, k,0);                     // distance in x-dir from centroid(i-1, j, k) to face(i,j,k)

          Real d0_jj    = 0.5 + ccent( i,jj, k,0);                     // distance in x-dir from centroid(i  ,jj, k) to face(i,jj,k)
          Real d1_jj    = 0.5 - ccent(ii,jj, k,0);                     // distance in x-dir from centroid(i-1,jj, k) to face(i,jj,k)

          Real d0_kk    = 0.5 + ccent( i, j,kk,0);                     // distance in x-dir from centroid(i  , j,kk) to face(i,j,kk)
          Real d1_kk    = 0.5 - ccent(ii, j,kk,0);                     // distance in x-dir from centroid(i-1, j,kk) to face(i,j,kk)

          Real d0_jj_kk = 0.5 + ccent( i,jj,kk,0);                     // distance in x-dir from centroid(i  ,jj,kk) to face(i,jj,kk)
          Real d1_jj_kk = 0.5 - ccent(ii,jj,kk,0);                     // distance in x-dir from centroid(i-1,jj,kk) to face(i,jj,kk)

          Real a0       = d1 / (d0 + d1);                             // coefficient of phi(i  , j, k,n)
          Real a1       = 1.0 - a0;                                   // coefficient of phi(i-1, j, k,n)

          Real a0_jj    = d1_jj / (d0_jj + d1_jj);                    // coefficient of phi(i  ,jj, k,n)
          Real a1_jj    = 1.0 - a0_jj;                                // coefficient of phi(i-1,jj, k,n)

          Real a0_kk    = d1_kk / (d0_kk + d1_kk);                    // coefficient of phi(i  , j,kk,n)
          Real a1_kk    = 1.0 - a0_kk;                                // coefficient of phi(i-1, j,kk,n)

          Real a0_jj_kk = d1_jj_kk / (d0_jj_kk + d1_jj_kk);           // coefficient of phi(i  ,jj,kk,n)
          Real a1_jj_kk = 1.0 - a0_jj_kk;                             // coefficient of phi(i-1,jj,kk,n)

          Real phi01       = a0       *   phi(i, j, k,n)  + a1       *   phi(ii, j, k,n); // interpolation in x-dir of phi  
          Real y01         = a0       * ccent(i, j, k,1)  + a1       * ccent(ii, j, k,1); // interpolation in x-dir of y-loc
          Real z01         = a0       * ccent(i, j, k,2)  + a1       * ccent(ii, j, k,2); // interpolation in x-dir of z-loc

          Real phi01_jj    = a0_jj    *   phi(i,jj, k,n)  + a1_jj    *   phi(ii,jj, k,n); // interpolation in x-dir of phi  
          Real y01_jj      = a0_jj    * ccent(i,jj, k,1)  + a1_jj    * ccent(ii,jj, k,1); // interpolation in x-dir of y-loc
          Real z01_jj      = a0_jj    * ccent(i,jj, k,2)  + a1_jj    * ccent(ii,jj, k,2); // interpolation in x-dir of z-loc

          Real phi01_kk    = a0_kk    *   phi(i, j,kk,n)  + a1_kk    *   phi(ii, j,kk,n); // interpolation in x-dir of phi  
          Real y01_kk      = a0_kk    * ccent(i, j,kk,1)  + a1_kk    * ccent(ii, j,kk,1); // interpolation in x-dir of y-loc
          Real z01_kk      = a0_kk    * ccent(i, j,kk,2)  + a1_kk    * ccent(ii, j,kk,2); // interpolation in x-dir of z-loc

          Real phi01_jj_kk = a0_jj_kk *   phi(i,jj,kk,n)  + a1_jj_kk *   phi(ii,jj,kk,n); // interpolation in x-dir of phi  
          Real y01_jj_kk   = a0_jj_kk * ccent(i,jj,kk,1)  + a1_jj_kk * ccent(ii,jj,kk,1); // interpolation in x-dir of y-loc
          Real z01_jj_kk   = a0_jj_kk * ccent(i,jj,kk,2)  + a1_jj_kk * ccent(ii,jj,kk,2); // interpolation in x-dir of z-loc
          
          // 2D interpolation on x-face from interpolated points at (y01,z01), (y01_jj,z01_jj), (y01_kk,z01_kk), (y01_jj_kk,z01_jj_kk)
          // to centroid of x-face(i,j,k) at (fcx(i,j,k,0), fcx(i,j,k,1))

          // We order these in order to form a set of four points on the x-face ... 
          // (x0,y0) : lower left   aka (jm,km)
          // (x1,y1) : lower right  aka (jp,km)
          // (x2,y2) : upper right  aka (jp,kp)
          // (x3,y3) : upper left   aka (jm,kp)

          int jm = std::min(j,jj);
          int jp = jm+1;
          int km = std::min(k,kk);
          int kp = km+1;


          Real y,z;

          // This is the location of the face centroid relative to the central node
          if (jm == j)
             y = -0.5 + fcx(i,j,k,0);  // (j,k) is in lower half of stencil so y < 0
          else
             y =  0.5 + fcx(i,j,k,0);  // (j,k) is in upper half of stencil so y > 0

          if (km == k)
             z = -0.5 + fcx(i,j,k,1);  // (j,k) is in lower half of stencil so z < 0
          else
             z =  0.5 + fcx(i,j,k,1);  // (j,k) is in upper half of stencil so z > 0

          if (jm == j && km == k) // (j,k) is lower left, (j+1,k+1) is upper right 
          {
              std::cout << "0:CALLING WITH (JM,KM) " << jm << " " << km << std::endl;
              Real yloc_0 = -0.5+y01;
              Real zloc_0 = -0.5+z01;
              Real yloc_1 =  0.5+y01_jj;
              Real zloc_1 = -0.5+z01_jj;
              Real yloc_2 =  0.5+y01_jj_kk;
              Real zloc_2 =  0.5+z01_jj_kk;
              Real yloc_3 = -0.5+y01_kk;
              Real zloc_3 =  0.5+z01_kk;
              phi_x(i,j,k,n) = EB_interp_2D_cell_centroid_to_face_centroid(y,z,phi01,phi01_jj,phi01_jj_kk, phi01_kk,
                                                                           yloc_0, zloc_0, yloc_1, zloc_1, yloc_2, zloc_2, yloc_3, zloc_3);
          }
          else if (jm == j-1 && km == k) // (j-1,k) is lower left, (j,k+1) is upper right 
          {
              std::cout << "1:CALLING WITH (JM,KM) " << jm << " " << km << std::endl;
              phi_x(i,j,k,n) = EB_interp_2D_cell_centroid_to_face_centroid(y,z,phi01_jj,phi01,phi01_kk, phi01_jj_kk,
                                                                           y01_jj,z01_jj, y01,z01, y01_kk,z01_kk, y01_jj_kk, z01_jj_kk);
          }
          else if (jm == j && km == k-1) // (j,k-1) is lower left, (j+1,k) is upper right 
          {
              std::cout << "2:CALLING WITH (JM,KM) " << jm << " " << km << std::endl;
              phi_x(i,j,k,n) = EB_interp_2D_cell_centroid_to_face_centroid(y,z,phi01_kk,phi01_jj_kk,phi01_jj, phi01,
                                                                           y01_kk,z01_kk, y01_jj_kk,z01_jj_kk, y01_jj,z01_jj, y01, z01);
          }
          else if (jm == j-1 && km == k-1) // (j-1,k-1) is lower left, (j,k) is upper right 
          {
              std::cout << "3:CALLING WITH (JM,KM) " << jm << " " << km << std::endl;
              phi_x(i,j,k,n) = EB_interp_2D_cell_centroid_to_face_centroid(y,z,phi01_jj_kk,phi01_kk,phi01, phi01_jj,
                                                                           y01_jj_kk,z01_jj_kk, y01_kk,z01_kk, y01,z01, y01_jj, z01_jj);
          }
          else 
            amrex::Abort("Bad option in interpolation from cell centroid to face centroid!");
      }
    });
  
#if 0
  //
  // ===================== Y =====================
  //
  amrex::Loop(vbx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apy(i,j,k) == 0)
        {
            phi_y(i,j,k,n) = 1e40;
        }
        else if ( (j == domlo.y) and (bc[n].lo(1) == BCType::ext_dir) )
        {
            phi_y(i,j,k,n) = phi(i,domlo.y-1,k,n);
        }
        else if ( (j == domhi.y+1) and (bc[n].hi(1) == BCType::ext_dir) )
        {
            phi_y(i,j,k,n) = phi(i,domhi.y+1,k,n);
        }
        else if (apy(i,j,k) == 1 && cvol(i,j,k) == 1 && cvol(i,j-1,k) == 1)
        {
            phi_y(i,j,k,n) = 0.5 * (phi(i,j  ,k,n) + phi(i,j-1,k,n));
        }
        else 
        {
          int jj = j - 1;
          int ii = (fcy(i,j,k,0) < 0.0) ? i - 1 : i + 1;
          int kk = (fcy(i,j,k,1) < 0.0) ? k - 1 : k + 1;

          Real d0 = 0.5 + ccent(i, j,k,1); // distance in y-dir from centroid(i,j  ) to face(i,j)
          Real d1 = 0.5 - ccent(i,jj,k,1); // distance in y-dir from centroid(i,j-1) to face(i,j)

          Real a0    = d1 / (d0 + d1);                              // coefficient of phi(i,j  ,k,n)
          Real a1    = d0 / (d0 + d1);                              // coefficient of phi(i,j-1,k,n)

          Real phi01 = a0 *   phi(i,j,k,n)  + a1 *   phi(i,jj,k,n); // interpolation in y-dir of phi  
          Real x01   = a0 * ccent(i,j,k,0)  + a1 * ccent(i,jj,k,0); // interpolation in y-dir of x-loc

          Real d0_ii = 0.5 + ccent(ii, j,k,1); // distance in y-dir from centroid(ii,j  ) to face(ii,j)
          Real d1_ii = 0.5 - ccent(ii,jj,k,1); // distance in y-dir from centroid(ii,j-1) to face(ii,j)

          Real a0_ii    = d1_ii / (d0_ii + d1_ii);                             // coefficient of phi(i  ,jj,k,n)
          Real a1_ii    = d0_ii / (d0_ii + d1_ii);                             // coefficient of phi(i-1,jj,k,n)
          Real phi01_ii = a0_ii *   phi(ii,j,k,n)  + a1_ii *   phi(ii,jj,k,n); // interpolation in y-dir of phi 
          Real x01_ii   = a0_ii * ccent(ii,j,k,0)  + a1_ii * ccent(ii,jj,k,0); // interpolation in y-dir of x-loc
          
          // 2D interpolation on x-face from interpolated points to centroid of x-face(i,j,k)
          phi_y(i,j,k,n) = ( (     fcy(i,j,k,0) - x01   ) * phi01_ii + 
                             (1. - fcy(i,j,k,0) + x01_ii) * phi01    ) / (1. - x01 + x01_ii);
        }
    });

  //
  // ===================== Z =====================
  //
  amrex::Loop(wbx, ncomp, [=] (int i, int j, int k, int n) noexcept
    {
        if (apz(i,j,k) == 0)
        {
            phi_z(i,j,k,n) = 1e40;
        }
        else if ( (k == domlo.z) and (bc[n].lo(2) == BCType::ext_dir) )
        {
            phi_z(i,j,k,n) = phi(i,j,domlo.z-1,n);
        }
        else if ( (k == domhi.z+1) and (bc[n].hi(2) == BCType::ext_dir) )
        {
            phi_z(i,j,k,n) = phi(i,j,domhi.z+1,n);
        }
        else if (apz(i,j,k) == 1 && cvol(i,j,k) == 1 && cvol(i,j,k-1) == 1)
        {
            phi_z(i,j,k,n) = 0.5 * (phi(i,j,k,n) + phi(i,j,k-1,n));
        }
        else 
        {
          int kk = k - 1;
          int ii = (fcz(i,j,k,0) < 0.0) ? i - 1 : i + 1;
          int jj = (fcz(i,j,k,1) < 0.0) ? j - 1 : j + 1;

          Real d0    = 0.5 + ccent(i,j,k ,2);  // distance in z-dir from centroid(i,j,k  ) to face(i,j,k)
          Real d1    = 0.5 - ccent(i,j,kk,2);  // distance in z-dir from centroid(i,j,k-1) to face(i,j,k)

          Real d0_ii = 0.5 + ccent(ii,j, k,2); // distance in z-dir from centroid(ii,j,k  ) to face(ii,j,k)
          Real d1_ii = 0.5 - ccent(ii,j,kk,2); // distance in z-dir from centroid(ii,j,k-1) to face(ii,j,k)

          Real a0    = d1 / (d0 + d1);                              // coefficient of phi(i,j,k  ,n)
          Real a1    = d0 / (d0 + d1);                              // coefficient of phi(i,j,k-1,n)

          Real a0_ii    = d1_ii / (d0_ii + d1_ii);                  // coefficient of phi(i ,j,kk,n)
          Real a1_ii    = d0_ii / (d0_ii + d1_ii);                  // coefficient of phi(ii,j,kk,n)

          Real phi01 = a0 *   phi(i,j,k,n)  + a1 *   phi(i,j,kk,n); // interpolation of phi   in z-direction
          Real x01   = a0 * ccent(i,j,k,0)  + a1 * ccent(i,j,kk,0); // interpolation of x-loc in z-direction
          Real y01   = a0 * ccent(i,j,k,1)  + a1 * ccent(i,j,kk,1); // interpolation of y-loc in z-direction

          Real phi01_ii = a0_ii *   phi(ii,j,k,n)  + a1_ii *   phi(ii,j,kk,n); // interpolation of phi   in z-direction
          Real x01_ii   = a0_ii * ccent(ii,j,k,0)  + a1_ii * ccent(ii,j,kk,0); // interpolation of x-loc in z-direction
          
          // 2D interpolation on x-face from interpolated points to centroid of x-face(i,j,k)
          phi_z(i,j,k,n) = ( (     fcz(i,j,k,0) - x01   ) * phi01_ii + 
                             (1. - fcz(i,j,k,0) + x01_ii) * phi01    ) / (1. - x01 + x01_ii);
        }
    });
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_interp_cc2face (Box const& box,
                        Array4<Real> const& phi,
                        Array4<Real> const& edg_x,
                        Array4<Real> const& edg_y,
                        Array4<Real> const& edg_z,
                        int ncomp,
                        const Box& domain,
                        const BCRec* bc) noexcept
{
  const Dim3 domlo = amrex::lbound(domain);
  const Dim3 domhi = amrex::ubound(domain);

  const Box ubx = amrex::surroundingNodes(box,0);
  const Box vbx = amrex::surroundingNodes(box,1);
  const Box wbx = amrex::surroundingNodes(box,2);

  //
  //===================== X =====================
  //
  amrex::Loop(ubx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if ( (i == domlo.x) and (bc[n].lo(0) == BCType::ext_dir) )
    {
      edg_x(i,j,k,n) = phi(domlo.x-1,j,k,n);
    }
    else if ( (i == domhi.x+1) and (bc[n].hi(0) == BCType::ext_dir) )
    {
      edg_x(i,j,k,n) = phi(domhi.x+1,j,k,n);
    }
    else
    {
      edg_x(i,j,k,n) = 0.5 * ( phi(i,j,k,n) + phi(i-1,j,k,n) );
    } 
  });

  //
  //===================== Y =====================
  //
  amrex::Loop(vbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if ( (j == domlo.y) and (bc[n].lo(1) == BCType::ext_dir) )
    {
      edg_y(i,j,k,n) = phi(i,domlo.y-1,k,n);
    }
    else if ( (j == domhi.y+1) and (bc[n].hi(1) == BCType::ext_dir) )
    {
      edg_y(i,j,k,n) = phi(i,domhi.y+1,k,n);
    }
    else
    {
      edg_y(i,j,k,n) = 0.5 * (phi(i,j  ,k,n) + phi(i,j-1,k,n));
    }
  });

  //
  //===================== Z =====================
  //
  amrex::Loop(wbx, ncomp, [=] (int i, int j, int k, int n) noexcept
  {
    if ( (k == domlo.z) and (bc[n].lo(2) == BCType::ext_dir) )
    {
      edg_z(i,j,k,n) = phi(i,j,domlo.z-1,n);
    }
    else if ( (k == domhi.z+1) and (bc[n].hi(2) == BCType::ext_dir) )
    {
      edg_z(i,j,k,n) = phi(i,j,domhi.z+1,n);
    }
    else
    {
      edg_z(i,j,k,n) = 0.5 * (phi(i,j  ,k,n) + phi(i,j,k-1,n));
    }
  });
}

}
#endif
