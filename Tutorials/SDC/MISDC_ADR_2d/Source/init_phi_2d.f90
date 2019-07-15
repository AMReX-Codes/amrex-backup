subroutine init_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi) bind(C, name="init_phi")
  !  Initialize the scalar field phi
  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 
  real(amrex_real), intent(in   ) :: prob_hi(2) 

  integer          :: i,j
  double precision :: x,y,tupi,t0,d, pi

  tupi=3.14159265358979323846d0*2d0
  pi=3.14159265358979323846d0
  d=0.1
  t0=0.0025d0/d

!print*, prob_lo, prob_hi, philo,phihi
  do j = philo(2), phihi(2)
     y = prob_lo(2) + (dble(j)+1.d0/2.d0 )* dx(2)

     do i = philo(1), phihi(1)
        x = prob_lo(1) + (dble(i)+1.d0/2.d0) * dx(1)

       ! phi(i,j) =sin(x*tupi)*sin(y*tupi)
       ! phi(i,j)=exp(-x*x/(4.0d0*d*t0))*exp(-y*y/(4.0d0*d*t0))
       ! phi(i,j)=cos(pi*(x+y))
        phi(i,j) = 0;
     end do
  end do
end subroutine init_phi

subroutine init_beta(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi,epsilon) bind(C, name="init_beta")
!  Initialize the scalar field phi
use amrex_fort_module, only : amrex_real

implicit none

integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
real(amrex_real), intent(in   ) :: dx(2)
real(amrex_real), intent(in   ) :: prob_lo(2)
real(amrex_real), intent(in   ) :: prob_hi(2)
real(amrex_real), intent(in   ) :: epsilon

integer          :: i,j
double precision :: x,y,pi

pi=3.14159265358979323846d0

do j = philo(2), phihi(2)
y = prob_lo(2) + (dble(j)+(1.d0/2.d0)) * dx(2)
do i = philo(1), phihi(1)
x = prob_lo(1) + (dble(i)+(1.d0/2.d0)) * dx(1)

phi(i,j)=1.d0+epsilon*sin(pi*(x+y))
end do
end do
end subroutine init_beta

subroutine err_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi,a,d,r,time) bind(C, name="err_phi")
  !  Subtract the exact solution from phi.  This will only work for special initial conditions
  !  We use the exact discretized form for diffusion and reaction, and exact translation for advection
  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 
  real(amrex_real), intent(in   ) :: prob_hi(2) 
  real(amrex_real), intent(in   ) :: a,d,r
  real(amrex_real), intent(in   ) :: time

  integer          :: i,j,kx,ky,nbox
  double precision :: x,y,sym,tupi, maxphi,xx,yy,t0,pi

  tupi=3.14159265358979323846d0*2d0
  pi=3.14159265358979323846d0

  !  Form the diffusion coefficient for the 2nd order Laplacian produces
!  sym=d*(-2.0d0+2.0d0*cos(tupi*dx(1)))/(dx(1)*dx(1))
!  sym=sym+d*(-2.0d0+2.0d0*cos(tupi*dx(2)))/(dx(2)*dx(2))
!  sym=    d*(-3.0d1+32.0d0*cos(tupi*dx(1))-2.0d0*cos(tupi*2.0d0*dx(1)))/(1.2d1*dx(1)*dx(1))
 ! sym=sym+d*(-3.0d1+32.0d0*cos(tupi*dx(2))-2.0d0*cos(tupi*2.0d0*dx(2)) )/(1.2d1*dx(2)*dx(2))
 ! sym=sym-r    !  Add reaction

 ! t0=0.0025d0/d
  !nbox = ceiling(sqrt(4.0*d*(t0+time)*37.0)/2.0)  !  Decide how many periodic images
  
  do j = philo(2), phihi(2)
     y = prob_lo(2) + (dble(j)+1.d0/2.d0) * dx(2) !+a*time
     do i = philo(1), phihi(1)
        x = prob_lo(1) + (dble(i)+1.d0/2.d0) * dx(1) !+ a*time
        !        phi(i,j) = phi(i,j)-sin(x*tupi)*sin(y*tupi)*exp(time*sym)

       ! do kx = -nbox,nbox
       ! do ky = -nbox,nbox
        !  xx = x + real(2*kx)
         ! yy = y + real(2*ky)
!print*, phi(i,j), (time**2)/2.d0
         ! phi(i,j) = phi(i,j) - t0/(t0+time)*exp(-xx*xx/(4.0d0*d*(t0+time)))*exp(-yy*yy/(4.0d0*d*(t0+time)))
           ! print*, 'phi', phi(i,j)
            ! phi(i,j) = phi(i,j) - exp(-2.d0*d*(pi**2)*time)*cos(pi*(x+y))
           ! print*, phi(i,j)
            phi(i,j) = phi(i,j) - (time**2)/2.d0
      ! end do
      ! end do

        
     end do
  end do
  
end subroutine err_phi
