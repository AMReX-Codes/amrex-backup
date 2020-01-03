subroutine init_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi, k_freq, Nprob) bind(C, name="init_phi")
  !  Initialize the scalar field phi
  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2), Nprob
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 
  real(amrex_real), intent(in   ) :: prob_hi(2) 
  real(amrex_real), intent(in   ) :: k_freq
  integer          :: i,j, i_quad, j_quad
  double precision :: x,y,tupi,t0,d, pi, x_quad, y_quad,om
double precision :: gauss_nodeFrac(0:2)
double precision :: gauss_weights(0:2)
gauss_nodeFrac = (/ (1.d0-(3.d0/5.d0)**(0.5d0))/2.d0,0.5d0,(1.d0+(3.d0/5.d0)**(0.5d0))/2.d0 /)
gauss_weights = (/ (5.d0/18.d0),(8.d0/18.d0),(5.d0/18.d0)/)
  tupi=3.14159265358979323846d0*2d0
  pi=3.14159265358979323846d0

om=k_freq*pi
!print*, prob_lo, prob_hi, philo,phihi
  do j = lo(2), hi(2)
     !y = prob_lo(2) + (dble(j)+1.d0/2.d0 )* dx(2)
     y = prob_lo(2) + dble(j)* dx(2)
     do i = lo(1), hi(1) ! Went philo - phihi before
        !x = prob_lo(1) + (dble(i)+1.d0/2.d0) * dx(1)
        x = prob_lo(1) + dble(i) * dx(1)

       ! phi(i,j) =sin(x*tupi)*sin(y*tupi)
       ! phi(i,j)=exp(-x*x/(4.0d0*d*t0))*exp(-y*y/(4.0d0*d*t0))
        phi(i,j) = 0.d0
        do j_quad = 0,2
        y_quad = y + dx(2)*gauss_nodeFrac(j_quad)
        !print*, y_quad
        do i_quad = 0,2
        x_quad = x + dx(1)*gauss_nodeFrac(i_quad)
        if ((Nprob .EQ. 1) .OR. (Nprob .EQ. 4)) then
            phi(i,j) = phi(i,j)+ gauss_weights(j_quad)*gauss_weights(i_quad)* &
                        sin(om*(x_quad))*sin(om*(y_quad))
        elseif (Nprob .EQ. 3) then
            phi(i,j) = phi(i,j)+ gauss_weights(j_quad)*gauss_weights(i_quad)* &
                cos(om*(x_quad+y_quad))
        endif


        end do
        end do
       ! phi(i,j) = 0

       ! phi(lo(1)+2,hi(2)-176) = 1

        !phi(i,j)=cos(pi*(x+y))
       ! phi(i,j) = 0;
     end do
  end do
end subroutine init_phi

subroutine init_beta(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi,epsilon, k_freq, Nprob) bind(C, name="init_beta")
!  Initialize the scalar field phi
use amrex_fort_module, only : amrex_real

implicit none

integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2), Nprob
real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
real(amrex_real), intent(in   ) :: dx(2)
real(amrex_real), intent(in   ) :: prob_lo(2)
real(amrex_real), intent(in   ) :: prob_hi(2)
real(amrex_real), intent(in   ) :: epsilon, k_freq

integer          :: i,j, i_quad, j_quad
double precision :: x,y,pi, x_quad, y_quad,om
double precision :: gauss_nodeFrac(0:4)
double precision :: gauss_weights(0:4)
!gauss_nodeFrac = (/ (1.d0-(3.d0/5.d0)**(0.5d0))/2.d0,0.5d0,(1.d0+(3.d0/5.d0)**(0.5d0))/2.d0 /)
!gauss_weights = (/ (5.d0/18.d0),(8.d0/18.d0),(5.d0/18.d0)/)
pi=3.14159265358979323846d0
 om=k_freq*pi
! Even higher order initialization!
gauss_nodeFrac = (/ (1.d0 - ((1.d0/3.d0)*dsqrt(5.d0+2.d0*dsqrt(10.d0/7.d0))))/2.d0, &
                    (1.d0 - ((1.d0/3.d0)*dsqrt(5.d0-2.d0*dsqrt(10.d0/7.d0))))/2.d0, &
                    0.5d0,&
                    (1.d0 + ((1.d0/3.d0)*dsqrt(5.d0-2.d0*dsqrt(10.d0/7.d0))))/2.d0, &
                    (1.d0 + ((1.d0/3.d0)*dsqrt(5.d0+2.d0*dsqrt(10.d0/7.d0))))/2.d0 /)

gauss_weights = (/(322.d0-13.d0*dsqrt(70.d0))/1800.d0 , &
                   (322.d0+13.d0*dsqrt(70.d0))/1800.d0, &
                  128.d0/450.d0, &
                  (322.d0+13.d0*dsqrt(70.d0))/1800.d0, &
                  (322.d0-13.d0*dsqrt(70.d0))/1800.d0/)




do j = philo(2), phihi(2)
!y = prob_lo(2) + (dble(j)+(1.d0/2.d0)) * dx(2)
y = prob_lo(2) + dble(j) * dx(2)
    do i = philo(1), phihi(1)
    !x = prob_lo(1) + (dble(i)+(1.d0/2.d0)) * dx(1)
    x = prob_lo(1) + dble(i) * dx(1)
    !phi(i,j)=1.d0+epsilon*sin(om*(x+y))
    phi(i,j) = 0.d0
        do j_quad = 0,4
        y_quad = y + dx(2)*gauss_nodeFrac(j_quad)
        !print*, y_quad
            do i_quad = 0,4
            x_quad = x + dx(1)*gauss_nodeFrac(i_quad)
            if ((Nprob .EQ. 1) .OR. (Nprob .EQ. 4)) then
                phi(i,j) = phi(i,j)+ gauss_weights(j_quad)*gauss_weights(i_quad)* &
                    (1.d0+epsilon*sin(om*(x_quad))*sin(om*(y_quad)))
            elseif (Nprob .EQ. 3) then
                phi(i,j) = phi(i,j)+ gauss_weights(j_quad)*gauss_weights(i_quad)* &
                    (1.d0+epsilon*sin(om*(x_quad+y_quad)))
            endif

            end do
        end do
    end do
end do
end subroutine init_beta

subroutine err_phi(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi,a,d,r,time, epsilon, k_freq, kappa, Nprob) bind(C, name="err_phi")
  !  Subtract the exact solution from phi.  This will only work for special initial conditions
  !  We use the exact discretized form for diffusion and reaction, and exact translation for advection
  use amrex_fort_module, only : amrex_real

  implicit none

  integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2), Nprob
  real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in   ) :: dx(2) 
  real(amrex_real), intent(in   ) :: prob_lo(2) 
  real(amrex_real), intent(in   ) :: prob_hi(2) 
  real(amrex_real), intent(in   ) :: a,d,r
  real(amrex_real), intent(in   ) :: time
real(amrex_real), intent(in   ) :: epsilon, k_freq, kappa

  integer          :: i,j,kx,ky,nbox, i_quad, j_quad
  double precision :: x,y,sym,tupi, maxphi,xx,yy,t0,pi,x_quad,y_quad,om
double precision :: gauss_nodeFrac(0:2)
double precision :: gauss_weights(0:2)
gauss_nodeFrac = (/ (1.d0-(3.d0/5.d0)**(0.5d0))/2.d0,0.5d0,(1.d0+(3.d0/5.d0)**(0.5d0))/2.d0 /)
gauss_weights = (/ (5.d0/18.d0),(8.d0/18.d0),(5.d0/18.d0)/)



  tupi=3.14159265358979323846d0*2d0
  pi=3.14159265358979323846d0
 om=k_freq*pi

  do j = philo(2), phihi(2)
    !y = prob_lo(2) + (dble(j)+1.d0/2.d0) * dx(2) !+a*time
    y = prob_lo(2) + dble(j) * dx(2)
     do i = philo(1), phihi(1)
        !x = prob_lo(1) + (dble(i)+1.d0/2.d0) * dx(1) !+ a*time
        x = prob_lo(1) + dble(i) * dx(1)
        !        phi(i,j) = phi(i,j)-sin(x*tupi)*sin(y*tupi)*exp(time*sym)

        !do kx = -nbox,nbox
        !do ky = -nbox,nbox
          !xx = x + real(2*kx)
          !yy = y + real(2*ky)
         ! phi(i,j) = phi(i,j) - t0/(t0+time)*exp(-xx*xx/(4.0d0*d*(t0+time)))*exp(-yy*yy/(4.0d0*d*(t0+time)))
           ! print*, 'phi', phi(i,j)

        !phi(i,j) = phi(i,j) - exp(-2.d0*d*(pi**2)*time)*cos(pi*(x+y))

            ! 5th order quadrature if necessary
            do j_quad = 0,2
            y_quad = y + dx(2)*gauss_nodeFrac(j_quad)
            !print*, y_quad
                do i_quad = 0,2
                x_quad = x + dx(1)*gauss_nodeFrac(i_quad)
                if ((Nprob .EQ. 1) .OR. (Nprob .EQ. 4)) then
                    phi(i,j) = phi(i,j)- gauss_weights(j_quad)*gauss_weights(i_quad)* &
                        exp(-kappa*time)*sin(om*(x_quad))*sin(om*(y_quad))
                elseif (Nprob .EQ. 3) then
                    phi(i,j) = phi(i,j)- gauss_weights(j_quad)*gauss_weights(i_quad)* &
                        exp(-kappa*time)*cos(om*(x_quad+y_quad))
                endif

                end do
            end do

       !end do
       !end do

        
     end do
  end do
  
end subroutine err_phi

subroutine print_multifab(phi, philo, phihi) bind(C, name="print_multifab")
!  Subtract the exact solution from phi.  This will only work for special initial conditions
!  We use the exact discretized form for diffusion and reaction, and exact translation for advection
use amrex_fort_module, only : amrex_real

implicit none

integer, intent(in) ::  philo(2), phihi(2)
real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))


integer :: i,j


!  Form the diffusion coefficient for the 2nd order Laplacian produces
!  sym=d*(-2.0d0+2.0d0*cos(tupi*dx(1)))/(dx(1)*dx(1))
!  sym=sym+d*(-2.0d0+2.0d0*cos(tupi*dx(2)))/(dx(2)*dx(2))
!  sym=    d*(-3.0d1+32.0d0*cos(tupi*dx(1))-2.0d0*cos(tupi*2.0d0*dx(1)))/(1.2d1*dx(1)*dx(1))
! sym=sym+d*(-3.0d1+32.0d0*cos(tupi*dx(2))-2.0d0*cos(tupi*2.0d0*dx(2)) )/(1.2d1*dx(2)*dx(2))
! sym=sym-r    !  Add reaction

! t0=0.0025d0/d
!nbox = ceiling(sqrt(4.0*d*(t0+time)*37.0)/2.0)  !  Decide how many periodic images

do j = philo(2), phihi(2)
!y = prob_lo(2) + (dble(j)+1.d0/2.d0) * dx(2) !+a*time

do i = philo(1), phihi(1)

 print*, 'index ', i , j , ' with val ', phi(i,j)

end do
end do

end subroutine print_multifab
