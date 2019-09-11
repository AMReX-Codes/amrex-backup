
subroutine SDC_feval_F(lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi,f, flo,fhi, &
                         dx,a,d,r,facex, facexlo, facexhi,facey, faceylo, faceyhi,prodx, prodxlo, prodxhi,prody, prodylo, prodyhi,n,time, epsilon, k_freq, kappa,Nprob) bind(C, name="SDC_feval_F")
 ! facex, facexlo, facexhi,facey, faceylo, faceyhi
  !  Compute the rhs terms
! Assumes you have 2 ghost cells for flux that way you can do product rule

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2),flo(2), fhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(inout) :: f   (flo(1):fhi(1),flo(2):fhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), intent(in)    :: a,d,r,time, epsilon, k_freq, kappa
  integer, intent(in) :: n,Nprob
 integer facexlo(2), facexhi(2)
real(amrex_real), intent(in) :: facex( facexlo(1): facexhi(1), facexlo(2): facexhi(2))
 integer faceylo(2), faceyhi(2)
real(amrex_real), intent(in) :: facey( faceylo(1): faceyhi(1), faceylo(2): faceyhi(2))
integer prodxlo(2), prodxhi(2)
real(amrex_real), intent(inout) :: prodx( prodxlo(1): prodxhi(1), prodxlo(2): prodxhi(2))
integer prodylo(2), prodyhi(2)
real(amrex_real), intent(inout) :: prody( prodylo(1): prodyhi(1), prodylo(2): prodyhi(2))

integer          :: i,j, i_quad, j_quad
double precision :: x,y,pi, x_quad, y_quad
real(amrex_real) :: phi_one, face_one

double precision :: gauss_nodeFrac(0:2)
double precision :: gauss_weights(0:2)
gauss_nodeFrac = (/ (1.d0-(3.d0/5.d0)**(0.5d0))/2.d0,0.5d0,(1.d0+(3.d0/5.d0)**(0.5d0))/2.d0 /)
gauss_weights = (/ (5.d0/18.d0),(8.d0/18.d0),(5.d0/18.d0)/)

pi=3.14159265358979323846d0

  select case(n)
     case (0)  !  Explicit term (here it is advection)
        ! x-fluxes
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              fluxx(i,j) = (-phi(i+1,j)+ 7.0d0*( phi(i,j) + phi(i-1,j) )-phi(i-2,j)) / 12.0d0
           end do
        end do
        
        ! y-fluxes
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
              fluxy(i,j) =(-phi(i,j+1)+7.0d0*( phi(i,j) + phi(i,j-1) )-phi(i,j-2)) / 12.0d0
           end do
        end do
        !print*, lo, hi, domlo, domhi
        !  Function value is divergence of flux
        do j = lo(2), hi(2)
            y = -1.d0 + dble(j) * dx(2)
        ! print*, y
           do i = lo(1), hi(1)
                x = -1.d0 + dble(i) * dx(1)
              !f(i,j) =  a*((fluxx(i+1,j  ) - fluxx(i,j))/dx(1) &
               !    + (fluxy(i  ,j+1) - fluxy(i,j))/dx(2))
                ! 5th order quadrature loop
                f(i,j) = 0.d0
                do j_quad = 0,2
                    y_quad = y + dx(2)*gauss_nodeFrac(j_quad)
                    !print*, y_quad
                do i_quad = 0,2
                    x_quad = x + dx(1)*gauss_nodeFrac(i_quad)
                    f(i,j) = f(i,j)+ gauss_weights(j_quad)*gauss_weights(i_quad)* &
                            (&
                            -kappa*exp(-kappa*time)*sin(k_freq*(x_quad))*sin(k_freq*(y_quad)) &
                            -d*k_freq*( &
                                -2.d0*k_freq*exp(-kappa*time)*sin(k_freq*(x_quad))*sin(k_freq*(y_quad)) &
                                + epsilon*k_freq*exp(-kappa*time)*(((sin(k_freq*y_quad))**2)*((cos(k_freq*x_quad))**2) + ((cos(k_freq*y_quad))**2)*((sin(k_freq*x_quad))**2)) &
                                -2.d0*epsilon*k_freq*exp(-kappa*time)*((sin(k_freq*x_quad))**2)*((sin(k_freq*y_quad))**2) &
                            )&
                            )


                end do
                end do
                !f(i,j) = 0


           end do
           
        end do
     case (1)  !  First implicit piece (here it is diffusion)
     ! We let beta vary in space via cell centered data bcc.

      !  do j = philo(2), phihi(2)
      !  do i = philo(1), phihi(1)
      !  ! fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
      !  print*, i, " x index ", j, " y index ", phi(i,j)
      !  end do
      !  end do


        ! x-fluxes
        do j = philo(2), phihi(2)
        !do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              ! fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
              fluxx(i,j) = ( -phi(i+1,j)  +15.0d0*(phi(i,j) - phi(i-1,j)) + phi(i-2,j) ) /(12.0d0*dx(1))
           end do
        end do
        
        ! y-fluxes
        do j = lo(2), hi(2)+1
           do i = philo(1), phihi(1)
           !do i = lo(1), hi(1)
              ! fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx(2)
              fluxy(i,j) = ( -phi(i,j+1)  +15.0d0*(phi(i,j) - phi(i,j-1)) + phi(i,j-2) ) /(12.0d0*dx(2))
           end do
        end do

       ! We now apply the product rule to find 4th order approximation
        ! x-flux products
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               ! no dx in one variables!
               phi_one = (-5.d0*fluxx(i,j+2)+34.d0*(fluxx(i,j+1)-fluxx(i,j-1))+5.d0*fluxx(i,j-2))/48.d0
               face_one =(-5.d0*facex(i,j+2)+34.d0*(facex(i,j+1)-facex(i,j-1))+5.d0*facex(i,j-2))/48.d0
               !prodx(i,j) = fluxx(i,j)*facex(i,j) !+ phi_one*face_one/12.d0 ! right added term changes from L^2,4,4 to L^4,4,4
               prodx(i,j) = fluxx(i,j)*facex(i,j) + phi_one*face_one/12.d0 ! right added term changes from L^2,4,4 to L^4,4,4
            end do
        end do

        ! y-fluxes
        do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
                ! no dx in one variables!
                phi_one = (-5.d0*fluxy(i+2,j)+34.d0*(fluxy(i+1,j)-fluxy(i-1,j))+5.d0*fluxy(i-2,j))/48.d0
                face_one =(-5.d0*facey(i+2,j)+34.d0*(facey(i+1,j)-facey(i-1,j))+5.d0*facey(i-2,j))/48.d0
                !prody(i,j) = fluxy(i,j)*facey(i,j) !+ phi_one*face_one/12.d0 ! right added term changes from L^2,4,4 to L^4,4,4
                prody(i,j) = fluxy(i,j)*facey(i,j) + phi_one*face_one/12.d0 ! right added term changes from L^2,4,4 to L^4,4,4
            end do
        end do

        !  Function value is divergence of flux
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              f(i,j) =  d*((prodx(i+1,j  ) - prodx(i,j))/dx(1) &
                   + (prody(i  ,j+1) - prody(i,j))/dx(2))
           end do
        end do










     case (2)  !  Second implicit piece (here it is reaction)
        do j = lo(2), hi(2)
            y = -1.d0 + dble(j) * dx(2)
           do i = lo(1), hi(1)
                x = -1.d0 + (dble(i)) * dx(1)
!              f(i,j) =  r*phi(i,j)*(1.0d0-phi(i,j))*(0.5d0-phi(i,j))
              f(i,j) =  -r*phi(i,j)
               ! f(i,j) = r*phi(i,j)*sin(pi*(x+y))
             end do
          end do
     case default
        print *, 'bad case in advance_2d'
     end select
     
   end subroutine SDC_feval_F
   

   subroutine SDC_fcomp_reaction_F (lo, hi, domlo, domhi, phi, philo, phihi, &
        rhs, rhslo, rhshi, &
        f, flo,fhi, dtq,n,Nprob) bind(C, name="SDC_fcomp_reaction_F")

  !  Solve for the reaction term
  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2)
  integer rhslo(2), rhshi(2)
  integer flo(2), fhi(2)
  real(amrex_real), intent(inout) :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in)    :: rhs  (rhslo(1):rhshi(1),rhslo(2):rhshi(2))  
  real(amrex_real), intent(in)     :: f   (flo(1):fhi(1),flo(2):fhi(2))
  real(amrex_real), intent(in)    :: dtq
  integer, intent(in)             :: n,Nprob
  
  ! local variables
  integer i,j
  real(amrex_real) c

  !  Function 
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        c = 1.0d0-dtq*(1.0d0-phi(i,j))*(0.5d0-phi(i,j))
        c = 1.0d0+dtq
        phi(i,j) =  rhs(i,j)/c
     end do
  end do

end subroutine SDC_fcomp_reaction_F

subroutine SDC_Lresid_F (lo, hi, domlo, domhi, phi, philo, phihi, &
        rhs, rhslo, rhshi, &
        res, reslo, reshi, &
        corr, corrlo, corrhi, &        
        dtq,dx,Nprob) bind(C, name="SDC_Lresid_F")

  !  Compute the residual for the diffusion operator
  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2)
  integer rhslo(2), rhshi(2)
  integer reslo(2), reshi(2)
  integer corrlo(2), corrhi(2)  
  real(amrex_real), intent(in) :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(in)    :: rhs  (rhslo(1):rhshi(1),rhslo(2):rhshi(2))  
  real(amrex_real), intent(inout)    :: res  (reslo(1):reshi(1),reslo(2):reshi(2))  
  real(amrex_real), intent(inout)    :: corr  (corrlo(1):corrhi(1),corrlo(2):corrhi(2))  
  real(amrex_real), intent(in)    :: dtq
  real(amrex_real), intent(in)    :: dx(2)
  integer, intent(in) :: Nprob
  ! local variables
  integer i,j
  real(amrex_real) Lap

  !  Function 
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
!        Lap=(phi(i-1,j)-2.0d0*phi(i,j)+phi(i+1,j))/(dx(1)*dx(1)) + (phi(i,j-1)-2.0d0*phi(i,j)+phi(i,j+1))/(dx(2)*dx(2))
        !  First compute Grad(phi) and then do multiplication
        Lap=(-(phi(i-2,j)+phi(i+2,j))+1.6d1*(phi(i+1,j)+phi(i-1,j))- 3.0d1*phi(i,j))/(1.2d1*dx(1)*dx(1)) 
        Lap=Lap+(-(phi(i,j-2)+phi(i,j+2))+1.6d1*(phi(i,j+1)+phi(i,j-1))- 3.0d1*phi(i,j))/(1.2d1*dx(2)*dx(2)) 
        res(i,j)=rhs(i,j)- (phi(i,j) -dtq*Lap)
        corr(i,j)=Lap     !0.0d0
     end do
  end do

end subroutine SDC_Lresid_F

subroutine cc_to_face_loc(lo, hi, cc_dat, cc_lo, cc_hi, face_dat, face_lo,face_hi, dir,Nprob) bind(C, name="cc_to_face_loc")
!  move cell centered data to be face centered
use amrex_fort_module, only : amrex_real

implicit none

integer, intent(in) :: lo(3), hi(3), cc_lo(3), cc_hi(3), face_lo(3), face_hi(3), dir
real(amrex_real), intent(inout) :: cc_dat(cc_lo(1):cc_hi(1),cc_lo(2):cc_hi(2),cc_lo(3):cc_hi(3))
real(amrex_real), intent(inout) :: face_dat(face_lo(1):face_hi(1),face_lo(2):face_hi(2),face_lo(3):face_hi(3))
integer, intent(in) :: Nprob

integer          :: i,j, k


!print *,'lo, hi, cc_lo, cc_hi, face_lo, face_hi=',lo, hi, cc_lo, cc_hi, face_lo, face_hi

!print *,'cc_dat(0,0,0)=',cc_dat(0,0,0)


! We condition based on direction
if (dir .EQ. 0) then
do i = lo(1),hi(1) + 1
do j = cc_lo(2),cc_hi(2)
do k = cc_lo(3),cc_hi(3)
face_dat(i,j,k) = (-1.d0*cc_dat(i-2,j,k) + 7.d0*cc_dat(i-1,j,k) &
+ 7.d0*cc_dat(i,j,k) - 1.d0*cc_dat(i+1,j,k))/(12.d0)

end do
end do
end do
elseif(dir .EQ. 1) then

do i = cc_lo(1),cc_hi(1)
do j = lo(2),hi(2) + 1
do k = cc_lo(3),cc_hi(3)

face_dat(i,j,k) = (-1.d0*cc_dat(i,j-2,k) + 7.d0*cc_dat(i,j-1,k) &
+ 7.d0*cc_dat(i,j,k) - 1.d0*cc_dat(i,j+1,k))/(12.d0)

end do
end do
end do

else
do i = cc_lo(1),cc_hi(1)
do j = cc_lo(2),cc_hi(2)
do k = lo(3),hi(3) + 1

face_dat(i,j,k) = (-1.d0*cc_dat(i,j,k-2) + 7.d0*cc_dat(i,j,k-1) &
+ 7.d0*cc_dat(i,j,k) - 1.d0*cc_dat(i,j,k+1))/(12.d0)

end do
end do
end do


endif






end subroutine cc_to_face_loc





subroutine fill_bdry_values(lo, hi, phi, philo, phihi, dx, prob_lo, prob_hi,time, epsilon, k_freq, kappa,Nprob)bind(C, name="fill_bdry_values")
!  Initialize the scalar field phi
use amrex_fort_module, only : amrex_real

implicit none

integer, intent(in) :: lo(2), hi(2), philo(2), phihi(2)
real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
real(amrex_real), intent(in   ) :: dx(2)
real(amrex_real), intent(in   ) :: prob_lo(2)
real(amrex_real), intent(in   ) :: prob_hi(2)
real(amrex_real), intent(in   ) :: k_freq, time, kappa, epsilon
integer, intent(in) :: Nprob

integer          :: i,j, i_quad, j_quad, face_index
double precision :: x,y,tupi,t0,d, pi, x_quad, y_quad
double precision :: gauss_nodeFrac(0:4)
double precision :: gauss_weights(0:4)
!gauss_nodeFrac = (/ (1.d0-(3.d0/5.d0)**(0.5d0))/2.d0,0.5d0,(1.d0+(3.d0/5.d0)**(0.5d0))/2.d0 /)
!gauss_weights = (/ (5.d0/18.d0),(8.d0/18.d0),(5.d0/18.d0)/)
tupi=3.14159265358979323846d0*2d0
pi=3.14159265358979323846d0
gauss_nodeFrac = (/ (1.d0 - ((1.d0/3.d0)*sqrt(5.d0+2.d0*sqrt(10.d0/7.d0))))/2.d0, &
(1.d0 - ((1.d0/3.d0)*sqrt(5.d0-2.d0*sqrt(10.d0/7.d0))))/2.d0, &
0.5d0,&
(1.d0 + ((1.d0/3.d0)*sqrt(5.d0-2.d0*sqrt(10.d0/7.d0))))/2.d0, &
(1.d0 + ((1.d0/3.d0)*sqrt(5.d0+2.d0*sqrt(10.d0/7.d0))))/2.d0 /)

gauss_weights = (/(322.d0-13.d0*sqrt(70.d0))/1800.d0 , &
(322.d0+13.d0*sqrt(70.d0))/1800.d0, &
128.d0/450.d0, &
(322.d0+13.d0*sqrt(70.d0))/1800.d0, &
(322.d0-13.d0*sqrt(70.d0))/1800.d0/)
 ! We fill ghost cells with Dir values

! x faces

do i = 0,1
if (i==0) then
    face_index = philo(1)
    x = prob_lo(1)
else
    face_index = phihi(1)
    x = prob_hi(1)
end if
    do j = lo(2), hi(2)
    !y = prob_lo(2) + (dble(j)+1.d0/2.d0 )* dx(2)
    y = prob_lo(2) + dble(j)* dx(2)

        phi(face_index,j) = 0.d0
            do j_quad = 0,4
            y_quad = y + dx(2)*gauss_nodeFrac(j_quad)
            !print*, y_quad

                phi(face_index,j) = phi(face_index,j)+ gauss_weights(j_quad)* &
                    cos(k_freq*x)*cos(k_freq*y_quad)

            end do
            !phi(face_index,j) = 0.d0  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end do
end do


! y faces
do j = 0, 1
!y = prob_lo(2) + (dble(j)+1.d0/2.d0 )* dx(2)
if (j==0) then
face_index = philo(2)
y = prob_lo(2)
else
face_index = phihi(2)
y = prob_hi(2)
end if
    do i = lo(1), hi(1)
    !x = prob_lo(1) + (dble(i)+1.d0/2.d0) * dx(1)
    x = prob_lo(1) + dble(i) * dx(1)

    phi(i,face_index) = 0.d0

            do i_quad = 0,4
            x_quad = x + dx(1)*gauss_nodeFrac(i_quad)
            phi(i,face_index) = phi(i,face_index)+ gauss_weights(i_quad)* &
                cos(k_freq*x_quad)*cos(k_freq*y)

            end do

          !  phi(i,face_index) = 0.d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !phi(i,j)=cos(pi*(x+y))
    ! phi(i,j) = 0;
    end do
end do



end subroutine fill_bdry_values
