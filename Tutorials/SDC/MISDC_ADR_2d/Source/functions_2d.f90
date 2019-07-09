
subroutine SDC_feval_F(lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi,f, flo,fhi, &
                         dx,a,d,r,facex, facexlo, facexhi,facey, faceylo, faceyhi,prodx, prodxlo, prodxhi,prody, prodylo, prodyhi,n) bind(C, name="SDC_feval_F")
 ! facex, facexlo, facexhi,facey, faceylo, faceyhi
  !  Compute the rhs terms
  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2),flo(2), fhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(inout) :: f   (flo(1):fhi(1),flo(2):fhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), intent(in)    :: a,d,r
  integer, intent(in) :: n
 integer facexlo(2), facexhi(2)
real(amrex_real), intent(in) :: facex( facexlo(1): facexhi(1), facexlo(2): facexhi(2))
 integer faceylo(2), faceyhi(2)
real(amrex_real), intent(in) :: facey( faceylo(1): faceyhi(1), faceylo(2): faceyhi(2))
integer prodxlo(2), prodxhi(2)
real(amrex_real), intent(inout) :: prodx( prodxlo(1): prodxhi(1), prodxlo(2): prodxhi(2))
integer prodylo(2), prodyhi(2)
real(amrex_real), intent(inout) :: prody( prodylo(1): prodyhi(1), prodylo(2): prodyhi(2))
  
  ! local variables
  integer i,j
real(amrex_real) :: phi_one, face_one
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

        !  Function value is divergence of flux
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              f(i,j) =  a*((fluxx(i+1,j  ) - fluxx(i,j))/dx(1) &
                   + (fluxy(i  ,j+1) - fluxy(i,j))/dx(2))
           end do
           
        end do
     case (1)  !  First implicit piece (here it is diffusion)
     ! We let beta vary in space via cell centered data bcc.





        ! x-fluxes
        do j = philo(2), phihi(2)
           do i = lo(1), hi(1)+1
              ! fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
              fluxx(i,j) = ( -phi(i+1,j)  +15.0d0*(phi(i,j) - phi(i-1,j)) + phi(i-2,j) ) /(12.0d0*dx(1))
           end do
        end do
        
        ! y-fluxes
        do j = lo(2), hi(2)+1
           do i = philo(1), phihi(1)
              ! fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx(2)
              fluxy(i,j) = ( -phi(i,j+1)  +15.0d0*(phi(i,j) - phi(i,j-1)) + phi(i,j-2) ) /(12.0d0*dx(2))
           end do
        end do

       ! We now apply the product rule to find 4th order approximation
        ! x-flux products
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)+1
               ! no dx in one variables!
               phi_one = (-5.d0*fluxx(i,j+2)+34.d0*(fluxx(i,j+1)-fluxx(i,j-1))+5*fluxx(i,j-2))/48.d0
               face_one = (-5.d0*facex(i,j+2)+34.d0*(facex(i,j+1)-facex(i,j-1))+5*facex(i,j-2))/48.d0
               prodx(i,j) = fluxx(i,j)*facex(i,j) + dx(1)*phi_one*face_one/12.d0
            end do
        end do

        ! y-fluxes
        do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)
                ! no dx in one variables!
                phi_one = (-5.d0*fluxy(i+2,j)+34.d0*(fluxy(i+1,j)-fluxy(i-1,j))+5*fluxy(i-2,j))/48.d0
                face_one = (-5.d0*facey(i+2,j)+34.d0*(facey(i+1,j)-facey(i-1,j))+5*facey(i-2,j))/48.d0
                prody(i,j) = fluxy(i,j)*facey(i,j) + dx(2)*phi_one*face_one/12.d0
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
           do i = lo(1), hi(1)
!              f(i,j) =  r*phi(i,j)*(1.0d0-phi(i,j))*(0.5d0-phi(i,j))
              f(i,j) =  -r*phi(i,j)
             end do
          end do
     case default
        print *, 'bad case in advance_2d'
     end select
     
   end subroutine SDC_feval_F
   

   subroutine SDC_fcomp_reaction_F (lo, hi, domlo, domhi, phi, philo, phihi, &
        rhs, rhslo, rhshi, &
        f, flo,fhi, dtq,n) bind(C, name="SDC_fcomp_reaction_F")

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
  integer, intent(in)             :: n
  
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
        dtq,dx) bind(C, name="SDC_Lresid_F")

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



