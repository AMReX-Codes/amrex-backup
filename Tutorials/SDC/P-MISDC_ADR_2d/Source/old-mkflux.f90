module mkflux_module
  ! MM_mc  module mkflux_module
  ! MM_mc  Make the edge average of the flux function to be used in the conservative
  ! MM_mc  finite-volume form of the nonlinear terms
  ! MM_mc  The main subroutine mkflux proceed by calling mv_vels_to_edges
  ! MM_mc  and then a mkflux_Nd routine to actually compute the flux functions

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use mv_vels_to_edges_module

  implicit none

  private

  public :: mkflux, mkflux_s

contains

  !------------------------------------------------------------------------------
  subroutine mkflux(mla,sold,uold,flux,umac,q,dx,the_bc_level)
  !------------------------------------------------------------------------------
    ! MM_sc subroutine mkflux(mla,sold,uold,flux,umac,q,dx,the_bc_level)
    ! MM_sc Compute the edge average of the flux functions
    ! MM_sc  sold and uold are the cell-average conserved quantities, and umac is the
    ! MM_sc  we also need these quantities at the opposite faces (u at top and v at sides)
    ! MM_sc  These could be computed here using phi if phi was really the gradient
    ! MM_sc  which came from the mac_project on uold to produce umac
   
    use ml_restriction_module, only: ml_edge_restriction_c

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)     ! cell average
    type(multifab) , intent(in   ) :: uold(:)     !  cell average
    type(multifab) , intent(inout) :: flux(:,:)   ! the thing we want to compute
    type(multifab) , intent(in   ) :: umac(:,:)   ! the div=0 edge ave. vel.
    type(multifab) , intent(in   ) :: q(:)        ! MM adding to flux
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    type(multifab), allocatable :: sedge(:,:)
    type(multifab), allocatable :: qedge(:,:)

    integer                  :: n,i,dm,ng,comp,ncomp,bccomp,nlevs
    integer                  :: ng_edge,ng_u
    integer                  :: lo(sold(1)%dim),hi(sold(1)%dim)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: uop(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: qepx(:,:,:,:)
    real(kind=dp_t), pointer :: qepy(:,:,:,:)
    real(kind=dp_t), pointer :: qepz(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpx(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpy(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpz(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)

    nlevs = mla%nlevel
    dm    = mla%dim

    ncomp = multifab_ncomp(sold(1))

    bccomp = 1

    allocate(sedge(nlevs,dm))
    allocate(qedge(nlevs,dm))

    ng_edge = 2

    do n = 1, nlevs
       do i = 1,dm
         call multifab_build_edge(sedge(n,i),mla%la(n),dm,ng_edge,i)
         call multifab_build_edge(qedge(n,i),mla%la(n), 1,ng_edge,i)
         call setval(sedge(n,i),ZERO,all=.true.)
         call setval(qedge(n,i),ZERO,all=.true.)
       end do
    enddo

    call mv_vels_to_edges(mla,sold,sedge,the_bc_level,all=.true.) ! all=.true. indicates all 4 edges
    call mv_vels_to_edges(mla,q   ,qedge,the_bc_level,all=.true.) ! all=.true. indicates all 4 edges

    ng   = sold(1)%ng
    ng_u = umac(1,1)%ng

    do n=1,nlevs
       do i = 1, nfabs(sold(n))
          sop    => dataptr(sold(n), i)
          uop    => dataptr(uold(n), i)
          sepx   => dataptr(sedge(n,1), i)
          sepy   => dataptr(sedge(n,2), i)
          qepx   => dataptr(qedge(n,1), i)
          qepy   => dataptr(qedge(n,2), i)
          fluxpx => dataptr(flux(n,1), i)
          fluxpy => dataptr(flux(n,2), i)
          ump    => dataptr(umac(n,1), i)
          vmp    => dataptr(umac(n,2), i)
          lo = lwb(get_box(uold(n), i))
          hi = upb(get_box(uold(n), i))

          select case (dm)
             case (2)
                call mkflux_2d(sop(:,:,1,:), &
                               sepx(:,:,1,:), sepy(:,:,1,:), &
                               qepx(:,:,1,1), qepy(:,:,1,1), ng_edge, &
                               fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), ng_u, &
                               lo, dx(n,:),  &
                               the_bc_level(n)%phys_bc_level_array(i,:,:), &
                               the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:bccomp+ncomp-1),&
                               ng)
             case (3)
                sepz   => dataptr(sedge(n,3), i)
                qepz   => dataptr(qedge(n,3), i)
                fluxpz => dataptr( flux(n,3), i)
                wmp    => dataptr( umac(n,3), i)
                call mkflux_3d(sop(:,:,:,:), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               qepx(:,:,:,1), qepy(:,:,:,1), qepz(:,:,:,1), ng_edge, &
                               fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_u, &
                               lo, dx(n,:),  &
                               the_bc_level(n)%phys_bc_level_array(i,:,:), &
                               the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:bccomp+ncomp-1),&
                               ng)
          end select
       end do
      enddo

    do n = 1, nlevs
       do i = 1,dm
          call multifab_destroy(sedge(n,i))
          call multifab_destroy(qedge(n,i))
       end do
    end do

    do n = nlevs,2,-1
       do comp = 1, ncomp
          do i = 1, dm
             call ml_edge_restriction_c(flux(n-1,i),comp,flux(n,i),comp,mla%mba%rr(n-1,:),i,1)
          enddo
       enddo
    enddo
    
  end subroutine mkflux

  !------------------------------------------------------------------------------
  subroutine mkflux_s(mla,sold,uold,flux,umac,dx,the_bc_level)
  !------------------------------------------------------------------------------
    ! MM_sc subroutine mkflux_s(mla,sold,uold,sedge,flux,umac,dx,the_bc_level)
    ! MM_sc Compute the edge average of the flux functions
    ! MM_sc  sold and uold are the cell-average conserved quantities, and umac is the
    ! MM_sc  divergence free velocity used for the update 
    ! MM_sc  We pass in umac which should be the mac-div=0 velocities, but 
    ! MM_sc  we also need these quantities at the opposite faces (u at top and v at sides)
    ! MM_sc  These could be computed here using phi if phi was really the gradient
    ! MM_sc  which came from the mac_project on uold to produce umac
   
    use ml_restriction_module, only: ml_edge_restriction_c
    use probin_module, only        : nscal

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)   ! cell average
    type(multifab) , intent(in   ) :: uold(:)   !  cell average
    type(multifab) , intent(inout) :: flux(:,:)   ! the thing we want to compute
    type(multifab) , intent(in   ) :: umac(:,:)   ! the div=0 edge ave. vel.
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    type(multifab), allocatable :: sedge(:,:)

    integer                  :: n,i,dm,ng,comp,ncomp,bccomp,nlevs
    integer                  :: ng_edge,ng_u
    integer                  :: lo(sold(1)%dim),hi(sold(1)%dim)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: uop(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpx(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpy(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpz(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)

    nlevs = mla%nlevel
    dm    = mla%dim

    ncomp = multifab_ncomp(sold(1))
    bccomp = dm+1

    allocate(sedge(nlevs,dm))

    ng_edge = 2

    do n = 1, nlevs
       do i = 1,dm
         call multifab_build_edge(sedge(n,i),mla%la(n),nscal,ng_edge,i)
         call setval(sedge(n,i),ZERO,all=.true.)
       end do
    enddo

    ! Setting all = .true. means put all components of sold onto all edges
    call mv_vels_to_edges(mla,sold,sedge,the_bc_level,all=.true.) 
    
    ng   = sold(1)%ng
    ng_u = umac(1,1)%ng

    do n=1,nlevs
       do i = 1, nfabs(sold(n))
          sop    => dataptr(sold(n), i)
          uop    => dataptr(uold(n), i)
          sepx   => dataptr(sedge(n,1), i)
          sepy   => dataptr(sedge(n,2), i)
          fluxpx => dataptr(flux(n,1), i)
          fluxpy => dataptr(flux(n,2), i)
          ump    => dataptr(umac(n,1), i)
          vmp    => dataptr(umac(n,2), i)
          lo = lwb(get_box(uold(n), i))
          hi = upb(get_box(uold(n), i))

          select case (dm)
             case (2)
                call mkflux_s_2d(sop(:,:,1,:), &
                                     sepx(:,:,1,:), sepy(:,:,1,:), ng_edge, &
                                     fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                                     ump(:,:,1,1), vmp(:,:,1,1), ng_u, &
                                     lo, dx(n,:), &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                     the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:bccomp+ncomp-1),&
                                     ng)
             case (3)
                sepz   => dataptr(sedge(n,3), i)
                fluxpz => dataptr(flux(n,3), i)
                wmp  => dataptr(umac(n,3), i)
                call mkflux_s_3d(sop(:,:,:,:), &
                                 sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), ng_edge,&
                                 fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                                 ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_u, &
                                 lo, dx(n,:), &
                                 the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                 the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:bccomp+ncomp-1),&
                                 ng)
          end select
       end do
      enddo

    do n = 1, nlevs
       do i = 1,dm
          call multifab_destroy(sedge(n,i))
       end do
    end do

    do n = nlevs,2,-1
       do comp = 1, ncomp
          do i = 1, dm
             call ml_edge_restriction_c(flux(n-1,i),comp,flux(n,i),comp,mla%mba%rr(n-1,:),i,1)
          enddo
       enddo
    enddo
    
  end subroutine mkflux_s

  !-------------------------------------------------------------------------
  subroutine mkflux_2d(s,sedgex,sedgey,qedgex,qedgey,ng_e,fluxx,fluxy,umac,vmac,ng_u,lo,dx, &
                       phys_bc,adv_bc,ng)
  !-------------------------------------------------------------------------
    ! MM_sc  subroutine mkflux_2d(s,sedgex,sedgey,fluxx,fluxy,umac,vmac,ng_u,lo,dx, &
    ! MM_sc                       phys_bc,adv_bc,ng)
    ! MM_sc  Compute the flux functions from the edge average quantitities

    use bc_module
    use bl_constants_module
    use probin_module, only: spatial_order

    integer, intent(in) :: lo(:),ng,ng_e,ng_u

    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: sedgex(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: sedgey(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: qedgex(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(in   ) :: qedgey(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(inout) ::  fluxx(lo(1)     :,lo(2)     :,:)
    real(kind=dp_t), intent(inout) ::  fluxy(lo(1)     :,lo(2)     :,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)-ng_u:,lo(2)-ng_u:)

    real(kind=dp_t),intent(in) :: dx(:)
    integer        ,intent(in) :: phys_bc(:,:)
    integer        ,intent(in) :: adv_bc(:,:,:)

    ! Local variables
    real(kind=dp_t) hx,hy,hxinv,hyinv
    real(kind=dp_t) u_slopey,v_slopex,su_slopey,sv_slopex
    real(kind=dp_t) s1,s2,facx,facy

    integer :: hi(2)
    integer :: i,j,is,js,ie,je,n
    integer :: ncomp

    ncomp = size(s,dim=3)

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    hx    = dx(1)
    hy    = dx(2)
    hxinv = 1.0d0/dx(1)
    hyinv = 1.0d0/dx(2)

    !  Coefficients for slope computation
    s1 = 34.d0/48.d0
    s2 = -5.d0/48.d0 

    facx = ((hx*hx)/TWELVE)
    facy = ((hy*hy)/TWELVE)

    if (spatial_order .eq. 2) then

    ! SECOND ORDER

       do j=js,je
          do i=is,ie+1
             do n = 1,ncomp
                fluxx(i,j,n) = sedgex(i,j,n)*umac(i,j) 
             enddo
          enddo
       enddo

       do j=js,je+1
          do i=is,ie
             do n = 1,ncomp
                fluxy(i,j,n) = sedgey(i,j,n)*vmac(i,j) 
             enddo
          enddo
       enddo
          
    else

    ! FOURTH ORDER

       ! loop over components with fourth-order stuff -- see KMM eqn. 68
       do j=js,je
          do i=is,ie+1
             u_slopey = (s1*(umac(i,j+1)-umac(i,j-1))&
                       + s2*(umac(i,j+2)-umac(i,j-2)) )*hyinv
             do n = 1,ncomp
                fluxx(i,j,n) = sedgex(i,j,n)*umac(i,j) 
                su_slopey = (s1*(sedgex(i,j+1,n)-sedgex(i,j-1,n))&
                           + s2*(sedgex(i,j+2,n)-sedgex(i,j-2,n)) )*hyinv
                fluxx(i,j,n) = fluxx(i,j,n) + facy*u_slopey*su_slopey
             enddo
          enddo
       enddo
          
       do j=js,je+1
          do i=is,ie
             v_slopex =(s1*(vmac(i+1,j)-vmac(i-1,j)) &
                      + s2*(vmac(i+2,j)-vmac(i-2,j)))*hxinv
             do n = 1,ncomp
                fluxy(i,j,n) = sedgey(i,j,n)*vmac(i,j) 
                sv_slopex =(s1*(sedgey(i+1,j,n)-sedgey(i-1,j,n)) &
                          + s2*(sedgey(i+2,j,n)-sedgey(i-2,j,n)))*hxinv
                fluxy(i,j,n) = fluxy(i,j,n) + facx*v_slopex*sv_slopex
             enddo
          enddo
       enddo ! end loop over components

    end if
    
    !  Add in the pressure term
    do j=js,je
       do i=is,ie+1
          fluxx(i,j,1) = fluxx(i,j,1) + qedgex(i,j)
       enddo
    enddo
    
    do j=js,je+1
       do i=is,ie
          fluxy(i,j,2) = fluxy(i,j,2) + qedgey(i,j)
       enddo
    enddo
       
  end subroutine mkflux_2d

  !-------------------------------------------------------------------------
  subroutine mkflux_3d(s,sedgex,sedgey,sedgez,qedgex,qedgey,qedgez,ng_e,&
                       fluxx,fluxy,fluxz,umac,vmac,wmac,ng_u,lo,dx, &
                       phys_bc,adv_bc,ng)
  !-------------------------------------------------------------------------
    ! MM_sc  subroutine mkflux_3d(s,sedgex,sedgey,sedgez,qedgex,qedgey,qedgez,ng_e,&
    ! MM_sc                       fluxx,fluxy,fluxz,&
    ! MM_sc                       umac,vmac,wmac,ng_u,lo,dx, &
    ! MM_sc                       phys_bc,adv_bc,ng)
    ! MM_sc  Compute the flux functions from the edge average quantitities

    use bc_module
    use bl_constants_module
    use probin_module, only: spatial_order

    integer, intent(in) :: lo(:),ng,ng_e,ng_u

    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng  :,lo(2)-ng  :,lo(3)-ng  :,:)
    real(kind=dp_t), intent(in   ) :: sedgex(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: sedgey(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: sedgez(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: qedgex(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: qedgey(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: qedgez(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) ::  fluxx(lo(1)     :,lo(2)     :,lo(3)     :,:)
    real(kind=dp_t), intent(inout) ::  fluxy(lo(1)     :,lo(2)     :,lo(3)     :,:)
    real(kind=dp_t), intent(inout) ::  fluxz(lo(1)     :,lo(2)     :,lo(3)     :,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) ::   wmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)

    real(kind=dp_t),intent(in) :: dx(:)
    integer        ,intent(in) :: phys_bc(:,:)
    integer        ,intent(in) :: adv_bc(:,:,:)

    ! Local variables
    real(kind=dp_t) hx,hy,hz,hxinv,hyinv,hzinv
    real(kind=dp_t)  u_slopey, u_slopez, v_slopex, v_slopez, w_slopex, w_slopey
    real(kind=dp_t) su_slopey,su_slopez,sv_slopex,sv_slopez,sw_slopex,sw_slopey
    real(kind=dp_t) s1,s2,facx,facy,facz

    integer :: hi(3)
    integer :: i,j,k,is,js,ks,ie,je,ke,n
    integer :: ncomp

    ncomp = size(s,dim=4)

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
    hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

    is = lo(1)
    js = lo(2)
    ks = lo(3)
    ie = hi(1)
    je = hi(2)
    ke = hi(3)

    hx    = dx(1)
    hy    = dx(2)
    hz    = dx(3)
    hxinv = 1.0d0/dx(1)
    hyinv = 1.0d0/dx(2)
    hzinv = 1.0d0/dx(3)

    !  Coefficients for slope computation
    s1 = 34.d0/48.d0
    s2 = -5.d0/48.d0 

    facx = ((hx*hx)/TWELVE)
    facy = ((hy*hy)/TWELVE)
    facz = ((hz*hz)/TWELVE)

    if (spatial_order .eq. 2) then

    ! SECOND ORDER

       do k=ks,ke
       do j=js,je
       do i=is,ie+1
          do n = 1,ncomp
             fluxx(i,j,k,n) = sedgex(i,j,k,n)*umac(i,j,k) 
          enddo
       enddo
       enddo
       enddo
       
       do k=ks,ke
       do j=js,je+1
       do i=is,ie
          do n = 1,ncomp
             fluxy(i,j,k,n) = sedgey(i,j,k,n)*vmac(i,j,k) 
          enddo
       enddo
       enddo 
       enddo 
       
       do k=ks,ke+1
       do j=js,je
       do i=is,ie
          do n = 1,ncomp
             fluxz(i,j,k,n) = sedgez(i,j,k,n)*wmac(i,j,k) 
          enddo
       enddo
       enddo 
       enddo 

    else

    ! FOURTH ORDER

       do k=ks,ke
       do j=js,je
       do i=is,ie+1
          u_slopey = (s1*(umac(i,j+1,k)-umac(i,j-1,k))&
                    + s2*(umac(i,j+2,k)-umac(i,j-2,k)) )*hyinv
          u_slopez = (s1*(umac(i,j,k+1)-umac(i,j,k-1))&
                    + s2*(umac(i,j,k+2)-umac(i,j,k-2)) )*hzinv
          do n = 1,ncomp
             fluxx(i,j,k,n) = sedgex(i,j,k,n)*umac(i,j,k) 
             su_slopey = (s1*(sedgex(i,j+1,k,n)-sedgex(i,j-1,k,n))&
                        + s2*(sedgex(i,j+2,k,n)-sedgex(i,j-2,k,n)) )*hyinv
             su_slopez = (s1*(sedgex(i,j,k+1,n)-sedgex(i,j,k-1,n))&
                        + s2*(sedgex(i,j,k+2,n)-sedgex(i,j,k-2,n)) )*hzinv
             fluxx(i,j,k,n) = fluxx(i,j,k,n) + facy*u_slopey*su_slopey &
                                             + facz*u_slopez*su_slopez
          enddo
       enddo
       enddo
       enddo
       
       do k=ks,ke
       do j=js,je+1
       do i=is,ie
          v_slopex =(s1*(vmac(i+1,j,k)-vmac(i-1,j,k)) &
                   + s2*(vmac(i+2,j,k)-vmac(i-2,j,k)))*hxinv
          v_slopez =(s1*(vmac(i,j,k+1)-vmac(i,j,k-1)) &
                   + s2*(vmac(i,j,k+2)-vmac(i,j,k-2)))*hzinv
          do n = 1,ncomp
             fluxy(i,j,k,n) = sedgey(i,j,k,n)*vmac(i,j,k) 
             sv_slopex =(s1*(sedgey(i+1,j,k,n)-sedgey(i-1,j,k,n)) &
                       + s2*(sedgey(i+2,j,k,n)-sedgey(i-2,j,k,n)))*hxinv
             sv_slopez =(s1*(sedgey(i,j,k+1,n)-sedgey(i,j,k-1,n)) &
                       + s2*(sedgey(i,j,k+2,n)-sedgey(i,j,k-2,n)))*hzinv
             fluxy(i,j,k,n) = fluxy(i,j,k,n) + facx*v_slopex*sv_slopex &
                                             + facz*v_slopez*sv_slopez
          enddo
       enddo
       enddo 
       enddo 
       
       do k=ks,ke+1
       do j=js,je
       do i=is,ie
          w_slopex =(s1*(wmac(i+1,j,k)-wmac(i-1,j,k)) &
                   + s2*(wmac(i+2,j,k)-wmac(i-2,j,k)))*hxinv
          w_slopey =(s1*(wmac(i,j+1,k)-wmac(i,j-1,k)) &
                   + s2*(wmac(i,j+2,k)-wmac(i,j-2,k)))*hyinv
          do n = 1,ncomp
             fluxz(i,j,k,n) = sedgez(i,j,k,n)*wmac(i,j,k) 
             sw_slopex =(s1*(sedgez(i+1,j,k,n)-sedgez(i-1,j,k,n)) &
                       + s2*(sedgez(i+2,j,k,n)-sedgez(i-2,j,k,n)))*hxinv
             sw_slopey =(s1*(sedgez(i,j+1,k,n)-sedgez(i,j-1,k,n)) &
                       + s2*(sedgez(i,j+2,k,n)-sedgez(i,j-2,k,n)))*hzinv
             fluxz(i,j,k,n) = fluxz(i,j,k,n) + facx*w_slopex*sw_slopex &
                                             + facz*w_slopey*sw_slopey
          enddo
       enddo
       enddo 
       enddo 

    end if
    
    !  Add in the pressure term
    do k=ks,ke
    do j=js,je
       do i=is,ie+1
          fluxx(i,j,k,1) = fluxx(i,j,k,1) + qedgex(i,j,k)
       enddo
    enddo
    enddo
    
    do k=ks,ke
    do j=js,je+1
       do i=is,ie
          fluxy(i,j,k,2) = fluxy(i,j,k,2) + qedgey(i,j,k)
       enddo
    enddo
    enddo
    
    do k=ks,ke+1
    do j=js,je
       do i=is,ie
          fluxz(i,j,k,3) = fluxz(i,j,k,3) + qedgez(i,j,k)
       enddo
    enddo
    enddo
       
  end subroutine mkflux_3d

  !-------------------------------------------------------------------------
  subroutine mkflux_s_2d(s,sedgex,sedgey,ng_e,fluxx,fluxy,umac,vmac,ng_u,lo,dx, &
                         phys_bc,adv_bc,ng)
  !-------------------------------------------------------------------------
    ! MM_sc  subroutine mkflux_s_2d(s,sedgex,sedgey,fluxx,fluxy,umac,vmac,lo,dx, &
    ! MM_sc                         phys_bc,adv_bc,ng)
    ! MM_sc  Compute the flux functions from the edge average quantitities

    use bc_module
    use bl_constants_module
    use probin_module, only: spatial_order

    integer, intent(in) :: lo(:),ng,ng_e,ng_u

    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t), intent(inout) ::  fluxx(lo(1)     :,lo(2)     :,:)
    real(kind=dp_t), intent(inout) ::  fluxy(lo(1)     :,lo(2)     :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-ng_u:,lo(2)-ng_u:)

    real(kind=dp_t),intent(in) :: dx(:)
    integer        ,intent(in) :: phys_bc(:,:)
    integer        ,intent(in) :: adv_bc(:,:,:)

    real(kind=dp_t) hx, hy, hxinv, hyinv
    real(kind=dp_t) s1, s2, facx, facy

    integer :: hi(2)
    integer :: i,j,is,js,ie,je,n
    integer :: ncomp

    real(kind=dp_t) :: u_slopey,v_slopex,su_slopey,sv_slopex

    ncomp = size(s,dim=3)

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    hx    = dx(1)
    hy    = dx(2)
    hxinv = 1.0d0/dx(1)
    hyinv = 1.0d0/dx(2)

    !  Coefficients for slope computation
    s1 = 34.d0/48.d0
    s2 = -5.d0/48.d0 

    facx = ((hx*hx)/TWELVE)
    facy = ((hy*hy)/TWELVE)

    if (spatial_order .eq. 2) then

    ! SECOND ORDER

       !  Flux at bottom of cell          
       do j=js,je
       do i=is,ie+1
         do n = 1,ncomp
             fluxx(i,j,n) = sedgex(i,j,n)*umac(i,j)
         end do
       enddo
       enddo
       
       !  Flux at top of cell          
       do j=js,je+1
       do i=is,ie
         do n = 1,ncomp
             fluxy(i,j,n) = sedgey(i,j,n)*vmac(i,j)
         end do
       enddo
       enddo

    else

    ! FOURTH ORDER

       !  Flux at bottom of cell          
       do j=js,je
       do i=is,ie+1
         u_slopey = (s1*(umac(i,j+1)-umac(i,j-1))&
                   + s2*(umac(i,j+2)-umac(i,j-2)) )*hyinv
         do n = 1,ncomp
             fluxx(i,j,n) = sedgex(i,j,n)*umac(i,j)
             su_slopey = (s1*(sedgex(i,j+1,n)-sedgex(i,j-1,n))&
                        + s2*(sedgex(i,j+2,n)-sedgex(i,j-2,n)) )*hyinv
             fluxx(i,j,n) = fluxx(i,j,n) + facy*u_slopey*su_slopey
         end do
       enddo
       enddo
       
       !  Flux at top of cell          
       do j=js,je+1
       do i=is,ie
         v_slopex =(s1*(vmac(i+1,j)-vmac(i-1,j)) &
                  + s2*(vmac(i+2,j)-vmac(i-2,j)))*hxinv
         do n = 1,ncomp
             fluxy(i,j,n) = sedgey(i,j,n)*vmac(i,j)
             sv_slopex =(s1*(sedgey(i+1,j,n)-sedgey(i-1,j,n)) &
                       + s2*(sedgey(i+2,j,n)-sedgey(i-2,j,n)))*hxinv
             fluxy(i,j,n) = fluxy(i,j,n) + facx*v_slopex*sv_slopex
         end do
       enddo
       enddo

    end if

  end subroutine mkflux_s_2d

  !-------------------------------------------------------------------------
  subroutine mkflux_s_3d(s,sedgex,sedgey,sedgez,ng_e,&
                         fluxx,fluxy,fluxz,umac,vmac,wmac,ng_u,lo,dx, &
                         phys_bc,adv_bc,ng)
  !-------------------------------------------------------------------------
    ! MM_sc  subroutine mkflux_s_3d(s,sedgex,sedgey,sedgez,fluxx,fluxy,fluxz,&
    ! MM_sc                         umac,vmac,wmac,ng_u,lo,dx, &
    ! MM_sc                         phys_bc,adv_bc,ng)
    ! MM_sc  Compute the flux functions from the edge average quantitities

    use bc_module
    use bl_constants_module
    use probin_module, only: spatial_order

    integer, intent(in) :: lo(:),ng,ng_e,ng_u

    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng  :,lo(2)-ng  :,lo(3)-ng  :,:)
    real(kind=dp_t), intent(in   ) :: sedgex(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: sedgey(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) :: sedgez(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(inout) ::  fluxx(lo(1)     :,lo(2)     :,lo(3)     :,:)
    real(kind=dp_t), intent(inout) ::  fluxy(lo(1)     :,lo(2)     :,lo(3)     :,:)
    real(kind=dp_t), intent(inout) ::  fluxz(lo(1)     :,lo(2)     :,lo(3)     :,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) ::   wmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)

    real(kind=dp_t),intent(in) :: dx(:)
    integer        ,intent(in) :: phys_bc(:,:)
    integer        ,intent(in) :: adv_bc(:,:,:)

    ! Local variables
    real(kind=dp_t) hx,hy,hz,hxinv,hyinv,hzinv
    real(kind=dp_t)  u_slopey, u_slopez, v_slopex, v_slopez, w_slopex, w_slopey
    real(kind=dp_t) su_slopey,su_slopez,sv_slopex,sv_slopez,sw_slopex,sw_slopey
    real(kind=dp_t) s1,s2,facx,facy,facz

    integer :: hi(3)
    integer :: i,j,k,is,js,ks,ie,je,ke,n
    integer :: ncomp

    ncomp = size(s,dim=4)

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
    hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)
    ks = lo(3)
    ke = hi(3)

    hx    = dx(1)
    hy    = dx(2)
    hz    = dx(3)
    hxinv = 1.0d0/dx(1)
    hyinv = 1.0d0/dx(2)
    hzinv = 1.0d0/dx(3)

    !  Coefficients for slope computation
    s1 = 34.d0/48.d0
    s2 = -5.d0/48.d0 

    facx = ((hx*hx)/TWELVE)
    facy = ((hy*hy)/TWELVE)
    facz = ((hz*hz)/TWELVE)

    if (spatial_order .eq. 2) then

    ! SECOND ORDER

       do k=ks,ke
       do j=js,je
       do i=is,ie+1
          do n = 1,ncomp
             fluxx(i,j,k,n) = sedgex(i,j,k,n)*umac(i,j,k) 
          enddo
       enddo
       enddo
       enddo
       
       do k=ks,ke
       do j=js,je+1
       do i=is,ie
          do n = 1,ncomp
             fluxy(i,j,k,n) = sedgey(i,j,k,n)*vmac(i,j,k) 
          enddo
       enddo
       enddo 
       enddo 
       
       do k=ks,ke+1
       do j=js,je
       do i=is,ie
          do n = 1,ncomp
             fluxz(i,j,k,n) = sedgez(i,j,k,n)*wmac(i,j,k) 
          enddo
       enddo
       enddo 
       enddo 

    else

    ! FOURTH ORDER

       do k=ks,ke
       do j=js,je
       do i=is,ie+1
          u_slopey = (s1*(umac(i,j+1,k)-umac(i,j-1,k))&
                    + s2*(umac(i,j+2,k)-umac(i,j-2,k)) )*hyinv
          u_slopez = (s1*(umac(i,j,k+1)-umac(i,j,k-1))&
                    + s2*(umac(i,j,k+2)-umac(i,j,k-2)) )*hzinv
          do n = 1,ncomp
             fluxx(i,j,k,n) = sedgex(i,j,k,n)*umac(i,j,k) 
             su_slopey = (s1*(sedgex(i,j+1,k,n)-sedgex(i,j-1,k,n))&
                        + s2*(sedgex(i,j+2,k,n)-sedgex(i,j-2,k,n)) )*hyinv
             su_slopez = (s1*(sedgex(i,j,k+1,n)-sedgex(i,j,k-1,n))&
                        + s2*(sedgex(i,j,k+2,n)-sedgex(i,j,k-2,n)) )*hzinv
             fluxx(i,j,k,n) = fluxx(i,j,k,n) + facy*u_slopey*su_slopey &
                                             + facz*u_slopez*su_slopez
          enddo
       enddo
       enddo
       enddo
       
       do k=ks,ke
       do j=js,je+1
       do i=is,ie
          v_slopex =(s1*(vmac(i+1,j,k)-vmac(i-1,j,k)) &
                   + s2*(vmac(i+2,j,k)-vmac(i-2,j,k)))*hxinv
          v_slopez =(s1*(vmac(i,j,k+1)-vmac(i,j,k-1)) &
                   + s2*(vmac(i,j,k+2)-vmac(i,j,k-2)))*hzinv
          do n = 1,ncomp
             fluxy(i,j,k,n) = sedgey(i,j,k,n)*vmac(i,j,k) 
             sv_slopex =(s1*(sedgey(i+1,j,k,n)-sedgey(i-1,j,k,n)) &
                       + s2*(sedgey(i+2,j,k,n)-sedgey(i-2,j,k,n)))*hxinv
             sv_slopez =(s1*(sedgey(i,j,k+1,n)-sedgey(i,j,k-1,n)) &
                       + s2*(sedgey(i,j,k+2,n)-sedgey(i,j,k-2,n)))*hzinv
             fluxy(i,j,k,n) = fluxy(i,j,k,n) + facx*v_slopex*sv_slopex &
                                             + facz*v_slopez*sv_slopez
          enddo
       enddo
       enddo 
       enddo 
       
       do k=ks,ke+1
       do j=js,je
       do i=is,ie
          w_slopex =(s1*(wmac(i+1,j,k)-wmac(i-1,j,k)) &
                   + s2*(wmac(i+2,j,k)-wmac(i-2,j,k)))*hxinv
          w_slopey =(s1*(wmac(i,j+1,k)-wmac(i,j-1,k)) &
                   + s2*(wmac(i,j+2,k)-wmac(i,j-2,k)))*hyinv
          do n = 1,ncomp
             fluxz(i,j,k,n) = sedgez(i,j,k,n)*wmac(i,j,k) 
             sw_slopex =(s1*(sedgez(i+1,j,k,n)-sedgez(i-1,j,k,n)) &
                       + s2*(sedgez(i+2,j,k,n)-sedgez(i-2,j,k,n)))*hxinv
             sw_slopey =(s1*(sedgez(i,j+1,k,n)-sedgez(i,j-1,k,n)) &
                       + s2*(sedgez(i,j+2,k,n)-sedgez(i,j-2,k,n)))*hzinv
             fluxz(i,j,k,n) = fluxz(i,j,k,n) + facx*w_slopex*sw_slopex &
                                             + facz*w_slopey*sw_slopey
          enddo
       enddo
       enddo 
       enddo 

    end if
    
  end subroutine mkflux_s_3d

end module mkflux_module
