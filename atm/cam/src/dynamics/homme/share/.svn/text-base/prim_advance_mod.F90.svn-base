#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
!#define _DBG_ !DBG
!
!
module prim_advance_mod
  use edge_mod, only : EdgeBuffer_t
  use kinds, only : real_kind, iulog

  implicit none
  private
  public :: prim_advance_exp, prim_advance_si, prim_advance_init, preq_robert3,&
       applyCAMforcing_dynamics, applyCAMforcing, applyCAMforcing_leapfrog, smooth_phis

  type (EdgeBuffer_t) :: edge1
  type (EdgeBuffer_t) :: edge2
  type (EdgeBuffer_t) :: edge3p1

  real (kind=real_kind) :: initialized_for_dt   = 0

  real (kind=real_kind), allocatable :: ur_weights(:)


contains

  subroutine prim_advance_init(integration)
    use edge_mod, only : initEdgeBuffer
    use dimensions_mod, only : nlev
    use control_mod, only : qsplit
    character(len=*)    , intent(in) :: integration 
    integer :: i

    call initEdgeBuffer(edge3p1,3*nlev+1)

    if(integration == 'semi_imp') then
       call initEdgeBuffer(edge1,nlev)
       call initEdgeBuffer(edge2,2*nlev)
    end if

    allocate(ur_weights(qsplit))
    ur_weights(:)=0.0d0

    if(mod(qsplit,2).NE.0)then
       ur_weights(1)=1.0d0/qsplit
       do i=3,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    else
       do i=2,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    endif

  end subroutine prim_advance_init




  subroutine prim_advance_exp(elem, deriv, hvcoord, Flt, hybrid,&
       dt, tl,  nets, nete, compute_diagnostics,n_Q)
    use bndry_mod, only : bndry_exchangev
    use control_mod, only : hypervis_order, prescribed_wind, qsplit, compute_mean_flux, tstep_type
    use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
    use dimensions_mod, only : nv, nlev, nlevp
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, initEdgeBuffer
    use element_mod, only : element_t
    use filter_mod, only : filter_t, preq_filter, prim_filter
    use hybvcoord_mod, only : hvcoord_t
    use hybrid_mod, only : hybrid_t
    use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use time_mod, only : TimeLevel_t, smooth
    use diffusion_mod, only :  prim_diffusion

#ifdef CAM
    use cam_control_mod, only : ideal_phys
#else
    use asp_tests_mod, only : asp_advection_vertical ! _EXTERNAL
#endif

    implicit none

    type (element_t), intent(inout), target   :: elem(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (filter_t)                   :: flt

    type (hybrid_t)    , intent(in):: hybrid

    real (kind=real_kind), intent(in) :: dt
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    integer              , intent(in) :: n_Q    ! use Q at this timelevel
    logical, intent(in)               :: compute_diagnostics

    ! =================
    ! Local
    ! =================
    real (kind=real_kind) ::  dt2, time, dt_vis, eta_ave_w
    real (kind=real_kind) ::  eta_dot_dpdn(nv,nv,nlevp)
    real (kind=real_kind) ::  dp(nv,nv)
    integer :: ie,nm1,n0,np1,nstep,method,qsplit_stage,nflux,k


    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep
    dt2 = 2*dt


!collect mean quantities:
! omega_P        used by CAM physics
! Udp            used by advection, when subcycling for consistenent advection 
! eta_dot_dpdn   used by advection when not subcycling
!
    if(tstep_type==1)then        ! forward in time scheme 
       method=1
       qsplit_stage = mod(nstep,qsplit)
       if (qsplit_stage==0) method=0
       eta_ave_w=ur_weights(qsplit_stage+1)
    else     ! pure leapfrog
       method=1  
       eta_ave_w = 1d0/qsplit
    endif
    if (nstep==0) method=0  ! RK2 for first timestep




    if (1==prescribed_wind) then
       ! fix dynamical variables, skip dynamics
       time=tl%nstep*dt

       do ie=nets,nete
          ! assume most fields are constant in time
          elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,n0)
          elem(ie)%state%lnps(:,:,np1) = elem(ie)%state%lnps(:,:,n0)
          elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,n0)
          elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,n0)
          elem(ie)%derived%div = 0
#ifndef CAM
	  !get eta_dot_dpdn at n0, ps_v is not used or calculated in this routine
          call asp_advection_vertical(time,hvcoord,elem(ie)%state%ps_v(:,:,n0),&
              eta_dot_dpdn)
#endif
          ! accumulate mean fluxes for advection
          elem(ie)%derived%eta_dot_dpdn(:,:,:) = &
               elem(ie)%derived%eta_dot_dpdn(:,:,:) + &
               eta_dot_dpdn(:,:,:)*eta_ave_w

          if(compute_mean_flux==1)then
             ! subcycling code uses a mean flux to advect tracers
             ! for consistency.  But we need the 
             ! mean flux Udp defined so that:
             ! dp(t+1)-dp(t) = -div ( Udp )  - d/dn ( eta dp )
             ! but this would require solving for Udp, 
             do k=1,nlev
                dp(:,:) =&
                      ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0) 

                elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+&
                     eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*dp(:,:)
                elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+&
                     eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*dp(:,:)
             enddo

          endif

       end do
       call t_stopf('prim_advance_exp')
       return
    endif



    ! ==================================
    ! Take timestep
    ! ==================================
    if (method==0) then ! RK2
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
       call compute_and_apply_rhs(np1,n0,n0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,n_Q,0d0)

       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
       call compute_and_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.,n_Q,eta_ave_w)

       dt_vis = dt                      
    endif
    if (method==1) then
       ! regular LF step
       call compute_and_apply_rhs(np1,nm1,n0,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics,n_Q,eta_ave_w)
       dt_vis = dt2  ! dt to use for time-split dissipation
    endif
 

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================
#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
       do ie = nets,nete
          elem(ie)%accum%DIFF(:,:,:,:)=elem(ie)%state%v(:,:,:,:,np1)
          elem(ie)%accum%DIFFT(:,:,:)=elem(ie)%state%T(:,:,:,np1)
       enddo
    endif
#endif



    ! note:time step computes u(t+1)= u(t*) + RHS. 
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    if (tstep_type==1) then  ! forward-in-time
       call advance_hypervis(edge3p1,elem,hvcoord,hybrid,deriv,np1,np1,np1,nets,nete,dt_vis)
    else ! leapfrog
       call advance_hypervis(edge3p1,elem,hvcoord,hybrid,deriv,nm1,n0,np1,nets,nete,dt_vis)
    endif
#ifdef ENERGY_DIAGNOSTICS
    if (compute_diagnostics) then
       do ie = nets,nete
          elem(ie)%accum%DIFF(:,:,:,:)=( elem(ie)%state%v(:,:,:,:,np1) -&
               elem(ie)%accum%DIFF(:,:,:,:) ) / dt_vis
          elem(ie)%accum%DIFFT(:,:,:)=( elem(ie)%state%T(:,:,:,np1) -&
               elem(ie)%accum%DIFFT(:,:,:) ) / dt_vis
       enddo
    endif
#endif
       ! rest of the code/diagnostics are using lnps, so compute:
    do ie = nets,nete
       elem(ie)%state%lnps(:,:,np1)= LOG(elem(ie)%state%ps_v(:,:,np1))
    end do


    call t_stopf('prim_advance_exp')
    end subroutine prim_advance_exp



subroutine prim_advance_si(elem, nets, nete, cg, blkjac, red, &
          refstate, hvcoord, deriv, flt, hybrid, tl, dt)
       use bndry_mod, only : bndry_exchangev
       use cg_mod, only : cg_t, cg_create
       use control_mod, only : filter_freq,debug_level, precon_method
       use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
       use dimensions_mod, only : nv, nlev, nlevp
       use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, initEdgeBuffer
       use element_mod, only : element_t
       use filter_mod, only : filter_t, preq_filter, prim_filter
       use hybvcoord_mod, only : hvcoord_t
       use hybrid_mod, only : hybrid_t
       use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
       use prim_si_ref_mod, only : ref_state_t, set_vert_struct_mat
       use reduction_mod, only : reductionbuffer_ordered_1d_t
       use solver_mod, only : pcg_solver, blkjac_t, blkjac_init
       use time_mod, only : TimeLevel_t, smooth
       use prim_si_mod, only : preq_impsys, preq_bldrhs, preq_backsub, preq_vertadv, &
            preq_omegap, preq_pressure
       use diffusion_mod, only :  prim_diffusion
       use physical_constants, only : kappa, rearth, rgas, cp, cpwater_vapor, rwater_vapor
       use physics_mod, only : virtual_temperature, virtual_specific_heat
#ifdef CAM
       use cam_control_mod, only : ideal_phys
#endif
       implicit none

       integer, intent(in)               :: nets,nete
       type (element_t), intent(inout), target :: elem(:)
       type (blkjac_t)                   :: blkjac(nets:nete)

       type (cg_t)                       :: cg

       type (ReductionBuffer_ordered_1d_t), intent(inout) :: red

       type (ref_state_t), intent(in), target :: refstate
       type (hvcoord_t), intent(in)      :: hvcoord
       type (derivative_t), intent(in)   :: deriv
       type (filter_t), intent(in)       :: flt
       type (hybrid_t), intent(in)       :: hybrid
       type (TimeLevel_t), intent(in)    :: tl
       real(kind=real_kind), intent(in)  :: dt
       real(kind=real_kind)              :: time_adv
       ! ==========================
       ! Local variables...
       ! ==========================

       real(kind=real_kind)                           :: ps0
       real(kind=real_kind)                           :: psref

       real(kind=real_kind), dimension(nv,nv)         :: ps
       real(kind=real_kind), dimension(nv,nv)         :: rps
       real(kind=real_kind), dimension(nv,nv,nlev)    :: rpmid
       real(kind=real_kind), dimension(nv,nv,nlev)    :: omegap
       real(kind=real_kind), dimension(nv,nv,nlev)    :: rpdel

       real(kind=real_kind) :: pintref(nlevp)
       real(kind=real_kind) :: pdelref(nlev)
       real(kind=real_kind) :: pmidref(nlev)
       real(kind=real_kind) :: rpdelref(nlev)
       real(kind=real_kind) :: rpmidref(nlev)

       real(kind=real_kind) :: pint(nv,nv,nlevp)
       real(kind=real_kind) :: pdel(nv,nv,nlev)
       real(kind=real_kind) :: pmid(nv,nv,nlev)

       real(kind=real_kind), dimension(nv,nv,nlevp) :: eta_dot_dp_deta
       real(kind=real_kind), dimension(nv,nv,nlev)  :: vgrad_ps

       real(kind=real_kind), dimension(nv,nv,nlev)   :: T_vadv
       real(kind=real_kind), dimension(nv,nv,2,nlev) :: v_vadv

       real(kind=real_kind), dimension(nv,nv)      :: HT
       real(kind=real_kind), dimension(nv,nv)      :: HrefT
       real(kind=real_kind), dimension(nv,nv)      :: HrefTm1

       real(kind=real_kind), dimension(nv,nv)      :: Gref0
       real(kind=real_kind), dimension(nv,nv)      :: Grefm1
       real(kind=real_kind), dimension(nv,nv)      :: E
       real(kind=real_kind), dimension(nv,nv)      :: Phi
       real(kind=real_kind), dimension(nv,nv)      :: dGref

       real(kind=real_kind), dimension(nv,nv,2)    :: vco
       real(kind=real_kind), dimension(nv,nv,2)    :: gradT
       real(kind=real_kind), dimension(nv,nv,2)    :: grad_Phi

       real(kind=real_kind), dimension(:,:), pointer  :: Emat
       real(kind=real_kind), dimension(:,:), pointer  :: Emat_inv
       real(kind=real_kind), dimension(:,:), pointer  :: Amat
       real(kind=real_kind), dimension(:,:), pointer  :: Amat_inv
       real(kind=real_kind), dimension(:), pointer    :: Lambda

       real(kind=real_kind), dimension(:), pointer    :: Tref
       real(kind=real_kind), dimension(:), pointer    :: RTref
       real(kind=real_kind), dimension(:), pointer    :: Pvec
       real(kind=real_kind), dimension(:,:), pointer  :: Href
       real(kind=real_kind), dimension(:,:), pointer  :: Tmat

       real(kind=real_kind) :: Vscript(nv,nv,2,nlev,nets:nete)
       real(kind=real_kind) :: Tscript(nv,nv,nlev,nets:nete)
       real(kind=real_kind) :: Pscript(nv,nv,nets:nete)

       real(kind=real_kind), dimension(nv,nv)      :: HrefTscript
       real(kind=real_kind), dimension(nv,nv)      :: suml
       real(kind=real_kind), dimension(nv,nv,2)    :: gVscript
       real(kind=real_kind), dimension(nv,nv,nlev) :: div_Vscript

       real(kind=real_kind) :: B(nv,nv,nlev,nets:nete)
       real(kind=real_kind) :: C(nv,nv,nlev,nets:nete)
       real(kind=real_kind) :: D(nv,nv,nlev,nets:nete)

       real(kind=real_kind) :: Gamma_ref(nv,nv,nlev,nets:nete)

       real(kind=real_kind) :: Gref(nv,nv,nlev,nets:nete)
       real(kind=real_kind) :: grad_dGref(nv,nv,2,nlev)
       real(kind=real_kind) :: grad_Gref(nv,nv,2,nlev)

       real(kind=real_kind) :: div(nv,nv)
       real(kind=real_kind) :: gv(nv,nv,2)

       real(kind=real_kind) :: dt2
       real(kind=real_kind) :: rpsref
       real(kind=real_kind) :: rdt
       real(kind=real_kind) :: hkk, hkl
       real(kind=real_kind) :: ddiv

       real(kind=real_kind) :: vgradT
       real(kind=real_kind) :: hybfac
       real(kind=real_kind) :: Crkk
       real(kind=real_kind) :: v1,v2
       real(kind=real_kind) :: term

       real(kind=real_kind) :: rdx,rdy
       real(kind=real_kind) :: dx,dy

       real(kind=real_kind) :: Vs1,Vs2
       real(kind=real_kind) :: glnps1, glnps2
       real(kind=real_kind) :: gGr1,gGr2

       real (kind=real_kind),allocatable :: solver_wts(:,:)  ! solver weights array for nonstaggered grid

       integer              :: nm1,n0,np1,nfilt
       integer              :: nstep
       integer              :: i,j,k,l,ie,kptr

       call t_startf('prim_advance_si')

       nm1   = tl%nm1
       n0    = tl%n0
       np1   = tl%np1
       nstep = tl%nstep


       if ( dt /= initialized_for_dt ) then
          if(hybrid%par%masterproc) print *,'Initializing semi-implicit matricies for dt=',dt

          !$OMP MASTER
          call set_vert_struct_mat(dt, refstate, hvcoord, hybrid%masterthread)
          !$OMP END MASTER

          allocate(solver_wts(nv*nv,nete-nets+1))
          do ie=nets,nete
             kptr=1
             do j=1,nv
                do i=1,nv
                   solver_wts(kptr,ie-nets+1) = elem(ie)%solver_wts(i,j)
                   kptr=kptr+1
                end do
             end do
          end do
          call cg_create(cg, nv*nv, nlev, nete-nets+1, hybrid, debug_level, solver_wts)
          deallocate(solver_wts)
          if (precon_method == "block_jacobi") then
             call blkjac_init(elem, deriv,refstate%Lambda,nets,nete,blkjac)
          end if
          initialized_for_dt = dt
       endif


       nfilt = tl%nm1     ! time level at which filter is applied (time level n)
       dt2   = 2.0_real_kind*dt
       rdt   = 1.0_real_kind/dt

       ps0      = hvcoord%ps0
       psref    = refstate%psr

       Emat     => refstate%Emat
       Emat_inv => refstate%Emat_inv
       Amat     => refstate%Amat
       Amat_inv => refstate%Amat_inv
       Lambda   => refstate%Lambda

       RTref    => refstate%RTref
       Tref     => refstate%Tref
       Href     => refstate%Href
       Tmat     => refstate%Tmat
       Pvec     => refstate%Pvec

       ! ============================================================
       ! If the time is right, apply a filter to the state variables
       ! ============================================================

       if (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0 ) then
          call preq_filter(elem, edge3p1, flt, cg%hybrid, nfilt, nets, nete)
       end if

       do ie = nets, nete

          rdx=2.0_real_kind/(elem(ie)%dx*rearth) ! strong derivative inverse x length
          rdy=2.0_real_kind/(elem(ie)%dy*rearth) ! strong derivative inverse y length

          elem(ie)%derived%grad_lnps(:,:,:) = gradient(elem(ie)%state%lnps(:,:,n0),deriv,rdx,rdy)

       end do

       ! ================================================
       ! boundary exchange grad_lnps
       ! ================================================

       do ie = nets, nete

          dx=0.5_real_kind*elem(ie)%dx/rearth ! weak derivative element x dimension
          dy=0.5_real_kind*elem(ie)%dy/rearth ! weak derivative element y dimension

          rdx=2.0_real_kind/(elem(ie)%dx*rearth) ! strong derivative inverse x length
          rdy=2.0_real_kind/(elem(ie)%dy*rearth) ! strong derivative inverse y length

          do k=1,nlevp
             pintref(k)  = hvcoord%hyai(k)*ps0 + hvcoord%hybi(k)*psref
          end do

          do k=1,nlev
             pmidref(k)  = hvcoord%hyam(k)*ps0 + hvcoord%hybm(k)*psref
             pdelref(k)  = pintref(k+1) - pintref(k)
             rpmidref(k) = 1.0_real_kind/pmidref(k)
             rpdelref(k) = 1.0_real_kind/pdelref(k)
          end do

          rpsref   = 1.0_real_kind/psref

          ps(:,:) = EXP(elem(ie)%state%lnps(:,:,n0))
          rps(:,:) = 1.0_real_kind/ps(:,:)

          call preq_pressure(ps0,ps,hvcoord%hyai,hvcoord%hybi,hvcoord%hyam,hvcoord%hybm,pint,pmid,pdel)

          rpmid = 1.0_real_kind/pmid
          rpdel = 1.0_real_kind/pdel

          do k=1,nlev
             do j=1,nv
                do i=1,nv
                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! Contravariant velocities

                   vco(i,j,1) = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
                   vco(i,j,2) = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2

                   vgrad_ps(i,j,k) = ps(i,j)*(vco(i,j,1)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        vco(i,j,2)*elem(ie)%derived%grad_lnps(i,j,2))

                end do
             end do
          end do

          call preq_omegap(elem(ie)%derived%div(:,:,:,n0),vgrad_ps,pdel,rpmid, &
               hvcoord%hybm,hvcoord%hybd,elem(ie)%derived%omega_p)

          Pscript(:,:,ie)        = 0.0_real_kind
          eta_dot_dp_deta(:,:,1) = 0.0_real_kind

          do k=1,nlev
             do j=1,nv
                do i=1,nv
                   eta_dot_dp_deta(i,j,k+1) = eta_dot_dp_deta(i,j,k) + &
                        vgrad_ps(i,j,k)*hvcoord%hybd(k) + elem(ie)%derived%div(i,j,k,n0)*pdel(i,j,k)
                   ddiv = elem(ie)%derived%div(i,j,k,n0) - 0.5_real_kind*elem(ie)%derived%div(i,j,k,nm1)
                   Pscript(i,j,ie) = Pscript(i,j,ie) + ddiv*pdelref(k)
                end do
             end do
          end do

          do j=1,nv
             do i=1,nv
                Pscript(i,j,ie) = elem(ie)%state%lnps(i,j,nm1) + &
                     dt2*( rpsref*Pscript(i,j,ie) - rps(i,j)*eta_dot_dp_deta(i,j,nlev+1) )
             end do
          end do

          do k=1,nlev-1
             do j=1,nv
                do i=1,nv
                   eta_dot_dp_deta(i,j,k+1) = hvcoord%hybi(k+1)*eta_dot_dp_deta(i,j,nlev+1) - &
                        eta_dot_dp_deta(i,j,k+1)
                end do
             end do
          end do

          eta_dot_dp_deta(:,:,nlev+1) = 0.0_real_kind

          call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
               eta_dot_dp_deta,rpdel,T_vadv,v_vadv)

          suml(:,:) = 0.0_real_kind

          do k=1,nlev

             gradT(:,:,:) = gradient(elem(ie)%state%T(:,:,k,n0),deriv,rdx,rdy)
             Crkk       = 0.5_real_kind

             do j=1,nv
                do i=1,nv
                   term = Crkk*(elem(ie)%derived%div(i,j,k,n0) - &
                        0.5_real_kind*elem(ie)%derived%div(i,j,k,nm1))*pdelref(k)
                   suml(i,j)  = suml(i,j) + term

                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! Contravariant velocities

                   vco(i,j,1) = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
                   vco(i,j,2) = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2

                   vgradT = vco(i,j,1)*gradT(i,j,1) + vco(i,j,2)*gradT(i,j,2)

                   Tscript(i,j,k,ie) = elem(ie)%state%T(i,j,k,nm1) &
                        + dt2*(- vgradT - T_vadv(i,j,k)           &
                        + kappa*(elem(ie)%state%T(i,j,k,n0)*elem(ie)%derived%omega_p(i,j,k) &
                        + Tref(k)*rpmidref(k)*suml(i,j)))
                   suml(i,j)  = suml(i,j) + term
                end do
             end do
          end do

          HrefT(:,:)   = 0.0_real_kind
          HrefTm1(:,:) = 0.0_real_kind
          HT(:,:)      = 0.0_real_kind

          do k=nlev,1,-1

             do j=1,nv
                do i=1,nv
                   hkl = rpmidref(k)*pdelref(k)
                   hkk = hkl*0.5_real_kind
                   Gref0(i,j)   = HrefT(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,n0)
                   HrefT(i,j)   = HrefT(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,n0)
                   Grefm1(i,j)  = HrefTm1(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,nm1)
                   HrefTm1(i,j) = HrefTm1(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,nm1)
                   hkl = rpmid(i,j,k)*pdel(i,j,k)
                   hkk = hkl*0.5_real_kind
                   Phi(i,j) = HT(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,n0)
                   HT(i,j)  = HT(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,n0)
                end do
             end do

             do j=1,nv
                do i=1,nv
                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! covariant velocity

                   vco(i,j,1) = elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(2,1,i,j)*v2
                   vco(i,j,2) = elem(ie)%D(1,2,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2

                   E(i,j) = 0.5_real_kind*( v1*v1 + v2*v2 )

                   Gref0(i,j)  =  Gref0(i,j)  + elem(ie)%state%phis(i,j) + RTref(k)*elem(ie)%state%lnps(i,j,n0)
                   Grefm1(i,j) =  Grefm1(i,j) + elem(ie)%state%phis(i,j) + RTref(k)*elem(ie)%state%lnps(i,j,nm1)

                   Phi(i,j)    =  Phi(i,j) + E(i,j) + elem(ie)%state%phis(i,j)
                   dGref(i,j)  =  -(Gref0(i,j)  - 0.5_real_kind*Grefm1(i,j))
                end do
             end do

             elem(ie)%derived%zeta(:,:,k) = vorticity(vco,deriv,rdx,rdy)
             grad_Phi(:,:,:)     = gradient(Phi,deriv,rdx,rdy)
             grad_dGref(:,:,:,k) = gradient_wk(dGref,deriv,dx,dy)

             do j=1,nv
                do i=1,nv

                   elem(ie)%derived%zeta(i,j,k) = elem(ie)%rmetdetp(i,j)*elem(ie)%derived%zeta(i,j,k)
                   hybfac =  hvcoord%hybm(k)*(ps(i,j)*rpmid(i,j,k))

                   glnps1 = elem(ie)%Dinv(1,1,i,j)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        elem(ie)%Dinv(2,1,i,j)*elem(ie)%derived%grad_lnps(i,j,2)
                   glnps2 = elem(ie)%Dinv(1,2,i,j)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        elem(ie)%Dinv(2,2,i,j)*elem(ie)%derived%grad_lnps(i,j,2)

                   v1 = elem(ie)%Dinv(1,1,i,j)*grad_Phi(i,j,1) + elem(ie)%Dinv(2,1,i,j)*grad_Phi(i,j,2)
                   v2 = elem(ie)%Dinv(1,2,i,j)*grad_Phi(i,j,1) + elem(ie)%Dinv(2,2,i,j)*grad_Phi(i,j,2)

                   Vscript(i,j,1,k,ie) = - v_vadv(i,j,1,k) &
                        + elem(ie)%state%v(i,j,2,k,n0) * (elem(ie)%fcor(i,j) + elem(ie)%derived%zeta(i,j,k)) &
                        - v1 - Rgas*hybfac*elem(ie)%state%T(i,j,k,n0)*glnps1

                   Vscript(i,j,2,k,ie) = - v_vadv(i,j,2,k) &
                        - elem(ie)%state%v(i,j,1,k,n0) * (elem(ie)%fcor(i,j) + elem(ie)%derived%zeta(i,j,k)) &
                        - v2 - Rgas*hybfac*elem(ie)%state%T(i,j,k,n0)*glnps2

                end do
             end do

          end do

          do k=1,nlev

             do j=1,nv
                do i=1,nv
                   Vs1 = elem(ie)%Dinv(1,1,i,j)*grad_dGref(i,j,1,k) + elem(ie)%Dinv(2,1,i,j)*grad_dGref(i,j,2,k)
                   Vs2 = elem(ie)%Dinv(1,2,i,j)*grad_dGref(i,j,1,k) + elem(ie)%Dinv(2,2,i,j)*grad_dGref(i,j,2,k)

                   Vscript(i,j,1,k,ie) = elem(ie)%mv(i,j)*Vscript(i,j,1,k,ie) + Vs1
                   Vscript(i,j,2,k,ie) = elem(ie)%mv(i,j)*Vscript(i,j,2,k,ie) + Vs2

                   Vscript(i,j,1,k,ie) = elem(ie)%mv(i,j)*elem(ie)%state%v(i,j,1,k,nm1) + dt2*Vscript(i,j,1,k,ie)
                   Vscript(i,j,2,k,ie) = elem(ie)%mv(i,j)*elem(ie)%state%v(i,j,2,k,nm1) + dt2*Vscript(i,j,2,k,ie)
                end do
             end do

          end do

          HrefTscript(:,:) = 0.0_real_kind

          do k=nlev,1,-1

             do j=1,nv
                do i=1,nv
                   hkl = rpmidref(k)*pdelref(k)
                   hkk = hkl*0.5_real_kind
                   B(i,j,k,ie)      = HrefTscript(i,j) + Rgas*hkk*Tscript(i,j,k,ie)
                   B(i,j,k,ie)      = B(i,j,k,ie) +  elem(ie)%state%phis(i,j) + RTref(k)*Pscript(i,j,ie)
                   HrefTscript(i,j) = HrefTscript(i,j) + Rgas*hkl*Tscript(i,j,k,ie)
                end do
             end do

          end do

          kptr=0
          call edgeVpack(edge2, Vscript(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

       end do

       call bndry_exchangeV(cg%hybrid,edge2)

       do ie = nets, nete

          rdx=2.0_real_kind/(elem(ie)%dx*rearth) ! strong derivative inverse x length
          rdy=2.0_real_kind/(elem(ie)%dy*rearth) ! strong derivative inverse y length

          kptr=0
          call edgeVunpack(edge2, Vscript(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
          do k=1,nlev

             do j=1,nv
                do i=1,nv
                   Vscript(i,j,1,k,ie) = elem(ie)%rmv(i,j)*Vscript(i,j,1,k,ie)
                   Vscript(i,j,2,k,ie) = elem(ie)%rmv(i,j)*Vscript(i,j,2,k,ie)
                end do
             end do

             do j=1,nv
                do i=1,nv

                   ! Contravariant Vscript

                   gVscript(i,j,1) = elem(ie)%Dinv(1,1,i,j)*Vscript(i,j,1,k,ie) + &
                        elem(ie)%Dinv(1,2,i,j)*Vscript(i,j,2,k,ie)
                   gVscript(i,j,2) = elem(ie)%Dinv(2,1,i,j)*Vscript(i,j,1,k,ie) + &
                        elem(ie)%Dinv(2,2,i,j)*Vscript(i,j,2,k,ie) 

                   gVscript(i,j,1) = elem(ie)%metdetp(i,j)*gVscript(i,j,1)
                   gVscript(i,j,2) = elem(ie)%metdetp(i,j)*gVscript(i,j,2)

                end do
             end do

             div_Vscript(:,:,k) = divergence(gVscript,deriv,rdx,rdy)

          end do

          do k=1,nlev

             do j=1,nv
                do i=1,nv
                   C(i,j,k,ie) = elem(ie)%metdet(i,j)*B(i,j,k,ie)
                end do
             end do

             do l=1,nlev
                do j=1,nv
                   do i=1,nv
                      C(i,j,k,ie) = C(i,j,k,ie) - dt*Amat(l,k)*div_Vscript(i,j,l)
                   end do
                end do
             end do

          end do

          ! ===============================================================
          !  Weight C (the RHS of the helmholtz problem) by the mass matrix
          ! ===============================================================

          do k=1,nlev
             do j=1,nv
                do i=1,nv
                   C(i,j,k,ie) = elem(ie)%mv(i,j)*C(i,j,k,ie)
                end do
             end do
          end do

          ! ===================================
          ! Pack C into the edge1 buffer
          ! ===================================

          kptr=0
          call edgeVpack(edge1,C(1,1,1,ie),nlev,kptr,elem(ie)%desc)

       end do

       ! ==================================
       ! boundary exchange C
       ! ==================================

       call bndry_exchangeV(cg%hybrid,edge1)

       do ie=nets,nete

          ! ===================================
          ! Unpack C from the edge1 buffer
          ! ===================================

          kptr=0
          call edgeVunpack(edge1, C(1,1,1,ie), nlev, kptr, elem(ie)%desc)

          ! ===============================================
          ! Complete global assembly by normalizing by rmv
          ! ===============================================

          do k=1,nlev

             do j=1,nv
                do i=1,nv
                   D(i,j,k,ie) = elem(ie)%rmv(i,j)*C(i,j,k,ie)
                end do
             end do

          end do

          do k=1,nlev

             do j=1,nv
                do i=1,nv
                   C(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,nv
                   do i=1,nv
                      C(i,j,k,ie) = C(i,j,k,ie) + Emat_inv(l,k)*D(i,j,l,ie)
                   end do
                end do
             end do

          end do

       end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
       ! ==========================================
       ! solve for Gamma_ref, given C as RHS input
       ! ==========================================

       Gamma_ref = pcg_solver(elem, &
            C,          &
            cg,         &
            red,        &
            edge1,      &
            edge2,      &   
            Lambda,     &   
            deriv,      &   
            nets,       & 
            nete,       &
            blkjac)


       ! ================================================================
       ! Backsubstitute Gamma_ref into semi-implicit system of equations
       ! to find prognostic variables at time level n+1
       ! ================================================================

       do ie = nets, nete

          dx=0.5_real_kind*elem(ie)%dx/rearth ! weak derivative element x dimension
          dy=0.5_real_kind*elem(ie)%dy/rearth ! weak derivative element y dimension

          do k=1,nlev

             do j=1,nv
                do i=1,nv
                   Gref(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,nv
                   do i=1,nv
                      Gref(i,j,k,ie) = Gref(i,j,k,ie) + Emat(l,k)*Gamma_ref(i,j,l,ie)
                   end do
                end do
             end do

             do j=1,nv  
                do i=1,nv
                   B(i,j,k,ie) = elem(ie)%mv(i,j) * dt * (B(i,j,k,ie) - Gref(i,j,k,ie)) 
                end do
             end do

          end do

          kptr=0
          call edgeVpack(edge1,B(:,:,:,ie),nlev,kptr,elem(ie)%desc)

       end do

       call bndry_exchangeV(cg%hybrid,edge1)

       do ie = nets, nete

          dx=0.5_real_kind*elem(ie)%dx/rearth ! weak derivative element x dimension
          dy=0.5_real_kind*elem(ie)%dy/rearth ! weak derivative element y dimension

          kptr=0
          call edgeVunpack(edge1, B(:,:,:,ie), nlev, kptr, elem(ie)%desc)
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
          do k=1,nlev

             do j=1,nv
                do i=1,nv
                   B(i,j,k,ie) = elem(ie)%rmv(i,j)*B(i,j,k,ie)
                end do
             end do

          end do

          do k=1,nlev

             do j=1,nv
                do i=1,nv
                   D(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,nv
                   do i=1,nv
                      D(i,j,k,ie) = D(i,j,k,ie) + Emat_inv(l,k)*B(i,j,l,ie)
                   end do
                end do
             end do

          end do

#if 1
          do k=1,nlev 

             do j=1,nv
                do i=1,nv
                   elem(ie)%derived%div(i,j,k,np1) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,nv
                   do i=1,nv
                      elem(ie)%derived%div(i,j,k,np1) = elem(ie)%derived%div(i,j,k,np1) + Emat(l,k)*D(i,j,l,ie)/Lambda(l)
                   end do
                end do
             end do

          end do
#endif

          do k=1,nlev

             grad_Gref(:,:,:,k)=gradient_wk(Gref(:,:,k,ie),deriv,dx,dy)

             do j=1,nv
                do i=1,nv
                   gGr1 = grad_Gref(i,j,1,k)
                   gGr2 = grad_Gref(i,j,2,k)
                   elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%Dinv(1,1,i,j)*gGr1 + elem(ie)%Dinv(2,1,i,j)*gGr2
                   elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%Dinv(1,2,i,j)*gGr1 + elem(ie)%Dinv(2,2,i,j)*gGr2
                end do
             end do

             do j=1,nv
                do i=1,nv
                   Pscript(i,j,ie) = Pscript(i,j,ie) - dt*Pvec(k)*elem(ie)%derived%div(i,j,k,np1)
                end do
             end do


             do l=1,nlev
                do j=1,nv
                   do i=1,nv
                      Tscript(i,j,k,ie) = Tscript(i,j,k,ie) - dt*Tmat(l,k)*elem(ie)%derived%div(i,j,l,np1)
                   end do
                end do
             end do

          end do

          do j=1,nv
             do i=1,nv
                Pscript(i,j,ie) = elem(ie)%mv(i,j)*Pscript(i,j,ie)
             end do
          end do
          do k=1,nlev
             do j=1,nv
                do i=1,nv
                   Tscript(i,j,k,ie) = elem(ie)%mv(i,j)*Tscript(i,j,k,ie)
                end do
             end do
          end do

          ! ===============================================
          ! Pack v at time level n+1 into the edge3p1 buffer
          ! ===============================================

          kptr=0
          call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,elem(ie)%desc)

          kptr=2*nlev
          call edgeVpack(edge3p1, Tscript(:,:,:,ie),nlev,kptr,elem(ie)%desc)

          kptr=3*nlev
          call edgeVpack(edge3p1, Pscript(:,:,ie),1,kptr,elem(ie)%desc)

       end do

       ! ======================================
       ! boundary exchange v at time level n+1
       ! ======================================

       call bndry_exchangeV(cg%hybrid,edge3p1)

       do ie=nets,nete

          ! ===================================
          ! Unpack v from the edge2 buffer
          ! ===================================

          kptr=0
          call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, elem(ie)%desc)

          kptr=2*nlev
          call edgeVunpack(edge3p1, Tscript(:,:,:,ie), nlev, kptr, elem(ie)%desc)

          kptr=3*nlev
          call edgeVunpack(edge3p1, Pscript(:,:,ie), 1, kptr, elem(ie)%desc)

          ! ==========================================================
          ! Complete global assembly by normalizing velocity by rmv
          ! Vscript = Vscript - dt*grad(Gref)
          ! ==========================================================

          do k=1,nlev

             do j=1,nv
                do i=1,nv
                   elem(ie)%state%v(i,j,1,k,np1) = Vscript(i,j,1,k,ie) + dt*elem(ie)%rmv(i,j)*elem(ie)%state%v(i,j,1,k,np1)
                   elem(ie)%state%v(i,j,2,k,np1) = Vscript(i,j,2,k,ie) + dt*elem(ie)%rmv(i,j)*elem(ie)%state%v(i,j,2,k,np1)
                end do
             end do

          end do

          do j=1,nv
             do i=1,nv
                elem(ie)%state%lnps(i,j,np1) = elem(ie)%rmv(i,j)*Pscript(i,j,ie)
             end do
          end do

          do k=1,nlev
             do j=1,nv
                do i=1,nv
                   elem(ie)%state%T(i,j,k,np1) = elem(ie)%rmv(i,j)*Tscript(i,j,k,ie)
                end do
             end do
          end do

       end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
#if 1
       call prim_diffusion(elem, nets,nete,np1,deriv,dt2,cg%hybrid)
#endif

       call t_stopf('prim_advance_si')
       end subroutine prim_advance_si


  subroutine preq_robert3(nm1,n0,np1,elem,hvcoord,nets,nete)
  use dimensions_mod, only : nv, nlev, qsize
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
  use time_mod, only: smooth
  use control_mod, only : TRACERADV_TOTAL_DIVERGENCE, tracer_advection_formulation, integration, tstep_type

  implicit none
  integer              , intent(in) :: nm1,n0,np1,nets,nete
  type (hvcoord_t), intent(in)      :: hvcoord
  type (element_t)     , intent(inout) :: elem(:)
  
  
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: dp
  logical :: filter_ps = .false.
  if (integration == "explicit") filter_ps = .true.

  call t_startf('preq_robert')
  do ie=nets,nete
     if (filter_ps) then
        elem(ie)%state%ps_v(:,:,n0) = elem(ie)%state%ps_v(:,:,n0) + smooth*(elem(ie)%state%ps_v(:,:,nm1) &
             - 2.0D0*elem(ie)%state%ps_v(:,:,n0)   + elem(ie)%state%ps_v(:,:,np1))
        elem(ie)%state%lnps(:,:,n0) = LOG(elem(ie)%state%ps_v(:,:,n0))
     else
        elem(ie)%state%lnps(:,:,n0) = elem(ie)%state%lnps(:,:,n0) + smooth*(elem(ie)%state%lnps(:,:,nm1) &
             - 2.0D0*elem(ie)%state%lnps(:,:,n0)   + elem(ie)%state%lnps(:,:,np1))
        elem(ie)%state%ps_v(:,:,n0) = EXP(elem(ie)%state%lnps(:,:,n0))
     endif
     
     elem(ie)%state%T(:,:,:,n0) = elem(ie)%state%T(:,:,:,n0) + smooth*(elem(ie)%state%T(:,:,:,nm1) &
          - 2.0D0*elem(ie)%state%T(:,:,:,n0)   + elem(ie)%state%T(:,:,:,np1))
     elem(ie)%state%v(:,:,:,:,n0) = elem(ie)%state%v(:,:,:,:,n0) + smooth*(elem(ie)%state%v(:,:,:,:,nm1) &
          - 2.0D0*elem(ie)%state%v(:,:,:,:,n0) + elem(ie)%state%v(:,:,:,:,np1))
     
     ! only apply Robert filter if running leapfrog for tracers
     if (tstep_type==0) then
        do q=1,qsize
           if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
              elem(ie)%state%Qdp(:,:,:,q,n0) = elem(ie)%state%Qdp(:,:,:,q,n0) + &
                   smooth*(elem(ie)%state%Qdp(:,:,:,q,nm1) &
                   - 2.0D0*elem(ie)%state%Qdp(:,:,:,q,n0)   + elem(ie)%state%Qdp(:,:,:,q,np1))
              ! Qdp(n0) and ps_v(n0) were filtered -update Q(n)
              do k=1,nlev
                 do j=1,nv
                    do i=1,nv
                       dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)
                       elem(ie)%state%Q(i,j,k,q,n0) = elem(ie)%state%Qdp(i,j,k,q,n0)/dp
                    enddo
                 enddo
              enddo
           else
              elem(ie)%state%Q(:,:,:,q,n0) = elem(ie)%state%Q(:,:,:,q,n0) + &
                   smooth*(elem(ie)%state%Q(:,:,:,q,nm1) &
                   - 2.0D0*elem(ie)%state%Q(:,:,:,q,n0)   + elem(ie)%state%Q(:,:,:,q,np1))
           endif
        enddo
     endif
  end do
  call t_stopf('preq_robert')
  
  end subroutine preq_robert3
  
  


  subroutine applyCAMforcing_leapfrog(elem,hvcoord,n0,np1,dt,nets,nete)
  use dimensions_mod, only : nv, nlev, qsize
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t
  use control_mod, only : TRACERADV_TOTAL_DIVERGENCE, tracer_advection_formulation,&
       tstep_type, moisture
  use time_mod, only : smooth
  
  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  real (kind=real_kind), intent(in) :: dt
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: n0,np1,nets,nete
  
  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: dt2,dt_q,v1,dp
  logical :: wet  
  
  wet = (moisture /= "dry")  
  
  dt2=2*dt
  dt_q=dt2
  if (tstep_type==1) dt_q=dt
  
  do ie=nets,nete
     ! apply forcing to Q
     if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
        
        elem(ie)%derived%FQps(:,:,np1)=0
        do q=1,qsize
           do k=1,nlev
              do j=1,nv	
                 do i=1,nv
                    v1 = dt_q*elem(ie)%derived%FQ(i,j,k,q,1)
                    if (elem(ie)%state%Qdp(i,j,k,q,np1) + v1 < 0 .and. v1<0) then
                       if (elem(ie)%state%Qdp(i,j,k,q,np1) < 0 ) then
                          v1=0  ! Q already negative, dont make it more so
                       else
                          v1 = -elem(ie)%state%Qdp(i,j,k,q,np1)
                       endif
                    endif
                    elem(ie)%state%Qdp(i,j,k,q,np1) = elem(ie)%state%Qdp(i,j,k,q,np1)+v1
                    if (q==1) then
                       elem(ie)%derived%FQps(i,j,np1)=elem(ie)%derived%FQps(i,j,np1)+v1/dt_q
                    endif
                 enddo
              enddo
           enddo
        enddo

        if (wet .and. qsize>0) then
           elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,np1) + dt_q*elem(ie)%derived%FQps(:,:,np1)
        endif
        ! Qdp(np1) and ps_v(np1) were updated by forcing - update Q(np1)
        ! ps_v(n0) may also have been changed if using Robert,
        ! but Q(n0) will be updated after robert filter
        ! so no need to do that now
        do q=1,qsize
           do k=1,nlev
              do j=1,nv	
                 do i=1,nv
                    dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                         ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,np1)
                    elem(ie)%state%Q(i,j,k,q,np1) = elem(ie)%state%Qdp(i,j,k,q,np1)/dp
                 enddo
              enddo
           enddo
        enddo
        
     else
        elem(ie)%state%Q(:,:,:,:,np1) = elem(ie)%state%Q(:,:,:,:,np1)+&
             dt_q*elem(ie)%derived%FQ(:,:,:,:,1)
     endif
     
     elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
          dt2*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
          dt2*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing_leapfrog



  subroutine applyCAMforcing(elem,hvcoord,np1,dt_q,nets,nete)
  use dimensions_mod, only : nv, nlev, qsize
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t
  use control_mod, only : TRACERADV_TOTAL_DIVERGENCE, tracer_advection_formulation,&
       moisture
  use time_mod, only : smooth
  
  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  real (kind=real_kind), intent(in) :: dt_q
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: np1,nets,nete
  
  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  logical :: wet  
  
  wet = (moisture /= "dry")  
  
  do ie=nets,nete
     ! apply forcing to Qdp
     elem(ie)%derived%FQps(:,:,1)=0
     do q=1,qsize
        do k=1,nlev
           do j=1,nv	
              do i=1,nv
                 v1 = dt_q*elem(ie)%derived%FQ(i,j,k,q,1)
                 if (elem(ie)%state%Qdp(i,j,k,q,np1) + v1 < 0 .and. v1<0) then
                    if (elem(ie)%state%Qdp(i,j,k,q,np1) < 0 ) then
                       v1=0  ! Q already negative, dont make it more so
                    else
                       v1 = -elem(ie)%state%Qdp(i,j,k,q,np1)
                    endif
                 endif
                 elem(ie)%state%Qdp(i,j,k,q,np1) = elem(ie)%state%Qdp(i,j,k,q,np1)+v1
                 if (q==1) then
                    elem(ie)%derived%FQps(i,j,1)=elem(ie)%derived%FQps(i,j,1)+v1/dt_q
                 endif
              enddo
           enddo
        enddo
     enddo
     
     if (wet .and. qsize>0) then
        ! to conserve dry mass in the precese of Q1 forcing:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,np1) + &
             dt_q*elem(ie)%derived%FQps(:,:,1)
     endif
     
     ! Qdp(np1) and ps_v(np1) were updated by forcing - update Q(np1)
     ! ps_v(n0) may also have been changed if using Robert,
     ! but Q(n0) will be updated after robert filter
     ! so no need to do that now
     do q=1,qsize
        do k=1,nlev
           do j=1,nv	
              do i=1,nv
                 dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,np1)
                 elem(ie)%state%Q(i,j,k,q,np1) = elem(ie)%state%Qdp(i,j,k,q,np1)/dp
              enddo
           enddo
        enddo
     enddo
     
     elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
          dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
          dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing



  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,dt_q,nets,nete)
  use dimensions_mod, only : nv, nlev, qsize
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t
  use time_mod, only : smooth
  
  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  real (kind=real_kind), intent(in) :: dt_q
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: np1,nets,nete
  
  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  logical :: wet  
  
  do ie=nets,nete
     elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
          dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
          dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing_dynamics



  subroutine advance_hypervis(edge3,elem,hvcoord,hybrid,deriv,nm1,n0,nt,nets,nete,dt2)
  !
  !  take one timestep of:  
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : nv, np, nlev
  use control_mod, only : nu, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, moisture, psurf_vis
  use hybrid_mod, only : hybrid_t
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk
  use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
  use physical_constants, only: Cp,Cpwater_vapor
!  use time_mod, only : TimeLevel_t
  implicit none
  
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
!  type (TimeLevel_t)   , intent(in) :: tl
  
  real (kind=real_kind) :: dt2
  integer :: nets,nete
  
  ! local
  real (kind=real_kind) :: nu_scale, dpdn,dpdn0, cp_star, nu_scale_top
  logical :: use_cp_star=.false.
  integer :: k,kptr,i,j,ie,ic,n0,nt,nm1
  real (kind=real_kind), dimension(nv,nv,2,nlev,nets:nete)      :: vtens   
  real (kind=real_kind), dimension(nv,nv,nlev,nets:nete)        :: ptens
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens	
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dXdp
  
  
! NOTE: PGI compiler bug: when using spheremv, rspheremv and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.  
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremv,rspheremv
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps
  
  real (kind=real_kind), dimension(nv,nv) :: lap_p
  real (kind=real_kind), dimension(nv,nv,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp

  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
  call t_startf('advance_hypervis')
  
!  nm1 = tl%nm1   ! heating term uses U,V at average of nt and nm1 levels
!  n0 = tl%n0     ! timelevel used for ps scaling.  use n0 for leapfrog.  
!  nt = tl%np1    ! apply viscosity to this timelevel  (np1)
  
  cp_star = cp
  !       "improvements" from this term do not seem to help TE balance
  !       if(moisture /= "dry")then
  !          use_cp_star=.true.
  !       endif
  
  
  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) stop 'ERROR: hypervis_order == 1 not coded for nu_p>0'
     do ic=1,hypervis_subcycle
        do ie=nets,nete
           
           do k=1,nlev
              lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie))
              lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
              ! advace in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,nv
                 do i=1,nv             
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremv(i,j)  +  dt*nu_s*lap_p(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremv(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremv(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo
           
           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,elem(ie)%desc)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,elem(ie)%desc)
        enddo
        
        call bndry_exchangeV(hybrid,edge3)
        
        do ie=nets,nete
           
           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, elem(ie)%desc)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, elem(ie)%desc)
           
           ! apply inverse mass matrix
           do k=1,nlev
              do j=1,nv
                 do i=1,nv             
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremv(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremv(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremv(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
     enddo  ! subcycle
  endif
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
        do ie=nets,nete
           nu_scale=1
           do k=1,nlev
              ! advace in time.  
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"
              
              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie))
                 lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2
              
              do j=1,nv
                 do i=1,nv
                    if (psurf_vis==0) then
                       ! normalize so as to conserve IE  (not needed when using p-surface viscosity)
                       ! scale velosity by 1/rho (normalized to be O(1))
                       ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                       dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)  ! nt ?
                       dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                       nu_scale = dpdn0/dpdn
                    endif

                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=nu_scale*(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=nu_scale*(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
                    else
                       utens_tmp=-nu_scale*nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu_scale*nu*vtens(i,j,2,k,ie)
                       ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
                    endif
                    
                    ptens(i,j,k,ie) = ptens_tmp  
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo
           
           pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
           kptr=0
           call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,elem(ie)%desc)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,elem(ie)%desc)
           kptr=3*nlev
           call edgeVpack(edge3,pstens(:,:,ie),1,kptr,elem(ie)%desc)
        enddo
        
        
        call bndry_exchangeV(hybrid,edge3)
        
        do ie=nets,nete
           
           kptr=0
           call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, elem(ie)%desc)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, elem(ie)%desc)
           kptr=3*nlev
           call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, elem(ie)%desc)

           if (psurf_vis == 1 ) then
              ! apply p-surface correction
              do k=1,nlev
                 p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(:,:,nt)
              enddo
              do k=1,nlev
                 if (k.eq.1) then
                    ! no correction needed
                 else if (k.eq.nlev) then
                    ! one-sided difference
                    dXdp = (elem(ie)%state%T(:,:,k,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
                        (p(:,:,k)-p(:,:,k-1)) 
                    ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
                 else
                    dXdp = (elem(ie)%state%T(:,:,k+1,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
                         (p(:,:,k+1)-p(:,:,k-1)) 
                    ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
                 endif
              enddo
           endif

           
           ! apply inverse mass matrix, accumulate tendencies
           do k=1,nlev
              do j=1,nv
                 do i=1,nv
                    ! note: Q has not been updated yet - so this code may be broken.  need
                    ! to use Q at n0
                    ! if(use_cp_star) then
                    !    cp_star = cp + (cpwater_vapor-cp)*elem(ie)%state%Q(i,j,k,1,n0)
                    !    cp_star = cp + (cpwater_vapor-cp)*&
                    !       ( elem(ie)%state%Q(i,j,k,1,nt) + elem(ie)%state%Q(i,j,k,1,nm1))/2
                    ! endif
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt) + &
                         dt*elem(ie)%rspheremv(i,j)*vtens(i,j,1,k,ie)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt) +  &
                         dt*elem(ie)%rspheremv(i,j)*vtens(i,j,2,k,ie) 
                    
                    v1=.5*(elem(ie)%state%v(i,j,1,k,nt)+elem(ie)%state%v(i,j,1,k,nm1))
                    v2=.5*(elem(ie)%state%v(i,j,2,k,nt)+elem(ie)%state%v(i,j,2,k,nm1))
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)     + &
                         dt*elem(ie)%rspheremv(i,j)*(cp*ptens(i,j,k,ie) - heating)/cp_star
                    
                 enddo
              enddo
           enddo
           elem(ie)%state%ps_v(:,:,nt)=elem(ie)%state%ps_v(:,:,nt) + dt*elem(ie)%rspheremv(:,:)*pstens(:,:,ie)

           
           
        enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
     enddo
  endif
  
  call t_stopf('advance_hypervis')
  
  end subroutine advance_hypervis
  
  


  subroutine compute_and_apply_rhs(np1,nm1,n0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics,n_Q,eta_ave_w)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine is normally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
  !
  ! Combining the RHS and DSS pack operation in one routine 
  ! allows us to fuse these two loops for more cache reuse
  !
  ! Combining the dt advance and DSS unpack operation in one routine 
  ! allows us to fuse these two loops for more cache reuse
  !
  ! note: for prescribed velocity case, velocity will be computed at
  ! "real_time", which should be the time of timelevel n0.  
  !
  ! note: Q(timelevel=n_Q) is used to compute Tv
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : nv, np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use control_mod, only : moisture, qsplit, compute_mean_flux
  use hybvcoord_mod, only : hvcoord_t

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv, preq_omega_ps, preq_hydrostatic


  implicit none
  integer :: np1,nm1,n0,nets,nete,n_Q
  real*8 :: dt2
  logical  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
  real (kind=real_kind), pointer, dimension(:,:)      :: ps         ! surface pressure for current tiime level
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi

  real (kind=real_kind), dimension(nv,nv,nlev)   :: omega_p       
  real (kind=real_kind), dimension(nv,nv,nlev)   :: T_v         
  real (kind=real_kind), dimension(nv,nv,nlev)   :: divdp
  real (kind=real_kind), dimension(nv,nv,nlev+1)   :: eta_dot_dpdn  ! half level vertical velocity on p-grid
  real (kind=real_kind), dimension(nv,nv)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(nv,nv,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind), dimension(nv,nv)      :: vgrad_T    ! v.grad(T)
  real (kind=real_kind), dimension(nv,nv)      :: Ephi       ! kinetic energy + PHI term
  real (kind=real_kind), dimension(nv,nv,2)      :: grad_ps    ! lat-lon coord version
  real (kind=real_kind), dimension(nv,nv,nlev)   :: vort       ! vorticity
  real (kind=real_kind), dimension(nv,nv,nlev)   :: p          ! pressure
  real (kind=real_kind), dimension(nv,nv,nlev)   :: dp         ! delta pressure
  real (kind=real_kind), dimension(nv,nv,nlev)   :: rdp        ! inverse of delta pressure
  real (kind=real_kind), dimension(nv,nv,nlev)   :: T_vadv     ! temperature vertical advection
  real (kind=real_kind), dimension(nv,nv,nlev)   :: vgrad_ps ! v.grad(lnps) 
  real (kind=real_kind), dimension(nv,nv,nlev+1) :: ph               ! half level pressures on p-grid
  real (kind=real_kind), dimension(nv,nv,2,nlev) :: v_vadv   ! velocity vertical advection
  real (kind=real_kind) ::  Kappa_star(nv,nv,nlev)
  real (kind=real_kind) ::  vtens1(nv,nv,nlev)
  real (kind=real_kind) ::  vtens2(nv,nv,nlev)
  real (kind=real_kind) ::  ttens(nv,nv,nlev)

  real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2
  real (kind=real_kind) ::  glnps1,glnps2,gpterm
  integer :: i,j,k,kptr,ie
  logical :: wet

  do ie=nets,nete
     ps => elem(ie)%state%ps_v(:,:,n0)
!     divdp => elem(ie)%derived%div(:,:,:,n0)
     phi => elem(ie)%derived%phi(:,:,:)
!     omega_p => elem(ie)%derived%omega_p(:,:,:)
     
     ! ==================================================
     ! compute pressure (p) on half levels from ps 
     ! using the hybrid coordinates relationship, i.e.
     ! e.g. equation (3.a.92) of the CCM-2 description, 
     ! (NCAR/TN-382+STR), June 1993, p. 24.
     ! ==================================================
     do k=1,nlev+1
        do j=1,nv
           do i=1,nv
              ph(i,j,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*ps(i,j)
           end do
        end do
     end do
     
     ! ============================
     ! compute p and delta p
     ! ============================
     do k=1,nlev
        do j=1,nv
           do i=1,nv
              p(i,j,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
              dp(i,j,k)  = ph(i,j,k+1) - ph(i,j,k)
              rdp(i,j,k) = 1.0D0/dp(i,j,k)
           end do
        end do
        ! ============================
        ! compute vgrad_lnps 
        ! ============================
        grad_ps = gradient_sphere(elem(ie)%state%ps_v(:,:,n0),deriv,elem(ie))
        do j=1,nv
           do i=1,nv
              v1 = elem(ie)%state%v(i,j,1,k,n0)
              v2 = elem(ie)%state%v(i,j,2,k,n0)
              vgrad_ps(i,j,k) = v1*grad_ps(i,j,1) + v2*grad_ps(i,j,2)
              vtemp(i,j,1) = v1*dp(i,j,k)
              vtemp(i,j,2) = v2*dp(i,j,k)
           end do
        end do


      
        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
        if(compute_mean_flux==1)then
        do j=1,nv
           do i=1,nv
              elem(ie)%derived%vn0(i,j,:,k)=elem(ie)%derived%vn0(i,j,:,k)+eta_ave_w*vtemp(i,j,:)
           end do
        end do        
        endif


        ! =========================================
        !
        ! Compute relative vorticity and divergence
        !
        ! =========================================
        divdp(:,:,k)=divergence_sphere(vtemp,deriv,elem(ie))
        vort(:,:,k)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie))

     enddo
     
     
     
     ! compute T_v for timelevel n0
     wet = ( moisture /= "dry")
     if(wet) then
        do k=1,nlev
           do j=1,nv
              do i=1,nv
                 Qt = elem(ie)%state%Q(i,j,k,1,n_Q) 
                 T_v(i,j,k) = Virtual_Temperature(elem(ie)%state%T(i,j,k,n0),Qt)
                 kappa_star(i,j,k) =  Rgas/Virtual_Specific_Heat(Qt)
              end do
           end do
        end do
     else
        do k=1,nlev
           do j=1,nv
              do i=1,nv
                 T_v(i,j,k) = elem(ie)%state%T(i,j,k,n0)
                 kappa_star(i,j,k) = kappa
              end do
           end do
        end do
     end if
     
     
     
     ! ====================================================
     ! Compute Hydrostatic equation, modeld after CCM-3
     ! ====================================================
     !call geopotential_t(p,dp,T_v,Rgas,phi)
     call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p,dp)
     
     ! ====================================================
     ! Compute omega_p according to CCM-3 
     ! ====================================================
     call preq_omega_ps(omega_p,hvcoord,ps,p,vgrad_ps,divdp)
     
     ! ==================================================
     ! zero partial sum for accumulating sum
     !    (div(v_k) + v_k.grad(lnps))*dsigma_k = div( v dp )
     ! used by eta_dot_dpdn and lnps tendency
     ! ==================================================
     sdot_sum=0
     
     ! ==================================================
     ! Compute eta_dot_dpdn 
     ! save sdot_sum as this is the -RHS of ps_v equation
     ! ==================================================
     do k=1,nlev
        ! ==================================================
        ! add this term to PS equation so we exactly conserve dry mass
        ! ==================================================
        do j=1,nv
           do i=1,nv
              sdot_sum(i,j) = sdot_sum(i,j) + divdp(i,j,k) 
              eta_dot_dpdn(i,j,k+1) = sdot_sum(i,j)      
           end do
        end do
     end do
     
     
     ! ===========================================================
     ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
     ! compute at interfaces:
     !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
     ! for reference: at mid layers we have:
     !    omega = v grad p  - integral_etatop^eta[ divdp ]
     ! ===========================================================
     do k=1,nlev-1
        eta_dot_dpdn(:,:,k+1) = hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
     end do
     
     eta_dot_dpdn(:,:,1     ) = 0.0D0
     eta_dot_dpdn(:,:,nlev+1) = 0.0D0



     ! accumulate mean fluxes:
     elem(ie)%derived%eta_dot_dpdn(:,:,:) = &
          elem(ie)%derived%eta_dot_dpdn(:,:,:) + eta_ave_w*eta_dot_dpdn(:,:,:)
     elem(ie)%derived%omega_p(:,:,:) = &
          elem(ie)%derived%omega_p(:,:,:) + eta_ave_w*omega_p(:,:,:)


     ! ===========================================================
     ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
     ! ==============================================
     call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
          eta_dot_dpdn,rdp,T_vadv,v_vadv)
     
     
     
     ! ==============================================
     ! Compute phi + kinetic energy term: 10*nv*nv Flops
     ! ==============================================
     do k=1,nlev   
        do j=1,nv
           do i=1,nv
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              E = 0.5D0*( v1*v1 + v2*v2 )
              Ephi(i,j)=E+phi(i,j,k)
           end do
        end do
        ! ================================================
        ! compute gradp term (ps/p)*(dp/dps)*T
        ! ================================================
        vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie))
        do j=1,nv
           do i=1,nv
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2) 
           end do
        end do
        
        
        ! vtemp = grad ( E + PHI )
        vtemp = gradient_sphere(Ephi(:,:),deriv,elem(ie))
        
        do j=1,nv
           do i=1,nv
              gpterm = hvcoord%hybm(k)*T_v(i,j,k)/p(i,j,k)
              glnps1 = Rgas*gpterm*grad_ps(i,j,1)
              glnps2 = Rgas*gpterm*grad_ps(i,j,2)
              
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              
              vtens1(i,j,k) =   - v_vadv(i,j,1,k)                           &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,1) - glnps1   
              
              vtens2(i,j,k) =   - v_vadv(i,j,2,k)                            &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,2) - glnps2   
              
              ttens(i,j,k)  = - T_vadv(i,j,k) - vgrad_T(i,j) + Kappa_star(i,j,k)*T_v(i,j,k)*omega_p(i,j,k)

           end do
        end do



     end do
       
#ifdef ENERGY_DIAGNOSTICS
     ! =========================================================
     !
     ! diagnostics
     ! recomputes some gradients that were not saved above
     ! uses:  sdot_sum(), eta_dot_dpdn(), grad_ps()
     ! grad_phi(), dp(), p(), T_vadv(), v_vadv(), divdp()
     ! =========================================================
     if (compute_diagnostics) then
        elem(ie)%accum%KEhorz1=0
        elem(ie)%accum%KEhorz2=0
        elem(ie)%accum%IEhorz1=0
        elem(ie)%accum%IEhorz2=0
        elem(ie)%accum%IEhorz1_wet=0
        elem(ie)%accum%IEhorz2_wet=0
        elem(ie)%accum%KEvert1=0
        elem(ie)%accum%KEvert2=0
        elem(ie)%accum%IEvert1=0
        elem(ie)%accum%IEvert2=0
        elem(ie)%accum%IEvert1_wet=0
        elem(ie)%accum%IEvert2_wet=0
        elem(ie)%accum%T1=0
        elem(ie)%accum%T2=0
        elem(ie)%accum%T2_s=0
        elem(ie)%accum%S1=0
        elem(ie)%accum%S1_wet=0
        elem(ie)%accum%S2=0
        
        do j=1,nv
           do i=1,nv
              elem(ie)%accum%S2(i,j) = elem(ie)%accum%S2(i,j) - &
                   sdot_sum(i,j)*elem(ie)%state%phis(i,j)
           enddo
        enddo
        
        do k=1,nlev
           ! vtemp = grad_E(:,:,k)
           do j=1,nv
              do i=1,nv
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 Ephi(i,j)=0.5D0*( v1*v1 + v2*v2 )
              enddo
           enddo
           vtemp = gradient_sphere(Ephi,deriv,elem(ie))
           do j=1,nv
              do i=1,nv
                 ! dp/dn u dot grad(E)
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 elem(ie)%accum%KEhorz2(i,j) = elem(ie)%accum%KEhorz2(i,j) + &
                      (v1*vtemp(i,j,1)  + v2*vtemp(i,j,2))*dp(i,j,k)
                 ! E div( u dp/dn )
                 elem(ie)%accum%KEhorz1(i,j) = elem(ie)%accum%KEhorz1(i,j) + Ephi(i,j)*divdp(i,j,k)
                 
                 ! Cp T div( u dp/dn)   ! dry horizontal advection component
                 elem(ie)%accum%IEhorz1(i,j) = elem(ie)%accum%IEhorz1(i,j) + Cp*elem(ie)%state%T(i,j,k,n0)*divdp(i,j,k)
                 
                 
              enddo
           enddo
           
           
           ! vtemp = grad_phi(:,:,k)
           vtemp = gradient_sphere(phi(:,:,k),deriv,elem(ie))
           do j=1,nv
              do i=1,nv
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 E = 0.5D0*( v1*v1 + v2*v2 )
                 ! NOTE:  Cp_star = Cp + (Cpv-Cp)*q   
                 ! advection terms can thus be broken into two components: dry and wet
                 ! dry components cancel exactly
                 ! wet components should cancel exactly
                 ! 
                 ! some diagnostics
                 ! e = eta_dot_dpdn()    
                 de =  eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k) 
                 ! Cp T de/dn, integral dn:
                 elem(ie)%accum%IEvert1(i,j)=elem(ie)%accum%IEvert1(i,j) + Cp*elem(ie)%state%T(i,j,k,n0)*de
                 ! E de/dn
                 elem(ie)%accum%KEvert1(i,j)=elem(ie)%accum%KEvert1(i,j) + E*de
                 ! Cp T_vadv dp/dn
                 elem(ie)%accum%IEvert2(i,j)=elem(ie)%accum%IEvert2(i,j) + Cp*T_vadv(i,j,k)*dp(i,j,k)
                 ! dp/dn V dot V_vadv
                 elem(ie)%accum%KEvert2(i,j)=elem(ie)%accum%KEvert2(i,j) + (v1*v_vadv(i,j,1,k) + v2*v_vadv(i,j,2,k)) *dp(i,j,k)
                 
                 ! IEvert1_wet():  (Cpv-Cp) T Qdp_vadv  (Q equation)
                 ! IEvert2_wet():  (Cpv-Cp) Qdp T_vadv   T equation
                 if (wet) then
                 elem(ie)%accum%IEvert2_wet(i,j)=elem(ie)%accum%IEvert2_wet(i,j) +&
                      (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1,n_Q)*T_vadv(i,j,k)*dp(i,j,k)
                 endif

                 gpterm = hvcoord%hybm(k)*T_v(i,j,k)/p(i,j,k)
                 elem(ie)%accum%T1(i,j) = elem(ie)%accum%T1(i,j) - &
                      Rgas*gpterm*(grad_ps(i,j,1)*v1 + grad_ps(i,j,2)*v2)*dp(i,j,k)
                 
                 elem(ie)%accum%T2(i,j) = elem(ie)%accum%T2(i,j) - &
                      (vtemp(i,j,1)*v1 + vtemp(i,j,2)*v2)*dp(i,j,k)
                 
                 ! S1 = < Cp_star dp/dn , RT omega_p/cp_star >  
                 elem(ie)%accum%S1(i,j) = elem(ie)%accum%S1(i,j) + &
                      Rgas*T_v(i,j,k)*omega_p(i,j,k)*dp(i,j,k)
                 
                 ! cp_star = cp + cp2
                 if (wet) then
                 cp2 = (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1,n_Q)
                 cp_ratio = cp2/(cp+cp2)
                 elem(ie)%accum%S1_wet(i,j) = elem(ie)%accum%S1_wet(i,j) + &
                      cp_ratio*(Rgas*T_v(i,j,k)*omega_p(i,j,k)*dp(i,j,k))
                 endif
                 
                 elem(ie)%accum%CONV(i,j,:,k)=-Rgas*gpterm*grad_ps(i,j,:)-vtemp(i,j,:)
              enddo
           enddo
           
           vtemp(:,:,:) = gradient_sphere(elem(ie)%state%phis(:,:),deriv,elem(ie))
           do j=1,nv
              do i=1,nv
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 elem(ie)%accum%T2_s(i,j) = elem(ie)%accum%T2_s(i,j) - &
                      (vtemp(i,j,1)*v1 + vtemp(i,j,2)*v2)*dp(i,j,k)
              enddo
           enddo
           
           vtemp(:,:,:)   = gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie))
           do j=1,nv
              do i=1,nv
                 v1     = elem(ie)%state%v(i,j,1,k,n0)
                 v2     = elem(ie)%state%v(i,j,2,k,n0)
                 
                 ! Cp dp/dn u dot gradT
                 elem(ie)%accum%IEhorz2(i,j) = elem(ie)%accum%IEhorz2(i,j) + &
                      Cp*(v1*vtemp(i,j,1) + v2*vtemp(i,j,2))*dp(i,j,k)
                 
                 if (wet) then
                 elem(ie)%accum%IEhorz2_wet(i,j) = elem(ie)%accum%IEhorz2_wet(i,j) + &
                      (Cpwater_vapor-Cp)*elem(ie)%state%Q(i,j,k,1,n_Q)*&
                      (v1*vtemp(i,j,1) + v2*vtemp(i,j,2))*dp(i,j,k)
                 endif
                 
              enddo
           enddo
           
        enddo
     endif
#endif
       
     ! =========================================================
     ! local element timestep, store in np1.  
     ! note that we allow np1=n0 or nm1
     ! apply mass matrix
     ! =========================================================
     if (dt2<0) then
        ! calling program just wanted DSS'd RHS, skip time advance
        do k=1,nlev
           elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremv(:,:)*vtens1(:,:,k) 
           elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremv(:,:)*vtens2(:,:,k) 
           elem(ie)%state%T(:,:,k,np1) = elem(ie)%spheremv(:,:)*ttens(:,:,k)
        enddo
        elem(ie)%state%ps_v(:,:,np1) = -elem(ie)%spheremv(:,:)*sdot_sum
     else
        do k=1,nlev
           elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremv(:,:)*( elem(ie)%state%v(:,:,1,k,nm1) + dt2*vtens1(:,:,k) )
           elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremv(:,:)*( elem(ie)%state%v(:,:,2,k,nm1) + dt2*vtens2(:,:,k) )
           elem(ie)%state%T(:,:,k,np1) = elem(ie)%spheremv(:,:)*(elem(ie)%state%T(:,:,k,nm1) + dt2*ttens(:,:,k))
        enddo
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%spheremv(:,:)*( elem(ie)%state%ps_v(:,:,nm1) - dt2*sdot_sum )
     endif

     
     ! =========================================================
     !
     ! Pack ps(np1), T, and v tendencies into comm buffer
     !
     ! =========================================================
     kptr=0
     call edgeVpack(edge3p1, elem(ie)%state%ps_v(:,:,np1),1,kptr,elem(ie)%desc)
     
     kptr=1
     call edgeVpack(edge3p1, elem(ie)%state%T(:,:,:,np1),nlev,kptr,elem(ie)%desc)
     
     kptr=nlev+1
     call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,elem(ie)%desc)
     
  end do
  
  ! =============================================================
    ! Insert communications here: for shared memory, just a single
  ! sync is required
  ! =============================================================
  
  call bndry_exchangeV(hybrid,edge3p1)
  
  do ie=nets,nete
     ! ===========================================================
     ! Unpack the edges for vgrad_T and v tendencies...
     ! ===========================================================
     kptr=0
     call edgeVunpack(edge3p1, elem(ie)%state%ps_v(:,:,np1), 1, kptr, elem(ie)%desc)
     
     kptr=1
     call edgeVunpack(edge3p1, elem(ie)%state%T(:,:,:,np1), nlev, kptr, elem(ie)%desc)
     
     kptr=nlev+1
     call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, elem(ie)%desc)
     
     
     ! ====================================================
     ! Scale tendencies by inverse mass matrix, Leapfog update 
     ! ====================================================
     do j=1,nv
        do i=1,nv
           elem(ie)%state%ps_v(i,j,np1) = elem(ie)%rspheremv(i,j)*elem(ie)%state%ps_v(i,j,np1)
        end do
     end do
     do k=1,nlev
        do j=1,nv
           do i=1,nv
              elem(ie)%state%T(i,j,k,np1)   = elem(ie)%rspheremv(i,j)*elem(ie)%state%T(i,j,k,np1)
              elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%rspheremv(i,j)*elem(ie)%state%v(i,j,1,k,np1)
              elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%rspheremv(i,j)*elem(ie)%state%v(i,j,2,k,np1)
           end do
        end do
     end do
  end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
  end subroutine compute_and_apply_rhs
  



  subroutine smooth_phis(phis,elem,hybrid,deriv,nets,nete,minf,numcycle)
  use dimensions_mod, only : nv, np, nlev
  use control_mod, only : smooth_phis_nudt
  use hybrid_mod, only : hybrid_t
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t , laplace_sphere_wk
  use viscosity_mod, only : biharmonic_wk
  use time_mod, only : TimeLevel_t
  implicit none
  
  real (kind=real_kind), dimension(nv,nv,nets:nete), intent(inout)   :: phis
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind), intent(in)   :: minf
  integer,               intent(in) :: numcycle
  
  integer :: nets,nete
  ! local 
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens	
  real (kind=real_kind) :: mx,mn
  integer :: nt,ie,ic,i,j,order,order_max

  do ic=1,numcycle

     pstens=phis

     ! order = 1   laplacian
     ! order = 2   grad^4 (need to add a negative sign)
     ! order = 3   grad^6
     ! order = 4   grad^8 (need to add a negative sign)  
     order_max = 2
     do order=1,order_max-1
        do ie=nets,nete
           pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie))
           call edgeVpack(edge3p1,pstens(:,:,ie),1,0,elem(ie)%desc)
        enddo
        call bndry_exchangeV(hybrid,edge3p1)
        do ie=nets,nete
           call edgeVunpack(edge3p1, pstens(:,:,ie), 1, 0, elem(ie)%desc)
           pstens(:,:,ie)=pstens(:,:,ie)*elem(ie)%rspheremv(:,:)
        enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
     enddo
     do ie=nets,nete
        pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie))
     enddo
     if (mod(order_max,2)==0) pstens=-pstens

     do ie=nets,nete
        !  ps(t+1) = ps(t) + Minv * DSS * M * RHS
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  RHS ]
        ! but output of biharminc_wk is of the form M*RHS.  rewrite as:
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  M*RHS/M ]
        ! so we can apply limiter to ps(t) +  (M*RHS)/M

        mn=minval(phis(:,:,ie))
        mx=maxval(phis(:,:,ie))
        phis(:,:,ie)=phis(:,:,ie) + &
           smooth_phis_nudt*pstens(:,:,ie)/elem(ie)%spheremv(:,:)


        ! remove new extrema.  could use conservative reconstruction from advection
        ! but no reason to conserve mean PHI.
        if (ic < numcycle/2) then
        do i=1,nv
        do j=1,nv
           if (phis(i,j,ie)>mx) phis(i,j,ie)=mx
           if (phis(i,j,ie)<mn) phis(i,j,ie)=mn
        enddo
        enddo
        endif

        ! user specified minimum 
        do i=1,nv
        do j=1,nv
           if (phis(i,j,ie)<minf) phis(i,j,ie)=minf
        enddo
        enddo

        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%spheremv(:,:)
        call edgeVpack(edge3p1,phis(:,:,ie),1,0,elem(ie)%desc)
     enddo
     call bndry_exchangeV(hybrid,edge3p1)
     do ie=nets,nete
        call edgeVunpack(edge3p1, phis(:,:,ie), 1, 0, elem(ie)%desc)
        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%rspheremv(:,:)
     enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
  enddo
  end subroutine smooth_phis








end module prim_advance_mod
