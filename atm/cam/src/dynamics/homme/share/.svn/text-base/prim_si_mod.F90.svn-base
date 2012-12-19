#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_si_mod
  implicit none
  private
  public :: preq_impsys
  public :: preq_omegap
  public :: preq_omega_ps
  public :: preq_omega_lnps
  public :: preq_hydrostatic, geopotential_t
  public :: preq_pressure
  public :: preq_vertadv
  public :: preq_backsub
  public :: preq_bldrhs
contains
	
! ==========================================================
! Implicit system for semi-implicit primitive equations.
! ==========================================================

  subroutine preq_impsys(dt,dx,dy,rdx,rdy,deriv,                &
       np1,n0,nm1,                            &
       Pscript,Tscript,Vscript,               &
       met,metdet,rmetdet,fcor,               &
       phis,lnps,T,v,                         &
       div, zeta, grad_lnps,                  &
       hyai, hybi, hyam, hybm, hybd, nprlev,  &
       ps0, psref,Tref,RTref,Href,Tmat,Pvec,  &
       Emat,Emat_inv,Amat_inv,Lambda,         &
       B, grad_dGref)
    use kinds,              only : real_kind
    use derivative_mod,     only : gradient, vorticity, gradient_wk, derivative_t
    use dimensions_mod,     only : nlev, nlevp, nv
    use physical_constants, only : rgas, kappa

    implicit none

    real(kind=real_kind),  intent(in) :: dt
    real(kind=real_kind),  intent(in) :: dx
    real(kind=real_kind),  intent(in) :: dy
    real(kind=real_kind),  intent(in) :: rdx
    real(kind=real_kind),  intent(in) :: rdy
    type (derivative_t),   intent(in) :: deriv

    integer, intent(in)               :: np1
    integer, intent(in)               :: n0
    integer, intent(in)               :: nm1
    integer, intent(in)               :: nprlev

    real(kind=real_kind), intent(in), dimension(nlevp)      :: hyai
    real(kind=real_kind), intent(in), dimension(nlevp)      :: hybi
    real(kind=real_kind), intent(in), dimension(nlev)       :: hyam
    real(kind=real_kind), intent(in), dimension(nlev)       :: hybm
    real(kind=real_kind), intent(in), dimension(nlev)       :: hybd    ! Hybrid dB coefficient on mid levels

    real(kind=real_kind), intent(in)                        :: ps0
    real(kind=real_kind), intent(in)                        :: psref
    real(kind=real_kind), intent(in), dimension(nlev)       :: Tref
    real(kind=real_kind), intent(in), dimension(nlev)       :: RTref
    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Href
    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Tmat
    real(kind=real_kind), intent(in), dimension(nlev)       :: Pvec

    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Emat
    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Emat_inv
    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Amat_inv
    real(kind=real_kind), intent(in), dimension(nlev)       :: Lambda

    real(kind=real_kind), intent(in), dimension(2,2,nv,nv)  :: met
    real(kind=real_kind), intent(in), dimension(nv,nv)      :: metdet
    real(kind=real_kind), intent(in), dimension(nv,nv)      :: rmetdet
    real(kind=real_kind), intent(in), dimension(nv,nv)      :: fcor

    real(kind=real_kind), intent(inout), dimension(nv,nv)          :: phis
    real(kind=real_kind), intent(inout), dimension(nv,nv,3)        :: lnps
    real(kind=real_kind), intent(inout), dimension(nv,nv,nlev,3)   :: T
    real(kind=real_kind), intent(inout), dimension(nv,nv,2,nlev,3) :: v

    real(kind=real_kind), intent(inout), dimension(nv,nv,nlev,3) :: div
    real(kind=real_kind), intent(inout), dimension(nv,nv,nlev)   :: zeta
    real(kind=real_kind), intent(inout), dimension(nv,nv,2)      :: grad_lnps


    real(kind=real_kind), intent(inout), dimension(nv,nv)        :: Pscript
    real(kind=real_kind), intent(inout), dimension(nv,nv,nlev)   :: Tscript
    real(kind=real_kind), intent(inout), dimension(nv,nv,2,nlev) :: Vscript
    real(kind=real_kind), intent(inout), dimension(nv,nv,nlev)   :: B
    real(kind=real_kind), intent(inout), dimension(nv,nv,2,nlev) :: grad_dGref

    ! ==============================
    ! Local Variables
    ! ==============================

    real(kind=real_kind), dimension(nv,nv)      :: ps                ! surface pressure
    real(kind=real_kind), dimension(nv,nv)      :: rps               ! 1/ps
    real(kind=real_kind), dimension(nv,nv,nlev) :: rpmid             ! 1./pmid
    real(kind=real_kind), dimension(nv,nv,nlev) :: omegap            ! omega/p
    real(kind=real_kind), dimension(nv,nv,nlev) :: rpdel             ! 1./pdel

    real(kind=real_kind), dimension(nv,nv,nlevp) :: eta_dot_dp_deta       ! eta dot * dp/deta at time level n
    real(kind=real_kind), dimension(nv,nv,nlev)  :: vgrad_ps              ! ps*(v.grad(lnps))

    real(kind=real_kind), dimension(nv,nv,nlev)   :: T_vadv
    real(kind=real_kind), dimension(nv,nv,2,nlev) :: v_vadv

    real(kind=real_kind), dimension(nv,nv)      :: HT
    real(kind=real_kind), dimension(nv,nv)      :: HrefT
    real(kind=real_kind), dimension(nv,nv)      :: HrefTm1

    real(kind=real_kind), dimension(nv,nv)      :: Gref
    real(kind=real_kind), dimension(nv,nv)      :: Grefm1
    real(kind=real_kind), dimension(nv,nv)      :: E
    real(kind=real_kind), dimension(nv,nv)      :: Phi
    real(kind=real_kind), dimension(nv,nv)      :: dGref
    real(kind=real_kind), dimension(nv,nv)      :: tmp
    real(kind=real_kind), dimension(nv,nv,2,nlev) :: vtmp
    real(kind=real_kind), dimension(nv,nv)      :: suml


    real(kind=real_kind), dimension(nv,nv,2)    :: vco
    real(kind=real_kind), dimension(nv,nv,2)    :: gradT
    real(kind=real_kind), dimension(nv,nv,2)    :: grad_Phi

    real(kind=real_kind) :: pintref(nlevp)
    real(kind=real_kind) :: pdelref(nlev)
    real(kind=real_kind) :: pmidref(nlev)
    real(kind=real_kind) :: rpdelref(nlev)
    real(kind=real_kind) :: rpmidref(nlev)

    real(kind=real_kind) :: pint(nv,nv,nlevp)
    real(kind=real_kind) :: pdel(nv,nv,nlev)
    real(kind=real_kind) :: pmid(nv,nv,nlev)

    real(kind=real_kind) :: rdt,dt2
    real(kind=real_kind) :: rpsref
    real(kind=real_kind) :: ddiv
    real(kind=real_kind) :: vgradT
    real(kind=real_kind) :: hybfac
    real(kind=real_kind) :: Crkk
    real(kind=real_kind) :: v1,v2
    real(kind=real_kind) :: term

    integer :: i,j,k,l

    rdt = 1.0_real_kind/dt
    dt2 = 2.0_real_kind*dt

    ! ==========================================================
    ! This procedure assumes grad(lnps) 
    ! has been calculated and assembled...
    ! ==========================================================

    ! ==========================================================
    ! Compute all reference pressure values
    ! ==========================================================  

    do k=1,nlevp
       pintref(k)  = hyai(k)*ps0 + hybi(k)*psref
    end do

    do k=1,nlev
       pmidref(k)  = hyam(k)*ps0 + hybm(k)*psref
       pdelref(k)  = pintref(k+1) - pintref(k)
       rpmidref(k) = 1.0_real_kind/pmidref(k)
       rpdelref(k) = 1.0_real_kind/pdelref(k)
    end do
    rpsref   = 1.0_real_kind/psref

    ! ==========================================================
    ! Compute all necessary pressure values
    ! ==========================================================  

    ps(:,:) = EXP(lnps(:,:,n0))
    rps(:,:) = 1.0_real_kind/ps(:,:)

    call preq_pressure(ps0,  ps,               &
         hyai, hybi, hyam, hybm, &
         pint, pmid, pdel)

    rpmid = 1.0_real_kind/pmid
    rpdel = 1.0_real_kind/pdel

#if 0
    ! print *, "ps=",MINVAL(ps(:,:)),MAXVAL(ps(:,:))
#endif

    ! ========================================================
    ! Compute omega/p..
    ! ========================================================

    do k=1,nlev
       do j=1,nv
          do i=1,nv
             vgrad_ps(i,j,k) = ps(i,j)*(v(i,j,1,k,n0)*grad_lnps(i,j,1) + v(i,j,2,k,n0)*grad_lnps(i,j,2))
          end do
       end do
    end do

    call preq_omegap(div(:,:,:,n0)   , vgrad_ps, pdel    ,rpmid    , &
         hybm    ,hybd   , omegap)

    ! ====================================================================
    ! Compute etadot*(dp/deta) and Pscript
    ! ====================================================================

    Pscript(:,:)           = 0.0_real_kind
    eta_dot_dp_deta(:,:,1) = 0.0_real_kind
    do k=1,nlev
       do j=1,nv
          do i=1,nv
             eta_dot_dp_deta(i,j,k+1) = eta_dot_dp_deta(i,j,k) + vgrad_ps(i,j,k)*hybd(k) + div(i,j,k,n0)*pdel(i,j,k)
             ddiv = div(i,j,k,n0) - 0.5_real_kind*div(i,j,k,nm1)
             Pscript(i,j) = Pscript(i,j) + ddiv*pdelref(k)
          end do
       end do
    end do

    do j=1,nv
       do i=1,nv
          Pscript(i,j) = lnps(i,j,nm1) + dt2*( rpsref*Pscript(i,j) - rps(i,j)*eta_dot_dp_deta(i,j,nlev+1) )
       end do
    end do

    do k=1,nlev-1
       do j=1,nv
          do i=1,nv
             eta_dot_dp_deta(i,j,k+1) = hybi(k+1)*eta_dot_dp_deta(i,j,nlev+1) - eta_dot_dp_deta(i,j,k+1)
          end do
       end do
    end do
    eta_dot_dp_deta(:,:,nlev+1) = 0.0_real_kind

    ! ====================================================================
    ! Compute the vertical advection terms needed by Vscript and Tscript
    ! Inputs are eta_dot_dp_deta and 1/pdel
    ! ====================================================================

    call preq_vertadv(T(:,:,:,n0), v(:,:,:,:,n0), eta_dot_dp_deta, rpdel, &
         T_vadv, v_vadv)

    suml(:,:) = 0.0_real_kind

    do k=1,nlev 

       gradT(:,:,:) = gradient(T(:,:,k,n0),deriv,rdx,rdy)
       Crkk       = 0.5_real_kind
          
       do j=1,nv
          do i=1,nv
             term = Crkk*(div(i,j,k,n0) - 0.5_real_kind*div(i,j,k,nm1))*pdelref(k)
             suml(i,j)  = suml(i,j) + term
             vgradT = v(i,j,1,k,n0)*gradT(i,j,1) + v(i,j,2,k,n0)*gradT(i,j,2)
             Tscript(i,j,k) = T(i,j,k,nm1) + dt2*(- vgradT - T_vadv(i,j,k)           &
                  + kappa*(T(i,j,k,n0)*omegap(i,j,k) &
                  + Tref(k)*rpmidref(k)*suml(i,j)))
             suml(i,j)  = suml(i,j) + term
          end do
       end do
    end do

    ! ===================================================
    ! Calculate B, which depends on Tscript and Vscript...
    ! Assumes Tscript and Pscript are already calculated.
    ! ===================================================

    do k=1,nlev

       do j=1,nv
          do i=1,nv
             HT(i,j)          = 0.0_real_kind
             HrefT(i,j)       = 0.0_real_kind
             HrefTm1(i,j)     = 0.0_real_kind
          end do
       end do

       do l=k,nlev
          do j=1,nv
             do i=1,nv
                HrefT(i,j)       = HrefT(i,j)       + Href(l,k)*T(i,j,l,n0)
                HrefTm1(i,j)     = HrefTm1(i,j)     + Href(l,k)*T(i,j,l,nm1)
             end do
          end do
       end do

       ! l=k

       do j=1,nv
          do i=1,nv
             HT(i,j) = (0.5_real_kind*rpmid(i,j,k))*pdel(i,j,k)*T(i,j,k,n0)
          end do
       end do

       ! l > k

       do l=k+1,nlev
          do j=1,nv
             do i=1,nv
                HT(i,j) = HT(i,j) + rpmid(i,j,l)*pdel(i,j,l)*T(i,j,l,n0)
             end do
          end do
       end do

#if 0
       print *, "RHT(",k,")=",MINVAL(rair*HT(:,:)),MAXVAL(rair*HT(:,:))
#endif

       ! ============================================
       ! Compute Phi dGref and B terms
       ! ============================================

       do j=1,nv
          do i=1,nv
             v1     = v(i,j,1,k,n0)
             v2     = v(i,j,2,k,n0)

             ! covariant velocity

             vco(i,j,1) = met(1,1,i,j)*v1 + met(1,2,i,j)*v2
             vco(i,j,2) = met(2,1,i,j)*v1 + met(2,2,i,j)*v2

             E(i,j) = 0.5_real_kind*( vco(i,j,1)*v1 + vco(i,j,2)*v2 )

             Gref(i,j)   =  phis(i,j) + HrefT(i,j)   + RTref(k)*lnps(i,j,n0)
             Grefm1(i,j) =  phis(i,j) + HrefTm1(i,j) + RTref(k)*lnps(i,j,nm1)

             Phi(i,j)    =  E(i,j) + phis(i,j) + Rgas*HT(i,j)
             dGref(i,j)  =  -(Gref(i,j)  - 0.5_real_kind*Grefm1(i,j))  ! minus sign flipped here by weak gradient
          end do
       end do

       zeta(:,:,k)         = vorticity(vco,deriv,rdx,rdy)
       grad_Phi(:,:,:)     = gradient(Phi,deriv,rdx,rdy) 
       grad_dGref(:,:,:,k) = gradient_wk(dGref,deriv,dx,dy)

#if 0
       print *, "grad_Phi()=",MINVAL(grad_Phi(:,:,:)),MAXVAL(grad_Phi(:,:,:))
#endif

       do j=1,nv
          do i=1,nv

             ! ====================================
             ! Construct covariant terms in Vscript
             ! ====================================

             v1 = - v_vadv(i,j,1,k)
             v2 = - v_vadv(i,j,2,k)

             vco(i,j,1) = met(1,1,i,j)*v1 + met(1,2,i,j)*v2
             vco(i,j,2) = met(2,1,i,j)*v1 + met(2,2,i,j)*v2

             zeta(i,j,k) = rmetdet(i,j)*zeta(i,j,k)
             hybfac =  hybm(k)*(ps(i,j)*rpmid(i,j,k))
#if 1
             Vscript(i,j,1,k) =      (   vco(i,j,1) +  v(i,j,2,k,n0)*metdet(i,j)*(fcor(i,j) + zeta(i,j,k)) &
                  - grad_Phi(i,j,1) - (Rgas*hybfac)*T(i,j,k,n0)*grad_lnps(i,j,1))

             Vscript(i,j,2,k) =      (   vco(i,j,2) - v(i,j,1,k,n0)*metdet(i,j)*(fcor(i,j) + zeta(i,j,k))  &
                  - grad_Phi(i,j,2) - (Rgas*hybfac)*T(i,j,k,n0)*grad_lnps(i,j,2))
#else
             Vscript(i,j,1,k) =      (   vco(i,j,1)                                                        &
                  - grad_Phi(i,j,1) - (Rgas*hybfac)*T(i,j,k,n0)*grad_lnps(i,j,1))

             Vscript(i,j,2,k) =      (   vco(i,j,2)                                                        &
                  - grad_Phi(i,j,2) - (Rgas*hybfac)*T(i,j,k,n0)*grad_lnps(i,j,2))
#endif
          end do
       end do

    end do

  end subroutine preq_impsys

  subroutine preq_bldrhs(dt,rdx,rdy,deriv,metdet,rmetdet,      &
       Href, Emat_inv, Lambda, RTref,        &
       phis, Pscript, Tscript, Vscript,      &
       B,C)
    use kinds,              only : real_kind
    use derivative_mod,     only : derivative_t, divergence
    use dimensions_mod,     only : nlev, nv

    implicit none

    real(kind=real_kind),  intent(in) :: dt
    real(kind=real_kind),  intent(in) :: rdx
    real(kind=real_kind),  intent(in) :: rdy
    type (derivative_t),   intent(in) :: deriv

    real(kind=real_kind), intent(in), dimension(nv,nv)        :: metdet
    real(kind=real_kind), intent(in), dimension(nv,nv)        :: rmetdet

    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Href
    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Emat_inv
    real(kind=real_kind), intent(in), dimension(nlev)       :: Lambda
    real(kind=real_kind), intent(in), dimension(nlev)       :: RTref

    real(kind=real_kind), intent(in), dimension(nv,nv)        :: phis
    real(kind=real_kind), intent(in), dimension(nv,nv)        :: Pscript
    real(kind=real_kind), intent(in), dimension(nv,nv,nlev)   :: Tscript
    real(kind=real_kind), intent(in), dimension(nv,nv,2,nlev) :: Vscript

    real(kind=real_kind), intent(inout), dimension(nv,nv,nlev)  :: B
    real(kind=real_kind), intent(out), dimension(nv,nv,nlev)  :: C        ! right hand side of Helmholtz problem

    ! Local Variables

    real(kind=real_kind), dimension(nv,nv,2)    :: gVscript
    real(kind=real_kind), dimension(nv,nv,nlev) :: div_Vscript
    real(kind=real_kind), dimension(nv,nv)      :: HrefTscript
    real(kind=real_kind) :: rdt

    integer :: i,j,k,l

    do k=1,nlev

       do j=1,nv
          do i=1,nv
             HrefTscript(i,j) = 0.0_real_kind
          end do
       end do

       do l=k,nlev
          do j=1,nv
             do i=1,nv
                HrefTscript(i,j) = HrefTscript(i,j) + Href(l,k)*Tscript(i,j,l)
             end do
          end do
       end do

       do j=1,nv
          do i=1,nv
             B(i,j,k)    =  phis(i,j) + HrefTscript(i,j) + RTref(k)*Pscript(i,j)
          end do
       end do

    end do

    do k=1,nlev
       do j=1,nv
          do i=1,nv
             gVscript(i,j,1) = metdet(i,j)*Vscript(i,j,1,k)
             gVscript(i,j,2) = metdet(i,j)*Vscript(i,j,2,k)
          end do
       end do
       div_Vscript(:,:,k) = divergence(gVscript(:,:,:),deriv,rdx,rdy)
    end do

    rdt = 1.0_real_kind/dt

    do k=1,nlev

       do j=1,nv
          do i=1,nv
             C(i,j,k) = 0.0_real_kind
          end do
       end do

       do l=1,nlev
          do j=1,nv
             do i=1,nv
                C(i,j,k) = C(i,j,k) + Emat_inv(l,k)*B(i,j,l)
             end do
          end do
       end do

       B(:,:,k) = C(:,:,k)

       do j=1,nv
          do i=1,nv
             C(i,j,k) = metdet(i,j)*C(i,j,k)
          end do
       end do

#if 1
       do l=1,nlev
          do j=1,nv
             do i=1,nv
                C(i,j,k) = C(i,j,k) - rdt*Lambda(k)*Emat_inv(l,k)*div_Vscript(i,j,l)
             end do
          end do
       end do
#endif

    end do

  end subroutine preq_bldrhs
  subroutine preq_backsub(dt,dx,dy,deriv,metdet,rmetdet,      &
       Tref,Href,RTref,Emat,Ainv,Lambda,Tmat,Pvec,     &
       Pscript,Tscript,Vscript,Gamma_ref,B, &
       phis,lnps,T,grad_Gref,D)
    use kinds,              only : real_kind
    use derivative_mod,     only : derivative_t, gradient_wk
    use dimensions_mod,     only : nlev, nv

    implicit none

    real(kind=real_kind),  intent(in) :: dt
    real(kind=real_kind),  intent(in) :: dx
    real(kind=real_kind),  intent(in) :: dy
    type (derivative_t),   intent(in) :: deriv

    real(kind=real_kind), intent(in), dimension(nv,nv)        :: rmetdet
    real(kind=real_kind), intent(in), dimension(nv,nv)        :: metdet

    real(kind=real_kind), intent(in), dimension(nlev)       :: Tref
    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Href
    real(kind=real_kind), intent(in), dimension(nlev)       :: RTref
    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Emat
    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Ainv
    real(kind=real_kind), intent(in), dimension(nlev)       :: Lambda
    real(kind=real_kind), intent(in), dimension(nlev,nlev)  :: Tmat
    real(kind=real_kind), intent(in), dimension(nlev)       :: Pvec

    real(kind=real_kind), intent(in), dimension(nv,nv)        :: Pscript
    real(kind=real_kind), intent(in), dimension(nv,nv,nlev)   :: Tscript
    real(kind=real_kind), intent(in), dimension(nv,nv,2,nlev) :: Vscript
    real(kind=real_kind), intent(inout), dimension(nv,nv,nlev)   :: Gamma_ref
    real(kind=real_kind), intent(in), dimension(nv,nv,nlev)   :: B

    real(kind=real_kind), intent(out), dimension(nv,nv)        :: phis
    real(kind=real_kind), intent(out), dimension(nv,nv)        :: lnps
    real(kind=real_kind), intent(out), dimension(nv,nv,nlev)   :: T
    real(kind=real_kind), intent(out), dimension(nv,nv,2,nlev) :: grad_Gref
    real(kind=real_kind), intent(out), dimension(nv,nv,nlev)   :: D

    ! Local Variables

    real(kind=real_kind), dimension(nv,nv)      :: PD
    real(kind=real_kind), dimension(nv,nv,nlev) :: Gref

    real(kind=real_kind), dimension(nv,nv)      :: HrefT
    real(kind=real_kind), dimension(nv,nv,nlev) :: Grefp1

    real(kind=real_kind) :: rdt
    integer :: i,j,k,l

    rdt = 1.0D0/dt

    do k=1,nlev
#if 0
       do j=1,nv
          do i=1,nv
             Gref(i,j,k)=B(i,j,k)  !!! RDL =0.0D0
          end do
       end do
#else
       do j=1,nv
          do i=1,nv
             Gref(i,j,k)=0.0_real_kind
          end do
       end do
       do l=1,nlev
          do j=1,nv
             do i=1,nv
                Gref(i,j,k) = Gref(i,j,k) + Emat(l,k)*Gamma_ref(i,j,l)
             end do
          end do
       end do
#endif
#if 0
       print *, "Gref(",k,")=",MINVAL(Gref(:,:,k)),MAXVAL(Gref(:,:,k))
       print *, "B(",k,")=",MINVAL(B(:,:,k)),MAXVAL(B(:,:,k))
       print *
#endif
    end do

    do j=1,nv
       do i=1,nv
          PD(i,j)=0.0_real_kind
       end do
    end do

    do k=1,nlev
       do j=1,nv
          do i=1,nv
             D(i,j,k)=0.0_real_kind
          end do
       end do
#if 1
       do l=1,nlev
          do j=1,nv
             do i=1,nv
                D(i,j,k) = D(i,j,k) + dt*Emat(l,k)*(B(i,j,l) - Gamma_ref(i,j,l))/Lambda(l)
             end do
          end do
       end do
#endif

       do j=1,nv
          do i=1,nv
             PD(i,j) = PD(i,j) + Pvec(k)*D(i,j,k)
          end do
       end do
    end do

#if 0
    print *, "D=",MINVAL(D),MAXVAL(D)
#endif

    do k=1,nlev
       do j=1,nv
          do i=1,nv
             T(i,j,k)=Tscript(i,j,k)
          end do
       end do
#if 1
       do l=1,nlev
          do j=1,nv
             do i=1,nv
                T(i,j,k) = T(i,j,k) -  dt*Tmat(l,k)*D(i,j,l)
             end do
          end do
       end do
#endif
       grad_Gref(:,:,:,k)=gradient_wk(Gamma_ref(:,:,k),deriv,dx,dy)
    end do

    do j=1,nv
       do i=1,nv
          lnps(i,j) = Pscript(i,j) - dt*PD(i,j)
       end do
    end do

  end subroutine preq_backsub


  subroutine preq_vertadv(T, v, eta_dot_dp_deta, rpdel, &
       T_vadv, v_vadv)
    use kinds,              only : real_kind
    use dimensions_mod,     only : nlev, nv, nlevp
    implicit none
    
    real (kind=real_kind), intent(in) :: T(nv,nv,nlev)
    real (kind=real_kind), intent(in) :: v(nv,nv,2,nlev)
    real (kind=real_kind), intent(in) :: eta_dot_dp_deta(nv,nv,nlevp)
    real (kind=real_kind), intent(in) :: rpdel(nv,nv,nlev)

    real (kind=real_kind), intent(out) :: T_vadv(nv,nv,nlev)
    real (kind=real_kind), intent(out) :: v_vadv(nv,nv,2,nlev)

    ! ========================
    ! Local Variables
    ! ========================

    integer :: i,j,k
    real (kind=real_kind) :: facp, facm

    ! ===========================================================
    ! Compute vertical advection of T and v from eq. (3.b.1)
    !
    ! k = 1 case:
    ! ===========================================================

    k=1
    do j=1,nv
       do i=1,nv 
          facp            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)
          T_vadv(i,j,k)   = facp*(T(i,j,k+1)- T(i,j,k))
          v_vadv(i,j,1,k) = facp*(v(i,j,1,k+1)- v(i,j,1,k))
          v_vadv(i,j,2,k) = facp*(v(i,j,2,k+1)- v(i,j,2,k))
       end do
    end do

    ! ===========================================================
    ! vertical advection
    !
    ! 1 < k < nlev case:
    ! ===========================================================

    do k=2,nlev-1

       do j=1,nv
          do i=1,nv
             facp            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)
             facm            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)
             T_vadv(i,j,k)   = facp*(T(i,j,k+1)- T(i,j,k)) + &
                  facm*(T(i,j,k)- T(i,j,k-1))
             v_vadv(i,j,1,k) = facp*(v(i,j,1,k+1)- v(i,j,1,k)) + &
                  facm*(v(i,j,1,k)- v(i,j,1,k-1))
             v_vadv(i,j,2,k) = facp*(v(i,j,2,k+1)- v(i,j,2,k)) + &
                  facm*(v(i,j,2,k)- v(i,j,2,k-1))
          end do
       end do

    end do

    ! ===========================================================
    ! vertical advection
    !
    ! k = nlev case:
    ! ===========================================================

    k=nlev

    do j=1,nv
       do i=1,nv
          facm            = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)
          T_vadv(i,j,k)   = facm*(T(i,j,k)- T(i,j,k-1))
          v_vadv(i,j,1,k) = facm*(v(i,j,1,k)- v(i,j,1,k-1))
          v_vadv(i,j,2,k) = facm*(v(i,j,2,k)- v(i,j,2,k-1))
       end do
    end do

  end subroutine preq_vertadv




!----------------------------------------------------------------------- 
! preq_omegap:

! Purpose: 
! Calculate (omega/p) needed for the Thermodynamics Equation
! 
! Method: 
! Simplified version in CAM2 for clarity
! 
! Author: Modified by Rich Loft for use in HOMME. 
! 
!-----------------------------------------------------------------------

  subroutine preq_omegap(div     ,vgrad_ps,pdel    ,rpmid, &
       hybm    ,hybd    ,omegap   )
    use kinds, only : real_kind
    use dimensions_mod, only : nv, nlev
    implicit none


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(in) :: div(nv,nv,nlev)      ! divergence
    real(kind=real_kind), intent(in) :: vgrad_ps(nv,nv,nlev) ! v.grad(ps)
    real(kind=real_kind), intent(in) :: pdel(nv,nv,nlev)     ! layer thicknesses (pressure)
    real(kind=real_kind), intent(in) :: rpmid(nv,nv,nlev)    ! 1./pmid
    real(kind=real_kind), intent(in) :: hybm(nlev)           ! Hybrid B coefficient on mid levels
    real(kind=real_kind), intent(in) :: hybd(nlev)           ! Hybrid dB coefficient on mid levels
    real(kind=real_kind), intent(out):: omegap(nv,nv,nlev)   ! vertical pressure velocity
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) term             ! one half of basic term in omega/p summation 
    real(kind=real_kind) Ckk              ! diagonal term of energy conversion matrix
    real(kind=real_kind) suml(nv,nv)      ! partial sum over l = (1, k-1)
    !-----------------------------------------------------------------------

    ! =========================
    ! Zero partial sum
    ! =========================

    do j=1,nv
       do i=1,nv
          suml(i,j)=0.0_real_kind
       end do
    end do

    ! =============================
    ! Compute omegap 
    ! =============================

    do k=1,nlev
       do j=1,nv
          do i=1,nv
             Ckk       = 0.5_real_kind
             term      = Ckk*(div(i,j,k)*pdel(i,j,k) + vgrad_ps(i,j,k)*hybd(k))
             suml(i,j) = suml(i,j) + term
             omegap(i,j,k) = rpmid(i,j,k)*(hybm(k)*vgrad_ps(i,j,k) - suml(i,j))
             suml(i,j) = suml(i,j) + term
          end do
       end do
    end do

  end subroutine preq_omegap



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  compute omega/p using ps, modeled after CCM3 formulas 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine preq_omega_ps(omega_p,hvcoord,ps,p,vgrad_ps,divdp)
    use kinds, only : real_kind
    use dimensions_mod, only : nv, nlev
    use hybvcoord_mod, only : hvcoord_t
    implicit none


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(in) :: divdp(nv,nv,nlev)      ! divergence
    real(kind=real_kind), intent(in) :: vgrad_ps(nv,nv,nlev) ! v.grad(ps)
    real(kind=real_kind), intent(in) :: p(nv,nv,nlev)     ! layer thicknesses (pressure)
    real(kind=real_kind), intent(in) :: ps(nv,nv)
    type (hvcoord_t),     intent(in) :: hvcoord
    real(kind=real_kind), intent(out):: omega_p(nv,nv,nlev)   ! vertical pressure velocity
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) term             ! one half of basic term in omega/p summation 
    real(kind=real_kind) Ckk,Ckl          ! diagonal term of energy conversion matrix
    real(kind=real_kind) suml(nv,nv)      ! partial sum over l = (1, k-1)
    !-----------------------------------------------------------------------
       do j=1,nv	
          do i=1,nv
             ckk = 0.5d0/p(i,j,1)
             term = divdp(i,j,1)
             omega_p(i,j,1) = hvcoord%hybm(1)*vgrad_ps(i,j,1)/p(i,j,1)
             omega_p(i,j,1) = omega_p(i,j,1) - ckk*term
             suml(i,j) = term
          end do
       end do

       do k=2,nlev-1

	  do j=1,nv
             do i=1,nv
                ckk = 0.5d0/p(i,j,k)
                ckl = 2*ckk
                term = divdp(i,j,k)
                omega_p(i,j,k) = hvcoord%hybm(k)*vgrad_ps(i,j,k)/p(i,j,k)
                omega_p(i,j,k) = omega_p(i,j,k) - ckl*suml(i,j) - ckk*term
                suml(i,j) = suml(i,j) + term

             end do
          end do

       end do

       do j=1,nv
          do i=1,nv
             ckk = 0.5d0/p(i,j,nlev)
             ckl = 2*ckk
             term = divdp(i,j,nlev)
             omega_p(i,j,nlev) = hvcoord%hybm(nlev)*vgrad_ps(i,j,nlev)/p(i,j,nlev)
             omega_p(i,j,nlev) = omega_p(i,j,nlev) - ckl*suml(i,j) - ckk*term
          end do
       end do

  end subroutine preq_omega_ps




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  compute omega/p using lnps 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine preq_omega_lnps(omega_p,hvcoord,ps,p,dp,vgrad_lnps,div)
    use kinds, only : real_kind
    use dimensions_mod, only : nv, nlev
    use hybvcoord_mod, only : hvcoord_t
    implicit none


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(in) :: div(nv,nv,nlev)      ! divergence
    real(kind=real_kind), intent(in) :: vgrad_lnps(nv,nv,nlev) ! v.grad(ps)
    real(kind=real_kind), intent(in) :: p(nv,nv,nlev)      ! pressure
    real(kind=real_kind), intent(in) :: dp(nv,nv,nlev)     ! dp/dn
    real(kind=real_kind), intent(in) :: ps(nv,nv)
    type (hvcoord_t),     intent(in) :: hvcoord
    real(kind=real_kind), intent(out):: omega_p(nv,nv,nlev)   ! vertical pressure velocity
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) term             ! one half of basic term in omega/p summation 
    real(kind=real_kind) Ckk,Ckl          ! diagonal term of energy conversion matrix
    real(kind=real_kind) suml(nv,nv)      ! partial sum over l = (1, k-1)
    !-----------------------------------------------------------------------
       do j=1,nv	
          do i=1,nv
             ckk = 0.5d0/p(i,j,1)
             term = div(i,j,1)*dp(i,j,1) + vgrad_lnps(i,j,1)*ps(i,j)*hvcoord%hybd(1)
             omega_p(i,j,1) = hvcoord%hybm(1)*(ps(i,j)/p(i,j,1))*vgrad_lnps(i,j,1)
             omega_p(i,j,1) = omega_p(i,j,1) - ckk*term
             suml(i,j) = term
          end do
       end do

       do k=2,nlev-1

	  do j=1,nv
             do i=1,nv
                ckk = 0.5d0/p(i,j,k)
                ckl = 2*ckk
                term = div(i,j,k)*dp(i,j,k) + vgrad_lnps(i,j,k)*ps(i,j)*hvcoord%hybd(k)
                omega_p(i,j,k) = hvcoord%hybm(k)*(ps(i,j)/p(i,j,k))*vgrad_lnps(i,j,k)
                omega_p(i,j,k) = omega_p(i,j,k) - ckl*suml(i,j) - ckk*term
                suml(i,j) = suml(i,j) + term

             end do
          end do

       end do

       do j=1,nv
          do i=1,nv
             ckk = 0.5d0/p(i,j,nlev)
             ckl = 2*ckk
             term = div(i,j,nlev)*dp(i,j,nlev) + vgrad_lnps(i,j,nlev)*ps(i,j)*hvcoord%hybd(nlev)
             omega_p(i,j,nlev) = hvcoord%hybm(nlev)*(ps(i,j)/p(i,j,nlev))*vgrad_lnps(i,j,nlev)
             omega_p(i,j,nlev) = omega_p(i,j,nlev) - ckl*suml(i,j) - ckk*term
          end do
       end do

  end subroutine preq_omega_lnps



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  CCM3 hydrostatic integral
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine preq_hydrostatic(phi,phis,T_v,p,dp)
    use kinds, only : real_kind
    use dimensions_mod, only : nv, nlev
    use physical_constants, only : rgas
!    use hybvcoord_mod, only : hvcoord_t
    implicit none


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(out) :: phi(nv,nv,nlev)     
    real(kind=real_kind), intent(in) :: phis(nv,nv)
    real(kind=real_kind), intent(in) :: T_v(nv,nv,nlev)
    real(kind=real_kind), intent(in) :: p(nv,nv,nlev)   
    real(kind=real_kind), intent(in) :: dp(nv,nv,nlev)  
 !   type (hvcoord_t),     intent(in) :: hvcoord
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) term             ! one half of basic term in omega/p summation 
    real(kind=real_kind) Hkk,Hkl          ! diagonal term of energy conversion matrix
    real(kind=real_kind), dimension(nv,nv,nlev) :: phii       ! Geopotential at interfaces
    !-----------------------------------------------------------------------
       do j=1,nv
          do i=1,nv
             hkk = dp(i,j,nlev)*0.5d0/p(i,j,nlev)
             hkl = 2*hkk
             phii(i,j,nlev)  = Rgas*T_v(i,j,nlev)*hkl
             phi(i,j,nlev) = phis(i,j) + Rgas*T_v(i,j,nlev)*hkk 
          end do
       end do

       do k=nlev-1,2,-1

          do j=1,nv	
             do i=1,nv
                ! hkk = dp*ckk
                hkk = dp(i,j,k)*0.5d0/p(i,j,k)
                hkl = 2*hkk
                phii(i,j,k) = phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkl
                phi(i,j,k) = phis(i,j) + phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkk
             end do
          end do

       end do

       do j=1,nv	
          do i=1,nv
             ! hkk = dp*ckk
             hkk = 0.5d0*dp(i,j,1)/p(i,j,1)
             phi(i,j,1) = phis(i,j) + phii(i,j,2) + Rgas*T_v(i,j,1)*hkk
          end do
       end do


end subroutine preq_hydrostatic



!
!  The hydrostatic routine from CAM physics.
!  (FV stuff removed)
!  t,q input changed to take t_v
!  removed gravit, so this routine returns PHI, not zm
subroutine geopotential_t(                                 &
       pmid   , pdel   ,  tv      , rair   ,  zm)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the geopotential height (above the surface) at the midpoints and 
! interfaces using the input temperatures and pressures.
!
!-----------------------------------------------------------------------
    use dimensions_mod,     only : nlev, nlevp, nv
    use kinds, only : real_kind
    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments



    real(real_kind), intent(in) :: pmid (nv*nv,nlev)    ! Midpoint pressures
    real(real_kind), intent(in) :: pdel (nv*nv,nlev)    ! layer thickness
    real(real_kind), intent(in) :: tv    (nv*nv,nlev)    ! temperature
    real(real_kind), intent(in) :: rair                 ! Gas constant for dry air
    ! real(real_kind), intent(in) :: gravit               ! Acceleration of gravity
    ! real(real_kind), intent(in) :: zvir                 ! rh2o/rair - 1

! Output arguments

    real(real_kind), intent(out) :: zm(nv*nv,nlev)      ! Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
    integer :: ncol=nv*nv             ! Number of longitudes

    integer  :: i,k                ! Lon, level indices
    real(real_kind) :: hkk(nv*nv)         ! diagonal element of hydrostatic matrix
    real(real_kind) :: hkl(nv*nv)         ! off-diagonal element
    real(real_kind) :: rog                ! Rair / gravit
    real(real_kind) :: zi(nv*nv,nlevp)     ! Height above surface at interfaces
!
!-----------------------------------------------------------------------
!
!    rog = rair/gravit
    rog = rair

! The surface height is zero by definition.
    do i = 1,ncol
       zi(i,nlevp) = 0.0_real_kind
    end do

! Compute zi, zm from bottom up. 
! Note, zi(i,k) is the interface above zm(i,k)
    do k = nlev, 1, -1
! First set hydrostatic elements consistent with dynamics
       do i = 1,ncol
          hkl(i) = pdel(i,k) / pmid(i,k)
          hkk(i) = 0.5_real_kind * hkl(i)
       end do

! Now compute tv, zm, zi
       do i = 1,ncol
          ! tvfac   = 1._r8 + zvir * q(i,k)
          ! tv      = t(i,k) * tvfac
          zm(i,k) = zi(i,k+1) + rog * tv(i,k) * hkk(i)
          zi(i,k) = zi(i,k+1) + rog * tv(i,k) * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_t





!----------------------------------------------------------------------- 
! preq_pressure:
!
! Purpose: 
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure. Originally plevs0!
! 
! Method: 
! 
! Author: B. Boville/ Adapted for HOMME by Rich Loft
! 
!-----------------------------------------------------------------------
!
! $Id: prim_si_mod.F90,v 2.10 2005/10/14 20:17:22 jedwards Exp $
! $Author: jedwards $
!
!-----------------------------------------------------------------------

  subroutine preq_pressure (ps0,  ps,               &
       hyai, hybi, hyam, hybm, &
       pint, pmid, pdel)
    use kinds, only : real_kind
    use dimensions_mod, only : nv, nlev, nlevp
    implicit none

    !-----------------------------------------------------------------------

    real(kind=real_kind), intent(in)  :: ps0                ! Hybrid coordinate reference pressure (pascals)
    real(kind=real_kind), intent(in)  :: ps(nv,nv)          ! Surface pressure (pascals)
    real(kind=real_kind), intent(in)  :: hyai(nlevp)        ! Hybrid interface A coefficients
    real(kind=real_kind), intent(in)  :: hybi(nlevp)        ! Hybrid interface B coefficients
    real(kind=real_kind), intent(in)  :: hyam(nlev)         ! Hybrid midpoint  A coefficients
    real(kind=real_kind), intent(in)  :: hybm(nlev)         ! Hybrid midpoint  B coefficients
    real(kind=real_kind), intent(out) :: pint(nv,nv,nlevp)  ! Pressure at model interfaces
    real(kind=real_kind), intent(out) :: pmid(nv,nv,nlev)   ! Pressure at model levels
    real(kind=real_kind), intent(out) :: pdel(nv,nv,nlev)   ! Layer thickness (pint(k+1) - pint(k))
    !-----------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k             ! Horizontal, level indices
    !-----------------------------------------------------------------------
    !
    ! Set interface pressures
    !
    do k=1,nlevp
       do j=1,nv
          do i=1,nv
             pint(i,j,k) = hyai(k)*ps0 + hybi(k)*ps(i,j)
          end do
       end do
    end do
    !
    ! Set midpoint pressures and layer thicknesses
    !
    do k=1,nlev
       do j=1,nv
          do i=1,nv
             pmid(i,j,k) = hyam(k)*ps0 + hybm(k)*ps(i,j)
             pdel(i,j,k) = pint(i,j,k+1) - pint(i,j,k)
          end do
       end do
    end do

  end subroutine preq_pressure




end module prim_si_mod
