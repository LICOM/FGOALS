#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#undef _DGEMV
module solver_mod
  use kinds, only : real_kind, int_kind
  use dimensions_mod, only : npsq, nlev
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  implicit none
  private

  character(len=8), private, parameter :: blkjac_storage = "inverse"
  !  character(len=8), private, parameter :: blkjac_storage = "LUfactor"

  type, public :: blkjac_t
     sequence
     real (kind=real_kind), dimension(npsq,npsq,nlev) :: E
     integer(kind=int_kind),     dimension(npsq,nlev) :: ipvt
  end type blkjac_t


  public  :: pcg_solver
  public  :: blkjac_init

  interface pcg_solver
     module procedure pcg_solver_stag
     module procedure pcg_solver_nonstag
  end interface

contains

  function pcg_solver_stag(elem,  & 
       rhs,      &
       cg,       &
       red,      &
       edge2,    &   
       lambdasq, &   
       deriv,    &   
       nets,     & 
       nete,     &
       blkjac) result(x) 
    use dimensions_mod, only : nlev, nv, np, nvsq
    use element_mod, only : element_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use cg_mod, only : cg_t, congrad
    use edge_mod, only : edgebuffer_t, edgevpack, edgerotate, edgevunpack
    use derivative_mod, only : derivative_t, gradient_wk, gradient, divergence, derivative_stag_t
    use control_mod, only : maxits, while_iter, tol, precon_method
    use physical_constants, only : rearth
    use bndry_mod, only : bndry_exchangeV
    use linear_algebra_mod, only : matvec
    use parallel_mod, only : syncmp, haltmp


    type(element_t), intent(in), target :: elem(:)
    integer, intent(in)  :: nets,nete
    real (kind=real_kind), intent(in) :: rhs(np,np,nlev,nets:nete) ! right hand side of operator
    type (cg_t)                       :: cg             ! conjugate gradient    (private)
    type (ReductionBuffer_ordered_1d_t)  :: red         ! CG reduction buffer   (shared memory)
    type (EdgeBuffer_t)               :: edge2          ! Laplacian edge buffer (shared memory)
    real (kind=real_kind)             :: lambdasq(nlev) ! Helmholtz lengthscale (private)
    type (derivative_t)          :: deriv          ! Staggered derivative struct     (private)
    type (blkjac_t)		      :: blkjac(nets:nete)
    real (kind=real_kind)             :: x(np,np,nlev,nets:nete)     ! solution (result)

#ifdef CAM
    call haltmp('semi-implicit method not yet supported in cam')
#else

    ! ===========
    ! Local
    ! ===========

    real (kind=real_kind),pointer :: metinv(:,:,:,:)
    real (kind=real_kind),pointer :: metdet(:,:)
    real (kind=real_kind),pointer :: rmv(:,:)
    real (kind=real_kind),pointer :: metdetp(:,:)
    real (kind=real_kind),pointer :: mp(:,:)

    real (kind=real_kind) :: gradp(nv,nv,2,nlev,nets:nete)
    real (kind=real_kind) :: div(np,np)
    real (kind=real_kind) :: p(np,np)
    real (kind=real_kind) :: r(npsq)
    real (kind=real_kind) :: z(npsq)

    real (kind=real_kind) :: gradp1
    real (kind=real_kind) :: gradp2
    real (kind=real_kind) :: dx,dy
    real (kind=real_kind) :: rdx,rdy
    real (kind=real_kind) :: lenscale

    integer :: ie
    integer :: i,j,k
    integer :: kptr
    integer :: iptr
    integer :: ieptr

    character(len=1) :: Trans
    real(kind=real_kind) :: one,zero
    integer(kind=int_kind) :: inc

    call t_startf('pcg_solver_stag')

    lenscale=rearth
#ifdef DEBUG
    print *,"lenscale=",lenscale
#endif

    ! ========================================
    ! Load the rhs in the cg struct...
    ! ========================================

#ifndef CAM

    !DBG print *,'pcg_solver: point #1'
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             end do
          end do
       end do
    end do
    !DBG print *,'pcg_solver: point #2'

    Trans='N'
    one=1.0D0
    zero=0.0D0
    inc=1

    ! cg%debug_level = 3
    do while (congrad(cg,red,maxits,tol))
       call t_startf('timer_helm')

!       print *,'CG inter = ',cg%iter
       do ie=nets,nete
          ieptr=ie-nets+1
          metinv => elem(ie)%metinv
          metdet => elem(ie)%metdet
          dx=0.5D0*elem(ie)%dx/lenscale    ! weak derivative element x dimension
          dy=0.5D0*elem(ie)%dy/lenscale    ! weak derivative element y dimension

          do k=1,nlev
             if (.not.cg%converged(k)) then

                ! ========================================
                ! Apply preconditioner: wrk2 = (M^-1) * wrk1 
                ! ========================================

                if (precon_method == "block_jacobi") then

                   if (blkjac_storage == "LUfactor") then
                      stop 'dgesl needs linpack'
!                      call dgesl(cg%state(ieptr)%r(:,k),cg%state(ieptr)%z(:,k),blkjac(ie)%E(:,:,k),blkjac(ie)%ipvt(:,k))
                   else if (blkjac_storage == "inverse") then
#ifdef _DGEMV
                      call dgemv(Trans,npsq,npsq,one,blkjac(ie)%E(:,:,k),npsq,cg%state(ieptr)%r(:,k),inc,zero,cg%state(ieptr)%z(:,k),inc)
#else
                      call matvec(cg%state(ieptr)%r(:,k),cg%state(ieptr)%z(:,k),blkjac(ie)%E(:,:,k),npsq)
#endif
                   end if

                else if (precon_method == "identity") then

                   do j=1,npsq
                      cg%state(ieptr)%z(j,k)=cg%state(ieptr)%r(j,k)
		   end do

                end if


                !JMD===========================================
                !JMD   2*nv*np*(nv + np) Flops 
   		!JMD  SR = (4*np*nv + 2*nv*nv + np*np)*Ld
   		!JMD  SUM(WS) = (6*np*nv + 2*nv*nv + np*np
                !JMD===========================================
#ifdef _WK_GRAD
#ifdef _NEWSTRUCT
                gradp(:,:,:,k,ie)=gradient_wk(cg%state(ieptr,k)%z(:),deriv,dx,dy)
#else
                !JPE temp solution get rid of reshape!
                gradp(:,:,:,k,ie)=gradient_wk(reshape(cg%state(ieptr)%z(:,k),(/np,np/)),deriv,dx,dy)
#endif
#else
#ifdef _NEWSTRUCT
                gradp(:,:,:,k,ie)=gradient(cg%state(ieptr,k)%z(:),deriv,dx,dy)
#else
                !JPE temp solution get rid of reshape!
                gradp(:,:,:,k,ie)=gradient(reshape(cg%state(ieptr)%z(:,k),(/np,np/)),deriv,dx,dy)
#endif
#endif

                ! =======================================
                ! rotate gradient to form contravariant
                !JMD  4*nv*nv Flops
                ! =======================================

                do j=1,nv
                   do i=1,nv
                      gradp1       = gradp(i,j,1,k,ie)
                      gradp2       = gradp(i,j,2,k,ie)
                      gradp(i,j,1,k,ie) = metdet(i,j)*(metinv(1,1,i,j)*gradp1 + &
                           metinv(1,2,i,j)*gradp2)
                      gradp(i,j,2,k,ie) = metdet(i,j)*(metinv(2,1,i,j)*gradp1 + &
                           metinv(2,2,i,j)*gradp2)
                   end do
                end do

             end if
          end do

          kptr=0
          call edgeVpack(edge2,gradp(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

          !DBG print *,'pcg_solver: point #11'
          kptr=0
          call edgerotate(edge2,2*nlev,kptr,elem(ie)%desc)

       end do
       while_iter = while_iter + 1 

       call bndry_exchangeV(cg%hybrid,edge2)

       do ie=nets,nete
          ieptr=ie-nets+1

          rdx=2.0/(elem(ie)%dx*lenscale)
          rdy=2.0/(elem(ie)%dy*lenscale)

          rmv     => elem(ie)%rmv
          metdetp => elem(ie)%metdetp
          mp      => elem(ie)%mp

          kptr=0
          call edgeVunpack(edge2, gradp(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

          do k=1,nlev
             if (.not.cg%converged(k)) then

                ! ====================
                ! 2*nv*nv Flops
                ! ====================

                do j=1,nv
                   do i=1,nv
                      gradp(i,j,1,k,ie) = rmv(i,j)*gradp(i,j,1,k,ie)
                      gradp(i,j,2,k,ie) = rmv(i,j)*gradp(i,j,2,k,ie)
                   end do
                end do

                ! ================================================
                ! Compute  Pseudo Laplacian(p), store in div
                !JMD   2*nv*np*(nv + np) Flops 
                ! ================================================

                div(:,:) = divergence(gradp(:,:,:,k,ie),deriv,rdx,rdy)

                ! ==========================================
                ! compute Helmholtz operator, store in wrk3
                !  4*np*np Flops
                ! ==========================================

                iptr=1
                do j=1,np
                   do i=1,np
                      cg%state(ieptr)%s(iptr,k) = mp(i,j)*(metdetp(i,j)*cg%state(ieptr)%z(iptr,k) + lambdasq(k)*div(i,j))
                      iptr=iptr+1
                   end do
                end do

             end if
          end do
       end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
       call t_stopf('timer_helm')

    end do  ! CG solver while loop


    ! ===============================
    ! Converged! unpack wrk3 into x
    ! ===============================

    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
#ifdef _NEWSTRUCT
                x(i,j,k,ie)=cg%state(ieptr,k)%x(iptr)
#else
                x(i,j,k,ie)=cg%state(ieptr)%x(iptr,k)
#endif
                iptr=iptr+1
             end do
          end do
       end do
    end do

#endif
    call t_startf('pcg_solver_stag')
! CAM
#endif 
  end function pcg_solver_stag

  ! ================================================
  ! pcg_solver_nonstag:
  !
  ! Preconditioned conjugate gradient solver on the
  ! Gauss-Lobatto nonstaggered grid (np = nv).
  ! 
  ! ================================================

  function pcg_solver_nonstag(elem,       &
       rhs,        &
       cg,         &
       red,        &
       edge1,      &
       edge2,      &   
       lambdasq,   &   
       deriv,      &   
       nets,       & 
       nete,       &
       blkjac) result(x) 
    use dimensions_mod, only : nlev, nv, np, nvsq
    use element_mod, only : element_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use cg_mod, only : cg_t, congrad
    use edge_mod, only : edgebuffer_t, edgevpack, edgerotate, edgevunpack
    use derivative_mod, only : derivative_t, gradient_wk, gradient, divergence
    use control_mod, only : maxits, while_iter, tol, precon_method
    use physical_constants, only : rearth
    use bndry_mod, only : bndry_exchangeV
    use linear_algebra_mod, only : matvec
    use parallel_mod, only : haltmp

    integer, intent(in)  :: nets,nete
    type(element_t), intent(in), target :: elem(:)
    real (kind=real_kind), intent(in) :: rhs(nv,nv,nlev,nets:nete) ! right hand side of operator
    type (cg_t)                       :: cg             ! conjugate gradient    (private)
    type (ReductionBuffer_ordered_1d_t)  :: red         ! CG reduction buffer   (shared memory)
    type (EdgeBuffer_t)               :: edge1          ! Laplacian divergence edge buffer (shared memory)
    type (EdgeBuffer_t)               :: edge2          ! Laplacian gradient edge buffer (shared memory)
    real (kind=real_kind), intent(in) :: lambdasq(nlev) ! Helmholtz lengthscale (private)
    type (derivative_t)               :: deriv          ! non staggered derivative struct     (private)
    type (blkjac_t)		      :: blkjac(nets:nete)

    real (kind=real_kind)             :: x(nv,nv,nlev,nets:nete)     ! solution (result)

    ! ===========
    ! Local
    ! ===========
#ifdef CAM
    call haltmp('semi-implicit method not yet supported in cam')
#else

    real (kind=real_kind),pointer :: metinv(:,:,:,:)
    real (kind=real_kind),pointer :: metdet(:,:)
    real (kind=real_kind),pointer :: metdetp(:,:)
    real (kind=real_kind),pointer :: rmetdet(:,:)
    real (kind=real_kind),pointer :: rmv(:,:)
    real (kind=real_kind),pointer :: mv(:,:)

    real (kind=real_kind) :: gradp(nv,nv,2,nlev,nets:nete)
    real (kind=real_kind) :: div(nv,nv,nlev,nets:nete)
    real (kind=real_kind) :: p(nv,nv)
    real (kind=real_kind) :: r(nvsq)
    real (kind=real_kind) :: z(nvsq)

    real (kind=real_kind) :: gradp1
    real (kind=real_kind) :: gradp2
    real (kind=real_kind) :: dx,dy
    real (kind=real_kind) :: rdx,rdy
    real (kind=real_kind) :: lenscale


    integer :: ie
    integer :: i,j,k
    integer :: kptr
    integer :: iptr
    integer :: ieptr


    lenscale=rearth

    ! ========================================
    ! Load the rhs in the cg struct...
    ! ========================================

    call t_startf('pcg_solver_nonstag')

    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,nv
             do i=1,nv
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             end do
          end do
       end do
    end do

    !cg%debug_level=1
    do while (congrad(cg,red,maxits,tol))
       call t_startf('timer_helm')

       do ie=nets,nete
          ieptr=ie-nets+1
          metinv => elem(ie)%metinv
          metdet => elem(ie)%metdet
          rmetdet  => elem(ie)%rmetdetp
          dx=0.5D0*elem(ie)%dx/lenscale    ! weak derivative element x dimension
          dy=0.5D0*elem(ie)%dy/lenscale    ! weak derivative element y dimension

          do k=1,nlev
             if (.not.cg%converged(k)) then

                ! ========================================
                ! Apply preconditioner: wrk2 = (M^-1) * wrk1 
                ! ========================================

                if (precon_method == "block_jacobi") then

                   if (blkjac_storage == "LUfactor") then
                      stop 'dgesl needs linpack'
!                      call dgesl(cg%state(ieptr)%r(:,k),cg%state(ieptr)%z(:,k),blkjac(ie)%E(:,:,k),blkjac(ie)%ipvt(:,k),nvsq)
                   else if (blkjac_storage == "inverse") then
                      call matvec(cg%state(ieptr)%r(:,k),cg%state(ieptr)%z(:,k),blkjac(ie)%E(:,:,k),nvsq)
                   end if

                   !                   iptr=1
                   !                   do j=1,nv
                   !                      do i=1,nv
                   !                         cg%wrk2(iptr,k,ieptr)=z(iptr)
                   !                         p(i,j) = cg%state(ieptr)%z(iptr,k)
                   !                         iptr=iptr+1
                   !                      end do
                   !                   end do

                else if (precon_method == "identity") then

                   iptr=1
                   do j=1,nv
                      do i=1,nv
                         cg%state(ieptr)%z(iptr,k) = cg%state(ieptr)%r(iptr,k)*rmetdet(i,j)
                         iptr=iptr+1
                      end do
                   end do


                end if

                !JMD===========================================
                !JMD   2*nv*nv*(nv + nv) Flops 
   		!JMD  SR = (4*nv*nv + 2*nv*nv + nv*nv)*Ld
   		!JMD  SUM(WS) = (6*nv*nv + 2*nv*nv + nv*nv
                !JMD===========================================

#ifdef _WK_GRAD
                gradp(:,:,:,k,ie)=gradient_wk(reshape(cg%state(ieptr)%z(:,k),(/np,np/)),deriv,dx,dy)
#else
                gradp(:,:,:,k,ie)=gradient(cg%state(ieptr)%z(:,k),deriv,dx,dy)
#endif


                ! =======================================
                ! rotate gradient to form contravariant
                !JMD  4*nv*nv Flops
                ! =======================================

                do j=1,nv
                   do i=1,nv
                      gradp1       = gradp(i,j,1,k,ie)
                      gradp2       = gradp(i,j,2,k,ie)
#if 1
                      gradp(i,j,1,k,ie) = metdet(i,j)*(metinv(1,1,i,j)*gradp1 + &
                           metinv(1,2,i,j)*gradp2)
                      gradp(i,j,2,k,ie) = metdet(i,j)*(metinv(2,1,i,j)*gradp1 + &
                           metinv(2,2,i,j)*gradp2)
#else
                      gradp(i,j,1,k,ie) = (metinv(1,1,i,j)*gradp1 + metinv(1,2,i,j)*gradp2)
                      gradp(i,j,2,k,ie) = (metinv(2,1,i,j)*gradp1 + metinv(2,2,i,j)*gradp2)
#endif
                   end do
                end do
             end if
          end do

          kptr=0
          call edgeVpack(edge2,gradp(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

          kptr=0
          call edgerotate(edge2,2*nlev,kptr,elem(ie)%desc)

       end do

       while_iter = while_iter + 1 

       call bndry_exchangeV(cg%hybrid,edge2)

       do ie=nets,nete
          ieptr=ie-nets+1

          rdx=2.0/(elem(ie)%dx*lenscale)
          rdy=2.0/(elem(ie)%dy*lenscale)

          rmv     => elem(ie)%rmv
          metdet  => elem(ie)%metdet
          metdetp => elem(ie)%metdetp
          mv      => elem(ie)%mv

          kptr=0
          call edgeVunpack(edge2, gradp(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
          do k=1,nlev
             if (.not.cg%converged(k)) then

                ! ====================
                ! 2*nv*nv Flops
                ! ====================

                do j=1,nv
                   do i=1,nv
                      gradp(i,j,1,k,ie) = rmv(i,j)*gradp(i,j,1,k,ie)
                      gradp(i,j,2,k,ie) = rmv(i,j)*gradp(i,j,2,k,ie)
                      !                     gradp(i,j,1,k,ie) = metdet(i,j)*rmv(i,j)*gradp(i,j,1,k,ie)
                      !                     gradp(i,j,2,k,ie) = metdet(i,j)*rmv(i,j)*gradp(i,j,2,k,ie)
                   end do
                end do

                ! ================================================
                ! Compute  Pseudo Laplacian(p), store in div
                !JMD   2*nv*np*(nv + np) Flops 
                ! ================================================

                div(:,:,k,ie) = divergence(gradp(:,:,:,k,ie),deriv,rdx,rdy)

                do j=1,nv
                   do i=1,nv
                      div(i,j,k,ie) = mv(i,j)*div(i,j,k,ie)
                   end do
                end do
             end if
          end do

          k=0
          call edgeVpack(edge1, div(1,1,1,ie), nlev, kptr, elem(ie)%desc)

       end do


       call bndry_exchangeV(cg%hybrid,edge1)


       ! ==========================================
       ! compute Helmholtz operator, store in wrk3
       !  4*np*np Flops
       ! ==========================================

       do ie=nets,nete
          ieptr=ie-nets+1

          rmv      => elem(ie)%rmv
          rmetdet  => elem(ie)%rmetdetp
          metdet   => elem(ie)%metdet
          metdetp  => elem(ie)%metdetp

          kptr=0
          call edgeVunpack(edge1, div(1,1,1,ie), nlev, kptr, elem(ie)%desc)

          do k=1,nlev
             if (.not.cg%converged(k)) then

                iptr=1
                do j=1,nv
                   do i=1,nv
                      cg%state(ieptr)%s(iptr,k) = metdetp(i,j)*cg%state(ieptr)%z(iptr,k)+lambdasq(k)*rmv(i,j)*div(i,j,k,ie)
                      iptr=iptr+1
                   end do
                end do

             end if
          end do
       end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
       call t_stopf('timer_helm')

    end do  ! CG solver while loop
!    print *,'CG inter = ',cg%iter
    ! ===============================
    ! Converged! unpack wrk3 into x
    ! ===============================

    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,nv
             do i=1,nv
                !                x(i,j,k,ie)=cg%wrk3(iptr,k,ieptr)
                x(i,j,k,ie)=cg%state(ieptr)%x(iptr,k)
                iptr=iptr+1
             end do
          end do
       end do
    end do
    call t_stopf('pcg_solver_nonstag')
#endif
  end function pcg_solver_nonstag

  subroutine blkjac_init(elem, deriv,lambdasq,nets,nete,blkjac)
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t, gradient_wk, gradient, divergence, derivative_stag_t
    use physical_constants, only : rearth
    use dimensions_mod, only : nlev, nv, np
    use parallel_mod, only : haltmp
    type(element_t), intent(in), target :: elem(:)
    type (derivative_t)                  :: deriv
    real (kind=real_kind), intent(in)    :: lambdasq(nlev)
    integer(kind=int_kind), intent(in) :: nets
    integer(kind=int_kind), intent(in) :: nete
    type (blkjac_t)      , intent(inout) :: blkjac(nets:nete)
#if 0
    !JMD    real (kind=real_kind)             :: E(npsq,npsq,nlev,nets:nete)
    !JMD    integer                           :: ipvt(npsq,nlev,nets:nete)
#endif

    ! ===============
    ! Local Variables
    ! ===============

    integer :: ie,kk
    integer :: i,j,k
    integer :: iptr
    integer :: ieptr
    integer :: lwork
    integer :: info

    real (kind=real_kind) :: p(npsq)
    real (kind=real_kind) :: z(npsq)
    real (kind=real_kind) :: gradp(nv,nv,2) ! velocity buffer
    real (kind=real_kind) :: div(np,np)

    real (kind=real_kind), pointer :: metdetp(:,:)
    real (kind=real_kind), pointer :: mp(:, :)

    real (kind=real_kind), pointer :: rmv(:,:)
    real (kind=real_kind), pointer :: mv(:,:)

    real (kind=real_kind), pointer :: metdet(:,:)
    real (kind=real_kind), pointer :: rmetdet(:,:)
    real (kind=real_kind), pointer :: metinv(:,:,:,:)

    real (kind=real_kind) :: rdx,rdy
    real (kind=real_kind) ::  dx, dy

    real (kind=real_kind) :: gradp1
    real (kind=real_kind) :: gradp2
    real (kind=real_kind) :: det(2)

    real (kind=real_kind) :: lenscale,lsq

    ! =================================
    ! Begin executable code
    ! =================================    
#ifndef CAM
    !DBG print *,'blkjac_init: point #0 nets,nete',nets,nete
    lenscale=rearth
    !DBG print *,'blkjac_init: point #0.1 nets,nete',nets,nete

    lsq=lenscale*lenscale

    ! ===========================
    ! zero out E matrix
    ! ===========================
    !DBG print *,'blkjac_init: point #0.2 nets,nete',nets,nete
    do ie=nets,nete
       blkjac(ie)%E = 0.0d0
       blkjac(ie)%ipvt = 0
    enddo
    !DBG print *,'blkjac_init: point #0.3 nets,nete',nets,nete

    do k=1,nlev

       !JMDprint *,'blkjac_init: point #1'
       do ie=nets,nete

          !JMD print *,'blkjac_init: point #2 ie,nets,nete',ie,nets,nete
          !JMD print *,'blkjac_init: point #3 ie,nets,nete',ie,nets,nete


          metdet  => elem(ie)%metdet
          metdetp => elem(ie)%metdetp
          metinv  => elem(ie)%metinv
          rmv     => elem(ie)%rmv
          mv     => elem(ie)%mv
          mp      => elem(ie)%mp          

          dx=0.5D0*elem(ie)%dx/lenscale ! weak derivative element x dimension
          dy=0.5D0*elem(ie)%dy/lenscale    ! weak derivative element y dimension
          rdx=2.0D0/(elem(ie)%dx*lenscale) ! strong derivative inverse x length 
          rdy=2.0D0/(elem(ie)%dy*lenscale) ! strong derivative inverse y length 

          do kk = 1, npsq     ! delta fn excitation index

             p(:) = 0.0d0
             p(kk)= 1.0d0 

             ! =========================================================
             ! Compute (weak form) pressure gradient(s) on velocity grid
             ! =========================================================

             gradp(:,:,:)= gradient_wk(reshape(p,(/np,np/)),deriv,dx,dy)

             do j=1,nv
                do i=1,nv
                   gradp1       = gradp(i,j,1)
                   gradp2       = gradp(i,j,2)
                   gradp(i,j,1) = metdet(i,j)*(metinv(1,1,i,j)*gradp1 + metinv(1,2,i,j)*gradp2)
                   gradp(i,j,2) = metdet(i,j)*(metinv(2,1,i,j)*gradp1 + metinv(2,2,i,j)*gradp2)
                end do
             end do

             ! ================================================
             ! Apply inverse velocity mass matrix Mu^{-1}
             ! ================================================

             do j=2,nv-1
                do i=2,nv-1
                   gradp(i,j,1) = gradp(i,j,1)*rmv(i,j)
                   gradp(i,j,2) = gradp(i,j,2)*rmv(i,j)
                end do
             end do

             ! ================================================
             ! Zero (Neumann) pressure gradients along the boundaries
             !
             ! zero north/south edge
             ! ================================================

             do j=1,nv
                gradp(1 ,j,1)=0.0D0
                gradp(1 ,j,2)=0.0D0
                gradp(nv,j,1)=0.0D0
                gradp(nv,j,2)=0.0D0
             end do

             ! ================================================
             ! zero east/west edge
             ! ================================================

             do i=1,nv
                gradp(i,1 ,1)=0.0D0
                gradp(i,1 ,2)=0.0D0
                gradp(i,nv,1)=0.0D0
                gradp(i,nv,2)=0.0D0
             end do

             ! ==================================================
             ! Compute divergence on pressure grid
             ! ==================================================

             div = divergence(gradp,deriv,rdx,rdy)

             ! ==================================================
             ! Compute Rt = Mp.(pp + dt^2.pmean.D.Mu^-1.D'.pp)
             ! ==================================================

             iptr=1
             do j=1,np
                do i=1,np
                   blkjac(ie)%E(iptr,kk,k) = metdet(i,j)*p(iptr) + lambdasq(k) * div(i,j)/(lsq*mv(i,j))
                   iptr=iptr+1
                end do
             end do

             !JMD print *,'blkjac_init: point #10'
          end do  ! end loop over kk (delta function exitation location)
          !DBG print *,'blkjac_init: point #11 ie,nets,nete',ie,nets,nete

          ! ======================================
          ! Lu factorize Block Jacobi matrix...
          ! ======================================
          !             print *, __FILE__,__LINE__,minval( E(:,:,k,ie)),maxval( E(:,:,k,ie))

          !          call dgefa(blkjac(ie)%E(1,1,k),npsq,npsq,blkjac(ie)%ipvt(1,k), info) !LINPACK
          call dgetrf(npsq,npsq,blkjac(ie)%E(1,1,k),npsq,blkjac(ie)%ipvt(1,k), info) !LAPACK

          !DBG print *,'blkjac_init: point #12 ie,nets,nete',ie,nets,nete

          if (info /= 0) then
             print *,"WARNING: singular Block Jacobi matrix detected during preconditioning setup"
          end if

          if (blkjac_storage == "inverse") then

             if (k==1 .and. ie == 1) then
                !debug                 print *,"storing Block Jacobi Inverse"
             end if

             ! =======================================
             ! invert E matrix on each processor
             ! =======================================


             !             call dgedi( blkjac(ie)%E(:,:,k), npsq, npsq, blkjac(ie)%ipvt(:,k), det, z(:), 1 ) !LINPACK
             lwork = npsq
             call dgetri(npsq, blkjac(ie)%E(1,1,k), npsq, blkjac(ie)%ipvt(1,k),z,lwork,info) !LAPACK


          else if (blkjac_storage == "LUfactor") then

             do i=1,npsq
                blkjac(ie)%E(i,i,k)=1.0D0/blkjac(ie)%E(i,i,k)
             end do

             if (k==1 .and. ie == 1) then
                print *,"LU Block Jacobi storage"
             end if

          else
             call haltmp("bad Block Jacobi storage option")
          end if

          !DBG print *,'blkjac_init: point #15'
       end do  ! end element loop
       !DBG print *,'blkjac_init: point #16'

    end do ! end level loop
    !DBG print *,'blkjac_init: point #17'
#endif
  end subroutine blkjac_init

end module solver_mod
