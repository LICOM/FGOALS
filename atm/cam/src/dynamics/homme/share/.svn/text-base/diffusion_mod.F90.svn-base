#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module diffusion_mod
  ! =======================
  use kinds,              only : real_kind
  ! =======================
  use dimensions_mod,     only : nlev, nv, qsize
  ! =======================
  use derivative_mod,     only : gradient, vorticity, derivative_t, divergence
  ! =======================
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, initedgebuffer
  ! =======================
  private
  public :: diffusion_init, prim_diffusion, scalar_diffusion
  type(EdgeBuffer_t)  :: edgeS1, edgeS2
  type (EdgeBuffer_t) :: edge3
  type (EdgeBuffer_t) :: edge4

contains
  subroutine diffusion_init()

       call initEdgeBuffer(edgeS1,qsize*nlev)
       call initEdgeBuffer(edgeS2,2*qsize*nlev)
       call initEdgeBuffer(edge3, 3*nlev)
       call initEdgeBuffer(edge4, 4*nlev)

  end subroutine diffusion_init

  subroutine prim_diffusion(elem, nets,nete,np1,deriv,dt2,hybrid)
    use hybrid_mod, only : hybrid_t
    use physical_constants, only : rearth
    use element_mod, only : element_t
    use bndry_mod, only : bndry_exchangev
    use control_mod, only : nu
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    implicit none
    type(element_t), intent(inout), target :: elem(:)
    integer, intent(in)                :: nets,nete
    integer, intent(in)                :: np1
    type(derivative_t), intent(in)     :: deriv
    real(kind=real_kind),intent(in)    :: dt2
    type(hybrid_t),intent(in)       :: hybrid

    ! ======================
    ! Local Variables
    ! ======================

    real(kind=real_kind) :: temp(nv,nv)

    real(kind=real_kind) :: grad_T_np1(nv,nv,2,nlev,nets:nete)
    real(kind=real_kind) :: zeta_np1(nv,nv,nlev,nets:nete)
    real(kind=real_kind) :: div_np1(nv,nv,nlev,nets:nete)
    real(kind=real_kind) :: lap_v_np1(nv,nv,2,nlev,nets:nete)
    real(kind=real_kind) :: lap_T_np1(nv,nv,nlev,nets:nete)

    real(kind=real_kind) :: grad_zeta(nv,nv,2)
    real(kind=real_kind) :: grad_div(nv,nv,2)
    real(kind=real_kind) :: grad_tmp(nv,nv,2)
    real(kind=real_kind) :: vco(nv,nv,2)
    real(kind=real_kind) :: gv(nv,nv,2)
    real(kind=real_kind) :: zeta_tmp(nv,nv)
    real(kind=real_kind) :: div_tmp(nv,nv)

    real(kind=real_kind) :: v1,v2
    real(kind=real_kind) :: rdx,rdy

    real (kind=real_kind), dimension(nv,nv)      :: rmetdetv   

    real(kind=real_kind), dimension(:,:), pointer :: mv
    real(kind=real_kind), dimension(:,:), pointer :: rmv
    real(kind=real_kind), dimension(:,:), pointer :: metdet
    real(kind=real_kind), dimension(:,:,:,:), pointer :: met
    real(kind=real_kind), dimension(:,:,:,:), pointer :: metinv
    real(kind=real_kind), dimension(:,:,:,:), pointer :: D
    real(kind=real_kind), dimension(:,:,:,:), pointer :: Dinv

    integer :: i,j,k,ie,kptr

    call t_startf('prim_diffusion')

    do ie=nets,nete

       rdx=2.0_real_kind/(elem(ie)%dx*rearth) ! strong derivative inverse x length
       rdy=2.0_real_kind/(elem(ie)%dy*rearth) ! strong derivative inverse y length

       mv       => elem(ie)%mv
       rmv      => elem(ie)%rmv
       met      => elem(ie)%met
       metinv   => elem(ie)%metinv
       metdet   => elem(ie)%metdet
       rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)

       do k=1,nlev   

          grad_tmp(:,:,:) = gradient(elem(ie)%state%T(:,:,k,np1),deriv,rdx,rdy)

          do j=1,nv
             do i=1,nv
                 grad_T_np1(i,j,1,k,ie) = mv(i,j)*(metinv(1,1,i,j)*grad_tmp(i,j,1)+metinv(1,2,i,j)*grad_tmp(i,j,2))
                 grad_T_np1(i,j,2,k,ie) = mv(i,j)*(metinv(2,1,i,j)*grad_tmp(i,j,1)+metinv(2,2,i,j)*grad_tmp(i,j,2))
                v1     = elem(ie)%state%v(i,j,1,k,np1)
                v2     = elem(ie)%state%v(i,j,2,k,np1)

                gv(i,j,1) = metdet(i,j)*(elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2)
                gv(i,j,2) = metdet(i,j)*(elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2)
                vco(i,j,1) = elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(2,1,i,j)*v2
                vco(i,j,2) = elem(ie)%D(1,2,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2
             end do
          end do

          div_tmp(:,:)  = divergence(gv,deriv,rdx,rdy)
          zeta_tmp(:,:) = vorticity(vco,deriv,rdx,rdy)

          do j=1,nv
             do i=1,nv
                div_np1(i,j,k,ie)  = mv(i,j)*rmetdetv(i,j)*div_tmp(i,j)
                zeta_np1(i,j,k,ie) = mv(i,j)*rmetdetv(i,j)*zeta_tmp(i,j)
             end do
          end do

       end do

       kptr=0
       call edgeVpack(edge4,grad_T_np1(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

       kptr=2*nlev
       call edgeVpack(edge4,zeta_np1(1,1,1,ie),nlev,kptr,elem(ie)%desc)

       kptr=3*nlev
       call edgeVpack(edge4,div_np1(1,1,1,ie),nlev,kptr,elem(ie)%desc)

       kptr=0
       call edgerotate(edge4,2*nlev,kptr,elem(ie)%desc)

    end do

    call bndry_exchangeV(hybrid,edge4)

    do ie=nets,nete

       rdx=2.0_real_kind/(elem(ie)%dx*rearth) ! strong derivative inverse x length
       rdy=2.0_real_kind/(elem(ie)%dy*rearth) ! strong derivative inverse y length

       mv       => elem(ie)%mv
       rmv      => elem(ie)%rmv
       metdet   => elem(ie)%metdet
       metinv   => elem(ie)%metinv
       rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)
       D        => elem(ie)%D
       Dinv     => elem(ie)%Dinv

       kptr=0
       call edgeVunpack(edge4, grad_T_np1(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

       kptr=2*nlev
       call edgeVunpack(edge4, zeta_np1(1,1,1,ie), nlev, kptr, elem(ie)%desc)

       kptr=3*nlev
       call edgeVunpack(edge4, div_np1(1,1,1,ie), nlev, kptr, elem(ie)%desc)
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
       ! ======================================
       ! compute Laplacian of T(n+1)
       ! ======================================

       do k=1,nlev
          do j=1,nv
             do i=1,nv

                ! ======================================================
                ! complete global assembly of gradT(n+1), zeta(n+1), div(n+1)
                ! ======================================================

!                v1 = rmv(i,j)*grad_T_np1(i,j,1,k,ie)
!                v2 = rmv(i,j)*grad_T_np1(i,j,2,k,ie)

                zeta_np1(i,j,k,ie) = rmv(i,j)*zeta_np1(i,j,k,ie)
                div_np1(i,j,k,ie)  = rmv(i,j)*div_np1(i,j,k,ie)
                elem(ie)%derived%zeta(i,j,k) = zeta_np1(i,j,k,ie)

                ! ==========================================
                ! Compute contravariant gradient(T(n+1))
                ! ==========================================

!                grad_T_np1(i,j,1,k,ie) = metdet(i,j)*(metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2)
!                grad_T_np1(i,j,2,k,ie) = metdet(i,j)*(metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2)
                grad_T_np1(i,j,1,k,ie) = metdet(i,j)*rmv(i,j)*grad_T_np1(i,j,1,k,ie)
                grad_T_np1(i,j,2,k,ie) = metdet(i,j)*rmv(i,j)*grad_T_np1(i,j,2,k,ie)

             end do
          end do

          div_tmp(:,:) = divergence(grad_T_np1(:,:,:,k,ie),deriv,rdx,rdy)

          do j=1,nv
             do i=1,nv
                lap_T_np1(i,j,k,ie)   = mv(i,j)*rmetdetv(i,j)*div_tmp(i,j)
             end do
          end do

          ! ==========================================
          ! compute covariant Laplacian of v(n+1)  
          ! ==========================================

          grad_div(:,:,:)  = gradient(div_np1(:,:,k,ie),deriv,rdx,rdy)
          grad_zeta(:,:,:) = gradient(zeta_np1(:,:,k,ie),deriv,rdx,rdy)

          do j=1,nv
             do i=1,nv

                v2 = Dinv(1,1,i,j)*grad_zeta(i,j,1) + Dinv(2,1,i,j)*grad_zeta(i,j,2)
                v1 = Dinv(1,2,i,j)*grad_zeta(i,j,1) + Dinv(2,2,i,j)*grad_zeta(i,j,2)

                v2 = -v2
                grad_zeta(i,j,1) = v1
                grad_zeta(i,j,2) = v2
                v1 = grad_div(i,j,1)
                v2 = grad_div(i,j,2)

                grad_div(i,j,1) = Dinv(1,1,i,j)*v1 + Dinv(2,1,i,j)*v2
                grad_div(i,j,2) = Dinv(1,2,i,j)*v1 + Dinv(2,2,i,j)*v2
                lap_v_np1(i,j,1,k,ie) = mv(i,j)*(grad_div(i,j,1) - grad_zeta(i,j,1))
                lap_v_np1(i,j,2,k,ie) = mv(i,j)*(grad_div(i,j,2) - grad_zeta(i,j,2))

             end do
          end do

       end do

       kptr=0
       call edgeVpack(edge3, lap_v_np1(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

       kptr=2*nlev
       call edgeVpack(edge3, lap_T_np1(1,1,1,ie),nlev,kptr,elem(ie)%desc)

    end do

    call bndry_exchangeV(hybrid,edge3)

    do ie=nets,nete

       rmv   => elem(ie)%rmv
       metdet   => elem(ie)%metdet
       rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)

       kptr=0
       call edgeVunpack(edge3, lap_v_np1(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

       kptr=2*nlev
       call edgeVunpack(edge3, lap_T_np1(1,1,1,ie), nlev, kptr, elem(ie)%desc)

       do k=1,nlev   

          do j=1,nv
             do i=1,nv
                lap_T_np1(i,j,k,ie)   = rmv(i,j)*lap_T_np1(i,j,k,ie)
                lap_v_np1(i,j,1,k,ie) = rmv(i,j)*lap_v_np1(i,j,1,k,ie)
                lap_v_np1(i,j,2,k,ie) = rmv(i,j)*lap_v_np1(i,j,2,k,ie)

                elem(ie)%state%T(i,j,k,np1)   = elem(ie)%state%T(i,j,k,np1)   + nu*dt2*lap_T_np1(i,j,k,ie)
                elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%state%v(i,j,1,k,np1) + nu*dt2*lap_v_np1(i,j,1,k,ie)
                elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%state%v(i,j,2,k,np1) + nu*dt2*lap_v_np1(i,j,2,k,ie)
             end do
          end do

       end do

    end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
    call t_stopf('prim_diffusion')

  end subroutine prim_diffusion

  subroutine scalar_diffusion(elem, nets,nete,np1,deriv,dt2,hybrid)
    use hybrid_mod, only : hybrid_t
    use physical_constants, only : rearth
    use element_mod, only : element_t
    use bndry_mod, only : bndry_exchangev
    use control_mod, only : nu_s
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    implicit none
    type(element_t), intent(inout), target :: elem(:)
    integer, intent(in)                :: nets,nete
    integer, intent(in)                :: np1
    type(derivative_t), intent(in)     :: deriv
    real(kind=real_kind),intent(in)    :: dt2
    type(hybrid_t),intent(in)       :: hybrid

    ! ======================
    ! Local Variables
    ! ======================
    real(kind=real_kind) :: temp(nv,nv)

    real(kind=real_kind) :: grad_Q_np1(nv,nv,2,nlev,qsize,nets:nete)
    real(kind=real_kind) :: lap_Q_np1(nv,nv,nlev,qsize,nets:nete)

    real(kind=real_kind) :: grad_tmp(nv,nv,2)
    real(kind=real_kind) :: div_tmp(nv,nv)

    real(kind=real_kind) :: v1,v2
    real(kind=real_kind) :: rdx,rdy

    real (kind=real_kind), dimension(nv,nv)      :: rmetdetv   

    real(kind=real_kind), dimension(:,:), pointer :: mv
    real(kind=real_kind), dimension(:,:), pointer :: rmv
    real(kind=real_kind), dimension(:,:), pointer :: metdet
    real(kind=real_kind), dimension(:,:,:,:), pointer :: met
    real(kind=real_kind), dimension(:,:,:,:), pointer :: metinv
    real(kind=real_kind), dimension(:,:,:,:), pointer :: D
    real(kind=real_kind), dimension(:,:,:,:), pointer :: Dinv

    integer :: i,j,k,ie, q

    call t_startf('scalar_diffusion')
    do ie=nets,nete
       do q=1,qsize

          rdx=2.0_real_kind/(elem(ie)%dx*rearth) ! strong derivative inverse x length
          rdy=2.0_real_kind/(elem(ie)%dy*rearth) ! strong derivative inverse y length

          mv       => elem(ie)%mv
          metinv   => elem(ie)%metinv
          metdet   => elem(ie)%metdet

          do k=1,nlev   
             
             grad_tmp(:,:,:) = gradient(elem(ie)%state%Q(:,:,k,q,np1),deriv,rdx,rdy)

             do j=1,nv
                do i=1,nv
                   !                grad_Q_np1(i,j,1,k,ie) = mv(i,j)*grad_tmp(i,j,1)
                   !                grad_Q_np1(i,j,2,k,ie) = mv(i,j)*grad_tmp(i,j,2)

                   v1 = mv(i,j)*grad_tmp(i,j,1)
                   v2 = mv(i,j)*grad_tmp(i,j,2)

                   grad_Q_np1(i,j,1,k,q,ie) = metdet(i,j)*(metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2)
                   grad_Q_np1(i,j,2,k,q,ie) = metdet(i,j)*(metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2)

                end do
             end do

          end do

       end do
       
       call edgeVpack(edgeS2,grad_Q_np1(:,:,:,:,:,ie),2*nlev*qsize,0,elem(ie)%desc)

       call edgerotate(edgeS2,2*nlev*qsize,0,elem(ie)%desc)

    end do

    call bndry_exchangeV(hybrid,edgeS2)

    do ie=nets,nete
       rdx=2.0_real_kind/(elem(ie)%dx*rearth) ! strong derivative inverse x length
       rdy=2.0_real_kind/(elem(ie)%dy*rearth) ! strong derivative inverse y length

       mv       => elem(ie)%mv
       rmv      => elem(ie)%rmv
       metdet   => elem(ie)%metdet
       metinv   => elem(ie)%metinv
       rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)
       D        => elem(ie)%D
       Dinv     => elem(ie)%Dinv
       
       call edgeVunpack(edgeS2, grad_Q_np1(:,:,:,:,:,ie), 2*nlev*qsize, 0, elem(ie)%desc)
#ifdef DEBUGOMP
!$OMP BARRIER
#endif

       do q=1,qsize

          do k=1,nlev
             do j=1,nv
                do i=1,nv

                   !                v1 = rmv(i,j)*grad_Q_np1(i,j,1,k,ie)
                   !                v2 = rmv(i,j)*grad_Q_np1(i,j,2,k,ie)

                   !                grad_Q_np1(i,j,1,k,ie) = metdet(i,j)*(metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2)
                   !                grad_Q_np1(i,j,2,k,ie) = metdet(i,j)*(metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2)

                   grad_Q_np1(i,j,1,k,q,ie) = rmv(i,j)*grad_Q_np1(i,j,1,k,q,ie)
                   grad_Q_np1(i,j,2,k,q,ie) = rmv(i,j)*grad_Q_np1(i,j,2,k,q,ie)

                end do
             end do

             div_tmp(:,:) = divergence(grad_Q_np1(:,:,:,k,q,ie),deriv,rdx,rdy)

             do j=1,nv
                do i=1,nv
                   lap_Q_np1(i,j,k,q,ie)   = mv(i,j)*rmetdetv(i,j)*div_tmp(i,j)
                end do
             end do

          end do

       end do

       call edgeVpack(edgeS1, lap_Q_np1(1,1,1,1,ie),nlev*qsize, 0,elem(ie)%desc)

    end do
    call bndry_exchangeV(hybrid,edgeS1)
    do ie=nets,nete
          
       rmv   => elem(ie)%rmv
       metdet   => elem(ie)%metdet
       rmetdetv(:,:)=1.0_real_kind/elem(ie)%metdet(:,:)

       
       call edgeVunpack(edgeS1, lap_Q_np1(1,1,1,1,ie), nlev*qsize, 0, elem(ie)%desc)

       do q=1,qsize
          do k=1,nlev   

             do j=1,nv
                do i=1,nv
                   lap_Q_np1(i,j,k,q,ie)   = rmv(i,j)*lap_Q_np1(i,j,k,q,ie)
                   elem(ie)%state%Q(i,j,k,q,np1)   = elem(ie)%state%Q(i,j,k,q,np1)   + nu_s*dt2*lap_Q_np1(i,j,k,q,ie)
                end do
             end do

          end do

       end do
    end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
    call t_stopf('scalar_diffusion')
  end subroutine scalar_diffusion

end module diffusion_mod
