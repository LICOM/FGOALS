#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module global_norms_mod
  use kinds, only : iulog
  implicit none
  private

  public :: l1_snorm
  public :: l2_snorm
  public :: linf_snorm

  public :: l1_vnorm
  public :: l2_vnorm
  public :: linf_vnorm

  public :: test_global_integral
  public :: global_integral

  private :: global_maximum

contains

  ! ================================
  ! global_integral:
  !
  ! eq 81 in Williamson, et. al. p 218
  ! for spectral elements
  !
  ! ================================
  ! --------------------------
  function global_integral(elem, h,hybrid,npts,nets,nete) result(I_sphere)
    use kinds,       only : real_kind
    use hybrid_mod,  only : hybrid_t
    use element_mod, only : element_t
    use dimensions_mod, only : nv, np, nelemd
    use physical_constants, only : dd_pi
#ifdef CAM
    use repro_sum_mod, only: repro_sum ! _EXTERNAL
#else
    use reduction_mod, only : red_sum, psum_mt
#endif
    use parallel_mod, only: global_shared_buf, global_shared_sum

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind) :: I_sphere

    real (kind=real_kind) :: I_priv
    real (kind=real_kind) :: I_shared
    common /gblintcom/I_shared

    ! Local variables

    integer :: ie,j,i
    real(kind=real_kind) :: I_tmp(1)

    real (kind=real_kind) :: da
#ifdef CAM
    real (kind=real_kind) :: J_tmp(nets:nete)
!
! This algorythm is independent of thread count and task count.
! This is a requirement of consistancy checking in cam.
!
    J_tmp = 0.0D0
    if(npts==np) then
       do ie=nets,nete
          do j=1,np
             do i=1,np
                da = elem(ie)%mp(i,j)*elem(ie)%metdetp(i,j)
                J_tmp(ie) = J_tmp(ie) + da*h(i,j,ie)
             end do
          end do
       end do       
    else if(npts==nv) then
       do ie=nets,nete
          do j=1,nv
             do i=1,nv
                da = elem(ie)%mv(i,j)*elem(ie)%metdet(i,j)
                J_tmp(ie) = J_tmp(ie) + da*h(i,j,ie)
             end do
          end do
       end do       
    end if
    do ie=nets,nete
      global_shared_buf(ie,1) = J_tmp(ie)
    enddo
!$OMP BARRIER
!$OMP MASTER
    call repro_sum(global_shared_buf, global_shared_sum, nelemd, nelemd, 1, commid=hybrid%par%comm)
!$OMP END MASTER
!$OMP BARRIER
    I_tmp = global_shared_sum(1)
    I_sphere = I_tmp(1)/(4.0D0*DD_PI)

#else
    I_shared = 0.0D0
    I_priv    = 0.0D0
    if (npts==np) then
       do ie=nets,nete
          do j=1,np
             do i=1,np
                da = elem(ie)%mp(i,j)*elem(ie)%metdetp(i,j)
                I_priv = I_priv + da*h(i,j,ie)
             end do
          end do
       end do
    else if (npts==nv) then
       do ie=nets,nete
          do j=1,nv
             do i=1,nv
                da = elem(ie)%mv(i,j)*elem(ie)%metdet(i,j)
                I_priv = I_priv + da*h(i,j,ie)
             end do
          end do
       end do
    end if

    I_tmp(1)=I_priv
    call psum_mt(red_sum,I_tmp,1,hybrid)
    !   I_sphere = red_sum%buf(1,1)/(4.0D0*DD_PI)
    I_sphere = red_sum%buf(1)/(4.0D0*DD_PI)
#endif

  end function global_integral

  ! ================================
  ! test_global_integral:
  !
  ! test that the global integral of 
  ! the area of the sphere is 1.
  !
  ! ================================

  subroutine test_global_integral(elem,hybrid,nets,nete,mindxout)
    use kinds,       only : real_kind
    use hybrid_mod,  only : hybrid_t
    use element_mod, only : element_t
    use dimensions_mod, only : np,nv,ne

    use reduction_mod, only : ParallelMin,ParallelMax,ParallelSum
    use quadrature_mod, only : gausslobatto, quadrature_t
    use physical_constants, only : rearth,dd_pi
    use control_mod, only : nu, hypervis_order, nu_top

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: nets,nete
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind),intent(out), optional :: mindxout

    real (kind=real_kind)             :: I_sphere
    ! Local variables

    real (kind=real_kind), allocatable :: h(:,:,:)
    real (kind=real_kind) :: min_area,max_area,area,min_len,mindx
    real (kind=real_kind) :: eq_len,eq_min,eq_ave,eq_max
    integer :: ie,T_equiv
    type (quadrature_t)    :: gv

    allocate(h(np,np,nets:nete))

    h(:,:,nets:nete)=1.0D0

    I_sphere = global_integral(elem, h(:,:,nets:nete),hybrid,np,nets,nete)

    min_area=1d99
    max_area=0
    do ie=nets,nete
       area=sum(elem(ie)%spheremv(:,:))
       min_area=min(min_area,area)
       max_area=max(max_area,area)
    enddo
    min_area=ParallelMin(min_area,hybrid)
    max_area=ParallelMax(max_area,hybrid)

    ! tot_area will be 4*pi (area of sphere radius 1). Convert to meters:
    min_len=sqrt(min_area)*rearth

    gv=gausslobatto(nv)
    mindx = .5*abs(gv%points(1)-gv%points(2)) * min_len  ! smallest grid spacing in domain

    eq_len=(2*dd_pi*rearth/(4*ne))/1000               ! average element length at equator, in km
    eq_min = .5*abs(gv%points(1)-gv%points(2)) * eq_len  ! min grid spacing at equator
    eq_ave = eq_len/(nv-1)                               ! average grid spacing at equator
    eq_max = .5*abs(gv%points(nv/2)-gv%points(nv/2 +1)) * eq_len  ! max grid spacing at equator
    T_equiv = (4*ne*(nv-1))/3d0                ! equivelent resolution of 2/3 dealiased spherical harmonics

    ! for an equation du/dt = i c u, leapfrog is stable for |c u dt| < 1
    ! Consider a gravity wave at the equator, c=340m/s  
    ! u = exp(i kmax x/ a ) with x = longitude,  and kmax =  pi a / dx, 
    ! u = exp(i pi x / dx ),   so du/dt = c du/dx becomes du/dt = i c pi/dx u
    ! stable for dt < dx/(c*pi)
    ! CAM 26 level AMIP simulation: max gravity wave speed 341.75 m/s
    if (hybrid%masterthread) then
       write(iulog,* )""
       write(iulog,* )"Running Global Integral Diagnostic..."
       write(iulog,*)"Area of unit sphere is",I_sphere
       write(iulog,*)"Should be 1.0 to round off..."
       write(iulog,'(a,f7.3)') 'Element area:  max/min',(max_area/min_area)
       write(iulog,'(a,3f8.2,a,i4,a)') 'Equatorial grid spacing: ave,min,max (km) = ',eq_ave,eq_min,eq_max,' (T',T_equiv,')'
       write(iulog,'(a,f10.2,a)') 'gravity wave dt:  (global min_dx/pi)    /342m/s: ',(mindx/dd_pi)/342,'s'
       write(iulog,'(a,f10.2,a)') 'gravity wave dt:  (equatorial min_dx/pi)/342m/s: ',(1000*eq_min/dd_pi)/342,'s'
       if (nu>0) then
          if (hypervis_order==1) write(iulog,'(a,f10.2,a)') 'viscosity dt: (min_dx/pi)**2/nu:   ',(mindx/dd_pi)**2/nu,'s'
          if (hypervis_order==2) then
             write(iulog,'(a,f10.2,a)') 'hyper viscosity dt: (min_dx/pi)**4/nu:   ',(mindx/dd_pi)**4/nu,'s'
             if (nv==4) then
                write(iulog,'(a,f10.2,a)') 'Observed hypervis dt: (Leapfrog, nv=4) ',1.25d23/(nu*ne**4.0),'s'
                write(iulog,'(a,f10.2,a)') 'Observed hypervis dt: (RK2, nv=4) ',2*1.25d23/(nu*ne**4.0),'s'
             endif
          endif
       endif
       if(nu_top>0) then
          write(iulog,'(a,f10.2,a)') 'TOP3 viscosity dt: (min_dx/pi)**2/nu_top: ',&
             (mindx/dd_pi)**2/nu_top,'s'
          write(iulog,'(a,f10.2,a)') 'nu_top(unscaled)dt: (Leapfrog, nv=4) ',sqrt(1.25d23)/(nu_top*ne**2.0),'s'
          write(iulog,'(a,f10.2,a)') 'nu_top(unscaled)dt: (RK2, nv=4) ',2*sqrt(1.25d23)/(nu_top*ne**2.0),'s'
       end if
    end if

    deallocate(h)
    
    if(present(mindxout))mindxout=mindx

  end subroutine test_global_integral



  ! ================================
  ! global_maximum:
  !
  ! Find global maximum on sphere
  !
  ! ================================

  function global_maximum(h,hybrid,npts,nets,nete) result(Max_sphere)
    use kinds, only : real_kind
    use hybrid_mod, only : hybrid_t
    use reduction_mod, only : red_max, pmax_mt

    integer              , intent(in) :: npts,nets,nete     
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind) :: Max_sphere

    ! Local variables

    real (kind=real_kind) :: redp(1)

    Max_sphere = MAXVAL(h(:,:,nets:nete))

    redp(1) = Max_sphere
    call pmax_mt(red_max,redp,1,hybrid)
    Max_sphere = red_max%buf(1)

  end function global_maximum
  ! ==========================================================
  ! l1_snorm:
  !
  ! computes the l1 norm per Williamson et al, p. 218 eq(8)
  ! for a scalar quantity
  ! ===========================================================

  function l1_snorm(elem, h,ht,hybrid,npts,nets,nete) result(l1)
    use kinds, only : real_kind
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in) :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: ht(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: l1     

    ! Local variables

    real (kind=real_kind) :: dhabs(npts,npts,nets:nete)
    real (kind=real_kind) :: htabs(npts,npts,nets:nete)
    real (kind=real_kind) :: dhabs_int
    real (kind=real_kind) :: htabs_int
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dhabs(i,j,ie) = ABS(h(i,j,ie)-ht(i,j,ie))
             htabs(i,j,ie) = ABS(ht(i,j,ie))
          end do
       end do
    end do

    dhabs_int = global_integral(elem, dhabs(:,:,nets:nete),hybrid,npts,nets,nete)
    htabs_int = global_integral(elem, htabs(:,:,nets:nete),hybrid,npts,nets,nete)

    l1 = dhabs_int/htabs_int

  end function l1_snorm

  ! ===========================================================
  ! l1_vnorm:
  !
  ! computes the l1 norm per Williamson et al, p. 218 eq(97),
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function l1_vnorm(elem, v,vt,hybrid,npts,nets,nete) result(l1)
    use kinds, only : real_kind
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in), target :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: l1     

    ! Local variables

    real (kind=real_kind), dimension(:,:,:,:), pointer :: met
    real (kind=real_kind) :: dvsq(npts,npts,nets:nete)
    real (kind=real_kind) :: vtsq(npts,npts,nets:nete)
    real (kind=real_kind) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: dv1,dv2
    real (kind=real_kind) :: vt1,vt2
    real (kind=real_kind) :: dvsq_int
    real (kind=real_kind) :: vtsq_int

    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met
       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(1,1,i,j)*dv1 + met(1,2,i,j)*dv2
             dvco(i,j,2) = met(2,1,i,j)*dv1 + met(2,2,i,j)*dv2

             vtco(i,j,1) = met(1,1,i,j)*vt1 + met(1,2,i,j)*vt2
             vtco(i,j,2) = met(2,1,i,j)*vt1 + met(2,2,i,j)*vt2

             dvsq(i,j,ie) = SQRT(dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2)
             vtsq(i,j,ie) = SQRT(vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2)

          end do
       end do
    end do

    dvsq_int = global_integral(elem, dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_int = global_integral(elem, vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    l1 = dvsq_int/vtsq_int

  end function l1_vnorm

  ! ==========================================================
  ! l2_snorm:
  !
  ! computes the l2 norm per Williamson et al, p. 218 eq(83)
  ! for a scalar quantity on the pressure grid.
  !
  ! ===========================================================

  function l2_snorm(elem, h,ht,hybrid,npts,nets,nete) result(l2)
    use kinds, only : real_kind
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t), intent(in) :: elem(:)	
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: ht(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: l2   

    ! Local variables

    real (kind=real_kind) :: dh2(npts,npts,nets:nete)
    real (kind=real_kind) :: ht2(npts,npts,nets:nete)
    real (kind=real_kind) :: dh2_int
    real (kind=real_kind) :: ht2_int
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dh2(i,j,ie)=(h(i,j,ie)-ht(i,j,ie))**2
             ht2(i,j,ie)=ht(i,j,ie)**2
          end do
       end do
    end do

    dh2_int = global_integral(elem,dh2(:,:,nets:nete),hybrid,npts,nets,nete)
    ht2_int = global_integral(elem,ht2(:,:,nets:nete),hybrid,npts,nets,nete)

    l2 = SQRT(dh2_int)/SQRT(ht2_int)

  end function l2_snorm

  ! ==========================================================
  ! l2_vnorm:
  !
  ! computes the l2 norm per Williamson et al, p. 219 eq(98)
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function l2_vnorm(elem, v,vt,hybrid,npts,nets,nete) result(l2)
    use kinds, only : real_kind
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t

    type(element_t)      , intent(in), target :: elem(:)
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: l2

    ! Local variables

    real (kind=real_kind), dimension(:,:,:,:), pointer :: met
    real (kind=real_kind) :: dvsq(npts,npts,nets:nete)
    real (kind=real_kind) :: vtsq(npts,npts,nets:nete)
    real (kind=real_kind) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: dv1,dv2
    real (kind=real_kind) :: vt1,vt2
    real (kind=real_kind) :: dvsq_int
    real (kind=real_kind) :: vtsq_int
    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met
       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(1,1,i,j)*dv1 + met(1,2,i,j)*dv2
             dvco(i,j,2) = met(2,1,i,j)*dv1 + met(2,2,i,j)*dv2

             vtco(i,j,1) = met(1,1,i,j)*vt1 + met(1,2,i,j)*vt2
             vtco(i,j,2) = met(2,1,i,j)*vt1 + met(2,2,i,j)*vt2

             dvsq(i,j,ie) = dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2
             vtsq(i,j,ie) = vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2

          end do
       end do
    end do

    dvsq_int = global_integral(elem, dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_int = global_integral(elem, vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    l2 = SQRT(dvsq_int)/SQRT(vtsq_int)

  end function l2_vnorm

  ! ==========================================================
  ! linf_snorm:
  !
  ! computes the l infinity norm per Williamson et al, p. 218 eq(84)
  ! for a scalar quantity on the pressure grid...
  !
  ! ===========================================================

  function linf_snorm(h,ht,hybrid,npts,nets,nete) result(linf)
    use kinds, only : real_kind
    use hybrid_mod, only : hybrid_t
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: h(npts,npts,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: ht(npts,npts,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: linf    

    ! Local variables

    real (kind=real_kind) :: dhabs(npts,npts,nets:nete)
    real (kind=real_kind) :: htabs(npts,npts,nets:nete)
    real (kind=real_kind) :: dhabs_max
    real (kind=real_kind) :: htabs_max
    integer i,j,ie

    do ie=nets,nete
       do j=1,npts
          do i=1,npts
             dhabs(i,j,ie)=ABS(h(i,j,ie)-ht(i,j,ie))
             htabs(i,j,ie)=ABS(ht(i,j,ie))
          end do
       end do
    end do

    dhabs_max = global_maximum(dhabs(:,:,nets:nete),hybrid,npts,nets,nete)
    htabs_max = global_maximum(htabs(:,:,nets:nete),hybrid,npts,nets,nete)

    linf = dhabs_max/htabs_max

  end function linf_snorm


  ! ==========================================================
  ! linf_vnorm:
  !
  ! computes the linf norm per Williamson et al, p. 218 eq(99),
  ! for a contravariant vector quantity on the velocity grid.
  !
  ! ===========================================================

  function linf_vnorm(elem,v,vt,hybrid,npts,nets,nete) result(linf)
    use kinds, only : real_kind
    use hybrid_mod, only : hybrid_t
    use element_mod, only : element_t

    type(element_t)      , intent(in), target :: elem(:) 
    integer              , intent(in) :: npts,nets,nete
    real (kind=real_kind), intent(in) :: v(npts,npts,2,nets:nete)  ! computed soln
    real (kind=real_kind), intent(in) :: vt(npts,npts,2,nets:nete) ! true soln
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind)             :: linf     

    ! Local variables

    real (kind=real_kind), dimension(:,:,:,:), pointer :: met
    real (kind=real_kind) :: dvsq(npts,npts,nets:nete)
    real (kind=real_kind) :: vtsq(npts,npts,nets:nete)
    real (kind=real_kind) :: dvco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: vtco(npts,npts,2)         ! covariant velocity
    real (kind=real_kind) :: dv1,dv2
    real (kind=real_kind) :: vt1,vt2
    real (kind=real_kind) :: dvsq_max
    real (kind=real_kind) :: vtsq_max
    integer i,j,ie

    do ie=nets,nete
       met => elem(ie)%met

       do j=1,npts
          do i=1,npts

             dv1     = v(i,j,1,ie)-vt(i,j,1,ie)
             dv2     = v(i,j,2,ie)-vt(i,j,2,ie)

             vt1     = vt(i,j,1,ie)
             vt2     = vt(i,j,2,ie)

             dvco(i,j,1) = met(1,1,i,j)*dv1 + met(1,2,i,j)*dv2
             dvco(i,j,2) = met(2,1,i,j)*dv1 + met(2,2,i,j)*dv2

             vtco(i,j,1) = met(1,1,i,j)*vt1 + met(1,2,i,j)*vt2
             vtco(i,j,2) = met(2,1,i,j)*vt1 + met(2,2,i,j)*vt2

             dvsq(i,j,ie) = SQRT(dvco(i,j,1)*dv1 + dvco(i,j,2)*dv2)
             vtsq(i,j,ie) = SQRT(vtco(i,j,1)*vt1 + vtco(i,j,2)*vt2)

          end do
       end do
    end do

    dvsq_max = global_maximum(dvsq(:,:,nets:nete),hybrid,npts,nets,nete)
    vtsq_max = global_maximum(vtsq(:,:,nets:nete),hybrid,npts,nets,nete)

    linf = dvsq_max/vtsq_max

  end function linf_vnorm

end module global_norms_mod
