#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module derivative_mod
  use kinds, only : real_kind, longdouble_kind
  use dimensions_mod, only : nv, np
  use quadrature_mod, only : quadrature_t, gauss, gausslobatto,legendre, jacobi

  ! needed for spherical differential operators:
  use physical_constants, only : rearth 
  use element_mod, only : element_t

implicit none
private

  type, public :: derivative_t
     sequence
     real (kind=real_kind) :: Dvv(nv,nv)
     real (kind=real_kind) :: Dvv_twt(nv,nv)
     real (kind=real_kind) :: Mvv_twt(nv,nv)
     real (kind=real_kind) :: vvtemp(nv,nv)
     real (kind=real_kind) :: vvtempt(nv,nv,2)
     real (kind=real_kind) :: M(nv,nv-2)
  end type derivative_t

  type, public :: derivative_stag_t
     sequence
     real (kind=real_kind) :: D(nv,np)
     real (kind=real_kind) :: M(nv,np)
     real (kind=real_kind) :: Dpv(np,nv)
     real (kind=real_kind) :: D_twt(np,nv)
     real (kind=real_kind) :: M_twt(np,nv)
     real (kind=real_kind) :: M_t(np,nv)
     real (kind=real_kind) :: vtemp(nv,np,2)
     real (kind=real_kind) :: vtempt(np,nv,2)
  end type derivative_stag_t

! ======================================
! Public Interfaces
! ======================================

  public :: derivinit
  public :: deriv_print

  public :: gradient
  public :: gradient_wk
  public :: vorticity
  public :: divergence

  public :: interpolate_v2p
  public :: interpolate_p2v

  interface divergence
      module procedure divergence_nonstag
      module procedure divergence_stag
  end interface

  interface gradient
      module procedure gradient_str_nonstag
      module procedure gradient_str_stag
  end interface

  interface gradient_wk
      module procedure gradient_wk_nonstag
      module procedure gradient_wk_stag
  end interface

  private :: dmatinit
  private :: v2pinit
  private :: p2vinit
  private :: dvvinit
  private :: dpvinit

! these routines compute spherical differential operators as opposed to
! the gnomonic coordinate operators above.  Vectors (input or output)
! are always expressed in lat-lon coordinates
  public  :: gradient_sphere
  public  :: vorticity_sphere
  public  :: divergence_sphere
  public  :: curl_sphere
  public  :: divergence_sphere_wk
  public  :: laplace_sphere_wk
  public  :: vlaplace_sphere_wk




contains

! ==========================================
! derivinit:
!
! Initialize the matrices for taking 
! derivatives and interpolating
! ==========================================

  subroutine derivinit(deriv,deriv_stag)
    type (derivative_t)      :: deriv
    type (derivative_stag_t), optional :: deriv_stag

    ! Local variables
    
    type (quadrature_t) :: gp   ! Quadrature points and weights on pressure grid
    type (quadrature_t) :: gv   ! Quadrature points and weights on velocity grid
    real (kind=longdouble_kind) :: dmat(nv,np)
    real (kind=longdouble_kind) :: dpv(np,nv)
    real (kind=longdouble_kind) :: v2p(nv,np)
    real (kind=longdouble_kind) :: p2v(np,nv)
    real (kind=longdouble_kind) :: dvv(nv,nv)
    real (kind=longdouble_kind) :: v2v(nv,nv)

    integer i,j

    ! ============================================
    ! initialize matrices in longdouble_kind precision
    ! and transfer results into real_kind
    ! floating point precision
    ! ============================================

    if (present(deriv_stag) .and. (np==nv-2)) then
       gp=gauss(np)
       call dmatinit(dmat)
       call v2pinit(v2p)
       call p2vinit(p2v)

       deriv_stag%D(:,:)     = dmat(:,:)
       deriv_stag%M(:,:)     = v2p(:,:)
       deriv_stag%M_t(:,:)   = p2v(:,:)

       ! =====================================
       ! Premultiply dmat,v2p to initialize
       ! weak gradient...
       ! =====================================

       do i=1,np
          do j=1,nv
             dmat(j,i) = dmat(j,i)*gp%weights(i)
             v2p(j,i)  = v2p(j,i)*gp%weights(i)
          end do
       end do

       deriv_stag%D_twt = TRANSPOSE(dmat)
       deriv_stag%M_twt = TRANSPOSE(v2p)

       call dpvinit(dpv)
       deriv_stag%Dpv(:,:) = dpv(:,:)

       deallocate(gp%points)
       deallocate(gp%weights)

    endif

    gv=gausslobatto(nv)

    call dvvinit(dvv)
    deriv%Dvv(:,:)   = dvv(:,:)


    v2v = 0.0D0
    do i=1,nv
       v2v(i,i) = gv%weights(i)
    end do

    do i=1,nv
       do j=1,nv
          dvv(j,i) = dvv(j,i)*gv%weights(i)
       end do
    end do

    deriv%Dvv_twt = TRANSPOSE(dvv)
    deriv%Mvv_twt = v2v

    deallocate(gv%points)
    deallocate(gv%weights)

  end subroutine derivinit

  subroutine deriv_print(deriv,deriv_stag)
    type (derivative_t) :: deriv
    type (derivative_stag_t), optional :: deriv_stag
    
    ! Local variables

    integer j
    print *,"Derivative Matrix Dvv"
    do j=1,nv
       write(6,*)deriv%Dvv(:,j)
    end do

    print *,"Weak Derivative Matrix Dvv_twt"
    do j=1,nv
       write(6,*)deriv%Dvv_twt(:,j)
    end do

    if (present(deriv_stag) .and. (np == nv-2)) then
       print *,"Staggered Derivative Matrix D"
       do j=1,np
          write(6,*)deriv_stag%D(:,j)
       end do

       print *,"Staggered Interpolation Matrix M"
       do j=1,np
          write(6,*)deriv_stag%M(:,j), SUM(deriv_stag%M(:,j))
       end do

       print *, "SUM of staggered M matrix columns (v->p) check..."
       do j=1,np
          write(6,*)j,SUM(deriv_stag%M(:,j))
       end do
       print *, "SUM of staggered M_t matrix columns (p->v) check ..."
       do j=1,nv
          write(6,*)j,SUM(deriv_stag%M_t(:,j))
       end do

    end if

  end subroutine deriv_print

! =======================================
! dmatinit:
!
! Compute rectangular v->p 
! derivative matrix (dmat)
! =======================================

  subroutine dmatinit(dmat)

    real (kind=longdouble_kind) :: dmat(nv,np)

    ! Local variables

    type (quadrature_t) :: gll
    type (quadrature_t) :: gs

    integer i,j
    real(kind=longdouble_kind)  fact,f1,f2
    real(kind=longdouble_kind)  func0,func1
    real(kind=longdouble_kind)  dis,c0,c1

    real(kind=longdouble_kind)  :: leg(nv,nv)
    real(kind=longdouble_kind)  ::  jac(0:nv-1)
    real(kind=longdouble_kind)  :: djac(0:nv-1)

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

    gll= gausslobatto(nv)
    gs = gauss(np)

    ! =============================================================
    ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
    ! =============================================================

    do i=1,nv
       leg(:,i) = legendre(gll%points(i),nv-1)
    end do

    ! ================================================================
    !  Derivatives of velocity cardinal functions on pressure grid
    !  d(i,j) = D(j,i) = D' (D-transpose) since D(i,j) = dh_j(x_i)/dx
    ! ================================================================

    fact = nv*(nv-1)

    do j=1,np
       call jacobi(nv-1,gs%points(j),c0,c0,jac(0:nv-1),djac(0:nv-1))
       func0 =  jac(nv-1)
       func1 = djac(nv-1)
       f1 = fact*func0
       f2 = (c1 - gs%points(j))*(c1 + gs%points(j)) * func1
       do i = 1, nv
          if ( gs%points(j) /= gll%points(i) ) then
             dis = gs%points(j) - gll%points(i)
             dmat(i,j) = func0 / ( leg(nv,i)*dis ) + f2 / (fact*leg(nv,i)*dis*dis)
!!! OTHER             dmat(i,j) = (1.0D0/(fact*leg(nv,i)*dis*dis))* (func0*fact*dis + f2)
          else
             dmat(i,j) = c0
          endif
       end do
    end do

    deallocate(gll%points)
    deallocate(gll%weights)

    deallocate(gs%points)
    deallocate(gs%weights)

  end subroutine dmatinit

! =======================================
! dpvinit:
!
! Compute rectangular p->v
! derivative matrix (dmat) 
! for strong gradients
! =======================================

  subroutine dpvinit(dmat)

    real (kind=longdouble_kind) :: dmat(np,nv)

    ! Local variables

    type (quadrature_t) :: gll
    type (quadrature_t) :: gs

    integer i,j
    real(kind=longdouble_kind)  dis,c0,c1

    real(kind=longdouble_kind)  :: legv(0:np,nv)
    real(kind=longdouble_kind)  :: dlegv(0:np,nv)

    real(kind=longdouble_kind)  :: leg(0:np)
    real(kind=longdouble_kind)  :: dleg(0:np)

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

    gll= gausslobatto(nv)
    gs = gauss(np)

    ! =============================================================
    ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
    ! =============================================================

    do i=1,nv
       call jacobi(np,gll%points(i),c0,c0,legv(0:np,i),dlegv(0:np,i))
    end do

    ! ================================================================
    !  Derivatives of velocity cardinal functions on pressure grid
    !  d(i,j) = D(j,i) = D' (D-transpose) since D(i,j) = dh_j(x_i)/dx
    ! ================================================================

    do j=1,np
       call jacobi(np,gs%points(j),c0,c0,leg(0:np),dleg(0:np))
       do i = 1, nv
          if ( gs%points(j) /= gll%points(i) ) then
             dis = gll%points(i) - gs%points(j)
             dmat(j,i) = dlegv(np,i)/( dleg(np)*dis ) -  legv(np,i)/ (dleg(np)*dis*dis)
          else
             dmat(j,i) = c0
          endif
       end do
    end do

    deallocate(gll%points)
    deallocate(gll%weights)

    deallocate(gs%points)
    deallocate(gs%weights)

  end subroutine dpvinit

! =======================================
! v2pinit:
!
! Compute interpolation matrix from 
! velocity to pressure grid (v2p)
! =======================================

  subroutine v2pinit(v2p)

    real(kind=longdouble_kind)  ::  v2p(nv,np)

    ! Local variables

    type (quadrature_t) :: gll
    type (quadrature_t) :: gs

    integer i,j
    real(kind=longdouble_kind)  fact,f1
    real(kind=longdouble_kind)  func0,func1

    real(kind=longdouble_kind)  :: leg(nv,nv)
    real(kind=longdouble_kind)  ::  jac(0:nv-1)
    real(kind=longdouble_kind)  :: djac(0:nv-1)
    real(kind=longdouble_kind)  :: c0,c1

    gll= gausslobatto(nv)
    gs = gauss(np)

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

    ! ==============================================================
    ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
    ! ==============================================================

    fact = -nv*(nv-1)
    do i=1,nv
       leg(:,i) = legendre(gll%points(i),nv-1)
       leg(nv,i) = fact * leg(nv,i)
    end do

    ! ===================================================
    !  Velocity cardinal functions on pressure grid
    ! ===================================================

    do j=1,np
       call jacobi(nv-1,gs%points(j),c0,c0,jac(0:nv-1),djac(0:nv-1))
       func0 = jac(nv-1)
       func1 = djac(nv-1)
       f1 = (c1 - gs%points(j)**2) * func1
       do i = 1, nv
          if ( gs%points(j) /= gll%points(i) ) then
             v2p(i,j) = f1 / ( leg(nv,i) * (gs%points(j)-gll%points(i)))
          else
             v2p(i,j) = c1
          endif
       end do
    end do

    deallocate(gll%points)
    deallocate(gll%weights)

    deallocate(gs%points)
    deallocate(gs%weights)

  end subroutine v2pinit

! =======================================
! p2vinit:
!
! Compute interpolation matrix from 
! pressure grid to velocity grid (p2v)
! Note: this is NOT (v2p)^T.
!
! =======================================

  subroutine p2vinit(p2v)

    real(kind=longdouble_kind)  ::  p2v(np,nv)

    ! Local variables

    type (quadrature_t) :: gll
    type (quadrature_t) :: gs

    integer i,j
    real(kind=longdouble_kind)  func0,func1

    real(kind=longdouble_kind)  :: leg(np+1,nv)
    real(kind=longdouble_kind)  ::  jac(0:np)
    real(kind=longdouble_kind)  :: djac(0:np)
    real(kind=longdouble_kind)  :: c0,c1

    gll= gausslobatto(nv)
    gs = gauss(np)

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

!############ ORIGINAL CODE FROM STEVE THOMAS #######################
!     Interpolation matrix from pressure to velocity grid
!     NOTE: Gauss-Legendre Cardinal Functions !!
!
!      do i = 1,npts
!         call LegAll(xvg(i),nptp,leg(1,i))
!      enddo
!
!     Pressure cardinal functions on velocity grid
!     NOTE: I(j,i) = h_j(x_i) = p2v(j,i)
!
!      do i = 1, nptp
!         call Jacobi(func0,func1,func2,xpg(i),nptp)
!         do j = 1, npts
!            if ( xvg(j) /= xpg(i) ) then
!               p2v(i,j) = leg(nptp+1,j) / ( func1 * (xvg(j)-xpg(i)))
!            else
!               p2v(i,j) = 1.0d0
!            endif
!         enddo
!      enddo
!####################################################################

    ! ==============================================================
    ! Compute Legendre polynomials on Gauss-Lobatto grid (velocity)
    ! ==============================================================

    do i=1,nv
       leg(:,i) = legendre(gll%points(i),np)
    end do

    ! ===================================================
    !  Pressure cardinal functions on velocity grid
    !  NOTE: I(j,i) = h_j(x_i) = p2v(j,i)
    ! ===================================================

    do i=1,np
       call jacobi(np,gs%points(i),c0,c0,jac(0:np),djac(0:np))
       func0 = jac(np)
       func1 = djac(np)
       do j = 1, nv
          if ( gll%points(j) /= gs%points(i) ) then
             p2v(i,j) = leg(np+1,j) / ( func1 * (gll%points(j)-gs%points(i)))
          else
             p2v(i,j) = c1
          endif
       end do
    end do

    deallocate(gll%points)
    deallocate(gll%weights)

    deallocate(gs%points)
    deallocate(gs%weights)

  end subroutine p2vinit

#if 0

! =======================================
! p2vinit:
!
! Compute interpolation matrix from 
! velocity to pressure grid (v2p)
! =======================================

  subroutine p2vinit(p2v)

    real(kind=longdouble_kind)  ::  p2v(np,nv)

    ! Local variables

    type (quadrature_t) :: gll
    type (quadrature_t) :: gs

    integer i,j
    real(kind=longdouble_kind)  func0,func1
    real(kind=longdouble_kind)  dis

    real(kind=longdouble_kind)  ::  jac(0:np)
    real(kind=longdouble_kind)  :: djac(0:np)
    real(kind=longdouble_kind)  :: c0,c1

    gll= gausslobatto(nv)
    gs = gauss(np)

    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind

    ! ===================================================
    !  Pressure cardinal functions on velocity grid
    ! ===================================================

    do j=1,nv
       call jacobi(np,gll%points(j),c0,c0,jac(0:np),djac(0:np))
       func0 = jac(np)
       func1 = djac(np)
       do i = 1, np
          if ( gs%points(i) /= gll%points(j) ) then
             p2v(i,j) = func0 / ( func1 * (gll%points(j) - gs%points(i)))
          else
             p2v(i,j) = c1
          endif
       end do
    end do

    deallocate(gll%points)
    deallocate(gll%weights)

    deallocate(gs%points)
    deallocate(gs%weights)

  end subroutine p2vinit
#endif

! =======================================
! dvvinit:
!
! Compute rectangular v->v
! derivative matrix (dvv)
! =======================================

  subroutine dvvinit(dvv)

    real(kind=longdouble_kind)  ::  dvv(nv,nv)

    ! Local variables

    type (quadrature_t)   :: gll
    real(kind=longdouble_kind)  :: leg(nv,nv)
    real(kind=longdouble_kind)  :: c0,c1,c4

    integer i,j

    gll= gausslobatto(nv)
    c0 = 0.0_longdouble_kind
    c1 = 1.0_longdouble_kind
    c4 = 4.0_longdouble_kind

    do i=1,nv
       leg(:,i) = legendre(gll%points(i),nv-1)
    end do

    dvv(:,:) = c0
    do j=1,nv
       do i=1,j-1
          dvv(j,i) = (c1/(gll%points(i)-gll%points(j)))*leg(nv,i)/leg(nv,j)
       end do
       dvv(j,j) = c0
       do i=j+1,nv
          dvv(j,i) = (c1/(gll%points(i)-gll%points(j)))*leg(nv,i)/leg(nv,j)
       end do
    end do


    dvv(nv,nv) = + nv*(nv-1)/c4
    dvv(1,1)   = - nv*(nv-1)/c4

    deallocate(gll%points)
    deallocate(gll%weights)

  end subroutine dvvinit

!  ================================================
!  divergence_stag: 
!
!  Compute divergence (maps v grid -> p grid)
!  ================================================

  function divergence_stag(v,deriv,rdx,rdy) result(div)

    real(kind=real_kind), intent(in) :: v(nv,nv,2)
    type (derivative_stag_t)         :: deriv
    real(kind=real_kind), intent(in) :: rdx
    real(kind=real_kind), intent(in) :: rdy

    real(kind=real_kind) :: div(np,np)

    ! Local

    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11

#ifdef DEBUG
    print *, "divergence_stag"
#endif
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    !JMD====================================
    !JMD  2*nv*nv*np Flops
    !JMD====================================
    do j=1,nv,2
       do l=1,np,2

          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do i=1,nv
             sumx00 = sumx00 + deriv%D(i,l  )*v(i,j  ,1)
             sumx01 = sumx01 + deriv%D(i,l+1)*v(i,j  ,1)
             sumx10 = sumx10 + deriv%D(i,l  )*v(i,j+1,1)
             sumx11 = sumx11 + deriv%D(i,l+1)*v(i,j+1,1)

             sumy00 = sumy00 + deriv%M(i,l  )*v(i,j  ,2)
             sumy01 = sumy01 + deriv%M(i,l+1)*v(i,j  ,2)
             sumy10 = sumy10 + deriv%M(i,l  )*v(i,j+1,2)
             sumy11 = sumy11 + deriv%M(i,l+1)*v(i,j+1,2)
          end do

          deriv%vtemp(j  ,l  ,1) = sumx00
          deriv%vtemp(j  ,l+1,1) = sumx01
          deriv%vtemp(j+1,l  ,1) = sumx10
          deriv%vtemp(j+1,l+1,1) = sumx11

          deriv%vtemp(j  ,l  ,2) = sumy00
          deriv%vtemp(j  ,l+1,2) = sumy01
          deriv%vtemp(j+1,l  ,2) = sumy10
          deriv%vtemp(j+1,l+1,2) = sumy11

       end do
    end do


    !JMD====================================
    !JMD  2*nv*np*np Flops
    !JMD====================================
    do j=1,np,2
       do i=1,np,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do l=1,nv
             sumx00 = sumx00 +  deriv%M(l,j  )*deriv%vtemp(l,i  ,1)
             sumx01 = sumx01 +  deriv%M(l,j+1)*deriv%vtemp(l,i  ,1)
             sumx10 = sumx10 +  deriv%M(l,j  )*deriv%vtemp(l,i+1,1)
             sumx11 = sumx11 +  deriv%M(l,j+1)*deriv%vtemp(l,i+1,1)

             sumy00 = sumy00 +  deriv%D(l,j  )*deriv%vtemp(l,i  ,2)
             sumy01 = sumy01 +  deriv%D(l,j+1)*deriv%vtemp(l,i  ,2)
             sumy10 = sumy10 +  deriv%D(l,j  )*deriv%vtemp(l,i+1,2)
             sumy11 = sumy11 +  deriv%D(l,j+1)*deriv%vtemp(l,i+1,2)
          end do

          div(i  ,j  ) = rdx*sumx00 + rdy*sumy00
          div(i  ,j+1) = rdx*sumx01 + rdy*sumy01
          div(i+1,j  ) = rdx*sumx10 + rdy*sumy10
          div(i+1,j+1) = rdx*sumx11 + rdy*sumy11

       end do
    end do
else
     do j=1,nv
        do l=1,np
 
           sumx00=0.0d0
           sumy00=0.0d0
           do i=1,nv
              sumx00 = sumx00 + deriv%D(i,l  )*v(i,j  ,1)
              sumy00 = sumy00 + deriv%M(i,l  )*v(i,j  ,2)
 	   enddo
          deriv%vtemp(j  ,l  ,1) = sumx00
          deriv%vtemp(j  ,l  ,2) = sumy00
       enddo
    enddo
    do j=1,np
       do i=1,np
          sumx00=0.0d0
	  sumy00=0.0d0
          do l=1,nv
             sumx00 = sumx00 +  deriv%M(l,j  )*deriv%vtemp(l,i  ,1)
	     sumy00 = sumy00 +  deriv%D(l,j  )*deriv%vtemp(l,i  ,2)
	  enddo
          div(i  ,j  ) = rdx*sumx00 + rdy*sumy00

	enddo
    enddo
endif

  end function divergence_stag

!  ================================================
!  divergence_nonstag: 
!
!  Compute divergence (maps v->v)
!  ================================================

  function divergence_nonstag(v,deriv,rdx,rdy) result(div)

    real(kind=real_kind), intent(in) :: v(nv,nv,2)
    type (derivative_t)              :: deriv
    real(kind=real_kind), intent(in) :: rdx
    real(kind=real_kind), intent(in) :: rdy

    real(kind=real_kind) :: div(nv,nv)

    ! Local

    integer i
    integer j
    integer l

    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind) ::  dudx00,dudx01
    real(kind=real_kind) ::  dudx10,dudx11

    real(kind=real_kind) ::  dvdy00,dvdy01
    real(kind=real_kind) ::  dvdy10,dvdy11

if(modulo(nv,2) .eq. 0 .and. UseUnroll) then
! this is just loop unrolling - a good compiler should do it for you jpe
       do j=1,nv,2
          do l=1,nv,2

             dudx00=0.0d0
             dudx01=0.0d0
             dudx10=0.0d0
             dudx11=0.0d0

             dvdy00=0.0d0
             dvdy01=0.0d0
             dvdy10=0.0d0
             dvdy11=0.0d0

             do i=1,nv
                
                dudx00 = dudx00 + deriv%Dvv(i,l  )*v(i,j  ,1)
                dudx01 = dudx01 + deriv%Dvv(i,l+1)*v(i,j  ,1)
                dudx10 = dudx10 + deriv%Dvv(i,l  )*v(i,j+1,1)
                dudx11 = dudx11 + deriv%Dvv(i,l+1)*v(i,j+1,1)
                
                dvdy00 = dvdy00 + deriv%Dvv(i,l  )*v(j  ,i,2)
                dvdy01 = dvdy01 + deriv%Dvv(i,l+1)*v(j  ,i,2)
                dvdy10 = dvdy10 + deriv%Dvv(i,l  )*v(j+1,i,2)
                dvdy11 = dvdy11 + deriv%Dvv(i,l+1)*v(j+1,i,2)

             end do

             div(l  ,j  ) = dudx00
             div(l+1,j  ) = dudx01
             div(l  ,j+1) = dudx10
             div(l+1,j+1) = dudx11

             deriv%vvtemp(j  ,l  ) = dvdy00
             deriv%vvtemp(j  ,l+1) = dvdy01
             deriv%vvtemp(j+1,l  ) = dvdy10
             deriv%vvtemp(j+1,l+1) = dvdy11

          end do
       end do
    else

       do j=1,nv
          do l=1,nv
             dudx00=0.0d0
             dvdy00=0.0d0

             do i=1,nv
                dudx00 = dudx00 + deriv%Dvv(i,l  )*v(i,j  ,1)
                dvdy00 = dvdy00 + deriv%Dvv(i,l  )*v(j  ,i,2)
             end do

             div(l  ,j  ) = dudx00
             deriv%vvtemp(j  ,l  ) = dvdy00


          end do
       end do
    end if
    do j=1,nv
       do i=1,nv
          div(i,j)=rdx*div(i,j)+rdy*deriv%vvtemp(i,j)
       end do
    end do

  end function divergence_nonstag

!  ================================================
!  gradient_wk_stag:
! 
!  Compute the weak form gradient:
!  maps scalar field on the pressure grid to the
!  velocity grid
!  ================================================

  function gradient_wk_stag(p,deriv,dx,dy) result(dp)

    type (derivative_stag_t)         :: deriv
    real(kind=real_kind), intent(in) :: p(np,np)
    real(kind=real_kind), intent(in) :: dx
    real(kind=real_kind), intent(in) :: dy

    real(kind=real_kind)             :: dp(nv,nv,2)

    ! Local
      
    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11

#ifdef DEBUG
    print *, "gradient_wk_stag"
#endif
    !JMD ================================
    !JMD 2*nv*np*np Flops 
    !JMD ================================

if(MODULO(np,2) == 0 .and. UseUnroll) then 


    do j=1,np,2
       do l=1,nv,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do i=1,np
             sumx00 = sumx00 + deriv%D_twt(i,l  )*p(i,j  )
             sumx01 = sumx01 + deriv%D_twt(i,l+1)*p(i,j  )
             sumx10 = sumx10 + deriv%D_twt(i,l  )*p(i,j+1)
             sumx11 = sumx11 + deriv%D_twt(i,l+1)*p(i,j+1)

             sumy00 = sumy00 + deriv%M_twt(i,l  )*p(i,j  )
             sumy01 = sumy01 + deriv%M_twt(i,l+1)*p(i,j  )
             sumy10 = sumy10 + deriv%M_twt(i,l  )*p(i,j+1)
             sumy11 = sumy11 + deriv%M_twt(i,l+1)*p(i,j+1)
          end do

          deriv%vtempt(j  ,l  ,1) = sumx00
          deriv%vtempt(j  ,l+1,1) = sumx01
          deriv%vtempt(j+1,l  ,1) = sumx10
          deriv%vtempt(j+1,l+1,1) = sumx11

          deriv%vtempt(j  ,l  ,2) = sumy00
          deriv%vtempt(j  ,l+1,2) = sumy01
          deriv%vtempt(j+1,l  ,2) = sumy10
          deriv%vtempt(j+1,l+1,2) = sumy11
       end do
    end do

	
    !JMD ================================
    !JMD 2*nv*nv*np Flops 
    !JMD ================================

    do j=1,nv,2
       do i=1,nv,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do l=1,np
             sumx00 = sumx00 +  deriv%M_twt(l,j  )*deriv%vtempt(l,i  ,1)
             sumx01 = sumx01 +  deriv%M_twt(l,j+1)*deriv%vtempt(l,i  ,1)
             sumx10 = sumx10 +  deriv%M_twt(l,j  )*deriv%vtempt(l,i+1,1)
             sumx11 = sumx11 +  deriv%M_twt(l,j+1)*deriv%vtempt(l,i+1,1)

             sumy00 = sumy00 +  deriv%D_twt(l,j  )*deriv%vtempt(l,i  ,2)
             sumy01 = sumy01 +  deriv%D_twt(l,j+1)*deriv%vtempt(l,i  ,2)
             sumy10 = sumy10 +  deriv%D_twt(l,j  )*deriv%vtempt(l,i+1,2)
             sumy11 = sumy11 +  deriv%D_twt(l,j+1)*deriv%vtempt(l,i+1,2)
          end do

          dp(i  ,j  ,1) = dy*sumx00
          dp(i  ,j+1,1) = dy*sumx01
          dp(i+1,j  ,1) = dy*sumx10
          dp(i+1,j+1,1) = dy*sumx11

          dp(i  ,j  ,2) = dx*sumy00
          dp(i  ,j+1,2) = dx*sumy01
          dp(i+1,j  ,2) = dx*sumy10
          dp(i+1,j+1,2) = dx*sumy11

       end do
    end do
else
    do j=1,np
       do l=1,nv
 	  sumx00=0.0d0
          sumy00=0.0d0
           do i=1,np
              sumx00 = sumx00 + deriv%D_twt(i,l  )*p(i,j  )
              sumy00 = sumy00 + deriv%M_twt(i,l  )*p(i,j  )
 	  enddo
           deriv%vtempt(j  ,l  ,1) = sumx00
           deriv%vtempt(j  ,l  ,2) = sumy00
        enddo
    enddo
    do j=1,nv
       do i=1,nv
          sumx00=0.0d0
	  sumy00=0.0d0
          do l=1,np
             sumx00 = sumx00 +  deriv%M_twt(l,j  )*deriv%vtempt(l,i  ,1)
	     sumy00 = sumy00 +  deriv%D_twt(l,j  )*deriv%vtempt(l,i  ,2)
	  enddo
	  dp(i  ,j  ,1) = dy*sumx00
          dp(i  ,j  ,2) = dx*sumy00
      enddo
    enddo
endif


  end function gradient_wk_stag

!  ================================================
!  gradient_wk_nonstag:
! 
!  Compute the weak form gradient:
!  maps scalar field on the Gauss-Lobatto grid to the
!  weak gradient on the Gauss-Lobbatto grid
!  ================================================

  function gradient_wk_nonstag(p,deriv,dx,dy) result(dp)

    type (derivative_t)         :: deriv
    real(kind=real_kind), intent(in) :: p(nv,nv)
    real(kind=real_kind), intent(in) :: dx
    real(kind=real_kind), intent(in) :: dy

    real(kind=real_kind)             :: dp(nv,nv,2)

    ! Local
      
    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11

    !JMD ================================
    !JMD 2*nv*nv*nv Flops 
    !JMD ================================

!   print *, "gradient_wk_nonstag"
    if(modulo(nv,2) .eq. 0 .and. UseUnroll) then
! this is just loop unrolling - a good compiler should do it for you jpe

       do j=1,nv,2
          do l=1,nv,2
             sumx00=0.0d0
             sumx01=0.0d0
             sumx10=0.0d0
             sumx11=0.0d0

             sumy00=0.0d0
             sumy01=0.0d0
             sumy10=0.0d0
             sumy11=0.0d0

             do i=1,nv
                sumx00 = sumx00 + deriv%Dvv_twt(i,l  )*p(i,j  )
                sumx01 = sumx01 + deriv%Dvv_twt(i,l+1)*p(i,j  )
                sumx10 = sumx10 + deriv%Dvv_twt(i,l  )*p(i,j+1)
                sumx11 = sumx11 + deriv%Dvv_twt(i,l+1)*p(i,j+1)

                sumy00 = sumy00 + deriv%Mvv_twt(i,l  )*p(i,j  )
                sumy01 = sumy01 + deriv%Mvv_twt(i,l+1)*p(i,j  )
                sumy10 = sumy10 + deriv%Mvv_twt(i,l  )*p(i,j+1)
                sumy11 = sumy11 + deriv%Mvv_twt(i,l+1)*p(i,j+1)
             end do

             deriv%vvtempt(j  ,l  ,1) = sumx00
             deriv%vvtempt(j  ,l+1,1) = sumx01
             deriv%vvtempt(j+1,l  ,1) = sumx10
             deriv%vvtempt(j+1,l+1,1) = sumx11

             deriv%vvtempt(j  ,l  ,2) = sumy00
             deriv%vvtempt(j  ,l+1,2) = sumy01
             deriv%vvtempt(j+1,l  ,2) = sumy10
             deriv%vvtempt(j+1,l+1,2) = sumy11

          end do
       end do
       ! vvtempt1 = p'*Dvv_twt
       ! vvtempt2 = p'*Mvv_twt
       ! dp1 = dy*Mvv_twt*vvtempt1' = dy*Mvv_twt*(p'*Dvv_twt)' = dy*Mvv_twt*Dvv_twt'*p
       ! dp2 = dx*Dvv_twt*vvtempt2' = dx*Dvv_twt*(p'*Mvv_twt)' = dx*Dvv_twt*Mvv_twt'*p
       !     New formulation 
       ! dp1 = dy*MvvDvvt*p
       ! dp2 = dx*DvvMvvt*p
       ! MvvDvvt = Mvv_twt*Dvv_twt'
       ! DvvMvvt = Dvv_twt*Mvv_twt'


       !JMD ================================
       !JMD 2*nv*nv*nv Flops 
       !JMD ================================

       do j=1,nv,2
          do i=1,nv,2
             sumx00=0.0d0
             sumx01=0.0d0
             sumx10=0.0d0
             sumx11=0.0d0
             
             sumy00=0.0d0
             sumy01=0.0d0
             sumy10=0.0d0
             sumy11=0.0d0
             
             do l=1,nv
                sumx00 = sumx00 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i  ,1)
                sumx01 = sumx01 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i  ,1)
                sumx10 = sumx10 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i+1,1)
                sumx11 = sumx11 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i+1,1)

                sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i  ,2)
                sumy01 = sumy01 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i  ,2)
                sumy10 = sumy10 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i+1,2)
                sumy11 = sumy11 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i+1,2)
             end do
             
             dp(i  ,j  ,1) = dy*sumx00
             dp(i  ,j+1,1) = dy*sumx01
             dp(i+1,j  ,1) = dy*sumx10
             dp(i+1,j+1,1) = dy*sumx11
             
             dp(i  ,j  ,2) = dx*sumy00
             dp(i  ,j+1,2) = dx*sumy01
             dp(i+1,j  ,2) = dx*sumy10
             dp(i+1,j+1,2) = dx*sumy11

          end do
       end do
    else

       do j=1,nv
          do l=1,nv
             sumx00=0.0d0

             sumy00=0.0d0

             do i=1,nv
                sumx00 = sumx00 + deriv%Dvv_twt(i,l  )*p(i,j  )

                sumy00 = sumy00 + deriv%Mvv_twt(i,l  )*p(i,j  )
             end do

             deriv%vvtempt(j  ,l  ,1) = sumx00

             deriv%vvtempt(j  ,l  ,2) = sumy00

          end do
       end do

       !JMD ================================
       !JMD 2*nv*nv*nv Flops 
       !JMD ================================

       do j=1,nv
          do i=1,nv
             sumx00=0.0d0
             
             sumy00=0.0d0
             
             do l=1,nv
                sumx00 = sumx00 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i  ,1)

                sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i  ,2)
             end do
             
             dp(i  ,j  ,1) = dy*sumx00
             
             dp(i  ,j  ,2) = dx*sumy00

          end do
       end do
    end if
  end function gradient_wk_nonstag

!  ================================================
!  gradient_str_stag:
! 
!  Compute the *strong* form gradient:
!  maps scalar field on the pressure grid to the
!  velocity grid
!  ================================================

  function gradient_str_stag(p,deriv,rdx,rdy) result(dp)

    type (derivative_stag_t)         :: deriv
    real(kind=real_kind), intent(in) :: p(np,np)
    real(kind=real_kind), intent(in) :: rdx
    real(kind=real_kind), intent(in) :: rdy

    real(kind=real_kind)             :: dp(nv,nv,2)

    ! Local
      
    integer i
    integer j
    integer l

    logical, parameter :: UseUnroll=.TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11

#ifdef DEBUG
    print *, "gradient_str_stag"
#endif
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,nv,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do i=1,np
             sumx00 = sumx00 + deriv%Dpv(i,l  )*p(i,j  )
             sumx01 = sumx01 + deriv%Dpv(i,l+1)*p(i,j  )
             sumx10 = sumx10 + deriv%Dpv(i,l  )*p(i,j+1)
             sumx11 = sumx11 + deriv%Dpv(i,l+1)*p(i,j+1)

             sumy00 = sumy00 + deriv%M_t(i,l  )*p(i,j  )
             sumy01 = sumy01 + deriv%M_t(i,l+1)*p(i,j  )
             sumy10 = sumy10 + deriv%M_t(i,l  )*p(i,j+1)
             sumy11 = sumy11 + deriv%M_t(i,l+1)*p(i,j+1)
          end do

          deriv%vtempt(j  ,l  ,1) = sumx00
          deriv%vtempt(j  ,l+1,1) = sumx01
          deriv%vtempt(j+1,l  ,1) = sumx10
          deriv%vtempt(j+1,l+1,1) = sumx11

          deriv%vtempt(j  ,l  ,2) = sumy00
          deriv%vtempt(j  ,l+1,2) = sumy01
          deriv%vtempt(j+1,l  ,2) = sumy10
          deriv%vtempt(j+1,l+1,2) = sumy11

       end do
    end do




    do j=1,nv,2
       do i=1,nv,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          sumy00=0.0d0
          sumy01=0.0d0
          sumy10=0.0d0
          sumy11=0.0d0

          do l=1,np
             sumx00 = sumx00 +  deriv%M_t(l,j  )*deriv%vtempt(l,i  ,1)
             sumx01 = sumx01 +  deriv%M_t(l,j+1)*deriv%vtempt(l,i  ,1)
             sumx10 = sumx10 +  deriv%M_t(l,j  )*deriv%vtempt(l,i+1,1)
             sumx11 = sumx11 +  deriv%M_t(l,j+1)*deriv%vtempt(l,i+1,1)

             sumy00 = sumy00 +  deriv%Dpv(l,j  )*deriv%vtempt(l,i  ,2)
             sumy01 = sumy01 +  deriv%Dpv(l,j+1)*deriv%vtempt(l,i  ,2)
             sumy10 = sumy10 +  deriv%Dpv(l,j  )*deriv%vtempt(l,i+1,2)
             sumy11 = sumy11 +  deriv%Dpv(l,j+1)*deriv%vtempt(l,i+1,2)
          end do

          dp(i  ,j  ,1) = rdx*sumx00
          dp(i  ,j+1,1) = rdx*sumx01
          dp(i+1,j  ,1) = rdx*sumx10
          dp(i+1,j+1,1) = rdx*sumx11

          dp(i  ,j  ,2) = rdy*sumy00
          dp(i  ,j+1,2) = rdy*sumy01
          dp(i+1,j  ,2) = rdy*sumy10
          dp(i+1,j+1,2) = rdy*sumy11

       end do
    end do
else
    do j=1,np
       do l=1,nv
          sumx00=0.0d0
          sumy00=0.0d0
          do i=1,np
             sumx00 = sumx00 + deriv%Dpv(i,l  )*p(i,j  )
             sumy00 = sumy00 + deriv%M_t(i,l  )*p(i,j  )
   	   enddo
	   deriv%vtempt(j  ,l  ,1) = sumx00
   	   deriv%vtempt(j  ,l  ,2) = sumy00
	enddo
    enddo
    do j=1,nv
       do i=1,nv
          sumx00=0.0d0
          sumy00=0.0d0
          do l=1,np
             sumx00 = sumx00 +  deriv%M_t(l,j  )*deriv%vtempt(l,i  ,1)
             sumy00 = sumy00 +  deriv%Dpv(l,j  )*deriv%vtempt(l,i  ,2)
	  enddo
          dp(i  ,j  ,1) = rdx*sumx00
	  dp(i  ,j  ,2) = rdy*sumy00
       enddo
    enddo
endif

  end function gradient_str_stag

!  ================================================
!  gradient_str_nonstag:
!
!  Compute the *strong* gradient on the velocity grid
!  of a scalar field on the velocity grid
!  ================================================

  function gradient_str_nonstag(s,deriv,rdx,rdy) result(ds)

    type (derivative_t)              :: deriv
    real(kind=real_kind), intent(in) :: s(nv,nv)
    real(kind=real_kind), intent(in) :: rdx
    real(kind=real_kind), intent(in) :: rdy

    real(kind=real_kind) :: ds(nv,nv,2)

    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind) ::  dsdx00,dsdx01
    real(kind=real_kind) ::  dsdx10,dsdx11

    real(kind=real_kind) ::  dsdy00,dsdy01
    real(kind=real_kind) ::  dsdy10,dsdy11
#ifdef DEBUG
    print *, "gradient_str_nonstag"
!   write(17) nv,s,deriv,rdx,rdy
#endif
    if(modulo(nv,2) .eq. 0 .and. UseUnroll) then
       do j=1,nv,2
          do l=1,nv,2
             dsdx00=0.0d0
             dsdx01=0.0d0
             dsdx10=0.0d0
             dsdx11=0.0d0

             dsdy00=0.0d0
             dsdy01=0.0d0
             dsdy10=0.0d0
             dsdy11=0.0d0

             do i=1,nv
                dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
                dsdx01 = dsdx01 + deriv%Dvv(i,l+1)*s(i,j  )
                dsdx10 = dsdx10 + deriv%Dvv(i,l  )*s(i,j+1)
                dsdx11 = dsdx11 + deriv%Dvv(i,l+1)*s(i,j+1)

                dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
                dsdy01 = dsdy01 + deriv%Dvv(i,l+1)*s(j  ,i)
                dsdy10 = dsdy10 + deriv%Dvv(i,l  )*s(j+1,i)
                dsdy11 = dsdy11 + deriv%Dvv(i,l+1)*s(j+1,i)
             end do
#ifdef DEBUG
             if(j.eq.3.and.l.eq.1) then
                print *, dsdx00,rdx,dsdx00*rdx
             endif
#endif
             ds(l  ,j  ,1) = rdx*dsdx00
             ds(l+1,j  ,1) = rdx*dsdx01
             ds(l  ,j+1,1) = rdx*dsdx10
             ds(l+1,j+1,1) = rdx*dsdx11

             ds(j  ,l  ,2) = rdy*dsdy00
             ds(j  ,l+1,2) = rdy*dsdy01
             ds(j+1,l  ,2) = rdy*dsdy10
             ds(j+1,l+1,2) = rdy*dsdy11

          end do

       end do
    else
       do j=1,nv
          do l=1,nv
             dsdx00=0.0d0

             dsdy00=0.0d0

             do i=1,nv
                dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )

                dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
             end do
             ds(l  ,j  ,1) = rdx*dsdx00

             ds(j  ,l  ,2) = rdy*dsdy00

          end do

       end do
    end if
  end function gradient_str_nonstag

!  ================================================
!  vorticity:
!
!  Compute the vorticity of the velocity field on the
!  velocity grid
!  ================================================

  function vorticity(v,deriv,rdx,rdy) result(vort)

    type (derivative_t)              :: deriv
    real(kind=real_kind), intent(in) :: v(nv,nv,2)
    real(kind=real_kind), intent(in) :: rdx
    real(kind=real_kind), intent(in) :: rdy

    real(kind=real_kind) :: vort(nv,nv)

    integer i
    integer j
    integer l
    
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind) ::  dvdx00,dvdx01
    real(kind=real_kind) ::  dvdx10,dvdx11

    real(kind=real_kind) ::  dudy00,dudy01
    real(kind=real_kind) ::  dudy10,dudy11

if(MODULO(nv,2) == 0 .and. UseUnroll) then 
    do j=1,nv,2
       do l=1,nv,2

          dudy00=0.0d0
          dudy01=0.0d0
          dudy10=0.0d0
          dudy11=0.0d0

          dvdx00=0.0d0
          dvdx01=0.0d0
          dvdx10=0.0d0
          dvdx11=0.0d0

          do i=1,nv

             dvdx00 = dvdx00 + deriv%Dvv(i,l  )*v(i,j  ,2)
             dvdx01 = dvdx01 + deriv%Dvv(i,l+1)*v(i,j  ,2)
             dvdx10 = dvdx10 + deriv%Dvv(i,l  )*v(i,j+1,2)
             dvdx11 = dvdx11 + deriv%Dvv(i,l+1)*v(i,j+1,2)

             dudy00 = dudy00 + deriv%Dvv(i,l  )*v(j  ,i,1)
             dudy01 = dudy01 + deriv%Dvv(i,l+1)*v(j  ,i,1)
             dudy10 = dudy10 + deriv%Dvv(i,l  )*v(j+1,i,1)
             dudy11 = dudy11 + deriv%Dvv(i,l+1)*v(j+1,i,1)

          end do

          vort(l  ,j  ) = dvdx00
          vort(l+1,j  ) = dvdx01
          vort(l  ,j+1) = dvdx10
          vort(l+1,j+1) = dvdx11

          deriv%vvtemp(j  ,l  ) = dudy00
          deriv%vvtemp(j  ,l+1) = dudy01
          deriv%vvtemp(j+1,l  ) = dudy10
          deriv%vvtemp(j+1,l+1) = dudy11

        end do
    end do
else
    do j=1,nv
       do l=1,nv

          dudy00=0.0d0
	  dvdx00=0.0d0

          do i=1,nv
             dvdx00 = dvdx00 + deriv%Dvv(i,l  )*v(i,j  ,2)
             dudy00 = dudy00 + deriv%Dvv(i,l  )*v(j  ,i,1)
	  enddo
 
	  vort(l  ,j  ) = dvdx00
	  deriv%vvtemp(j  ,l  ) = dudy00
	enddo
     enddo

endif

    do j=1,nv
       do i=1,nv
          vort(i,j)=rdx*vort(i,j)-rdy*deriv%vvtemp(i,j)
       end do
    end do

  end function vorticity

!  ================================================
!  interpolate_v2p:
!
!  Interpolate a scalar quantity from velocity to 
!  pressure grid
!  ================================================

  function interpolate_v2p(v,deriv) result(p)

    real(kind=real_kind), intent(in) :: v(nv,nv)
    type (derivative_stag_t)         :: deriv

    real(kind=real_kind) :: p(np,np)

    ! Local

    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll=.TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumx10,sumx11

if(MODULO(nv,2) == 0 .and. UseUnroll) then 
    do j=1,nv,2
       do l=1,np,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          do i=1,nv
             sumx00 = sumx00 + deriv%M(i,l  )*v(i,j  )
             sumx01 = sumx01 + deriv%M(i,l+1)*v(i,j  )
             sumx10 = sumx10 + deriv%M(i,l  )*v(i,j+1)
             sumx11 = sumx11 + deriv%M(i,l+1)*v(i,j+1)
          end do

          deriv%vtemp(j  ,l  ,1) = sumx00
          deriv%vtemp(j  ,l+1,1) = sumx01
          deriv%vtemp(j+1,l  ,1) = sumx10
          deriv%vtemp(j+1,l+1,1) = sumx11

       end do
    end do

    do j=1,np,2
       do i=1,np,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          do l=1,nv
             sumx00 = sumx00 + deriv%M(l,j  )*deriv%vtemp(l,i  ,1)
             sumx01 = sumx01 + deriv%M(l,j+1)*deriv%vtemp(l,i  ,1)
             sumx10 = sumx10 + deriv%M(l,j  )*deriv%vtemp(l,i+1,1)
             sumx11 = sumx11 + deriv%M(l,j+1)*deriv%vtemp(l,i+1,1)
          end do

          p(i  ,j  ) = sumx00
          p(i  ,j+1) = sumx01
          p(i+1,j  ) = sumx10
          p(i+1,j+1) = sumx11

       end do
     end do
else
    do j=1,nv
       do l=1,np
          sumx00=0.0d0
          do i=1,nv
             sumx00 = sumx00 + deriv%M(i,l  )*v(i,j  )
	  enddo
	  deriv%vtemp(j  ,l  ,1) = sumx00
	enddo
    enddo
    do j=1,np
       do i=1,np
          sumx00=0.0d0
          do l=1,nv
             sumx00 = sumx00 + deriv%M(l,j  )*deriv%vtemp(l,i  ,1)
	  enddo
          p(i  ,j  ) = sumx00
       enddo
    enddo
endif

  end function interpolate_v2p

!  ================================================
!  interpolate_p2v:
!
!  Interpolate a scalar quantity from pressure to 
!  velocity grid
!
!  ================================================

  function interpolate_p2v(p,deriv) result(v)

    real(kind=real_kind), intent(in) :: p(np,np)
    type (derivative_stag_t)         :: deriv

    real(kind=real_kind) :: v(nv,nv)

    ! Local

    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumx10,sumx11

if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,nv,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          do i=1,np
             sumx00 = sumx00 + deriv%M_t(i,l  )*p(i,j  )
             sumx01 = sumx01 + deriv%M_t(i,l+1)*p(i,j  )
             sumx10 = sumx10 + deriv%M_t(i,l  )*p(i,j+1)
             sumx11 = sumx11 + deriv%M_t(i,l+1)*p(i,j+1)
          end do

          deriv%vtempt(j  ,l  ,1) = sumx00
          deriv%vtempt(j  ,l+1,1) = sumx01
          deriv%vtempt(j+1,l  ,1) = sumx10
          deriv%vtempt(j+1,l+1,1) = sumx11

       end do
    end do

    do j=1,nv,2
       do i=1,nv,2
          sumx00=0.0d0
          sumx01=0.0d0
          sumx10=0.0d0
          sumx11=0.0d0

          do l=1,np
             sumx00 = sumx00 + deriv%M_t(l,j  )*deriv%vtempt(l,i  ,1)
             sumx01 = sumx01 + deriv%M_t(l,j+1)*deriv%vtempt(l,i  ,1)
             sumx10 = sumx10 + deriv%M_t(l,j  )*deriv%vtempt(l,i+1,1)
             sumx11 = sumx11 + deriv%M_t(l,j+1)*deriv%vtempt(l,i+1,1)
          end do

          v(i  ,j  ) = sumx00
          v(i  ,j+1) = sumx01
          v(i+1,j  ) = sumx10
          v(i+1,j+1) = sumx11

       end do
     end do
else
    do j=1,np
       do l=1,nv
          sumx00=0.0d0
          do i=1,np
             sumx00 = sumx00 + deriv%M_t(i,l  )*p(i,j  )
	  end do
	  deriv%vtempt(j  ,l  ,1) = sumx00
       enddo
    enddo
    do j=1,nv
       do i=1,nv
          sumx00=0.0d0
          do l=1,np
             sumx00 = sumx00 + deriv%M_t(l,j  )*deriv%vtempt(l,i  ,1)
	  enddo
	  v(i  ,j  ) = sumx00
       enddo
    enddo

endif

  end function interpolate_p2v







  function gradient_sphere(s,deriv,elem) result(ds)
!
!   input s:  scalar
!   output  ds: spherical gradient of s, lat-lon coordinates
!

    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind), intent(in) :: s(nv,nv)

    real(kind=real_kind) :: ds(nv,nv,2)

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dsdx00
    real(kind=real_kind) ::  dsdy00
    real(kind=real_kind) ::  rdx,rdy
    real(kind=real_kind) ::  v1(nv,nv),v2(nv,nv)

    rdx=2.0D0/(elem%dx*rearth) ! strong derivative inverse x length
    rdy=2.0D0/(elem%dy*rearth) ! strong derivative inverse y length

    do j=1,nv
       do l=1,nv
          dsdx00=0.0d0
          dsdy00=0.0d0
          do i=1,nv
             dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
             dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
          end do
          v1(l  ,j  ) = rdx*dsdx00
          v2(j  ,l  ) = rdy*dsdy00
       end do
    end do
    ! convert covarient to latlon
    do j=1,nv
       do i=1,nv
          ds(i,j,1)=elem%Dinv(1,1,i,j)*v1(i,j) + elem%Dinv(2,1,i,j)*v2(i,j) 
          ds(i,j,2)=elem%Dinv(1,2,i,j)*v1(i,j) + elem%Dinv(2,2,i,j)*v2(i,j) 
       enddo
    enddo
 
    end function gradient_sphere



  function curl_sphere(s,deriv,elem) result(ds)
!
!   input s:  scalar  (assumed to be  s khat)
!   output  curl(s khat) in lat-lon coordinates
! 
!   This subroutine can be used to compute divergence free velocity fields,
!   since div(ds)=0
!
!    first compute:  
!    curl(s khat) =  ( ds/dy, -ds/dx ) in contra-variant coordinates
!    then map to lat-lon
!
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind), intent(in) :: s(nv,nv)

    real(kind=real_kind) :: ds(nv,nv,2)

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dsdx00
    real(kind=real_kind) ::  dsdy00
    real(kind=real_kind) ::  rdx,rdy
    real(kind=real_kind) ::  v1(nv,nv),v2(nv,nv)

    rdx=2.0D0/(elem%dx*rearth) ! strong derivative inverse x length
    rdy=2.0D0/(elem%dy*rearth) ! strong derivative inverse y length

    do j=1,nv
       do l=1,nv
          dsdx00=0.0d0
          dsdy00=0.0d0
          do i=1,nv
             dsdx00 = dsdx00 + deriv%Dvv(i,l  )*s(i,j  )
             dsdy00 = dsdy00 + deriv%Dvv(i,l  )*s(j  ,i)
          end do
          v2(l  ,j  ) = -rdx*dsdx00
          v1(j  ,l  ) = rdy*dsdy00 
       end do
    end do
    ! convert contra -> latlon
    do j=1,nv
       do i=1,nv
          ds(i,j,1)=elem%rmetdetp(i,j)* (elem%D(1,1,i,j)*v1(i,j) + elem%D(1,2,i,j)*v2(i,j) )
          ds(i,j,2)=elem%rmetdetp(i,j)* (elem%D(2,1,i,j)*v1(i,j) + elem%D(2,2,i,j)*v2(i,j) ) 
       enddo
    enddo
 
    end function curl_sphere




  function divergence_sphere_wk(v,deriv,elem) result(div)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v, integrated by parts
!
!   Computes  -< grad(psi) dot v > 
!   (the integrated by parts version of < psi div(v) > )
!
!   note: after DSS, divergence_sphere() and divergence_sphere_wk() 
!   are identical to roundoff, as theory predicts.
!
    real(kind=real_kind), intent(in) :: v(nv,nv,2)  ! in lat-lon coordinates
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: div(nv,nv)

    ! Local

    integer i,j,m,n

    real(kind=real_kind) ::  vtemp(nv,nv,2)
    real(kind=real_kind) ::  ggtemp(nv,nv,2)
    real(kind=real_kind) ::  gtemp(nv,nv,2)
    real(kind=real_kind) ::  psi(nv,nv)
    real(kind=real_kind) :: rdx,rdy,xtmp

    rdx=2.0D0/(elem%dx*rearth) ! strong derivative inverse x length
    rdy=2.0D0/(elem%dy*rearth) ! strong derivative inverse y length


    ! latlon- > contra
    do j=1,nv
       do i=1,nv
          vtemp(i,j,1)=(elem%Dinv(1,1,i,j)*v(i,j,1) + elem%Dinv(1,2,i,j)*v(i,j,2))
          vtemp(i,j,2)=(elem%Dinv(2,1,i,j)*v(i,j,1) + elem%Dinv(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    do n=1,nv
       do m=1,nv

          div(m,n)=0
          do j=1,nv
             div(m,n)=div(m,n)-rdx*elem%spheremv(j,n)*vtemp(j,n,1)*deriv%Dvv(m,j) &
                              -rdy*elem%spheremv(m,j)*vtemp(m,j,2)*deriv%Dvv(n,j)
          enddo

#if 0
! debug the above formula using the N^4 slow formulation:
          psi=0
          psi(m,n)=1
          ggtemp=gradient_sphere(psi,deriv,elem)
          ! latlon -> covarient
          do j=1,nv
             do i=1,nv
                gtemp(i,j,1)=(elem%D(1,1,i,j)*ggtemp(i,j,1) + elem%D(2,1,i,j)*ggtemp(i,j,2))
                gtemp(i,j,2)=(elem%D(1,2,i,j)*ggtemp(i,j,1) + elem%D(2,2,i,j)*ggtemp(i,j,2))
             enddo
          enddo
! grad(psi) dot v:
          xtmp=0
          do j=1,nv
          do i=1,nv
             xtmp=xtmp-elem%spheremv(i,j)*(vtemp(i,j,1)*gtemp(i,j,1)+vtemp(i,j,2)*gtemp(i,j,2))
          enddo
          enddo
          if (abs(xtmp-div(m,n)) > 3e-17) then
             print *,m,n,xtmp,div(m,n),xtmp-div(m,n)
          endif
#endif          
       end do
    end do
    
  end function divergence_sphere_wk






    

  function vorticity_sphere(v,deriv,elem) result(vort)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  spherical vorticity of v
!

    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind), intent(in) :: v(nv,nv,2)

    real(kind=real_kind) :: vort(nv,nv)

    integer i
    integer j
    integer l
    
    real(kind=real_kind) ::  dvdx00
    real(kind=real_kind) ::  dudy00
    real(kind=real_kind) ::  vco(nv,nv,2)
    real(kind=real_kind) ::  vtemp(nv,nv)
    real(kind=real_kind) :: rdx
    real(kind=real_kind) :: rdy

    rdx=2.0D0/(elem%dx*rearth) ! strong derivative inverse x length
    rdy=2.0D0/(elem%dy*rearth) ! strong derivative inverse y length

    ! convert to covariant form
    do j=1,nv
       do i=1,nv
          vco(i,j,1)=(elem%D(1,1,i,j)*v(i,j,1) + elem%D(2,1,i,j)*v(i,j,2))
          vco(i,j,2)=(elem%D(1,2,i,j)*v(i,j,1) + elem%D(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    do j=1,nv
       do l=1,nv

          dudy00=0.0d0
	  dvdx00=0.0d0

          do i=1,nv
             dvdx00 = dvdx00 + deriv%Dvv(i,l  )*vco(i,j  ,2)
             dudy00 = dudy00 + deriv%Dvv(i,l  )*vco(j  ,i,1)
	  enddo
 
	  vort(l  ,j  ) = dvdx00
	  vtemp(j  ,l  ) = dudy00
	enddo
     enddo

    do j=1,nv
       do i=1,nv
          vort(i,j)=elem%rmetdetp(i,j)*(rdx*vort(i,j)-rdy*vtemp(i,j))
       end do
    end do

  end function vorticity_sphere



  function divergence_sphere(v,deriv,elem) result(div)
!
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v
!


    real(kind=real_kind), intent(in) :: v(nv,nv,2)  ! in lat-lon coordinates
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: div(nv,nv)

    ! Local

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dudx00
    real(kind=real_kind) ::  dvdy00
    real(kind=real_kind) ::  gv(nv,nv,2),vvtemp(nv,nv)
    real(kind=real_kind) :: rdx
    real(kind=real_kind) :: rdy

    rdx=2.0D0/(elem%dx*rearth) ! strong derivative inverse x length
    rdy=2.0D0/(elem%dy*rearth) ! strong derivative inverse y length

    ! convert to contra variant form and multiply by g
    do j=1,nv
       do i=1,nv
          gv(i,j,1)=elem%metdet(i,j)*(elem%Dinv(1,1,i,j)*v(i,j,1) + elem%Dinv(1,2,i,j)*v(i,j,2))
          gv(i,j,2)=elem%metdet(i,j)*(elem%Dinv(2,1,i,j)*v(i,j,1) + elem%Dinv(2,2,i,j)*v(i,j,2))
       enddo
    enddo

    ! compute d/dx and d/dy         
    do j=1,nv
       do l=1,nv
          dudx00=0.0d0
          dvdy00=0.0d0
          do i=1,nv
             dudx00 = dudx00 + deriv%Dvv(i,l  )*gv(i,j  ,1)
             dvdy00 = dvdy00 + deriv%Dvv(i,l  )*gv(j  ,i,2)
          end do
          div(l  ,j  ) = dudx00
          vvtemp(j  ,l  ) = dvdy00
       end do
    end do

    do j=1,nv
       do i=1,nv
          div(i,j)=elem%rmetdetp(i,j)*(rdx*div(i,j)+rdy*vvtemp(i,j))
       end do
    end do
    
  end function divergence_sphere


  function laplace_sphere_wk(s,deriv,elem) result(laplace)
!
!   input:  s = scalar
!   ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
!     note: for this form of the operator, grad(s) does not need to be made C0
!            
    real(kind=real_kind), intent(in) :: s(nv,nv) 
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: laplace(nv,nv)

    ! Local
    real(kind=real_kind) :: grads(nv,nv,2)

    grads=gradient_sphere(s,deriv,elem)
    ! note: divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
    ! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().  
    laplace=divergence_sphere_wk(grads,deriv,elem)

  end function laplace_sphere_wk




  function vlaplace_sphere_wk(v,deriv,elem) result(laplace)
!
!   input:  v = vector in lat-lon coordinates
!   ouput:  weak laplacian of v, in lat-lon coordinates
!   note: integrals must be performed in contra variant coordinates,
!         convert to lat-lon after computing integral
!
!   laplace(v) =  grad(div) -  curl(vor) 
!   weak form   < PHI , grad(div) > - < PHI, curl(vor*khat) >    
!   by parts:   -< div(PHI) div >   - < vor curl(PHI) >      
!             
!
! used vector identity: div(F cross G)=G dot curl F - F dot curl G
!                  OR:    < G, curl F> = < F, curl G >
!
!  NOTE: for the equation:   < PHI, LAPLACE > =   -< div(PHI) div >   - < vor curl(PHI) >      
!  if LAPLACE is in covarient, we test with PHI = (1,0) and (0,1) in contra-variant
!  if LAPLACE is in contra-varient, we test with PHI = (1,0) and (0,1) in co-varient
!
!  compute with two different test functions:
!  (then transform back to lat-lon (since we output in lat-lon))
!  test function 1:  contra:  (phi,0)        covarient:   (met11 phi, met21 phi)
!  test function 2:  contra:  (0,phi)        covariant:   (met12 phi, met22 phi)
!  
!  div acts on contra components
!   < div(PHI) div >  =  <  1/g (g phi)_x div >  = <  phi_x div >    (test 1)
!                        <  1/g (g phi)_y div >  = <  phi_y div >    (test 2)
!
!
!  curl  acts on co-variant
!   < curl(PHI) vor >  =   <  1/g [-(met11 phi )_y + (met21 phi)_x  ] vor > 
!                          <  1/g [-(met12 phi )_y + (met22 phi)_x  ] vor >
!                      =   <  1/g [-met11 phi_y + met21 phi_x  ] vor > 
!                          <  1/g [-met12 phi_y + met22 phi_x  ] vor >
!             
!
!  curl acting on co-variant test functions:
!  test function 1:  co:  (phi,0)        curl(PHI) = -phi_y
!  test function 2:  co:  (0,phi)        curl(PHI) = phi_x
!   < curl(PHI) vor >  =  <  1/g -phi_y vor  >  = <  1/g -phy_y vor  >    (test 1)
!                         <  1/g phi_x vor >    = <  1/g phi_x vor >      (test 2)
!
!                         
!   
!
! NOTE:  SEM can compute (g11 phi)_y in two ways:
!        project, than derivative:     met11 phi_y   (we are using this formula)
!        expand first:                 met11 phi_y  + met11_y phi
!
!
! compare to weak divergence:  < grad(PHI), v >
!  if v is in contra,  grad(PHI) in covariant = (phi_x, phi_y)
!  < phi_x v1  + phi_y v2 >   THUS:  < phi_x div > loop should look like d/dx
!  loop in divergence_sphere_wk()
!
! NOTE: dont forget < u v > = integral g*u*v = sum mv()*metdet()*u*v  
!        with g=metdet(), spheremv=mv()*metdet()
!
!  < phy_y div > = sum spheremv * phy_y * div
!  < 1/g F >     = sum mv * F
!


    real(kind=real_kind), intent(in) :: v(nv,nv,2) 
    type (derivative_t)              :: deriv
    type (element_t)                 :: elem
    real(kind=real_kind) :: laplace(nv,nv,2)

    ! Local

    integer l,m,n
    real(kind=real_kind) :: vor(nv,nv),div(nv,nv)
    real(kind=real_kind) :: v1,v2,div1,div2,vor1,vor2,phi_x,phi_y,rdx,rdy

    rdx=2.0D0/(elem%dx*rearth) ! cube face -> ref element metric term
    rdy=2.0D0/(elem%dy*rearth) ! cube face -> ref element metric term

    div=divergence_sphere(v,deriv,elem)
    vor=vorticity_sphere(v,deriv,elem)


#undef DIVCONTRA_VORCO

    do n=1,nv
       do m=1,nv

          div1=0; div2=0;
          vor1=0; vor2=0; 

          do l=1,nv
             phi_x=deriv%Dvv(m,l)*rdx   ! (m,n) cardinal function, d/dx  
                                        ! phi_x(i,j) = 0  for j<>n.  so treat this as phi_x(i,n)

             phi_y=deriv%Dvv(n,l)*rdy   ! (m,n) cardinal function, d/dy
                                        ! phi_y(i,j) = 0  for i<>m.  so treat this as phi_x(m,j)

             div1=div1 + elem%spheremv(l,n)*div(l,n)*phi_x
             div2=div2 + elem%spheremv(m,l)*div(m,l)*phi_y
             
#ifdef DIVCONTRA_VORCO
             vor1=vor1 - elem%mv(m,l)*vor(m,l)*phi_y
             vor2=vor2 + elem%mv(l,n)*vor(l,n)*phi_x
#else
             vor1=vor1 - elem%mv(m,l)*vor(m,l)*elem%met(1,1,m,l)*phi_y &
                         + elem%mv(l,n)*vor(l,n)*elem%met(2,1,l,n)*phi_x

             vor2=vor2 - elem%mv(m,l)*vor(m,l)*elem%met(1,2,m,l)*phi_y &
                         + elem%mv(l,n)*vor(l,n)*elem%met(2,2,l,n)*phi_x

#endif

          enddo
#ifdef DIVCONTRA_VORCO
          v1=-div1
          v2=-div2

          !  (v1,v2) = divergence componet tested against contra-variant, so result is CO-variant
          laplace(m,n,1)=elem%Dinv(1,1,m,n)*v1 + elem%Dinv(2,1,m,n)*v2   ! co->latlon
          laplace(m,n,2)=elem%Dinv(1,2,m,n)*v1 + elem%Dinv(2,2,m,n)*v2   ! co->latlon

          v1=-vor1
          v2=-vor2
          !  (v1,v2) = vorticity component tested against co-variant.  so result is CONTRA 
          laplace(m,n,1)=laplace(m,n,1) + elem%D(1,1,m,n)*v1 + elem%D(1,2,m,n)*v2   ! contra->latlon
          laplace(m,n,2)=laplace(m,n,2) + elem%D(2,1,m,n)*v1 + elem%D(2,2,m,n)*v2   ! contra->latlon
#else
          v1=-( div1 + vor1 )
          v2=-( div2 + vor2 )

          !  (v1,v2) = RHS tested agains contra-variant delta functions, so result is CO-varient
          laplace(m,n,1)=elem%Dinv(1,1,m,n)*v1 + elem%Dinv(2,1,m,n)*v2   ! co->latlon
          laplace(m,n,2)=elem%Dinv(1,2,m,n)*v1 + elem%Dinv(2,2,m,n)*v2   ! co->latlon
#endif
          ! add in correction so we dont damp rigid rotation
#define UNDAMPRR
#ifdef UNDAMPRR
          laplace(m,n,1)=laplace(m,n,1) + 2*elem%spheremv(m,n)*v(m,n,1)/(rearth**2)
          laplace(m,n,2)=laplace(m,n,2) + 2*elem%spheremv(m,n)*v(m,n,2)/(rearth**2)
#endif
       enddo
    enddo
  end function vlaplace_sphere_wk



end module derivative_mod








