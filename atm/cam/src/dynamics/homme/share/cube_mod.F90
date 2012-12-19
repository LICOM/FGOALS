#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _BEGIN_FACE 1
#define _END_FACE   4
#undef _FACE_6
#undef _FACE_5

module cube_mod
  use kinds, only : real_kind
  use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t
  use physical_constants, only : dd_pi
  implicit none
  private

  integer,public, parameter :: nfaces = 6          ! number of faces on the cube
  integer,public, parameter :: nInnerElemEdge = 8  ! number of edges for an interior element
  integer,public, parameter :: nCornerElemEdge = 4 ! number of corner elements

  type, public :: face_t
     sequence
     type (spherical_polar_t) :: sphere0       ! tangent point of face on sphere
     type (spherical_polar_t) :: sw            ! sw corner of face on sphere
     type (spherical_polar_t) :: se            ! se corner of face on sphere
     type (spherical_polar_t) :: ne            ! ne corner of face on sphere
     type (spherical_polar_t) :: nw            ! nw corner of face on sphere
     type (cartesian3D_t)     :: P0
     type (cartesian3D_t)     :: X0
     type (cartesian3D_t)     :: Y0
     integer                  :: number
     integer                  :: padding       ! padd the struct
  end type face_t

  type, public :: cube_face_coord_t
     sequence
     real(real_kind) :: x             ! x coordinate
     real(real_kind) :: y             ! y coordinate
     type (face_t), pointer :: face     ! face
  end type cube_face_coord_t

  ! ==========================================
  ! Public Interfaces
  ! ==========================================

  public :: CubeTopology

  ! ===============================================
  ! problem domain size: cube face (equal angular)
  ! ===============================================

  real (kind=real_kind), public, parameter :: cube_xstart = -0.25D0*DD_PI
  real (kind=real_kind), public, parameter :: cube_xend   =  0.25D0*DD_PI

  real (kind=real_kind), public, parameter :: cube_ystart = -0.25D0*DD_PI
  real (kind=real_kind), public, parameter :: cube_yend   =  0.25D0*DD_PI

  ! Rotate the North Pole:  used for JW baroclinic test case
  ! Settings this only changes Coriolis.  
  ! User must also rotate initial condition
  real (kind=real_kind), public :: rotate_grid = 0

  integer, private, parameter :: nface = 6

  ! ===============================
  ! Public methods for cube
  ! ===============================

  public  :: cube_init,cube_init_atomic
  public  :: convert_gbl_index
  public  :: cube_assemble
  public  :: vmap
  public  :: covariant_rot
  public  :: contravariant_rot

  public  :: CubeEdgeCount
  public  :: CubeElemCount
  public  :: CubeSetupEdgeIndex
  public  :: rotation_init_atomic

  ! ===============================
  ! Private methods
  ! ===============================

  private :: coordinates,coordinates_atomic
  private :: metric,metric_atomic
  private :: coreolis_init, coreolis_init_atomic
  private :: solver_weights
  private :: GetLatticeSpacing


contains

  ! =======================================
  ! cube_init:
  !
  ! Initialize element descriptors for
  ! cube sphere case...
  ! =======================================

  subroutine cube_init(elem)
    use element_mod, only : element_t

    type (element_t) :: elem(:)

    call coordinates(elem)

    call metric(elem)

    call coreolis_init(elem)

    elem(:)%desc%use_rotation=0

    call solver_weights(elem)

  end subroutine cube_init

  ! =======================================
  !  cube_init_atomic:
  !
  ! Initialize element descriptors for 
  ! cube sphere case for each element ... 
  ! =======================================
  subroutine cube_init_atomic(elem,alpha_in)
    use element_mod, only : element_t
    type (element_t) :: elem
    real (kind=real_kind),optional :: alpha_in
    real (kind=real_kind)          :: alpha=1

    if(present(alpha_in)) alpha=alpha_in
    
    call coordinates_atomic(elem)

    call metric_atomic(elem,alpha)

    call coreolis_init_atomic(elem)
    elem%desc%use_rotation= 0
    call solver_weights_atomic(elem)

  end subroutine cube_init_atomic
  ! =======================================
  ! coordinates:
  !
  ! Initialize element coordinates for
  ! cube-sphere case. 
  !
  ! =======================================

  subroutine coordinates(elem)
    use element_mod, only : element_t, element_coordinates
    use coordinate_systems_mod, only : cartesian2d_t, project
    use dimensions_mod, only : nv, ne, np
    use quadrature_mod, only : quadrature_t, gauss, gausslobatto

    type (element_t) :: elem(:)

    ! Local variables

    integer ii
    integer ie,je,face_no
    real (kind=real_kind)  :: dx,dy
    real (kind=real_kind)  :: vspacing
    real (kind=real_kind)  :: pspacing
    type (cartesian2D_t)   :: start,end
    type (cartesian2D_t)   :: tanv(nv,nv)
    type (cartesian2D_t)   :: tanp(np,np)
    type (quadrature_t)    :: gv
    type (quadrature_t)    :: gp

    dx = (cube_xend-cube_xstart)/ne
    dy = (cube_yend-cube_ystart)/ne

    gv=gausslobatto(nv)
    if (np==nv) then
       gp=gausslobatto(nv)
    else if (np==nv-2) then
       gp=gauss(np)
    endif

    do ii=1,SIZE(elem)

       ! =========================================
       ! compute cube face coordinates of element
       ! =========================================

       call convert_gbl_index(elem(ii)%vertex%number,ie,je,face_no)

       start%x=cube_xstart+ie*dx
       start%y=cube_ystart+je*dy

       end%x  =start%x + dx
       end%y  =start%y + dy

       elem(ii)%dx=dx
       elem(ii)%dy=dy

       elem(ii)%cartv=element_coordinates(start,end,gv%points)
       elem(ii)%cartp=element_coordinates(start,end,gp%points)

       tanv(:,:)%x = TAN(elem(ii)%cartv(:,:)%x)
       tanv(:,:)%y = TAN(elem(ii)%cartv(:,:)%y)
       call project(elem(ii)%spherev,tanv,face_no)

       tanp(:,:)%x = TAN(elem(ii)%cartp(:,:)%x)
       tanp(:,:)%y = TAN(elem(ii)%cartp(:,:)%y)
       call project(elem(ii)%spherep,tanp,face_no)

       if (ii == 1) then
          call GetLatticeSpacing(elem(ii)%spherev,vspacing,elem(ii)%spherep,pspacing)       
       end if

    end do

  end subroutine coordinates
  ! =======================================
  ! coordinates_atomic:
  !
  ! Initialize element coordinates for
  ! cube-sphere case ... (atomic) 
  !
  ! =======================================

  subroutine coordinates_atomic(elem)
    use element_mod, only : element_t, element_coordinates
    use coordinate_systems_mod, only : cartesian2d_t, project
    use dimensions_mod, only : nv, ne, np
    use quadrature_mod, only : quadrature_t, gauss, gausslobatto
    
    type (element_t) :: elem

    ! Local variables

    integer ii
    integer ie,je,face_no
    real (kind=real_kind)  :: dx,dy
    real (kind=real_kind)  :: vspacing
    real (kind=real_kind)  :: pspacing
    type (cartesian2D_t)   :: start,end
    type (cartesian2D_t)   :: tanv(nv,nv)
    type (cartesian2D_t)   :: tanp(np,np)
    type (quadrature_t)    :: gv
    type (quadrature_t)    :: gp

    dx = (cube_xend-cube_xstart)/ne
    dy = (cube_yend-cube_ystart)/ne

    gv=gausslobatto(nv)
    if (np==nv) then
       gp=gausslobatto(nv)
    else if (np==nv-2) then
       gp=gauss(np)
    endif


    ! =========================================
    ! compute cube face coordinates of element
    ! =========================================

    call convert_gbl_index(elem%vertex%number,ie,je,face_no)

    start%x=cube_xstart+ie*dx
    start%y=cube_ystart+je*dy

    end%x  =start%x + dx
    end%y  =start%y + dy

    elem%dx=dx
    elem%dy=dy

    elem%cartv=element_coordinates(start,end,gv%points)
    elem%cartp=element_coordinates(start,end,gp%points)

    tanv(:,:)%x = TAN(elem%cartv(:,:)%x)
    tanv(:,:)%y = TAN(elem%cartv(:,:)%y)
    call project(elem%spherev,tanv,face_no)

    tanp(:,:)%x = TAN(elem%cartp(:,:)%x)
    tanp(:,:)%y = TAN(elem%cartp(:,:)%y)
    call project(elem%spherep,tanp,face_no)


  end subroutine coordinates_atomic

  ! =========================================
  ! metric:
  !
  ! Initialize cube-sphere metric terms:
  ! equal angular elements
  ! =========================================

  subroutine metric(elem)
    use element_mod, only : element_t
    use dimensions_mod, only : np, nv
    type (element_t) :: elem(:)

    ! Local variables

    integer ii,ie,je,face_no
    integer i,j
    integer iptr

    real (kind=real_kind) :: r         ! distance from origin for point on cube tangent to unit sphere

    real (kind=real_kind) :: const     
    real (kind=real_kind) :: detD      ! determinant of vector field mapping matrix.  

    real (kind=real_kind) :: x1        ! 1st cube face coordinate
    real (kind=real_kind) :: x2        ! 2nd cube face coordinate

    do ii=1,SIZE(elem)

       call convert_gbl_index(elem(ii)%vertex%number,ie,je,face_no)

       ! ==============================================
       ! Initialize rgdet=1/gdet on pressure grid.
       !       gdet = SQRT(ABS(DET(gij)))
       !
       ! ==============================================

       iptr=1
       do j=1,np
          do i=1,np
             x1=elem(ii)%cartp(i,j)%x
             x2=elem(ii)%cartp(i,j)%y
             r=SQRT(1.0D0 + TAN(x1)**2 + TAN(x2)**2)        
             elem(ii)%rmetdetp(i,j)  = (r**3)*(COS(x1)**2)*(COS(x2)**2)
             !             elem(ii)%metdetp(i,j)   = 1.0D0/((r**3)*(COS(x1)**2)*(COS(x2)**2))
             elem(ii)%metdetp(i,j)   = 1.0D0/elem(ii)%rmetdetp(i,j)
             iptr=iptr+1
          end do
       end do

       ! ===============================================
       !
       ! Initialize equal angular metric tensor on each 
       ! on velocity grid for unit sphere.
       !
       ! Initialize gdet = SQRT(ABS(DET(gij)))
       !
       ! These quantities are the same on every face
       ! of the cube.
       !
       ! =================================================

       do j=1,nv
          do i=1,nv
             x1=elem(ii)%cartv(i,j)%x
             x2=elem(ii)%cartv(i,j)%y
             r=SQRT(1.0D0 + TAN(x1)**2 + TAN(x2)**2)        

             ! ==============================
             ! Metric term: g_{ij}
             ! ==============================

             const=1.0D0/((r**4)*(COS(x1)**2)*(COS(x2)**2))

             elem(ii)%met(1,1,i,j) =   const*(1.0D0 + TAN(x1)**2)
             elem(ii)%met(1,2,i,j) = - const*TAN(x1)*TAN(x2)
             elem(ii)%met(2,1,i,j) = - const*TAN(x1)*TAN(x2)
             elem(ii)%met(2,2,i,j) =   const*(1.0D0 + TAN(x2)**2)

             ! ==============================
             ! metdet = SQRT(DET(g_{ij}))
             ! ==============================

             elem(ii)%metdet(i,j)= 1.0D0/((r**3)*(COS(x1)**2)*(COS(x2)**2))

             ! ==============================
             ! Inverse of metric term: g^{ij}
             ! ==============================

             const=(r**2)*(COS(x1)**2)*(COS(x2)**2)

             elem(ii)%metinv(1,1,i,j) =   const*(1.0D0 + TAN(x2)**2)
             elem(ii)%metinv(1,2,i,j) =   const*TAN(x1)*TAN(x2)
             elem(ii)%metinv(2,1,i,j) =   const*TAN(x1)*TAN(x2)
             elem(ii)%metinv(2,2,i,j) =   const*(1.0D0 + TAN(x1)**2)

          end do
       end do

       ! ==============================================
       ! Initialize differential mapping operator
       ! to and from vector fields on the sphere to 
       ! contravariant vector fields on the cube
       ! i.e. dM/dx^i in Sadourney (1972) and it's 
       ! inverse
       ! ==============================================

       do j=1,nv
          do i=1,nv
             x1=elem(ii)%cartv(i,j)%x
             x2=elem(ii)%cartv(i,j)%y

             call vmap(elem(ii)%D(:,:,i,j),x1,x2,face_no)

             ! compute D^-1...
             ! compute determinant of D mapping matrix... if not zero compute inverse

             detD = elem(ii)%D(1,1,i,j)*elem(ii)%D(2,2,i,j) - elem(ii)%D(1,2,i,j)*elem(ii)%D(2,1,i,j)      

             elem(ii)%Dinv(1,1,i,j) =  elem(ii)%D(2,2,i,j)/detD
             elem(ii)%Dinv(1,2,i,j) = -elem(ii)%D(1,2,i,j)/detD
             elem(ii)%Dinv(2,1,i,j) = -elem(ii)%D(2,1,i,j)/detD
             elem(ii)%Dinv(2,2,i,j) =  elem(ii)%D(1,1,i,j)/detD

          end do
       end do

    end do

  end subroutine metric
  ! =========================================
  ! metric_atomic:
  !
  ! Initialize cube-sphere metric terms:
  ! equal angular elements (atomic)
  ! initialize:  
  !         metdetp, rmetdetp  (analytic)    = detD, 1/detD
  !         met                (analytic)    D'D or DD' ?                   
  !         metdet             (analytic)    = detD
  !         metinv             (analytic)    Dinv'Dinv  or Dinv Dinv' ?     
  !         D     (from subroutine vmap)
  !         Dinv  (computed directly from D)
  ! 
  ! so if we want to tweak the mapping by a factor alpha (so he weights add up to 4pi, for example)
  ! we take:
  !    NEW       OLD     
  !       D = sqrt(alpha) D  and then rederive all quantities.  
  !    detD = alpha detD
  !    
  ! where alpha = 4pi/SEMarea, SEMarea = global sum elem(ie)%mv(i,j)*elem(ie)%metdet(i,j)
  ! 
  ! =========================================

  subroutine metric_atomic(elem,alpha)
    use element_mod, only : element_t
    use dimensions_mod, only : np, nv

    type (element_t) :: elem
    real(kind=real_kind) :: alpha

    ! Local variables

    integer ii,ie,je,face_no
    integer i,j
    integer iptr

    real (kind=real_kind) :: r         ! distance from origin for point on cube tangent to unit sphere

    real (kind=real_kind) :: const     
    real (kind=real_kind) :: detD      ! determinant of vector field mapping matrix.  

    real (kind=real_kind) :: x1        ! 1st cube face coordinate
    real (kind=real_kind) :: x2        ! 2nd cube face coordinate


    call convert_gbl_index(elem%vertex%number,ie,je,face_no)

    ! ==============================================
    ! Initialize rgdet=1/gdet on pressure grid.
    !       gdet = SQRT(ABS(DET(gij)))
    !
    ! ==============================================

    do j=1,np
       do i=1,np
          x1=elem%cartp(i,j)%x
          x2=elem%cartp(i,j)%y
          r=SQRT(1.0D0 + TAN(x1)**2 + TAN(x2)**2)        
          elem%rmetdetp(i,j)  = (r**3)*(COS(x1)**2)*(COS(x2)**2)
          !             elem%metdetp(i,j) = 1.0D0/((r**3)*(COS(x1)**2)*(COS(x2)**2))
          elem%metdetp(i,j) = 1.0D0/elem%rmetdetp(i,j)
       end do
    end do

    ! ===============================================
    !
    ! Initialize equal angular metric tensor on each 
    ! on velocity grid for unit sphere.
    !
    ! Initialize gdet = SQRT(ABS(DET(gij)))
    !
    ! These quantities are the same on every face
    ! of the cube.
    !
    ! =================================================

    do j=1,nv
       do i=1,nv
          x1=elem%cartv(i,j)%x
          x2=elem%cartv(i,j)%y
          r=SQRT(1.0D0 + TAN(x1)**2 + TAN(x2)**2)        

          ! ==============================
          ! Metric term: g_{ij}
          ! ==============================

          const=1.0D0/((r**4)*(COS(x1)**2)*(COS(x2)**2))

          elem%met(1,1,i,j) =   const*(1.0D0 + TAN(x1)**2)
          elem%met(1,2,i,j) = - const*TAN(x1)*TAN(x2)
          elem%met(2,1,i,j) = - const*TAN(x1)*TAN(x2)
          elem%met(2,2,i,j) =   const*(1.0D0 + TAN(x2)**2)

          ! ==============================
          ! metdet = SQRT(DET(g_{ij}))
          ! ==============================

          elem%metdet(i,j)= 1.0D0/((r**3)*(COS(x1)**2)*(COS(x2)**2))

          ! ==============================
          ! Inverse of metric term: g^{ij}
          ! ==============================

          const=(r**2)*(COS(x1)**2)*(COS(x2)**2)

          elem%metinv(1,1,i,j) =   const*(1.0D0 + TAN(x2)**2)
          elem%metinv(1,2,i,j) =   const*TAN(x1)*TAN(x2)
          elem%metinv(2,1,i,j) =   const*TAN(x1)*TAN(x2)
          elem%metinv(2,2,i,j) =   const*(1.0D0 + TAN(x1)**2)

       end do
    end do

    ! ==============================================
    ! Initialize differential mapping operator
    ! to and from vector fields on the sphere to 
    ! contravariant vector fields on the cube
    ! i.e. dM/dx^i in Sadourney (1972) and it's 
    ! inverse
    ! ==============================================

    do j=1,nv
       do i=1,nv
          x1=elem%cartv(i,j)%x
          x2=elem%cartv(i,j)%y

          call vmap(elem%D(1,1,i,j),x1,x2,face_no)

          ! compute D^-1...
          ! compute determinant of D mapping matrix... if not zero compute inverse

          detD = elem%D(1,1,i,j)*elem%D(2,2,i,j) - elem%D(1,2,i,j)*elem%D(2,1,i,j)      

          elem%Dinv(1,1,i,j) =  elem%D(2,2,i,j)/detD
          elem%Dinv(1,2,i,j) = -elem%D(1,2,i,j)/detD
          elem%Dinv(2,1,i,j) = -elem%D(2,1,i,j)/detD
          elem%Dinv(2,2,i,j) =  elem%D(1,1,i,j)/detD

       end do
    end do

    ! mt: better might be to compute all these quantities directly from D
    ! for consistency?
    elem%D = elem%D * sqrt(alpha) 
    elem%Dinv = elem%Dinv / sqrt(alpha) 
    elem%metdetp = elem%metdetp * alpha
    elem%metdet = elem%metdet * alpha
    elem%rmetdetp = elem%rmetdetp / alpha
    elem%met = elem%met * alpha
    elem%metinv = elem%metinv / alpha

  end subroutine metric_atomic

  ! =======================================
  ! solver_weights:
  !
  ! For nonstaggered GaussLobatto elements,
  ! compute weights for redundant points in 
  ! cg solver.
  !
  ! =======================================

  subroutine solver_weights(elem)
    use element_mod, only : element_t
    use dimensions_mod, only : nv, ne

    type (element_t) :: elem(:)

    ! Local variables

    integer ii,i,j
    integer ie,je,face_no

    do ii=1,SIZE(elem)

       ! =========================================
       ! compute cube face coordinates of element
       ! =========================================

       call convert_gbl_index(elem(ii)%vertex%number,ie,je,face_no)

       if ((ie == 0) .and. (je==0)) then
          elem(ii)%solver_wts(1,1) = 1.0_real_kind/3.0_real_kind
       else
          elem(ii)%solver_wts(1,1) = 0.25_real_kind
       end if

       if ((ie == 0) .and. (je==ne-1)) then
          elem(ii)%solver_wts(1,nv) = 1.0_real_kind/3.0_real_kind
       else
          elem(ii)%solver_wts(1,nv) = 0.25_real_kind
       end if

       if ((ie == ne-1) .and. (je==0)) then
          elem(ii)%solver_wts(nv,1) = 1.0_real_kind/3.0_real_kind
       else
          elem(ii)%solver_wts(nv,1) = 0.25_real_kind
       end if

       if ((ie == ne-1) .and. (je==ne-1)) then
          elem(ii)%solver_wts(nv,nv) = 1.0_real_kind/3.0_real_kind
       else
          elem(ii)%solver_wts(nv,nv) = 0.25_real_kind
       end if

       do i=2,nv-1
          elem(ii)%solver_wts(i,1) = 0.5_real_kind
          elem(ii)%solver_wts(i,nv)= 0.5_real_kind
       end do

       do j=2,nv-1
          elem(ii)%solver_wts(1,j) = 0.5_real_kind
          elem(ii)%solver_wts(nv,j) = 0.5_real_kind
          do i=2,nv-1
             elem(ii)%solver_wts(i,j) = 1.0_real_kind
          end do
       end do
    end do
  end subroutine solver_weights
  ! =======================================
  ! solver_weights:
  !
  ! For nonstaggered GaussLobatto elements,
  ! compute weights for redundant points in 
  ! cg solver.
  !
  ! =======================================

  subroutine solver_weights_atomic(elem)
    use element_mod, only : element_t
    use dimensions_mod, only : nv, ne

    type (element_t) :: elem

    ! Local variables

    integer :: i, j, ie,je,face_no
    ! =========================================
    ! compute cube face coordinates of element
    ! =========================================

    call convert_gbl_index(elem%vertex%number,ie,je,face_no)

    if ((ie == 0) .and. (je==0)) then
       elem%solver_wts(1,1) = 1.0_real_kind/3.0_real_kind
    else
       elem%solver_wts(1,1) = 0.25_real_kind
    end if

    if ((ie == 0) .and. (je==ne-1)) then
       elem%solver_wts(1,nv) = 1.0_real_kind/3.0_real_kind
    else
       elem%solver_wts(1,nv) = 0.25_real_kind
    end if

    if ((ie == ne-1) .and. (je==0)) then
       elem%solver_wts(nv,1) = 1.0_real_kind/3.0_real_kind
    else
       elem%solver_wts(nv,1) = 0.25_real_kind
    end if

    if ((ie == ne-1) .and. (je==ne-1)) then
       elem%solver_wts(nv,nv) = 1.0_real_kind/3.0_real_kind
    else
       elem%solver_wts(nv,nv) = 0.25_real_kind
    end if

    do i=2,nv-1
       elem%solver_wts(i,1) = 0.5_real_kind
       elem%solver_wts(i,nv)= 0.5_real_kind
    end do

    do j=2,nv-1
       elem%solver_wts(1,j) = 0.5_real_kind
       elem%solver_wts(nv,j) = 0.5_real_kind
       do i=2,nv-1
          elem%solver_wts(i,j) = 1.0_real_kind
       end do
    end do
  end subroutine solver_weights_atomic

#if 1
  ! ========================================
  ! covariant_rot:
  !
  ! 2 x 2 matrix multiply:  Db^T * Da^-T
  ! for edge rotations: maps face a to face b
  !
  ! ========================================

  function covariant_rot(Da,Db) result(R)

    real (kind=real_kind) :: Da(2,2)
    real (kind=real_kind) :: Db(2,2)
    real (kind=real_kind) :: R(2,2)

    real (kind=real_kind) :: detDa

    detDa = Da(2,2)*Da(1,1) - Da(1,2)*Da(2,1)

    R(1,1)=(Da(2,2)*Db(1,1) - Da(1,2)*Db(2,1))/detDa
    R(1,2)=(Da(1,1)*Db(2,1) - Da(2,1)*Db(1,1))/detDa
    R(2,1)=(Da(2,2)*Db(1,2) - Da(1,2)*Db(2,2))/detDa
    R(2,2)=(Da(1,1)*Db(2,2) - Da(2,1)*Db(1,2))/detDa

  end function covariant_rot
#else

  ! ========================================
  ! covariant_rot:
  !
  ! 2 x 2 matrix multiply:  Db * Da^-1
  ! for edge rotations: maps face a to face b
  !
  ! ========================================

  function covariant_rot(Da,Db) result(R)

    real (kind=real_kind) :: Da(2,2)
    real (kind=real_kind) :: Db(2,2)
    real (kind=real_kind) :: R(2,2)

    real (kind=real_kind) :: detDa

    detDa = Da(2,2)*Da(1,1) - Da(1,2)*Da(2,1)

    R(1,1)=(Da(2,2)*Db(1,1) - Da(2,1)*Db(1,2))/detDa
    R(1,2)=(Da(1,1)*Db(1,2) - Da(1,2)*Db(1,1))/detDa
    R(2,1)=(Da(2,2)*Db(2,1) - Da(2,1)*Db(2,2))/detDa
    R(2,2)=(Da(1,1)*Db(2,2) - Da(1,2)*Db(2,1))/detDa

  end function covariant_rot

#endif

  ! ========================================
  ! contravariant_rot:
  !
  ! 2 x 2 matrix multiply:  Db^-1 * Da
  ! that maps a contravariant vector field
  ! from an edge of cube face a to a contiguous 
  ! edge of cube face b.
  !
  ! ========================================

  function contravariant_rot(Da,Db) result(R)

    real (kind=real_kind) :: Da(2,2)
    real (kind=real_kind) :: Db(2,2)
    real (kind=real_kind) :: R(2,2)

    real (kind=real_kind) :: detDb

    detDb = Db(2,2)*Db(1,1) - Db(1,2)*Db(2,1)

    R(1,1)=(Da(1,1)*Db(2,2) - Da(2,1)*Db(1,2))/detDb
    R(1,2)=(Da(1,2)*Db(2,2) - Da(2,2)*Db(1,2))/detDb
    R(2,1)=(Da(2,1)*Db(1,1) - Da(1,1)*Db(2,1))/detDb
    R(2,2)=(Da(2,2)*Db(1,1) - Da(1,2)*Db(2,1))/detDb

  end function contravariant_rot

  ! ========================================================
  ! vmap:
  !
  ! Initialize mapping that tranforms contravariant 
  ! vector fields on the cube onto vector fields on
  ! the sphere. This follows Taylor's D matrix 
  !
  !       | cos(theta)dlambda/dx1  cos(theta)dlambda/dx2 |
  !   D = |                                              |
  !       |     dtheta/dx1              dtheta/dx2       |
  !
  ! ========================================================

  subroutine vmap(D, x1, x2, face_no) 
    use coordinate_systems_mod, only : dist_threshold
    real (kind=real_kind), intent(inout)  :: D(2,2)
    real (kind=real_kind), intent(in)     :: x1
    real (kind=real_kind), intent(in)     :: x2
    integer              , intent(in)     :: face_no

    ! Local variables

    real (kind=real_kind) :: poledist  ! SQRT(TAN(x1)**2 +TAN(x2)**2)
    real (kind=real_kind) :: r         ! distance from cube point to center of sphere

    real (kind=real_kind) :: D11
    real (kind=real_kind) :: D12
    real (kind=real_kind) :: D21
    real (kind=real_kind) :: D22

    r=SQRT(1.0D0 + TAN(x1)**2 + TAN(x2)**2)

    if (face_no >= 1 .and. face_no <= 4) then

       D11 = 1.0D0/(r*COS(x1))
       D12 = 0.0D0
       D21 = -TAN(x1)*TAN(x2)/(COS(x1)*r*r)        
       D22 = 1.0D0/(r*r*COS(x1)*COS(x2)*COS(x2))

       D(1,1) =  D11
       D(1,2) =  D12
       D(2,1) =  D21
       D(2,2) =  D22


    else if (face_no==6) then
       poledist=SQRT( TAN(x1)**2 + TAN(x2)**2)
       if ( poledist <= DIST_THRESHOLD ) then

          ! we set the D transform to the identity matrix 
          ! which works ONLY for swtc1, phi starting at 
          ! 3*PI/2... assumes lon at pole == 0

          D(1,1) =  1.0D0
          D(1,2) =  0.0D0
          D(2,1) =  0.0D0
          D(2,2) =  1.0D0

       else

          D11 = -TAN(x2)/(poledist*COS(x1)*COS(x1)*r)
          D12 =  TAN(x1)/(poledist*COS(x2)*COS(x2)*r)
          D21 = -TAN(x1)/(poledist*COS(x1)*COS(x1)*r*r)
          D22 = -TAN(x2)/(poledist*COS(x2)*COS(x2)*r*r)

          D(1,1) =  D11
          D(1,2) =  D12
          D(2,1) =  D21
          D(2,2) =  D22

       end if
    else if (face_no==5) then
       poledist=SQRT( TAN(x1)**2 + TAN(x2)**2)
       if ( poledist <= DIST_THRESHOLD ) then

          ! we set the D transform to the identity matrix 
          ! which works ONLY for swtc1, phi starting at 
          ! 3*PI/2... assumes lon at pole == 0, i.e. very specific

          D(1,1) =  1.0D0
          D(1,2) =  0.0D0
          D(2,1) =  0.0D0
          D(2,2) =  1.0D0

       else

          D11 =  TAN(x2)/(poledist*COS(x1)*COS(x1)*r)
          D12 = -TAN(x1)/(poledist*COS(x2)*COS(x2)*r)
          D21 =  TAN(x1)/(poledist*COS(x1)*COS(x1)*r*r)
          D22 =  TAN(x2)/(poledist*COS(x2)*COS(x2)*r*r)

          D(1,1) =  D11
          D(1,2) =  D12
          D(2,1) =  D21
          D(2,2) =  D22

       end if
    end if

  end subroutine vmap

  ! ========================================
  ! coreolis_init:
  !
  ! Initialize coreolis term ...
  !
  ! ========================================

  subroutine coreolis_init(elem)
    use element_mod, only : element_t
    use dimensions_mod, only : nv
    use physical_constants, only : omega
    type (element_t) :: elem(:)

    ! Local variables

    integer                  :: i,j
    integer                  :: ii

    do ii=1,SIZE(elem)
       call coreolis_init_atomic(elem(ii))
    end do

  end subroutine coreolis_init

  ! ========================================
  ! coreolis_init_atomic:
  !
  ! Initialize coreolis term ...
  !
  ! ========================================

  subroutine coreolis_init_atomic(elem)
    use element_mod, only : element_t
    use dimensions_mod, only : nv
    use physical_constants, only : omega

    type (element_t) :: elem

    ! Local variables

    integer                  :: i,j
    real (kind=real_kind) :: lat,lon,rangle

    rangle = rotate_grid*DD_PI/180
    do j=1,nv
       do i=1,nv
             if ( rotate_grid /= 0) then
                lat = elem%spherev(i,j)%lat
                lon = elem%spherev(i,j)%lon
             	elem%fcor(i,j)= 2*omega* &
                     (-cos(lon)*cos(lat)*sin(rangle) + sin(lat)*cos(rangle))
             else
                elem%fcor(i,j) = 2.0D0*omega*SIN(elem%spherev(i,j)%lat)
             endif
       end do
    end do

  end subroutine coreolis_init_atomic

  ! =========================================
  ! rotation_init_atomic:
  !
  ! Initialize cube rotation terms resulting
  ! from changing cube face coordinate systems
  !
  ! =========================================

  subroutine rotation_init_atomic(elem, rot_type)
    use element_mod, only : element_t
    use dimensions_mod, only : nv
    use control_mod, only : north, south, east, west, neast, seast, swest, nwest

    type (element_t) :: elem
    character(len=*) rot_type

    ! =======================================
    ! Local variables
    ! =======================================

    integer :: ie,je
    integer :: myface_no        ! current element face number
    integer :: nbrface_no       ! neighbor element face number
    integer :: inbr
    integer :: nrot,irot
    integer :: ii,i,j
    integer :: ir,jr
    integer :: nbrnum

    real (kind=real_kind) :: Dloc(2,2,nv)
    real (kind=real_kind) :: Drem(2,2,nv)
    real (kind=real_kind) :: x1,x2


    call convert_gbl_index(elem%vertex%number,ie,je,myface_no)

    nrot   = 0

    elem%FaceNum = myface_no

    do inbr=1,8
       if (elem%vertex%nbrs(inbr) /= -1 ) then
          call convert_gbl_index(elem%vertex%nbrs(inbr),ie,je,nbrface_no)
          if (myface_no /= nbrface_no) nrot=nrot+1
       end if
    end do

    if(associated(elem%desc%rot)) then
       if (size(elem%desc%rot) > 0) then
          !         deallocate(elem%desc%rot)
          NULLIFY(elem%desc%rot)
       endif
    end if

    ! =====================================================
    ! If there are neighbors on other cube faces, allocate 
    ! an array of rotation matrix structs.
    ! =====================================================

    if (nrot > 0) then
       allocate(elem%desc%rot(nrot))
       elem%desc%use_rotation=1
       irot=0          
       do inbr=1,8

          call convert_gbl_index(elem%vertex%nbrs(inbr),ie,je,nbrface_no)

          ! The cube edge (myface_no,nbrface_no) and inbr defines 
          ! a unique rotation given by (D^-1) on myface_no x (D) on nbrface_no

          if (myface_no /= nbrface_no .and. elem%vertex%nbrs(inbr) /= -1 ) then           

             irot=irot+1

             if (inbr <= 4) then      
                allocate(elem%desc%rot(irot)%R(2,2,nv))  ! edge
             else                     
                allocate(elem%desc%rot(irot)%R(2,2,1 ))   ! corner
             end if

             ! must compute Dloc on my face, Drem on neighbor face, 
             ! for each point on edge or corner.

             ! ==================================== 
             ! Equatorial belt east/west neighbors
             ! ==================================== 

             if (nbrface_no <= 4 .and. myface_no <= 4) then

                if (inbr == west) then
                   do j=1,nv
                      x1 = elem%cartv(1,j)%x
                      x2 = elem%cartv(1,j)%y
                      call Vmap(Dloc(1,1,j), x1,x2,myface_no)
                      call Vmap(Drem(1,1,j),-x1,x2,nbrface_no)
                   end do
                else if (inbr == east) then
                   do j=1,nv
                      x1 = elem%cartv(nv,j)%x
                      x2 = elem%cartv(nv,j)%y
                      call Vmap(Dloc(1,1,j), x1,x2,myface_no)
                      call Vmap(Drem(1,1,j),-x1,x2,nbrface_no)
                   end do
                else if (inbr == swest ) then
                   x1 = elem%cartv(1,1)%x
                   x2 = elem%cartv(1,1)%y
                   call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                   call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                else if (inbr == nwest ) then
                   x1 = elem%cartv(1,nv)%x
                   x2 = elem%cartv(1,nv)%y
                   call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                   call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                else if (inbr == seast ) then
                   x1 = elem%cartv(nv,1)%x
                   x2 = elem%cartv(nv,1)%y
                   call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                   call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                else if (inbr == neast ) then
                   x1 = elem%cartv(nv,nv)%x
                   x2 = elem%cartv(nv,nv)%y
                   call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                   call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                end if

             end if

             ! Northern Neighbors of Equatorial Belt

             if ( myface_no <= 4 .and. nbrface_no == 6 ) then
                if (inbr == north) then
                   do i=1,nv
                      ir=nv+1-i
                      x1 = elem%cartv(i,nv)%x
                      x2 = elem%cartv(i,nv)%y
                      if ( myface_no == 1) then
                         call Vmap(Dloc(1,1,i), x1,x2,myface_no)
                         call Vmap(Drem(1,1,i),x1,-x2,nbrface_no)
                      end if
                      if ( myface_no == 2) then
                         call Vmap(Dloc(1,1,i),x1,x2,myface_no)
                         call Vmap(Drem(1,1,i),x2,x1,nbrface_no)
                         nbrnum = elem%vertex%nbrs(inbr)

                      end if
                      if ( myface_no == 3) then
                         call Vmap(Dloc(1,1,ir), x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                      end if
                      if ( myface_no == 4) then
                         call Vmap(Dloc(1,1,ir), x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x2,-x1,nbrface_no)
                      end if
                   end do
                else if (inbr == nwest) then
                   x1 = elem%cartv(1,nv)%x
                   x2 = elem%cartv(1,nv)%y
                   call Vmap(Dloc(1,1,1), x1,x2,myface_no)
                   if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   if ( myface_no == 2) call Vmap(Drem(1,1,1),x2, x1,nbrface_no)
                   if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   if ( myface_no == 4) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                else if (inbr == neast) then
                   x1 = elem%cartv(nv,nv)%x
                   x2 = elem%cartv(nv,nv)%y
                   call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                   if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   if ( myface_no == 2) call Vmap(Drem(1,1,1),x2, x1,nbrface_no)
                   if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   if ( myface_no == 4) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                end if

             end if

             ! Southern Neighbors of Equatorial Belt

             if ( myface_no <= 4 .and. nbrface_no == 5 ) then
                if (inbr == south) then
                   do i=1,nv
                      ir=nv+1-i
                      x1 = elem%cartv(i,1)%x
                      x2 = elem%cartv(i,1)%y
                      if ( myface_no == 1) then
                         call Vmap(Dloc(1,1,i), x1, x2,myface_no)
                         call Vmap(Drem(1,1,i), x1,-x2,nbrface_no)
                      end if
                      if ( myface_no == 2) then
                         call Vmap(Dloc(1,1,ir),x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x2,-x1,nbrface_no)
                      end if
                      if ( myface_no == 3) then
                         call Vmap(Dloc(1,1,ir), x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                      end if
                      if ( myface_no == 4) then
                         call Vmap(Dloc(1,1,i), x1,x2,myface_no)
                         call Vmap(Drem(1,1,i), x2,x1,nbrface_no)
                      end if
                   end do
                else if (inbr == swest) then
                   x1 = elem%cartv(1,1)%x
                   x2 = elem%cartv(1,1)%y
                   call Vmap(Dloc(1,1,1),x1,x2,myface_no)


                   if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   if ( myface_no == 2) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   if ( myface_no == 4) call Vmap(Drem(1,1,1),x2,x1,nbrface_no)

                else if (inbr == seast) then
                   x1 = elem%cartv(nv,1)%x
                   x2 = elem%cartv(nv,1)%y
                   call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                   if ( myface_no == 1) call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   if ( myface_no == 2) call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   if ( myface_no == 3) call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   if ( myface_no == 4) call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                end if

             end if

             ! Neighbors of Northern Capping Face Number 6

             if ( myface_no == 6 ) then
                if (nbrface_no == 1) then
                   if (inbr == south) then
                      do i=1,nv
                         x1 = elem%cartv(i,1)%x
                         x2 = elem%cartv(i,1)%y
                         call Vmap(Dloc(1,1,i),x1,x2,myface_no)
                         call Vmap(Drem(1,1,i),x1,-x2,nbrface_no)
                      end do
                   else if (inbr == swest) then
                      x1 = elem%cartv(1,1)%x
                      x2 = elem%cartv(1,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   else if (inbr == seast) then
                      x1 = elem%cartv(nv,1)%x
                      x2 = elem%cartv(nv,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   end if
                else if (nbrface_no == 2) then
                   if (inbr == east) then
                      do j=1,nv
                         x1 = elem%cartv(nv,j)%x
                         x2 = elem%cartv(nv,j)%y
                         call Vmap(Dloc(1,1,j),x1,x2,myface_no)
                         call Vmap(Drem(1,1,j),x2,x1,nbrface_no)
                      end do
                   else if (inbr == seast) then
                      x1 = elem%cartv(nv,1)%x
                      x2 = elem%cartv(nv,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                   else if (inbr == neast) then
                      x1 = elem%cartv(nv,nv)%x
                      x2 = elem%cartv(nv,nv)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                   end if
                else if (nbrface_no == 3) then
                   if (inbr == north) then
                      do i=1,nv
                         ir =nv+1-i
                         x1 = elem%cartv(i,nv)%x
                         x2 = elem%cartv(i,nv)%y
                         call Vmap(Dloc(1,1,ir),x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                      end do
                   else if (inbr == nwest) then
                      x1 = elem%cartv(1,nv)%x
                      x2 = elem%cartv(1,nv)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   else if (inbr == neast) then
                      x1 = elem%cartv(nv,nv)%x
                      x2 = elem%cartv(nv,nv)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   end if
                else if (nbrface_no == 4) then
                   if (inbr == west) then
                      do j=1,nv
                         jr=nv+1-j
                         x1 = elem%cartv(1,j)%x
                         x2 = elem%cartv(1,j)%y
                         call Vmap(Dloc(1,1,jr), x1, x2,myface_no )
                         call Vmap(Drem(1,1,jr),-x2,-x1,nbrface_no)
                      end do
                   else if (inbr == swest) then
                      x1 = elem%cartv(1,1)%x
                      x2 = elem%cartv(1,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   else if (inbr == nwest) then
                      x1 = elem%cartv(1,nv)%x
                      x2 = elem%cartv(1,nv)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   end if
                end if
             end if

             ! Neighbors of South Capping Face Number 5

             if ( myface_no == 5 ) then
                if (nbrface_no == 1) then
                   if (inbr == north) then
                      do i=1,nv
                         x1 = elem%cartv(i,nv)%x
                         x2 = elem%cartv(i,nv)%y
                         call Vmap(Dloc(1,1,i),x1,x2,myface_no)
                         call Vmap(Drem(1,1,i),x1,-x2,nbrface_no)
                      end do
                   else if (inbr == nwest) then
                      x1 = elem%cartv(1,nv)%x
                      x2 = elem%cartv(1,nv)%y
                      call Vmap(Dloc(:,:,1),x1,x2,myface_no)
                      call Vmap(Drem(:,:,1),x1,-x2,nbrface_no)
                   else if (inbr == neast) then
                      x1 = elem%cartv(nv,nv)%x
                      x2 = elem%cartv(nv,nv)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x1,-x2,nbrface_no)
                   end if
                else if (nbrface_no == 2) then
                   if (inbr == east) then
                      do j=1,nv
                         jr=nv+1-j
                         x1 = elem%cartv(nv,j)%x
                         x2 = elem%cartv(nv,j)%y
                         call Vmap(Dloc(1,1,jr),x1,  x2,myface_no)
                         call Vmap(Drem(1,1,jr),-x2,-x1,nbrface_no)
                      end do
                   else if (inbr == seast) then
                      x1 = elem%cartv(nv,1)%x
                      x2 = elem%cartv(nv,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   else if (inbr == neast) then
                      x1 = elem%cartv(nv,nv)%x
                      x2 = elem%cartv(nv,nv)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x2,-x1,nbrface_no)
                   end if
                else if (nbrface_no == 3) then
                   if (inbr == south) then
                      do i=1,nv
                         ir=nv+1-i
                         x1 = elem%cartv(i,1)%x
                         x2 = elem%cartv(i,1)%y
                         call Vmap(Dloc(1,1,ir),x1,x2,myface_no)
                         call Vmap(Drem(1,1,ir),-x1,x2,nbrface_no)
                      end do
                   else if (inbr == swest) then
                      x1 = elem%cartv(1,1)%x
                      x2 = elem%cartv(1,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   else if (inbr == seast) then
                      x1 = elem%cartv(nv,1)%x
                      x2 = elem%cartv(nv,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),-x1,x2,nbrface_no)
                   end if
                else if (nbrface_no == 4) then
                   if (inbr == west) then
                      do j=1,nv
                         x1 = elem%cartv(1,j)%x
                         x2 = elem%cartv(1,j)%y
                         call Vmap(Dloc(1,1,j),x1,x2,myface_no)
                         call Vmap(Drem(1,1,j),x2,x1,nbrface_no)
                      end do
                   else if (inbr == swest) then
                      x1 = elem%cartv(1,1)%x
                      x2 = elem%cartv(1,1)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                   else if (inbr == nwest) then
                      x1 = elem%cartv(1,nv)%x
                      x2 = elem%cartv(1,nv)%y
                      call Vmap(Dloc(1,1,1),x1,x2,myface_no)
                      call Vmap(Drem(1,1,1),x2,x1,nbrface_no)
                   end if
                end if
             end if

             elem%desc%rot(irot)%nbr = inbr
             if (rot_type == "covariant") then
                do i=1,SIZE(elem%desc%rot(irot)%R(:,:,:),3)
                   elem%desc%rot(irot)%R(:,:,i)=covariant_rot(Dloc(:,:,i),Drem(:,:,i))
                end do
             else if (rot_type == "contravariant") then
                do i=1,SIZE(elem%desc%rot(irot)%R(:,:,:),3)
                   elem%desc%rot(irot)%R(:,:,i)=contravariant_rot(Dloc(:,:,i),Drem(:,:,i))
                end do
             end if

          endif

       end do
    end if

  end subroutine rotation_init_atomic

  ! ================================================
  ! convert_gbl_index:
  !
  ! Convert global element index to cube index
  ! ================================================

  subroutine convert_gbl_index(number,ie,je,face_no)
    use dimensions_mod, only : ne
    integer, intent(in)  :: number
    integer, intent(out) :: ie,je,face_no

    face_no=((number-1)/(ne*ne))+1
    ie=MODULO(number-1,ne)
    je=(number-1)/ne - (face_no-1)*ne

  end subroutine convert_gbl_index

#if 0
  ! ================================
  ! cube_topology:
  !
  ! ================================

  subroutine cube_topology(graph,ne,nv) 
    type (vertex_t),intent(out) :: graph(6*ne*ne)
    integer        ,intent(in)  :: ne
    integer        ,intent(in)  :: nv

    integer :: ii,jj,ir,k,l
    integer :: icube(ne,ne,6)

    integer edge,corner,rev

    edge  = nv
    corner= 1
    rev   = 0

    do k=1,6 
       do jj=1,ne
          do ii=1,ne
             icube(ii,jj,k)=ii+(jj-1)*ne+(k-1)*ne*ne
             graph(icube(ii,jj,k))%number=icube(ii,jj,k)
          end do
       end do
    end do

    ! =============================================
    ! Initialize the diagonal neighbors at cube
    ! face corners to be null...
    ! =============================================

    do k=1,6
       call nbrinit_sub(graph(icube(1 ,1 ,k))%nbr(swest),-1,neast,corner,rev)
       call nbrinit_sub(graph(icube(1 ,ne,k))%nbr(nwest),-1,seast,corner,rev)
       call nbrinit_sub(graph(icube(ne,1 ,k))%nbr(seast),-1,nwest,corner,rev)
       call nbrinit_sub(graph(icube(ne,ne,k))%nbr(neast),-1,swest,corner,rev)
    end do

    ! ======================
    ! cube face interiors
    ! ======================

    do k=1,6

       ! ===================================
       ! west neighbor
       ! ===================================

       do jj=1,ne
          do ii=2,ne
             call nbrinit_sub(graph(icube(ii,jj,k))%nbr(west),icube(ii-1,jj,k),east,edge,rev)  
          end do
       end do

       ! ===================================
       ! east neighbor
       ! ===================================

       do jj=1,ne
          do ii=1,ne-1
             call nbrinit_sub(graph(icube(ii,jj,k))%nbr(east),icube(ii+1,jj,k),west,edge,rev)
          end do
       end do

       ! ===================================
       ! south neighbor
       ! ===================================

       do jj=2,ne
          do ii=1,ne
             call nbrinit_sub(graph(icube(ii,jj,k))%nbr(south),icube(ii,jj-1,k),north,edge,rev)
          end do
       end do

       ! ===================================
       ! north neighbor
       ! ===================================

       do jj=1,ne-1
          do ii=1,ne
             call nbrinit_sub(graph(icube(ii,jj,k))%nbr(north),icube(ii,jj+1,k),south,edge,rev)
          end do
       end do

       ! ===================================
       ! south west neighbor
       ! ===================================

       do jj=2,ne
          do ii=2,ne
             call nbrinit_sub(graph(icube(ii,jj,k))%nbr(swest),icube(ii-1,jj-1,k),neast,corner,rev)
          end do
       end do

       ! ===================================
       ! south east neighbor
       ! ===================================

       do jj=2,ne
          do ii=1,ne-1
             call nbrinit_sub(graph(icube(ii,jj,k))%nbr(seast),icube(ii+1,jj-1,k),nwest,corner,rev)
          end do
       end do

       ! ===================================
       ! north west neighbor
       ! ===================================

       do jj=1,ne-1
          do ii=2,ne
             call nbrinit_sub(graph(icube(ii,jj,k))%nbr(nwest),icube(ii-1,jj+1,k),seast,corner,rev)
          end do
       end do

       ! ===================================
       ! north east neighbor
       ! ===================================

       do jj=1,ne-1
          do ii=1,ne-1
             call nbrinit_sub(graph(icube(ii,jj,k))%nbr(neast),icube(ii+1,jj+1,k),swest,corner,rev)
          end do
       end do

    end do

    ! ===================================
    ! cube faces 1-4 west/east edges
    ! ===================================

    do k=1,4

       ! ====================
       ! west neighbor
       ! ====================

       do jj=1,ne
          call nbrinit_sub(graph(icube(1 ,jj,k))%nbr(west),icube(ne,jj,MODULO(2+k,4)+1),east,edge,rev)
       end do

       ! ====================
       ! south west neighbor
       ! ====================

       do jj=2,ne
          call nbrinit_sub(graph(icube(1 ,jj,k))%nbr(swest),icube(ne,jj-1,MODULO(2+k,4)+1),neast,corner,rev)
       end do


       ! ====================
       ! north west neighbor
       ! ====================

       do jj=1,ne-1
          call nbrinit_sub(graph(icube(1 ,jj,k))%nbr(nwest),icube(ne,jj+1,MODULO(2+k,4)+1),seast,corner,rev)
       end do

       ! ====================
       ! east neighbor
       ! ====================

       do jj=1,ne
          call nbrinit_sub(graph(icube(ne,jj,k))%nbr(east),icube(1 ,jj,MODULO(k  ,4)+1),west,edge,rev)
       end do

       ! ====================
       ! south east neighbor
       ! ====================

       do jj=2,ne
          call nbrinit_sub(graph(icube(ne,jj,k))%nbr(seast),icube(1 ,jj-1,MODULO(k  ,4)+1),nwest,corner,rev)
       end do

       ! ====================
       ! north east neighbor
       ! ====================

       do jj=1,ne-1
          call nbrinit_sub(graph(icube(ne,jj,k))%nbr(neast),icube(1 ,jj+1,MODULO(k  ,4)+1),swest,corner,rev)
       end do

    end do

    ! =============================
    ! 1 north edge - 6 south edge
    ! =============================

    do ii=1,ne
       call nbrinit_sub(graph(icube(ii,ne,1))%nbr(north),icube(ii,1 ,6),south,edge,rev)
       call nbrinit_sub(graph(icube(ii,1 ,6))%nbr(south),icube(ii,ne,1),north,edge,rev)
    end do

    do ii=2,ne
       call nbrinit_sub(graph(icube(ii  ,ne,1))%nbr(nwest),icube(ii-1,1 ,6),seast,corner,rev)
       call nbrinit_sub(graph(icube(ii-1,1 ,6))%nbr(seast),icube(ii  ,ne,1),nwest,corner,rev)
    end do

    do ii=1,ne-1

       call nbrinit_sub(graph(icube(ii  ,ne,1))%nbr(neast),icube(ii+1,1 ,6),swest,corner,rev)
       call nbrinit_sub(graph(icube(ii+1,1 ,6))%nbr(swest),icube(ii  ,ne,1),neast,corner,rev)
    end do

    ! =============================
    ! 2 north edge - 6 east edge
    ! =============================

    do ii=1,ne
       call nbrinit_sub(graph(icube(ii,ne,2))%nbr(north),icube(ne,ii,6),east,edge,rev)
       call nbrinit_sub(graph(icube(ne,ii,6))%nbr(east), icube(ii,ne,2),north,edge,rev)
    end do

    do ii=2,ne
       call nbrinit_sub(graph(icube(ii,ne  ,2))%nbr(nwest),icube(ne,ii-1,6),neast,corner,rev)
       call nbrinit_sub(graph(icube(ne,ii-1,6))%nbr(neast),icube(ii,ne  ,2),nwest,corner,rev)
    end do

    do ii=1,ne-1
       call nbrinit_sub(graph(icube(ii,ne  ,2))%nbr(neast),icube(ne,ii+1,6),seast,corner,rev)
       call nbrinit_sub(graph(icube(ne,ii+1,6))%nbr(seast),icube(ii,ne  ,2),neast,corner,rev)
    end do

    ! =============================
    ! 3 north edge - 6 north edge
    ! =============================

    do ii=1,ne
       ir=ne+1-ii
       call nbrinit_sub(graph(icube(ii,ne,3))%nbr(north),icube(ir,ne,6),north,edge,1-rev)
       call nbrinit_sub(graph(icube(ii,ne,6))%nbr(north),icube(ir,ne,3),north,edge,1-rev)
    end do

    do ii=2,ne
       ir=ne+1-(ii-1)
       call nbrinit_sub(graph(icube(ii,ne,3))%nbr(nwest),icube(ir,ne,6),nwest,corner,rev)
       call nbrinit_sub(graph(icube(ii,ne,6))%nbr(nwest),icube(ir,ne,3),nwest,corner,rev)
    end do

    do ii=1,ne-1
       ir=ne+1-(ii+1)
       call nbrinit_sub(graph(icube(ii,ne,3))%nbr(neast),icube(ir,ne,6),neast,corner,rev)
       call nbrinit_sub(graph(icube(ii,ne,6))%nbr(neast),icube(ir,ne,3),neast,corner,rev)
    end do

    ! =============================
    ! 4 north edge - 6 west edge
    ! =============================

    do ii=1,ne
       ir=ne+1-ii
       call nbrinit_sub(graph(icube(ii,ne,4))%nbr(north),icube(1 ,ir,6),west ,edge,1-rev)
       call nbrinit_sub(graph(icube(1 ,ii,6))%nbr(west) ,icube(ir,ne,4),north,edge,1-rev)
    end do

    do ii=2,ne
       ir=ne+1-(ii-1)
       call nbrinit_sub(graph(icube(ii,ne,4))%nbr(nwest),icube(1 ,ir,6),swest,corner,rev)
       call nbrinit_sub(graph(icube(1 ,ii,6))%nbr(swest),icube(ir,ne,4),nwest,corner,rev)
    end do

    do ii=1,ne-1
       ir=ne+1-(ii+1)
       call nbrinit_sub(graph(icube(ii,ne,4))%nbr(neast),icube(1 ,ir,6),nwest,corner,rev)
       call nbrinit_sub(graph(icube(1 ,ii,6))%nbr(nwest),icube(ir,ne,4),neast,corner,rev)
    end do

    ! =============================
    ! 1 south edge - 5 north edge 
    ! =============================

    do ii=1,ne
       call nbrinit_sub(graph(icube(ii,1 ,1))%nbr(south),icube(ii,ne,5),north,edge,rev)
       call nbrinit_sub(graph(icube(ii,ne,5))%nbr(north),icube(ii,1 ,1),south,edge,rev)
    end do

    do ii=2,ne
       call nbrinit_sub(graph(icube(ii  ,1 ,1))%nbr(swest),icube(ii-1,ne,5),neast,corner,rev)
       call nbrinit_sub(graph(icube(ii-1,ne,5))%nbr(neast),icube(ii  ,1 ,1),swest,corner,rev)
    end do

    do ii=1,ne-1
       call nbrinit_sub(graph(icube(ii  ,1 ,1))%nbr(seast),icube(ii+1,ne,5),nwest,corner,rev)
       call nbrinit_sub(graph(icube(ii+1,ne,5))%nbr(nwest),icube(ii  ,1 ,1),seast,corner,rev)
    end do

    ! =============================
    ! 2 south edge - 5 east edge 
    ! =============================

    do ii=1,ne
       ir = ne+1-ii
       call nbrinit_sub(graph(icube(ii,1 ,2))%nbr(south),icube(ne,ir,5),east,edge,1-rev)
       call nbrinit_sub(graph(icube(ne,ii,5))%nbr(east) ,icube(ir,1 ,2),south,edge,1-rev)
    end do

    do ii=2,ne
       ir = ne+1-(ii-1)
       call nbrinit_sub(graph(icube(ii,1 ,2))%nbr(swest),icube(ne,ir,5),seast,corner,rev)
       call nbrinit_sub(graph(icube(ne,ii,5))%nbr(seast),icube(ir,1 ,2),swest,corner,rev)
    end do

    do ii=1,ne-1
       ir = ne+1-(ii+1)
       call nbrinit_sub(graph(icube(ii,1 ,2))%nbr(seast),icube(ne,ir,5),neast,corner,rev)
       call nbrinit_sub(graph(icube(ne,ii,5))%nbr(neast),icube(ir,1 ,2),seast,corner,rev)
    end do

    ! =============================
    ! 3 south edge - 5 south edge 
    ! =============================

    do ii=1,ne
       ir=ne+1-ii
       call nbrinit_sub(graph(icube(ii,1 ,3))%nbr(south),icube(ir,1 ,5),south,edge,1-rev)
       call nbrinit_sub(graph(icube(ii,1 ,5))%nbr(south),icube(ir,1 ,3),south,edge,1-rev)
    end do

    do ii=2,ne
       ir=ne+1-(ii-1)
       call nbrinit_sub(graph(icube(ii,1 ,3))%nbr(swest),icube(ir,1 ,5),swest,corner,rev)
       call nbrinit_sub(graph(icube(ii,1 ,5))%nbr(swest),icube(ir,1 ,3),swest,corner,rev)
    end do

    do ii=1,ne-1
       ir=ne+1-(ii+1)
       call nbrinit_sub(graph(icube(ii,1 ,3))%nbr(seast),icube(ir,1 ,5),seast,corner,rev)
       call nbrinit_sub(graph(icube(ii,1 ,5))%nbr(seast),icube(ir,1 ,3),seast,corner,rev)
    end do

    ! =============================
    ! 4 south edge - 5 west edge 
    ! =============================

    do ii=1,ne
       call nbrinit_sub(graph(icube(ii,1 ,4))%nbr(south),icube(1 ,ii,5),west,edge,rev)
       call nbrinit_sub(graph(icube(1 ,ii,5))%nbr(west) ,icube(ii,1 ,4),south,edge,rev)  
    end do

    do ii=2,ne
       call nbrinit_sub(graph(icube(ii,1 ,4))%nbr(swest),icube(1 ,ii-1,5),nwest,corner,rev)
       call nbrinit_sub(graph(icube(1 ,ii-1,5))%nbr(nwest),icube(ii,1   ,4),swest,corner,rev)
    end do

    do ii=1,ne-1
       call nbrinit_sub(graph(icube(ii,1 ,4))%nbr(seast),icube(1 ,ii+1,5),swest,corner,rev)
       call nbrinit_sub(graph(icube(1 ,ii+1,5))%nbr(swest),icube(ii,1   ,4),seast,corner,rev)
    end do


  end subroutine cube_topology
#endif

  subroutine CubeTopology(GridEdge, GridVertex)
    use params_mod, only : RECURSIVE, SFCURVE
    use control_mod, only: partmethod
    use gridgraph_mod, only : GridEdge_t, GridVertex_t, initgridedge
    use dimensions_mod, only : nv, np, ne
    use spacecurve_mod, only :  IsFactorable, genspacecurve
    use control_mod, only : north, south, east, west, neast, seast, swest, nwest
    use parallel_mod, only : abortmp
    !-----------------------
    implicit none

    type (GridEdge_t),   intent(out),target     :: GridEdge(:)
    type (GridVertex_t), intent(out),target     :: GridVertex(:)


    integer,allocatable       :: Mesh(:,:)
    integer,allocatable       :: Mesh2(:,:),Mesh2_map(:,:,:),sfcij(:,:)
    type (GridVertex_t),allocatable        :: GridElem(:,:,:)
    integer                   :: i,j,k,number,irev,ne2,i2,j2,sfc_index
    integer                   :: EdgeWgtV,EdgeWgtP,CornerWgt
    integer                   :: ielem,nedge

    integer                   :: offset, ierr

#if 0
    ! ===================================
    ! character strings
    ! ===================================

    character(len=6)              :: charne,charnpart
    character(len=80)             :: cornerfile
    character(len=80)             :: labelfile
    character(len=80)             :: partfile

    !=================================================
    ! Setup the Output Files if necessary
    !=================================================
    if(OutputFiles) then
       ! =============================================
       ! Create 3 files... each file called
       !  xxx{ne}.cube
       !
       ! labelfile: xxx=label
       !    label each corner of each element with a
       !    unique integer (1,... 4*number_of_elements)
       !
       ! cornerfile: xxx=corner
       !    3-D cartesian coordinates of element corners
       !
       ! vertexfile: xxx=vertex
       !    METIS element connectivity graph
       !
       ! =============================================

       write(charne,'(i6)')ne

       labelfile  = "cube_"//TRIM(ADJUSTL(charne))//"x"//TRIM(ADJUSTL(charne))//".con"
       cornerfile = "cube_"//TRIM(ADJUSTL(charne))//"x"//TRIM(ADJUSTL(charne))//".xyz"

       write(charnpart,'(i6)')npart
       partfile   = "cube_"//TRIM(ADJUSTL(charne))//"x"//TRIM(ADJUSTL(charne))//".part."//TRIM(ADJUSTL(charnpart))

       open(7 ,file=labelfile ,form="formatted",status="unknown")
       open(8 ,file=cornerfile,form="formatted",status="unknown")
       open(10 ,file=partfile,form="formatted",status="unknown")
    endif
#endif
    allocate(GridElem(ne,ne,nfaces),stat=ierr)
    if(ierr/=0) then
       call abortmp('error in allocation of GridElem structure')
    end if

    number=1
    EdgeWgtV   = nv
    EdgeWgtP   = np
    CornerWgt = 1
    do k=1,nfaces
       do j=1,ne
          do i=1,ne

             ! ====================================
             ! Number elements
             ! ====================================


             !   Adding the corner 'edges'
             GridElem(i,j,k)%degree=8

#if 0
             allocate(GridElem(i,j,k)%nbrs(GridElem(i,j,k)%degree))
             allocate(GridElem(i,j,k)%wgt(GridElem(i,j,k)%degree))
#endif

             ! Do some initalization here
             GridElem(i,j,k)%nbrs(:) = -1
             GridElem(i,j,k)%wgtV(:)=0
             GridElem(i,j,k)%wgtP(:)=0
             GridElem(i,j,k)%SpaceCurve=0
             GridElem(i,j,k)%number=number 
             number=number+1

#if 0
             call set_GridVertex_number(GridElem(i,j,k),number)
             if(OutputFiles) then
                write(7,*)(number-1)*4+1,(number-1)*4+2,(number-1)*4+3,(number-1)*4+4
             endif

             ! ====================================
             ! Initialize location on cube/sphere
             ! ====================================

             start%x = -1.0D0 + (2.0D0/ne)*(i-1)
             start%y = -1.0D0 + (2.0D0/ne)*(j-1)

             end%x = -1.0D0 + (2.0D0/ne)*(i)
             end%y = -1.0D0 + (2.0D0/ne)*(j)


             call set_GridVertex_location(GridElem(i,j,k),start,end,face(k))

             ! =======================================
             ! Output data for plotting quadrilateral
             ! elements on spheres
             ! =======================================
             corners(1)=start
             corners(1)%x=start%x  ! GridElem(i,j,k)%sw%x
             corners(1)%y=start%y  ! GridElem(i,j,k)%sw%y
             corners(1)%face=>face(k)

             corners(2)%x=end%x    ! GridElem(i,j,k)%se%x
             corners(2)%y=start%y  !GridElem(i,j,k)%se%y
             corners(2)%face=>face(k)

             corners(3)%x=end%x  ! GridElem(i,j,k)%ne%x
             corners(3)%y=end%y  ! GridElem(i,j,k)%ne%y
             corners(3)%face=>face(k)

             corners(4)%x=start%x  ! GridElem(i,j,k)%nw%x
             corners(4)%y=end%y    ! GridElem(i,j,k)%nw%y
             corners(4)%face=>face(k)


             call project_old(corners,cart_corners)

             if(OutputFiles) then
                write(8,*)cart_corners(1)%x,cart_corners(1)%y,cart_corners(1)%z
                write(8,*)cart_corners(2)%x,cart_corners(2)%y,cart_corners(2)%z
                write(8,*)cart_corners(3)%x,cart_corners(3)%y,cart_corners(3)%z
                write(8,*)cart_corners(4)%x,cart_corners(4)%y,cart_corners(4)%z
             endif

#endif
          end do
       end do
    end do

    !    print *,'CubeTopology: Ne, IsFactorable, IsLoadBalanced : ',ne,IsFactorable(ne),IsLoadBalanced(nelem,npart)

    allocate(Mesh(ne,ne))
    if(IsFactorable(ne)) then
       call GenspaceCurve(Mesh)
       !      call PrintCurve(Mesh) 
    else
       ! find the smallest ne2 which is a power of 2 and ne2>ne
       ne2=2**ceiling( log(real(ne))/log(2d0) )
       if (ne2<ne) call abortmp('Fatel SFC error')

       allocate(Mesh2(ne2,ne2))
       allocate(Mesh2_map(ne2,ne2,2))
       allocate(sfcij(0:ne2*ne2,2))

       call GenspaceCurve(Mesh2)  ! SFC partition for ne2

       ! associate every element on the ne x ne mesh (Mesh)
       ! with its closest element on the ne2 x ne2 mesh (Mesh2)
       ! Store this as a map from Mesh2 -> Mesh in Mesh2_map.
       ! elements in Mesh2 which are not mapped get assigned a value of 0
       Mesh2_map=0
       do j=1,ne
          do i=1,ne
             ! map this element to an (i2,j2) element
             ! [ (i-.5)/ne , (j-.5)/ne ]  = [ (i2-.5)/ne2 , (j2-.5)/ne2 ]
             i2=nint( ((i-.5)/ne)*ne2 + .5 )
             j2=nint( ((j-.5)/ne)*ne2 + .5 )
             if (i2<1) i2=1
             if (i2>ne2) i2=ne2
             if (j2<1) j2=1
             if (j2>ne2) j2=ne2
             Mesh2_map(i2,j2,1)=i
             Mesh2_map(i2,j2,2)=j
          enddo
       enddo

       ! create a reverse index array for Mesh2
       ! k = Mesh2(i,j) 
       ! (i,j) = (sfcij(k,1),sfci(k,2)) 
       do j=1,ne2
          do i=1,ne2
             k=Mesh2(i,j)
             sfcij(k,1)=i
             sfcij(k,2)=j
          enddo
       enddo

       ! generate a SFC for Mesh with the same ordering as the 
       ! elements in Mesh2 which map to Mesh.
       sfc_index=0
       do k=0,ne2*ne2-1
          i2=sfcij(k,1)
          j2=sfcij(k,2)
          i=Mesh2_map(i2,j2,1)
          j=Mesh2_map(i2,j2,2)
          if (i/=0) then
             ! (i2,j2) element maps to (i,j) element
             Mesh(i,j)=sfc_index
             sfc_index=sfc_index+1
          endif
       enddo
#if 0
       print *,'SFC Mapping to non powers of 2,3 used.  Mesh:'  
       do j=1,ne
          write(*,'(99i3)') (Mesh(i,j),i=1,ne)
       enddo
       call PrintCurve(Mesh2) 
#endif
       deallocate(Mesh2)
       deallocate(Mesh2_map)
       deallocate(sfcij)
    endif


    ! -------------------------------------------
    !  Setup the space-filling curve for face 1
    ! -------------------------------------------
    offset=0
    do j=1,ne
       do i=1,ne
          GridElem(i,j,1)%SpaceCurve = offset + Mesh(i,ne-j+1)
       enddo
    enddo

    ! -------------------------------------------
    !  Setup the space-filling curve for face 2
    ! -------------------------------------------
    offset = offset + ne*ne
    do j=1,ne
       do i=1,ne
          GridElem(i,j,2)%SpaceCurve = offset + Mesh(i,ne-j+1)
       enddo
    enddo

    ! -------------------------------------------
    !  Setup the space-filling curve for face 6
    ! -------------------------------------------
    offset = offset + ne*ne
    do j=1,ne
       do i=1,ne
          GridElem(i,j,6)%SpaceCurve = offset + Mesh(ne-i+1,ne-j+1)
       enddo
    enddo

    ! -------------------------------------------
    !  Setup the space-filling curve for face 4
    ! -------------------------------------------
    offset = offset + ne*ne
    do j=1,ne
       do i=1,ne
          GridElem(i,j,4)%SpaceCurve = offset + Mesh(ne-j+1,i)
       enddo
    enddo

    ! -------------------------------------------
    !  Setup the space-filling curve for face 5
    ! -------------------------------------------
    offset = offset + ne*ne
    do j=1,ne
       do i=1,ne
          GridElem(i,j,5)%SpaceCurve = offset + Mesh(i,j)
       enddo
    enddo


    ! -------------------------------------------
    !  Setup the space-filling curve for face 3
    ! -------------------------------------------
    offset = offset + ne*ne
    do j=1,ne
       do i=1,ne
          GridElem(i,j,3)%SpaceCurve = offset + Mesh(i,j)
       enddo
    enddo

    ! ==================
    ! face interiors
    ! ==================
    do k=1,6
       ! setup  SOUTH, WEST, SW neighbors
       do j=2,ne
          do i=2,ne
             GridElem(i,j,k)%nbrs(west)  = GridElem(i-1,j,k)%number
             GridElem(i,j,k)%wgtV(west)   = EdgeWgtV
             GridElem(i,j,k)%wgtP(west)   = EdgeWgtP
             GridElem(i,j,k)%nbrs(south) = GridElem(i,j-1,k)%number
             GridElem(i,j,k)%wgtV(south)  = EdgeWgtV
             GridElem(i,j,k)%wgtP(south)  = EdgeWgtP
             GridElem(i,j,k)%nbrs(swest)    = GridElem(i-1,j-1,k)%number
             GridElem(i,j,k)%wgtV(swest)     = CornerWgt
             GridElem(i,j,k)%wgtP(swest)     = CornerWgt
          end do
       end do

       !  setup EAST, NORTH, NE neighbors
       do j=1,ne-1
          do i=1,ne-1
             GridElem(i,j,k)%nbrs(east)   = GridElem(i+1,j,k)%number
             GridElem(i,j,k)%wgtV(east)    = EdgeWgtV
             GridElem(i,j,k)%wgtP(east)    = EdgeWgtP
             GridElem(i,j,k)%nbrs(north)  = GridElem(i,j+1,k)%number
             GridElem(i,j,k)%wgtV(north)   = EdgeWgtV
             GridElem(i,j,k)%wgtP(north)   = EdgeWgtP
             GridElem(i,j,k)%nbrs(neast)     = GridElem(i+1,j+1,k)%number
             GridElem(i,j,k)%wgtV(neast)      = CornerWgt
             GridElem(i,j,k)%wgtP(neast)      = CornerWgt
          end do
       end do

       ! Setup the remaining SOUTH, EAST, and SE neighbors
       do j=2,ne
          do i=1,ne-1
             GridElem(i,j,k)%nbrs(south)  = GridElem(i,j-1,k)%number
             GridElem(i,j,k)%wgtV(south)   = EdgeWgtV
             GridElem(i,j,k)%wgtP(south)   = EdgeWgtP
             GridElem(i,j,k)%nbrs(east)   = GridElem(i+1,j,k)%number
             GridElem(i,j,k)%wgtV(east)    = EdgeWgtV
             GridElem(i,j,k)%wgtP(east)    = EdgeWgtP
             GridElem(i,j,k)%nbrs(seast)     = GridElem(i+1,j-1,k)%number
             GridElem(i,j,k)%wgtV(seast)      = CornerWgt
             GridElem(i,j,k)%wgtP(seast)      = CornerWgt
          enddo
       enddo

       ! Setup the remaining NORTH, WEST, and NW neighbors
       do j=1,ne-1
          do i=2,ne
             GridElem(i,j,k)%nbrs(north)  = GridElem(i,j+1,k)%number
             GridElem(i,j,k)%wgtV(north)   = EdgeWgtV
             GridElem(i,j,k)%wgtP(north)   = EdgeWgtP
             GridElem(i,j,k)%nbrs(west)   = GridElem(i-1,j,k)%number
             GridElem(i,j,k)%wgtV(west)    = EdgeWgtV
             GridElem(i,j,k)%wgtP(west)    = EdgeWgtP
             GridElem(i,j,k)%nbrs(nwest)     = GridElem(i-1,j+1,k)%number
             GridElem(i,j,k)%wgtV(nwest)      = CornerWgt
             GridElem(i,j,k)%wgtP(nwest)      = CornerWgt
          enddo
       enddo
    end do

    ! ======================
    ! west/east "belt" edges
    ! ======================

    do k=1,4
       do j=1,ne
          GridElem(1 ,j,k)%nbrs(west) = GridElem(ne,j,MODULO(2+k,4)+1)%number
          GridElem(1 ,j,k)%wgtV(west)  = EdgeWgtV
          GridElem(1 ,j,k)%wgtP(west)  = EdgeWgtP
          GridElem(ne,j,k)%nbrs(east) = GridElem(1 ,j,MODULO(k  ,4)+1)%number
          GridElem(ne,j,k)%wgtV(east)  = EdgeWgtV
          GridElem(ne,j,k)%wgtP(east)  = EdgeWgtP

          !  Special rules for corner 'edges'
          if( j .ne. 1) then
             GridElem(1 ,j,k)%nbrs(swest)   = GridElem(ne,j-1,MODULO(2+k,4)+1)%number
             GridElem(1 ,j,k)%wgtV(swest)    = CornerWgt
             GridElem(1 ,j,k)%wgtP(swest)    = CornerWgt
             GridElem(ne,j,k)%nbrs(seast)   = GridElem(1 ,j-1,MODULO(k  ,4)+1)%number
             GridElem(ne,j,k)%wgtV(seast)    = CornerWgt
             GridElem(ne,j,k)%wgtP(seast)    = CornerWgt
          endif
          if( j .ne. ne) then
             GridElem(1 ,j,k)%nbrs(nwest)   = GridElem(ne,j+1,MODULO(2+k,4)+1)%number
             GridElem(1 ,j,k)%wgtV(nwest)    = CornerWgt
             GridElem(1 ,j,k)%wgtP(nwest)    = CornerWgt
             GridElem(ne,j,k)%nbrs(neast)   = GridElem(1 ,j+1,MODULO(k  ,4)+1)%number
             GridElem(ne,j,k)%wgtV(neast)    = CornerWgt
             GridElem(ne,j,k)%wgtP(neast)    = CornerWgt
          endif
       end do
    end do


    ! ==================================
    ! south edge of 1 / north edge of 5
    ! ==================================

    do i=1,ne
       GridElem(i,1 ,1)%nbrs(south) = GridElem(i,ne,5)%number
       GridElem(i,1 ,1)%wgtV(south)  = EdgeWgtV
       GridElem(i,1 ,1)%wgtP(south)  = EdgeWgtP
       GridElem(i,ne,5)%nbrs(north) = GridElem(i,1 ,1)%number
       GridElem(i,ne,5)%wgtV(north)  = EdgeWgtV
       GridElem(i,ne,5)%wgtP(north)  = EdgeWgtP

       !  Special rules for corner 'edges'
       if( i .ne. 1) then
          GridElem(i,1 ,1)%nbrs(swest)    = GridElem(i-1,ne,5)%number
          GridElem(i,1 ,1)%wgtV(swest)     = CornerWgt
          GridElem(i,1 ,1)%wgtP(swest)     = CornerWgt
          GridElem(i,ne,5)%nbrs(nwest)    = GridElem(i-1,1 ,1)%number
          GridElem(i,ne,5)%wgtV(nwest)     = CornerWgt
          GridElem(i,ne,5)%wgtP(nwest)     = CornerWgt
       endif
       if( i .ne. ne) then
          GridElem(i,1 ,1)%nbrs(seast)    = GridElem(i+1,ne,5)%number
          GridElem(i,1 ,1)%wgtV(seast)     = CornerWgt
          GridElem(i,1 ,1)%wgtP(seast)     = CornerWgt
          GridElem(i,ne,5)%nbrs(neast)    = GridElem(i+1,1 ,1)%number
          GridElem(i,ne,5)%wgtV(neast)     = CornerWgt
          GridElem(i,ne,5)%wgtP(neast)     = CornerWgt
       endif

    end do

    ! ==================================
    ! south edge of 2 / east edge of 5
    ! ==================================

    do i=1,ne
       irev=ne+1-i
       GridElem(i,1 ,2)%nbrs(south) = GridElem(ne,irev,5)%number
       GridElem(i,1 ,2)%wgtV(south)  = EdgeWgtV
       GridElem(i,1 ,2)%wgtP(south)  = EdgeWgtP
       GridElem(ne,i,5)%nbrs(east) = GridElem(irev,1 ,2)%number
       GridElem(ne,i,5)%wgtV(east)  = EdgeWgtV
       GridElem(ne,i,5)%wgtP(east)  = EdgeWgtP

       !  Special rules for corner 'edges'
       if( i .ne. 1) then
          GridElem(i,1 ,2)%nbrs(swest)   = GridElem(ne,irev+1,5)%number
          GridElem(i,1 ,2)%wgtV(swest)    = CornerWgt
          GridElem(i,1 ,2)%wgtP(swest)    = CornerWgt
          GridElem(ne,i,5)%nbrs(seast)   = GridElem(irev+1,1 ,2)%number
          GridElem(ne,i,5)%wgtV(seast)    = CornerWgt
          GridElem(ne,i,5)%wgtP(seast)    = CornerWgt
       endif
       if(i .ne. ne) then
          GridElem(i,1 ,2)%nbrs(seast)   = GridElem(ne,irev-1,5)%number
          GridElem(i,1 ,2)%wgtV(seast)    = CornerWgt
          GridElem(i,1 ,2)%wgtP(seast)    = CornerWgt
          GridElem(ne,i,5)%nbrs(neast)   = GridElem(irev-1,1 ,2)%number
          GridElem(ne,i,5)%wgtV(neast)    = CornerWgt
          GridElem(ne,i,5)%wgtP(neast)    = CornerWgt
       endif
    enddo
    ! ==================================
    ! south edge of 3 / south edge of 5
    ! ==================================

    do i=1,ne
       irev=ne+1-i
       GridElem(i,1,3)%nbrs(south) = GridElem(irev,1,5)%number
       GridElem(i,1,3)%wgtV(south)  = EdgeWgtV
       GridElem(i,1,3)%wgtP(south)  = EdgeWgtP
       GridElem(i,1,5)%nbrs(south) = GridElem(irev,1,3)%number
       GridElem(i,1,5)%wgtV(south)  = EdgeWgtV
       GridElem(i,1,5)%wgtP(south)  = EdgeWgtP

       !  Special rules for corner 'edges'
       if( i .ne. 1) then
          GridElem(i,1,3)%nbrs(swest) = GridElem(irev+1,1,5)%number
          GridElem(i,1,3)%wgtV(swest)  = CornerWgt
          GridElem(i,1,3)%wgtP(swest)  = CornerWgt
          GridElem(i,1,5)%nbrs(swest) = GridElem(irev+1,1,3)%number
          GridElem(i,1,5)%wgtV(swest)  = CornerWgt
          GridElem(i,1,5)%wgtP(swest)  = CornerWgt
       endif
       if(i .ne. ne) then
          GridElem(i,1,3)%nbrs(seast)    = GridElem(irev-1,1,5)%number
          GridElem(i,1,3)%wgtV(seast)     = CornerWgt
          GridElem(i,1,3)%wgtP(seast)     = CornerWgt
          GridElem(i,1,5)%nbrs(seast)    = GridElem(irev-1,1,3)%number
          GridElem(i,1,5)%wgtV(seast)     = CornerWgt
          GridElem(i,1,5)%wgtP(seast)     = CornerWgt
       endif
    end do

    ! ==================================
    ! south edge of 4 / west edge of 5
    ! ==================================

    do i=1,ne
       irev=ne+1-i
       GridElem(i,1,4)%nbrs(south) = GridElem(1,i,5)%number
       GridElem(i,1,4)%wgtV(south)  = EdgeWgtV
       GridElem(i,1,4)%wgtP(south)  = EdgeWgtP
       GridElem(1,i,5)%nbrs(west)  = GridElem(i,1,4)%number
       GridElem(1,i,5)%wgtV(west)   = EdgeWgtV
       GridElem(1,i,5)%wgtP(west)   = EdgeWgtP
       !  Special rules for corner 'edges'
       if( i .ne. 1) then
          GridElem(i,1,4)%nbrs(swest)    = GridElem(1,i-1,5)%number
          GridElem(i,1,4)%wgtV(swest)     = CornerWgt
          GridElem(i,1,4)%wgtP(swest)     = CornerWgt
          GridElem(1,i,5)%nbrs(swest)  = GridElem(i-1,1,4)%number
          GridElem(1,i,5)%wgtV(swest)   = CornerWgt
          GridElem(1,i,5)%wgtP(swest)   = CornerWgt
       endif
       if( i .ne. ne) then
          GridElem(i,1,4)%nbrs(seast) = GridElem(1,i+1,5)%number
          GridElem(i,1,4)%wgtV(seast)  = CornerWgt
          GridElem(i,1,4)%wgtP(seast)  = CornerWgt
          GridElem(1,i,5)%nbrs(nwest)  = GridElem(i+1,1,4)%number
          GridElem(1,i,5)%wgtV(nwest)   = CornerWgt
          GridElem(1,i,5)%wgtP(nwest)   = CornerWgt
       endif
    end do

    ! ==================================
    ! north edge of 1 / south edge of 6
    ! ==================================

    do i=1,ne
       GridElem(i,ne,1)%nbrs(north) = GridElem(i,1 ,6)%number
       GridElem(i,ne,1)%wgtV(north)  = EdgeWgtV
       GridElem(i,ne,1)%wgtP(north)  = EdgeWgtP
       GridElem(i,1 ,6)%nbrs(south) = GridElem(i,ne,1)%number
       GridElem(i,1 ,6)%wgtV(south)  = EdgeWgtV
       GridElem(i,1 ,6)%wgtP(south)  = EdgeWgtP
       !  Special rules for corner 'edges'
       if( i .ne. 1) then
          GridElem(i,ne,1)%nbrs(nwest) = GridElem(i-1,1 ,6)%number
          GridElem(i,ne,1)%wgtV(nwest)  = CornerWgt
          GridElem(i,ne,1)%wgtP(nwest)  = CornerWgt
          GridElem(i,1 ,6)%nbrs(swest) = GridElem(i-1,ne,1)%number
          GridElem(i,1 ,6)%wgtV(swest)  = CornerWgt
          GridElem(i,1 ,6)%wgtP(swest)  = CornerWgt
       endif
       if( i .ne. ne) then
          GridElem(i,ne,1)%nbrs(neast) = GridElem(i+1,1 ,6)%number
          GridElem(i,ne,1)%wgtV(neast)  = CornerWgt
          GridElem(i,ne,1)%wgtP(neast)  = CornerWgt
          GridElem(i,1 ,6)%nbrs(seast) = GridElem(i+1,ne,1)%number
          GridElem(i,1 ,6)%wgtV(seast)  = CornerWgt
          GridElem(i,1 ,6)%wgtP(seast)  = CornerWgt
       endif
    end do

    ! ==================================
    ! north edge of 2 / east edge of 6
    ! ==================================

    do i=1,ne
       GridElem(i,ne,2)%nbrs(north) = GridElem(ne,i,6)%number
       GridElem(i,ne,2)%wgtV(north)  = EdgeWgtV
       GridElem(i,ne,2)%wgtP(north)  = EdgeWgtP
       GridElem(ne,i,6)%nbrs(east)  = GridElem(i,ne,2)%number
       GridElem(ne,i,6)%wgtV(east)   = EdgeWgtV
       GridElem(ne,i,6)%wgtP(east)   = EdgeWgtP
       !  Special rules for corner 'edges'
       if( i .ne. 1) then
          GridElem(i,ne,2)%nbrs(nwest)    = GridElem(ne,i-1,6)%number
          GridElem(i,ne,2)%wgtV(nwest)     = CornerWgt
          GridElem(i,ne,2)%wgtP(nwest)     = CornerWgt
          GridElem(ne,i,6)%nbrs(seast)    = GridElem(i-1,ne,2)%number
          GridElem(ne,i,6)%wgtV(seast)     = CornerWgt
          GridElem(ne,i,6)%wgtP(seast)     = CornerWgt
       endif
       if( i .ne. ne) then
          GridElem(i,ne,2)%nbrs(neast) = GridElem(ne,i+1,6)%number
          GridElem(i,ne,2)%wgtV(neast)   = CornerWgt
          GridElem(i,ne,2)%wgtP(neast)   = CornerWgt
          GridElem(ne,i,6)%nbrs(neast)  = GridElem(i+1,ne,2)%number
          GridElem(ne,i,6)%wgtV(neast)   = CornerWgt
          GridElem(ne,i,6)%wgtP(neast)   = CornerWgt
       endif
    end do

    ! ===================================
    ! north edge of 3 / north edge of 6
    ! ===================================

    do i=1,ne
       irev=ne+1-i
       GridElem(i,ne,3)%nbrs(north) = GridElem(irev,ne,6)%number
       GridElem(i,ne,3)%wgtV(north)  = EdgeWgtV
       GridElem(i,ne,3)%wgtP(north)  = EdgeWgtP
       GridElem(i,ne,6)%nbrs(north) = GridElem(irev,ne,3)%number
       GridElem(i,ne,6)%wgtV(north)  = EdgeWgtV
       GridElem(i,ne,6)%wgtP(north)  = EdgeWgtP
       !  Special rules for corner 'edges'
       if( i .ne. 1) then
          GridElem(i,ne,3)%nbrs(nwest) = GridElem(irev+1,ne,6)%number
          GridElem(i,ne,3)%wgtV(nwest)  = CornerWgt
          GridElem(i,ne,3)%wgtP(nwest)  = CornerWgt
          GridElem(i,ne,6)%nbrs(nwest) = GridElem(irev+1,ne,3)%number
          GridElem(i,ne,6)%wgtV(nwest)  = CornerWgt
          GridElem(i,ne,6)%wgtP(nwest)  = CornerWgt
       endif
       if( i .ne. ne) then
          GridElem(i,ne,3)%nbrs(neast) = GridElem(irev-1,ne,6)%number
          GridElem(i,ne,3)%wgtV(neast)  = CornerWgt
          GridElem(i,ne,3)%wgtP(neast)  = CornerWgt
          GridElem(i,ne,6)%nbrs(neast) = GridElem(irev-1,ne,3)%number
          GridElem(i,ne,6)%wgtV(neast)  = CornerWgt
          GridElem(i,ne,6)%wgtP(neast)  = CornerWgt
       endif
    end do

    ! ===================================
    ! north edge of 4 / west edge of 6
    ! ===================================

    do i=1,ne
       irev=ne+1-i
       GridElem(i,ne,4)%nbrs(north) = GridElem(1,irev,6)%number
       GridElem(i,ne,4)%wgtV(north)  = EdgeWgtV
       GridElem(i,ne,4)%wgtP(north)  = EdgeWgtP
       GridElem(1,i,6)%nbrs(west)   = GridElem(irev,ne,4)%number
       GridElem(1,i,6)%wgtV(west)    = EdgeWgtV
       GridElem(1,i,6)%wgtP(west)    = EdgeWgtP
       !  Special rules for corner 'edges'
       if( i .ne. 1) then
          GridElem(i,ne,4)%nbrs(nwest) = GridElem(1,irev+1,6)%number
          GridElem(i,ne,4)%wgtV(nwest)  = CornerWgt
          GridElem(i,ne,4)%wgtP(nwest)  = CornerWgt
          GridElem(1,i,6)%nbrs(swest)   = GridElem(irev+1,ne,4)%number
          GridElem(1,i,6)%wgtV(swest)    = CornerWgt
          GridElem(1,i,6)%wgtP(swest)    = CornerWgt
       endif
       if( i .ne. ne) then
          GridElem(i,ne,4)%nbrs(neast) = GridElem(1,irev-1,6)%number
          GridElem(i,ne,4)%wgtV(neast)  = CornerWgt
          GridElem(i,ne,4)%wgtP(neast)  = CornerWgt
          GridElem(1,i,6)%nbrs(nwest)   = GridElem(irev-1,ne,4)%number
          GridElem(1,i,6)%wgtV(nwest)    = CornerWgt
          GridElem(1,i,6)%wgtP(nwest)    = CornerWgt
       endif
    end do
    

    ielem = 1                       ! Element counter
    do k=1,6
       do j=1,ne
          do i=1,ne
             GridVertex(ielem)%nbrs       = GridElem(i,j,k)%nbrs
             GridVertex(ielem)%wgtV       = GridElem(i,j,k)%wgtV
             GridVertex(ielem)%wgtP       = GridElem(i,j,k)%wgtP
             GridVertex(ielem)%degree     = GridElem(i,j,k)%degree
             GridVertex(ielem)%number     = GridElem(i,j,k)%number
             GridVertex(ielem)%partition  = 0
             GridVertex(ielem)%SpaceCurve = GridElem(i,j,k)%SpaceCurve
             ielem=ielem+1
          end do
       end do
    end do

#if 0
    if(OutputFiles) then
       close(7)
       close(8)
    endif
#endif

    ! =======================================
    ! Generate cube graph...
    ! =======================================

#if 0
    if(OutputFiles) then
       write(9,*)nelem,2*nelem      ! METIS requires this first line
    endif
#endif

    ! ============================================
    !  Setup the Grid edges (topology independent)
    ! ============================================
    call initgridedge(GridEdge,GridVertex)

    ! ============================================
    !  Setup the Grid edge Indirect addresses
    !          (topology dependent)
    ! ============================================
    nedge = SIZE(GridEdge)
    do i=1,nedge
       call CubeSetupEdgeIndex(GridEdge(i))
    enddo

  end subroutine CubeTopology

  ! =======================================
  ! cube_assemble:
  !
  ! Assemble the cube field element by element
  ! this routine is assumed to be single 
  ! threaded...
  ! =======================================

  function cube_assemble(gbl,fld,elem,par,nelemd,nelem,ielem) result(ierr)
    use element_mod, only : element_t
#ifdef _MPI
    use parallel_mod, only : parallel_t, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_STATUS_SIZE, MPI_REAL8,MPI_TAG
#else
    use parallel_mod, only : parallel_t
#endif
    real (kind=real_kind) :: gbl(:,:,:,:)    ! global output field 
    real (kind=real_kind) :: fld(:,:,:)      ! local model field  
    type (element_t)      :: elem            ! element to assemble 
    type (parallel_t)     :: par             ! parallel structure 
    integer               :: nelemd          ! number of elements on the node
    integer               :: nelem           ! number of elements on the node
    integer               :: ielem           ! local element ctr 
    integer               :: ierr            ! returned error code

    ! Local variables

    integer :: ie,je,face_no
    integer :: ibase,jbase
    integer :: i,j,k
    integer :: elem_number

    integer :: ne1,ne2    ! element dimensions
    integer :: n1,n2      ! gbl face dimensions
    integer :: nface      ! number of faces (must be 6)
    integer :: nlyr       ! number of layers

#if defined(_MPI)
    integer :: ectr       ! global element counter
    integer tag
    integer :: count      ! w/o "::", triggers PGI 3.1 F90 bug 
    integer pe
    integer status(MPI_STATUS_SIZE)
    integer mpi_err
#endif      

    ne1   = SIZE(fld,1)
    ne2   = SIZE(fld,2)
    nlyr  = SIZE(fld,3)

    n1    = SIZE(gbl,1)
    n2    = SIZE(gbl,2)
    nface = SIZE(gbl,3)

    ! =========================
    ! Enforce certain rules...
    ! =========================

    ierr=0

    if (MODULO(n1,ne1) /= 0) then
       ierr=-1
       return
    end if

    if (MODULO(n2,ne2) /= 0) then 
       ierr=-2
       return
    end if

    if (nface /= 6) then
       ierr=-3
       return
    end if

    ! =========================================================
    ! Perform global assembly procedure element by element ...
    ! =========================================================

    if (par%rank==par%root) then

       if (ielem<=nelemd) then
          elem_number = elem%vertex%number

          call convert_gbl_index(elem_number,ie,je,face_no)

          ibase=ie*ne1
          jbase=je*ne2

          do k=1,nlyr
             do j=1,ne2
                do i=1,ne1
                   gbl(i+ibase,j+jbase,face_no,k)=fld(i,j,k)
                end do
             end do
          end do
       end if

#if defined(_MPI)
       if (ielem==nelemd) then
          ectr=nelemd
          do while(ectr<nelem)
             pe    = MPI_ANY_SOURCE
             tag   = MPI_ANY_TAG
             count = ne1*ne2*nlyr
             call MPI_RECV(fld(1,1,1),   &
                  count,        &
                  MPI_REAL8,    &
                  pe,           &
                  tag,          &  
                  par%comm,     &
                  status,       &
                  mpi_err) 

             elem_number = status(MPI_TAG)
             call convert_gbl_index(elem_number,ie,je,face_no)

             ibase=ie*ne1
             jbase=je*ne2

             do k=1,nlyr
                do j=1,ne2
                   do i=1,ne1
                      gbl(i+ibase,j+jbase,face_no,k)=fld(i,j,k)
                   end do
                end do
             end do

             ectr=ectr+1
          end do
       end if

    else

       pe    = par%root
       tag   = elem%vertex%number
       count = ne1*ne2*nlyr
       call MPI_SEND(fld(1,1,1),    &
            count,         &
            MPI_REAL8,     &
            pe,            &
            tag,           &
            par%comm,      &
            mpi_err)
#endif
    end if

  end function cube_assemble

  ! ===================================================================
  ! CubeEdgeCount:
  !
  !  Determine the number of Grid Edges
  !
  ! ===================================================================

  function CubeEdgeCount(ne)  result(nelem)

    implicit none
    integer, intent(in)         :: ne
    integer                     :: nelem

    nelem = nfaces*(ne*ne*nInnerElemEdge - nCornerElemEdge)

  end function CubeEdgeCount

  ! ===================================================================
  ! CubeElemCount:
  !
  !  Determine the number of Grid Elem
  !
  ! ===================================================================

  function CubeElemCount(ne)  result(nelem)

    implicit none
    integer, intent(in)         :: ne
    integer                     :: nelem

    nelem = nfaces*ne*ne

  end function CubeElemCount

  subroutine CubeSetupEdgeIndex(Edge)
    use gridgraph_mod, only : gridedge_t
    use dimensions_mod, only : nv, np
    use control_mod, only : north, south, east, west, neast, seast, swest, nwest
    type (GridEdge_t),target           :: Edge

    integer                            :: nv0,np0,sFace,dFace
    logical                            :: reverse
    integer,allocatable                :: forwardV(:), forwardP(:)
    integer,allocatable                :: backwardV(:), backwardP(:)
    integer                            :: i,ii

    ii=Edge%tail_face
    nv0 = Edge%tail%wgtV(ii)
    np0 = Edge%tail%wgtP(ii)

#ifdef TESTGRID
    allocate(forwardV(nv0))
    allocate(backwardV(nv0))

    allocate(forwardP(np0))
    allocate(backwardP(np0))

    do i=1,nv0
       forwardV(i)  = i
       backwardV(i) = nv0-i+1
    enddo

    do i=1,np0
       forwardP(i)  = i
       backwardP(i) = np0-i+1
    enddo
#endif

    sFace = Edge%tail_face
    dFace = Edge%head_face
    ! Do not reverse the indices
    reverse=.FALSE.

    ! Under special conditions use index reversal
    if(   (SFace == south .AND. dFace == east)  &
         .OR. (sFace == east  .AND. dFace == south) &
         .OR. (sFace == south .AND. dFace == south) &
         .OR. (sFace == north .AND. dFace == north) &
         .OR. (sFace == north .AND. dFace == west)  &
         .OR. (sFace == west  .AND. dFace == north)   ) then
       reverse=.TRUE.
       Edge%reverse=.TRUE.
    endif

#ifdef TESTGRID
    !  Setup the destination indices
    select case(dFace)
    case(east)
       Edge%HeadIndex%ixV=nv
       Edge%HeadIndex%iyV=forwardV

       Edge%HeadIndex%ixP=np
       Edge%HeadIndex%iyP=forwardP
    case(west)
       Edge%HeadIndex%ixV=1
       Edge%HeadIndex%iyV=forwardV

       Edge%HeadIndex%ixP=1
       Edge%HeadIndex%iyP=forwardP
    case(north)
       Edge%HeadIndex%ixV=forwardV
       Edge%HeadIndex%iyV=nv

       Edge%HeadIndex%ixP=forwardP
       Edge%HeadIndex%iyP=np
    case(south)
       Edge%HeadIndex%ixV=forwardV
       Edge%HeadIndex%iyV=1

       Edge%HeadIndex%ixP=forwardP
       Edge%HeadIndex%iyP=1
    case(swest)
       Edge%HeadIndex%ixV=1
       Edge%HeadIndex%iyV=1

       Edge%HeadIndex%ixP=1
       Edge%HeadIndex%iyP=1
    case(seast)
       Edge%HeadIndex%ixV=nv
       Edge%HeadIndex%iyV=1

       Edge%HeadIndex%ixP=np
       Edge%HeadIndex%iyP=1
    case(nwest)
       Edge%HeadIndex%ixV=1
       Edge%HeadIndex%iyV=nv

       Edge%HeadIndex%ixP=1
       Edge%HeadIndex%iyP=np
    case(neast)
       Edge%HeadIndex%ixV=nv
       Edge%HeadIndex%iyV=nv

       Edge%HeadIndex%ixP=np
       Edge%HeadIndex%iyP=np
    case default
       write (*,*) 'SetupEdgeIndex: Error in dFace select statement'
    end select

    ! Setup the source indices
    select case(sFace)
    case(north)
       Edge%TailIndex%ixV=forwardV
       if(reverse) Edge%TailIndex%ixV=backwardV
       Edge%TailIndex%iyV=nv

       Edge%TailIndex%ixP=forwardP
       if(reverse) Edge%TailIndex%ixP=backwardP
       Edge%TailIndex%iyP=np
    case(south)
       Edge%TailIndex%ixV=forwardV
       if(reverse) Edge%TailIndex%ixV=backwardV
       Edge%TailIndex%iyV=1

       Edge%TailIndex%ixP=forwardP
       if(reverse) Edge%TailIndex%ixP=backwardP
       Edge%TailIndex%iyP=1
    case(east)
       Edge%TailIndex%ixV=nv
       Edge%TailIndex%iyV=forwardV
       if(reverse) Edge%TailIndex%iyV=backwardV

       Edge%TailIndex%ixP=np
       Edge%TailIndex%iyP=forwardP
       if(reverse) Edge%TailIndex%iyP=backwardP
    case(west)
       Edge%TailIndex%ixV=1
       Edge%TailIndex%iyV=forwardV
       if(reverse) Edge%TailIndex%iyV=backwardV

       Edge%TailIndex%ixP=1
       Edge%TailIndex%iyP=forwardP
       if(reverse) Edge%TailIndex%iyP=backwardP
    case(swest)
       Edge%TailIndex%ixV=1
       Edge%TailIndex%iyV=1

       Edge%TailIndex%ixP=1
       Edge%TailIndex%iyP=1
    case(seast)
       Edge%TailIndex%ixV=nv
       Edge%TailIndex%iyV=1

       Edge%TailIndex%ixP=np
       Edge%TailIndex%iyP=1
    case(nwest)
       Edge%TailIndex%ixV=1
       Edge%TailIndex%iyV=nv

       Edge%TailIndex%ixP=1
       Edge%TailIndex%iyP=np
    case(neast)
       Edge%TailIndex%ixV=nv
       Edge%TailIndex%iyV=nv

       Edge%TailIndex%ixP=np
       Edge%TailIndex%iyP=np
    case default
       write (*,*) 'SetupEdgeIndex: Error in sFace select statement'
    end select

    deallocate(forwardV)
    deallocate(forwardP)
    deallocate(backwardV)
    deallocate(backwardP)
#endif

  end subroutine CubeSetupEdgeIndex

  subroutine GetLatticeSpacing(spherev,dxv,spherep,dxp)
    use physical_constants, only : rearth
    use dimensions_mod, only : nv, np

    type (spherical_polar_t), intent(in) :: spherev(nv,nv)
    real (kind=real_kind)                :: dxv
    type (spherical_polar_t), intent(in) :: spherep(np,np)
    real (kind=real_kind)                :: dxp

    real (kind=real_kind) xcorner,ycorner,zcorner
    real (kind=real_kind) x,y,z
    real (kind=real_kind) chord
    real (kind=real_kind) theta

    xcorner=COS(spherev(1,1)%lat)*COS(spherev(1,1)%lon)
    ycorner=COS(spherev(1,1)%lat)*SIN(spherev(1,1)%lon)
    zcorner=SIN(spherev(1,1)%lat)

    x=COS(spherev(2,1)%lat)*COS(spherev(2,1)%lon)
    y=COS(spherev(2,1)%lat)*SIN(spherev(2,1)%lon)
    z=SIN(spherev(2,1)%lat)

    chord = SQRT( (xcorner-x)**2 + &
         (ycorner-y)**2 + &
         (zcorner-z)**2 )

    theta = 2.0D0*ASIN(0.50D0*chord)

    dxv   = theta*rearth

    x=COS(spherep(1,1)%lat)*COS(spherep(1,1)%lon)
    y=COS(spherep(1,1)%lat)*SIN(spherep(1,1)%lon)
    z=SIN(spherep(1,1)%lat)

    chord = SQRT( (xcorner-x)**2 + &
         (ycorner-y)**2 + &
         (zcorner-z)**2 )

    theta = 2.0D0*ASIN(0.50D0*chord)

    dxp   = theta*rearth

  end subroutine GetLatticeSpacing

end module cube_mod



