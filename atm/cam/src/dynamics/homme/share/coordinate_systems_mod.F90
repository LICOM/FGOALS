#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module coordinate_systems_mod
  use kinds, only : real_kind
  implicit none
  private

  real (kind=real_kind), public, parameter :: DIST_THRESHOLD= 1.0D-9
  real (kind=real_kind), parameter :: one=1.0D0, two=2.0D0

  type, public :: cartesian2D_t
     sequence
     real(real_kind) :: x             ! x coordinate
     real(real_kind) :: y             ! y coordinate
  end type cartesian2D_t

  type, public :: cartesian3D_t
     sequence
     real(real_kind) :: x             ! x coordinate
     real(real_kind) :: y             ! y coordinate
     real(real_kind) :: z             ! z coordinate
  end type cartesian3D_t

  type, public :: spherical_polar_t
     sequence
     real(real_kind) :: r             ! radius
     real(real_kind) :: lon           ! longitude
     real(real_kind) :: lat           ! latitude
  end type spherical_polar_t


  interface assignment ( = )
     module procedure copy_cart2d
  end interface

  interface operator( == )
     module procedure eq_cart2d
  end interface

  interface distance
     module procedure distance_cart2D
     module procedure distance_cart2D_v
     module procedure distance_cart3D
     module procedure distance_cart3D_v
  end interface

  interface change_coordinates
     module procedure spherical_to_cart_v
     module procedure spherical_to_cart
     module procedure cart_to_spherical_v
     module procedure cart_to_spherical
  end interface


  ! ==========================================
  ! Public Interfaces
  ! ==========================================

  public :: distance
  public :: change_coordinates
  public :: project
  public :: cart2cubedsphere
  ! note: cant make these next two an interface since they only differ by return arg
  public :: cubedsphere2cart
  public :: cubedsphere2cartB
  public :: cart2face
contains

  ! ============================================
  ! copy_cart2d:
  !
  ! Overload assignment operator for cartesian2D_t
  ! ============================================

  subroutine copy_cart2d(cart2,cart1)
    type(cartesian2D_t), intent(out) :: cart2
    type(cartesian2D_t), intent(in)  :: cart1
    cart2%x=cart1%x
    cart2%y=cart1%y
  end subroutine copy_cart2d

  ! ============================================
  ! eq_cart2d:
  !
  ! Overload == operator for cartesian2D_t
  ! ============================================

  function eq_cart2d(cart2,cart1) result(is_same)
    type(cartesian2D_t), intent(in)  :: cart2
    type(cartesian2D_t), intent(in)  :: cart1

    logical :: is_same    

    if (distance(cart1,cart2)<DIST_THRESHOLD) then
       is_same=.true.
    else
       is_same=.false.
    end if

  end function eq_cart2d

  ! ===================================================
  ! distance_cart2D  : scalar version
  ! distance_cart2D_v: vector version
  !
  ! computes distance between cartesian 2D coordinates
  ! ===================================================

  function distance_cart2D(cart1,cart2) result(dist)

    type(cartesian2D_t) :: cart1
    type(cartesian2D_t), optional :: cart2
    real(real_kind)   :: dist

    if (present(cart2)) then
       dist = SQRT((cart1%x-cart2%x)**2 + &
            (cart1%y-cart2%y)**2   )
    else
       dist = SQRT(cart1%x*cart1%x + &
            cart1%y*cart1%y   )
    end if

  end function distance_cart2D

  function distance_cart2D_v(cart1,cart2) result(dist)

    type(cartesian2D_t) :: cart1(:)
    type(cartesian2D_t), optional :: cart2(:)
    real(real_kind)   :: dist(SIZE(cart1))

    integer             :: i

    if (present(cart2)) then
       do i=1,SIZE(cart1)
          dist(i) = SQRT((cart1(i)%x-cart2(i)%x)**2 + &
               (cart1(i)%y-cart2(i)%y)**2   )
       end do
    else
       do i=1,SIZE(cart1)
          dist(i) = SQRT(cart1(i)%x*cart1(i)%x + &
               cart1(i)%y*cart1(i)%y   )
       end do
    end if

  end function distance_cart2D_v


  ! ===================================================
  ! distance_cart3D  : scalar version
  ! distance_cart3D_v: vector version
  ! ===================================================

  function distance_cart3D(cart1,cart2) result(dist)

    type(cartesian3D_t)          :: cart1
    type(cartesian3D_t),optional :: cart2
    real(real_kind)   :: dist

    if (present(cart2)) then
       dist = SQRT((cart1%x-cart2%x)**2 + &
            (cart1%y-cart2%y)**2 + &
            (cart1%z-cart2%z)**2   )
    else
       dist = SQRT(cart1%x*cart1%x + &
            cart1%y*cart1%y + &
            cart1%z*cart1%z   )
    end if
  end function distance_cart3D

  function distance_cart3D_v(cart1,cart2) result(dist)

    type(cartesian3D_t)          :: cart1(:)
    type(cartesian3D_t),optional :: cart2(:)
    real(real_kind)   :: dist(SIZE(cart1))

    integer             :: i

    if (present(cart2)) then
       do i=1,SIZE(cart1)
          dist(i) = SQRT((cart1(i)%x-cart2(i)%x)**2 + &
               (cart1(i)%y-cart2(i)%y)**2 + &
               (cart1(i)%z-cart2(i)%z)**2   )
       end do
    else
       do i=1,SIZE(cart1)
          dist(i) = SQRT(cart1(i)%x*cart1(i)%x + &
               cart1(i)%y*cart1(i)%y + &
               cart1(i)%z*cart1(i)%z   )
       end do
    end if
  end function distance_cart3D_v

  ! ===================================================================
  ! spherical_to_cart_v:
  ! converts spherical polar {lon,lat}  to 3D cartesian {x,y,z}
  ! on unit sphere
  ! ===================================================================

  function spherical_to_cart_v(sphere) result (cart)

    type(spherical_polar_t), intent(in) :: sphere(:)
    type(cartesian3D_t)                 :: cart(SIZE(sphere))

    integer                 :: i

    do i=1,SIZE(sphere)
       cart(i)%x=sphere(i)%r*COS(sphere(i)%lat)*COS(sphere(i)%lon)
       cart(i)%y=sphere(i)%r*COS(sphere(i)%lat)*SIN(sphere(i)%lon)
       cart(i)%z=sphere(i)%r*SIN(sphere(i)%lat)
    end do

  end function spherical_to_cart_v

  ! ===================================================================
  ! spherical_to_cart:
  ! converts spherical polar {lon,lat}  to 3D cartesian {x,y,z}
  ! on unit sphere
  ! ===================================================================

  function spherical_to_cart(sphere) result (cart)

    type(spherical_polar_t), intent(in) :: sphere
    type(cartesian3D_t)                 :: cart

    cart%x=sphere%r*COS(sphere%lat)*COS(sphere%lon)
    cart%y=sphere%r*COS(sphere%lat)*SIN(sphere%lon)
    cart%z=sphere%r*SIN(sphere%lat)

  end function spherical_to_cart

  ! ==========================================================================
  ! cart_to_spherical:
  !
  ! converts 3D cartesian {x,y,z} to spherical polar {lon,lat} 
  ! on unit sphere
  ! ==========================================================================

  function cart_to_spherical_v(cart) result (sphere)
    use physical_constants, only : dd_pi
    type(cartesian3D_t), intent(in) :: cart(:)
    type(spherical_polar_t)         :: sphere(SIZE(cart))

    integer                 :: i

    sphere(:)%r=distance(cart(:))
    sphere(:)%lat=ASIN(cart(:)%z/sphere(:)%r)
    sphere(:)%lon=0

    ! ==========================================================
    ! enforce three facts:
    !
    ! 1) lon at poles is defined to be zero
    ! 
    ! 2) Grid points must be separated by about .01 Meter (on earth)
    !    from pole to be considered "not the pole".
    !
    ! 3) range of lon is { 0<= lon < 2*pi }
    !
    ! ==========================================================

    do i=1,SIZE(cart)

       if (distance(cart(i)) >= DIST_THRESHOLD) then
          sphere(i)%lon=ATAN2(cart(i)%y,cart(i)%x)
          if (sphere(i)%lon<0) then
             sphere(i)%lon=sphere(i)%lon+2*DD_PI
          end if
       end if

    end do

  end function cart_to_spherical_v

  ! scalar version

  function cart_to_spherical(cart) result (sphere)
    use physical_constants, only : dd_pi
    type(cartesian3D_t), intent(in) :: cart         
    type(spherical_polar_t)         :: sphere

    sphere%r=distance(cart)
    sphere%lat=ASIN(cart%z/sphere%r)
    sphere%lon=0

    ! ==========================================================
    ! enforce three facts:
    !
    ! 1) lon at poles is defined to be zero
    ! 
    ! 2) Grid points must be separated by about .01 Meter (on earth)
    !    from pole to be considered "not the pole".
    !
    ! 3) range of lon is { 0<= lon < 2*pi }
    !
    ! ==========================================================

    if (distance(cart) >= DIST_THRESHOLD) then
       sphere%lon=ATAN2(cart%y,cart%x)
       if (sphere%lon<0) then
          sphere%lon=sphere%lon+2*DD_PI
       end if
    end if

  end function cart_to_spherical

  subroutine project(sphere,cart,face_no)         
    use physical_constants, only : dd_pi

    type (spherical_polar_t)             :: sphere(:,:)
    type (cartesian2d_t), intent(in)     :: cart(:,:)   ! assumed to be cartesian coordinates of cube
    ! face, circumscribed to sphere
    integer, intent(in)                  :: face_no

    ! Local

    integer npts
    integer i,j
    real (kind=real_kind) :: r

    npts=SIZE(cart,1)

    do j=1,npts
       do i=1,npts
          sphere(i,j)%r=one     ! r=1 ... unit sphere
          r=SQRT( one + (cart(i,j)%x)**2 + (cart(i,j)%y)**2)
          if (face_no == 1) then
             sphere(i,j)%lat=ASIN((cart(i,j)%y)/r)
             sphere(i,j)%lon=ATAN2(cart(i,j)%x,one)
          else if (face_no == 2) then
             sphere(i,j)%lat=ASIN((cart(i,j)%y)/r)
             sphere(i,j)%lon=ATAN2(one,-cart(i,j)%x)
          else if (face_no == 3) then
             sphere(i,j)%lat=ASIN((cart(i,j)%y)/r)
             sphere(i,j)%lon=ATAN2(-cart(i,j)%x,-one)
          else if (face_no == 4) then
             sphere(i,j)%lat=ASIN((cart(i,j)%y)/r)
             sphere(i,j)%lon=ATAN2(-one,cart(i,j)%x)
          else if (face_no == 5) then
             if (ABS(cart(i,j)%y) > DIST_THRESHOLD .or. ABS(cart(i,j)%x) > DIST_THRESHOLD ) then
                sphere(i,j)%lon=ATAN2(cart(i,j)%x, cart(i,j)%y )
             else
                sphere(i,j)%lon= 0.0D0     ! longitude is meaningless at south pole set to 0.0
             end if
             sphere(i,j)%lat=ASIN(-one/r)
          else if (face_no == 6) then
             if (ABS(cart(i,j)%y) > DIST_THRESHOLD .or. ABS(cart(i,j)%x) > DIST_THRESHOLD ) then
                sphere(i,j)%lon = ATAN2(cart(i,j)%x, -cart(i,j)%y)
             else
                sphere(i,j)%lon= 0.0D0     ! longitude is meaningless at north pole set to 0.0
             end if
             sphere(i,j)%lat=ASIN(one/r)
          end if

          if (sphere(i,j)%lon < 0.0D0) then
             sphere(i,j)%lon=sphere(i,j)%lon + two*DD_PI
          end if

       end do
    end do

  end subroutine project

  ! takes a point on a face of the cube and projects it 
  ! into R^3
  function cubedsphere2cart(cartin,face_no) result(xyz)
    type (cartesian2d_t), intent(in)    :: cartin   ! assumed to be cartesian coordinates of cube
    integer, intent(in)                 :: face_no

    real (kind=real_kind)               :: xyz(3)
    type(cartesian3D_t)                 :: cart

    cart = spherical_to_cart(projectpoint(cartin,face_no))

    xyz(1) = cart%x
    xyz(2) = cart%y
    xyz(3) = cart%z

  end function cubedsphere2cart

  ! takes a point on a face of the cube and projects it 
  ! into R^3
  function cubedsphere2cartB(cartin,face_no) result(cart)
    type (cartesian2d_t), intent(in)    :: cartin   ! assumed to be cartesian coordinates of cube
    integer, intent(in)                 :: face_no

    type(cartesian3D_t)                 :: cart

    cart = spherical_to_cart(projectpoint(cartin,face_no))

  end function cubedsphere2cartB


  function cart2cubedsphere(cart3D,face_no) result(cart)

    type(cartesian3D_t),intent(in) :: cart3d
    integer, intent(in)            :: face_no
    type (cartesian2d_t)           :: cart   

    ! Local
    real (kind=real_kind) :: xony,yonx,yonz,zony,xonz,zonx

    if (face_no == 1) then
       yonx = cart3D%y/cart3D%x
       zonx = cart3D%z/cart3D%x
       cart%x = ATAN(yonx)
       cart%y = ATAN(zonx)
    else if (face_no == 2) then
       xony = cart3D%x/cart3D%y
       zony = cart3D%z/cart3D%y
       cart%x = ATAN(-xony)
       cart%y = ATAN(zony)
    else if (face_no == 3) then
       yonx = cart3D%y/cart3D%x
       zonx = cart3D%z/cart3D%x
       cart%x = ATAN(yonx)
       cart%y = ATAN(-zonx)
    else if (face_no == 4) then
       xony = cart3D%x/cart3D%y
       zony = cart3D%z/cart3D%y
       cart%x = ATAN(-xony)
       cart%y = ATAN(-zony)
    else if (face_no == 6) then
       yonz  = cart3D%y/cart3D%z
       xonz  = cart3D%x/cart3D%z
       cart%x = ATAN(yonz)
       cart%y = ATAN(-xonz)
    else if (face_no == 5) then
       yonz  = cart3D%y/cart3D%z
       xonz  = cart3D%x/cart3D%z
       cart%x = ATAN(-yonz)
       cart%y = ATAN(-xonz)
    end if

  end function cart2cubedsphere



  function cart2face(cart3D) result(face_no)

    type(cartesian3D_t),intent(in) :: cart3d
    integer :: face_no

    real(real_kind) :: x,y,z
    x=cart3d%x
    y=cart3d%y
    z=cart3d%z

    if (y<x .and. y>-x) then      ! x>0
      if (z>x) then
         face_no=6  ! north pole
      else if (z<-x) then
         face_no=5 ! south pole
      else
         face_no = 1
      endif
   else if (y>x .and. y<-x) then  ! x<0
      if (z>-x) then
         face_no=6 ! north pole
      else if (z<x) then
         face_no=5 ! south pole
      else 
         face_no=3
      endif
   else if (y>x .and. y>-x) then  ! y>0
      if (z>y) then
         face_no=6 ! north pole
      else if (z<-y) then
        face_no = 5 ! south pole
      else 
         face_no=2
      endif
   else if (y<x .and. y<-x) then  ! y<0
      if (z>-y) then
         face_no=6 ! north pole
      else if (z<y) then
         face_no=5 ! south pole
      else 
         face_no=4
      endif
   else
      ! y = x.  point is on cube edge, or on face 5 or 6:
      if (z>x) then
         face_no=6  ! north pole
      else if (z<-x) then
         face_no=5 ! south pole
      else
         print *,'x=',x
         print *,'y=',y
         print *,'z=',z
         stop 'ERROR: cart2face: probably called with point on cube edge'
      endif
   endif
   
   end function cart2face


  function projectpoint(cartin,face_no) result(sphere)         
    use physical_constants, only : dd_pi

    type (spherical_polar_t)             :: sphere
    type (cartesian2d_t), intent(in)     :: cartin   ! assumed to be cartesian coordinates of cube
    type (cartesian2d_t)                 :: cart   ! assumed to be cartesian coordinates of cube
    ! face, circumscribed to sphere
    integer, intent(in)                  :: face_no

    ! Local
    integer i,j
    real (kind=real_kind) :: r

    !ASC  This is X and Y and not xhi eta ...

    cart%x = TAN(cartin%x)
    cart%y = TAN(cartin%y)

    sphere%r=one     ! r=1 ... unit sphere
    r=SQRT( one + (cart%x)**2 + (cart%y)**2)
    if (face_no == 1) then
       sphere%lat=ASIN((cart%y)/r)
       sphere%lon=ATAN2(cart%x,one)
    else if (face_no == 2) then
       sphere%lat=ASIN((cart%y)/r)
       sphere%lon=ATAN2(one,-cart%x)
    else if (face_no == 3) then
       sphere%lat=ASIN((cart%y)/r)
       sphere%lon=ATAN2(-cart%x,-one)
    else if (face_no == 4) then
       sphere%lat=ASIN((cart%y)/r)
       sphere%lon=ATAN2(-one,cart%x)
    else if (face_no == 5) then
       if (ABS(cart%y) > DIST_THRESHOLD .or. ABS(cart%x) > DIST_THRESHOLD ) then
          sphere%lon=ATAN2(cart%x, cart%y )
       else
          sphere%lon= 0.0D0     ! longitude is meaningless at south pole set to 0.0
       end if
       sphere%lat=ASIN(-one/r)
    else if (face_no == 6) then
       if (ABS(cart%y) > DIST_THRESHOLD .or. ABS(cart%x) > DIST_THRESHOLD ) then
          sphere%lon = ATAN2(cart%x, -cart%y)
       else
          sphere%lon= 0.0D0     ! longitude is meaningless at north pole set to 0.0
       end if
       sphere%lat=ASIN(one/r)
    end if

    if (sphere%lon < 0.0D0) then
       sphere%lon=sphere%lon + two*DD_PI
    end if

  end function projectpoint

end module coordinate_systems_mod
