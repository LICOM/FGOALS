#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module interpolate_mod
  use kinds, only : real_kind, iulog
  use element_mod, only : element_t
  use dimensions_mod, only : nv, np,ne
  use quadrature_mod, only : quadrature_t, legendre, quad_norm
  use coordinate_systems_mod, only : spherical_polar_t, cartesian2d_t
  use physical_constants, only : DD_PI
  use quadrature_mod, only : quadrature_t, gauss, gausslobatto
  use parallel_mod, only : abortmp
  use cube_mod, only : vmap, convert_gbl_index
  use hybrid_mod,  only : hybrid_t
  use reduction_mod, only : ParallelMin,ParallelMax,ParallelSum


  implicit none
  private
  integer, parameter, public :: MAX_VECVARS=25
  character(len=10), public :: vector_uvars(MAX_VECVARS), vector_vvars(MAX_VECVARS)

  type, public :: interpolate_t
     sequence
     real (kind=real_kind), dimension(:,:), pointer :: Imat  ! P_k(xj)*wj/gamma(k)
     real (kind=real_kind), dimension(:)  , pointer :: rk    ! 1/k
     real (kind=real_kind), dimension(:)  , pointer :: vtemp ! temp results
     real (kind=real_kind), dimension(:)  , pointer :: glp   ! GLL pts (nair)
  end type interpolate_t
  
  type, public :: interpdata_t
     ! Output Interpolation points.  Used to output data on lat-lon (or other grid)
     ! with native element interpolation.  Each element keeps a list of points from the
     ! interpolation grid that are in this element
     type (cartesian2D_t),pointer,dimension(:):: interp_xy      ! element coordinate
     type (cartesian2D_t),pointer,dimension(:):: interp_cube    ! coordinates on face of cube
     integer, pointer,dimension(:)            :: ilat,ilon   ! position of interpolation point in lat-lon grid
     integer                                  :: n_interp
     integer                                  :: nlat  
     integer                                  :: nlon  
     logical                                  :: first_entry = .TRUE.
  end type interpdata_t

  real (kind=real_kind), private :: delta  = 1.0D-9  ! move tiny bit off center to 
  ! avoid landing on element edges
  public :: interp_init
  public :: setup_latlon_interp
  public :: interpolate_scalar
  public :: interpolate_vector
  public :: set_interp_parameter
  public :: get_interp_parameter
  public :: get_interp_gweight
  public :: get_interp_lat
  public :: get_interp_lon
  public :: var_is_vector_uvar, var_is_vector_vvar

!OG
  public :: interpolate_2d
  public :: interpolate_create 
  public :: cube_facepoint


  interface interpolate_scalar
     module procedure interpolate_scalar2d
     module procedure interpolate_scalar3d
  end interface
  interface interpolate_vector
     module procedure interpolate_vector2d
     module procedure interpolate_vector3d
  end interface

  type (interpolate_t), target ::  interp_v,interp_p
 
  ! store the  lat-lon grid
  ! gridtype = 1       equally spaced, including poles (FV scalars output grid)
  ! gridtype = 2       Gauss grid (CAM Eulerian)
  ! gridtype = 3       equally spaced, no poles (FV staggered velocity)
  ! Seven possible history files, last one is inithist and should be native grid
  logical, public, save :: interpolate_analysis(8) = (/.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)
  integer :: nlat,nlon
  real (kind=real_kind), pointer,dimension(:)   :: lat(:),lon(:),gweight(:)
  integer :: gridtype = 1        ! 
  integer :: itype = 1           ! 0 = native high order 
                                 ! 1 = bilinear 

  
  

contains
  subroutine set_interp_parameter(parm_name, value)
    character*(*), intent(in) :: parm_name
    character(len=80) :: msg
    integer :: value,power
    real (kind=real_kind) :: value_target 

    if(parm_name .eq. 'itype') then
       itype=value
    else if(parm_name .eq. 'nlon') then
       nlon=value 
    else if(parm_name .eq. 'nlat') then
       nlat=value
    else if(parm_name.eq. 'gridtype') then
       gridtype=value
    else if(parm_name.eq. 'auto') then
       ! compute recommended nlat,nlon which has slightly higher
       ! resolution than the specifed number of points around equator given in "value"
       ! computed recommended lat-lon grid.  
       ! nlon > peq   peq = points around equator cubed sphere grid 
       ! take nlon power of 2, and at most 1 power of 3
       value_target=value*1.25
       power = nint(.5 +  log( value_target)/log(2d0) )
       power = max(power,7) ! min grid: 64x128
       if ( 3*2**(power-2) > value_target) then
          nlon=3*2**(power-2)   ! use 1 power of 3
       else
          nlon=2**power
       endif
       nlat=nlon/2
       if (gridtype==1) nlat=nlat+1
    else
       write(msg,*) 'Did not recognize parameter named ',parm_name,' in interpolate_mod:set_interp_parameter'
       call abortmp(msg)
    end if
  end subroutine set_interp_parameter
  function get_interp_parameter(parm_name) result(value)
    character*(*), intent(in) :: parm_name
    integer :: value
    character(len=80) :: msg
    if(parm_name .eq. 'itype') then
       value=itype
    else if(parm_name .eq. 'nlon') then
       value=nlon
    else if(parm_name .eq. 'nlat') then
       value=nlat
    else if(parm_name.eq. 'gridtype') then
       value=gridtype
    else
       write(msg,*) 'Did not recognize parameter named ',parm_name,' in interpolate_mod:get_interp_parameter'
       value=-1
       call abortmp(msg)
    end if
    return
  end function get_interp_parameter
  function get_interp_gweight() result(gw)
    real(kind=real_kind) :: gw(nlat)
    gw=gweight
    return
  end function get_interp_gweight
  function get_interp_lat() result(thislat)
    use physical_constants, only : DD_PI
    real(kind=real_kind) :: thislat(nlat)
    thislat=lat*180.0D0/DD_PI
    return
  end function get_interp_lat
  function get_interp_lon() result(thislon)
    use physical_constants, only : DD_PI
    real(kind=real_kind) :: thislon(nlon)
    thislon=lon*180.0D0/DD_PI
    return
  end function get_interp_lon

  subroutine interpolate_create(gquad,interp)
    type (quadrature_t) , intent(in)   :: gquad
    type (interpolate_t), intent(out)  :: interp


    ! Local variables

    integer k,j
    integer npts
    real (kind=real_kind), dimension(:), allocatable :: gamma
    real (kind=real_kind), dimension(:), allocatable :: leg

    npts = size(gquad%points)

    allocate(interp%Imat(npts,npts))
    allocate(interp%rk(npts))
    allocate(interp%vtemp(npts))
    allocate(interp%glp(npts))
    allocate(gamma(npts))
    allocate(leg(npts))

    gamma = quad_norm(gquad,npts)

    do k=1,npts
       interp%rk(k) = 1.0D0/k
       interp%glp(k) = gquad%points(k)    !nair
    end do

    do j=1,npts
       leg=legendre(gquad%points(j),npts-1)
       do k=1,npts
          interp%Imat(j,k)=leg(k)*gquad%weights(j)/gamma(k)
       end do
    end do

    deallocate(gamma)
    deallocate(leg)

  end subroutine interpolate_create

  function interpolate_2d(cart, f, interp, npts, fillvalue) result(fxy)
    integer, intent(in)               :: npts
    type (cartesian2D_t), intent(in)  :: cart
    real (kind=real_kind), intent(in) :: f(npts,npts)
    type (interpolate_t)              :: interp
    real (kind=real_kind)             :: fxy     ! value of f interpolated to (x,y)
    real (kind=real_kind), intent(in), optional :: fillvalue
    ! local variables

    real (kind=real_kind)             :: tmp_1,tmp_2
    real (kind=real_kind)             :: fk0,fk1
    real (kind=real_kind)             :: pk

    integer                           :: l,j,k

    if(present(fillvalue)) then
       if (any(f==fillvalue)) then
          fxy = fillvalue
          return
       endif
    endif


    do l=1,npts,2

       ! Compute Pk(cart%x) for Legendre order 0

       pk = 1.0D0

       fk0=0.0D0
       fk1=0.0D0
       do j=1,npts
          fk0 = fk0 + interp%Imat(j,1)*f(j,l  )
          fk1 = fk1 + interp%Imat(j,1)*f(j,l+1)
       end do
       interp%vtemp(l  ) = pk*fk0
       interp%vtemp(l+1) = pk*fk1

       ! Compute Pk(cart%x) for Legendre order 1

       tmp_2 = pk
       pk    = cart%x

       fk0=0.0D0
       fk1=0.0D0
       do j=1,npts
          fk0 = fk0 + interp%Imat(j,2)*f(j,l  )
          fk1 = fk1 + interp%Imat(j,2)*f(j,l+1)
       end do
       interp%vtemp(l  ) = interp%vtemp(l  ) + pk*fk0
       interp%vtemp(l+1) = interp%vtemp(l+1) + pk*fk1

       ! Compute Pk(cart%x) for Legendre order 2 to npts-1

       do k = 2,npts-1

          tmp_1  = tmp_2
          tmp_2  = pk
          pk = ( (2*k-1)*cart%x*tmp_2 - (k-1)*tmp_1 )*interp%rk(k)

          fk0=0.0D0
          fk1=0.0D0
          do j=1,npts
             fk0 = fk0 + interp%Imat(j,k+1)*f(j,l  )
             fk1 = fk1 + interp%Imat(j,k+1)*f(j,l+1)
          end do
          interp%vtemp(l  ) = interp%vtemp(l  ) + pk*fk0
          interp%vtemp(l+1) = interp%vtemp(l+1) + pk*fk1

       end do

    end do

    ! Compute Pk(cart%y) for Legendre order 0

    pk = 1.0

    fk0 = 0.0D0
    do j=1,npts
       fk0 = fk0 + interp%Imat(j,1)*interp%vtemp(j)
    end do
    fxy = pk*fk0

    ! Compute Pk(cart%y) for Legendre order 1

    tmp_2 = pk
    pk    = cart%y

    fk0=0.0D0
    do j=1,npts
       fk0 = fk0 + interp%Imat(j,2)*interp%vtemp(j)
    end do
    fxy = fxy + pk*fk0

    ! Compute Pk(cart%y) for Legendre order 2, npts-1

    do k = 2,npts-1
       tmp_1  = tmp_2
       tmp_2  = pk
       pk = ( (2*k-1)*cart%y*tmp_2 - (k-1)*tmp_1 )*interp%rk(k)

       fk0 = 0.0D0
       do j=1,npts
          fk0 = fk0 + interp%Imat(j,k+1)*interp%vtemp(j)
       end do

       fxy = fxy + pk*fk0

    end do

  end function interpolate_2d

  !===============================
  !(Nair) Bilinear interpolation for every GLL grid cell
  !===============================

  function interpol_bilinear(cart, f, interp, npts, fillvalue) result(fxy)

    integer, intent(in)               :: npts
    type (cartesian2D_t), intent(in)  :: cart
    real (kind=real_kind), intent(in) :: f(npts,npts)
    type (interpolate_t)              :: interp
    real (kind=real_kind)             :: fxy     ! value of f interpolated to (x,y)
    real (kind=real_kind), intent(in), optional :: fillvalue
    ! local variables

    real (kind=real_kind)             :: xoy(npts)
    real (kind=real_kind)             :: p,q,xp,yp ,y4(4)

    integer                           :: l,j,k, ii, jj, na,nb,nm

    xp = cart%x
    yp = cart%y

    xoy(:) = interp%glp(:)

    ! Search index along "x"  (bisection method)

    na = 1
    nb = npts
    do
       if  ((nb-na) <=  1)  exit
       nm = (nb + na)/2
       if (xp  >  xoy(nm)) then
          na = nm
       else
          nb = nm
       endif
    enddo
    ii = na

    ! Search index along "y"

    na = 1
    nb = npts
    do
       if  ((nb-na) <=  1)  exit
       nm = (nb + na)/2
       if (yp  >  xoy(nm)) then
          na = nm
       else
          nb = nm
       endif
    enddo
    jj = na

    ! GLL cell containing (xp,yp)

    y4(1) = f(ii,jj)
    y4(2) = f(ii+1,jj)
    y4(3) = f(ii+1,jj+1)
    y4(4) = f(ii,jj+1)

    if(present(fillvalue)) then
       if (any(y4==fillvalue)) then
          fxy = fillvalue
          return
       endif
    endif
       p = (xp - xoy(ii))/(xoy(ii+1) - xoy(ii))
       q = (yp - xoy(jj))/(xoy(jj+1) - xoy(jj))

       fxy = (1.0D0 - p)*(1.0D0 - q)* y4(1) + p*(1.0D0 - q) * y4(2)   &
            + p*q* y4(3) + (1.0D0 - p)*q * y4(4)

  end function interpol_bilinear




  !===============================================
  !  Extracting integer part of the division (a/b) 
  !===============================================

  function integer_part(a, b) result(ii)

    real (kind=real_kind), intent(in) :: a,b
    integer                           :: ii     ! value of f interpolated to (x,y)

    ! local variables

    real (kind=real_kind)             :: a1,b1,tol

    tol = 1.0D-12

    a1 = abs(a)
    b1 = abs((b - a1)/b)

    if (b1 < tol )  then
       if (a < 0.0D0 )  then
          ii = -1
       else
          ii =  1
       endif
    else
       ii = int(a/b)
    endif

  end function integer_part





  !================================================
  !  (Nair) Cube face index and local coordinates
  !================================================

  subroutine cube_facepoint(sphere,ne,cart,cube,ie,je,face_no)
    type (spherical_polar_t), intent (in) :: sphere
    integer             , intent(in)      :: ne
    type (cartesian2D_t), intent(out)     :: cart,cube
    integer             , intent(out)     :: ie,je,face_no

    real (kind=real_kind) :: x,y,z
    real (kind=real_kind) :: lat,lon
    real (kind=real_kind) :: x1,x2
    real (kind=real_kind) :: dx,eps
    real (kind=real_kind) :: pi,twopi, pi2, pi4
    real (kind=real_kind) :: xp,yp, pm, xl,yt

    integer               :: ix,iy,iz 


    lat = sphere%lat
    lon = sphere%lon

    x = COS(lat)*COS(lon)
    y = COS(lat)*SIN(lon)
    z = SIN(lat)  

    dx = (0.5D0*DD_PI)/ne
    eps = 0.0D0 

    xl = lon
    yt = lat

    pi = DD_PI
    twopi = 2.0D0 * pi
    pi2   = pi * 0.5D0
    pi4   = pi * 0.25D0

    pm = max(abs(x),abs(y),abs(z))

    ix =  integer_part(x,pm)
    iy =  integer_part(y,pm)
    iz =  integer_part(z,pm)

    ! For Homme grids

    if (iz  ==  1) then
       face_no    = 6
    elseif (iz  == -1) then
       face_no    = 5
    elseif ((ix == 1).and.(iy /= 1)) then
       face_no    = 1
    elseif ((ix == -1).and.(iy /= -1)) then
       face_no    = 3
    elseif ((iy == 1).and.(ix /= -1)) then
       face_no    = 2
    elseif ((iy == -1).and.(ix /=  1)) then
       face_no    = 4
    endif

    ! ip = face_no   

    if (face_no  == 1) then

       if (xl >= (twopi - pi4)) xp = xl - twopi
       if (xl <=  pi4) xp = xl

       yp = atan(tan(yt)/cos(xp))

    elseif  (face_no == 2) then
       xp = xl - pi2
       yp = atan(tan(yt)/cos(xp))

    elseif  (face_no == 3) then
       xp = xl - pi
       yp = atan(tan(yt)/cos(xp))

    elseif  (face_no == 4) then
       xp = xl - pi *1.5D0         
       yp = atan(tan(yt)/cos(xp))

    elseif  (face_no == 6) then
       xp = atan(sin(xl)/tan(yt))
       yp = atan(-cos(xl)/tan(yt))

    elseif  (face_no == 5) then
       xp = atan(-sin(xl)/tan(yt))
       yp = atan(-cos(xl)/tan(yt))

    endif

    ! coordinates on the cube:
    cube%x = xp
    cube%y = yp
    

    x1 = xp + 0.25D0*DD_PI
    x2 = yp + 0.25D0*DD_PI

!    ie = INT(ABS(x1-eps)/dx) + 1
!    je = INT(ABS(x2-eps)/dx) + 1

    
    ie = INT(ABS(x1)/dx) 
    je = INT(ABS(x2)/dx)
    ! if we are exactly on an element edge, we can put the point in
    ! either the ie or ie+1 element, EXCEPT if ie==ne.  
    if ( ABS(x1) < ne*dx  ) ie=ie+1
    if ( ABS(x2) < ne*dx  ) je=je+1

    if (ie>ne .or. je>ne) then
       write(iulog,*)'ERROR: ',ie,je,ne
       write(iulog,*)'lat,lon=',lat,lon
       write(iulog,*)'x,y,z=',x,y,z
       write(iulog,*)'face no=',face_no
       write(iulog,*)x1,x2,x1/dx,x2/dx
       call abortmp('interpolate_mod: bad argument')
    endif

! bug fix MT 1/2009.  This was creating a plotting error at
! the row of elements in iface=2 at 50 degrees (NE=16 128x256 lat/lon grid)
! For point on element edge, we can have ie=2, but x1=dx 
! but if ie>1, we must execute this statement.
! The only time we can skip this statement is if ie=1, but then
! the statement has no effect, so lets never skip it:
!    if (x1 > dx ) then 
       x1 = x1 - dble(ie-1)*dx  
!    endif

    x1 = 2.0D0*(x1/dx)-1.0D0

!    if (x2 > dx ) then    ! removed MT 1/2009, see above
       x2 = x2 - dble(je-1)*dx  
!    endif

    x2 = 2.0D0*(x2/dx)-1.0D0

    ! coordinates within an element [-1,1]
    cart%x = x1
    cart%y = x2

  end subroutine cube_facepoint

  subroutine interp_init()
    type (quadrature_t)   :: gp

    gp = gausslobatto(nv)
    call interpolate_create(gp,interp_v)
    gp = gausslobatto(np)
    call interpolate_create(gp,interp_p)
  end subroutine interp_init



  subroutine setup_latlon_interp(elem,interpdata,hybrid, nets,nete)
    !
    ! initialize interpolation data structures to interpolate to a lat-lon grid
    ! 
    ! 


    implicit none
    type (element_t)     , intent(in), target :: elem(:)
    type (interpdata_t)  , intent(out) :: interpdata(:)
    type (hybrid_t)      , intent(in)            :: hybrid 
    integer, intent(in) :: nets,nete


    ! local
    integer i,j,ii,count,count_total,n_interp,ie,je,face_no,count_max
    integer ie2,je2,face_no2,ngrid,plat
    real (kind=real_kind)    ::  countx
    real (kind=real_kind)    ::  dp,latdeg(nlat+1),clat(nlat+1),w(nlat+1),w_staggered(nlat)
    real (kind=real_kind)    ::  clat_staggered(nlat),latdeg_st(nlat)

    type (spherical_polar_t) :: sphere
    type (cartesian2D_t) :: cart,cube
    type (quadrature_t)   :: gp


    logical,save                                  :: first_time=.TRUE.

    if (hybrid%masterthread) then
       write(iulog,'(a,i4,a,i4,a)') 'Initializing ',nlat,' x ',nlon,' lat-lon interpolation grid: '
    endif

    do ii=nets,nete
       interpdata(ii)%n_interp=0  ! reset counter
    enddo

    if(first_time)then
       NULLIFY(lat)
       NULLIFY(gweight)
       NULLIFY(lon)
       first_time=.FALSE.
    end if

    if (associated(lat))then
       if(size(lat)>0)deallocate(lat)
    endif
    if (associated(gweight))then
       if(size(gweight)>0)deallocate(gweight)
    endif

    if (associated(lon))then
       if(size(lon)>0)deallocate(lon)
    endif

    allocate(lat(nlat))
    allocate(gweight(nlat))
    allocate(lon(nlon))
    call interp_init()
    gweight=0
    do i=1,nlon
       lon(i)=2*dd_pi*(i-1)/nlon
    enddo
    if (gridtype==1) then
       do j=1,nlat
          lat(j) = -dd_pi/2 + dd_pi*(j-1)/(nlat-1)
       end do
       plat=nlat
    endif
    if (gridtype==2) then
       gp=gauss(nlat)
       do j=1,nlat
          lat(j) = asin(gp%points(j))
          gweight(j) = gp%weights(j)
       end do
    endif
    if (gridtype==3) then
       do j=1,nlat
          lat(j) = -dd_pi/2 + dd_pi*(j-.5d0)/nlat
       end do
       plat=nlat+1
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! compute weighs for equally spaced grids, based on CAM-FV formulas
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if (gridtype==1 .or. gridtype==3) then
       ! gridtype=1    plat=nlat    gweight(1:nlat)=w(1:plat)
       ! gridtype=3    plat=nlat+1  gweight(1:nlat)=w_staggered(1:plat-1)
       
       ! L-R dynamics uses a regular latitude distribution (not gausian).
       ! The algorithm below is a bastardized version of LSM: map.F.
       dp = 180d0/(plat-1)
       do j = 1, plat
          latdeg(j) = -90d0 + (j-1)*dp
          clat(j) = latdeg(j)*dd_pi/180d0
       end do
       
       ! Calculate latitudes for the staggered grid
       
       do j = 1, plat-1
          clat_staggered(j) = (clat(j) + clat(j+1)) / 2
          latdeg_st     (j) = clat_staggered(j)*180d0/dd_pi
       end do
       
       ! Weights are defined as cos(phi)*(delta-phi)
       ! For a sanity check, the sum of w across all lats should be 2, or 1 across
       ! half of the latitudes.
       
       do j = 2, plat-1
          w(j) = sin(clat_staggered(j)) - sin(clat_staggered(j-1))
       end do
       w(1) = sin(clat_staggered(1)) + 1
       w(plat) = w(1)
       
       ! with nlat=2048, this error was 4e-16
       if (abs(sum(w(1:plat)) - 2) > 1e-8) then
          write(iulog,*) 'interpolate_mod: w weights do not sum to 2. sum=',sum(w(1:plat))
          call abortmp('interpolate_mod: weights do not sum to 2.')
       end if
       
       dp = dd_pi / (plat-1)
       do j = 1, plat-1
          w_staggered(j) = sin(clat(j+1)) - sin(clat(j))
       end do
       
       
       if (abs(sum(w_staggered(1:plat-1)) - 2) > 1e-8) then
          write(iulog,*) 'interpolate_mod: staggered weights do not sum to 2. sum=',sum(w_staggered(1:plat-1))
          call abortmp('interpolate_mod: weights do not sum to 2.')
       end if
       
       if (gridtype==1) then
          gweight(1:nlat)=w(1:plat)
       endif
       if (gridtype==3) then
          gweight(1:nlat)=w_staggered(1:plat-1)
       endif
    endif
           
 

    ! go through once, counting the number of points on each element
    sphere%r=1
    do j=1,nlat
       do i=1,nlon
          sphere%lat=lat(j)
          sphere%lon=lon(i)
          call cube_facepoint(sphere,ne,cart,cube,ie,je,face_no)
          
          ! the sphere point belongs to the element (ie,je) on face = face_no.
          ! do I own this element?
          do ii=nets,nete             
             call convert_gbl_index(elem(ii)%vertex%number,ie2,je2,face_no2)
             ie2=ie2+1
             je2=je2+1
             if (face_no == face_no2 .and. ie==ie2 .and. je==je2 ) then
                interpdata(ii)%n_interp = interpdata(ii)%n_interp + 1
                exit
             endif
          enddo
       enddo
       if (hybrid%masterthread) then
          if (mod(j,64).eq.1) then
             print *,'finished latitude ',j,' of ',nlat
          endif
       endif
    enddo
    ! check if every point in interpolation grid was claimed by an element:
    countx=sum(interpdata(nets:nete)%n_interp)
    count_total = ParallelSum(countx,hybrid)
    if (count_total /= nlat*nlon ) then
       if(hybrid%masterthread) then
          write(iulog,*)__FILE__,__LINE__,count_total, nlat, nlon
       end if
       call abortmp('Error setting up interpolation grid count_total<>nlat*nlon') 
    endif

    countx=maxval(interpdata(nets:nete)%n_interp)
    count_max = ParallelMax(countx,hybrid)

    if (hybrid%masterthread) then
       write(iulog,'(a,f8.1)') 'Average number of interpolation points per element: ',count_total/real(6*ne*ne)
       write(iulog,'(a,i6)') 'Maximum number of interpolation points on any element: ',count_max
    endif

    ! allocate storage
    do ii=nets,nete
       ngrid = interpdata(ii)%n_interp
       if(interpdata(ii)%first_entry)then
          NULLIFY(interpdata(ii)%interp_xy)
          NULLIFY(interpdata(ii)%interp_cube)
          NULLIFY(interpdata(ii)%ilat)
          NULLIFY(interpdata(ii)%ilon)

          interpdata(ii)%first_entry=.FALSE.
       endif
       if(associated(interpdata(ii)%interp_xy))then
          if(size(interpdata(ii)%interp_xy)>0)deallocate(interpdata(ii)%interp_xy)
       endif
       if(associated(interpdata(ii)%interp_cube))then
          if(size(interpdata(ii)%interp_cube)>0)deallocate(interpdata(ii)%interp_cube)
       endif
       if(associated(interpdata(ii)%ilat))then
          if(size(interpdata(ii)%ilat)>0)deallocate(interpdata(ii)%ilat)
       endif

       if (associated(interpdata(ii)%ilon))then
          if(size(interpdata(ii)%ilon)>0)deallocate(interpdata(ii)%ilon)
       endif
       allocate(interpdata(ii)%interp_xy( ngrid ) )
       allocate(interpdata(ii)%interp_cube( ngrid ) )
       allocate(interpdata(ii)%ilat( ngrid ) )
       allocate(interpdata(ii)%ilon( ngrid ) )
       interpdata(ii)%n_interp=0  ! reset counter
    enddo

    ! now go through the list again, adding the coordinates
    ! if this turns out to be slow, then it can be done in the loop above
    ! but we have to allocate and possibly resize the interp_xy() array.  
    do j=1,nlat
       do i=1,nlon
          sphere%lat=lat(j)
          sphere%lon=lon(i)
          call cube_facepoint(sphere,ne,cart,cube,ie,je,face_no)

          ! the sphere point belongs to the element (ie,je) on face = face_no.
          ! do I own this element?
          ! the sphere point belongs to the element (ie,je) on face = face_no.
          ! do I own this element?
          do ii=nets,nete
             call convert_gbl_index(elem(ii)%vertex%number,ie2,je2,face_no2)
             ie2=ie2+1
             je2=je2+1
             if (face_no == face_no2 .and. ie==ie2 .and. je==je2 ) then
                ngrid = interpdata(ii)%n_interp + 1
                interpdata(ii)%n_interp = ngrid
                interpdata(ii)%interp_xy( ngrid ) = cart
                interpdata(ii)%interp_cube( ngrid ) = cube
                interpdata(ii)%ilon( ngrid ) = i
                interpdata(ii)%ilat( ngrid ) = j
             endif
          enddo
       enddo
    enddo

  end subroutine setup_latlon_interp



  ! =======================================
  ! interpolate_scalar
  !
  ! Interpolate a scalar field given in an element (fld_cube) to the points in 
  ! interpdata%interp_xy(i), i=1 .. interpdata%n_interp.  
  !
  ! Note that it is possible the given element contains none of the interpolation points
  ! =======================================
  subroutine interpolate_scalar2d(interpdata,fld_cube,npts,fld, fillvalue)
    integer                  ::  npts
    real (kind=real_kind)    ::  fld_cube(npts,npts) ! cube field
    real (kind=real_kind)    ::  fld(:)          ! field at new grid lat,lon coordinates
    type (interpdata_t)         ::  interpdata
    real (kind=real_kind), intent(in), optional :: fillvalue
    ! Local variables
    type (interpolate_t), pointer  ::  interp          ! interpolation structure

    integer :: ne

    integer :: i

    type (cartesian2D_t) :: cart

    if (npts==nv) then
       interp => interp_v
    else if (npts==np) then
       interp => interp_p
    else
       call abortmp('Error in interpolate_scalar(): must be called with p or v grid data') 
    endif

       ! Choice for Native (high-order) or Bilinear interpolations
    if(present(fillvalue)) then
       if (itype == 0) then
          do i=1,interpdata%n_interp
             fld(i)=interpolate_2d(interpdata%interp_xy(i),fld_cube,interp,npts,fillvalue)
          end do
       elseif (itype == 1) then
          do i=1,interpdata%n_interp
             fld(i)=interpol_bilinear(interpdata%interp_xy(i),fld_cube,interp,npts,fillvalue)
          end do
       end if
    else
       if (itype == 0) then
          do i=1,interpdata%n_interp
             fld(i)=interpolate_2d(interpdata%interp_xy(i),fld_cube,interp,npts)
          end do
       elseif (itype == 1) then
          do i=1,interpdata%n_interp
             fld(i)=interpol_bilinear(interpdata%interp_xy(i),fld_cube,interp,npts)
          end do
       end if
    endif


  end subroutine interpolate_scalar2d
  subroutine interpolate_scalar3d(interpdata,fld_cube,npts,nlev,fld, fillvalue)
    integer , intent(in)                 ::  npts, nlev
    real (kind=real_kind)    ::  fld_cube(npts,npts,nlev) ! cube field
    real (kind=real_kind)    ::  fld(:,:)          ! field at new grid lat,lon coordinates
    type (interpdata_t)         ::  interpdata
    real (kind=real_kind), intent(in), optional :: fillvalue
    ! Local variables
    type (interpolate_t), pointer  ::  interp          ! interpolation structure

    integer :: ne

    integer :: i, k

    type (cartesian2D_t) :: cart

    if (npts==nv) then
       interp => interp_v
    else if (npts==np) then
       interp => interp_v
    else
       call abortmp('Error in interpolate_scalar(): must be called with p or v grid data') 
    endif

    ! Choice for Native (high-order) or Bilinear interpolations
    if(present(fillvalue)) then
       if (itype == 0) then
          do k=1,nlev
             do i=1,interpdata%n_interp
                fld(i,k)=interpolate_2d(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts,fillvalue)
             end do
          end do
       elseif (itype == 1) then
          do k=1,nlev
             do i=1,interpdata%n_interp
                fld(i,k)=interpol_bilinear(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts,fillvalue)
             end do
          end do
       endif
    else
       if (itype == 0) then
          do k=1,nlev
             do i=1,interpdata%n_interp
                fld(i,k)=interpolate_2d(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts)
             end do
          end do
       elseif (itype == 1) then
          do k=1,nlev
             do i=1,interpdata%n_interp
                fld(i,k)=interpol_bilinear(interpdata%interp_xy(i),fld_cube(:,:,k),interp,npts)
             end do
          end do
       else 
          write(iulog,*) itype
          call abortmp("wrong interpolation type")
       endif
    endif
  end subroutine interpolate_scalar3d



  ! =======================================
  ! interpolate_vector
  !
  ! Interpolate a vector field given in an element (fld_cube)
  ! to the points in interpdata%interp_xy(i), i=1 .. interpdata%n_interp.  
  !
  ! input_coords = 0    fld_cube given in lat-lon
  ! input_coords = 1    fld_cube given in contravariant
  !
  ! Note that it is possible the given element contains none of the interpolation points
  ! =======================================
  subroutine interpolate_vector2d(interpdata,Dinv,face_no, fld_cube,npts,fld,input_coords, fillvalue)
    integer                  ::  npts
    real (kind=real_kind)    ::  fld_cube(npts,npts,2) ! vector field
    real (kind=real_kind)    ::  fld(:,:)          ! field at new grid lat,lon coordinates
    type (interpdata_t)      ::  interpdata
    integer, intent(in)      ::  face_no
    real (kind=real_kind), intent(in)    ::  Dinv(:,:,:,:)   ! derivative of gnomonic mapping
    real (kind=real_kind), intent(in), optional :: fillvalue
    integer                  ::  input_coords

    ! Local variables
    real (kind=real_kind)    ::  fld_contra(npts,npts,2) ! vector field
    type (interpolate_t), pointer  ::  interp          ! interpolation structure

    real (kind=real_kind)    ::  v1,v2
    real (kind=real_kind)    ::  D(2,2)   ! derivative of gnomonic mapping
    integer :: ie,je


    integer :: i,j

    type (cartesian2D_t) :: cart

    if(present(fillvalue)) then
       if (any(fld_cube==fillvalue)) then
          fld = fillvalue
          return
       end if
    end if
    
    if (input_coords==0 ) then
       ! convert to contra
       do j=1,npts
          do i=1,npts
             ! latlon->contra
             fld_contra(i,j,1) = Dinv(1,1,i,j)*fld_cube(i,j,1) + Dinv(1,2,i,j)*fld_cube(i,j,2)
             fld_contra(i,j,2) = Dinv(2,1,i,j)*fld_cube(i,j,1) + Dinv(2,2,i,j)*fld_cube(i,j,2)
          enddo
       enddo
    else
       fld_contra=fld_cube
    endif


    if (npts==nv) then
       interp => interp_v
    else if (npts==np) then
       call abortmp('Error in interpolate_vector(): input must be on velocity grid')
    endif


       ! Choice for Native (high-order) or Bilinear interpolations

    if (itype == 0) then
       do i=1,interpdata%n_interp
          fld(i,1)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,1),interp,npts)
          fld(i,2)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,2),interp,npts)
       end do
    elseif (itype == 1) then
       do i=1,interpdata%n_interp
          fld(i,1)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,1),interp,npts)
          fld(i,2)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,2),interp,npts)
       end do
    else 
       write(iulog,*) itype
       call abortmp("wrong interpolation type")
    endif
    do i=1,interpdata%n_interp
       ! convert fld from contra->latlon

       ! compute D(:,:) at the point elem%interp_xy(i)
       ! note that vmap takes cube-face coordinates, not element local coordinates
       call vmap(D,interpdata%interp_cube(i)%x,interpdata%interp_cube(i)%y,face_no)

       v1 = fld(i,1)
       v2 = fld(i,2)

       fld(i,1)=D(1,1)*v1 + D(1,2)*v2  
       fld(i,2)=D(2,1)*v1 + D(2,2)*v2  
    end do

  end subroutine interpolate_vector2d
  ! =======================================
  ! interpolate_vector
  !
  ! Interpolate a vector field given in an element (fld_cube)
  ! to the points in interpdata%interp_xy(i), i=1 .. interpdata%n_interp.  
  !
  ! input_coords = 0    fld_cube given in lat-lon
  ! input_coords = 1    fld_cube given in contravariant
  !
  ! Note that it is possible the given element contains none of the interpolation points
  ! =======================================
  subroutine interpolate_vector3d(interpdata,Dinv,face_no, fld_cube,npts,nlev,fld,input_coords, fillvalue)
    type (interpdata_t),intent(in)       ::  interpdata
    real (kind=real_kind), intent(in)    ::  Dinv(:,:,:,:)   ! derivative of gnomonic mapping
    integer, intent(in)                  ::  npts, nlev, face_no
    real (kind=real_kind), intent(in)    ::  fld_cube(npts,npts,2,nlev) ! vector field
    real (kind=real_kind), intent(out)   ::  fld(:,:,:)          ! field at new grid lat,lon coordinates
    real (kind=real_kind), intent(in),optional :: fillvalue
    integer, intent(in)                  ::  input_coords

    ! Local variables
    real (kind=real_kind)    ::  fld_contra(npts,npts,2,nlev) ! vector field
    type (interpolate_t), pointer  ::  interp          ! interpolation structure

    real (kind=real_kind)    ::  v1,v2
    real (kind=real_kind)    ::  D(2,2)   ! derivative of gnomonic mapping
    integer :: ie,je


    integer :: i,j,k

    type (cartesian2D_t) :: cart
    if(present(fillvalue)) then
       if (any(fld_cube==fillvalue)) then
          fld = fillvalue
          return
       end if
    end if
    if (input_coords==0 ) then
       ! convert to contra
       do k=1,nlev
          do j=1,npts
             do i=1,npts
                ! latlon->contra
                fld_contra(i,j,1,k) = Dinv(1,1,i,j)*fld_cube(i,j,1,k) + Dinv(1,2,i,j)*fld_cube(i,j,2,k)
                fld_contra(i,j,2,k) = Dinv(2,1,i,j)*fld_cube(i,j,1,k) + Dinv(2,2,i,j)*fld_cube(i,j,2,k)
             enddo
          enddo
       end do
    else
       fld_contra=fld_cube
    endif

    if (npts==nv) then
       interp => interp_v
    else if (npts==np) then
       call abortmp('Error in interpolate_vector(): input must be on velocity grid')
    endif


       ! Choice for Native (high-order) or Bilinear interpolations

    if (itype == 0) then
       do k=1,nlev
          do i=1,interpdata%n_interp
             fld(i,k,1)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp,npts)
             fld(i,k,2)=interpolate_2d(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp,npts)
          end do
       end do
    elseif (itype == 1) then
       do k=1,nlev
          do i=1,interpdata%n_interp
             fld(i,k,1)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,1,k),interp,npts)
             fld(i,k,2)=interpol_bilinear(interpdata%interp_xy(i),fld_contra(:,:,2,k),interp,npts)
          end do
       end do
    else 
       call abortmp("wrong interpolation type")
    endif


    do i=1,interpdata%n_interp
       ! compute D(:,:) at the point elem%interp_cube(i)
       ! note that vmap takes cube-face coordinates, not element local coordinates
       call vmap(D(:,:),interpdata%interp_cube(i)%x,interpdata%interp_cube(i)%y,face_no)
       do k=1,nlev
          ! convert fld from contra->latlon
          v1 = fld(i,k,1)
          v2 = fld(i,k,2)

          fld(i,k,1)=D(1,1)*v1 + D(1,2)*v2  
          fld(i,k,2)=D(2,1)*v1 + D(2,2)*v2  
       end do
    end do
  end subroutine interpolate_vector3d
  function var_is_vector_uvar(name) 
    character(len=*), intent(in) :: name
    integer :: i, var_is_vector_uvar, null_index
    
    var_is_vector_uvar=0
    null_index=0
    do i=1,MAX_VECVARS
       if(trim(vector_uvars(i)).eq. '') then
          null_index=i
          exit
       endif
       if(trim(vector_uvars(i)).eq.name) then
          var_is_vector_uvar=i
          exit
       end if
    end do
#if 0
DISABLED: breaks in many cases, like UV
    if (var_is_vector_uvar==0) then
       ! default rules: if variable starts with U and was not found, add it to the list:
       if (name(1:1).eq.'U') then 
          if (null_index==0) then
             call abortmp("Error: MAX_VECVARS too small")
          endif
          vector_uvars(null_index)=name
          vector_vvars(null_index)=name
          vector_vvars(null_index)(1:1)='V'
          var_is_vector_uvar=null_index
       endif
    endif
#endif    
  end function var_is_vector_uvar


  function var_is_vector_vvar(name) 
    character(len=*), intent(in) :: name
    integer :: i, var_is_vector_vvar, null_index
    
    var_is_vector_vvar=0
    null_index=0
    do i=1,MAX_VECVARS
       if(trim(vector_vvars(i)).eq. '') then
          null_index=i
          exit
       endif
       if(trim(vector_vvars(i)).eq.name) then
          var_is_vector_vvar=i
          exit
       end if
    end do
#if 0
DISABLED: breaks in many cases, like UV
    if (var_is_vector_vvar==0) then
       ! default rules: if variable starts with V and was not found, add it to the list:
       if (name(1:1).eq.'V') then 
          if (null_index==0) then
             call abortmp("Error: MAX_VECVARS too small")
          endif
          vector_uvars(null_index)=name
          vector_uvars(null_index)(1:1)='U'
          vector_vvars(null_index)=name
          var_is_vector_vvar=null_index
       endif
    endif
#endif
    

  end function var_is_vector_vvar



end module interpolate_mod







