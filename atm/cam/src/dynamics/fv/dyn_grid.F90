module dyn_grid
!----------------------------------------------------------------------- 
! 
! Purpose: Definition of dynamics computational grid.
!
! Method: Variables are private; interface routines used to extract
!         information for use in user code. Global column index range
!         defined using full (unreduced) grid. 
! 
! Entry points:
!      get_block_bounds_d  get first and last indices in global 
!                          block ordering
!      get_block_gcol_d      get column indices for given block
!      get_block_gcol_cnt_d  get number of columns in given block
!      get_block_lvl_cnt_d get number of vertical levels in column
!      get_block_levels_d  get vertical levels in column
!      get_gcol_block_d      get global block indices and local columns 
!                            index for given global column index
!      get_gcol_block_cnt_d  get number of blocks containing data
!                            from a given global column index
!      get_block_owner_d   get process "owning" given block
!      get_horiz_grid_d      get horizontal grid coordinates
!      get_horiz_grid_dim_d  get horizontal dimensions of dynamics grid
!
! Author: John Drake and Patrick Worley
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plev
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog

   integer, private   :: ngcols_d = 0     ! number of dynamics columns
! WS 2006.04.12:  moved ptimelevels here from prognostics, which is gone
   integer, parameter :: ptimelevels = 2  ! number of time levels in the dycore

   real(r8), parameter ::  D0_5                    =   0.5_r8
   real(r8), parameter ::  D90_0                   =  90.0_r8
   real(r8), parameter ::  D360_0                  = 360.0_r8
contains
  subroutine get_block_ldof_d( nlev, ldof)
    use pmgrid, only : beglonxy, endlonxy, beglatxy, endlatxy
    integer, intent(in) :: nlev
    integer, intent(out) :: ldof(:)

    integer:: i,j,k, lcnt
#if 0
    lcnt=(endlatxy-beglatxy+1)*nlev*(endlonxy-beglonxy+1)
    allocate(ldof(lcnt))
    lcnt=0
    ldof(:)=0	
    do k=1,nlev
       do j=beglatxy,endlatxy
          do i=beglonxy, endlonxy
             lcnt=lcnt+1
             if(j.eq.1.and.hdim2_d==plat-1) then
                ldof(lcnt)=0
             else if(hdim2_d>1) then
                if(fileorder.eq.'xyz') then
                   ldof(lcnt)=i+(j-(plat-hdim2_d+1))*hdim1_d+(k-1)*hdim1_d*hdim2_d
                else 
                   ldof(lcnt)=i+(j-(plat-hdim2_d+1))*hdim1_d*nlev+(k-1)*hdim1_d
                end if
             else  ! lon, lev decomp used for history nacs
                ldof(lcnt)=i+(k-1)*hdim1_d
             end if
          end do
       end do
    end do
#endif


  end subroutine get_block_ldof_d




!========================================================================
!
   Subroutine get_block_bounds_d(block_first,block_last)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return first and last indices used in global block ordering
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: plat, spmd_on, nprxy_x, nprxy_y

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(out) :: block_first  ! first (global) index used for blocks
   integer, intent(out) :: block_last   ! last (global) index used for blocks

!-----------------------------------------------------------------------
!  latitude slice block
   block_first = 1
   if (spmd_on .eq. 1) then
! Assume 1 block per subdomain
      block_last  = nprxy_x*nprxy_y
   else
      block_last  = plat
   endif

   return
   end subroutine get_block_bounds_d
!
!========================================================================
!
   subroutine get_block_gcol_d(blockid,size,cdex)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return list of dynamics column indices in given block
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use rgrid,    only: nlon
   use pmgrid,   only: spmd_on, plon, nprxy_x
   use spmd_utils, only: iam
#if ( defined SPMD )
   use spmd_dyn, only: lonrangexy, latrangexy, npes_xy
#endif

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid      ! global block id
   integer, intent(in) :: size         ! array size

   integer, intent(out):: cdex(size)   ! global column indices
!---------------------------Local workspace-----------------------------
!
    integer i,j                        ! block coordinates
    integer blksiz                     ! block size
    integer k,l                        ! loop indices
    integer n                          ! column index
!-----------------------------------------------------------------------

   if (spmd_on .eq. 1) then
      j = (blockid-1) / nprxy_x + 1
      i = blockid - (j-1) * nprxy_x
#if ( defined SPMD )
      blksiz = (lonrangexy(2,i)-lonrangexy(1,i)+1) *       &
               (latrangexy(2,j)-latrangexy(1,j)+1)
      if (size < blksiz) then
         write(iulog,*)'GET_BLOCK_GCOL_D: array not large enough (', &
                             size,' < ',blksiz,' ) '
         call endrun
      else
         n = 0
         do k=latrangexy(1,j),latrangexy(2,j)
            do l=lonrangexy(1,i),lonrangexy(2,i)
               n = n + 1
               cdex(n) = l + (k-1)*plon
            enddo
         enddo
      endif
#endif
   else
      if (size < nlon(blockid)) then
         write(iulog,*)'GET_BLOCK_GCOL_D: array not large enough (', &
                             size,' < ',nlon(blockid),' ) '
         call endrun
      else
         n = (blockid-1)*plon
         do i = 1,nlon(blockid)
            n = n + 1
            cdex(i) = n
         enddo
      endif
   endif
!
   return
   end subroutine get_block_gcol_d
!
!========================================================================
!
   integer function get_block_gcol_cnt_d(blockid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return number of dynamics columns in indicated block
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use rgrid, only: nlon
   use pmgrid, only: spmd_on, nprxy_x
#if ( defined SPMD )
   use spmd_dyn, only: lonrangexy, latrangexy
#endif

   implicit none

!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id

!---------------------------Local workspace-----------------------------
   integer i, j

!-----------------------------------------------------------------------
   if (spmd_on .eq. 1) then
      j = (blockid-1) / nprxy_x + 1
      i = blockid - (j-1) * nprxy_x
#if ( defined SPMD )
      get_block_gcol_cnt_d = (lonrangexy(2,i)-lonrangexy(1,i)+1) *       &
         (latrangexy(2,j)-latrangexy(1,j)+1)
#endif
   else
      get_block_gcol_cnt_d = nlon(blockid)
   endif

   return
   end function get_block_gcol_cnt_d

!
!========================================================================
!
   integer function get_block_lvl_cnt_d(blockid,bcid)

!-----------------------------------------------------------------------
!
!
! Purpose: Return number of levels in indicated column. If column
!          includes surface fields, then it is defined to also
!          include level 0.
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid    ! column index within block

!-----------------------------------------------------------------------
!  latitude slice block
   get_block_lvl_cnt_d = plev + 1

   return
   end function get_block_lvl_cnt_d
!
!========================================================================
!
   subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

!-----------------------------------------------------------------------
!
!
! Purpose: Return level indices in indicated column. If column
!          includes surface fields, then it is defined to also
!          include level 0.
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid    ! column index within block
   integer, intent(in) :: lvlsiz   ! dimension of levels array

   integer, intent(out) :: levels(lvlsiz) ! levels indices for block

!---------------------------Local workspace-----------------------------
!
    integer k                      ! loop index
!-----------------------------------------------------------------------
!  latitude slice block
   if (lvlsiz < plev + 1) then
      write(iulog,*)'GET_BLOCK_LEVELS_D: levels array not large enough (', &
                          lvlsiz,' < ',plev + 1,' ) '
      call endrun
   else
      do k=0,plev
         levels(k+1) = k
      enddo
      do k=plev+2,lvlsiz
         levels(k) = -1
      enddo
   endif

   return
   end subroutine get_block_levels_d

!
!========================================================================
!
   subroutine get_gcol_block_d(gcol,cnt,blockid,bcid,localblockid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return global block index and local column index
!          for global column index
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: plon, spmd_on, nprxy_x, nprxy_y
#if ( defined SPMD )
   use spmd_dyn, only: lonrangexy, latrangexy
#endif

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: gcol     ! global column index
   integer, intent(in) :: cnt      ! size of blockid and bcid arrays

   integer, intent(out) :: blockid(cnt) ! block index
   integer, intent(out) :: bcid(cnt)    ! column index within block
   integer, intent(out), optional :: localblockid(cnt)
!---------------------------Local workspace-----------------------------
!
    integer i,j,ii,jj                   ! loop indices
    integer glon, glat                  ! global longitude and latitude
                                        ! indices
    integer ddlon                       ! number of longitudes in block
!-----------------------------------------------------------------------

!  lon/lat block
   if (cnt < 1) then
      write(iulog,*)'GET_GCOL_BLOCK_D: arrays not large enough (', cnt,' < ',1,' ) '
      call endrun
   else
      if (spmd_on .eq. 1) then
! Determine global latitude and longitude coordinate indices from
! global column index
         glat = (gcol-1)/plon + 1
         glon = gcol - ((glat-1)*plon)

! Determine block coordinates (ii,jj), where ii ranges from 1 to
! nprxy_x and jj ranges from 1 to nprxy_y.
#if ( defined SPMD )
         ii=0
         do i=1,nprxy_x
            if (lonrangexy(1,i) .le. glon .and. glon .le. lonrangexy(2,i)) ii=i
         enddo
         jj=0
         do j=1,nprxy_y
            if (latrangexy(1,j) .le. glat .and. glat .le. latrangexy(2,j)) jj=j
         enddo
         if (ii .eq. 0 .or. jj .eq. 0) then
            write(iulog,*)'GET_GCOL_BLOCK_D: could not find block indices for (', &
                      glon,',',glat,' ) '
            call endrun
         endif

! Global block index
         blockid(1) = (jj-1)*nprxy_x+ii

! Local coordinates in block
         j = glat-latrangexy(1,jj)+1
         i = glon-lonrangexy(1,ii)+1
         ddlon = lonrangexy(2,ii)-lonrangexy(1,ii)+1

! Local column index in block
         bcid(1) = (j-1)*ddlon+i
!
#endif
      else
         glat = (gcol-1)/plon + 1
         glon = gcol - ((glat-1)*plon)

         blockid(1) = glat
         bcid(1)    = glon
      endif
!
      do j=2,cnt
         blockid(j) = -1
         bcid(j)    = -1
      enddo
!
   endif
!

   return
   end subroutine get_gcol_block_d

   subroutine get_gcol_lat(gcid, lat)
     use pmgrid,     only: plat, plon
     use commap,     only: clat
     integer, intent(in) :: gcid(:)
     real(r8), intent(out) :: lat(:)
     integer :: k, glen, ilat 
     
     glen = size(gcid)

     do k=1,glen
        ilat = (gcid(k)-1)/plon + 1 
        lat(k) = CLAT(ilat)
     end do
     
   end subroutine get_gcol_lat

   subroutine get_gcol_lon(gcid, lon)
     use pmgrid,     only: plat, plon
     use commap,     only: clon
     integer, intent(in) :: gcid(:)
     real(r8), intent(out) :: lon(:)
     
     integer :: k, glen, ilat
     
     glen = size(gcid)
        

     do k=1,glen
        ilat = (gcid(k)-1)/plon + 1
        lon(k) = CLON(gcid(k) - ((ilat-1)*plon),1)
     end do
     
   end subroutine get_gcol_lon







!
!========================================================================
!
   integer function get_gcol_block_cnt_d(gcol)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return number of blocks contain data for the vertical column
!          with the given global column index
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: gcol     ! global column index
!-----------------------------------------------------------------------
!  lon/lat block
   get_gcol_block_cnt_d = 1

   return
   end function get_gcol_block_cnt_d

!
!========================================================================
!
   integer function get_block_owner_d(blockid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return id of processor that "owns" the indicated block
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: spmd_on
#if ( defined SPMD )
   use spmd_dyn, only: proc
#endif

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id

!-----------------------------------------------------------------------
!  latitude slice block
#if (defined SPMD)
   if (spmd_on .eq. 1) then
      get_block_owner_d = blockid - 1
   else
      get_block_owner_d = proc(blockid)
   endif
#else
   get_block_owner_d = 0
#endif

   return
   end function get_block_owner_d

!
!========================================================================
!
   subroutine get_horiz_grid_dim_d(hdim1_d,hdim2_d)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Returns declared horizontal dimensions of computational grid.
!          Note that global column ordering is assumed to be compatible
!          with the first dimension major ordering of the 2D array.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid,     only: plat, plon

!------------------------------Arguments--------------------------------
   integer, intent(out) :: hdim1_d           ! first horizontal dimension
   integer, intent(out) :: hdim2_d           ! second horizontal dimension
!-----------------------------------------------------------------------
   if (ngcols_d == 0) then
      ngcols_d = plat*plon
   endif
   hdim1_d = plon
   hdim2_d = plat

   return
   end subroutine get_horiz_grid_dim_d
!
!========================================================================
!
   subroutine get_horiz_grid_d(size,clat_d_out,clon_d_out, &
        area_d_out, wght_d_out)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return latitude and longitude (in radians), column surface
!          area (in radians squared) and surface integration weights
!          for global column indices that will be passed to/from physics
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid,     only: plat, plon
   use rgrid,      only: nlon
   use commap,     only: clat, clon, latdeg, londeg, w
   use physconst,  only: pi, spval

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in)   :: size             ! array sizes

   real(r8), intent(out), optional :: clat_d_out(size) ! column latitudes
   real(r8), intent(out), optional :: clon_d_out(size) ! column longitudes

   real(r8), intent(out), optional :: area_d_out(size) ! column surface 
                                                       !  area
   real(r8), intent(out), optional :: wght_d_out(size) ! column integration
                                                       !  weight
!---------------------------Local workspace-----------------------------
!
    integer i,j                 ! loop indices
    integer n                   ! column index
    real(r8) :: ns_vert(2,plon) ! latitude grid vertices
    real(r8) :: ew_vert(2,plon) ! longitude grid vertices
    real(r8) :: del_theta       ! difference in latitude at a grid point
    real(r8) :: del_phi         ! difference in longitude at a grid point
    real(r8),parameter :: degtorad=pi/180.0_r8 ! convert degrees to radians
    
!-----------------------------------------------------------------------
!   if (size < ngcols_d) then
!      write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
!                          size,' < ',ngcols_d,' ) '
!      call endrun
!   else
    if(present(clon_d_out)) then
       if(size == ngcols_d) then
          n = 0
          do j = 1,plat
             do i = 1,nlon(j)
                n = n + 1
                clon_d_out(n) = clon(i,j)
             enddo
          enddo
       else if(size == plon) then
          clon_d_out(:) = clon(:,1)
       else
          write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
               size,' < ',ngcols_d,' ) '
          call endrun
       end if
    end if
    if(present(clat_d_out)) then
       if(size == ngcols_d) then
          n = 0
          do j = 1,plat
             do i = 1,nlon(j)
                n = n + 1
                clat_d_out(n) = clat(j)
             enddo
          enddo
       else if(size == plat) then
          clat_d_out(:) = clat(:)
       else
          write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
               size,' < ',ngcols_d,' ) '
          call endrun
       end if
    end if
    if(size==plat .and. present(wght_d_out)) then
       wght_d_out(:) = w(:)
       return
    endif
      
    if ( ( present(area_d_out) ) .or. ( present(wght_d_out) ) ) then
       if(size==plat .and.present(area_d_out)) then
          call endrun('size argument to get_horiz_grid_d too small')
       end if
       n = 0
       do j = 1,plat

          ! First, determine vertices of each grid point. 
          ! Verticies are ordered as follows: 
          ! ns_vert: 1=lower left, 2 = upper left
          ! ew_vert: 1=lower left, 2 = lower right

          ! Latitude vertices
          ns_vert(:,:) = spval
          if (j .eq. 1) then
             ns_vert(1,:nlon(1))    = -D90_0 + (latdeg(1) - latdeg(2))*D0_5
          else
             ns_vert(1,:nlon(j))    = (latdeg(j) + latdeg(j-1) )*D0_5
          endif
          
          if (j .eq. plat) then
             ns_vert(2,:nlon(plat)) =  D90_0 + (latdeg(plat) - latdeg(plat-1))*D0_5
          else
             ns_vert(2,:nlon(j))    = (latdeg(j) + latdeg(j+1) )*D0_5
          endif

          ! Longitude vertices
          ew_vert(:,:) = spval
          ew_vert(1,1)          = (londeg(1,j) - D360_0 + londeg(nlon(j),j))*D0_5
          ew_vert(1,2:nlon(j))  = (londeg(1:nlon(j)-1,j)+ londeg(2:nlon(j),j))*D0_5
          ew_vert(2,:nlon(j)-1) = ew_vert(1,2:nlon(j))
          ew_vert(2,nlon(j))    = (londeg(nlon(j),j) + (D360_0 + londeg(1,j)))*D0_5

          do i = 1,nlon(j)
             n = n + 1
             
             if (j .eq. 1) then
                del_phi = -sin( latdeg(j)*degtorad )    + sin( ns_vert(2,i)*degtorad )
             else if (j .eq. plat) then
                del_phi =  sin( latdeg(j)*degtorad )    - sin( ns_vert(1,i)*degtorad )
             else
                del_phi =  sin( ns_vert(2,i)*degtorad ) - sin( ns_vert(1,i)*degtorad )
             end if
             
             del_theta = ( ew_vert(2,i) - ew_vert(1,i) )*degtorad
             
             if ( present(area_d_out) ) area_d_out(n) = del_theta*del_phi
             if (present(wght_d_out) ) wght_d_out(n) = del_theta*del_phi
          enddo

       enddo

    endif

   

   return
   end subroutine get_horiz_grid_d


!#######################################################################
   function get_dyn_grid_parm_real2d(name) result(rval)
     use commap, only : londeg, londeg_st, clon
     character(len=*), intent(in) :: name
     real(r8), pointer :: rval(:,:)

     if(name.eq.'clon') then
        rval => clon
     else if(name.eq.'londeg') then
        rval => londeg
     else if(name.eq.'londeg_st') then
        rval => londeg_st
     else
        nullify(rval)
     end if
   end function get_dyn_grid_parm_real2d

!#######################################################################
   function get_dyn_grid_parm_real1d(name) result(rval)
     use commap, only : latdeg, latdeg_st, clat_staggered, w_staggered, clat, w
     character(len=*), intent(in) :: name
     real(r8), pointer :: rval(:)

     if(name.eq.'clat') then
        rval => clat
     else if(name.eq.'latdeg') then
        rval => latdeg
     else if(name.eq.'latdeg_st') then
        rval => latdeg_st
     else if(name.eq.'clatdeg_staggered') then
        rval => latdeg_st
     else if(name.eq.'w') then
        rval => w
     else if(name.eq.'w_staggered') then
        rval => w_staggered
     else
        nullify(rval)
     end if
   end function get_dyn_grid_parm_real1d



   integer function get_dyn_grid_parm(name) result(ival)
     use pmgrid, only : beglat, endlat, plat, plon, plev, plevp, &
          splon, beglev, endlev, endlevp, &
          beglonxy, endlonxy, beglatxy, endlatxy
     character(len=*), intent(in) :: name

     if(name.eq.'splon') then
        ival = splon
     else if(name.eq.'beglev') then
        ival = beglev
     else if(name.eq.'endlev') then
        ival = endlev
     else if(name.eq.'endlevp') then
        ival = endlevp
     else if(name.eq.'beglonxy') then
        ival = beglonxy
     else if(name.eq.'endlonxy') then
        ival = endlonxy
     else if(name.eq.'beglatxy') then
        ival = beglatxy
     else if(name.eq.'endlatxy') then
        ival = endlatxy
     else if(name.eq.'beglat') then
        ival = beglat
     else if(name.eq.'endlat') then
        ival = endlat
     else if(name.eq.'plat') then
        ival = plat
     else if(name.eq.'plon') then
        ival = plon
     else if(name.eq.'plev') then
        ival = plev
     else if(name.eq.'plevp') then
        ival = plevp
     else	
        ival = -1
     end if
     return
   end function get_dyn_grid_parm


end module dyn_grid
