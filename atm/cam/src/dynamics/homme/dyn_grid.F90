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
  !      get_block_bounds_d    get first and last indices in global 
  !                            block ordering
  !      get_block_gcol_d      get column indices for given block
  !      get_block_gcol_cnt_d  get number of columns in given block
  !      get_block_lvl_cnt_d   get number of vertical levels in column
  !      get_block_levels_d    get vertical levels in column
  !      get_gcol_block_d      get global block indices and local columns 
  !                            index for given global column index
  !      get_gcol_block_cnt_d  get number of blocks containing data
  !                            from a given global column index
  !      get_block_owner_d     get process "owning" given block
  !      get_horiz_grid_d      get horizontal grid coordinates
  !      get_horiz_grid_dim_d  get horizontal dimensions of dynamics grid
  !
  ! Author: Jim Edwards and Patrick Worley
  ! 
  !-----------------------------------------------------------------------
  use element_mod, only : element_t
  use cam_logfile, only : iulog
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private

  integer           :: ngcols_d = 0     ! number of dynamics columns
  integer, public, parameter :: ptimelevels=2
  logical :: gblocks_need_initialized=.true.


  type block_global_data
     integer :: UniquePtOffset
     integer :: NumUniqueP
     integer :: LocalID
     integer :: Owner
  end type block_global_data


  type(block_global_data), allocatable :: gblocks(:)
  type(element_t), public, pointer :: elem(:)
  public :: get_block_owner_d, get_gcol_block_d, get_gcol_block_cnt_d, &
       get_block_gcol_cnt_d, get_horiz_grid_dim_d, get_block_levels_d, &
       get_block_gcol_d, get_block_bounds_d, get_horiz_grid_d, &
       get_block_lvl_cnt_d, set_horiz_grid_cnt_d, get_dyn_grid_parm, &
       get_dyn_grid_parm_real2d, get_dyn_grid_parm_real1d, get_block_ldof_d
  integer, pointer :: fdofP_local(:,:,:)

contains
  !========================================================================
  !
  Subroutine get_block_bounds_d(block_first,block_last)
    use dimensions_mod, only: nelem
    !----------------------------------------------------------------------- 
    ! 
    !                          
    ! Purpose: Return first and last indices used in global block ordering
    ! 
    ! Method: 
    ! 
    ! Author: Jim Edwards
    ! 
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    integer, intent(out) :: block_first  ! first (global) index used for blocks
    integer, intent(out) :: block_last   ! last (global) index used for blocks

    !-----------------------------------------------------------------------
    block_first = 1
    block_last = nelem

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
    ! Author: Jim Edwards
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !------------------------------Arguments--------------------------------
    integer, intent(in) :: blockid      ! global block id
    integer, intent(in) :: size         ! array size

    integer, intent(out):: cdex(size)   ! global column indices
    !
    integer :: ic
    if(gblocks_need_initialized) call gblocks_init()
    do ic=1,size
       cdex(ic)=gblocks(blockid)%UniquePtOffset+ic-1
    end do

    return
  end subroutine get_block_gcol_d

  subroutine get_block_ldof_d(nlev, ldof)
    use dof_mod, only : CreateMetaData
    use dimensions_mod, only : npsq, np, nelemd
    use parallel_mod,   only : par
    use element_mod, only : index_t

    integer, intent(in) :: nlev
    integer, pointer, intent(out) :: ldof(:)
    integer :: i, j, k, ie, ii
    type (index_t), pointer  :: idx     
    integer :: mydof(npsq)

    logical, save :: firstcall = .true.


    if(firstcall) then
       allocate(fdofp_local(np,np,nelemd))
       fdofp_local=0
       do ie=1,nelemd
          idx => elem(ie)%idxP
          do ii=1,idx%NumUniquePts
             i=idx%ia(ii)
             j=idx%ja(ii)          
             fdofp_local(i,j,ie) = (idx%UniquePtoffset+ii-1)
          end do
       end do
       firstcall=.false.
!    print *,__FILE__,__LINE__,maxval(fdofp_local),count(fdofp_local/=0)
    end if

    allocate(ldof(nelemd*npsq*nlev))       
    ldof = 0
    do ie=1,nelemd
       do j=1,np
          do i=1,np
             if(fdofp_local(i,j,ie)>0) then
                do k=0,nlev-1
                   ldof(i+(j-1)*np+k*npsq+(ie-1)*npsq*nlev ) = fdofP_local(i,j,ie) + ngcols_d*k
                end do
             end if
          end do
       end do
    end do
!    print *,__FILE__,__LINE__,maxval(ldof),count(ldof/=0)

  end subroutine get_block_ldof_d




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
    ! Author: Jim Edwards
    ! 
    !-----------------------------------------------------------------------
    integer, intent(in) :: blockid
    integer :: ie

    if(gblocks_need_initialized) call gblocks_init()
    get_block_gcol_cnt_d=gblocks(blockid)%NumUniqueP

    return
  end function get_block_gcol_cnt_d

  !
  !========================================================================
  !
  integer function get_block_lvl_cnt_d(blockid,bcid)
    use pmgrid, only: plevp

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

    get_block_lvl_cnt_d = plevp

    return
  end function get_block_lvl_cnt_d
  !
  !========================================================================
  !
  subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)
    use pmgrid, only: plev
    use abortutils, only: endrun

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
    use dimensions_mod, only : nelemd, nelem
    use kinds, only: int_kind
    use abortutils, only: endrun
    !----------------------------------------------------------------------- 
    ! 
    !                          
    ! Purpose: Return global block index and local column index
    !          for global column index
    !           element array is naturally ordered
    ! 
    ! Method: 
    ! 
    ! Author: Jim Edwards
    ! Modified: replaced linear search with best guess and binary search   
    !           (Pat Worley, 7/2/09)
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !------------------------------Arguments--------------------------------
    integer, intent(in) :: gcol     ! global column index
    integer, intent(in) :: cnt      ! size of blockid and bcid arrays

    integer, intent(out) :: blockid(cnt) ! block index
    integer, intent(out) :: bcid(cnt)    ! column index within block
    integer, intent(out), optional :: localblockid(cnt)
    integer :: sb, eb, ie, high, low
    logical :: found
    integer, save :: iedex_save = 1
    !---------------------------Local workspace-----------------------------
    !
    if(gblocks_need_initialized) call gblocks_init()
    found = .false.
    low = 1
    high = nelem

    ! check whether previous found element is the same here
    if (.not. found) then
       ie = iedex_save
       sb = gblocks(ie)%UniquePtOffset
       if (gcol >= sb) then
          eb = sb+gblocks(ie)%NumUniqueP
          if (gcol < eb) then
             found = .true.
          else
             low = ie
          endif
       else
          high = ie
       endif
    endif

    ! check whether next element  is the one wanted
    if ((.not. found) .and. &
        ((low .eq. iedex_save) .or. (iedex_save .eq. nelem))) then
       ie = iedex_save + 1
       if (ie > nelem) ie = 1

       sb = gblocks(ie)%UniquePtOffset
       if (gcol >= sb) then
          eb = sb+gblocks(ie)%NumUniqueP
          if (gcol < eb) then
             found = .true.
          else
             low = ie
          endif
       else
          high = ie
       endif
    endif

    ! otherwise, use a binary search to find element
    if (.not. found) then
       ! (start with a sanity check)
       ie = low
       sb = gblocks(ie)%UniquePtOffset

       ie = high
       eb = gblocks(ie)%UniquePtOffset + gblocks(ie)%NumUniqueP

       if ((gcol < sb) .or.  (gcol >= eb)) then
          do ie=1,nelemd
             write(iulog,*) __LINE__,ie,elem(ie)%idxP%UniquePtOffset,elem(ie)%idxP%NumUniquePts
          end do
          call endrun()
       end if

       do while (.not. found)
          ie = low + (high-low)/2;

          sb = gblocks(ie)%UniquePtOffset
          if (gcol >= sb) then
             eb = sb+gblocks(ie)%NumUniqueP
             if (gcol < eb) then
                found = .true.
             else
                low = ie+1
             endif
          else
             high = ie-1
          endif
       enddo
    endif

    blockid(1)=ie
    if(present(localblockid)) localblockid(1)=gblocks(ie)%LocalID
    bcid(1)=gcol-sb+1
    iedex_save = ie

 return
end subroutine get_gcol_block_d

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
 get_gcol_block_cnt_d = 1

 return
end function get_gcol_block_cnt_d

!
!========================================================================
!
integer function get_block_owner_d(blockid)
 use abortutils, only: endrun
 !----------------------------------------------------------------------- 
 ! 
 !                          
 ! Purpose: Return id of processor that "owns" the indicated block
 ! 
 ! Method: 
 ! 
 ! Author: Jim Edwards
 ! 
 !-----------------------------------------------------------------------

 implicit none
 !------------------------------Arguments--------------------------------
 integer, intent(in) :: blockid  ! global block id

 !-----------------------------------------------------------------------
 if(gblocks_need_initialized) call gblocks_init()

 if(gblocks(blockid)%Owner>-1) then
    get_block_owner_d=gblocks(blockid)%Owner
 else
    call endrun('Block owner not assigned in gblocks_init')
 end if

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
  !          For non-lon/lat grids, declare grid to be one-dimensional,
  !          i.e., (ngcols_d x 1)
  ! 
  ! Method: 
  ! 
  ! Author: Patrick Worley
  ! 
  !-----------------------------------------------------------------------

  !------------------------------Arguments--------------------------------
 integer, intent(out) :: hdim1_d           ! first horizontal dimension
 integer, intent(out), optional :: hdim2_d           ! second horizontal dimension
 !-----------------------------------------------------------------------
 hdim1_d = ngcols_d
 if(present(hdim2_d)) hdim2_d = 1

 return
end subroutine get_horiz_grid_dim_d
!
!========================================================================
!
subroutine set_horiz_grid_cnt_d(NumUniqueCols)
 integer, intent(in) :: NumUniqueCols
 !----------------------------------------------------------------------- 
 ! 
 !                          
 ! Purpose: Set number of columns in the dynamics computational grid
 ! 
 ! Method: 
 ! 
 ! Author: Jim Edwards
 ! 
 !-----------------------------------------------------------------------
 ngcols_d = NumUniqueCols


 return
end subroutine set_horiz_grid_cnt_d
!
!========================================================================
!
subroutine get_horiz_grid_d(nxy,clat_d_out,clon_d_out,area_d_out, &
    wght_d_out)

 use shr_kind_mod, only: r8 => shr_kind_r8
 use physconst, only : pi
 use dof_mod, only : UniqueCoords, UniquePoints
 use dimensions_mod, only: nelem, nelemd, nelemdmax, np
 use commap, only : clon, clat
#ifdef NO_MPI2
 use spmd_utils, only : iam, npes, mpi_integer, mpir8, mpicom
#else
 use spmd_utils, only : iam, npes, mpi_integer, mpir8, mpicom, mpi_in_place, mpi_datatype_null
#endif

 use abortutils, only: endrun
 !----------------------------------------------------------------------- 
 ! 
 !                          
 ! Purpose: Return latitude and longitude (in radians), column surface
 !          area (in radians squared) and surface integration weights
 !          for global column indices that will be passed to/from physics
 ! 
 ! Method: 
 ! 
 ! Author: Jim Edwards
 ! 
 !-----------------------------------------------------------------------
 integer, intent(in)   :: nxy                     ! array sizes

 real(r8), intent(out), optional :: clat_d_out(:) ! column latitudes
 real(r8), intent(out), optional :: clon_d_out(:) ! column longitudes
 real(r8), intent(out), target, optional :: area_d_out(:) ! column surface 
 !  area
 real(r8), intent(out), target, optional :: wght_d_out(:) ! column integration
 !  weight

 integer :: ie, i, sb, eb, p, j

 real (r8), dimension(np,np) :: areaw
 real(r8), pointer :: area_d(:)
 real(r8) :: scale_factor
 integer :: ierr
 integer :: rdispls(npes), recvcounts(npes), gid(npes)
#ifdef NO_MPI2
 integer :: ibuf
 real(r8) :: rbuf(nxy)
#endif
 logical :: compute_area 

 if (nxy < ngcols_d) then
    write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
         nxy,' < ',ngcols_d,' ) '
    call endrun
 endif

 compute_area = .false.
 if ( present(area_d_out) ) then
    if (size(area_d_out) .ne. nxy) then
       call endrun('bad area_d_out array size in dyn_grid')
    end if
    area_d => area_d_out
    compute_area = .true.
 else if ( present(wght_d_out) ) then
    if (size(wght_d_out) .ne. nxy) then
       call endrun('bad wght_d_out array size in dyn_grid')
    end if
    area_d => wght_d_out
    compute_area = .true.
 endif

 if(gblocks_need_initialized) call gblocks_init()
 ! could improve this
 if(.not. associated(clat)) then	
    allocate(clat(nxy), clon(nxy,1))
 end if

 clon(:,1) = -iam 
 do ie=1,nelemdmax        
    if(ie<=nelemd) then
       rdispls(iam+1)=elem(ie)%idxp%UniquePtOffset-1
       eb = rdispls(iam+1) + elem(ie)%idxp%NumUniquePts
       recvcounts(iam+1)=elem(ie)%idxP%NumUniquePts
       call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep, &
            clat(rdispls(iam+1)+1:eb), &
            clon(rdispls(iam+1)+1:eb,1)) 

       if (compute_area) then
          areaw = 1D0/elem(ie)%rspheremv(:,:)         
          call UniquePoints(elem(ie)%idxP, areaw, area_d(rdispls(iam+1)+1:eb))
       endif
    else     
       rdispls(iam+1) = 0
       recvcounts(iam+1) = 0
    endif
#ifdef NO_MPI2
    ibuf=rdispls(iam+1)
    call mpi_allgather(ibuf, 1, mpi_integer, rdispls, &
         1, mpi_integer, mpicom, ierr)

    ibuf=recvcounts(iam+1)
    call mpi_allgather(ibuf, 1, mpi_integer, recvcounts, &
         1, mpi_integer, mpicom, ierr)

    sb=rdispls(iam+1)+1
    eb =rdispls(iam+1) + recvcounts(iam+1)

    rbuf(1:recvcounts(iam+1)) = clat(sb:eb)  ! whats going to happen if end=0?
    call mpi_allgatherv(rbuf,recvcounts(iam+1),mpir8, clat, &
         recvcounts(:), rdispls(:), mpir8, mpicom, ierr)

    rbuf(1:recvcounts(iam+1)) = clon(sb:eb,1)
    call mpi_allgatherv(rbuf,recvcounts(iam+1),mpir8, clon, &
         recvcounts(:), rdispls(:), mpir8, mpicom, ierr)

    if (compute_area) then
       rbuf(1:recvcounts(iam+1)) = area_d(sb:eb)
       call mpi_allgatherv(rbuf,recvcounts(iam+1),mpir8, area_d, &
         recvcounts(:), rdispls(:), mpir8, mpicom, ierr)
    endif
#else
    call mpi_allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rdispls, &
         1, mpi_integer, mpicom, ierr)

    call mpi_allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvcounts, &
         1, mpi_integer, mpicom, ierr)

    call mpi_allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL, clat, &
         recvcounts(:), rdispls(:), mpir8, mpicom, ierr)

    call mpi_allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL, clon, &
         recvcounts(:), rdispls(:), mpir8, mpicom, ierr)
    if (compute_area) then
       call mpi_allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL, area_d, &
         recvcounts(:), rdispls(:), mpir8, mpicom, ierr)
    endif
#endif
 end do
 call mpi_barrier(mpicom,ierr)

 if ( present(clat_d_out) ) then
    if (size(clat_d_out) .ne. nxy) then
       call endrun('bad clat_d_out array size in dyn_grid')
    end if
    clat_d_out(:) = clat(:)
 endif

 if ( present(clon_d_out) ) then
    if (size(clon_d_out) .ne. nxy) then
       call endrun('bad clon_d_out array size in dyn_grid')
    end if
    clon_d_out(:) = clon(:,1)
 endif

 !if one of area_d_out  or wght_d_out was present, than it was computed
 !above.  if they were *both* present, then do this: 
 if ( present(area_d_out) .and. present(wght_d_out) ) then
    wght_d_out(:) = area_d_out(:)
 endif


 return
end subroutine get_horiz_grid_d


subroutine get_gcol_lat(gcid, lat)
  use shr_kind_mod, only: r8 => shr_kind_r8
  use commap, only : clat
  integer, intent(in) :: gcid(:)
  real(r8), intent(out) :: lat(:)
  integer :: k, glen

  glen = size(gcid)

  do k=1,glen
     lat(k) = clat(gcid(k))
  end do

end subroutine get_gcol_lat

subroutine get_gcol_lon(gcid, lon)
  use commap, only : clon
  integer, intent(in) :: gcid(:)
  real(r8), intent(out) :: lon(:)
  integer :: k, glen

  glen = size(gcid)

  do k=1,glen
     lon(k) = clon(gcid(k),1)
  end do

end subroutine get_gcol_lon



subroutine gblocks_init()
 use dimensions_mod, only : nelem, nelemd, nelemdmax
 use spmd_utils, only : iam, npes, mpi_integer, mpicom
 use abortutils, only: endrun
#ifndef NO_MPI2
 use spmd_utils, only : mpi_in_place, mpi_datatype_null
#endif
 integer :: ie, p
 integer :: ibuf
 integer :: ierr
 integer :: rdispls(npes), recvcounts(npes), gid(npes), lid(npes)


 if(.not.allocated(gblocks)) then
    allocate(gblocks(nelem))
    do ie=1,nelem
       gblocks(ie)%Owner=-1
       gblocks(ie)%UniquePtOffset=-1
       gblocks(ie)%NumUniqueP    =-1
       gblocks(ie)%LocalID       =-1
    end do
 end if

 do ie=1,nelemdmax        
    if(ie<=nelemd) then
       rdispls(iam+1)=elem(ie)%idxP%UniquePtOffset-1
       gid(iam+1)=elem(ie)%GlobalID
       lid(iam+1)=ie
       recvcounts(iam+1)=elem(ie)%idxP%NumUniquePts
    else     
       rdispls(iam+1) = 0
       recvcounts(iam+1) = 0
       gid(iam+1)=0
    endif
#ifdef NO_MPI2
    ibuf=lid(iam+1)
    call mpi_allgather(ibuf, 1, mpi_integer, lid, &
         1, mpi_integer, mpicom, ierr)

    ibuf=gid(iam+1)
    call mpi_allgather(ibuf,1,mpi_integer, gid, &
         1, mpi_integer, mpicom, ierr)

    ibuf=rdispls(iam+1)
    call mpi_allgather(ibuf, 1, mpi_integer, rdispls, &
         1, mpi_integer, mpicom, ierr)

    ibuf=recvcounts(iam+1)
    call mpi_allgather(ibuf, 1, mpi_integer, recvcounts, &
         1, mpi_integer, mpicom, ierr)
#else
    call mpi_allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, lid, &
         1, mpi_integer, mpicom, ierr)

    call mpi_allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, gid, &
         1, mpi_integer, mpicom, ierr)

    call mpi_allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rdispls, &
         1, mpi_integer, mpicom, ierr)

    call mpi_allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, recvcounts, &
         1, mpi_integer, mpicom, ierr)
#endif
    do p=1,npes
       if(gid(p)>0) then
          gblocks(gid(p))%UniquePtOffset=rdispls(p)+1
          gblocks(gid(p))%NumUniqueP=recvcounts(p)
          gblocks(gid(p))%LocalID=lid(p)
          gblocks(gid(p))%Owner=p-1
       end if
    end do
 end do
 gblocks_need_initialized=.false.

end subroutine gblocks_init



!#######################################################################
   function get_dyn_grid_parm_real2d(name) result(rval)
     use commap, only : londeg, clon
     character(len=*), intent(in) :: name
     real(r8), pointer :: rval(:,:)

     if(name.eq.'clon') then
        rval => clon
     else if(name.eq.'londeg') then
        rval => londeg
     else
        nullify(rval)
     end if
   end function get_dyn_grid_parm_real2d

!#######################################################################
   function get_dyn_grid_parm_real1d(name) result(rval)
     use commap, only : latdeg, clat, w
     character(len=*), intent(in) :: name
     real(r8), pointer :: rval(:)

     if(name.eq.'clat') then
        rval => clat
     else if(name.eq.'latdeg') then
        rval => latdeg
     else if(name.eq.'w') then
        rval => w
     else
        nullify(rval)
     end if
   end function get_dyn_grid_parm_real1d




integer function get_dyn_grid_parm(name) result(ival)
 use pmgrid, only : beglat, endlat, plat, plon, plev, plevp
 use dimensions_mod, only: ne, np, nelemd, npsq
 use interpolate_mod, only : get_interp_parameter
 character(len=*), intent(in) :: name


 if(name.eq.'ne') then
    ival = ne
 else if(name.eq.'np') then
    ival = np
 else if(name.eq.'npsq') then
    ival = npsq
 else if(name.eq.'nelemd') then
    ival = nelemd
 else if(name.eq.'beglat') then
    ival = beglat
 else if(name.eq.'endlat') then
    ival = endlat
 else if(name.eq.'plat') then
    ival = plat
 else if(name.eq.'plon') then
    ival = ngcols_d
 else if(name.eq.'plev') then
    ival = plev
 else if(name.eq.'plevp') then
    ival = plevp
 else if(name.eq.'nlon') then
    ival = get_interp_parameter('nlon')
 else if(name.eq.'nlat') then
    ival = get_interp_parameter('nlat')
 else	
    ival = -1
 end if
end function get_dyn_grid_parm

!#######################################################################
end module dyn_grid
