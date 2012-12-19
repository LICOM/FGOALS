module phys_buffer

!----------------------------------------------------------------------- 
! 
! Purpose: 
!   Implement a physics buffer to hold arrays that must persist
!   across timesteps or between calls to different physics packages.
!
!   Current implementation only supports one buffer.
!
! Author: B. Eaton
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use infnan,       only: nan
   use spmd_utils,   only: masterproc
   use ppgrid,       only: pcols, begchunk, endchunk
   use phys_grid,    only: physgrid_set, write_field_from_chunk, read_chunk_from_field, &
                           get_ncols_p
   use dyn_grid,     only: ptimelevels
   use string_utils, only: to_upper
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog
#ifdef SPMD
   use mpishorthand, only: mpicom, mpiint
#endif
   use pio,          only: file_desc_t, var_desc_t, io_desc_t, &
                           pio_setdebuglevel, pio_seterrorhandling, &
                           pio_double, pio_int, pio_noerr, &
                           pio_inq_dimid, pio_inq_varid, &
                           pio_def_dim, pio_def_var, &
                           pio_put_var, pio_get_var, &
                           pio_write_darray, pio_read_darray, &
                           pio_freedecomp, &
                           pio_bcast_error, pio_internal_error

   implicit none
   private
   save

! Public methods 

   public ::&
      pbuf_defaultopts,   &! set namelist defaults
      pbuf_setopts,       &! set user specified namelist values
      pbuf_init,          &! initialize physics buffer
      pbuf_add,           &! add field to physics buffer
      pbuf_print,         &! print summary of fields in physics buffer
      pbuf_get_fld_idx,   &! get index of specified field in the physics buffer
      pbuf_get_fld_name,  &! get name of specified field given index
      pbuf_old_tim_idx,   &! return the index for the oldest time
      pbuf_next_tim_idx,  &! return the index for the next time
      pbuf_update_tim_idx,&! update the index for the oldest time
      pbuf_allocate,      &! allocate memory for physics buffer fields
      pbuf_deallocate,    &! deallocate memory for physics buffer fields
      pbuf_setval,        &! set value for a field in the buffer
      pbuf_write_restart, &! write physics buffer to restart file
      pbuf_read_restart    ! read physics buffer from restart file
   public :: pbuf_init_restart
   type(var_desc_t) :: timeidx_desc

! Public types and data

   type, public :: pbuf_fld
      character(len=16)                       :: name
      character(len=16)                       :: scope
      integer                                 :: fdim, mdim, ldim
      real(r8), pointer, dimension(:,:,:,:,:) :: fld_ptr
      type(var_desc_t)                        :: vardesc
      type(var_desc_t), pointer               :: vardesc_list(:)   ! used to seperate variables of > 4MB (required for netcdf)
      type(io_desc_t), pointer                :: iodesc
   end type pbuf_fld


   integer, public, parameter :: pbuf_size_max=1000
   type(pbuf_fld), public, dimension(pbuf_size_max) :: pbuf

   integer, public :: pbuf_times  ! number of time levels in physics buffer (dycore dependent)

! Private module data

   integer :: pbuf_size = 0
   integer :: old_time_idx = 1
   logical :: global_allocate_all = .true.  ! allocate all buffers as global
   integer, parameter :: maxvarsize = 536870910
!=========================================================================================
contains
!=========================================================================================

subroutine pbuf_defaultopts(pbuf_global_allocate_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   logical, intent(out), optional :: pbuf_global_allocate_out
!-----------------------------------------------------------------------
   if ( present(pbuf_global_allocate_out) ) then
      pbuf_global_allocate_out = global_allocate_all
   endif
end subroutine pbuf_defaultopts

!=========================================================================================

subroutine pbuf_setopts(pbuf_global_allocate_in)
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
!-----------------------------------------------------------------------

   logical, intent(in), optional :: pbuf_global_allocate_in
!-----------------------------------------------------------------------
   if ( present(pbuf_global_allocate_in) ) then
      global_allocate_all = pbuf_global_allocate_in
   endif
end subroutine pbuf_setopts

!=========================================================================================

subroutine pbuf_init()

! Initialize physics buffer.

   implicit none
!-----------------------------------------------------------------------------------------
   integer :: i
   pbuf_times = ptimelevels - 1

   do i=1,pbuf_size_max
      nullify(pbuf(i)%fld_ptr)
      nullify(pbuf(i)%iodesc)
   end do


end subroutine pbuf_init

!=========================================================================================

subroutine pbuf_add(name, scope, fdim, mdim, ldim, index)

! Add a field to the physics buffer

   implicit none

   character(len=*), intent(in)  :: name   ! field name (case insensitive)
   character(len=*), intent(in)  :: scope  ! field scope, either 'global' or 'physpkg'
   integer,          intent(in)  :: fdim   ! first generic field dimension
   integer,          intent(in)  :: mdim   ! middle generic field dimension
   integer,          intent(in)  :: ldim   ! last generic field dimension
   integer,          intent(out) :: index  ! index in the physics buffer

! Local variables
   character(len=*), parameter :: sub = 'pbuf_add'
   integer :: i
   character(len=len(name)) :: uname       ! =to_upper(name)
!-----------------------------------------------------------------------------------------

   if ( pbuf_size >= pbuf_size_max ) then
      call endrun (sub//': max number physics buffer fields exceeded. Increase pbuf_size_max in phys_buffer.F90')
   end if

   do i = 1, pbuf_size
      if ( pbuf(i)%name == name ) then
         call endrun (sub//': ERROR: field name '//name//' is already in use.')
      end if
   end do

   if ( scope /= 'global' .and. scope /= 'physpkg' ) then
      call endrun (sub//': scope must be set to global or physpkg.  Current value is: '//scope)
   end if

   pbuf_size = pbuf_size + 1
   index = pbuf_size
   pbuf(index)%name = name
   pbuf(index)%scope = scope
   pbuf(index)%fdim = fdim
   pbuf(index)%mdim = mdim
   pbuf(index)%ldim = ldim

end subroutine pbuf_add

!=========================================================================================

subroutine pbuf_print()

   ! Print summary of fields in physics buffer

   ! Local variables
   integer :: i
   !-----------------------------------------------------------------------------------------

   write(iulog,*) ' Fields in physics buffer'
   write(iulog,*) ' name, scope, fdim, mdim, ldim'
   do i = 1, pbuf_size
      write(iulog,*) pbuf(i)%name, pbuf(i)%scope, pbuf(i)%fdim, pbuf(i)%mdim, pbuf(i)%ldim
   end do

end subroutine pbuf_print

!=========================================================================================

function pbuf_get_fld_name(idx)
  character(len=16)  :: pbuf_get_fld_name
  integer, intent(in) :: idx
  pbuf_get_fld_name = pbuf(idx)%name
  return
end function pbuf_get_fld_name

!=========================================================================================
function pbuf_get_fld_idx(name, failcode)

! Get index of specified field in the physics buffer.  String matching is case insensitive.
! Call endrun if name not found

   implicit none

   character(len=*), intent(in)  :: name   ! field name 
   integer, intent(in), optional :: failcode

! Return value
   integer :: pbuf_get_fld_idx

! Local variables
   integer :: i
   character(len=len(name)) :: Uname       ! =to_upper(name)
!-----------------------------------------------------------------------------------------

!
!  Search for specified field in physics buffer, assuming that case of
!  argument "name" matches definition in pbuf structure.
!

   do i = 1, pbuf_size
      if ( pbuf(i)%name == name ) then
         pbuf_get_fld_idx = i
         return
      end if
   end do

   if ( present(failcode)) then
      pbuf_get_fld_idx = failcode	
      return
   else
      call endrun ('PBUF_GET_FLD_IDX: index not found for '//name)
   endif

end function pbuf_get_fld_idx

!=========================================================================================

function pbuf_old_tim_idx()

! Return index of oldest time sample in the physics buffer.

   implicit none

! Return value
   integer :: pbuf_old_tim_idx
!-----------------------------------------------------------------------------------------

   pbuf_old_tim_idx = old_time_idx

end function pbuf_old_tim_idx

!=========================================================================================

function pbuf_next_tim_idx(idx)

! Return index of next time sample in the physics buffer.

   implicit none

   integer, intent(in) :: idx

! Return value
   integer :: pbuf_next_tim_idx
!-----------------------------------------------------------------------------------------

   pbuf_next_tim_idx = mod(idx, pbuf_times) + 1

end function pbuf_next_tim_idx

!=========================================================================================

subroutine pbuf_update_tim_idx()

! Update index of old time sample in the physics buffer.

   implicit none
!-----------------------------------------------------------------------------------------

   old_time_idx = mod(old_time_idx, pbuf_times) + 1

end subroutine pbuf_update_tim_idx

!=========================================================================================

subroutine pbuf_allocate(scope)

! Allocate storage for fields in the physics buffer with the specified scope.
! If global_allocate_all=.true. then storage for both global and physpkg scope 
! is allocated just once, when scope='global'.
!
! N.B. This routine must be called after phys_grid_init because that's
!      where begchunk and endchunk are set

   implicit none

   character(len=*), intent(in)  :: scope

! Local variables
   character(len=*), parameter :: sub = 'pbuf_allocate'
   integer :: i, fdim, mdim, ldim, istat
   logical :: allocate_all
!-----------------------------------------------------------------------------------------

   if ( .not. physgrid_set ) then
      call endrun (sub//': ERROR: called before physics grid initialized.')
   end if

   ! allocate_all is used to force allocation of all fields at same time as allocation
   ! for global scope.
   allocate_all = .false.
   if ( global_allocate_all ) then
      if ( scope == 'global' ) then
         allocate_all = .true.
      else
         return
      end if
   end if

   do i = 1, pbuf_size
      if ( pbuf(i)%scope == scope  .or.  allocate_all ) then
         fdim = pbuf(i)%fdim
         mdim = pbuf(i)%mdim
         ldim = pbuf(i)%ldim

         allocate(pbuf(i)%fld_ptr(fdim,pcols,mdim,begchunk:endchunk,ldim), stat=istat)
         if ( istat /= 0 ) then
            call endrun (sub//': ERROR: allocate failed for '//pbuf(i)%name)
         end if
         pbuf(i)%fld_ptr = nan
         
      end if
   end do

end subroutine pbuf_allocate

!=========================================================================================

subroutine pbuf_deallocate(scope)

! Deallocate storage for fields in the physics buffer with the specified scope.
! If global_allocate_all=.true. then storage for both global and physpkg scope 
! is deallocated just once, when scope='global'.

   implicit none

   character(len=*), intent(in)  :: scope

! Local variables
   character(len=*), parameter :: sub = 'pbuf_deallocate'
   integer :: i, fdim, mdim, ldim
   logical :: deallocate_all
!-----------------------------------------------------------------------------------------

   ! deallocate_all is used to force allocation of all fields at same time as allocation
   ! for global scope.
   deallocate_all = .false.
   if ( global_allocate_all ) then
      if ( scope == 'global' ) then
         deallocate_all = .true.
      else
         return
      end if
   end if

   do i = 1, pbuf_size
      if ( pbuf(i)%scope == scope   .or.  deallocate_all ) then
         if (associated(pbuf(i)%fld_ptr)) then
            deallocate(pbuf(i)%fld_ptr)
         else
            call endrun (sub//': ERROR: '//pbuf(i)%name//' is not allocated')
         end if
      end if
   end do

end subroutine pbuf_deallocate

!=========================================================================================

subroutine pbuf_setval(name, value)

! Set a value for a field in the physics buffer.

   implicit none

   character(len=*), intent(in)  :: name   ! field name 
   real(r8),         intent(in)  :: value

! Local variables
   character(len=*), parameter :: sub = 'pbuf_setval'
   integer :: idx, lchnk, ncols
   integer :: failcode = -1
!-----------------------------------------------------------------------------------------

   ! get field index
   idx = pbuf_get_fld_idx(name, failcode)
   ! check return value
   if (idx == failcode) then
      call endrun (sub//': ERROR: name not found:'//name)
   end if

   ! check that field pointer is associated -- if so then set the value
   if (associated(pbuf(idx)%fld_ptr)) then
      do lchnk = begchunk, endchunk
         ncols = get_ncols_p(lchnk)
         pbuf(idx)%fld_ptr(:,1:ncols,:,lchnk,:) = value
      end do
   else
      call endrun (sub//': ERROR: field '//name//' is not allocated')
   end if

end subroutine pbuf_setval

!=========================================================================================

subroutine pbuf_init_restart(file, hdimids, hdim1_d, hdim2_d, pver_id, pverp_id)

  use ppgrid, only : pver, pverp
  type(file_desc_t), intent(inout) :: file
  integer, intent(in) :: hdimids(:)
  integer, intent(in) :: hdim1_d, hdim2_d, pver_id, pverp_id

  integer :: ncdimid, londimid, latdimid
  integer :: dimids(5), ndims, dimlens(5)
  integer :: i, j, ierr, hdimcnt
  integer :: totalsize
  integer, pointer :: ldof(:)
  character(len=20) :: dimname, tmpname


  hdimcnt = size(hdimids)

  dimlens(1)=hdim1_d
  dimlens(2)=hdim2_d

  do i = 1, pbuf_size
     nullify(pbuf(i)%vardesc_list)

     ndims=hdimcnt
     dimids(1:ndims) = hdimids
     totalsize = 1
     do j = 1, ndims
        totalsize = totalsize*dimlens(j)
     enddo

     if ( pbuf(i)%scope == 'global' ) then
        if(pbuf(i)%fdim>1) then
           dimname =  trim(pbuf(I)%name)//'_fdim'
           ndims=ndims+1
           ierr = pio_def_dim(file, dimname, pbuf(i)%fdim, dimids(ndims))
           dimlens(ndims)=pbuf(i)%fdim
           totalsize = totalsize*dimlens(ndims)
        end if
!
! Assume that if the mdim is the same as pver or pverp then it is pver or pverp.  
! This is almost always the case and should not hurt if it isnt
        if(pbuf(i)%mdim>1) then
           ndims=ndims+1
	   if(pbuf(i)%mdim==pver) then
              dimids(ndims) = pver_id
              dimlens(ndims)=pver
           else if(pbuf(i)%mdim==pverp) then
              dimids(ndims) = pverp_id
              dimlens(ndims)=pverp
           else
              dimname =  trim(pbuf(I)%name)//'_mdim'
              ierr = pio_def_dim(file, trim(dimname), pbuf(i)%mdim, dimids(ndims))
              dimlens(ndims)=pbuf(i)%mdim
           end if
           totalsize = totalsize*dimlens(ndims)
        end if
        if(pbuf(i)%ldim>1) then
           totalsize=totalsize*pbuf(i)%ldim

! max bytes per field is 4GB for netcdf, 8 bytes per element for double 
           if(totalsize < maxvarsize) then
              dimname =  trim(pbuf(I)%name)//'_ldim'
              ndims=ndims+1
              ierr = pio_def_dim(file, dimname,  pbuf(i)%ldim, dimids(ndims))
              dimlens(ndims)=pbuf(i)%ldim
           else
              allocate(pbuf(i)%vardesc_list(pbuf(i)%ldim))
           end if
        end if
        if(totalsize < maxvarsize) then
           ierr = pio_def_var(file, trim(pbuf(i)%name), pio_double, dimids(1:ndims), pbuf(i)%vardesc)
        else
           if(pbuf(i)%ldim==1) then
              call endrun('size of var exceeds netcdf limitations')
           end if
           do j=1,pbuf(i)%ldim
              write(tmpname,'(a,i3.3)') trim(pbuf(i)%name),j
              ierr = pio_def_var(file, tmpname, pio_double, dimids(1:ndims), pbuf(i)%vardesc_list(j))
           end do
        end if
     end if
  end do
  ierr = pio_def_var(File, 'pbuf_time_idx', pio_int, timeidx_desc)

end subroutine pbuf_init_restart


subroutine pbuf_write_restart(File)

! write physics buffer to restart file
  use ppgrid,        only: pver, pverp
  use cam_pio_utils, only: get_phys_decomp
   
  type(file_desc_t), intent(inout) :: file
   
! Local variables
  character(len=*), parameter :: sub = 'pbuf_write_restart'
  integer :: i, j, ierr
  real(r8) :: mold(1), null(0)  ! required for transfer statement
!-----------------------------------------------------------------------------------------

  do i = 1, pbuf_size
     if ( pbuf(i)%scope == 'global' ) then
        if(associated(pbuf(i)%vardesc_list)) then
           call get_phys_decomp(pbuf(i)%iodesc, pbuf(i)%fdim, pbuf(i)%mdim, 1,pio_double)
           do j=1,pbuf(i)%ldim
              call pio_write_darray(file, pbuf(i)%vardesc_list(j), pbuf(i)%iodesc, pbuf(i)%fld_ptr(:,:,:,:,j), ierr)
           end do
        else
           call get_phys_decomp(pbuf(i)%iodesc, pbuf(i)%fdim, pbuf(i)%mdim, pbuf(i)%ldim,pio_double)
           call pio_write_darray(file, pbuf(i)%vardesc, pbuf(i)%iodesc, pbuf(i)%fld_ptr, ierr)
        end if
     end if
  end do

  ierr = pio_put_var(File, timeidx_desc, (/old_time_idx/))
  

end subroutine pbuf_write_restart

!=========================================================================================

subroutine pbuf_read_restart(File)
  use ppgrid,        only: pver, pverp
  use cam_pio_utils, only: get_phys_decomp, fillvalue
  
! write physics buffer to restart file

   type(file_desc_t), intent(inout) :: File
   type(var_desc_t) :: vardesc

! Local variables
   character(len=*), parameter :: sub = 'pbuf_read_restart'
   integer :: i, ierr
   real(r8), allocatable :: tmpfld(:)
   integer :: fdim, mdim, ldim, j
   character(len=20) :: tmpname
!-----------------------------------------------------------------------------------------
   ierr = pio_inq_varid(File, 'pbuf_time_idx', vardesc)
   ierr = pio_get_var(File, vardesc, old_time_idx) 

   do i = 1, pbuf_size
      if ( pbuf(i)%scope == 'global' ) then

         call get_phys_decomp(pbuf(i)%iodesc, pbuf(i)%fdim, pbuf(i)%mdim, pbuf(i)%ldim,pio_double)

         allocate(tmpfld(size(pbuf(i)%fld_ptr)))
	 tmpfld(:) = fillvalue
         fdim = pbuf(i)%fdim
         mdim = pbuf(i)%mdim
         ldim = pbuf(i)%ldim

         call pio_seterrorhandling(File, pio_bcast_error)
         ierr = pio_inq_varid(File, trim(pbuf(i)%name), vardesc)
         if(ierr==PIO_NOERR) then
            call pio_seterrorhandling(File, pio_internal_error)
            call pio_read_darray(File, vardesc, pbuf(i)%iodesc, tmpfld, ierr)
            pbuf(i)%fld_ptr = reshape(tmpfld,(/fdim, pcols, mdim, endchunk-begchunk+1, ldim/))
         else if (ldim>1) then
            write(tmpname,'(a,i3.3)') trim(pbuf(i)%name),1
            ierr = pio_inq_varid(File, tmpname, vardesc)
            if(ierr==PIO_NOERR) then
               call pio_seterrorhandling(File, pio_internal_error)
               call pio_read_darray(File, vardesc, pbuf(i)%iodesc, tmpfld, ierr)
               pbuf(i)%fld_ptr(:,:,:,:,1) = reshape(tmpfld,(/fdim, pcols, mdim, endchunk-begchunk+1/))
               do j=2,ldim
                  write(tmpname,'(a,i3.3)') trim(pbuf(i)%name),j
                  ierr = pio_inq_varid(File, tmpname, vardesc)
                  call pio_read_darray(File, vardesc, pbuf(i)%iodesc, tmpfld, ierr)
                  pbuf(i)%fld_ptr(:,:,:,:,j) = reshape(tmpfld,(/fdim, pcols, mdim, endchunk-begchunk+1/))
               end do
            end if
            if(ierr/=PIO_NOERR) then
               write(iulog,*) __FILE__,__LINE__,'error getting ', pbuf(i)%name, file%fh
            end if
         end if
         deallocate(tmpfld)
      end if
   end do

end subroutine pbuf_read_restart


end module phys_buffer
