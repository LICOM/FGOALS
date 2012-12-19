!-------------------------------------------------------------------
! manages reading and interpolation of prescribed aerosol deposition fluxes
! Created by: Francis Vitt
!-------------------------------------------------------------------
module aerodep_flx

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile
  use cam_logfile,  only : iulog

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: aerodep_flx_init
  public :: aerodep_flx_adv
  public :: aerodep_flx_readnl

  logical :: has_aerodep_flx = .false.
  integer, parameter, public :: N_DEPFLX = 14
  integer :: number_flds

  character(len=256) :: filename = ' '
  character(len=256) :: filelist = ' '
  character(len=256) :: datapath = ' '
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  character(len=16)  :: specifier(N_DEPFLX) = ' '
  
  character(len=8), parameter :: flx_names(N_DEPFLX) = (/ &
       'BCDEPWET', 'BCPHODRY', 'BCPHIDRY',  &
       'OCDEPWET', 'OCPHODRY', 'OCPHIDRY',  &
       'DSTX01DD', 'DSTX02DD', 'DSTX03DD', 'DSTX04DD', &
       'DSTX01WD', 'DSTX02WD', 'DSTX03WD', 'DSTX04WD'  /)

  integer:: index_map(N_DEPFLX)
  integer :: ibcphiwet,ibcphidry,ibcphodry
  integer :: iocphiwet,iocphidry,iocphodry

  integer :: idstdry1,idstdry2,idstdry3,idstdry4
  integer :: idstwet1,idstwet2,idstwet3,idstwet4

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine aerodep_flx_init()
    
    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk

    implicit none

    integer :: ndx, istat, i

    if ( has_aerodep_flx ) then
       if ( masterproc ) then
          write(iulog,*) 'aero dep fluxes are prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    if ( len_trim(specifier(1)) == 0 ) specifier = flx_names

    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype )

    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       has_aerodep_flx = .false.
       if (masterproc) then
          write(iulog,*) 'aerodep_flx_init: no aerosol deposition fluxes have been specified'
       endif
       return
    end if

    index_map(:) = -1

    do i = 1,number_flds
       ndx = get_ndx( fields(i)%fldnam )

       if (ndx < 1) then
          call endrun('aerodep_flx_init: '//trim(fields(i)%fldnam)//' is not one of the named aerosol fields in pbuf')
       endif
       index_map(ndx) = i

       call addfld(fields(i)%fldnam,'kg/m2/sec', 1, 'I', 'prescribed aero dep', phys_decomp )
    enddo

    ibcphiwet = index_map(1)
    ibcphodry = index_map(2)
    ibcphidry = index_map(3)

    iocphiwet = index_map(4)
    iocphodry = index_map(5)
    iocphidry = index_map(6)

    idstdry1 = index_map(7)
    idstdry2 = index_map(8)
    idstdry3 = index_map(9)
    idstdry4 = index_map(10)

    idstwet1 = index_map(11)
    idstwet2 = index_map(12)
    idstwet3 = index_map(13)
    idstwet4 = index_map(14)

  end subroutine aerodep_flx_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine aerodep_flx_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'aerodep_flx_readnl'

   character(len=16)  :: aerodep_flx_specifier(N_DEPFLX)
   character(len=256) :: aerodep_flx_file
   character(len=256) :: aerodep_flx_filelist
   character(len=256) :: aerodep_flx_datapath
   character(len=32)  :: aerodep_flx_type
   logical            :: aerodep_flx_rmfile
   integer            :: aerodep_flx_cycle_yr
   integer            :: aerodep_flx_fixed_ymd
   integer            :: aerodep_flx_fixed_tod

   namelist /aerodep_flx_nl/ &
      aerodep_flx_specifier, &
      aerodep_flx_file,      &
      aerodep_flx_filelist,  &
      aerodep_flx_datapath,  &
      aerodep_flx_type,      &
      aerodep_flx_rmfile,    &
      aerodep_flx_cycle_yr,  &
      aerodep_flx_fixed_ymd, &
      aerodep_flx_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   aerodep_flx_specifier= specifier
   aerodep_flx_file     = filename
   aerodep_flx_filelist = filelist
   aerodep_flx_datapath = datapath
   aerodep_flx_type     = datatype
   aerodep_flx_rmfile   = rmv_file
   aerodep_flx_cycle_yr = cycle_yr
   aerodep_flx_fixed_ymd= fixed_ymd
   aerodep_flx_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'aerodep_flx_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, aerodep_flx_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(aerodep_flx_specifier,len(aerodep_flx_specifier(1))*N_DEPFLX,     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_file,     len(aerodep_flx_file),     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_filelist, len(aerodep_flx_filelist), mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_datapath, len(aerodep_flx_datapath), mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_type,     len(aerodep_flx_type),     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(aerodep_flx_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(aerodep_flx_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(aerodep_flx_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier  = aerodep_flx_specifier
   filename   = aerodep_flx_file
   filelist   = aerodep_flx_filelist
   datapath   = aerodep_flx_datapath
   datatype   = aerodep_flx_type
   rmv_file   = aerodep_flx_rmfile
   cycle_yr   = aerodep_flx_cycle_yr
   fixed_ymd  = aerodep_flx_fixed_ymd
   fixed_tod  = aerodep_flx_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_aerodep_flx = .true.

end subroutine aerodep_flx_readnl


!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine aerodep_flx_set( cam_out, ncol, lchnk )
    use camsrfexch_types, only : cam_out_t     

    type(cam_out_t),     intent(inout) :: cam_out
    integer,             intent(in)    :: ncol, lchnk
    
    if( .not. has_aerodep_flx ) return

    call set_fluxes( cam_out%bcphiwet, ibcphiwet, ncol, lchnk )
    call set_fluxes( cam_out%bcphidry, ibcphidry, ncol, lchnk )
    call set_fluxes( cam_out%bcphodry, ibcphodry, ncol, lchnk )

    call set_fluxes( cam_out%ocphiwet, iocphiwet, ncol, lchnk )
    call set_fluxes( cam_out%ocphidry, iocphidry, ncol, lchnk )
    call set_fluxes( cam_out%ocphodry, iocphodry, ncol, lchnk )

    call set_fluxes( cam_out%dstdry1, idstdry1, ncol, lchnk )
    call set_fluxes( cam_out%dstdry2, idstdry2, ncol, lchnk )
    call set_fluxes( cam_out%dstdry3, idstdry3, ncol, lchnk )
    call set_fluxes( cam_out%dstdry4, idstdry4, ncol, lchnk )

    call set_fluxes( cam_out%dstwet1, idstwet1, ncol, lchnk )
    call set_fluxes( cam_out%dstwet2, idstwet2, ncol, lchnk )
    call set_fluxes( cam_out%dstwet3, idstwet3, ncol, lchnk )
    call set_fluxes( cam_out%dstwet4, idstwet4, ncol, lchnk )

  end subroutine aerodep_flx_set

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine set_fluxes( fluxes, fld_indx, ncol, lchnk )
    use cam_history,  only : outfld

    real(r8), intent(inout) :: fluxes(:)
    integer,  intent(in)    :: fld_indx, ncol, lchnk

    integer :: i

    if (fld_indx<1) return

    do i = 1,ncol
       fluxes(i) = max(fields(fld_indx)%data(i,1,lchnk), 0._r8)
    enddo

    call outfld(fields(fld_indx)%fldnam, fluxes(:ncol), ncol, lchnk )

  endsubroutine set_fluxes

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine aerodep_flx_adv( state, cam_out )

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use camsrfexch_types, only : cam_out_t

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)

    integer :: c, ncol
    
    if( .not. has_aerodep_flx ) return

    call advance_trcdata( fields, file, state  )

!$OMP PARALLEL DO PRIVATE (C, NCOL)
    do c = begchunk, endchunk
       ncol = state(c)%ncol
       call aerodep_flx_set( cam_out(c), ncol, c )
    enddo

  end subroutine aerodep_flx_adv

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  integer function get_ndx( name )

    implicit none
    character(len=*), intent(in) :: name

    integer :: i

    get_ndx = 0
    do i = 1,N_DEPFLX
      if ( trim(name) == trim(flx_names(i)) ) then
        get_ndx = i
        return
      endif
    enddo

  end function get_ndx

end module aerodep_flx
