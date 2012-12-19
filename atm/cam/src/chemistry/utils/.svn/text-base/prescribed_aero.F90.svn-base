!-------------------------------------------------------------------
! manages reading and interpolation of prescribed aerosols
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_aero

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

  public :: prescribed_aero_init
  public :: prescribed_aero_adv
  public :: write_prescribed_aero_restart
  public :: read_prescribed_aero_restart
  public :: has_prescribed_aero
  public :: prescribed_aero_register
  public :: init_prescribed_aero_restart
  public :: prescribed_aero_readnl

  logical :: has_prescribed_aero = .false.
  integer, parameter, public :: N_AERO = 13
  integer :: number_flds

  character(len=256) :: filename = ' '
  character(len=256) :: filelist = ' '
  character(len=256) :: datapath = ' '
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  character(len=16)  :: specifier(N_AERO) = ''

  real(r8), parameter :: molmass(N_AERO) = (/ 96.0635986_r8, &
                          12.0109997_r8, 12.0109997_r8, 12.0109997_r8, 12.0109997_r8, &
                          58.4424667_r8, 58.4424667_r8, 58.4424667_r8, 58.4424667_r8, &
                          135.064041_r8, 135.064041_r8, 135.064041_r8, 135.064041_r8 /)

  character(len=8)    :: aero_names(N_AERO) = (/ 'sulf    ', &
                          'bcar1   ',    'bcar2   ',    'ocar1   ',    'ocar2   ', &
                          'sslt1   ',    'sslt2   ',    'sslt3   ',    'sslt4   ', &
                          'dust1   ',    'dust2   ',    'dust3   ',    'dust4   ' /)

  integer :: index_map(N_AERO)

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_aero_register()
    use ppgrid,         only: pver
    use phys_buffer,    only: pbuf_add
    integer :: i,idx

    if (has_prescribed_aero) then
       do i = 1,N_AERO
          call pbuf_add(aero_names(i),'physpkg',1,pver,1,idx)
       enddo
    endif

  endsubroutine prescribed_aero_register
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_aero_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk

    implicit none

    integer :: ndx, istat, i
    
    if ( has_prescribed_aero ) then
       if ( masterproc ) then
          write(iulog,*) 'aero is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    file%in_pbuf = .true.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype )
        
    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'There are no prescribed aerosols'
          write(iulog,*) ' '
       endif
       return
    end if

    do i = 1,number_flds
       ndx = get_ndx( fields(i)%fldnam )
       index_map(i) = ndx

       if (ndx < 1) then
          call endrun('prescribed_aero_init: '//trim(fields(i)%fldnam)//' is not one of the named aerosol fields in pbuf')
       endif
       call addfld(aero_names(i),'kg/kg', pver, 'I', 'prescribed aero', phys_decomp )
    enddo


  end subroutine prescribed_aero_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_aero_readnl'

   character(len=16)  :: prescribed_aero_specifier(N_AERO)
   character(len=256) :: prescribed_aero_file
   character(len=256) :: prescribed_aero_filelist
   character(len=256) :: prescribed_aero_datapath
   character(len=32)  :: prescribed_aero_type
   logical            :: prescribed_aero_rmfile
   integer            :: prescribed_aero_cycle_yr
   integer            :: prescribed_aero_fixed_ymd
   integer            :: prescribed_aero_fixed_tod

   namelist /prescribed_aero_nl/ &
      prescribed_aero_specifier, &
      prescribed_aero_file,      &
      prescribed_aero_filelist,  &
      prescribed_aero_datapath,  &
      prescribed_aero_type,      &
      prescribed_aero_rmfile,    &
      prescribed_aero_cycle_yr,  &
      prescribed_aero_fixed_ymd, &
      prescribed_aero_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_aero_specifier= specifier
   prescribed_aero_file     = filename
   prescribed_aero_filelist = filelist
   prescribed_aero_datapath = datapath
   prescribed_aero_type     = datatype
   prescribed_aero_rmfile   = rmv_file
   prescribed_aero_cycle_yr = cycle_yr
   prescribed_aero_fixed_ymd= fixed_ymd
   prescribed_aero_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_aero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_aero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_aero_specifier,len(prescribed_aero_specifier(1))*N_AERO,     mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_file,     len(prescribed_aero_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_filelist, len(prescribed_aero_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_datapath, len(prescribed_aero_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_type,     len(prescribed_aero_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_aero_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_aero_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_aero_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_aero_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier  = prescribed_aero_specifier
   filename   = prescribed_aero_file
   filelist   = prescribed_aero_filelist
   datapath   = prescribed_aero_datapath
   datatype   = prescribed_aero_type
   rmv_file   = prescribed_aero_rmfile
   cycle_yr   = prescribed_aero_cycle_yr
   fixed_ymd  = prescribed_aero_fixed_ymd
   fixed_tod  = prescribed_aero_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_aero = .true.

end subroutine prescribed_aero_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_aero_adv( state )

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : amass => mwdry       ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz                ! J/K/molecule

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    integer :: ind,c,ncol,i,caseid,chnk_offset
    real(r8) :: to_mmr(pcols,pver)

    if( .not. has_prescribed_aero ) return

    call advance_trcdata( fields, file, state  )
    
    ! set the tracer fields with the correct units
    do i = 1,number_flds

       select case ( to_lower(trim(fields(i)%units(:GLC(fields(i)%units)))) )
       case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
          caseid = 1
       case ('kg/kg','mmr')
          caseid = 2
       case ('mol/mol','mole/mole','vmr','fraction')
          caseid = 3
       case default
          print*, 'prescribed_aero_adv: units = ',trim(fields(i)%units) ,' are not recognized'
          call endrun('prescribed_aero_adv: units are not recognized')
       end select

       ind = index_map(i)
       chnk_offset = fields(i)%chnk_offset

!$OMP PARALLEL DO PRIVATE (C, NCOL, TO_MMR)
       do c = begchunk,endchunk
          ncol = state(c)%ncol

          if (caseid == 1) then
             to_mmr(:ncol,:) = (molmass(ind)*1.e6_r8*boltz*state(c)%t(:ncol,:))/(amass*state(c)%pmiddry(:ncol,:))
          elseif (caseid == 2) then
             to_mmr(:ncol,:) = 1._r8
          else
             to_mmr(:ncol,:) = molmass(ind)/amass
          endif

          fields(i)%data(:ncol,:,c+chnk_offset) = to_mmr(:ncol,:) * fields(i)%data(:ncol,:,c+chnk_offset) 
          call outfld( fields(i)%fldnam, fields(i)%data(:ncol,:,c+chnk_offset), ncol, state(c)%lchnk )
       enddo
    enddo

  end subroutine prescribed_aero_adv

!-------------------------------------------------------------------


!-------------------------------------------------------------------
  subroutine init_prescribed_aero_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'prescribed_aero', piofile, file )

  end subroutine init_prescribed_aero_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_aero_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_aero_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_prescribed_aero_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_aero', piofile, file )

  end subroutine read_prescribed_aero_restart
!-------------------------------------------------------------------
  integer function get_ndx( name )

    implicit none
    character(len=*), intent(in) :: name

    integer :: i

    get_ndx = 0
    do i = 1,N_AERO
      if ( trim(name) == trim(aero_names(i)) ) then
        get_ndx = i
        return
      endif
    enddo

  end function get_ndx

end module prescribed_aero
