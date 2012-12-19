#define _FILE 'control/readinitial.F90 '
module readinitial

  public :: read_initial

contains

  subroutine read_initial(ncid)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Ensure that requisite netcdf variables are on the initial dataset.
    !          Set base day and date info using the "current" values from it.
    ! 
    ! Method: Issue proper netcdf wrapper calls.  Broadcast to slaves if SPMD
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use pmgrid,         only: plev, plat, plon, plevp
    use pspect
    use rgrid
    use abortutils,     only: endrun
    use hycoef,         only: read_restart_hycoef
    use error_messages, only: handle_ncerr
    use scamMod,        only: single_column   
    use cam_logfile,    only: iulog
    use pio, only : pio_inq_dimlen, pio_inq_varid,  pio_inq_dimid, &
         var_desc_t, pio_internal_error, pio_noerr, pio_bcast_error, &
         pio_seterrorhandling, file_desc_t
    use dycore, only : dycore_is
    !-----------------------------------------------------------------------
    implicit none
    !------------------------------Parameters-------------------------------
    !
    ! Arguments
    !
    type(file_desc_t), intent(inout) :: ncid  ! History file unit
    !                                
    ! Local variables
    !                                
    integer :: j
    integer :: lonid            !------------------------------------------------------
    integer :: levid            ! 
    integer :: latid            ! 
    type(var_desc_t) :: ntrkid           ! 
    type(var_desc_t) :: ntrmid           ! 
    type(var_desc_t) :: ntrnid           ! 
    type(var_desc_t) :: ncdateid         ! PIO variable and dimension ids for variable of that
    type(var_desc_t) :: ncsecid          ! name with "id" tacked on to the end
    type(var_desc_t) :: hyaiid           ! 
    type(var_desc_t) :: hybiid           ! 
    type(var_desc_t) :: hyamid           ! 
    type(var_desc_t) :: hybmid           ! 
    integer :: ilev             ! 
    integer :: ilevid           ! 
    integer :: rlonid           ! 
    type(var_desc_t) :: nlonid           ! 
    type(var_desc_t) :: wnummaxid        !------------------------------------------------------

    integer :: mlon             ! longitude dimension length from dataset
    integer :: mlev             ! level dimension length from dataset
    integer :: morec            ! latitude dimension length from dataset

    integer :: slonid           ! staggered longitude dimension length from dataset
    integer :: slatid           ! staggered latitude dimension length from dataset

    integer :: ncollen
    integer :: ierr
    logical :: isncol
    !
    !-----------------------------------------------------------------------
    !
    isncol = .false.

    !  Tells pio to communicate any error to all tasks so that it can be handled locally
    !  otherwise and by default pio error are handled internally.

    call pio_seterrorhandling( ncid, PIO_BCAST_ERROR)
    ierr = PIO_inq_dimid( ncid, 'ncol', lonid)
    call pio_seterrorhandling( ncid, PIO_INTERNAL_ERROR)

    if(ierr==PIO_NOERR) then
       isncol=.true.
       ierr = PIO_inq_dimlen (ncid, lonid , ncollen)
    else

       !
       ! Get and check dimension/date info
       !
       ierr = PIO_inq_dimid (ncid, 'lon' , lonid)
       ierr = PIO_inq_dimid (ncid, 'lat' , latid)

       if(dycore_is('LR')) then
          ierr = PIO_inq_dimid (ncid, 'slon' , slonid)
          ierr = PIO_inq_dimid (ncid, 'slat' , slatid)
       end if
       ierr = PIO_inq_dimlen (ncid, lonid , mlon)
       ierr = PIO_inq_dimlen (ncid, latid , morec)
    endif
    ierr = PIO_inq_dimid (ncid, 'lev' , levid)
    ierr = PIO_inq_dimid (ncid, 'ilev', ilevid)
    !
    ierr = PIO_inq_dimlen (ncid, levid , mlev)
    ierr = PIO_inq_dimlen (ncid, ilevid, ilev)
    !
    ierr = PIO_inq_varid (ncid, 'ntrk'   , ntrkid)
    ierr = pio_inq_varid (ncid, 'ntrm'   , ntrmid)
    ierr = pio_inq_varid (ncid, 'ntrn'   , ntrnid)
    ierr = pio_inq_varid (ncid, 'date'   , ncdateid)
    ierr = pio_inq_varid (ncid, 'datesec', ncsecid)
    !
    ! Reduced grid is not supported
    !
    wnummax(:) = ptrm
    nlon   (:) = plon


    if(isncol) then
       if(mlev /= plev ) then
          write(iulog,*)'READINITIAL: model parameters do not match initial dataset parameters'
          write(iulog,*)'Model Parameters:   plev = ',plev
          write(iulog,*)'Dataset Parameters: dlev = ',mlev
          call endrun(_FILE)
       end if
    else if (single_column) then
       if (mlev /= plev) then
          write(iulog,*)'READINITIAL: model parameters for single column run do not match initial dataset parameters'
          write(iulog,*)'Model Parameters:   plev = ',plev,' plon = ',plon,' plat = ',plat
          write(iulog,*)'Dataset Parameters: dlev = ',mlev,' dlon = ',mlon,' dlat = ',morec
          call endrun(_FILE)
       endif
    else if (mlev /= plev.or.mlon /= plon.or.morec /= plat) then
       write(iulog,*)'READINITIAL: model parameters do not match initial dataset parameters'
       write(iulog,*)'Model Parameters:   plev = ',plev,' plon = ',plon,' plat = ',plat
       write(iulog,*)'Dataset Parameters: dlev = ',mlev,' dlon = ',mlon,' dlat = ',morec
       call endrun(_FILE)
    end if

    call read_restart_hycoef(ncid)

    return
  end subroutine read_initial



end module readinitial
