!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: ice_data
!
! !DESCRIPTION:	Module to handle dealing with the ICE data.
!
! Public interfaces:
!
!	iceini -- Initialization and reading of dataset.
!	iceint -- Interpolate dataset ICE to current time.
!
!----------------------------------------------------------------------- 

module ice_data

#if ( ! defined COUP_SOM )
!
! USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_scam_mod, only: shr_scam_GetCloseLatLon
  use pmgrid,         only: plon, plat
  use ppgrid,         only: pcols, begchunk, endchunk
  use phys_grid,      only: scatter_field_to_chunk, get_ncols_p,get_lat_all_p
  use commap,         only: clat, clon
  use abortutils,     only: endrun
  use wrap_nf
  use cam_logfile,    only: iulog

  use ice_constants
  use ice_kinds_mod, only: log_kind
  use ice_dh,        only: prognostic_icesnow
  use ice_spmd

  use scamMod,       only :scmlon,scmlat,isrestart,single_column

  implicit none
!----------------------------------------------------------------------- 
! PUBLIC: Make default data and interfaces private
!----------------------------------------------------------------------- 
!
! ! PUBLIC MEMBER FUNCTIONS:
!
  public iceini   ! Initialization
  public iceint   ! Time interpolation of ICE data
  logical (kind=log_kind), parameter :: snowice_climatology = .true.
!
! ! PUBLIC DATA:
!
  logical, public :: icecyc            ! true => cycle ice fraction dataset

!===============================================================================
!EOP
!===============================================================================
!----------------------------------------------------------------------- 
! PRIVATE: Everthing else is private to this module
!----------------------------------------------------------------------- 
  private   ! By default all data is private to this module
  integer, parameter :: toticesz=2000
  real(r8), parameter :: daysperyear = 365.0_dbl_kind  ! Number of days in a year

  real(r8), allocatable, dimension(:,:,:) :: &
      icebdy         ! ICE values on boundary dataset (pcols,begchunk:endchunk,2)
  real(r8), allocatable, dimension(:,:) :: &
      ice            ! Interpolated model ice values (pcols,begchunk:endchunk)

  real(r8) :: cdayicem   ! Calendar day for prv. month ICE values read in
  real(r8) :: cdayicep   ! Calendar day for nxt. month ICE values read in
      
  integer :: nm,np   ! Array indices for prv., nxt month ice data
  integer :: nm1
  integer :: nmshift,npshift ! Array indices for prv., nxt month ice data
  integer :: iceid   ! netcdf id for ice variable
  integer :: lonsiz  ! size of longitude dimension on ice dataset
  integer :: levsiz  ! size of level dimension on ice dataset
  integer :: latsiz  ! size of latitude dimension on ice dataset
  integer :: timesiz ! size of time dimension on ice dataset
  integer :: np1     ! current forward time index of ice dataset
  integer :: date_ice(toticesz)! Date on ice dataset (YYYYMMDD)
  integer :: sec_ice(toticesz) ! seconds of date on ice dataset (0-86399)

  integer :: closelatidx,closelonidx
  real(r8):: closelat,closelon
 
  real(r8):: snwbdynh(2)   ! snow height boundary values nrthrn hmsphr
  real(r8):: snwbdysh(2)   ! snow height boundary values sthrn hmsphr
  real(r8):: snwcnh(12)   ! mean snow cover (m) first of month (nrthrn hmsphr)
  real(r8):: snwcsh(12)   ! mean snow cover (m) first of month (sthrn hmsphr)
  data snwcnh   /  .23_dbl_kind,  .25_dbl_kind,  .27_dbl_kind,  .29_dbl_kind,  .33_dbl_kind,  .18_dbl_kind, &
                 0._dbl_kind,   0._dbl_kind,  .02_dbl_kind,  .12_dbl_kind,  .18_dbl_kind,  .21_dbl_kind / 
  data snwcsh   /  0._dbl_kind,   0._dbl_kind,  .02_dbl_kind,  .12_dbl_kind,  .18_dbl_kind,  .21_dbl_kind, & 
                 .23_dbl_kind,  .25_dbl_kind,  .27_dbl_kind,  .29_dbl_kind,  .33_dbl_kind,  .18_dbl_kind/

  integer :: ncid_sst ! netcdf file handle for sst file

  save
!===============================================================================
CONTAINS
!===============================================================================

!======================================================================
! PUBLIC ROUTINES: Following routines are publically accessable
!======================================================================
!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: iceini
!
! !DESCRIPTION:
!
! Initialize the procedure for specifying sea surface temperatures
! Do initial read of time-varying ice boundary dataset, reading two
! consecutive months on either side of the current model date.
!
! Method: 
! 
! Author: L.Bath
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
subroutine iceini(bndtvs)
!
! !USES:
!
  use rgrid,            only: nlon, fullgrid
  use error_messages,   only: alloc_err, handle_ncerr
  use ioFileMod,        only: getfil
  use ice_time_manager, only: get_curr_date, get_curr_calday, get_step_size, &
                              is_first_step
  use ice_constants
!
! EOP

!
!---------------------------Input Arguments-----------------------------

  character(len=*), intent(in) :: bndtvs       !  Input character string

!
!---------------------------Local variables-----------------------------
  integer dtime                 ! timestep size [seconds]
  integer dateid                ! netcdf id for date variable
  integer secid                 ! netcdf id for seconds variable
  integer londimid              ! netcdf id for longitude variable
  integer latdimid              ! netcdf id for latitude variable
  integer lonid                 ! netcdf id for longitude variable
  integer latid                 ! netcdf id for latitude variable
  integer timeid                ! netcdf id for time variable
  integer nlonid                ! netcdf id for nlon variable (rgrid)
  integer cnt3(3)               ! array of counts for each dimension
  integer strt3(3)              ! array of starting indices
  integer n                     ! indices
  integer nlon_ice(plat)        ! number of lons per lat on bdy dataset
  integer i                     ! index into chunk
  integer j                     ! latitude index
  integer k
  integer ncol
  integer istat                 ! error return
  integer lchnk           ! chunk to process
  integer  :: yr, mon, day      ! components of a date
  integer  :: ncdate            ! current date in integer format [yyyymmdd]
  integer  :: ncsec             ! current time of day [seconds]
  real(r8) calday               ! calendar day (includes yr if no cycling)
  real(r8) caldayloc            ! calendar day (includes yr if no cycling)
  real(r8) xvar(plon,plat,2)    ! work space 
  real(r8) xvartmp    ! work space 
  integer  :: ret
  character(len=256) :: locfn   ! netcdf local filename to open
!-----------------------------------------------------------------------
  snwbdynh(:) = 0.0_dbl_kind
  snwbdysh(:) = 0.0_dbl_kind

  !
  ! Obtain time-variant sst datatset
  !
  if (masterproc) then
     call getfil(bndtvs, locfn)
     call wrap_open(locfn, 0, ncid_sst)
     write(iulog,*)'ICEINI: NCOPN returns id ',ncid_sst,' for file ',trim(locfn)
  endif
!
  if (.not. prognostic_icesnow) then
     if (.not.icecyc) then
        call endrun ('ice_data: Snowice climatology option not valid for icecyc = .false.')
     end if
  end if
! initialize ice constants	

  call init_constants

!
! Initialize time indices
!
  nm = 1
  np = 2
!
! Allocate space for data.
!
  allocate( ice(pcols,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'iceini', 'ice', &
       pcols*(endchunk-begchunk+1) )

  allocate( icebdy(pcols,begchunk:endchunk,2), stat=istat )
  call alloc_err( istat, 'iceini', 'icebdy', &
       pcols*(endchunk-begchunk+1)*2 )

  if (masterproc) then
!
! Use year information only if not cycling ice dataset
!
     if (is_first_step()) then
        dtime = get_step_size()
        dtime = -dtime
        calday = get_curr_calday(offset=dtime)
        call get_curr_date(yr, mon, day, ncsec,offset=dtime)
     else
        calday = get_curr_calday()
        call get_curr_date(yr, mon, day, ncsec)
     endif
     if (icecyc) then
        caldayloc = calday
     else
        caldayloc = calday + yr*daysperyear
     end if
     ncdate = yr*10000 + mon*100 + day
!
! Get and check dimension info
!
     call wrap_inq_dimid( ncid_sst, 'lon', londimid   )
     call wrap_inq_dimid( ncid_sst, 'time', timeid  )
     call wrap_inq_dimid( ncid_sst, 'lat', latdimid   )

     if (.not.single_column) then
        call wrap_inq_dimlen( ncid_sst, londimid, lonsiz   )
        if (lonsiz /= plon) then
           write(iulog,*)'ICEINI: lonsiz=',lonsiz,' must = plon=',plon
           call endrun
        end if

        call wrap_inq_dimlen( ncid_sst, latdimid, latsiz   )
        
        if (latsiz /= plat) then
           write(iulog,*)'ICEINI: latsiz=',latsiz,' must = plat=',plat
           call endrun
        end if
        
        call wrap_inq_dimlen( ncid_sst, timeid, timesiz   )
!
! Check to make sure space allocated for time variables is sufficient
!
        if (timesiz>toticesz) then
           write(iulog,*)'ICEINI:  Allocated space for ice data is insufficient.'
           write(iulog,*)'Please increase parameter toticesz to',timesiz,' and recompile.'
           call endrun
        end if
!
! Check to ensure reduced or not grid of dataset matches that of model
!
        if (fullgrid) then
           call wrap_inq_varid( ncid_sst, 'lon', lonid   )
        else
           call wrap_inq_varid (ncid_sst, 'nlon', nlonid)
           call wrap_get_var_int (ncid_sst, nlonid, nlon_ice)
           do j=1,plat
              if (nlon_ice(j) /= nlon(j)) then
                 write(iulog,*)'ICEINI: model grid does not match dataset grid'
                 call endrun
              end if
           end do
        end if
     else
        call wrap_inq_dimlen( ncid_sst, londimid, lonsiz   )
        call wrap_inq_dimlen( ncid_sst, latdimid, latsiz   )
        call wrap_inq_dimlen( ncid_sst, timeid, timesiz   )
        call wrap_inq_varid( ncid_sst, 'lon', lonid   )
     endif
     call wrap_inq_varid( ncid_sst, 'date', dateid   )
     call wrap_inq_varid( ncid_sst, 'datesec', secid   )
     call wrap_inq_varid( ncid_sst, 'ice_cov', iceid   )
     call wrap_inq_varid( ncid_sst, 'lat', latid   )
!
! Retrieve entire date and sec variables.
!
     call wrap_get_var_int (ncid_sst,dateid,date_ice)
     call wrap_get_var_int (ncid_sst,secid,sec_ice)
     if (icecyc) then
        if (timesiz<12) then 
           write(iulog,*)'ICEINI: ERROR' 
           write(iulog,*)'When cycling ice, ice data set must have 12' 
           write(iulog,*)'consecutive months of data starting with Jan'
           write(iulog,*)'Current dataset has only ',timesiz,' months'
           call endrun
        end if
        do n = 1,12
           if (mod(date_ice(n),10000)/100/=n) then
              write(iulog,*)'ICEINI: ERROR' 
              write(iulog,*)'When cycling ice, ice data set must have 12' 
              write(iulog,*)'consecutive months of data starting with Jan'
              write(iulog,*)'Month ',n,' of ice data set is out of order'
              call endrun
           end if
        end do
     end if

     if (single_column) then
        call shr_scam_getCloseLatLon(ncid_sst,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
        strt3(1) = closelonidx
        strt3(2) = closelatidx
        strt3(3) = 1
        cnt3(1)  = 1
        cnt3(2)  = 1
        cnt3(3)  = 1
     else
        strt3(1) = 1
        strt3(2) = 1
        strt3(3) = 1
        cnt3(1)  = lonsiz
        cnt3(2)  = latsiz
        cnt3(3)  = 1
     endif
!
! Special code for interpolation between December and January
!
     if (icecyc) then
        n = 12
        np1 = 1
        call bnddyi(date_ice(n  ), sec_ice(n  ), cdayicem)
        call bnddyi(date_ice(np1), sec_ice(np1), cdayicep)
        if (caldayloc<=cdayicep .or. caldayloc>cdayicem) then
           strt3(3) = n
           call wrap_get_vara_realx (ncid_sst,iceid,strt3,cnt3,xvar(1,1,nm))
           snwbdynh(nm)=snwcnh(n)
           snwbdysh(nm)=snwcsh(n)
           strt3(3) = np1                                      
           call wrap_get_vara_realx (ncid_sst,iceid,strt3,cnt3,xvar(1,1,np))
           snwbdynh(np)=snwcnh(np1)
           snwbdysh(np)=snwcsh(np1)
           goto 10
        end if
     end if
!
! Normal interpolation between consecutive time slices.
!
     do n=1,timesiz-1
        np1 = n + 1
        call bnddyi(date_ice(n  ), sec_ice(n  ), cdayicem)
        call bnddyi(date_ice(np1), sec_ice(np1), cdayicep)
        if (.not.icecyc) then
           yr = date_ice(n)/10000
           cdayicem = cdayicem + yr*daysperyear
           yr = date_ice(np1)/10000
           cdayicep = cdayicep + yr*daysperyear
        end if
        if (caldayloc>cdayicem .and. caldayloc<=cdayicep) then
           strt3(3) = n
           call wrap_get_vara_realx (ncid_sst,iceid,strt3,cnt3,xvar(1,1,nm))
           if (.not. prognostic_icesnow) then
              snwbdynh(nm)=snwcnh(n)
              snwbdysh(nm)=snwcsh(n)
           end if
           strt3(3) = np1                                      
           call wrap_get_vara_realx (ncid_sst,iceid,strt3,cnt3,xvar(1,1,np))
           if (.not. prognostic_icesnow) then
              snwbdynh(np)=snwcnh(np1)
              snwbdysh(np)=snwcsh(np1)
           end if
           goto 10
        end if
     end do
     write(iulog,*)'ICEINI: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
     call endrun
10   continue
     write(iulog,*)'ICEINI: Read ice data for dates ',date_ice(n),sec_ice(n), &
          ' and ',date_ice(np1),sec_ice(np1)

     call mpi_bcast( timesiz , 1       , MPI_INTEGER, 0, mpicom, ret )
     call mpi_bcast( date_ice, toticesz, MPI_INTEGER, 0, mpicom, ret )
     call mpi_bcast( sec_ice , toticesz, MPI_INTEGER, 0, mpicom, ret )
     call mpi_bcast( cdayicem, 1       , MPI_REAL8  , 0, mpicom, ret )
     call mpi_bcast( cdayicep, 1       , MPI_REAL8  , 0, mpicom, ret )
     call mpi_bcast( np1     , 1       , MPI_INTEGER, 0, mpicom, ret )
     call mpi_bcast( snwbdynh, 2       , MPI_REAL8  , 0, mpicom, ret )
     call mpi_bcast( snwbdysh, 2       , MPI_REAL8  , 0, mpicom, ret )
  else
     call mpi_bcast( timesiz , 1       , MPI_INTEGER, 0, mpicom, ret )
     call mpi_bcast( date_ice, toticesz, MPI_INTEGER, 0, mpicom, ret )
     call mpi_bcast( sec_ice , toticesz, MPI_INTEGER, 0, mpicom, ret )
     call mpi_bcast( cdayicem, 1       , MPI_REAL8  , 0, mpicom, ret )
     call mpi_bcast( cdayicep, 1       , MPI_REAL8  , 0, mpicom, ret )
     call mpi_bcast( np1     , 1       , MPI_INTEGER, 0, mpicom, ret )
     call mpi_bcast( snwbdynh, 2       , MPI_REAL8  , 0, mpicom, ret )
     call mpi_bcast( snwbdysh, 2       , MPI_REAL8  , 0, mpicom, ret )
  end if

  call scatter_field_to_chunk(1,1,2,plon,xvar,icebdy)
!
! Uncomment the following to get fractional model to give bit for bit 
! with non fractional
!
!!$  do lchnk=begchunk,endchunk
!!$     ncol = get_ncols_p(lchnk)
!!$     do i=1,ncol
!!$        do k=1,2
!!$           if (icebdy(i,lchnk,k)<=tsice) then
!!$              icebdy(i,lchnk,k) = 1.0
!!$           else
!!$              icebdy(i,lchnk,k) = 0.
!!$           end if
!!$        end do
!!$     end do
!!$  end do

  return
end subroutine iceini

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: iceint
!
! !DESCRIPTION:
!
! Otherwise, time interpolate ICE's to current time, reading in new monthly data if
! necessary.
!
! Method: 
! 
! Author: L.Bath
! 
!-----------------------------------------------------------------------
!
! !INTERFACE:
!
subroutine iceint( aice, snowhice, frac, prev_timestep)
!
! !USES:
!
  use rgrid,             only: nlon
  use ice_time_manager, only: is_first_step, get_curr_date, get_curr_calday, get_step_size
  use ice_types,         only: frac_t
!
! !PARAMETERS:
!
  real(r8), intent(inout) :: aice(pcols,begchunk:endchunk)      ! Sea-ice areal fraction over ocean-cell
  real(r8), intent(inout) :: snowhice(pcols,begchunk:endchunk)  ! Snow height over ice
  type(frac_t), intent(in):: frac(begchunk:endchunk)            ! land fraction
  logical, intent(in), optional :: prev_timestep                ! If using previous timestep, set to true
!
! EOP
!
!---------------------------Local variables-----------------------------
  integer cnt3(3)        ! array of counts for each dimension
  integer strt3(3)       ! array of starting indices
  integer i,j,lchnk      ! indices
  integer ncol           ! number of columns in current chunk
  integer ntmp           ! temporary
  integer :: dtime       ! timestep size [seconds]
  real(r8) fact1, fact2  ! time interpolation factors
  integer :: yr, mon, day! components of a date
  integer :: ncdate      ! current date in integer format [yyyymmdd]
  integer :: ncsec       ! current time of day [seconds]
  real(r8) :: calday     ! current calendar day
  real(r8) caldayloc     ! calendar day (includes yr if no cycling)
  real(r8) deltat        ! time (days) between interpolating ice data
!
! Aqua planet variables
!
  real(r8) pi            ! 3.14159...
  real(r8) pio180        ! pi/180.
  real(r8) tmp           ! temporary
  real(r8) tmp1          ! temporary
  real(r8) t0_max        ! max reference temperature
  real(r8) t0_min        ! min reference temperature
  real(r8) t0_max6       ! max asymmetric reference temperature for option 6
  real(r8) t0_max7       ! max asymmetric reference temperature for option 7
  real(r8) maxlat        ! cutoff latitude poleward of which ICE = 0 deg C
  real(r8) shift         ! number of degrees peak ICE is shifted off equator
  real(r8) shift9        ! number of degrees peak ICE is shifted off equator for opt. 9
  real(r8) shift10       ! number of degrees peak ICE is shifted off equator for opt. 10
  real(r8) latcen        ! center of asymmetric ICE forcing
  real(r8) latrad6       ! radius of asymmetric ICE forcing for option 6
  real(r8) latrad8       ! radius of asymmetric ICE forcing for option 8
  real(r8) loncen        ! center of asymmetric ICE forcing
  real(r8) lonrad        ! radius of asymmetric ICE forcing
  real(r8) xvar(plon,plat,2)    ! work space 
  integer  ice_option    ! option of analytical ICE algorithm
  integer :: lats(pcols)           ! array of latitude indices
  logical :: previous              ! If using previous timestep, set to true
  integer ret,n
  real(r8) xvartmp    ! work space 
!
!-----------------------------------------------------------------------
!
! SPMD: Master does all the work.  Sends needed info to slaves
!
!
! Use year information only if a multiyear dataset
!
     if ( .not. present(prev_timestep) ) then
        previous = .false.
     else
        previous = prev_timestep
     end if

     if (previous .and. is_first_step()) then
        dtime = get_step_size()
        dtime = -dtime
        calday = get_curr_calday(offset=dtime)
        call get_curr_date(yr, mon, day, ncsec,offset=dtime)
     else
        calday = get_curr_calday()
        call get_curr_date(yr, mon, day, ncsec)
     endif
     if (icecyc) then
        caldayloc = calday
     else
        caldayloc = calday + yr*daysperyear
     end if

     if (masterproc) then
        if (single_column) then
           call shr_scam_getCloseLatLon(ncid_sst,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
           strt3(1) = closelonidx
           strt3(2) = closelatidx
           strt3(3) = 1
           cnt3(1)  = 1
           cnt3(2)  = 1
           cnt3(3)  = 1
        else
           strt3(1) = 1
           strt3(2) = 1
           strt3(3) = 1
           cnt3(1)  = lonsiz
           cnt3(2)  = latsiz
           cnt3(3)  = 1
        endif
     endif
!
! If model time is past current forward ice timeslice, read in the next
! timeslice for time interpolation.  Messy logic is for icecyc = .true. 
! interpolation between December and January (np1==1).  Note that 
! np1 is never 1 when icecyc is .false.
!
     if (caldayloc > cdayicep .and. .not. (np1==1 .and. caldayloc>cdayicem)) then
        if (icecyc) then
           np1 = mod(np1,12) + 1
        else
           np1 = np1 + 1
        end if
        if (np1>timesiz) then
           call endrun ('ICEINT: Attempt to read past end of ICE dataset')
        end if
        cdayicem = cdayicep
        call bnddyi(date_ice(np1), sec_ice(np1), cdayicep)

        if (.not.icecyc) then
           yr = date_ice(np1)/10000
           cdayicep = cdayicep + yr*daysperyear
        end if
        if (np1==1 .or. caldayloc<=cdayicep) then
           ntmp = nm
           nm = np
           np = ntmp
           if (masterproc) then
              strt3(3) = np1
              call wrap_get_vara_realx (ncid_sst,iceid,strt3,cnt3,xvar(1,1,np))
              write(iulog,*)'ICEINT: Read ice for date (yyyymmdd) ',date_ice(np1), &
                   ' sec ',sec_ice(np1)
           endif
           if (.not. prognostic_icesnow) then
              snwbdynh(np)=snwcnh(np1)
              snwbdysh(np)=snwcsh(np1)
           end if
           call scatter_field_to_chunk(1,1,1,plon,xvar(1,1,np),icebdy(1,begchunk,np))
        else
           if (masterproc) then
              write(iulog,*)'ICEINT: Input ice for date',date_ice(np1), &
                   ' sec ',sec_ice(np1), 'does not exceed model date',ncdate,&
                   ' sec ',ncsec,' Stopping.'
           endif
           call endrun
        end if
     end if
!
! Time interpolation.  Account for December-January interpolation if
! cycling ice dataset.  Again note that np1 is never 1 when icecyc is false
!
     if (np1==1) then                    ! Dec-Jan interpolation
        deltat = cdayicep + daysperyear - cdayicem
        if (caldayloc>cdayicep) then      ! We're in December
           fact1 = (cdayicep + daysperyear - caldayloc)/deltat
           fact2 = (caldayloc - cdayicem)/deltat
        else                                ! We're in January
           fact1 = (cdayicep - caldayloc)/deltat
           fact2 = (caldayloc + daysperyear - cdayicem)/deltat
        end if
     else
        deltat = cdayicep - cdayicem
        fact1 = (cdayicep - caldayloc)/deltat
        fact2 = (caldayloc - cdayicem)/deltat
     end if
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
     if (abs(fact1+fact2-1._dbl_kind)>1.e-6_dbl_kind .or. &
             fact1>1.000001_dbl_kind .or. fact1<-1.e-6_dbl_kind .or. &
             fact2>1.000001_dbl_kind .or. fact2<-1.e-6_dbl_kind) then 
        if (masterproc) then
           write(iulog,*)'ICEINT: Bad fact1 and/or fact2=',fact1,fact2
        endif
        call endrun
     end if

     do lchnk=begchunk,endchunk
        ncol = get_ncols_p(lchnk)
        call get_lat_all_p(lchnk, ncol, lats)
        do i=1,ncol
              aice(i,lchnk) = icebdy(i,lchnk,nm)*fact1 + icebdy(i,lchnk,np)*fact2
              aice(i,lchnk)=min(aice(i,lchnk),1._dbl_kind)
              aice(i,lchnk)=max(aice(i,lchnk),0._dbl_kind)

              if (frac(lchnk)%land(i) >= 1._dbl_kind) then
                aice(i,lchnk) = 0._dbl_kind
              end if
!
! If we are not prognosing snow then use a snow climatology.
! Compute snow height for given snow climatology phase shift snow
! height 6 months for southern hemisphere this only works when
! using ice climatology dataset
!
              if (.not. prognostic_icesnow) then
                 if (clat(lats(i)).gt.0) then
                    snowhice(i,lchnk)=snwbdynh(nm)*rhos/rhofresh * fact1 &
                         + snwbdynh(np)*rhos/rhofresh * fact2
                 else 
                    snowhice(i,lchnk)=snwbdysh(nm)*rhos/rhofresh*fact1 + &
                         snwbdysh(np)*rhos/rhofresh*fact2
                 end if
              end if
        end do
     end do

     return
   end subroutine iceint

#endif
end module ice_data
