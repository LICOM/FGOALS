module ice_comp
!-----------------------------------------------------------------------
!
! Module interface to initial and run methods for CAM-CSIM
! (CAM version of Community Sea-Ice Model)
!
!-----------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use shr_sys_mod,     only: shr_sys_flush, shr_sys_abort

  use esmf_mod,        only: ESMF_Clock
  use seq_timemgr_mod, only: seq_timemgr_EClockGetData
  use seq_infodata_mod,only: seq_infodata_type,seq_infodata_GetData
  use cam_logfile,     only: iulog

  use cam_control_mod, only: nsrest
  use commap,          only: clat, w
  use rgrid,           only: nlon
  use pmgrid,          only: plon, plat, plev   	
  use ppgrid,          only: pcols, begchunk,endchunk
  use phys_grid,       only: get_rlon_all_p, get_rlat_all_p, &
                             get_lon_all_p, get_lat_all_p, get_ncols_p, &
                             read_chunk_from_field, write_field_from_chunk, &
                             gather_chunk_to_field, scatter_field_to_chunk
  use physconst,       only: stebol
  use cam_history,     only: outfld
  use units,           only: getunit, freeunit
  use ioFileMod,       only: opnfil, getfil
  use infnan,          only: inf, uninit_r8
  use startup_initialconds, only: initial_file_get_id
  use wrap_nf

  use ice_time_manager,only: get_step_size, get_nstep, get_curr_calday, &
                             is_end_curr_day, is_first_step, &
                             timemgr_write_restart, timemgr_read_restart, &
                             timemgr_restart, timemgr_init
  use ice_constants,   only: ni, saltz, rLfs, Tffresh, rhofresh, Lfus, rhos, emissivity_ice, Lvap, &
                             init_constants, plevmx, tsnam, TfrezK
  use ice_types
  use ice_spmd         
  use ice_dh,          only: prognostic_icesnow, reset_csim_iceprops, energ
#ifndef COUP_SOM
  use ice_data,        only: iceini, iceint, icecyc
#endif
  use ice_filenames,   only: interpret_filename_spec
  use perf_mod

  implicit none
  private                  ! By default make everything private to this module
  save
!
! Public methods
!
  public ice_init          ! Initialization method
  public ice_run           ! Run method
  public ice_final         ! Finalization method
  public ice_write_restart ! Write restart method
  public ice_alloc         ! Allocation method
!
! Private data
!
  integer, parameter  :: nlen = 256                ! Length of character strings
  integer             :: nrf = -1                  ! Logical unit number for ice restart dataset
  integer             :: nrpf = -1                 ! logical unit number for ocn restart pointer file
  character(len=nlen) :: fname                     ! surface restart filename
  character(len=nlen) :: pname                     ! surface restart full pathname
  character(len=nlen) :: csim_branch_file = ' '    ! full pathname of restart file to branch from (nsrest=3)
  character(len=nlen) :: rest_pfile = './rpointer.ice' ! restart pointer file contains name of most recently
  character(len=nlen) :: rsfilename_spec_ice = '%c.camice.r.%y-%m-%d-%s' ! ice restarts
  character(len=nlen) :: focndomain                ! ocean domain file
  character(len=nlen) :: bndtvs                    ! sst file
  integer             :: ncid_dom                  ! netcdf file handle for ocn domain file

  type(frac_t), public, allocatable :: frac(:) 
  real(r8), allocatable :: landfrac_glob(:,:)

  real(r8), allocatable:: tair(:,:)
  real(r8), allocatable:: tssub(:,:,:)     ! Ice surface/subsurface temperatures

  real(r8), allocatable:: tsice(:,:)       ! Ice surface temperature
  real(r8), allocatable:: tsice_rad(:,:)   ! Equivalent LW_up temperature
  real(r8), allocatable:: previcefrac(:,:) ! previous fraction of total grid cell covered by seaice
  real(r8), allocatable:: Focn(:,:)        ! ocean-ice heat flux for basal and lateral melt (<0)
  real(r8), allocatable:: aice(:,:)        ! CSIM ice fraction

  ! the following are used by ice_srf - not by any ocean module
  real(r8), public, allocatable:: sicthk(:,:)      ! cam sea-ice thickness (m)
  real(r8), public, allocatable:: snowhice(:,:)    ! snow depth (liquid water) over ice

  ! ------Ice namelist variables -------
  ! prognostic_icesnow = .T,  
  !   prognostic snow over ice, currently limited to 0.5m.
  !   if this is false then a snow climatology is used (default .T.)
  ! logical reset_csim_iceprops = .F.,
  !   if true => resets the csim ice properties to base state
  !   No Snow Cover, TSICE and TS1-4 are all set to
  !   freezing. Default is false.
  !   The csim is sensitive to imbalances between the
  !   surface temperature and ice temperatures. When
  !   using an initial conditions dataset interpolated
  !   from a different resolution you may have to set this
  !   to true to get csim to run.  If set to true you will
  !   have to allow time for the ice to "spin-up".
  ! ice_conschk_frq
  !   energy conservation check frequency in CSIM4 code
  ! icecyc
  !   true => cycle ice fraction dataset
  !--------------------------------------
  integer :: ice_conschk_frq

!======================================================================= 
contains
!======================================================================= 

subroutine ice_init( ice_mpicom, ice_in, ice_out, EClock, infodata )

  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: 
  ! CAM sea ice surface fluxes initialization method
  !
  !----------------------------------------------------------------------- 
  !
  ! Arguments
  !
  integer, intent(in)      :: ice_mpicom
  type(ice_in_t) , pointer :: ice_in(:)
  type(ice_out_t), pointer :: ice_out(:)
  type(ESMF_Clock), intent(in) :: EClock  ! Synchronization clock
  type(seq_infodata_type),intent(in) :: infodata
  !
  ! Local variables
  !
  integer  :: ncol                                 ! Number of columns
  integer  :: i                                    ! Column loop index
  integer  :: c                                    ! Chunk loop index
  integer  :: k                                    ! Level index
  real(r8) :: tsice_tmp(pcols,begchunk:endchunk)   ! Surface ice temp at start
  integer  :: unitn, ntspdy                        ! input/output unit
  integer  :: start_ymd                            ! Start date (YYYYMMDD)
  integer  :: start_tod                            ! Start time of day (sec)
  integer  :: ref_ymd                              ! Reference date (YYYYMMDD)
  integer  :: ref_tod                              ! Reference time of day (sec)
  integer  :: stop_ymd                             ! Stop date (YYYYMMDD)
  integer  :: stop_tod                             ! Stop time of day (sec)
  integer  :: dtime                                ! Time-step
  logical  :: log_print                            ! Flag to print out log information or not
  logical  :: perpetual_run                        ! If in perpetual mode or not
  integer  :: perpetual_ymd                        ! Perpetual date (YYYYMMDD)
  integer  :: ierr                                 ! Error code
  character(len=256) :: calendar                   ! Calendar type
  !
#ifdef COUP_SOM
  logical  :: icecyc                               ! needed for code to compile in SOM mode, but not used
#endif
   	
  namelist /csim_inparm/ prognostic_icesnow, reset_csim_iceprops, &
	                 ice_conschk_frq, icecyc, rest_pfile, csim_branch_file, &
	                 bndtvs, focndomain
  !----------------------------------------------------------------------- 
  !
  ! Allocate dynamic memory
  ! 
  call ice_types_alloc( ice_in, ice_out )
  call ice_alloc( )
  !
  ! Initialize ocn MPI communicator 
  !
  call ice_spmd_init( ice_mpicom )
  !
  ! Read in ice namelist
  !	
  icecyc              = .true.    ! false => do not cycle ice/sst dataset (assume multiyear)
  prognostic_icesnow  = .true.    ! snow falls on ice by default but it is limited to 0.5 meter.
  reset_csim_iceprops = .false.   ! use initial condition info unless need to reset ice properties in csim
  ice_conschk_frq     = 0
  
  if (masterproc) then
     unitn = getunit()
     write(iulog,*) 'Read in cam-csim namelist from ice_in'
     open( unitn, file='ice_in', status='old' )
     ierr = 1
     do while ( ierr /= 0 )
        read(unitn, csim_inparm, iostat=ierr)
        if (ierr < 0) then
           call shr_sys_abort( 'ice_comp encountered end-of-file on namelist read' )
        endif
     end do
     call freeunit( unitn )
  end if
  !
  ! Initialize time manager.
  !	
  call seq_timemgr_EClockGetData( EClock, start_ymd=start_ymd,       &
       start_tod=start_tod, ref_ymd=ref_ymd, &
       ref_tod=ref_tod, stop_ymd=stop_ymd,   &
       stop_tod=stop_tod, dtime=dtime,       &
       calendar=calendar )

  call seq_infodata_GetData(infodata, &
       perpetual=perpetual_run, perpetual_ymd=perpetual_ymd)

  if ( nsrest == 0 )then
     call timemgr_init( calendar_in=calendar, start_ymd=start_ymd, &
          start_tod=start_tod, ref_ymd=ref_ymd,      &
          ref_tod=ref_tod, stop_ymd=stop_ymd,        &
          stop_tod=stop_tod, dtime_in=dtime,         &
          perpetual_run=perpetual_run,               &
          perpetual_ymd=perpetual_ymd )
  end if
  !
  ! Read initial/restart data (including time manager restart info)
  !
  if (nsrest == 0) then
     call ice_read_inidat( ice_out )
  else
     call ice_read_restart( ice_out, stop_ymd, stop_tod )
  end if
  !
  ! Determine consistency checks with namelist
  !
  if (masterproc) then
     dtime = get_step_size()
     ntspdy = nint(86400._r8/dtime) 
     if (ice_conschk_frq < 0) then
        ice_conschk_frq = -ice_conschk_frq*ntspdy
     end if
     if (ice_conschk_frq > 0) then
        write(iulog,*)'ICE global energy checking will be done every ',ice_conschk_frq,' timesteps'
     end if

     if ( reset_csim_iceprops) then
        write(iulog,*)'CSIM ICE properties being reset to a new base state'
     end if
     if (prognostic_icesnow) then
        write(iulog,*)'Snow will accumulate to a maximum over sea-ice'
     else
        write(iulog,*)'Snow over sea-ice will be set to a climatology'
     end if

     if (icecyc) then
        write(iulog,*)'ICE dataset will be reused for each model year'
     else
        write(iulog,*)'ICE dataset will not be cycled'
     end if
  end if
  call mpi_bcast(icecyc             ,1,MPI_LOGICAL,0,ice_mpicom,ierr)
  call mpi_bcast(reset_csim_iceprops,1,MPI_LOGICAL,0,ice_mpicom,ierr)
  call mpi_bcast(prognostic_icesnow ,1,MPI_LOGICAL,0,ice_mpicom,ierr)
  call mpi_bcast(ice_conschk_frq    ,1,MPI_INTEGER,0,ice_mpicom,ierr)
  
#if ( defined COUP_SOM )

  !--------------------------------------------------------------------------------
  ! Sea-ice model for SOM
  !--------------------------------------------------------------------------------
  !
  ! Initialize freeze/melt potential and csim constants
  !
  if (is_first_step()) then
     do c=begchunk,endchunk
        ncol = get_ncols_p(c)
        do i=1,ncol
           Focn(i,c) = 0._r8
        end do
     end do
  end if
  call init_constants ()

#else

  !--------------------------------------------------------------------------------
  ! Sea-ice model for DOM
  !--------------------------------------------------------------------------------
  !
  ! on a restart the ice values come from the restart data sets
  !
  call iceini(bndtvs)
  if ( is_first_step() ) then
     do c = begchunk,endchunk
        ncol = get_ncols_p(c)
        tsice_tmp(:ncol,c) = tsice(:ncol,c)
     end do
     call iceint( aice, snowhice, frac, prev_timestep=.true.)
     do c = begchunk,endchunk
        ncol = get_ncols_p(c)
        do i = 1, ncol
           ice_out(c)%aice(i)     = aice(i,c)
           ice_out(c)%areafrac(i) = ice_fraction( aice(i,c),frac(c)%land(i) )
        end do
        tsice(:ncol,c) = tsice_tmp(:ncol,c)
     end do
     
  end if
  do c=begchunk,endchunk
     ncol = get_ncols_p(c)
     do i=1,ncol
        previcefrac(i,c) = ice_out(c)%areafrac(i) 
     end do
  end do
  !
  ! Set constant sea ice thickness for each hemisphere. For new ice and
  ! for the reset_csim_ice_props=.T. namelist variable set
  ! initial temperature profile to freezing
  !
  ! if reset_csim_iceprops set the previous ice fraction to 0 to make
  ! all ice appear as though it is new sea ice. This will effectively
  ! reset all ice properties to a the base state of new ice.  This is
  ! only provided on the first timestep of a run.
  !
  if ( is_first_step().and. reset_csim_iceprops ) then
     write(iulog,*)'**********Resetting previcefrac to 0********'
     previcefrac(:,:)=0._r8
  end if
  if ( is_first_step() ) then
     do c=begchunk,endchunk
        call newiceproperties ( c, ice_out(c), frac(c)%land, ice_in(c) )
     end do
  end if

#endif
  !
  ! Initialize sea-ice albedos at NSTEP = 0. 
  !
  if (is_first_step()) then
     do c = begchunk,endchunk
        ncol = get_ncols_p(c)
        call albice(c,ncol, &
             tair(1,c), snowhice(1,c), &
             ice_out(c)%asdir, ice_out(c)%aldir, &
             ice_out(c)%asdif, ice_out(c)%aldif, &
             ice_out(c)%areafrac, sicthk(1,c) )
        !
        ! Set TS for ice
        !
        do i = 1, ncol
           if ( ice_out(c)%areafrac(i) > 0._r8 )then
#ifdef COUP_SOM
              ice_out(c)%ts(i)   = tsice(i,c)
#else
              ice_out(c)%ts(i)   = tsice_rad(i,c)
#endif
              ice_out(c)%lwup(i) = stebol*ice_out(c)%ts(i)**4
           end if
        end do
     end do
  end if

end subroutine ice_init

!
!----------------------------------------------------------------------- 
!

subroutine ice_run ( ice_in, ice_out )

  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: 
  ! CAM sea ice surface fluxes run method
  !
  !----------------------------------------------------------------------- 
  !
  ! Input/Output arguments
  !
  type(ice_in_t) , intent(in)    :: ice_in(begchunk:endchunk)
  type(ice_out_t), intent(inout) :: ice_out(begchunk:endchunk)

  !---------------------------Local variables-----------------------------
  integer :: dtime            ! timestep size [seconds]
  real(r8):: rtime            ! calendar day for next timestep

  integer :: ncol           ! number of columns in chunk
  integer :: c              ! chunk index
  integer :: i              ! temporary variables
  integer :: k              ! ice depth index

  real(r8) :: snowfall(pcols,begchunk:endchunk)  ! total snowfall rate
  real(r8) :: fsns(pcols,begchunk:endchunk)      ! SW absorbed in ice
  real(r8) :: evap(pcols,begchunk:endchunk)      ! evaporative flux off snow and ice
  real(r8) :: aiceinit(pcols,begchunk:endchunk)
  real(r8) :: gmsie_in          ! global mean sea ice/snow internal energy on input
  real(r8) :: gmsie_out         ! global mean sea ice/snow internal energy on output
  real(r8) :: sier              ! sea ice energy rate
  real(r8) :: fluxsum           ! F_ice - F_ocn - F_frzmlt
  real(r8) :: F_ice
  real(r8) :: F_ocn
  real(r8) :: F_frzmlt
  real(r8) :: errterm           ! error term in global energy calculations
  real(r8) :: deltae(pcols,begchunk:endchunk)    ! change in energy
  real(r8) :: deltaaice(pcols,begchunk:endchunk) ! change in aice
  real(r8) :: tmp(pcols,begchunk:endchunk)       ! temporary field

  logical :: ice_conschk = .true.        ! whether to apply global energy conservation check

  !-----------------------------------------------------------------------
  !
  ! Exit from run method if sea-ice not in this grid
  !

  !
  ! Initial condition file output, get ice-area from dataset
  !
   call t_startf ('ice_write_inidat')
   call ice_write_inidat( ice_out )
   call t_stopf ('ice_write_inidat')
   !
#ifndef COUP_SOM
   call t_startf ('iceint')
   call iceint ( aice, snowhice, frac )
   call t_stopf ('iceint')
#endif
  !
  ! Time step 
  ! 
  call t_startf ('ice_st')
  dtime = get_step_size()
  rtime = dtime
  !
  ! set up snowfall here so it doesn't have to be private in the omp call
  !
  do c = begchunk,endchunk
     ncol = get_ncols_p(c)
     do i = 1,ncol
#ifndef COUP_SOM
        !
        ! Update previous aice fraction
        !
        previcefrac(i,c) = ice_out(c)%areafrac(i)
        ice_out(c)%areafrac(i) = ice_fraction( aice(i,c),frac(c)%land(i) )
#endif
	if (prognostic_icesnow) then
           snowfall(i,c) = ice_in(c)%snow(i)
        else
           snowfall(i,c) = 0._r8	
        end if
     end do
  end do
  call t_stopf ('ice_st')
  !
  ! check for new ice and give it some default values for sicthk and snowhice
  ! and tssub
  !
#ifndef COUP_SOM
  call t_startf('new_ice')
  do c = begchunk,endchunk
     call newiceproperties ( c, ice_out(c), frac(c)%land, ice_in(c) )
  end do
  call t_stopf('new_ice')
#endif
  !
  ! BPB Compute initial sea ice/snow internal energy and globally average over
  ! BPB ocean.  gmsie_in = global mean sea ice/snow internal energy on input
  !
#ifdef COUP_SOM
  call t_startf('sea_ice_energy')
  ice_conschk = .false.
  if (ice_conschk_frq > 0) then
     ice_conschk = mod (get_nstep(), ice_conschk_frq) == 0
  end if
  if (ice_conschk) then
     gmsie_in = sea_ice_energy (ice_in, dtime, 1, deltae, deltaaice, snowfall)
  end if
  call t_stopf('sea_ice_energy')
#endif

  call t_startf('sea_ice')
!$OMP PARALLEL DO PRIVATE (C, NCOL, I)
  do c = begchunk,endchunk
     ncol = get_ncols_p(c)

     ! Sea ice surface fluxes and temperatures

     call seaice (c, ncol, rtime, aice(1,c), tsice(1,c), &
                  sicthk(1,c), snowhice(1,c), ice_in(c)%ubot, &
                     ice_in(c)%vbot, ice_in(c)%tbot, &
                  ice_in(c)%qbot, ice_in(c)%thbot, ice_in(c)%zbot, &
                     ice_in(c)%pbot ,ice_in(c)%flwds, &
                  ice_in(c)%sols, ice_in(c)%soll, ice_in(c)%solsd, &
                     ice_in(c)%solld, ice_out(c)%asdir, &
                  ice_out(c)%aldir, ice_out(c)%asdif, ice_out(c)%aldif, &
    	             snowfall(1,c), &
                  ice_in(c)%tsocn, ice_in(c)%frzmlt, Focn(1,c),  &
                  tssub(1,1,c), ice_out(c)%cflx, ice_out(c)%wsx, ice_out(c)%wsy, &
                     ice_out(c)%ts, ice_out(c)%shf, &
         	  ice_out(c)%lhf, ice_out(c)%lwup, ice_out(c)%tref, fsns(1,c), evap(1,c), &
                  aiceinit(1,c), ice_in )

     do i = 1,ncol
        previcefrac(i,c)       = ice_out(c)%areafrac(i)
        ice_out(c)%areafrac(i) = ice_fraction( aice(i,c),frac(c)%land(i) )
	ice_out(c)%sicthk(i)   = sicthk(i,c)
	ice_out(c)%focn(i)     = Focn(i,c)
	ice_out(c)%aice(i)     = aice(i,c)
        ice_out(c)%fswabs(i)   = fsns(i,c)
     end do
     !
     ! save off tair for restart
     !
     do i = 1,ncol
        tair(i,c) = ice_in(c)%tbot(i) 
     end do
     !
     !    Calculate Tsice_rad
     !
     do i=1,ncol
        if(ice_out(c)%areafrac(i) > 0._r8) then
           tsice_rad(i,c) = sqrt(sqrt(-ice_out(c)%lwup(i)/stebol))
        else
           tsice_rad(i,c) = TfrezK
        endif
     end do
     call ice_output( c )
  end do
  !
  ! Albedos 
  ! Note the total albedo here that is returned to the atmosphere 
  ! model is based on a weighted sum of the albedo over ice and ocean
  ! using fractional areas from this time step. The absorbed shortwave over
  ! sea ice in the next step uses ice albedos that are saved at there
  ! present value but with a NEW fractional area that is input prior to 
  ! the next time through the sea ice model.  Hence
  ! there is a time step mismatch in the absorbed solar over sea ice. 
  ! CCSM would not allow such a thing, but here we are specifying sst, 
  ! over the ocean fraction anyway so it doesn't really matter. 
  ! ---save off ice albedos for sea ice routine per email Bitz.
  ! I should note that I made one change to the "physics" from John's
  ! original fracice implementation. John had the absorbed solar by the
  ! sea ice equal to the gridcell average.  This is pretty far off when
  ! the sea ice fraction is small. I realize that it is standard practise
  ! in many models, but it doesn't have to be.  Therefore I have compute a
  ! special srfrad over ice and I send the ice albedos to the restart file.
  !
  ! 
!$OMP PARALLEL DO PRIVATE (C, NCOL, I)
  do c = begchunk,endchunk
     ncol = get_ncols_p(c)
     call albice(c,ncol, &
                 ice_in(c)%tbot,snowhice(1,c), &
                 ice_out(c)%asdir, ice_out(c)%aldir, &
                 ice_out(c)%asdif, ice_out(c)%aldif, &
                 ice_out(c)%areafrac, sicthk(1,c) )
  end do

  call t_stopf('sea_ice')

#ifdef COUP_SOM

  if (ice_conschk) then
     call t_startf('ice_conschk')

     !BPB Check sea ice energy conservation: sier = sea ice energy rate

     gmsie_out = sea_ice_energy (ice_in, dtime, 2, deltae, deltaaice, snowfall)

     sier = (gmsie_out - gmsie_in) / dtime

     !BPB Compute global mean F_ice, F_frzmlt and F_ocn
     !BPB Then compute output sea ice/snow internal energy and globally
     !BPB average over ocean.  gmsie_out = global mean sea ice/snow
     !BPB internal energy on output

     call gmean_ice (ice_out, ice_in, fsns, snowfall, F_ice, &
                     F_ocn, F_frzmlt, dtime, deltae, evap, aiceinit)

     fluxsum = F_ice - F_ocn - F_frzmlt
     errterm = abs (sier/fluxsum - 1._r8)
     if (masterproc) then
        write(iulog,'(a,1p,4e12.4)') 'ICE: sier,fluxsum,diff,errterm=', &
                                 sier, fluxsum, fluxsum-sier, errterm
                   
        write(iulog,'(a,1p,5e12.4)') '        gmsie_in,gmsie_out,F_ice,F_ocn,F_frzmlt=', &
                                 gmsie_in, gmsie_out, F_ice, F_ocn, F_frzmlt
     end if
     call t_stopf('ice_conschk')
  end if

#endif

   if (is_end_curr_day ()) then
      do c = begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = ice_out(c)%areafrac(i) 
         end do
      end do
      call print_coverage ('icefrac', ' million km^2', tmp, 1.e-12_r8)

      do c = begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = ice_out(c)%areafrac(i)*sicthk(i,c)
         end do
      end do
      call print_coverage ('icevol ', ' 10^13m^3', tmp, 1.e-13_r8)

      do c = begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            tmp(i,c) = ice_out(c)%areafrac(i)*snowhice(i,c)
         end do
      end do
      call print_coverage ('snowvol', ' 10^13m^3', tmp, 1.e-13_r8)
   end if

end subroutine ice_run

subroutine ice_final( ice_in, ice_out )

  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: 
  ! CAM sea ice finalization
  !
  !----------------------------------------------------------------------- 
  
  type(ice_in_t) , pointer :: ice_in(:)
  type(ice_out_t), pointer :: ice_out(:)

  deallocate (ice_in)
  deallocate (ice_out)
  deallocate (sicthk)
  deallocate (snowhice)
  deallocate (previcefrac)
  deallocate (tsice_rad)
  deallocate (tsice)
  deallocate (Focn)
  deallocate (aice)
  deallocate (tssub)

end subroutine ice_final
!
!----------------------------------------------------------------------- 
!

subroutine ice_output( c )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Output Ice data to history files.
!
!----------------------------------------------------------------------- 
!
! Input/Output arguments
!
   integer, intent(in) :: c       ! Chunk index

!---------------------------Local variables-----------------------------
   integer :: k      ! Layer index

   !
   ! history Outfield calls
   !
   do k=1,plevmx
      call outfld(tsnam(k), tssub(1,k,c), pcols, c)
   end do
   call outfld('SICTHK',   sicthk   (1,c),  pcols, c)
   call outfld('TSICE',    tsice    (1,c),  pcols, c)
   call outfld('SNOWHICE', snowhice (1,c) , pcols, c)

end subroutine ice_output

!
!----------------------------------------------------------------------- 
!

subroutine ice_write_inidat( ice_out )
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: 
  ! Output Ice initial condition data to history files.
  !
  !----------------------------------------------------------------------- 
  type(ice_out_t), intent(in) :: ice_out(begchunk:endchunk)
  !---------------------------Local variables-----------------------------
  integer :: k       ! Layer index
  integer :: c       ! Chunk index
  !----------------------------------------------------------------------- 
  !
  do c = begchunk, endchunk
     do k=1,plevmx
        call outfld(trim(tsnam(k))//'&IC', tssub(1,k,c), pcols, c)
     end do
     call outfld('TSICERAD&IC'     , tsice_rad(1,  c)   , pcols, c)
     call outfld('TSICE&IC   '     , tsice    (1,  c)   , pcols, c)
     call outfld('SNOWHICE&IC'     , snowhice (1,  c)   , pcols, c)
     call outfld('ICEFRAC&IC '     , ice_out(c)%areafrac, pcols, c)
     call outfld('SICTHK&IC  '     , sicthk   (1,  c)   , pcols, c)
  end do
  
end subroutine ice_write_inidat

!
!----------------------------------------------------------------------- 
!

subroutine ice_read_inidat( ice_out)
  !-----------------------------------------------------------------------
  !
  ! Purpose: Read in CAM-CSIM sea-ice variables from initial datafile.
  !
  !-----------------------------------------------------------------------
  use ncdio_atm
  use scamMod,         only: scmlon,scmlat,single_column
  use netcdf, only : nf90_inq_varid, nf90_get_var
  use pio, only : file_desc_t
  use shr_scam_mod,    only: shr_scam_GetCloseLatLon
  ! 
  ! Arguments
  !
  type(ice_out_t), intent(inout) :: ice_out(begchunk:endchunk)
  ! 
  ! Local variables
  !
  logical :: readvar                ! inquiry:  true => variable exists on file
  integer :: m, c                   ! Indices
  integer :: ncols                  ! number of columns
  integer :: ncid_dom               ! netcdf file id
  type(file_desc_t), pointer :: ncid_ini               ! netcdf file id
  real(r8), pointer :: arr2d(:,:)   ! temporary 2D array
  character(len=256) :: locfn       ! netcdf local filename to open
  character(len=16)  :: fieldname   ! field name
  real(r8), allocatable :: t3_tmp(:,:,:)
  real(r8), allocatable :: t2_tmp(:,:)
  integer :: closelatidx
  integer :: closelonidx
  integer :: fracid,rcode
  real(r8) :: closelat,closelon
  !-----------------------------------------------------------------------
  allocate ( arr2d(1:pcols,begchunk:endchunk) )
  !
  ! read in land fraction (compute global land fraction on master)
  !
  if (masterproc) then
     call getfil(focndomain, locfn)
     call wrap_open(locfn, 0, ncid_dom)
  endif
  if (single_column) then
     call shr_scam_getCloseLatLon(ncid_dom,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
     rcode = nf90_inq_varid(ncid_dom, 'frac', fracid)
     rcode = nf90_get_var(ncid_dom,fracid,arr2d,start=(/closelonidx,closelatidx/),count=(/1,1/))
  else
     fieldname = 'frac'
     call read_domain (fieldname, ncid_dom, 'ni', 'nj', 'xc','yc',&
          1, pcols, begchunk, endchunk, arr2d, readvar)
     if(.not. readvar) call shr_sys_abort('som: error in reading frac')
  end if
  if (masterproc) then
     call wrap_close(ncid_dom)
  end if
  do c = begchunk, endchunk
     ncols = get_ncols_p(c)
     ! first convert from ocn fraction to land fraction
     arr2d(:ncols,c) = 1._r8 - arr2d(:ncols,c)	
     frac(c)%land(:ncols) = arr2d(:ncols,c)
  end do
  call gather_chunk_to_field(1, 1, 1, plon, arr2d, landfrac_glob)
  !
  ! Read initial data
  !
  ncid_ini => initial_file_get_id()
  fieldname = 'TSICE'
  call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
       tsice   , readvar, grid_map='phys')
  if(.not. readvar) call shr_sys_abort('ice_comp.F90')

  fieldname = 'SNOWHICE'
  call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
       snowhice, readvar, grid_map='phys')
  if(.not. readvar) call shr_sys_abort('ice_comp.F90')

  ! read TS1, TS2, TS3, TS4 in "plevmx" loop
  do m = 1,plevmx
     fieldname = tsnam(m)
     call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
          arr2d, readvar, grid_map='phys')
     if(.not. readvar) call shr_sys_abort('ice_comp.F90')
     do c = begchunk,endchunk
        tssub(:,m,c) = arr2d(:,c)
     end do
  end do

#if ( defined COUP_SOM )
  fieldname = 'SICTHK'
  call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
       sicthk  , readvar, grid_map='phys')
  if(.not. readvar) call shr_sys_abort('ice_comp.F90')
  do c = begchunk,endchunk
     ncols = get_ncols_p(c)
     ice_out(c)%sicthk(:ncols) = sicthk(:ncols,c)
  end do

  fieldname = 'ICEFRAC'
  call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
       arr2d, readvar, grid_map='phys')
  do c = begchunk,endchunk
     ncols = get_ncols_p(c)
     ice_out(c)%areafrac(:ncols) = arr2d(:ncols,c)
  end do
  if(.not. readvar) call shr_sys_abort('ice_comp.F90')

  do c = begchunk,endchunk
     ncols = get_ncols_p(c)
     where ((ice_out(c)%areafrac(:ncols) + frac(c)%land(:ncols)) > 1.0_r8)
        ice_out(c)%areafrac(:ncols) = 1._r8 - frac(c)%land(:ncols)
     end where
     where (frac(c)%land(:ncols) < 1._r8)
        aice(:ncols,c) = ice_out(c)%areafrac(:ncols) /(1._r8 - frac(c)%land(:ncols))
     elsewhere
        aice(:ncols,c) = 0._r8
     end where
     ice_out(c)%aice(:ncols) = aice(:ncols,c)
  end do
#endif

  ! --------------------------------------------------------------------
  ! Read optional ice fields (if not found, will be arbitrarily initialized)
  ! --------------------------------------------------------------------

  fieldname = 'TSICERAD'
  call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
       tsice_rad, readvar, grid_map='phys')
  if(.not. readvar) then
     do c = begchunk, endchunk
        ncols = get_ncols_p(c)
        tsice_rad(:ncols,c) = tsice(:ncols,c)
     end do
     if (masterproc) write(iulog,*) trim(fieldname), ' initialized with TSICE.'
  end if

  ! --------------------------------------------------------------------
  ! Read bottom level atmosphere temperature
  ! --------------------------------------------------------------------

  fieldname = 'TBOT'
  call infld(fieldname, ncid_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
       arr2d, readvar, grid_map='phys')


  if(readvar) then
     do c = begchunk,endchunk
        tair(:,c) = arr2d(:,c)
     end do
     deallocate ( arr2d )
  else	
     deallocate ( arr2d )
     if (masterproc) then
        write(iulog,*) trim(fieldname), ' initialized with lowest level of T'
     end if

     fieldname = 'T'
     allocate ( t3_tmp(pcols,plev,begchunk:endchunk) )
     call infld(fieldname, ncid_ini, 'lon', 'lev', 'lat', 1, pcols, 1, plev, &
	begchunk, endchunk, t3_tmp, readvar, grid_map='phys')
     if(.not. readvar) call shr_sys_abort('ice_comp.F90')
     do c = begchunk,endchunk
	ncols =  get_ncols_p(c)
        tair(:ncols,c) = t3_tmp(:ncols,plev,c)
     end do

     deallocate ( t3_tmp )
     
  end if



end subroutine ice_read_inidat
!
!----------------------------------------------------------------------- 
!
subroutine ice_alloc( )

  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: 
  ! CAM sea ice module data allocation
  !
  !----------------------------------------------------------------------- 
  
  !---------------------------------------------
  integer :: c	
  !---------------------------------------------

  allocate (tair(pcols,begchunk:endchunk))
  allocate (frac(begchunk:endchunk))
  allocate (landfrac_glob(plon, plat))

  allocate (sicthk(pcols,begchunk:endchunk))
  allocate (previcefrac(pcols,begchunk:endchunk))
  allocate (aice(pcols,begchunk:endchunk))
  allocate (snowhice(pcols,begchunk:endchunk))
  allocate (tsice_rad(pcols,begchunk:endchunk))
  allocate (tsice(pcols,begchunk:endchunk))
  allocate (Focn(pcols,begchunk:endchunk))
  allocate (tssub(pcols,plevmx,begchunk:endchunk))

  sicthk   (:,:)     = inf
  snowhice (:,:)     = 0._r8
  previcefrac  (:,:) = uninit_r8
  tssub(:,:,:)       = 0._r8
  tsice_rad (:,:)    = inf
  tsice(:,:)         = inf
  Focn(:,:)          = 0._r8
  aice (:,:)         = inf

  do c = begchunk,endchunk
     frac(c)%land(:) = 0._r8
  end do
  landfrac_glob(:,:) = 0._r8

end subroutine ice_alloc

!
!----------------------------------------------------------------------- 
!

subroutine ice_write_restart( ice_out, &
                              yr_spec, mon_spec, day_spec, sec_spec )

  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: 
  ! Write out sea-ice restart data to given file
  !
  !----------------------------------------------------------------------- 
  !
  ! Input arguments
  !
  type(ice_out_t), intent(in) :: ice_out(begchunk:endchunk)
  integer        , intent(in) :: yr_spec         ! Simulation year
  integer        , intent(in) :: mon_spec        ! Simulation month
  integer        , intent(in) :: day_spec        ! Simulation day
  integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
  !
  ! Local variables
  !
  integer  :: i                                           ! Loop index
  integer  :: ncol                                        ! Number of columns
  real(r8) :: tmpfield3d(pcols,plevmx,begchunk:endchunk)  ! Temporary field
  real(r8) :: tmpfield2d(pcols,begchunk:endchunk)         ! Temporary field
  integer  :: ioerr                                       ! Error status for IO
  !-----------------------------------------------------------------------

  ! Determine and open surface restart dataset

  if (masterproc) then
     fname = interpret_filename_spec( rsfilename_spec_ice, &
          yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
     if ( nrf == -1 ) nrf = getunit()
     call opnfil(fname, nrf, 'u')
  endif
  
  do i=begchunk,endchunk
     tmpfield2d(:,i) = ice_out(i)%asdir(:)
  end do
  call write_field_from_chunk(nrf,1,1,1,tmpfield2d)

  do i=begchunk,endchunk
     tmpfield2d(:,i) = ice_out(i)%asdif(:)
  end do
  call write_field_from_chunk(nrf,1,1,1,tmpfield2d)

  do i=begchunk,endchunk
     tmpfield2d(:,i) = ice_out(i)%aldir(:)
  end do
  call write_field_from_chunk(nrf,1,1,1,tmpfield2d)

  do i=begchunk,endchunk
     tmpfield2d(:,i) = ice_out(i)%aldif(:)
  end do
  call write_field_from_chunk(nrf,1,1,1,tmpfield2d)

  call write_field_from_chunk(nrf,1,1,1,tsice)
  call write_field_from_chunk(nrf,1,1,1,tair)
  call write_field_from_chunk(nrf,1,1,1,aice)
  call write_field_from_chunk(nrf,1,1,1,sicthk)
  call write_field_from_chunk(nrf,1,1,1,snowhice)
  call write_field_from_chunk(nrf,1,1,1,Focn)

  do i=begchunk,endchunk
     ncol = get_ncols_p(i)
     tmpfield2d(:ncol,i)   = ice_out(i)%areafrac(:ncol)
     tmpfield3d(:ncol,:,i) = tssub(:ncol,:,i)
  end do
  call write_field_from_chunk(nrf,1,1,1,tmpfield2d)
  call write_field_from_chunk(nrf,1,plevmx,1,tmpfield3d)
  
  do i=begchunk,endchunk
     tmpfield2d(:,i) = ice_out(i)%lwup(:)
  end do
  call write_field_from_chunk(nrf,1,1,1,tmpfield2d)

  do i=begchunk,endchunk
     tmpfield2d(:,i) = frac(i)%land(:)
  end do
  call write_field_from_chunk(nrf,1,1,1,tmpfield2d)

  ! write time manager restart
  
  if (masterproc) then
     call timemgr_write_restart(nrf)
  end if
  
  if (masterproc) then
     close (nrf)

     if ( nrpf == -1 ) nrpf = getunit()
     call opnfil(rest_pfile, nrpf, 'f')
     rewind nrpf
     write (nrpf,'(a)') trim(fname)
     close(nrpf)
     write(iulog,*)'(ice_write_restart): successfully wrote local restart pointer file ',trim(rest_pfile)
  end if

end subroutine ice_write_restart

!
!----------------------------------------------------------------------- 
!

subroutine ice_read_restart( ice_out, stop_ymd, stop_tod )

  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: 
  ! Read in sea-ice restart data to given file
  !
  !----------------------------------------------------------------------- 
  !
  ! Input arguments
  !
  type(ice_out_t), intent(inout) :: ice_out(begchunk:endchunk)
  integer        , intent(in)    :: stop_ymd      ! Stop date (YYYYMMDD)
  integer        , intent(in)    :: stop_tod      ! Stop time of day (sec)
  !
  ! Local variables
  !
  integer  :: i                                          ! Loop index
  integer  :: ncol                                       ! Number of columns
  real(r8) :: tmpfield3d(pcols,plevmx,begchunk:endchunk) ! Temporary field
  real(r8) :: tmpfield2d(pcols,begchunk:endchunk)        ! Temporary field
  integer  :: ioerr                                      ! Error status for IO
  !-----------------------------------------------------------------------

  ! Determine and open surface restart dataset
  ! restart run =>nsrest=1) and branch run=>nsrest=3.  
  ! Only read the restart pointer file for a restart run.
  
  if (masterproc) then
     if (nsrest == 1) then
        nrpf = getunit()
        call opnfil (rest_pfile, nrpf, 'f', status="old")
        read (nrpf,'(a)') pname
     else
        if ( trim(csim_branch_file) == '' )then
           call shr_sys_abort ('ocn_read_restart: csim_branch_file is empty')
        end if
        if ( csim_branch_file(1:1) /= '/' )then
           call shr_sys_abort ('ocn_read_restart: csim_branch_file is not an absolute pathname')
        end if
        if ( len_trim(csim_branch_file) > nlen )then
           call shr_sys_abort ('ocn_read_restart: csim_branch_file is too long :'//csim_branch_file)
        end if
        pname = trim(csim_branch_file)
     end if
     call getfil(pname, fname)
     nrf = getunit()
     call opnfil(fname, nrf, 'u')
  endif
  
  call read_chunk_from_field(nrf,1,1,1,tmpfield2d)
  do i = begchunk,endchunk
     ncol = get_ncols_p(i)
     ice_out(i)%asdir(:ncol) = tmpfield2d(:ncol,i)
  end do
  
  call read_chunk_from_field(nrf,1,1,1,tmpfield2d)
  do i = begchunk,endchunk
     ncol = get_ncols_p(i)
     ice_out(i)%asdif(:ncol) = tmpfield2d(:ncol,i)
  end do
  
  call read_chunk_from_field(nrf,1,1,1,tmpfield2d)
  do i = begchunk,endchunk
     ncol = get_ncols_p(i)
     ice_out(i)%aldir(:ncol) = tmpfield2d(:ncol,i)
  end do
  
  call read_chunk_from_field(nrf,1,1,1,tmpfield2d)
  do i = begchunk,endchunk
     ncol = get_ncols_p(i)
     ice_out(i)%aldif(:ncol) = tmpfield2d(:ncol,i)
  end do
  
  call read_chunk_from_field(nrf,1,1,1,tsice)
  call read_chunk_from_field(nrf,1,1,1,tair)
  call read_chunk_from_field(nrf,1,1,1,aice)
  call read_chunk_from_field(nrf,1,1,1,sicthk)
  call read_chunk_from_field(nrf,1,1,1,snowhice)
  call read_chunk_from_field(nrf,1,1,1,Focn)
  
  do i=begchunk,endchunk
     ncol = get_ncols_p(i)
     ice_out(i)%sicthk(:ncol) = sicthk(:ncol,i)
     ice_out(i)%focn(:ncol)   = Focn(:ncol,i)
     ice_out(i)%aice(:ncol)   = aice(:ncol,i)
  end do

  call read_chunk_from_field(nrf,1,1,1,tmpfield2d)
  call read_chunk_from_field(nrf,1,plevmx,1,tmpfield3d)
  do i=begchunk,endchunk
     ncol = get_ncols_p(i)
     ice_out(i)%areafrac(:ncol) = tmpfield2d(:ncol,i)
     tssub(:ncol,:,i)           = tmpfield3d(:ncol,:,i)
  end do
  
  call read_chunk_from_field(nrf,1,1,1,tmpfield2d)
  do i = begchunk,endchunk
     ncol = get_ncols_p(i)
     ice_out(i)%lwup(:ncol) = tmpfield2d(:ncol,i)
  end do
  
  call read_chunk_from_field(nrf,1,1,1,tmpfield2d)
  do i = begchunk,endchunk
     ncol = get_ncols_p(i)
     frac(i)%land(:ncol) = tmpfield2d(:ncol,i)
  end do
  call gather_chunk_to_field(1, 1, 1, plon, tmpfield2d, landfrac_glob)
 
  ! Restart the time manager.
  
  if (masterproc) then
     call timemgr_read_restart(nrf)
  end if
  call timemgr_restart( stop_ymd=stop_ymd, stop_tod=stop_tod )
  
  if (masterproc) then
     close(nrf)
  end if
  
end subroutine ice_read_restart

!
!----------------------------------------------------------------------- 
!
#ifndef COUP_SOM

subroutine newiceproperties( lchnk, ice_out, landfrac, ice_in )
!----------------------------------------------------------------------- 
! 
! Purpose: 
!	Set the initial sea-ice surface variables associated with ICE.
!	This routine is private to this module.
!
!-----------------------------------------------------------------------
!---------------------------Arguments ----------------------------------
  integer        , intent(in)    :: lchnk   ! Chunk index
  type(ice_out_t), intent(inout) :: ice_out
  real(r8)       , intent(in)    :: landfrac(pcols)
  type(ice_in_t) , intent(in)    :: ice_in
!---------------------------Local variables-----------------------------
  integer i,k              ! loop indices
  integer ncol             ! number of columns in current chunk
  integer :: lats(pcols)   ! array of latitude indices
!-----------------------------------------------------------------------
!
! Set initial ice surface values
!
  ncol = get_ncols_p(lchnk)
  call get_lat_all_p(lchnk, ncol, lats)
  do i=1,ncol
    if (ice_out%areafrac(i) > 0._r8) then
        if (clat(lats(i)).gt.0) then
           sicthk(i,lchnk) = 2.0_r8
        else
           sicthk(i,lchnk) = 1.0_r8
        endif
        if ( previcefrac(i,lchnk)==0._r8) then !newice
           snowhice(i,lchnk) = 0._r8
           tsice(i,lchnk) = Tfrezk
           do k=1,plevmx
              tssub(i,k,lchnk)=Tfrezk
           end do
        end if
     else !tssub for non ice boxes filled with open ocean temp
        snowhice(i,lchnk) = 0._r8
        sicthk(i,lchnk) = 0._r8
        tsice(i,lchnk) = 0._r8
        if ( landfrac(i) < 1._r8)then
           do k=1,plevmx
              tssub(i,k,lchnk)=ice_in%tsocn(i)+Tffresh
           end do
        else
           do k=1,plevmx
              tssub(i,k,lchnk)=Tfrezk
           end do
        end if
     end if
  end do
  return
end subroutine newiceproperties
#endif

!
!----------------------------------------------------------------------- 
!
   real(r8) function sea_ice_energy (ice_in, dtime, inout, deltae, deltaaice, snowfall)
     !----------------------------------------------------------------------- 
     !
     ! Arguments
     !
      type(ice_in_t), intent(in) :: ice_in(begchunk:endchunk)
      integer, intent(in) :: dtime
      integer, intent(in) :: inout
      real(r8), intent(inout) :: deltae(pcols,begchunk:endchunk)
      real(r8), intent(inout) :: deltaaice(pcols,begchunk:endchunk)
      real(r8), intent(in) :: snowfall(pcols,begchunk:endchunk)
      !
      ! Local workspace
      !
      real(r8) :: totenerg(pcols,begchunk:endchunk)
      real(r8) :: totenerg_field(plon,plat)
      real(r8) :: qi
      real(r8) :: wt
      real(r8) :: mean
      real(r8) :: tmp
      real(r8) :: wght
      real(r8) :: erate(pcols,begchunk:endchunk)
      real(r8) :: hs
      integer :: i, j, k, c, ncol
      character(len=8) :: varname
      integer :: ret
      !----------------------------------------------------------------------- 

      !      call t_startf ('sea_ice_energy')
      if (inout == 1) then
         varname = 'EICEIN'
      else
         varname = 'EICEOUT'
      end if
!$OMP PARALLEL DO PRIVATE (C, I, K, NCOL, QI, HS)
      do c = begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            totenerg(i,c) = 0._r8
            if (aice(i,c) > 0._r8) then
               qi = 0._r8
               do k=1,ni
                  qi = qi + energ (tssub(i,k,c)-Tffresh, saltz(k))
               end do
               qi = qi/ni
               if (inout == 1) then
                  hs = snowhice(i,c) + snowfall(i,c)*dtime
               else
                  hs = snowhice(i,c)
               end if
               totenerg(i,c) = (-rLfs*hs*rhofresh/rhos + qi*sicthk(i,c))*aice(i,c)
               ! if (inout == 1) then
               !    write(iulog,*)'i,c,snowterm2a=',i,c,(-rLfs*hs*rhofresh/rhos)*aice(i,c)
               !    write(iulog,*)'i,c,iceterm2a=',i,c,qi*sicthk(i,c)*aice(i,c)
               ! else
               !    write(iulog,*)'i,c,snowterm2b=',i,c,(-rLfs*hs*rhofresh/rhos)*aice(i,c)
               !    write(iulog,*)'i,c,iceterm2b=',i,c,qi*sicthk(i,c)*aice(i,c)
               ! end if
            end if
         end do
         if (inout == 1) then
            do i=1,ncol
               deltae(i,c)    = totenerg(i,c)
               deltaaice(i,c) = aice(i,c)
            end do
         else
            do i=1,ncol
               deltae(i,c)    = totenerg(i,c) - deltae(i,c)
               deltaaice(i,c) = aice(i,c) - deltaaice(i,c)
               erate(i,c)     = deltae(i,c)/dtime
            end do
            call outfld ('DELTAICE', deltaaice(1,c), pcols, c)
            call outfld ('NRGICE  ', totenerg(1,c), pcols, c)
            call outfld ('IIERATE ', erate(1,c), pcols, c)
         end if
      end do

      call gather_chunk_to_field (1, 1, 1, plon, totenerg, totenerg_field)

      if (masterproc) then
         mean = 0._r8
         wght = 0._r8
         do j=1,plat
            wt = w(j)/nlon(j)
            do i=1,nlon(j)
               if (landfrac_glob(i,j) < 1._r8) then
                  tmp = wt*(1._r8 - landfrac_glob(i,j))
                  mean = mean + tmp*totenerg_field(i,j)
                  wght = wght + tmp
               end if
            end do
         end do
         mean = mean/wght
      end if
      call mpi_bcast( mean, 1, MPI_REAL8, 0, mpicom, ret )
      sea_ice_energy = mean
      !      call t_stopf ('sea_ice_energy')

   end function sea_ice_energy

!
!----------------------------------------------------------------------- 
!

   subroutine gmean_ice (ice_out, ice_in, fsns, snowfall, F_ice, &
                         F_ocn, F_frzmlt, dtime, deltae, evap, aiceinit)

     !----------------------------------------------------------------------- 
      !
      ! Arguments
      !
      type(ice_out_t), intent(in) :: ice_out(begchunk:endchunk)
      type(ice_in_t) , intent(in) :: ice_in(begchunk:endchunk)

      real(r8), intent(in) :: fsns(pcols,begchunk:endchunk)
      real(r8), intent(in) :: snowfall(pcols,begchunk:endchunk)
      real(r8), intent(out):: F_ice
      real(r8), intent(out):: F_ocn
      real(r8), intent(out):: F_frzmlt
      integer , intent(in) :: dtime
      real(r8), intent(in) :: deltae(pcols,begchunk:endchunk)
      real(r8), intent(in) :: evap(pcols,begchunk:endchunk)
      real(r8), intent(in) :: aiceinit(pcols,begchunk:endchunk)
      !
      ! Local workspace
      !
      real(r8) :: imbalance(pcols,begchunk:endchunk)
      real(r8) :: F_ice_chunk(pcols,begchunk:endchunk)
      real(r8) :: F_ocn_chunk(pcols,begchunk:endchunk)
      real(r8) :: F_frzmlt_chunk(pcols,begchunk:endchunk)
      real(r8) :: F_ice_field(plon,plat)
      real(r8) :: F_ocn_field(plon,plat)
      real(r8) :: F_frzmlt_field(plon,plat)
      real(r8) :: wght
      real(r8) :: wt
      real(r8) :: tmp
      real(r8) :: netflux
      integer :: i, c, j, ncol
      integer :: ret
      !----------------------------------------------------------------------- 

      call t_startf ('gmean_ice')

!$OMP PARALLEL DO PRIVATE (C, I, NCOL, NETFLUX)
      do c = begchunk,endchunk
         ncol = get_ncols_p(c)
         do i=1,ncol
            F_ice_chunk(i,c) = 0._r8
            F_ocn_chunk(i,c) = 0._r8
            F_frzmlt_chunk(i,c) = 0._r8
            if (frac(c)%land(i) < 1._r8) then
               F_ocn_chunk(i,c) = Focn(i,c)
               F_frzmlt_chunk(i,c) = max (0._r8, ice_in(c)%frzmlt(i))
            end if
            if (aiceinit(i,c) > 0._r8) then
               !BPB note: CAM2 sign conventions on fluxes
               F_ice_chunk(i,c) = aiceinit(i,c)* &
                      (fsns(i,c) + ice_in(c)%flwds(i) - &
                       ice_out(c)%lwup(i) + evap(i,c)*Lvap - ice_out(c)%shf(i))
               ! write(iulog,*)'i,c,aiceinit2=',i,c,aiceinit(i,c)
               ! write(iulog,*)'i,c,aice2=',i,c,aice(i,c)
               ! write(iulog,*)'i,c,fswabs2=',i,c,fsns(i,c)
               ! write(iulog,*)'i,c,lwsum2=',i,c,ice_in(c)%flwds(i)-ice_out(c)%lwup(i)
               ! write(iulog,*)'i,c,shflx2=',i,c,ice_out(c)%shf(i)
               ! write(iulog,*)'i,c,evaplvap2=',i,c,evap(i,c)*Lvap
               ! write(iulog,*)'i,c,focn2=',i,c,F_ocn_chunk(i,c)
               ! write(iulog,*)'i,c,frzmlt2=',i,c,F_frzmlt_chunk(i,c)
            endif
            netflux = F_ice_chunk(i,c) - F_ocn_chunk(i,c) - F_frzmlt_chunk(i,c)
            imbalance(i,c) = deltae(i,c)/dtime - netflux
         end do
         call outfld ('F_ICE',    F_ice_chunk(1,c),    pcols, c)
         call outfld ('F_OCN',    F_ocn_chunk(1,c),    pcols, c)
         call outfld ('FRZMLTMX', F_frzmlt_chunk(1,c), pcols, c)
         call outfld ('IMBAL   ', imbalance(1,c), pcols, c)
      end do

      call gather_chunk_to_field (1, 1, 1, plon, F_ice_chunk, F_ice_field)
      call gather_chunk_to_field (1, 1, 1, plon, F_ocn_chunk, F_ocn_field)
      call gather_chunk_to_field (1, 1, 1, plon, F_frzmlt_chunk, F_frzmlt_field)

      if (masterproc) then
         F_ice    = 0._r8
         F_ocn    = 0._r8
         F_frzmlt = 0._r8
         wght     = 0._r8
         do j=1,plat
            wt = w(j)/nlon(j)
            do i=1,nlon(j)
               if (landfrac_glob(i,j) < 1._r8) then
                  tmp = wt*(1._r8 - landfrac_glob(i,j))
                  F_ice    = F_ice    + tmp*F_ice_field(i,j)
                  F_ocn    = F_ocn    + tmp*F_ocn_field(i,j)
                  F_frzmlt = F_frzmlt + tmp*F_frzmlt_field(i,j)
                  wght     = wght     + tmp
               end if
            end do
         end do
         F_ice    = F_ice/wght
         F_ocn    = F_ocn/wght
         F_frzmlt = F_frzmlt/wght
      end if
      call mpi_bcast (F_ice   , 1, MPI_REAL8, 0, mpicom, ret)
      call mpi_bcast (F_ocn   , 1, MPI_REAL8, 0, mpicom, ret)
      call mpi_bcast (F_frzmlt, 1, MPI_REAL8, 0, mpicom, ret)

      call t_stopf ('gmean_ice')

   end subroutine gmean_ice
!
!----------------------------------------------------------------------- 
!
    real(r8) function ice_fraction( aice, landfrac )
         real(r8), intent(in) :: aice       ! Areal fraction of sea-ice over open-ocean
         real(r8), intent(in) :: landfrac   ! Areal fraction of grid-square covered by land
         ice_fraction = aice*(1._r8 - landfrac)
    end function ice_fraction
!
!----------------------------------------------------------------------- 
!

    subroutine read_domain(varname, ncid , dimlonnam, dimlatnam, lonnam, latnam, &
                       dim1b, dim1e, dim2b, dim2e, field, readvar)
      !-----------------------------------------------------------------------
      !
      ! Arguments
      !
      character(len=*), intent(in)  :: varname     ! variable name
      integer         , intent(in)  :: ncid        ! input unit
      character(len=*), intent(in)  :: dimlonnam   ! name of longitude dimension of field on file
      character(len=*), intent(in)  :: dimlatnam   ! name of latitude  dimension of field on file
      character(len=*), intent(in)  :: lonnam      ! name of longitude variable  of field on file
      character(len=*), intent(in)  :: latnam      ! name of latitude  variable  of field on file
      integer         , intent(in)  :: dim1b       ! start of first  dimension of array to be returned
      integer         , intent(in)  :: dim1e       ! end   of first  dimension of array to be returned
      integer         , intent(in)  :: dim2b       ! start of second dimension of array to be returned
      integer         , intent(in)  :: dim2e       ! end   of second dimension of array to be returned
      real(r8)        , intent(out) :: field(dim1b:dim1e,dim2b:dim2e) ! array to be returned (decomposed or global)
      logical         , intent(out) :: readvar     ! true => variable is on initial dataset
      !
      ! local variables
      !
      integer :: i,j                      ! index
      integer :: ier                      ! error status
      integer :: varid                    ! variable id
      integer :: dimlon, dimlat           ! lon, lat, lev dimension lengths
      integer :: tmptype
      integer :: ndims                    ! number of dimensions
      integer :: dims(NF_MAX_VAR_DIMS)    ! variable shape
      integer :: londimid, latdimid       ! Dimension ID's
      integer :: strt(3)                  ! start lon, lat, time indices for netcdf 2-d
      integer :: cnt (3)                  ! lon, lat, time counts for netcdf 2-d
      data strt/3*1/                      ! 
      data cnt /1,1,1/                    ! 2-d arrs
      real(r8), pointer :: tmp(:,:)       ! input data
      logical :: readvar_tmp              ! if true, variable is on tape
      character(len=NF_MAX_NAME) tmpname
      character(len=32) :: subname='read_domain' ! subroutine name
      !-----------------------------------------------------------------------
      !
      if (masterproc) then
         call check_var(ncid, varname, varid, readvar_tmp)
         if (readvar_tmp) then
            call wrap_inq_dimid  (ncid, dimlonnam, londimid)
            call wrap_inq_dimlen (ncid, londimid , dimlon)
            call wrap_inq_dimid  (ncid, dimlatnam, latdimid)
            call wrap_inq_dimlen (ncid, latdimid , dimlat)

            ! Check order of dimensions in variable
            call wrap_inq_var (ncid, varid, tmpname, tmptype, ndims, dims , ier)
            if (dims(1) /= londimid .or. dims(2) /= latdimid .or. ndims > 3) then
               write(iulog,*) trim(subname), ' Error: Bad number of dims or ordering while reading field ', trim(varname)
               call endrun()
            end if

            ! Allocate memory and read variable
            cnt(1) = dimlon
            cnt(2) = dimlat
            allocate ( tmp(dimlon,dimlat) )
            call wrap_get_vara_realx (ncid, varid, strt, cnt, tmp)
         end if  ! end of readvar_tmp
      end if  ! end masterproc
      call mpi_bcast(readvar_tmp, 1, MPI_LOGICAL, 0, mpicom, ier)
      if (ier /= 0) then
         write(iulog,*) trim(subname),' Error:  broadcast error while reading ', trim(varname)
         call endrun()
      end if
      if (readvar_tmp) then
         call mpi_bcast(dimlon, 1, MPI_INTEGER, 0, mpicom, ier)
         if (ier /= 0) then
            write(iulog,*) trim(subname),' Error:  broadcast error while reading ', trim(varname)
            call endrun()
         end if
         call mpi_bcast(dimlat, 1, MPI_INTEGER, 0, mpicom, ier)
         if (ier /= 0) then
            write(iulog,*) trim(subname),' Error:  broadcast error while reading ', trim(varname)
            call endrun()
         end if
         if (.not. masterproc) then
            allocate ( tmp(dimlon,dimlat) )
         end if
         call scatter_field_to_chunk(1,1,1,dimlon,tmp,field)
         deallocate (tmp)
      end if
      readvar = readvar_tmp

    end subroutine read_domain
!
!----------------------------------------------------------------------- 
!
    subroutine check_var(ncid, varname, varid, readvar)
      !-----------------------------------------------------------------------
      !
      ! Arguments
      !
      integer, intent(in)          :: ncid
      character(len=*), intent(in) :: varname
      integer, intent(out)         :: varid
      logical, intent(out)         :: readvar 
      !
      ! Local Variables
      !
      integer :: ret     ! return value
      !-----------------------------------------------------------------------
      readvar = .true.
      ret = nf_inq_varid (ncid, varname, varid)
      if (ret/=NF_NOERR) then
         if (masterproc) then
            write(iulog,*)'CHECK_VAR Warning:  variable ',trim(varname),' is not on initial dataset'
         end if
         readvar = .false.
      end if
    end subroutine check_var

end module ice_comp





