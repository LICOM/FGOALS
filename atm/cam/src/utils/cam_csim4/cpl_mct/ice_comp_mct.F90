module ice_comp_mct

  use shr_kind_mod,      only: r8 => shr_kind_r8
  use shr_sys_mod,       only: shr_sys_abort
  use physconst,         only: tmelt,pi
  use cam_logfile,       only: iulog

  use mct_mod
  use esmf_mod,          only: ESMF_Clock
  use seq_flds_mod
  use seq_flds_indices
  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod

  use ppgrid,            only: pcols, begchunk,endchunk
  use phys_grid   ,      only: get_ncols_p,get_rlat_all_p,get_rlon_all_p,get_area_all_p,&
                               ngcols, get_gcol_p

  use ice_types,         only: ice_in_t, ice_out_t, ice_types_alloc 
  use ice_comp,          only: ice_init, ice_alloc, ice_run, ice_write_restart, frac  
  use ice_comp,          only: ice_init, ice_alloc, ice_run, frac  
  use ice_spmd,          only: masterproc, iam
  use ice_time_manager,  only: get_curr_date, advance_timestep, get_nstep
  use perf_mod

  implicit none
  private
  save

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: ice_init_mct         ! Initialization method
  public :: ice_run_mct          ! Run method
  public :: ice_final_mct        ! Finalization method

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  private :: ice_SetgsMap_mct
  private :: ice_export_mct
  private :: ice_import_mct
  private :: ice_domain_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(ice_in_t) , pointer :: ice_in(:) 
  type(ice_out_t), pointer :: ice_out(:) 

!===============================================================
contains
!===============================================================

  subroutine ice_init_mct( EClock, cdata_i, x2i_i, i2x_i, NLFilename )

    !----------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock)         , intent(in)    :: EClock
    type(seq_cdata)          , intent(inout) :: cdata_i
    type(mct_aVect)          , intent(inout) :: x2i_i, i2x_i
    character(len=*), optional  , intent(in) :: NLFilename ! Namelist filename
    !
    ! Local variables
    !
    integer                               :: ICEID	
    integer                               :: mpicom_ice
    type(mct_gsMap)             , pointer :: gsMap_ice
    type(mct_gGrid)             , pointer :: dom_i
    type(seq_infodata_type)     , pointer :: infodata   ! Input init object
    integer                               :: lsize
    !----------------------------------------------------------
    !
    ! Set cdata pointers
    !
    call seq_cdata_setptrs(cdata_i, ID=ICEID, mpicom=mpicom_ice, &
         gsMap=gsMap_ice, dom=dom_i, infodata=infodata)
    ! 
    ! Initialize ice model
    !
    call ice_init( mpicom_ice, ice_in, ice_out, EClock, infodata )
    !
    ! Initialize ice gsMap
    !
    call ice_SetgsMap_mct( mpicom_ice, ICEID, gsMap_ice ) 	
    lsize = mct_gsMap_lsize(gsMap_ice, mpicom_ice)
    !
    ! Initialize mct ice domain (needs ice initialization info)
    !
    call ice_domain_mct( lsize, gsMap_ice, dom_i )
    !
    ! Inialize mct attribute vectors
    !
    call mct_aVect_init(x2i_i, rList=seq_flds_x2i_fields, lsize=lsize)
    call mct_aVect_zero(x2i_i)

    call mct_aVect_init(i2x_i, rList=seq_flds_i2x_fields, lsize=lsize)
    call mct_aVect_zero(i2x_i)
    !
    ! Create mct ice export state
    !
    call ice_export_mct( ice_out, i2x_i )
    call seq_infodata_PutData( infodata, ice_prognostic=.true.)

  end subroutine ice_init_mct
  
!==========================================================================

  subroutine ice_run_mct( EClock, cdata_i, x2i_i, i2x_i)

    !----------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock)            , intent(in)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_i
    type(mct_aVect)             , intent(inout) :: x2i_i
    type(mct_aVect)             , intent(inout) :: i2x_i
    !
    ! Local variables
    !
    logical :: rstwr           ! .true. ==> write a restart file
    integer :: ymd             ! Current date (YYYYMMDD)
    integer :: yr              ! Current year
    integer :: mon             ! Current month
    integer :: day             ! Current day
    integer :: tod             ! Current time of day (sec)
    integer :: ymd_sync        ! Current year of sync clock
    integer :: tod_sync        ! Time of day of sync clock
    integer :: stop_ymd        ! stop time (YYYYMMDD)
    integer :: stop_tod        ! stop time (sec)	
    integer :: yr_sync         ! Sync current year
    integer :: mon_sync        ! Sync current month
    integer :: day_sync        ! Sync current day
    integer :: RestartNextYMD  
    integer :: RestartNextTOD  
    character(len=*), parameter :: SubName = "ice_run_mct"
    !----------------------------------------------------------

    call t_startf ('ice_import')
    call ice_import_mct( x2i_i, ice_in )
    call t_stopf ('ice_import')

    call t_startf ('ice_run')
    call ice_run( ice_in, ice_out )
    call t_stopf ('ice_run')

    call t_startf ('ice_export')
    call ice_Export_mct (ice_out, i2x_i )
    call t_stopf ('ice_export')

    ! Determine restart and last time step flags and write restart if appropriate

    rstwr = seq_timemgr_RestartAlarmIsOn(EClock)

    if (rstwr) then
       call t_startf('ice_write_restart')
       call seq_timemgr_EClockGetData(EClock, curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync, curr_tod=tod_sync)
       call ice_write_restart( ice_out, &
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync)
       call t_stopf('ice_write_restart')
    endif
    
    ! Advance ocn timestep

    call advance_timestep()

    ! Check that internal clock in sync with master clock

    call get_curr_date(yr, mon, day, tod )
    ymd = yr*10000 + mon*100 + day

    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' camice ymd=',ymd     ,'  camice tod= ',tod
       write(iulog,*)'   sync ymd=',ymd_sync,'    sync tod= ',tod_sync
       call shr_sys_abort( SubName//":: Internal sea-ice clock not in sync with Sync Clock")
    end if

  end subroutine ice_run_mct

!==========================================================================

  subroutine ice_final_mct( )
    
     ! fill this in
  end subroutine ice_final_mct
  
!==========================================================================

  subroutine ice_export_mct(ice_out, i2x_i ) 

    !-------------------------------------------------------------------
    implicit none
    type(ice_out_t), intent(in)    :: ice_out(begchunk:endchunk) 
    type(mct_aVect), intent(inout) :: i2x_i

    integer :: i,c,ig         ! indices
    integer :: ncols          ! number of columns
    !-----------------------------------------------------------------------

    ig=1
    do c=begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          i2x_i%rAttr(index_i2x_Si_ifrac  ,ig) = ice_out(c)%aice(i)
          i2x_i%rAttr(index_i2x_Si_t      ,ig) = ice_out(c)%ts(i)      
          i2x_i%rAttr(index_i2x_Si_tref   ,ig) = ice_out(c)%tref(i)    
          i2x_i%rAttr(index_i2x_Si_avsdr  ,ig) = ice_out(c)%asdir(i)   
          i2x_i%rAttr(index_i2x_Si_avsdf  ,ig) = ice_out(c)%asdif(i)   
          i2x_i%rAttr(index_i2x_Si_anidr  ,ig) = ice_out(c)%aldir(i)   
          i2x_i%rAttr(index_i2x_Si_anidf  ,ig) = ice_out(c)%aldif(i)   
          i2x_i%rAttr(index_i2x_Faii_lwup ,ig) = ice_out(c)%lwup (i)   
          i2x_i%rAttr(index_i2x_Faii_lat  ,ig) = ice_out(c)%lhf(i)     
          i2x_i%rAttr(index_i2x_Faii_sen  ,ig) = ice_out(c)%shf(i)     
          i2x_i%rAttr(index_i2x_Faii_evap ,ig) = ice_out(c)%cflx(i)  
          i2x_i%rAttr(index_i2x_Faii_taux ,ig) = ice_out(c)%wsx(i)     
          i2x_i%rAttr(index_i2x_Faii_tauy ,ig) = ice_out(c)%wsy(i)     
          i2x_i%rAttr(index_i2x_Faii_swnet,ig) = ice_out(c)%fswabs(i)     
          if (ice_out(c)%aice(i) > 0._r8) then
             ! need to divide by aice since the merge to the ocean multiplies this back in
             i2x_i%rAttr(index_i2x_Fioi_melth,ig) = ice_out(c)%focn(i)/ice_out(c)%aice(i) 
          else
             i2x_i%rAttr(index_i2x_Fioi_melth,ig) = 0._r8
          end if
          ig=ig+1
       end do
    end do

  end subroutine ice_export_mct

!==========================================================================

  subroutine ice_import_mct( x2i_i, ice_in )

    !-----------------------------------------------------------------------
    implicit none
    type(mct_aVect), intent(inout) :: x2i_i
    type(ice_in_t) , intent(inout) :: ice_in(begchunk:endchunk)

    integer  :: i,c,ig        ! indices
    integer  :: ncols         ! number of columns
    !-----------------------------------------------------------------------

    ig=1
    do c=begchunk,endchunk
       ncols= get_ncols_p(c)
       do i=1,ncols
          ice_in(c)%tsocn(i)  =  x2i_i%rAttr(index_x2i_So_t      ,ig) - tmelt     
          ice_in(c)%tbot(i)   =  x2i_i%rAttr(index_x2i_Sa_tbot   ,ig)    
          ice_in(c)%zbot(i)   =  x2i_i%rAttr(index_x2i_Sa_z      ,ig)       
          ice_in(c)%ubot(i)   =  x2i_i%rAttr(index_x2i_Sa_u      ,ig)       
          ice_in(c)%vbot(i)   =  x2i_i%rAttr(index_x2i_Sa_v      ,ig)       
          ice_in(c)%pbot(i)   =  x2i_i%rAttr(index_x2i_Sa_pbot   ,ig)       
          ice_in(c)%qbot(i)   =  x2i_i%rAttr(index_x2i_Sa_shum   ,ig)    
          ice_in(c)%thbot(i)  =  x2i_i%rAttr(index_x2i_Sa_ptem   ,ig)    
          ice_in(c)%flwds(i)  =  x2i_i%rAttr(index_x2i_Faxa_lwdn ,ig)  
          ice_in(c)%snow(i)   =  x2i_i%rAttr(index_x2i_Faxa_snow ,ig)/1000._r8 
          ice_in(c)%soll(i)   =  x2i_i%rAttr(index_x2i_Faxa_swndr,ig)      
          ice_in(c)%sols(i)   =  x2i_i%rAttr(index_x2i_Faxa_swvdr,ig)      
          ice_in(c)%solld(i)  =  x2i_i%rAttr(index_x2i_Faxa_swndf,ig)     
          ice_in(c)%solsd(i)  =  x2i_i%rAttr(index_x2i_Faxa_swvdf,ig)     
          ice_in(c)%frzmlt(i) =  x2i_i%rAttr(index_x2i_Fioo_q    ,ig) 
          ig=ig+1
       end do
    end do

  end subroutine ice_import_mct

!==========================================================================

  subroutine ice_SetgsMap_mct( mpicom_ice, ICEID, gsMap_ice )
    use phys_grid, only : get_nlcols_p
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)  :: mpicom_ice
    integer        , intent(in)  :: ICEID
    type(mct_gsMap), intent(out) :: gsMap_ice
    !
    ! Local variables
    !
    integer :: i,j,sizebuf,n,c,ncols, nlcols  ! indices
    integer :: ier                    ! error status
    integer, allocatable :: gindex(:)
    !-------------------------------------------------------------------

    ! Determine global seg map

    sizebuf=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i = 1,ncols
          sizebuf = sizebuf+1
       end do
    end do

    allocate(gindex(sizebuf))

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i = 1,ncols
          n=n+1
          gindex(n) = get_gcol_p(c,i)
       end do
    end do
    nlcols = get_nlcols_p()
    call mct_gsMap_init( gsMap_ice, gindex, mpicom_ice, ICEID, nlcols, ngcols )

    deallocate(gindex)

  end subroutine ice_SetgsMap_mct

!===============================================================================

  subroutine ice_domain_mct( lsize, gsMap_i, dom_i )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_i
    type(mct_ggrid), intent(inout) :: dom_i     
    !
    ! Local Variables
    !
    integer  :: n,j,i,c,ncols         ! indices	
    real(r8) :: lats(pcols)           ! array of global latitude indices
    real(r8) :: lons(pcols)           ! array of global longitude indices
    real(r8) :: area(ngcols)          ! area in radians squared for each grid point
    real(r8), pointer :: data(:)      ! temporary
    integer , pointer :: idata(:)     ! temporary
    real(r8), parameter :: radtodeg = 180.0_r8/pi
    !-------------------------------------------------------------------
    !
    ! Initialize domain type
    !
    call mct_gGrid_init( GGrid=dom_i, CoordChars=trim(seq_flds_dom_coord), &
      OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsMap_i, iam, idata)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_i,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_i,"mask" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"frac" ,data,lsize) 
    !
    ! Fill in correct values for domain components
    !
    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats)
       do i=1,ncols
          n = n+1
          data(n) = lats(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_i,"lat",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlon_all_p(c, ncols, lons)
       do i=1,ncols
          n = n+1
          data(n) = lons(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_i,"lon",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_area_all_p(c, ncols, area)
       do i=1,ncols
          n = n+1
          data(n) = area(i) 
       end do
    end do
    call mct_gGrid_importRAttr(dom_i,"area",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          n = n+1
          if (frac(c)%land(i) < 1._r8) then
             data(n) = 1._r8 ! mask
          else
             data(n) = 0._r8
          end if
       end do
    end do
    call mct_gGrid_importRAttr(dom_i,"mask",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          n = n+1
          data(n) = 1._r8 - frac(c)%land(i)
       end do
    end do
    call mct_gGrid_importRAttr(dom_i,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine ice_domain_mct

end module ice_comp_mct
