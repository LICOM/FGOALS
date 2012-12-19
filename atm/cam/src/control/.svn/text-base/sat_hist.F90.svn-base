!-------------------------------------------------------------------------------
! Outputs history field columns as specified by a satellite track data file
!
! Created by Francis Vitt -- 17 Sep 2010
!-------------------------------------------------------------------------------
module sat_hist

  use perf_mod,      only: t_startf, t_stopf
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use cam_logfile,   only: iulog
  use ppgrid,        only: pcols, pver, begchunk, endchunk
  use cam_history_support, only: fieldname_lenp2, max_string_len, ptapes
  use spmd_utils,    only: masterproc, iam
  use abortutils,    only: endrun

  use pio,           only: file_desc_t, iosystem_desc_t, iosystem_desc_t, var_desc_t, io_desc_t
  use pio,           only: pio_openfile, pio_redef, pio_enddef, pio_inq_dimid, pio_inq_varid, pio_seterrorhandling, pio_def_var
  use pio,           only: pio_inq_dimlen, pio_get_att, pio_put_att, pio_get_var, pio_put_var, pio_write_darray
  use pio,           only: pio_real, pio_int, pio_double, pio_copy_att
  use pio,           only: PIO_WRITE,PIO_NOWRITE, PIO_NOERR, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR, PIO_Rearr_box, PIO_GLOBAL
  use spmd_utils,    only: mpicom
#ifdef SPMD
  use mpishorthand,  only: mpichar, mpiint
#endif
  
  implicit none

  private

  public :: sat_hist_readnl
  public :: sat_hist_init
  public :: sat_hist_write
  public :: sat_hist_define
  public :: is_satfile

  character(len=max_string_len)  :: sathist_track_infile
  type(file_desc_t), save :: infile


  integer :: half_step
  logical :: has_sat_hist = .false.

  real(r8), allocatable :: obs_lats(:)
  real(r8), allocatable :: obs_lons(:)

  logical  :: doy_format
  real(r8) :: first_datetime
  real(r8) :: last_datetime
  real(r8) :: previous_datetime
  integer  :: time_ndx
  integer  :: t_buffer_size
  integer, allocatable :: date_buffer(:), time_buffer(:)
  integer :: sat_tape_num=ptapes-1


      ! input file
      integer :: n_profiles
      integer :: time_vid, date_vid, lat_vid, lon_vid, instr_vid, orbit_vid, prof_vid, zenith_vid
      
      integer :: in_julian_vid
      integer :: in_localtime_vid
      integer :: in_doy_vid
      integer :: in_occ_type_vid

      integer :: in_start_col


  ! output file
      type(var_desc_t) :: out_latid, out_lonid, out_instrid, out_zenithid, out_orbid, out_profid
      type(var_desc_t) :: out_instr_lat_vid, out_instr_lon_vid
      type(var_desc_t) :: out_obs_date_vid, out_obs_time_vid
      type(var_desc_t) :: out_julian_vid
      type(var_desc_t) :: out_localtime_vid
      type(var_desc_t) :: out_doy_vid
      type(var_desc_t) :: out_occ_type_vid

contains
  
!-------------------------------------------------------------------------------

  logical function is_satfile (file_index)
    integer, intent(in) :: file_index ! index of file in question
    is_satfile = file_index == sat_tape_num
  end function is_satfile

!-------------------------------------------------------------------------------
  subroutine sat_hist_readnl(nlfile, hfilename_spec, mfilt, fincl, nhtfrq, avgflag_pertape)
    
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use cam_history_support, only : pflds
    implicit none
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    character(len=*), intent(inout) :: hfilename_spec(:)
    character(len=*), intent(inout) :: fincl(:,:)
    character(len=1), intent(inout) :: avgflag_pertape(:)
    integer,          intent(inout) :: mfilt(:), nhtfrq(:)
    
    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'sat_hist_readnl'
    integer :: f, fcnt

    character(len=fieldname_lenp2) :: sathist_fincl(pflds)
    character(len=max_string_len)  :: sathist_hfilename_spec
    integer :: sathist_mfilt, sat_tape_num

    namelist /satellite_options_nl/ sathist_track_infile, sathist_hfilename_spec, sathist_fincl, sathist_mfilt

    ! set defaults

    sathist_track_infile = ' '
    sathist_hfilename_spec = '%c.cam2.hs.%y-%m-%d-%s.nc'
    sathist_fincl(:) = ' '
    sathist_mfilt = 100000

    !read namelist options

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'satellite_options_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, satellite_options_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! broadcast the options to all MPI tasks
    call mpibcast(sathist_track_infile,   len(sathist_track_infile),   mpichar, 0, mpicom)
    call mpibcast(sathist_hfilename_spec, len(sathist_hfilename_spec), mpichar, 0, mpicom)
    call mpibcast(sathist_fincl,          pflds*len(sathist_fincl(1)), mpichar, 0, mpicom)
    call mpibcast(sathist_mfilt,          1,                       mpiint,  0, mpicom)
#endif

    has_sat_hist = len_trim(sathist_track_infile) > 0

    if (.not.has_sat_hist) return

     sat_tape_num=ptapes-1
     hfilename_spec(sat_tape_num) = sathist_hfilename_spec
     mfilt(sat_tape_num) = sathist_mfilt
     fcnt=0
     do f=1, pflds
        fincl(f,sat_tape_num) = sathist_fincl(f)
        if(len_trim(sathist_fincl(f)) > 0) then
           fcnt=fcnt+1
        end if
     enddo
     
     nhtfrq(sat_tape_num) = 1
     avgflag_pertape(sat_tape_num) = 'I'

     if(masterproc) then
        write(iulog,*) 'sathist_track_infile: ',trim(sathist_track_infile)
        write(iulog,*) 'sathist_hfilename_spec: ',trim(sathist_hfilename_spec)
        write(iulog,*) 'sathist_fincl: ',(trim(sathist_fincl(f))//' ', f=1,fcnt)
        write(iulog,*) 'max columns per file sathist_mfilt: ',sathist_mfilt
     end if

   end subroutine sat_hist_readnl

  
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine sat_hist_init
    use cam_pio_utils, only: cam_pio_openfile
    use ioFileMod,     only: getfil
    use spmd_utils,    only: npes
    use time_manager,  only: get_step_size
    use string_utils,  only: to_lower, GLC

    implicit none

    character(len=max_string_len)  :: locfn       ! Local filename
    integer :: ierr, dimid, i
    integer :: num_iotasks, io_stride

    character(len=128) :: date_format

    if (.not.has_sat_hist) return

    call getfil (sathist_track_infile, locfn)
    call cam_pio_openfile(infile, locfn, PIO_NOWRITE)

    ierr = pio_inq_dimid(infile,'profs',dimid)
    ierr = pio_inq_dimlen(infile, dimid, n_profiles)

    ierr = pio_inq_varid( infile, 'time', time_vid )
    ierr = pio_inq_varid( infile, 'date', date_vid )

    ierr = pio_get_att( infile, date_vid, 'long_name', date_format)
    date_format = to_lower(trim( date_format(:GLC(date_format))))

    if ( index( date_format, 'yyyymmdd') > 0 ) then
       doy_format = .false.
    else if  ( index( date_format, 'yyyyddd') > 0 ) then
       doy_format = .true.
    else
       call endrun('sat_hist_init: date_format not recognized : '//trim(date_format))
    endif

    ierr = pio_inq_varid( infile, 'lat', lat_vid )
    ierr = pio_inq_varid( infile, 'lon', lon_vid )

    call pio_seterrorhandling(infile, PIO_BCAST_ERROR)
    ierr = pio_inq_varid( infile, 'instr_num', instr_vid )
    if(ierr/=PIO_NOERR) instr_vid=-1

    ierr = pio_inq_varid( infile, 'orbit_num', orbit_vid )
    if(ierr/=PIO_NOERR) orbit_vid=-1

    ierr = pio_inq_varid( infile, 'prof_num',  prof_vid )
    if(ierr/=PIO_NOERR) prof_vid=-1

    ierr = pio_inq_varid( infile, 'instr_sza', zenith_vid )
    if(ierr/=PIO_NOERR) zenith_vid=-1

    ierr = pio_inq_varid( infile, 'julian', in_julian_vid )
    if(ierr/=PIO_NOERR) in_julian_vid=-1

    ierr = pio_inq_varid( infile, 'local_time', in_localtime_vid )
    if(ierr/=PIO_NOERR) in_localtime_vid=-1

    ierr = pio_inq_varid( infile, 'doy', in_doy_vid )
    if(ierr/=PIO_NOERR) in_doy_vid=-1

    ierr = pio_inq_varid( infile, 'occ_type', in_occ_type_vid )
    if(ierr/=PIO_NOERR) in_occ_type_vid=-1

    call pio_seterrorhandling(infile, PIO_INTERNAL_ERROR)

    call read_datetime( first_datetime, 1 )
    call read_datetime( last_datetime, n_profiles )
    previous_datetime = first_datetime
    t_buffer_size = min(1000,n_profiles)
    allocate( date_buffer(t_buffer_size), time_buffer(t_buffer_size) )

    if ( last_datetime<first_datetime ) then
       call endrun('sat_hist_init: satellite track file has invalid date time info')
    endif

    time_ndx = 1
    half_step = get_step_size()*0.5_r8

    num_iotasks = 1
    io_stride = npes

  end subroutine sat_hist_init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine read_datetime( datetime, index )

    real(r8), intent( out ) :: datetime
    integer,  intent( in )  :: index

    integer :: ierr
    integer :: cnt(1)
    integer :: start(1)
    integer :: date(1), time(1)

    cnt = (/ 1 /)
    start = (/index/)

    ierr = pio_get_var( infile, time_vid, start, cnt, time )
    ierr = pio_get_var( infile, date_vid, start, cnt, date )
    
    datetime = convert_date_time( date(1),time(1) )

  end subroutine read_datetime

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine read_buffered_datetime( datetime, index )

    real(r8), intent( out ) :: datetime
    integer,  intent( in )  :: index

    integer :: ii

    integer :: ierr
    integer :: cnt
    integer :: start
    integer :: date, time
    
    ii = mod( index, t_buffer_size )
    if ( ii == 0 ) then
       ii = t_buffer_size
    endif
    if ( ii == 1 ) then
       start = index
       if ( start+t_buffer_size-1 <= n_profiles ) then
          cnt = t_buffer_size 
       else
          cnt = n_profiles-start+1
       endif
       ierr = pio_get_var( infile, time_vid, (/ start /), (/ cnt /), time_buffer(1:cnt) )
       ierr = pio_get_var( infile, date_vid, (/ start /), (/ cnt /), date_buffer(1:cnt) )
    endif
    time = time_buffer(ii)
    date = date_buffer(ii)
    datetime = convert_date_time( date,time )

  end subroutine read_buffered_datetime

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  function convert_date_time( date,time )
    use time_manager, only: set_time_float_from_date

    integer, intent(in) :: date,time
    real(r8) :: convert_date_time

    real(r8) :: datetime
    integer :: yr, doy, mon, dom

    if ( doy_format ) then
       yr = date/1000
       doy = date - yr*1000
       call set_time_float_from_date( datetime, yr, 1, doy, time )
    else 
       yr = date/10000
       mon = (date - yr*10000)/100
       dom = date - yr*10000 - mon*100
       call set_time_float_from_date( datetime, yr, mon, dom, time )
    endif
    convert_date_time = datetime

  end function convert_date_time
!-------------------------------------------------------------------------------
  subroutine sat_hist_define(outfile)
    use pio, only : pio_inquire
    type(file_desc_t), intent(inout) :: outfile

    integer :: coldim
    integer :: ierr
    
    ierr = pio_inquire(outfile, unlimitedDimId=coldim)

    call pio_seterrorhandling(outfile, PIO_BCAST_ERROR)
    ierr = define_var( 'instr_lat', coldim, infile, lat_vid,  outfile, out_instr_lat_vid )
    ierr = define_var( 'instr_lon', coldim, infile, lon_vid,  outfile, out_instr_lon_vid )
    ierr = define_var( 'obs_time', coldim, infile, time_vid,  outfile, out_obs_time_vid )
    ierr = define_var( 'obs_date', coldim, infile, date_vid,  outfile, out_obs_date_vid )

    if (orbit_vid>0) then
       ierr = define_var( 'orbit_num', coldim, infile, orbit_vid,  outfile, out_orbid )
    endif
    if (prof_vid>0) then
       ierr = define_var( 'prof_num', coldim, infile, prof_vid,  outfile, out_profid )
    endif
    if (instr_vid>0) then
       ierr = define_var( 'instr_num', coldim, infile, instr_vid,  outfile, out_instrid )
    endif
    if (zenith_vid>0) then
       ierr = define_var( 'instr_sza', coldim, infile, zenith_vid,  outfile, out_zenithid )
    endif
    if (in_occ_type_vid>0) then
       ierr = define_var( 'occ_type', coldim, infile, in_occ_type_vid,  outfile, out_occ_type_vid )
    endif
    if (in_julian_vid>0) then
       ierr = define_var( 'julian', coldim, infile, in_julian_vid,  outfile, out_julian_vid )
    endif
    if (in_localtime_vid>0) then
       ierr = define_var( 'local_time', coldim, infile, in_localtime_vid,  outfile, out_localtime_vid )
    endif
    if (in_doy_vid>0) then
       ierr = define_var( 'doy', coldim, infile, in_doy_vid,  outfile, out_doy_vid )
    endif

    call pio_seterrorhandling(outfile, PIO_INTERNAL_ERROR)
    ierr=pio_put_att (outfile, PIO_GLOBAL, 'satellite_track_file', sathist_track_infile)
  end subroutine sat_hist_define


!-------------------------------------------------------------------------------
  subroutine sat_hist_write( tape , nflds, nfils)

    use ppgrid,   only : pcols, begchunk, endchunk
    use dyn_grid, only : get_dyn_grid_parm
    use cam_pio_utils, only: phys_decomp, dyn_decomp
    use cam_history_support, only : active_entry
    use pio, only : pio_file_is_open
    implicit none
    type(active_entry) :: tape
    integer, intent(in) :: nflds
    integer, intent(inout) :: nfils

    integer :: t, f, i, ncols    
    integer :: ierr

    integer, allocatable :: col_ndxs(:)
    integer, allocatable :: chk_ndxs(:)
    integer, allocatable :: fdyn_ndxs(:)
    integer, allocatable :: ldyn_ndxs(:)
    integer, allocatable :: phs_owners(:)
    integer, allocatable :: dyn_owners(:)
    real(r8),allocatable :: mlats(:)
    real(r8),allocatable :: mlons(:)

    integer, pointer :: dof(:)

    integer :: coldim

    integer :: io_type


    if (.not.has_sat_hist) return

    call read_next_position( ncols )

    if ( ncols < 1 ) return

    call t_startf ('sat_hist_write')


    allocate( col_ndxs(ncols) )
    allocate( chk_ndxs(ncols) )
    allocate( fdyn_ndxs(ncols) )
    allocate( ldyn_ndxs(ncols) )
    allocate( phs_owners(ncols) )
    allocate( dyn_owners(ncols) )
    allocate( mlats(ncols) )
    allocate( mlons(ncols) )

    call get_indices( obs_lats, obs_lons, ncols, col_ndxs, chk_ndxs, fdyn_ndxs, ldyn_ndxs, phs_owners, dyn_owners, mlats, mlons )


    if ( .not. pio_file_is_open(tape%File) ) then
       call endrun('sat file not open')
    endif

    ierr = pio_inq_dimid(tape%File,'ncol',coldim )
    
    ierr = pio_inq_varid(tape%File, 'lat', out_latid )
    ierr = pio_inq_varid(tape%File, 'lon', out_lonid )


    call write_record_coord( tape, mlats, mlons, ncols, nfils )

    do f=1,nflds

       select case (tape%hlist(f)%field%decomp_type)
       case (phys_decomp)
          call dump_columns(tape%File, tape%hlist(f), ncols, nfils, col_ndxs, chk_ndxs, phs_owners )
       case (dyn_decomp)
          call dump_columns(tape%File, tape%hlist(f), ncols, nfils, fdyn_ndxs, ldyn_ndxs, dyn_owners )
       end select

    enddo

    deallocate( col_ndxs, chk_ndxs, fdyn_ndxs, ldyn_ndxs, phs_owners, dyn_owners )
    deallocate( mlons, mlats )
    deallocate( obs_lons, obs_lats )

    nfils = nfils + ncols

    call t_stopf ('sat_hist_write')

  end subroutine sat_hist_write

!-------------------------------------------------------------------------------
  subroutine dump_columns( File, hitem, ncols, nfils, fdims, ldims, owners  )
    use cam_history_support,  only: field_info, hentry
    use pionfwrite_mod, only: write_nf
    use cam_pio_utils, only : fillvalue
    use pio,            only: pio_initdecomp, pio_freedecomp, pio_setframe, pio_offset, pio_iam_iotask

    type(File_desc_t),intent(inout)  :: File
    type(hentry),     intent(in)     :: hitem
    integer,          intent(in)     :: ncols
    integer,          intent(in)     :: nfils
    integer,          intent(in)     :: fdims(:)
    integer,          intent(in)     :: ldims(:)
    integer,          intent(in)     :: owners(:)

    type(field_info) :: field
    type(var_desc_t) :: vardesc
    type(iosystem_desc_t), pointer :: sat_iosystem
    type(io_desc_t) :: iodesc
    integer :: t, ierr
    integer :: dimlens(2)

    real(r8), allocatable :: buf(:)
    integer,  allocatable :: dof(:)
    integer :: i,k, cnt

    call t_startf ('dump_columns')

    sat_iosystem => File%iosystem
    field = hitem%field
    vardesc = hitem%varid(1)

    dimlens(1) = field%numlev
    dimlens(2) = ncols
    allocate( buf( dimlens(1)*dimlens(2) ) )
    allocate( dof( dimlens(1)*dimlens(2) ) )

    vardesc%rec = -1
    cnt = 0
    
    buf = fillvalue
    dof = 0

    do i = 1,ncols
       do k = 1,field%numlev
          cnt = cnt+1
          if ( iam == owners(i) ) then
             buf(cnt) = hitem%hbuf( fdims(i), k, ldims(i) )
             dof(cnt) = cnt
          endif
       enddo
    enddo

    call pio_setframe(vardesc, int(-1,kind=PIO_OFFSET))

    if ( field%numlev>1 ) then
       call pio_initdecomp(sat_iosystem, pio_double, dimlens, dof, iodesc )
       if(pio_iam_iotask(sat_iosystem)) &
            iodesc%start(2)=iodesc%start(2)+nfils-1
    else
       call pio_initdecomp(sat_iosystem, pio_double, dimlens(2:2), dof, iodesc )
       if(pio_iam_iotask(sat_iosystem)) &
            iodesc%start(1)=iodesc%start(1)+nfils-1
    endif

    call pio_write_darray(File, vardesc, iodesc, buf, ierr, fillval=fillvalue)

    call pio_freedecomp(sat_iosystem, iodesc)

    deallocate( buf )
    deallocate( dof )
    

    call t_stopf ('dump_columns')

  end subroutine dump_columns

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine read_next_position( ncols )
    use time_manager, only: get_curr_date, get_prev_date
    use time_manager, only: set_time_float_from_date
   
    implicit none

    integer,  intent(out) :: ncols

    integer :: ierr
    integer :: yr, mon, day, tod
    real(r8) :: begdatetime, enddatetime
    integer :: beg_ndx, end_ndx, i

    real(r8) :: datetime

    call get_curr_date(yr, mon, day, tod)
    call set_time_float_from_date(begdatetime, yr, mon, day, tod-half_step)
    call set_time_float_from_date(enddatetime, yr, mon, day, tod+half_step)

    ncols = 0

    if ( first_datetime > enddatetime ) then
       if (masterproc) write(iulog,'(a)') 'sat_hist->read_next_position: all of the satellite date times are after the time window', first_datetime, enddatetime
       return
    endif
    if ( last_datetime < begdatetime ) then
       if (masterproc) write(iulog,'(a)') 'sat_hist->read_next_position: all of the satellite date times are before the time window', begdatetime, last_datetime
       return
    endif

    call t_startf ('read_next_position')

    beg_ndx = -99
    end_ndx = -99

    bnds_loop: do i = time_ndx,n_profiles

       call read_buffered_datetime( datetime, i )
       if ( datetime<previous_datetime ) then
          call endrun('sat_hist::read_next_position: datetimes are not valid')
       endif
       previous_datetime = datetime

       if ( datetime>begdatetime .and. beg_ndx<0 ) beg_ndx = i
       if ( datetime>enddatetime ) exit bnds_loop
       end_ndx = i

    enddo bnds_loop

    if (beg_ndx == -99 .and. end_ndx== -99) then
       if (masterproc) write(iulog,'(a)')  'sat_hist->read_next_position: must be beyond last position -- returning.'
       return
    endif

    if (end_ndx>0) time_ndx = end_ndx+1

    ncols = end_ndx-beg_ndx+1

    if (ncols > 0) then
       allocate( obs_lats(ncols), obs_lons(ncols) )
       in_start_col = beg_ndx

       ierr = pio_get_var( infile, lat_vid, (/beg_ndx/), (/ncols/), obs_lats )
       ierr = pio_get_var( infile, lon_vid, (/beg_ndx/), (/ncols/), obs_lons )

    endif

    call t_stopf ('read_next_position')
  end subroutine read_next_position
  
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine write_record_coord( tape, mod_lats, mod_lons, ncols, nfils )

    use time_manager,  only: get_nstep, get_curr_date, get_curr_time
    use cam_history_support, only : active_entry
    implicit none
    type(active_entry), intent(inout) :: tape

    integer,  intent(in) :: ncols
    real(r8), intent(in) :: mod_lats(ncols)
    real(r8), intent(in) :: mod_lons(ncols)
    integer,  intent(in) ::  nfils

    integer :: t, ierr
    integer :: yr, mon, day      ! year, month, and day components of a date
    integer :: nstep             ! current timestep number
    integer :: ncdate            ! current date in integer format [yyyymmdd]
    integer :: ncsec             ! current time of day [seconds]
    integer :: ndcur             ! day component of current time
    integer :: nscur             ! seconds component of current time
    real(r8) :: time             ! current time
    integer :: itmp(ncols)
    real(r8):: rtmp(ncols)

    call t_startf ('write_record_coord')

    nstep = get_nstep()
    call get_curr_date(yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day
    call get_curr_time(ndcur, nscur)


    time = ndcur + nscur/86400._r8

    itmp(:) = ncdate
    ierr = pio_put_var(tape%File, tape%dateid,(/nfils/), (/ncols/),itmp)
    itmp(:) = ncsec
    ierr = pio_put_var(tape%File, tape%datesecid,(/nfils/),(/ncols/),itmp)
    rtmp(:) = time
    ierr = pio_put_var(tape%File, tape%timeid, (/nfils/),(/ncols/),rtmp)

    ! output model column coordinates
    ierr = pio_put_var(tape%File, out_latid, (/nfils/),(/ncols/), mod_lats)
    ierr = pio_put_var(tape%File, out_lonid, (/nfils/),(/ncols/), mod_lons)
    
    ! output instrument location
    ierr = pio_put_var(tape%File, out_instr_lat_vid, (/nfils/),(/ncols/), obs_lats)
    ierr = pio_put_var(tape%File, out_instr_lon_vid, (/nfils/),(/ncols/), obs_lons)
    
    ierr = copy_data( infile, date_vid, tape%File, out_obs_date_vid, in_start_col, nfils, ncols )
    ierr = copy_data( infile, time_vid, tape%File, out_obs_time_vid, in_start_col, nfils, ncols )
    
    ! output observation identifiers
    if (instr_vid>0) then
       ierr = copy_data( infile, instr_vid, tape%File, out_instrid, in_start_col, nfils, ncols )
    endif
    if (orbit_vid>0) then
       ierr = copy_data( infile, orbit_vid, tape%File, out_orbid, in_start_col, nfils, ncols )
    endif
    if (prof_vid>0) then
       ierr = copy_data( infile, prof_vid, tape%File, out_profid, in_start_col, nfils, ncols )
    endif
    if (zenith_vid>0) then
       ierr = copy_data( infile, zenith_vid, tape%File, out_zenithid, in_start_col, nfils, ncols )
    endif
    if (in_julian_vid>0) then
       ierr = copy_data( infile, in_julian_vid, tape%File, out_julian_vid, in_start_col, nfils, ncols )
    endif
    if (in_occ_type_vid>0) then
       ierr = copy_data( infile, in_occ_type_vid, tape%File, out_occ_type_vid, in_start_col, nfils, ncols )
    endif
    if (in_localtime_vid>0) then
       ierr = copy_data( infile, in_localtime_vid, tape%File, out_localtime_vid, in_start_col, nfils, ncols )
    endif
    if (in_doy_vid>0) then
       ierr = copy_data( infile, in_doy_vid, tape%File, out_doy_vid, in_start_col, nfils, ncols )
    endif

    call t_stopf ('write_record_coord')
  end subroutine write_record_coord

!-------------------------------------------------------------------------------
! This returns the lat/lon information (and corresponding MPI task number (owner)) 
! of the global model grid column nearest to the input satellite coordinate (lat,lon)
!-------------------------------------------------------------------------------
  subroutine dyn_find_col( lat, lon, owner, first_ndx, last_ndx, rlat, rlon )

    use dyn_grid, only : get_dyn_grid_parm,get_horiz_grid_d,get_gcol_block_d,get_block_owner_d
    use dycore,           only: dycore_is
    use physconst,        only: pi 

    implicit none

    real(r8), intent(in)  :: lat, lon
    real(r8), optional, intent(out) :: rlat, rlon
    integer,  optional, intent(out) :: owner, first_ndx, last_ndx 

    integer :: i

    real(r8), allocatable :: clat_d(:), clon_d(:)
    real(r8) :: dist            ! the distance (in radians**2 from lat, lon)
    real(r8) :: distmin         ! the distance (in radians**2 from closest column)
    integer :: icol, ngcols
    integer :: blockid(1), bcid(1), lclblockid(1)
    integer :: plon, plat  
    real(r8), parameter :: rad2deg = 180._r8/pi




    call t_startf ('dyn_find_col')

    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')

    ngcols = plon*plat

    allocate( clat_d(1:ngcols) )
    allocate( clon_d(1:ngcols) )

    call get_horiz_grid_d(ngcols, clat_d_out=clat_d, clon_d_out=clon_d)

    clat_d(:) = rad2deg*clat_d(:)
    clon_d(:) = rad2deg*clon_d(:)

    icol    = -999
    distmin = 1.e10

    do i = 1,ngcols

       dist = (lat-clat_d(i))**2 + (lon-clon_d(i))**2
       if (dist < distmin ) then
          distmin = dist
          icol = i
       endif

    enddo

    if (present(rlat)) rlat = clat_d(icol)
    if (present(rlon)) rlon = clon_d(icol)
    
    call  get_gcol_block_d( icol, 1, blockid, bcid, lclblockid )

    if (present(first_ndx).and.present(last_ndx)) then
       if(dycore_is('UNSTRUCTURED')) then
          first_ndx = bcid(1)
          last_ndx = lclblockid(1)
       else
          last_ndx = (icol-1)/plon + 1 
          first_ndx = icol - (last_ndx-1)*plon
       endif
    endif

    if (present(owner)) then 
       owner = get_block_owner_d(blockid(1))
    endif

    deallocate( clat_d )
    deallocate( clon_d )

    call t_stopf ('dyn_find_col')
  end subroutine dyn_find_col

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine get_indices( lats, lons, ncols, col_ndxs, chk_ndxs, fdyn_ndxs, ldyn_ndxs, phs_owners, dyn_owners, mlats,mlons )

    use phys_grid,     only : phys_grid_find_col

    real(r8), intent(in)  :: lats(ncols)
    real(r8), intent(in)  :: lons(ncols)
    integer,  intent(in)  :: ncols
    integer,  intent(out) :: col_ndxs(ncols)
    integer,  intent(out) :: chk_ndxs(ncols)
    integer,  intent(out) :: fdyn_ndxs(ncols)
    integer,  intent(out) :: ldyn_ndxs(ncols)
    integer,  intent(out) :: phs_owners(ncols)
    integer,  intent(out) :: dyn_owners(ncols)
    real(r8), intent(out) :: mlats(ncols)
    real(r8), intent(out) :: mlons(ncols)
    
    integer :: i, ndx
    integer :: ichk,icol,idyn1,idyn2
    real(r8) :: rlat, rlon
    real(r8) :: lat, lon

    call t_startf ('get_indices')

    col_ndxs = -1
    chk_ndxs = -1
    fdyn_ndxs = -1
    ldyn_ndxs = -1
    phs_owners = -1
    dyn_owners = -1

    ndx = 0
    do i = 1,ncols

       lat = lats(i)
       lon = lons(i)

       if ( lon >= 360._r8) then
         lon = lon-360._r8
       endif
       if ( lon < 0._r8) then
         lon = lon+360._r8
       endif
       if (lat<-90._r8 .or. lat>90._r8) then
          write(iulog,*) 'sat_hist::get_indices lat = ',lat
          call endrun('sat_hist::get_indices : lat must be between -90 and 90 degrees (-90<=lat<=90)')
       endif
       if (lon<0._r8 .or. lon>=360._r8) then
          write(iulog,*) 'sat_hist::get_indices lon = ',lon
          call endrun('sat_hist::get_indices : lon must be between 0 and 360 degrees (0<=lon<360)')
       endif
       
       call phys_grid_find_col( lat, lon, phs_owners(i), ichk, icol )
       call dyn_find_col( lat, lon, owner=dyn_owners(i), first_ndx=idyn1, last_ndx=idyn2, rlat=rlat, rlon=rlon )  

       ndx = ndx+1         
       chk_ndxs(ndx) = ichk
       col_ndxs(ndx) = icol
       fdyn_ndxs(ndx) = idyn1
       ldyn_ndxs(ndx) = idyn2
       mlats(ndx) = rlat
       mlons(ndx) = rlon

    enddo

    call t_stopf ('get_indices')
  end subroutine get_indices

!-------------------------------------------------------------------------------
! utility function
!-------------------------------------------------------------------------------
  integer function define_var( var_name, coldim, infile, in_vid, outfile, out_id ) result(res)

    use pio, only: pio_inq_vartype

    character(len=*), intent(in) :: var_name
    integer,          intent(in) :: coldim
    type(File_desc_t),intent(inout) :: infile
    type(File_desc_t),intent(inout) :: outfile
    integer,          intent(in) :: in_vid
    type(var_desc_t), intent(out):: out_id

    integer :: type

    res = pio_inq_varid( outfile, var_name, out_id )
    if(res/=PIO_NOERR) then

       res = pio_inq_vartype( infile, in_vid, type )

       res = pio_def_var ( outfile, var_name, type, (/coldim/), out_id )

       res = copy_att( infile, in_vid, 'long_name', outfile, out_id )
       res = copy_att( infile, in_vid, 'units',     outfile, out_id )

    endif

  end function define_var

!-------------------------------------------------------------------------------
! utility function
!-------------------------------------------------------------------------------
  integer function copy_data( infile, in_vid, outfile, out_id, instart, outstart, ncols ) result(res)

    type(File_desc_t),intent(in) :: infile
    type(File_desc_t),intent(inout) :: outfile
    integer,          intent(in) :: in_vid
    type(var_desc_t), intent(in) :: out_id
    integer,          intent(in) :: instart, outstart, ncols

    real(r8), allocatable :: data(:)

    allocate( data(ncols) )

    res = pio_get_var( infile,  in_vid, (/instart/),  (/ncols/), data )
    res = pio_put_var( outfile, out_id, (/outstart/), (/ncols/), data )

    deallocate(data)

  end function copy_data

!-------------------------------------------------------------------------------
! utility function
! -- should be able to use pio_copy_att which does not seem to work
!-------------------------------------------------------------------------------
  integer function copy_att( infile, in_vid, att_name, outfile, out_id ) result(res)

    type(File_desc_t),intent(inout) :: infile
    type(File_desc_t),intent(inout) :: outfile
    character(len=*), intent(in) :: att_name
    integer,          intent(in) :: in_vid
    type(var_desc_t), intent(in) :: out_id

    character(len=128) :: att
    

    res = pio_get_att( infile, in_vid, trim(att_name), att )
    if (res==PIO_NOERR) then
       res = pio_put_att ( outfile, out_id, trim(att_name), trim(att))
    endif


  end function copy_att







end module sat_hist
