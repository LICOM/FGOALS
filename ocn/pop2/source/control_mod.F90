!=======================================================================
! CVS $Id: control_mod.F90,v 1.2 2001/05/18 20:19:42 mvertens Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/ocn/docn5/control_mod.F90,v $
! CVS $Name: ccsm2_0_beta58 $
!=======================================================================

  module control_mod

  use shr_kind_mod     ! precision
  implicit none

  !----- input namelist parms --------------------------------------
  character(len=16) :: rest_type    ! run type: initial, cont, regen, branch
  character(len=64) :: rest_pfile   ! IC data pointer file
  character(len=64) :: rest_bfile   ! IC data file for branch runs
  character(len=16) :: rest_freq    ! restart file frequency
  integer           :: rest_date    ! start date for initial or branch runs
  integer           :: rest_nday    ! restart file interval    (nday option
  integer           :: rest_odate   ! restart file offset date (nday option)
  integer           :: rest_lag     ! restart date: lag by 1 sync interval?
  character(len=16) :: case_name    ! case name
  character(len=64) :: case_desc    ! case description
  character(len=256):: data_dir     ! data: directory for input files
  character(len=256):: domain_file  ! domain data file
  integer           :: domain_nx    ! nx for domain fabrication option
  integer           :: domain_ny    ! nx for domain fabrication option
  character(len=256):: data_file    ! sst data file
  character(len=16) :: data_form    ! either "annual" or "multiyear"
  character(len=16) :: data_sstname ! name of sst field on input data file
  character(len=16) :: data_lonname ! name of longitude coordinate on input data file
  character(len=16) :: data_latname ! name of latitude  coordinate on input data file
  character(len=64) :: mss_dir      ! MSS file directory
  character(len=16) :: mss_pass     ! MSS file password
  character(len=16) :: mss_opts     ! MSS file mswrite options
  integer           :: mss_rtpd     ! MSS file MS retention period
  integer           :: mss_rmlf     ! T => remove local file after mswrite
  integer           :: info_dbug    ! dbug level: 0=lowest, 3=highest
  integer           :: info_time    ! T => print extra timing info
  real(SHR_KIND_R8)          :: info_sleep   ! simulate time to advance an active model 
  integer           :: ncpl         ! number of msg pairs per day

  !----- timing info -----------------------------------------------
  integer       :: eday         ! elapsed days, 0 => 1 Jan year 0
  integer       :: eday0        ! value of eday on startup
  integer       :: cdate        ! coded model date (yyyymmdd)
  integer       :: year         ! model year
!  integer       :: month        ! model month
  integer       :: day          ! model day
  integer       :: sec          ! model elapsed seconds during model day

  !----- control flags ---------------------------------------------
  integer       :: stop_now     ! T => stop model now
  integer       :: stop_eod     ! T => stop model at end of day
  integer       :: rest_now     ! T => create restart data now
  integer       :: rest_eod     ! T => create restart data at end of day
  integer       :: hist_now     ! T => create history data now
  integer       :: hist_eod     ! T => create history data at end of day
  integer       :: hist_tav     ! T => create monthly avg data
  integer       :: diag_now     ! T => create diagnostic data now
  integer       :: diag_eod     ! T => create diagnostic data now

  integer, private :: ierr         ! error code
  SAVE

end MODULE control_mod
