module date_expected
  use shr_kind_mod
  use shr_date_mod
  use shr_const_mod

  implicit none

  private

  public date_exp_data
  public date_exp_init
  public is_date_expected
  public is_prev_cdate_expected

  type date_exp_data
      integer(SHR_KIND_IN) :: y           ! calendar year
      integer(SHR_KIND_IN) :: m           ! calendar month
      integer(SHR_KIND_IN) :: d           ! calendar day
      integer(SHR_KIND_IN) :: s           ! elapsed seconds in current day
      integer(SHR_KIND_IN) :: cDate       ! coded calendar date (yymmdd)
      integer(SHR_KIND_IN) :: eDay        ! elsapsed days relative to calendar's reference date
      integer(SHR_KIND_IN) :: stepInDay   ! current time-step in current day
      integer(SHR_KIND_IN) :: stepsPerDay ! number of time-steps per day
  end type date_exp_data

  contains

  subroutine date_exp_init( date_exp, year, month, day, sec, stepsperday )
     implicit none

     type(date_exp_data), intent(OUT) :: date_exp
     integer, intent(IN) :: year, month, day, sec, stepsperday
 
     integer :: cdate, eday
     !--- this is the noleap calendar ---
     integer(SHR_KIND_IN),parameter :: dsm(12) = &     ! elapsed Days on Start of Month
   &                       (/ 0,31,59, 90,120,151, 181,212,243, 273,304,334/)

     date_exp%y    = year
     date_exp%m    = month
     date_exp%d    = day
     date_exp%s    = sec

     eday = year*365 + dsm(month) + day-1
     date_exp%eday = eday

     cdate = year*10000 + month*100 + day
     date_exp%cdate       = cdate
     date_exp%stepinday   = nint( sec * stepsperday / SHR_CONST_CDAY )
     date_exp%stepsperday = stepsperday

  end subroutine date_exp_init

  logical function is_date_expected( date, date_exp )
     implicit none

     type(shr_date),      intent(IN) :: date
     type(date_exp_data), intent(IN) :: date_exp

     integer :: year, month, day, sec, cdate, eday, stepsperday, stepinday

     is_date_expected = .true.
     call shr_date_getYMD( date, year, month, day, sec )
     if ( date_exp%y     /= year  ) is_date_expected = .false.
     if ( date_exp%m     /= month ) is_date_expected = .false.
     if ( date_exp%d     /= day   ) is_date_expected = .false.
     if ( date_exp%s     /= sec   ) is_date_expected = .false.
     call shr_date_getCDate( date, cdate, sec )
     if ( date_exp%cdate /= cdate ) is_date_expected = .false.
     call shr_date_getEDay( date, eday, sec )
     if ( date_exp%eday  /= eday  ) is_date_expected = .false.
     stepsperday = shr_date_getStepsPerDay( date )
     if ( date_exp%stepsperday  /= stepsperday  ) is_date_expected = .false.
     stepinday = shr_date_getStepInDay( date )
     if ( date_exp%stepinday  /= stepinday  ) is_date_expected = .false.
     
  end function is_date_expected

  logical function is_prev_cdate_expected( date, prev_date_exp )
     implicit none
     type(shr_date),      intent(IN) :: date
     type(date_exp_data), intent(IN) :: prev_date_exp

     integer :: cdate, sec

     is_prev_cdate_expected = .true.

     call shr_date_getCDate( date, cdate, sec, previous=.true. )
     if ( prev_date_exp%cdate /= cdate ) is_prev_cdate_expected = .false.
     if ( prev_date_exp%s     /= sec   ) is_prev_cdate_expected = .false.

  end function is_prev_cdate_expected

end module date_expected

program test_shr_date

  use shr_kind_mod, only : r8 => SHR_KIND_R8
  use shr_date_mod
  use test_mod
  use date_expected

  implicit none

  type(shr_date) :: date, next_date
  type(date_exp_data) :: date_exp, prev_date_exp
  real(r8) :: calday, calexp
  integer, parameter :: spday=48, nsteps = spday*366, dtime = 86400/spday
  integer :: i

  call test_init( 18 + nsteps )

  call date_exp_init( date_exp, year=2004, month=12, day=25, sec=3600, stepsperday=spday )
  date = shr_date_initYMD( y=2004, m=12, d=25, ns=spday, sec=3600)
  call test_is( is_date_expected( date, date_exp ), "test initYMD" )
  date = shr_date_initEDay( eday=2004*365+334+25-1, ns=spday )
  call shr_date_adv1step( date )
  call shr_date_adv1step( date )
  call test_is( is_date_expected( date, date_exp ), "test initEDay" )
  date = shr_date_initCDate( cdate=20041225, ns=spday, sec=3600 )
  call test_is( is_date_expected( date, date_exp ), "test initCDate" )

  calday = shr_date_getJulian( date )
  calexp = 334._r8+25._r8+(2._r8/real(spday,r8))
  write(*,*) 'calday = ', calday, ' expected = ', calexp
  call test_close( calday, calexp, 1.0e-14_r8, "test getJulian" )
  date = shr_date_initYMD( y=2004, m=1, d=1, ns=spday, sec=0)
  calday = shr_date_getJulian( date )
  calexp = 1.0
  write(*,*) 'calday = ', calday, ' expected = ', calexp
  call test_close( calday, calexp, 1.0e-14_r8, "test getJulian" )
  do i = 1, nsteps
     call shr_date_adv1step( date )
     calexp = calday + 1._r8 / real(spday,r8)
     if ( calexp >= 366.0_r8 ) calexp = calexp - 365.0_r8
     calday = shr_date_getJulian( date )
     write(*,*) '# ', i, 'calday = ', calday, ' expected = ', calexp
     call test_close( calday, calexp, 1.0e-10_r8, "test getJulian loop" )
  end do

  date = shr_date_initYMD( y=2004, m=12, d=25, ns=spday, sec=3600)
  prev_date_exp = date_exp
  call date_exp_init( date_exp, year=2004, month=12, day=25, sec=5400, stepsperday=spday)
  call shr_date_adv1step( date )
  call test_is( is_date_expected( date, date_exp ),  "test adv1step" )
  call test_is( .not. is_date_expected( date, prev_date_exp ), "test not same if advanced" )
  call test_is( is_prev_cdate_expected( date, prev_date_exp ), "test previous cdate" )

  date = shr_date_initCDate( cdate=20050101, sec=0, ns=spday)
  call date_exp_init( prev_date_exp, year=2004, month=12, day=31, sec=84600, &
                      stepsperday=spday)
  call test_is( is_prev_cdate_expected( date, prev_date_exp ), "test previous cdate for beginning of year/day" )
  call shr_date_adv1step( date )
  call date_exp_init( prev_date_exp, year=2005, month=1, day=1, sec=0, &
                      stepsperday=spday)
  call test_is( is_prev_cdate_expected( date, prev_date_exp ), "test previous cdate for 1-step after beginning of year/day" )

  next_date = date
  call test_is( date==next_date, "test date comparision works" )
  call shr_date_adv1step( next_date )
  call test_is( .not. date==next_date, "test date comparision is false if advanced" )
  call test_is( date<next_date, "test date lessthan works" )
  call test_is( next_date>date, "test date greaterhan works" )
  call test_is( .not. date>next_date, "test date greaterthan works if dates switched" )
  call test_is( .not. next_date<date, "test date lessthan works if dates switched" )

  date = shr_date_initYMD( y=2004, m=12, d=31, ns=spday, sec=86400-dtime)
  calday = shr_date_getJulian( date, shift=dtime )
  calexp = 1.0_r8
  call test_close( calday, calexp, 1.0e-10_r8, "test getJulian for end of year with an advance to next step" )

  date = shr_date_initYMD( y=2004, m=1, d=1, ns=spday, sec=25*dtime)
  calday = shr_date_getJulian( date, shift=2*dtime )
  calexp = 1.5625_r8
  call test_close( calday, calexp, 1.0e-10_r8, "test getJulian for a hour shift starting at an odd time-step" )

  call test_final( )

end program test_shr_date

