!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_cal_mod -- calendar module, relates elapsed days to calendar date.
!
! !DESCRIPTION:
!   These calendar routines do conversions between...
!   \begin{itemize}
!   \item the integer number of elapsed days 
!   \item the integers year, month, day (three inter-related integers)
!   \item the integer coded calendar date (yyyymmdd)
!   \end{itemize}
!   Possible uses include: a calling routine can increment the elapsed days 
!   integer and use this module to determine what the corresponding calendar 
!   date is;  this module can be used to determine how many days apart two
!   arbitrary calendar dates are.
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - created initial version, taken from cpl5
!
! !REMARKS:
!   Following are some internal assumptions.  These assumptions are somewhat
!   arbitrary -- they were chosen because they result in the simplest code given
!   the requirements of this module.  These assumptions can be relaxed as 
!   necessary: 
!   o the valid range of years is [-999,9999]
!   o elapsed days = 0 <=> January 1st, year 0000
!   o all years have 365 days (no leap years)
!     This module is hard-coded to implement a 365-day calendar, ie. there
!     are no leap years.  This module can be modified to implement a calendar
!     with leap years if this becomes desireable.  This would make the internal
!     logic of this module more complex, but would not require any change to the
!     module API or the calling code because the module API hides these details
!     from all external routines.
!
! !INTERFACE: ------------------------------------------------------------------

module shr_cal_mod

! !USES:

   use shr_kind_mod   ! kinds
   use shr_const_mod  ! constants
   use shr_sys_mod    ! system
   use shr_string_mod, only: shr_string_toLower
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit
   use esmf_mod

   implicit none

   private ! except

! !PUBLIC TYPES: 
 
   ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_cal_calendarName       ! checks and returns valid cal name
   public :: shr_cal_numDaysinMonth     ! number of days in a month
   public :: shr_cal_numDaysinYear      ! number of days in a year
   public :: shr_cal_elapsDaysStrtMonth ! elapsed days on start of month
   public :: shr_cal_eday2date  ! converts elapsed days to coded-date
   public :: shr_cal_eday2ymd   ! converts elapsed days to yr,month,day
   public :: shr_cal_date2ymd   ! converts coded-date   to yr,month,day
   public :: shr_cal_date2eday  ! converts coded-date   to elapsed days
   public :: shr_cal_date2julian! converts coded-date,sec to julian days
   public :: shr_cal_ymd2julian ! converts yr,month,day,sec to julian days
   public :: shr_cal_ymd2date   ! converts yr,month,day to coded-date
   public :: shr_cal_ymd2eday   ! converts yr,month,day to elapsed days
   public :: shr_cal_advDate    ! advance date/secs real seconds
   public :: shr_cal_advDateInt ! advance date/secs integer seconds
   public :: shr_cal_validDate  ! logical function: is coded-date valid?
   public :: shr_cal_validYMD   ! logical function: are yr,month,day valid?
   public :: shr_cal_validHMS   ! logical function: are hr, min, sec valid?
   public :: shr_cal_getDebug   ! get internal debug level
   public :: shr_cal_setDebug   ! set internal debug level

! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   integer(SHR_KIND_IN),parameter,public :: shr_cal_calMaxLen = 64
   character(len=*),parameter,public :: &
      shr_cal_noleap    = 'NO_LEAP', &
      shr_cal_gregorian = 'GREGORIAN'

   !--- this is the noleap calendar ---
   integer(SHR_KIND_IN),parameter :: dpy = 365
   integer(SHR_KIND_IN),parameter :: dpm(12) = &     ! Days Per Month
   &                    (/31,28,31, 30, 31, 30,  31, 31, 30,  31, 30, 31/)


   !--- trigger internal debug output ---
   integer(SHR_KIND_IN) :: debug = 0

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_calendarName - check calendar name and translate
!
! !DESCRIPTION:
!    Check the validity of the calendar name and translate to standard naming
!
! !REVISION HISTORY:
!     2010-oct-7 - T Craig - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_cal_calendarName(calendar,trap)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)         ,intent(in) :: calendar  ! calendar type
   logical     ,optional,intent(in) :: trap
   character(len=shr_cal_calMaxLen) :: shr_cal_calendarName

!EOP

   character(len=shr_cal_calMaxLen) :: lcal
   logical :: ltrap
   character(len=shr_cal_calMaxLen) :: lowercalendar
   character(*),parameter :: subName = "(shr_cal_calendarName)"

   ltrap = .true.
   if (present(trap)) then
      ltrap = trap
   endif

   lcal = ' '
   lowercalendar = trim(shr_string_toLower(trim(calendar)))
   lcal = trim(calendar)

   selectcase(trim(lowercalendar))

   case ('noleap')
      lcal = trim(shr_cal_noleap)
   case ('no_leap')
      lcal = trim(shr_cal_noleap)
   case ('365_day')
      lcal = trim(shr_cal_noleap)
   case ('365day')
      lcal = trim(shr_cal_noleap)
   case ('gregorian')
      lcal = trim(shr_cal_gregorian)
   case ('standard')
      lcal = trim(shr_cal_gregorian)
   case ('proleptic_gregorian')
      lcal = trim(shr_cal_gregorian)
   case default
      if (ltrap) then
         write(s_logunit,*) trim(subname),' : ERROR calendar not supported : ',trim(calendar)
         call shr_sys_abort(trim(subname)//' : calendar not supported : '//trim(calendar))
      endif
   end select

   shr_cal_calendarName = trim(lcal)

   end function shr_cal_calendarName

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_numDaysInMonth - return the number of days in a month.
!
! !DESCRIPTION:
!    Deturn the number of days in a month.
!
! !REVISION HISTORY:
!     2002-sep-18 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

integer function shr_cal_numDaysInMonth(year,month,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year,month  ! calendar year,month
   character(*)        ,intent(in)  :: calendar  ! calendar type

!EOP

   type(ESMF_time) :: ltime1,ltime2
   integer(SHR_KIND_IN) :: eday1,eday2
   character(len=shr_cal_calMaxLen) :: lcalendar
   integer :: rc
   character(*),parameter :: subName = "(shr_cal_numDaysInMonth)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lcalendar = shr_cal_calendarName(calendar)

   if (trim(lcalendar) == trim(shr_cal_gregorian)) then
#ifdef USE_ESMF_LIB 
      call ESMF_TimeSet(ltime1,yy=year,mm=month,dd=1,s=0,calendarType=ESMF_CAL_GREGORIAN,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      if (month < 12) then
         call ESMF_TimeSet(ltime2,yy=year,mm=month+1,dd=1,s=0,calendarType=ESMF_CAL_GREGORIAN,rc=rc)
         if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      else
         call ESMF_TimeSet(ltime2,yy=year+1,mm=1,dd=1,s=0,calendarType=ESMF_CAL_GREGORIAN,rc=rc)
         if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      endif
      call ESMF_TimeGet(ltime1,d=eday1,rc=rc)
      call ESMF_TimeGet(ltime2,d=eday2,rc=rc)
      shr_cal_numDaysInMonth = eday2-eday1
#else
      call shr_sys_abort(trim(subname)//' ERROR: shr_cal gregorian requires USE_ESMF_LIB')
#endif
   else
      shr_cal_numDaysInMonth = dpm(month)
   endif

end function shr_cal_numDaysInMonth

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_numDaysInYear - return the number of days in a year.
!
! !DESCRIPTION:
!    Deturn the number of days in a year.
!
! !REVISION HISTORY:
!     2002-sep-18 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

integer function shr_cal_numDaysInYear(year,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year      ! calendar year
   character(*)        ,intent(in)  :: calendar  ! calendar type

!EOP

   type(ESMF_time) :: ltime1,ltime2
   integer(SHR_KIND_IN) :: eday1,eday2
   character(len=shr_cal_calMaxLen) :: lcalendar
   integer :: rc
   character(*),parameter :: subName = "(shr_cal_numDaysInYear)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lcalendar = shr_cal_calendarName(calendar)

   if (trim(lcalendar) == trim(shr_cal_gregorian)) then
#ifdef USE_ESMF_LIB 
      call ESMF_TimeSet(ltime1,yy=year,mm=1,dd=1,s=0,calendarType=ESMF_CAL_GREGORIAN,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_TimeSet(ltime2,yy=year+1,mm=1,dd=1,s=0,calendarType=ESMF_CAL_GREGORIAN,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_TimeGet(ltime1,d=eday1,rc=rc)
      call ESMF_TimeGet(ltime2,d=eday2,rc=rc)
      shr_cal_numDaysInYear = eday2-eday1
#else
      call shr_sys_abort(trim(subname)//' ERROR: shr_cal gregorian requires USE_ESMF_LIB')
#endif
   else
      shr_cal_numDaysInYear = dpy
   endif

end function shr_cal_numDaysInYear

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_elapsDaysStrtMonth - return the number of elapsed days
!            at start of month
!
! !DESCRIPTION:
!    Return the number of elapsed days at start of a month.
!
! !REVISION HISTORY:
!     2002-Oct-29 - R. Jacob - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

integer function shr_cal_elapsDaysStrtMonth(year,month,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year,month  ! calendar year,month
   character(*)        ,intent(in)  :: calendar  ! calendar type

!EOP

   integer :: rc
   integer :: k
   character(*),parameter :: subName = "(shr_cal_elapsDaysStrtMonth)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   shr_cal_elapsDaysStrtMonth = 0
   do k = 1,month-1
      shr_cal_elapsDaysStrtMonth = shr_cal_elapsDaysStrtMonth + &
                                   shr_cal_numDaysInMonth(year,k,calendar)
   enddo

end function shr_cal_elapsDaysStrtMonth

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_eday2date - converts elapsed days to coded-date
!
! !DESCRIPTION:
!     Converts elapsed days to coded-date.
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_eday2date(eday,date,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in)  :: eday  ! number of elapsed days
   integer(SHR_KIND_IN),intent(out) :: date  ! coded (yyyymmdd) calendar date
   character(*)        ,intent(in)  :: calendar  ! calendar type

!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: yy,mm,dd
   character(*),parameter :: subName = "(shr_cal_eday2date)"

!-------------------------------------------------------------------------------
! ASSUMPTIONS:
!   this calendar has a year zero (but no day or month zero)
!-------------------------------------------------------------------------------

   if (debug > 1) write(s_logunit,*) trim(subname),'_a ',eday

   call shr_cal_eday2ymd(eday,yy,mm,dd,calendar)
   call shr_cal_ymd2date(yy,mm,dd,date)

   if (debug > 1) write(s_logunit,*) trim(subname),'_b ',date

end subroutine shr_cal_eday2date

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_eday2ymd - converts elapsed days to year/month/day.
!
! !DESCRIPTION:
!     Converts elapsed days to year/month/day.
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_eday2ymd (eday,year,month,day,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in)  :: eday             ! elapsed days
   integer(SHR_KIND_IN),intent(out) :: year,month,day   ! calendar year,month,day
   character(*),        intent(in)  :: calendar         ! calendar type

!EOP

   !--- local ---
   type(ESMF_Time) :: ltime
   integer(SHR_KIND_IN) :: k,tday,rc
   character(len=shr_cal_calMaxLen) :: lcalendar
   character(*),parameter :: subName = "(shr_cal_eday2ymd)"

!-------------------------------------------------------------------------------
! ASSUMPTIONS:
!   this calendar has a year zero (but no day or month zero)
!-------------------------------------------------------------------------------

   if (debug > 1) write(s_logunit,*) trim(subname),'_a ',eday

   lcalendar = shr_cal_calendarName(calendar)

   if (trim(lcalendar) == trim(shr_cal_gregorian)) then
#ifdef USE_ESMF_LIB 
      call ESMF_TimeSet(ltime,d=eday,s=0,calendarType=ESMF_CAL_GREGORIAN,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_TimeGet(ltime,yy=year,mm=month,dd=day,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
#else
      call shr_sys_abort(trim(subname)//' ERROR: shr_cal gregorian requires USE_ESMF_LIB')
#endif
   else

      if (eday < 0) then
         year = -( abs(eday+1)/dpy + 1)
      else
         year = eday/dpy       ! calendar year (note: Fortran truncation)
      endif

      tday = eday - year*dpy

      if (tday < 0 .or. tday > dpy-1) then
         write(s_logunit,*) trim(subname),' ERROR: tday ',eday,tday,year,dpy
         call shr_sys_abort(trim(subname)//' ERROR: tday error')
      endif

      day  = mod(tday,dpy) + 1  ! elapsed days within current year
      month = -1
      k = 0
      do while (month < 1 .and. k < 12)
         k = k + 1
         IF (day .le. dpm(k)) then
           month = k                    ! calendar month
         else
           day = day - dpm(k)           ! calendar day
         endif
      end do   ! while

   endif

   if (.not. shr_cal_validYMD(year,month,day,calendar)) then
      if (s_loglev > 0) write(s_logunit,*) trim(subname)," ERROR: invalid date = ", &
         eday,year,month,day
   endif

   if (debug > 1) write(s_logunit,*) trim(subname),'_b ',year,month,day

end subroutine shr_cal_eday2ymd 

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_date2ymd - converts coded-date to year/month/day.
!
! !DESCRIPTION:
!     Converts coded-date (yyyymmdd) to year/month/day.
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_date2ymd (date,year,month,day)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in)  :: date             ! coded-date (yyyymmdd)
   integer(SHR_KIND_IN),intent(out) :: year,month,day   ! calendar year,month,day

!EOP

   integer(SHR_KIND_IN) :: tdate   ! temporary date
   character(*),parameter :: subName = "(shr_cal_date2ymd)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (debug > 1) write(s_logunit,*) trim(subname),'_a ',date

   tdate = abs(date)
   year =int(     tdate       /10000)
   if (date < 0) year = -year
   month=int( mod(tdate,10000)/  100)
   day  =     mod(tdate,  100) 

   if (debug > 1) write(s_logunit,*) trim(subname),'_b ',year,month,day

end subroutine shr_cal_date2ymd 

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_date2eday - converts coded-date to elapsed days
!
! !DESCRIPTION:
!     Converts coded-date to elapsed days
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_date2eday(date,eday,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: date            ! coded (yyyymmdd) calendar date
   integer(SHR_KIND_IN),intent(out) :: eday            ! number of elapsed days
   character(*)        ,intent(in)  :: calendar        ! calendar type

!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: year,month,day
   character(*),parameter :: subName = "(shr_cal_date2eday)"

!-------------------------------------------------------------------------------
! NOTE:
!   elapsed days since yy-mm-dd = 00-01-01, with 0 elapsed seconds
!-------------------------------------------------------------------------------

   if (debug > 1) write(s_logunit,*) trim(subname),'_a ',date

   call shr_cal_date2ymd(date,year,month,day)
   call shr_cal_ymd2eday(year,month,day,eday,calendar)

   if (debug > 1) write(s_logunit,*) trim(subname),'_b ',eday

end subroutine shr_cal_date2eday

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_date2julian - converts coded-date to julian day of year
!
! !DESCRIPTION:
!     Converts coded-date to julian day of year
!
! !REVISION HISTORY:
!     2009-Oct-23 - T. Craig - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_date2julian(date,sec,jday,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: date            ! coded (yyyymmdd) calendar date
   integer(SHR_KIND_IN),intent(in ) :: sec             ! seconds
   real   (SHR_KIND_R8),intent(out) :: jday            ! julian day of year
   character(*)        ,intent(in)  :: calendar        ! calendar type

!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: year,month,day
   character(*),parameter :: subName = "(shr_cal_date2julian)"

!-------------------------------------------------------------------------------
! NOTE:
!   julian day of year since yy-01-01-00000
!-------------------------------------------------------------------------------

   if (debug > 1) write(s_logunit,*) trim(subname),'_a ',date,sec

   call shr_cal_date2ymd(date,year,month,day)
   call shr_cal_ymd2julian(year,month,day,sec,jday,calendar)

   if (debug > 1) write(s_logunit,*) trim(subname),'_b ',jday

end subroutine shr_cal_date2julian

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_ymd2julian - converts y,m,d,s to julian day of year
!
! !DESCRIPTION:
!     Converts y,m,d,s to julian day of year
!
! !REVISION HISTORY:
!     2009-Oct-23 - T. Craig - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_ymd2julian(year,month,day,sec,jday,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year            ! year
   integer(SHR_KIND_IN),intent(in ) :: month           ! month
   integer(SHR_KIND_IN),intent(in ) :: day             ! day
   integer(SHR_KIND_IN),intent(in ) :: sec             ! seconds
   real   (SHR_KIND_R8),intent(out) :: jday            ! julian day of year
   character(*)        ,intent(in)  :: calendar        ! calendar type

!EOP

   !--- local ---
   character(*),parameter :: subName = "(shr_cal_ymd2julian)"

!-------------------------------------------------------------------------------
! NOTE:
!   julian day of year since yy-01-01-000000
!-------------------------------------------------------------------------------

   if (debug > 1) write(s_logunit,*) trim(subname),'_a ',year,month,day,sec

   if (.not. shr_cal_validYMD(year,month,day,calendar)) then
      write(s_logunit,*) trim(subname),' ERROR: invalid ymd',year,month,day
      call shr_sys_abort(trim(subname)//' ERROR: invalid ymd')
   endif

   jday = shr_cal_elapsDaysStrtMonth(year,month,calendar) + day + sec/SHR_CONST_CDAY

   if (debug > 1) write(s_logunit,*) trim(subname),'_b ',jday

end subroutine shr_cal_ymd2julian

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_ymd2date - converts year, month, day to coded-date
!
! !DESCRIPTION:
!     Converts  year, month, day to coded-date
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_ymd2date(year,month,day,date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year,month,day  ! calendar year,month,day
   integer(SHR_KIND_IN),intent(out) :: date            ! coded (yyyymmdd) calendar date

!EOP

   !--- local ---
   character(*),parameter :: subName = "(shr_cal_ymd2date)"

!-------------------------------------------------------------------------------
! NOTE:
!   this calendar has a year zero (but no day or month zero)
!-------------------------------------------------------------------------------

   if (debug > 1) write(s_logunit,*) trim(subname),'_a ',year,month,day

   date = abs(year)*10000 + month*100 + day  ! coded calendar date
   if (year < 0) date = -date

   if (debug > 1) write(s_logunit,*) trim(subname),'_b ',date

end subroutine shr_cal_ymd2date

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_ymd2eday - converts year, month, day to elapsed days
!
! !DESCRIPTION:
!     Converts  year, month, day to elapsed days
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine  shr_cal_ymd2eday(year,month,day,eday,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year,month,day  ! calendar year,month,day
   integer(SHR_KIND_IN),intent(out) :: eday            ! number of elapsed days
   character(*)        ,intent(in)  :: calendar        ! calendar type

!EOP

   !--- local ---
   type(ESMF_Time) :: ltime
   integer(SHR_KIND_IN) :: eday0
   character(len=shr_cal_calMaxLen) :: lcalendar
   integer :: rc
   character(*),parameter :: subName = "(shr_cal_ymd2eday)"

!-------------------------------------------------------------------------------
! NOTE:
!   elapsed days since yy-mm-dd = 00-01-01, with 0 elapsed seconds
!-------------------------------------------------------------------------------

   if (debug > 1) write(s_logunit,*) trim(subname),'_a ',year,month,day

   lcalendar = shr_cal_calendarName(calendar)

   if (.not. shr_cal_validYMD(year,month,day,calendar)) then
      write(s_logunit,*) trim(subname),' ERROR: invalid ymd',year,month,day,trim(calendar)
      call shr_sys_abort(trim(subname)//' ERROR: invalid ymd')
   endif

   if (trim(lcalendar) == trim(shr_cal_gregorian)) then
#ifdef USE_ESMF_LIB 
      call ESMF_TimeSet(ltime,yy=year,mm=month,dd=day,s=0,calendarType=ESMF_CAL_GREGORIAN,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_TimeGet(ltime,d=eday,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
#else
      call shr_sys_abort(trim(subname)//' ERROR: shr_cal gregorian requires USE_ESMF_LIB')
#endif
   else
      eday = year*dpy + shr_cal_elapsDaysStrtMonth(year,month,calendar) + day - 1
   endif

   if (debug > 1) write(s_logunit,*) trim(subname),'_b ',eday

end subroutine  shr_cal_ymd2eday

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_advDate - advances a date and seconds with a delta time
!
! !DESCRIPTION:
!     Advances a date and seconds with a delta time
!
! !REVISION HISTORY:
!     2009-Jun-09 - T. Craig - allows delta < 0
!     2005-Jun-10 - B. Kauffman - bug fix, simplified algorithm
!     2005-May-15 - T. Craig - initial version 
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_advDate(delta,units,dateIN,secIN,dateOUT,secOUT,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   real   (SHR_KIND_R8) ,intent(in)  :: delta     ! time increment
   character(*)         ,intent(in)  :: units     ! units of increment
   integer(SHR_KIND_IN) ,intent(in)  :: dateIN    ! base date, yyyymmdd
   real   (SHR_KIND_R8) ,intent(in)  :: secIN     ! base seconds
   integer(SHR_KIND_IN) ,intent(out) :: dateOUT   ! new date, yyyymmdd
   real   (SHR_KIND_R8) ,intent(out) :: secOUT    ! new seconds
   character(*)         ,intent(in)  :: calendar  ! calendar type

!EOP

   !--- local ---
   real   (SHR_KIND_R8)   :: dSec    ! delta-sec: advance date this many seconds
   integer(SHR_KIND_IN)   :: dDay    ! advance this many whole days
   real   (SHR_KIND_R8)   :: rSec    ! advance this many "remainder seconds"
   integer(SHR_KIND_IN)   :: eDay    ! date as elapsed dates from ref date
   integer(SHR_KIND_IN)   :: dayadjust ! time adjustment to support negative delta
   integer(SHR_KIND_IN)   :: n       ! counter
   logical                :: found   ! check for valid type

   !--- formats ---
   character(*),parameter :: subName = "(shr_cal_advDate)"
   character(*),parameter :: F00 = "('(shr_cal_advDate) ',a,i5)"
   character(*),parameter :: F02 = "('(shr_cal_advDate) ',a,i8.8,f10.3)"
   
!-------------------------------------------------------------------------------
! NOTE:
!-------------------------------------------------------------------------------

   !--- calculate delta-time in seconds ---
   if     (trim(units) == 'days'   ) then
      dSec = delta * SHR_CONST_CDAY
   elseif (trim(units) == 'hours'  ) then
      dSec = delta * 3600.0_SHR_KIND_R8
   elseif (trim(units) == 'minutes') then
      dSec = delta *   60.0_SHR_KIND_R8
   elseif (trim(units) == 'seconds') then
      dSec = delta *    1.0_SHR_KIND_R8
   else
      call shr_sys_abort(trim(subname)//' ERROR: unrecognized time units '//trim(units))
   endif

   !--- starting from dateIN and zero seconds...      ---
   !--- advance (secIN + dSec) seconds to arrive at output date ---

   call shr_cal_date2eday(dateIN,eDay,calendar) ! starting from this eDay...
   dSec = dSec + secIN                          ! advance this many seconds

   !--- if dSec lt 0 then add dayadjust to dSec to make it positive
   !--- and subtract dayadjust from dDay to correct it.
   if (dSec < 0.0_SHR_KIND_R8) then
      dayadjust = int(abs(dSec)/SHR_CONST_CDAY) + 1
      dSec = dSec + dayadjust*SHR_CONST_CDAY
      dDay = int(dSec/SHR_CONST_CDAY) - dayadjust
      rSec = mod(dSec,SHR_CONST_CDAY)
   else
      dDay = int(dSec/SHR_CONST_CDAY)     ! advance this many whole days...
      rSec = mod(dSec,SHR_CONST_CDAY)     ! ...plus this many secs (less than a day)
   endif

   call shr_cal_eday2date(eDay+dDay,dateOUT,calendar)
   secOUT = rSec    

   if (debug>0) then
      if (s_loglev > 0) write(s_logunit,*) subName," units,delta,calendar=",trim(units),delta,' ',trim(calendar)
      if (s_loglev > 0) write(s_logunit,F02) "dateIN ,secIN =",dateIN ,secIN
      if (s_loglev > 0) write(s_logunit,F02) "dateOUT,secOUT=",dateOUT,secOUT
   end if

end subroutine shr_cal_advDate

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_advDateInt - advances a date and seconds with a delta time
!
! !DESCRIPTION:
!     Advances a date and seconds with a delta time
!
! !REVISION HISTORY:
!     2009-???-?? - ?? - replicated from shr_cal_advDate()
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_advDateInt(delta,units,dateIN,secIN,dateOUT,secOUT,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN) ,intent(in)  :: delta     ! time increment
   character(*)         ,intent(in)  :: units     ! units of increment
   integer(SHR_KIND_IN) ,intent(in)  :: dateIN    ! base date, yyyymmdd
   integer(SHR_KIND_IN) ,intent(in)  :: secIN     ! base seconds
   integer(SHR_KIND_IN) ,intent(out) :: dateOUT   ! new date, yyyymmdd
   integer(SHR_KIND_IN) ,intent(out) :: secOUT    ! new seconds
   character(*)         ,intent(in)  :: calendar  ! calendar type

!EOP

   !--- local ---
   integer(SHR_KIND_I8)   :: dSec    ! delta-sec: advance date this many seconds
   integer(SHR_KIND_IN)   :: dDay    ! advance this many whole days
   integer(SHR_KIND_IN)   :: eDay    ! date as elapsed dates from ref date
   integer(SHR_KIND_IN)   :: dayadjust ! time adjustment to support negative delta
   integer(SHR_KIND_IN)   :: n       ! counter
   integer(SHR_KIND_I8)   :: spd     ! seconds per day
   logical                :: found   ! check for valid type

   !--- formats ---
   character(*),parameter :: subName = "(shr_cal_advDateInt)"
   character(*),parameter :: F00 = "('(shr_cal_advDateInt) ',a,i5)"
   character(*),parameter :: F02 = "('(shr_cal_advDateInt) ',a,i8.8,f10.3)"
   
!-------------------------------------------------------------------------------
! NOTE:
!-------------------------------------------------------------------------------

   spd = nint(SHR_CONST_CDAY)

   !--- calculate delta-time in seconds ---
   if     (trim(units) == 'days'   ) then
      dSec = delta * spd
   elseif (trim(units) == 'hours'  ) then
      dSec = delta * 3600
   elseif (trim(units) == 'minutes') then
      dSec = delta *   60
   elseif (trim(units) == 'seconds') then
      dSec = delta *    1
   else
      call shr_sys_abort(trim(subname)//' ERROR: unrecognized time units '//trim(units))
   endif

   !--- starting from dateIN and zero seconds...      ---
   !--- advance (secIN + dSec) seconds to arrive at output date ---

   call shr_cal_date2eday(dateIN,eDay,calendar) ! starting from this eDay...
   dSec = dSec + secIN                          ! advance this many seconds

   !--- if dSec lt 0 then add dayadjust to dSec to make it positive
   !--- and subtract dayadjust from dDay to correct it.
   if (dSec < 0 ) then
      dayadjust = (abs(dSec)/spd) + 1
      dSec = dSec + dayadjust*spd
      dDay =    (dSec/spd) - dayadjust
      secOUT = mod(dSec,spd)
   else
      dDay =    (dSec/spd)     ! advance this many whole days...
      secOUT = mod(dSec,spd)     ! ...plus this many secs (less than a day)
   endif

   call shr_cal_eday2date(eDay+dDay,dateOUT,calendar)

   if (debug>0) then
      if (s_loglev > 0) write(s_logunit,*) subName," units,delta,calendar=",trim(units),delta,' ',trim(calendar)
      if (s_loglev > 0) write(s_logunit,F02) "dateIN ,secIN =",dateIN ,secIN
      if (s_loglev > 0) write(s_logunit,F02) "dateOUT,secOUT=",dateOUT,secOUT
   end if

end subroutine shr_cal_advDateInt

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_validDate - determines if coded-date is a valid date
!
! !DESCRIPTION:
!    Determines if the given coded-date is a valid date.
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

logical function shr_cal_validDate(date,calendar) 

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: date      ! coded (yyyymmdd) calendar date
   character(len=*)    ,intent(in ) :: calendar  ! calendar

!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: year,month,day
   integer(SHR_KIND_IN) :: tdate
   character(*),parameter :: subName = "(shr_cal_validDate)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call shr_cal_date2ymd(date,year,month,day)
   shr_cal_validDate = shr_cal_validYMD(year,month,day,calendar)

end function shr_cal_validDate

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_validYMD - determines if year, month, day is a valid date
!
! !DESCRIPTION:
!    Determines if the given year, month, and day indicate a valid date.
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

logical function shr_cal_validYMD(year,month,day,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year,month,day  ! year,month,day
   character(len=*)    ,intent(in ) :: calendar        ! calendar

!EOP

   !--- local ---

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   shr_cal_validYMD = .true.
   if (year  < -999) shr_cal_validYMD = .false.
   if (year  > 9999) shr_cal_validYMD = .false.
   if (month <    1) shr_cal_validYMD = .false.
   if (month >   12) shr_cal_validYMD = .false.
   if (day   <    1) shr_cal_validYMD = .false.
   if (day   > shr_cal_numDaysInMonth(year,month,calendar)) &
                     shr_cal_validYMD = .false.

end function shr_cal_validYMD

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_validHMS - determines if hour, min, sec is valid
!
! !DESCRIPTION:
!    Determines if the given hour, min, sec is valid
!
! !REVISION HISTORY:
!     2005-May-15 - T. Craig - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

logical function shr_cal_validHMS(hr,min,sec)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: hr, min, sec   ! hour, minute, second

!EOP

   !--- local ---
   character(*),parameter :: subName = "(shr_cal_validHMS)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   shr_cal_validHMS = .true.
   if (hr  <    0) shr_cal_validHMS = .false.
   if (hr  >   23) shr_cal_validHMS = .false.
   if (min <    0) shr_cal_validHMS = .false.
   if (min >   59) shr_cal_validHMS = .false.
   if (sec <    0) shr_cal_validHMS = .false.
   if (sec >   60) shr_cal_validHMS = .false.

end function shr_cal_validHMS

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_setDebug -- Set local shr_cal debug level
!
! !DESCRIPTION:
!    Set local shr\_cal debug level, 0 = production
!    \newline
!    General Usage: call shr\_cal\_setDebug(2)
!
! !REVISION HISTORY:
!     2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_cal_setDebug(level)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(in) :: level

!EOP

  !--- formats ---
  character(*),parameter :: subName =   "(shr_cal_setDebug) "
  character(*),parameter :: F00     = "('(shr_cal_setDebug) ',a) "
  character(*),parameter :: F01     = "('(shr_cal_setDebug) ',a,i4) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  debug = level
  if (s_loglev > 0) write(s_logunit,F01) "debug level reset to ",level

end subroutine shr_cal_setDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_getDebug -- get shr_cal internal debug level
!
! !DESCRIPTION:
!    Get shr_cal internal debug level, 0 = production
!    \newline
!    General Usage: call shr\_cal\_getDebug(level)
!
! !REVISION HISTORY:
!     2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_cal_getDebug(level)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(out) :: level

!EOP

  !--- formats ---
  character(*),parameter :: subName =   "(shr_cal_getDebug) "
  character(*),parameter :: F00     = "('(shr_cal_getDebug) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  level = debug

end subroutine shr_cal_getDebug

!===============================================================================
!===============================================================================

end module shr_cal_mod
