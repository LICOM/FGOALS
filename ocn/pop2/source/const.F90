!  CVS: $Id: const.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ================
      SUBROUTINE CONST
!     ================
!-----------------------------------------------------------------------
!
! Purpose: Set up some constants and control parameter.
!
! Author: Yongqiang Yu and Hailong Liu, Dec. 31, 2002
!
!
!-----------------------------------------------------------------------

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use diag_mod
use pmix_mod
#ifdef SPMD
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
#endif
#ifdef COUP
use shr_sys_mod
!use shr_msg_mod by linpf
use shr_msg_mod !,only:shr_sys_flush,shr_msg_stdio
#endif


      IMPLICIT NONE
!lhl090729
#include <netcdf.inc>
      integer :: ncid,iret
      real(r8) :: tmpy(s_jmt)
!lhl090729
      INTEGER :: IDTB,IDTC,IDTS
      REAL(r8)    :: AG,ALFA

      namelist /namctl/ AFB1,AFC1,AFT1,IDTB,IDTC,IDTS,AMV,AHV,NUMBER, &
                        NSTART,IO_HIST,IO_REST,klv,DLAM,AM_TRO,AM_EXT,&
                        diag_msf,diag_bsf,diag_budget,diag_mth,rest_freq,hist_freq,out_dir
!-------------------------------------------------------
!     Set up the default value for namelist.
!-------------------------------------------------------
!
!      diag_msf=.false.
!      diag_bsf=.false.
!      diag_budget=.false.
!      diag_mth=.false.
!
      IDTB=30   ; IDTC=1800 ; IDTS=7200
      AFB1=0.025D0; AFC1=0.43D0 ; AFT1=0.43D0
!-------------------------------------------------------
!     Set up the constants for calendar month.
!-------------------------------------------------------
      NMONTH=RESHAPE((/31,28,31,30,31,30,31,31,30,31,30,31/),(/12/))
     NNMONTH=RESHAPE((/15,46,74,105,135,166,196,227,258,288,319,349/),(/12/))
      ABMON =RESHAPE((/'Jan','Feb','Mar','Apr','May','Jun', &
                   'Jul','Aug','Sep','Oct','Nov','Dec'/),(/12/))
      ABMON1=RESHAPE((/'jan','feb','mar','apr','may','jun', &
                    'jul','aug','sep','oct','nov','dec'/),(/12/))
!-------------------------------------------------------
!     PHYSICAL CONSTANTS
!-------------------------------------------------------
!     G     Acceleration of gravity in m/sec**2
!     CP    Specific heat capacity of sea water in J/kg/K
!     D0    Density of sea water in kg/m**3
!     AG    Ekman bias angle
!     C0F    friction coefficient
!     TBICE the frozen point of seawater in C
!     KARMAN the Karman number for Smagrinsky horizontal voscosity
!            diffusion

      G = 9.806D0

      CP = 3996.0D0
      D0 = 1026.0D0
      AG = 10.0D0*3.1415926D0/180.0D0
      C0F = 2.6D-3
      TBICE = -1.8D0
#if ( defined SMAG)
!      KARMAN= 0.6D0
      KARMAN = 0.14D0
      RR = 0.5D0
#endif

      SAG = SIN (AG)
      CAG = COS (AG)
      OD0 = 1.0D0/ D0
      OD0CP = 1.0D0/ (D0* CP)

!-------------------------------------------------------
!     DIFFUSION & VISCOSITY
!-------------------------------------------------------
!     AM    laternal viscosity coeffcient in m**2/sec
!     AH    laternal diffusion coeffcient in m**2/sec
!     AMV   vertical viscosity coeffcient in m**2/sec
!     AHV   vertical diffusion coeffcient in m**2/sec
!     AHICE diffusion coeffcient between ice & water

!move AM  = 2.0E+5
!YU
!     AM_TRO  = 2.0E+3
!     AM_EXT  = 2.0E+5
!     DLAM    = 1.0
!YU
      AMV = 1.0D-3
!lhl  AH  = 2.0D+3
      AHV = 0.3D-4

!
#ifdef SPMD
!-------------------------------------------------------
!     Read namelist from the external file.
!-------------------------------------------------------
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.0"
      close(111)
      end if
      if (mytid==0) then
	      open(11,file='ocn.parm',form='formatted',status='OLD')
	      rewind 11
	      read(11,namctl)
	      close(11)
              write (112,namctl)
              close(112)
      end if
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.1", mpi_comm_ocn
      close(111)
      end if
!     call mpi_barrier(mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.2"
      close(111)
      end if
      call mpi_bcast(afb1,1,MPI_PR,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.3"
      close(111)
      end if
      call mpi_bcast(afc1,1,MPI_PR,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.4"
      close(111)
      end if
      call mpi_bcast(aft1,1,MPI_PR,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.5"
      close(111)
      end if
      call mpi_bcast(amv,1,MPI_PR,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.6"
      close(111)
      end if
      call mpi_bcast(ahv,1,MPI_PR,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.7"
      close(111)
      end if
      call mpi_bcast(idtb,1,mpi_integer,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.8"
      close(111)
      end if
      call mpi_bcast(idtc,1,mpi_integer,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.9"
      close(111)
      end if
      call mpi_bcast(idts,1,mpi_integer,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.10"
      close(111)
      end if
      call mpi_bcast(number,1,mpi_integer,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.11"
      close(111)
      end if
      call mpi_bcast(nstart,1,mpi_integer,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.12"
      close(111)
      end if
      call mpi_bcast(io_hist,1,mpi_integer,0,mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.13"
      close(111)
      end if
      call mpi_bcast(io_rest,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(klv,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(dlam,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(am_tro,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(am_ext,1,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(hist_freq,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(rest_freq,1,mpi_integer,0,mpi_comm_ocn,ierr)
!      call mpi_barrier(mpi_comm_ocn,ierr)
      if (mytid ==0 ) then
      write(111,*) "OK----2    0.14"
      close(111)
      end if
#else
      open(11,file='ocn.parm',form='formatted')
      read(11,namctl)
#endif


      AHICE = AHV

!-------------------------------------------------------
!     wndmix = min value for mixing at surface to simulate
!              high freq wind mixing. (m**2/sec)
!     fricmx = maximum mixing (m**2/sec)
!     diff_cbt_back = background "diff_cbt"(m**2/sec)
!     visc_cbu_back = background"visc_cbu" (m**2/sec)
!     diff_cbt_limit = largest "diff_cbt"  (m**2/sec)
!     visc_cbu_limit = largest "visc_cbu"  (m**2/sec)
!-------------------------------------------------------

      wndmix = 10.0d-4
      fricmx = 50.0d-4
      diff_cbt_back = 1.0d-4
      visc_cbu_back = 1.0d-4

      visc_cbu_limit = fricmx
      diff_cbt_limit = fricmx

!lhl1204
! 1 in 2x2 model, baroclinic time step is 7200s at first and 9600s later,
!   and annual mean forcina, and from the begining
!      dfricmx = 1000.0d-4  !OK
!   more larger parameter is tested!
!   the time step of baroclinic and thermal have been reduced
!   to 1/4 (3600) and 1/2 (14400) respectively.
!   it is same as when the viscosity decreased.
!      dfricmx = 2000.0d-4  !OK
! 2 this parameter can also use when the viscosity decreased to 1e+4
!   but the time step of baroclinic and thermal have been reduced to
!   1/4 (3600) and 1/2 (14400) respectively. the former is close to the
!   time step which used in Canoto's 2001 paper.
!      dfricmx = 1000.0d-4  !OK
! in 2x2 model, baroclinic time step is 7200s, seasonal forcing, from the first experiment
!      dfricmx = 500.0d-4  !OK
       dfricmx = 100.0d-4 !OK
!       dfricmx = 2000.0d-4 !OK
!      dfricmx = 5000.0d-4 !OK
      dwndmix =   10.0d-4
!lhl1204


      DTB = FLOAT (IDTB)
      DTC = FLOAT (IDTC)
      DTS = FLOAT (IDTS)

      DTB2 = 2.0D0* DTB
      DTC2 = 2.0D0* DTC

!-------------------------------------------------------
!     NSS   steps for thermohaline procoss per day
!     NCC   baroclinic steps within one thermohaline step
!     NBB   barotropic steps within one baroclinic step
!-------------------------------------------------------
      NSS = 24*3600/ IDTS
#if ( defined SYNCH)
      NCC = IDTS / IDTC
#else
      NCC = 1
#endif
      NBB = IDTC / IDTB

      ONBB = 1.0D0/ FLOAT (NBB +1)
      ONCC = 1.0D0/ FLOAT (NCC +1)
      ONBC = 1.0D0/ FLOAT (NBB * NCC +1)

      AFB2 = 1.0D0-2.0D0* AFB1
      AFC2 = 1.0D0-2.0D0* AFC1
      AFT2 = 1.0D0-2.0D0* AFT1


!-------------------------------------------------------
!     RESTORING TIME SCALE FOR SALINITY
!-------------------------------------------------------
      GAMMA = 1.0D0/ (30.0D0*86400.0D0)
!YU
!-------------------------------------------------------
!     FOR FOURIOUR FILTERING
!-------------------------------------------------------

!      ALFA = 3.1416D0*2.0D0/ FLOAT (IMT -2)
      ALFA =  4.0* ATAN (1.0)*2.0D0/ FLOAT (IMT -2)
!$OMP PARALLEL DO PRIVATE (I)
      DO I = 1,IMT
         CF1 (I)= COS (ALFA * FLOAT (I -1))
         CF2 (I)= COS (2.0D0* ALFA * FLOAT (I -1))
         SF1 (I)= SIN (ALFA * FLOAT (I -1))
         SF2 (I)= SIN (2.0D0* ALFA * FLOAT (I -1))
      END DO


!-------------------------------------------------------
!     Read reference temperature, salinity and
!     coefficients for UNESCO formula.
!-------------------------------------------------------
      if (mytid==0) then
      open(33,file='dncoef.h1',form='formatted')
      read(33,*)
      read(33,*)to
      read(33,*)
      write(113,*) to
      close(113)
      read(33,*)so
      write(114,*) so
      close(114)
!lhl1204   add a arry to read reference potential density
!
      read(33,*)
      read(33,*)po
      write(115,*) po
      close(115)
!      read(33,'(8(f7.4,1x))')po
!lhl1204
      do i=1,km
         read(33,*)
         read(33,*)(c(i,k),k=1,9)
      end do
      close(33)
      endif
#ifdef SPMD
!     call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_bcast(to,km,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(so,km,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(po,km,MPI_PR,0,mpi_comm_ocn,ierr)
      call mpi_bcast(c,km*9,MPI_PR,0,mpi_comm_ocn,ierr)
!     call mpi_barrier(mpi_comm_ocn,ierr)
#endif
!YU
!-------------------------------------------------------
!     the units of density departures is [g/cm**3] in
!     GFDL's MOM, and [kg/m**3] in present model
!-------------------------------------------------------

!$OMP PARALLEL DO PRIVATE (M,K)
      DO M = 1,9
         DO K = 1,KM
            C (K,M)= C (K,M)*1000.0D0
         END DO
      END DO

!YU
      if ((imt_global-2)/=nint(360./DLAM)) then

!lhl0711      if ((imt-2)>=nint(360./DLAM)) then
         write(6,*)"Error in DLAM ! (Const)"
         write(6,*)imt-2,nint(360./DLAM)
         stop
      end if
      if (mytid ==0 ) then
      write(6,*)"AM_EXT,AM_TRO",AM_EXT,AM_TRO
      endif
!YU
!lhl090729
#ifdef FRC_CORE
      if (mytid.eq.0) then
      iret=nf_open("t_10.db.clim.daymean.15JUNE2009.nc",nf_nowrite,ncid)
      call check_err (iret)
      iret=nf_get_var_double(ncid,3,s_lon)
      call check_err (iret)
      iret=nf_get_var_double(ncid,2,s_lat)
      call check_err (iret)
      iret=nf_close(ncid)
      call check_err (iret)
      do j=1,s_jmt
      tmpy(j)=s_lat(s_jmt-j+1)
      enddo
      do j=1,s_jmt
      s_lat(j)=tmpy(j)
      enddo
      endif

#ifdef SPMD
      call mpi_barrier(mpi_comm_ocn,ierr)
      call mpi_bcast(s_lon,s_imt,mpi_double_precision,0,mpi_comm_ocn,ierr)
      call mpi_bcast(s_lat,s_jmt,mpi_double_precision,0,mpi_comm_ocn,ierr)
      call mpi_barrier(mpi_comm_ocn,ierr)
#endif
#endif
!lhl090729
      RETURN
      END SUBROUTINE CONST



