!  CVS: $Id: rdriver.F90,v 1.7 2003/08/25 07:47:52 lhl Exp $
!     ======================
      SUBROUTINE RDRIVER
!     ======================
 
#include <def-undef.h>
use param_mod
use pconst_mod
use forc_mod
use work_mod
use dyn_mod, only : buffer
#ifdef SPMD
use msg_mod
#endif

#ifdef COUP
use shr_sys_mod
#endif 

      IMPLICIT NONE
#include <netcdf.inc>
!
!     Define Variables.
      integer*4   :: ncid1, iret
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
 
!      REAL    :: WCOE (JMT),ABC
 
      allocate(su3(imt,jmt,12),sv3(imt,jmt,12),psa3(imt,jmt,12),tsa3(imt,jmt,12),qar3(imt,jmt,12),uva3(imt,jmt,12))
      allocate(swv3(imt,jmt,12),cld3(imt,jmt,12),sss3(imt,jmt,12),sst3(imt,jmt,12),nswv3(imt,jmt,12),dqdt3(imt,jmt,12),chloro3(imt,jmt,12))
      allocate(seaice3(imt,jmt,12),runoff3(imt,jmt,12))
      allocate(wspd3(imt,jmt,12),wspdu3(imt,jmt,12),wspdv3(imt,jmt,12),lwv3(imt,jmt,12),rain3(imt,jmt,12),snow3(imt,jmt,12))
!
      allocate(buffer_real4(imt_global,jmt_global))
      allocate(buffer(imt_global,jmt_global))
!
!
      if (mytid==0)then
      write(6,*)"Beginning------RDRIVER! "
#ifdef COUP
      call shr_sys_flush(6)
#endif
      endif 

!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               WKA (I,J,K)= 0.0
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     SU3  : Sea surface zonal wind stress         (N/M**2)
!     SV3  : Sea surface meridional wind stress    (N/M**2)
!     PSA3 : Sea surface air pressure              (Pa)
!     SWV3 : Total net downward solar radiation    (W/M**2)
!    NSWV3 : None Solar flux                       (Wm-2)
!    DQDT3 : Dq/Dt                                 (WK-1m-2)
!     SST3 : Sea surface temperature               (Celsius)
!     SSS3 : Sea surface salinity                  (psu)
!   chloro3:chlorophll concentration               (mg m-3)
!-----------------------------------------------------------------------
 
!     READ FORCING FIELD
#ifdef SPMD
!
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      if(mytid==0)then
      iret=nf_open('MODEL.FRC',nf_nowrite,ncid1)
      call check_err (iret)
      end if
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      do k=1, 12
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=k ; count(4)=1

      if (mytid ==0 ) then
         iret=nf_get_vara_real(ncid1,   5,start,count,buffer_real4)
         call check_err (iret)
      end if
      call global_distribute_real(buffer_real4,swv3(1,1,k))
!
      if (mytid == 0) then
         iret=nf_get_vara_real(ncid1,   6,start,count,buffer_real4)
         call check_err (iret)
      end if
      call global_distribute_real(buffer_real4,nswv3(1,1,k))
!
      if (mytid == 0) then
          iret=nf_get_vara_real(ncid1,   7,start,count,buffer_real4)
          call check_err (iret)
      end if
      call global_distribute_real(buffer_real4,dqdt3(1,1,k))
!
      if (mytid == 0 ) then
          iret=nf_get_vara_real(ncid1,   8,start,count,buffer_real4)
          call check_err (iret)
      end if
      call global_distribute_real(buffer_real4,su3(1,1,k))
!
      if (mytid == 0 ) then
          iret=nf_get_vara_real(ncid1,   9,start,count,buffer_real4)
          call check_err (iret)
      end if
      call global_distribute_real(buffer_real4,sv3(1,1,k))
!
      if (mytid == 0) then
          iret=nf_get_vara_real(ncid1,  10,start,count,buffer_real4)
          call check_err (iret)
      end if
      call global_distribute_real(buffer_real4,sst3(1,1,k))
!
      if (mytid == 0 ) then
         iret=nf_get_vara_real(ncid1,  11,start,count,buffer_real4)
         call check_err (iret)
      end if
      call global_distribute_real(buffer_real4,sss3(1,1,k))
!
     end do 
      if (mytid == 0 ) then
         iret = nf_close (ncid1)
         call check_err (iret)
      end if

!===============================================
!input seaice
#ifdef FRC_CORE
      if (mytid == 0 ) then
      iret=nf_open('seaice.db.NSIDC.1979-2006.clim.monthmean.modelgrid.01x01.nc',nf_nowrite,ncid1)
      call check_err (iret)
      end if

      do k=1, 12
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=k ; count(4)=1
      if (mytid == 0 ) then
      iret=nf_get_vara_double(ncid1, 5,start,count,buffer)
      call check_err (iret)
      end if
      call global_distribute(buffer,seaice3(1,1,k))
      end do 

      if (mytid == 0 ) then
      iret = nf_close (ncid1)
      call check_err (iret)
      end if

!===============================================
!input runoff
      if (mytid == 0 ) then
      iret=nf_open('runoff.db.clim.modelgrid.01x01.nc',nf_nowrite,ncid1)
      call check_err (iret)
      end if

      do k=1, 1
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=k ; count(4)=1
      if (mytid == 0 ) then
      iret=nf_get_vara_double(ncid1, 4,start,count,buffer)
      call check_err (iret)
      end if
      call global_distribute(buffer,runoff3(1,1,k))
      end do 

      if (mytid == 0 ) then
      iret = nf_close (ncid1)
      call check_err (iret)
      end if
#endif


!===============================================
!input the chlorophyll concentration
!===============================================
#if (defined SOLARCHLORO)
      if (mytid == 0)  then
      iret=nf_open('MODEL_CHLFRC',nf_nowrite,ncid1)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,      5,start,count,chloro3_io)
!     swv3_io must be change
      call check_err (iret)
      iret = nf_close (ncid1)
      call check_err (iret)
      end if
#endif
!==============================================
!
!        write(*,'(i4,11f8.2)') mytid,((su3(i,j,1),i=190,200),j=10,20)
!
#else
!Yu
#if (!defined CDFIN)
      OPEN (90,FILE ='MODEL.FRC',STATUS ='OLD',FORM ='UNFORMATTED')
      READ (90) SWV3_io,NSWV3_io,DQDT3_io,SU3_io,SV3_io,SST3_io,SSS3_io
      CLOSE (90)
      OPEN (91,FILE ='MODEL_CHLFRC',STATUS ='OLD',FORM ='UNFORMATTED')
      READ (91) chloro3_io
      CLOSE (91)
!
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('MODEL.FRC',nf_nowrite,ncid1)
      call check_err (iret)

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=1 ; count(4)=12

      iret=nf_get_vara_real(ncid1,   5,start,count,swv3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,   6,start,count,nswv3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,   7,start,count,dqdt3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,   8,start,count,su3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,   9,start,count,sv3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,  10,start,count,sst3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,  11,start,count,sss3_io)
      call check_err (iret)
!
      iret = nf_close (ncid1)
      call check_err (iret)

 
#if (defined SOLARCHLORO)
      iret=nf_open('MODEL_CHLFRC',nf_nowrite,ncid1)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,      5,start,count,chloro3_io)
      call check_err (iret)
      iret = nf_close (ncid1)
      call check_err (iret)
#endif
!
#endif

      SWV3=SWV3_io
      NSWV3=NSWV3_io
      DQDT3=DQDT3_io
      SU3=SU3_io
      SV3=SV3_io
      SST3=SST3_io
      SSS3=SSS3_io
      chloro3=chloro3_io
#endif
!Yu
      if (mytid==0) then
#ifdef COUP
      call shr_sys_flush(6)
#endif
      end if
 
!-----------------------------------------------------------------------
!     land grids of the forcing fields assigned to 0 
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
                SWV3(I,J,K)= SWV3(I,J,K)*VIT(I,J,1)
               NSWV3(I,J,K)=NSWV3(I,J,K)*VIT(I,J,1)
               DQDT3(I,J,K)=DQDT3(I,J,K)*VIT(I,J,1)
                 SU3(I,J,K)=  SU3(I,J,K)*VIV(I,J,1)
                 SV3(I,J,K)=  SV3(I,J,K)*VIV(I,J,1)
                SST3(I,J,K)= SST3(I,J,K)*VIT(I,J,1)
                SSS3(I,J,K)= SSS3(I,J,K)*VIT(I,J,1)
                seaice3(I,J,K)= seaice3(I,J,K)*VIT(I,J,1)
                runoff3(I,J,K)= runoff3(I,J,K)*VIT(I,J,1)
#if (defined SOLARCHLORO)
        chloro3(I,J,K)= chloro3(I,J,K)*VIT(I,J,1) 
#endif
                END DO
           END DO
        END DO

!-----------------------------------------------------------------------
!     salinity = (psu-35)*0.001
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
               SSS3 (I,J,K) = (SSS3 (I,J,K) -35.0D0)*0.001D0
            END DO
         END DO
      END DO
 
!-----------------------------------------------------------------------
!     reverse VS (southward is positive)
!     notice: the former program VS is reversed during preparing forcing field
!-----------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE (K,J,I)
      DO K = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
               SV3 (I,J,K) = -SV3 (I,J,K)
            END DO
         END DO
      END DO
 
 
!-----------------------------------------------------------------------
!     CALCULATING THE ANNUAL MEAN FORCING FIELD 
!-----------------------------------------------------------------------
 
#if (defined FRC_ANN)
 
      DO M = 1,12
!$OMP PARALLEL DO PRIVATE (J,I)
         DO J = 1,JMT
            DO I = 1,IMT
               WKA (I,J,1)= WKA (I,J,1) + SU3 (I,J,M)
               WKA (I,J,2)= WKA (I,J,2) + SV3 (I,J,M)
               WKA (I,J,3)= WKA (I,J,3) + SSS3 (I,J,M)
               WKA (I,J,4)= WKA (I,J,4) + SWV3 (I,J,M)
               WKA (I,J,5)= WKA (I,J,5) + SST3 (I,J,M)
               WKA (I,J,6)= WKA (I,J,6) + NSWV3 (I,J,M)
               WKA (I,J,7)= WKA (I,J,7) + DQDT3 (I,J,M)
#if (defined SOLARCHLORO)
               WKA (I,J,8)= WKA (I,J,8) + chloro3 (I,J,M)
#endif
            END DO
         END DO
      END DO
 
!$OMP PARALLEL DO PRIVATE (M,J,I)
      DO M = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
               SU3 (I,J,M) = WKA (I,J,1)/12.0D0
               SV3 (I,J,M) = WKA (I,J,2)/12.0D0
              SSS3 (I,J,M) = WKA (I,J,3)/12.0D0
              SWV3 (I,J,M) = WKA (I,J,4)/12.0D0
              SST3 (I,J,M) = WKA (I,J,5)/12.0D0
             NSWV3 (I,J,M) = WKA (I,J,6)/12.0D0
             DQDT3 (I,J,M) = WKA (I,J,7)/12.0D0
#if (defined SOLARCHLORO)
             chloro3 (I,J,M) = WKA (I,J,8)/12.0D0
#endif
            END DO
         END DO
      END DO
 
#endif
!
      if (mytid==0)then
      write(6,*)"END-----------RDRIVER!"
#ifdef COUP
      call shr_sys_flush(6)
#endif
      endif 
 
    deallocate(buffer_real4,buffer)

      RETURN
      END SUBROUTINE RDRIVER
 
