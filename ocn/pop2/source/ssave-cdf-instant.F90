!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine SSAVEINS
#include <def-undef.h>
!     =================
!     output in NETcdf format
!     written by liu hai long 2001 jun
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
!LPF 20120815
use buf_mod, only:t_cpl,s_cpl,u_cpl,v_cpl,dhdx,dhdy,q      
!LPF 20120815

      character (len=18) :: fname
      integer :: klevel
!
         number_day = iday !+1 !LPF 20120816
         number_month = month
         if ( number_day  > imd ) then
             number_day = 1
             number_month = month + 1
         end if
         nwmf= iyfm

      if(mytid==0) write(*,*)'in instant,iday=,imd=,rest_freq',iday,imd,rest_freq
!     if ( 1 .or. mod(iday,rest_freq) == 0 .or. iday == imd) then ! lihuimin, FOR DEBUG
     if ( mod(iday,rest_freq) == 0 .or. iday == imd) then
!
       if ( mytid == 0 ) then
         fname(1:8)='fort.22.'
         fname(13:13)='-'
         fname(16:16)='-'
         write(fname(14:15),'(i2.2)')mon0
         write(fname(9:12),'(i4.4)')nwmf
         write(fname(17:18),'(i2.2)') number_day
!
     if(mytid==0) write(6,*)'in instant',iday,imd,rest_freq
         if ( number_day ==1 .and. mon0 < 12) then
             write(fname(14:15),'(i2.2)')mon0+1
         end if
!
         if ( number_day ==1 .and. mon0 == 12) then
             write(fname(14:15),'(i2.2)')mon0-11
             write(fname(9:12),'(i4.4)')nwmf+1
         end if

         open (17, file="rpointer.ocn", form='formatted')
         write(17,'(a18)') fname
         close(17)
          write(*,*)'fname=',fname
         open(22,file=trim(out_dir)//fname,form='unformatted')
!         write(6,*) "in SSAVEcdf, month=",mon0,"day=",number_day
       end if
#ifdef SPMD
!
          allocate(buffer(imt_global,jmt_global))
!         write(*,*) 'allocate sucess'

         call local_to_global_4d_double(h0,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
!         open(92,file='z0_99.dat',form='unformatted')
!         WRITE (92)buffer
!          close(92) 
         end if
!           write(*,*)'finish h0'

         do klevel=1,km
         call local_to_global_4d_double(u(1,1,klevel),buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         end do
!         write(*,*)'finish U'
         
         do klevel=1,km
         call local_to_global_4d_double(v(1,1,klevel),buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         end do
!         write(*,*)'finish V'
!
         do klevel=1,km
         call local_to_global_4d_double(at(1,1,klevel,1),buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         end do
!         write(*,*)'finish at1'
!
         do klevel=1,km
         call local_to_global_4d_double(at(1,1,klevel,2),buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         end do
!         write(*,*)'finish at2'
!lhl20110728
         do klevel=1,km
         call local_to_global_4d_double(ws(1,1,klevel),buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
         end do
!         write(*,*)'finish at2'
         
         call local_to_global_4d_double(su,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish su'

         call local_to_global_4d_double(sv,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if

!         write(*,*)'finish sv'

         call local_to_global_4d_double(swv,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish swv'
         call local_to_global_4d_double(lwv,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish lwv'
         call local_to_global_4d_double(sshf,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish sshf'
         call local_to_global_4d_double(lthf,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish lthf'
         call local_to_global_4d_double(fresh,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish fresh'
!lhl20110728
         if (mytid==0) then
             write(22) number_month, number_day
         end if
!
         if(mytid==0) write(*,*)'finish number'
!lhl20120731
#ifdef COUP
         call local_to_global_4d_double(t_cpl,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish t_cpl'
         call local_to_global_4d_double(s_cpl,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish s_cpl'
         call local_to_global_4d_double(u_cpl,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish u_cpl'
         call local_to_global_4d_double(v_cpl,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish v_cpl'
         call local_to_global_4d_double(dhdx,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish dhdx'
         call local_to_global_4d_double(dhdy,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish dhdy'
         call local_to_global_4d_double(q,buffer,1,1)
         if (mytid==0) then
         WRITE (22)buffer
         end if
!         write(*,*)'finish q'
#endif

      if (mytid==0) then
         close(22)
      end if
      if(mytid==0) write(*,*)'ok  fort.22'
!lhl20120731
!
         deallocate( buffer)
#else
        write(22) h0
!
        do k =1, km
           write(22) ((u(i,j,k),i=1,imt_global),j=1,jmt_global)
        end do
!
        do k =1, km
           write(22) ((v(i,j,k),i=1,imt_global),j=1,jmt_global)
        end do
!
        do k =1, km
           write(22) ((at(i,j,k,1),i=1,imt_global),j=1,jmt_global)
        end do
!
        do k =1, km
           write(22) ((at(i,j,k,2),i=1,imt_global),j=1,jmt_global)
        end do
!lhl20110728
        do k =1, km
           write(22) ((ws(i,j,k),i=1,imt_global),j=1,jmt_global)
        end do
!lhl20110728
!
        write(22) number_month, number_day
!lhl20120731
#ifdef COUP
        write(22) t_cpl
        write(22) s_cpl
        write(22) u_cpl
        write(22) v_cpl
        write(22) dhdx
        write(22) dhdy
        write(22) q
#endif
!lhl20120731
        close(22)
!
#endif
      end if
!

!$OMP PARALLEL DO PRIVATE (k,j,i)
      DO k = 1,km
         DO j = 1,jmt ! Dec. 5, 2002, Yongqiang Yu
            DO i = 1,imt
               up (i,j,k) = u (i,j,k)
               vp (i,j,k) = v (i,j,k)
               utf (i,j,k) = u (i,j,k)
               vtf (i,j,k) = v (i,j,k)
               atb (i,j,k,1) = at (i,j,k,1)
               atb (i,j,k,2) = at (i,j,k,2)
            END DO
         END DO
      END DO

         CALL VINTEG (U,UB)
         CALL VINTEG (V,VB)

!$OMP PARALLEL DO PRIVATE (j,i)
      DO j = 1,jmt ! Dec. 5, 2002, Yongqiang Yu
         DO i = 1,imt
            h0p (i,j)= h0 (i,j)
            ubp (i,j)= ub (i,j)
            vbp (i,j)= vb (i,j)
            h0f (i,j)= h0 (i,j)
            h0bf (i,j)= h0 (i,j)
         END DO
      END DO
!
     ISB = 0
     ISC = 0
     IST = 0
!
!     write(*,*)'ok instant'
      return
      end

#if (defined NETCDF) || (defined ALL)
      SUBROUTINE check_err (iret)
#include <netcdf.inc>
      INTEGER :: iret
      IF (iret /= NF_NOERR) THEN
         PRINT *, nf_strerror (iret)
         STOP
      END IF
      END SUBROUTINE check_err
#endif

