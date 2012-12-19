!-----------------------------------------------------------------------------
!   Processing some variables from flux coupler.
!-----------------------------------------------------------------------------
!
!        By Yongqiang Yu , 16 Apr., 1999
!
!
!
      SUBROUTINE post_cpl

#include <def-undef.h>


use param_mod
use pconst_mod
use tracer_mod
use forc_mod
use buf_mod
use control_mod
use shr_sys_mod
use work_mod,only : wkj
use output_mod,only:spval
!
      implicit none
      real(r8),dimension(:,:),allocatable::tmp_su,tmp_sv
!
         allocate(tmp_su(imt,jmt))
         allocate(tmp_sv(imt,jmt))
!
        tsf=0.0D0
        ssf=0.0D0
        swv=0.0D0
        tmp_su=0.0D0
        tmp_sv=0.0D0

!$OMP PARALLEL DO PRIVATE (i,j)
        ! TODO
        do j=1,jmt
        do i= 2,imt-1 ! 1,imt !LPF 20120822
           TSF(i,j) = (lat1(i,j)+sen(i,j)  +lwup(i,j)+lwdn(i,j )&
                      +netsw(i,j)+melth(i,j)) *OD0CP  ! net heat flux
           SWV(i,j) =   netsw(i,j)                               ! net solar radiation
           NSWV(i,j) = lat1(i,j)+sen(i,j)+lwup(i,j)+lwdn(i,j)&
                      +melth(i,j)                                ! none solar radiation !for BUOY
           SSF(i,j) =  -(prec(i,j)+evap(i,j)+&
                         meltw(i,j)+roff(i,j)) &
                        *34.7*1.0e-3/DZP(1)*OD0                                     ! P+E+melting !linpf 25->DZP(1)
           tmp_su(i,j)= taux(i,j)
           tmp_sv(i,j) = -tauy(i,j)
!           if(vit(i,j,1)<0.5.and.(swv(i,j)>10..or.abs(nswv(i,j)).gt.10..or.abs(tsf(i,j)).gt.10)) &
!                write(*,*)'error onland-heat,vit=,swv=',i,j,vit(i,j,1),swv(i,j),nswv(i,j),tsf(i,j)
!           if(vit(i,j,1)<0.5.and.(abs(tmp_su(i,j))>10..or.abs(tmp_sv(i,j)).gt.10.)) &
!                write(*,*)'error onland-wind,vit=,tauxy=',i,j,vit(i,j,1),tmp_su(i,j),tmp_sv(i,j)
!           if(vit(i,j,1)>0.5.and.(swv(i,j)>=2000..or.abs(nswv(i,j)).ge.2000..or.abs(tsf(i,j)).gt.2000)) &
!                write(*,*)'error onocean-heat,vit=,swv=',i,j,vit(i,j,1),swv(i,j),nswv(i,j),tsf(i,j)
!           if(vit(i,j,1)>0.5.and.(abs(tmp_su(i,j))>10..or.abs(tmp_sv(i,j)).gt.10.)) &
!                write(*,*)'error onocean-wind,vit=,tauxy=',i,j,vit(i,j,1),tmp_su(i,j),tmp_sv(i,j)
!            if(vit(i,j,1)>0.5.and.swv(i,j)==0.0) &
!                 write(*,*)'vit=,swv=',vit(i,j,1),swv(i,j)
!            if(vit(i,j,1)>0.5.and.swv(i,j)==0.0) &
!             write(*,*)'i,j,vit=,swv=',i,j,&
!             90.-wkj(j)/torad,&
!              vit(i,j,1),swv(i,j)
!               if((roff(i,j)*(1-vit(i,j,1))).ne.0) write(*,*)'roff',roff(i,j),i,j
        end do
!LPF 20120819
#ifdef TEST_BUG
             if(mytid==0) then
              if(90-wkj(j)/torad<50.0.and.90.0-wkj(j)/torad>45) then
                 write(*,*) '500B1,j,swv(i,j)',j,90-wkj(j)/torad,(i,i=1,imt)
                 write(*,*) '500B2', (vit(i,j,1),i=1,imt)
                 write(*,*) '500B3', (swv(i,j),i=1,imt)
               endif

             if(90-wkj(j)/torad<=-20.0.and.90.0-wkj(j)/torad>=-22.0) then
               write(*,*) '-200B1,j,swv(i,j)',j,90-wkj(j)/torad,(i,i=1,imt)
               write(*,*) '-200B2',(vit(i,j,1),i=1,imt)
               write(*,*) '-200B3',(swv(i,j),i=1,imt)
             endif
             endif
#endif
!LPF 20120819
        end do
        
        ! lihuimin, 2012.7.17, add 1,1
        call exchange_2d(tsf,1,1)
        call exchange_2d(swv,1,1)
        call exchange_2d(ssf,1,1)
        call exchange_2d(tmp_su,1,1)
        call exchange_2d(tmp_sv,1,1)
        call exchange_2d(duu10n,1,1)
 !LPF 20120904
        call exchange_2d(roff,1,1)
#ifdef TEST_BUG
            do j=1,jmt
             if(mytid==0) then
              if(90-wkj(j)/torad<50.0.and.90.0-wkj(j)/torad>45) then
                 write(*,*) '500F1,j,swv(i,j)',j,90-wkj(j)/torad,(i,i=1,imt)
                 write(*,*) '500F2', (vit(i,j,1),i=1,imt)
                 write(*,*) '500F3', (swv(i,j),i=1,imt)
               endif

             if(90-wkj(j)/torad<=-20.0.and.90.0-wkj(j)/torad>=-22.0) then
               write(*,*) '-200F1,j,swv(i,j)',j,90-wkj(j)/torad,(i,i=1,imt)
               write(*,*) '-200F2',(vit(i,j,1),i=1,imt)
               write(*,*) '-200F3',(swv(i,j),i=1,imt)
              endif
             endif
            enddo
#endif

!linpf 2012 july27
!for ocn output
            lthf = lat1 !latent flux
            sshf = sen !sensible flux
            lwv = lwup + lwdn   !long wave flux
            fresh = ssf
            runoff=roff

!!linpf 2012Jul26 !2012Jul28
        where(vit(:,:,1)<0.5) tsf=spval
        where(vit(:,:,1)<0.5) swv=spval
        where(vit(:,:,1)<0.5) nswv=spval
        where(vit(:,:,1)<0.5) ssf=spval
        where(vit(:,:,1)<0.5) lwv=spval
        where(vit(:,:,1)<0.5) sshf=spval
        where(vit(:,:,1)<0.5) fresh=spval
        where(vit(:,:,1)<0.5) lthf=spval
!
#ifdef USE_OCN_CARBON
        call exchange_2d(pco2)
#endif         
!

! lihuimin, 2012.7.23, ft. lpf
!!$OMP PARALLEL DO PRIVATE (j)
!        do j=1,jmt
!           tsf(imt-1,j)=tsf(1,j)
!           tsf(imt,j)=tsf(2,j)
!           swv(imt-1,j)=swv(1,j)
!           swv(imt,j)=swv(2,j)
!           ssf(imt-1,j)=ssf(1,j)
!           ssf(imt,j)=ssf(2,j)
!           tmp_su(imt-1,j)=tmp_su(1,j)
!           tmp_su(imt,j)=tmp_su(2,j)
!           tmp_sv(imt-1,j)=tmp_sv(1,j)
!           tmp_sv(imt,j)=tmp_sv(2,j)
!           duu10n(imt-1,j)=duu10n(1,j)
!           duu10n(imt,j)=duu10n(2,j)
!        end do

!
! surface stress
!$OMP PARALLEL DO PRIVATE (i,j)
        do j=2,jmt-1
        do i=2,imt
           SU(i,j) = ( tmp_su(i-1,j)+tmp_su(i,j)+tmp_su(i-1,j+1)+tmp_su(i,j+1))*0.25
        end do
        end do
!
!$OMP PARALLEL DO PRIVATE (i,j)
        do j=2,jmt-1
        do i=2,imt
            SV(i,j) = ( tmp_sv(i-1,j)+tmp_sv(i,j)+tmp_sv(i-1,j+1)+tmp_sv(i,j+1))*0.25 !linpf
        end do
        end do
!
! lihuimin, 2012.7.23, ft. lpf
!!$OMP PARALLEL DO PRIVATE (j)
!        do j= 1,jmt
!           SU (1,j) =  SU(imt-1,j)
!           SV (1,j) =  SV(imt-1,j)
!        end do

        ! lihuimin, 2012.7.23, ft. lpf
        call exchange_2d(SU,1,1)
        call exchange_2d(SV,1,1)
!!linpf 2012Jul29
        where(viv(:,:,1)<0.5) su=spval
        where(viv(:,:,1)<0.5) sv=spval
!linpf 2012Jul29
!calculate USTAR

       DO J = 2,jmt-1
         DO I = 2,imt-1
          USTAR(I,J)=sqrt(sqrt(taux(i,j)*taux(i,j)+tauy(i,j)*tauy(i,j))*OD0)*vit(i,j,1) 
         END DO
      END DO        
      
       do i=1,imt
        do j=1,jmt
          if (USTAR(I,J)*vit(i,j,1)>100) then
           write(700+mytid,*)taux(i,j),tauy(i,j),USTAR(I,J),SWV(i,j),NSWV(i,j),i,j,vit(i,j,1),OD0
          endif
       end do
       end do
!
!$OMP PARALLEL DO PRIVATE (i,j)
       do i=1,imt
       do j=1,jmt
          licomqice (i,j)= 0.0D0
       end do
       end do

         deallocate (tmp_su,tmp_sv) !LPF 20120818
        return
        end
