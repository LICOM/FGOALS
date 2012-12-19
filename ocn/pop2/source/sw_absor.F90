!===========================================  
!short wave absorption  based on efolding or
!==========================================      
	  SUBROUTINE sw_absor
#include <def-undef.h>  
use precision_mod
use param_mod   !define km  imt ,jmt 
use pmix_mod    !define pen
use pconst_mod
use sw_mod
use forc_mod


#if (defined SOLARCHLORO)
!=========================
! chl for table look-up
!==========================
         

       chloc=(/ 0.001D0,0.005D0, 0.01D0, 0.02D0,& 
        .03D0,  .05D0,  .10D0,  .15D0, & 
        .20D0,  .25D0,  .30D0,  .35D0, &  
        .40D0,  .45D0,  .50D0,  .60D0, &  
        .70D0,  .80D0,  .90D0, 1.00D0, &
       1.50D0, 2.00D0, 2.50D0, 3.00D0, &
       4.00D0, 5.00D0, 6.00D0, 7.00D0, &
       8.00D0, 9.00D0,10.00D0 /)
       
         A_11=(/0.44210D0, 0.44510D0, 0.44880D0,  0.45630D0, &                                               
       0.46220D0, 0.47150D0, 0.48770D0,  0.49930D0, &                                          
       0.50840D0, 0.51590D0, 0.52230D0,  0.52780D0, &                                               
       0.53260D0, 0.53690D0, 0.54080D0,  0.54740D0, &                                              
       0.55290D0, 0.55760D0, 0.56150D0,  0.56490D0, &                                                 
       0.57570D0, 0.58020D0, 0.58080D0,  0.57880D0, &                                                  
       0.56965D0, 0.55638D0, 0.54091D0,  0.52442D0, &                                                
       0.50766D0, 0.49110D0, 0.47505D0 /)
       
        A_12=( /0.29810D0, 0.29630D0, 0.29400D0, 0.28940D0, &                                                   
       0.28580D0, 0.28000D0, 0.27030D0, 0.26280D0, &                                                   
       0.25710D0, 0.25230D0, 0.24810D0, 0.24440D0, &                                                   
       0.24110D0, 0.23820D0, 0.23560D0, 0.23090D0, &                                                   
       0.22690D0, 0.22350D0, 0.22060D0, 0.21810D0, &                                                  
       0.21060D0, 0.20890D0, 0.21130D0, 0.21670D0, &                                                   
       0.23357D0, 0.25504D0, 0.27829D0, 0.30274D0, &                                                 
       0.32698D0, 0.35056D0, 0.37303D0 /)
       
!       data
        B_11=(/0.02870D0, 0.03010D0, 0.03190D0, 0.03550D0, &                                                  
       0.03840D0, 0.04340D0, 0.05320D0, 0.06120D0, &                                                   
       0.06810D0, 0.07430D0, 0.08000D0, 0.08530D0, &                                                  
       0.09020D0, 0.09490D0, 0.09930D0, 0.10770D0, &                                                   
       0.11540D0, 0.12270D0, 0.12940D0, 0.13590D0, &                                                 
       0.16400D0, 0.18760D0, 0.20820D0, 0.22640D0, &                                                   
       0.25808D0, 0.28498D0, 0.30844D0, 0.32932D0, &                                                
       0.34817D0, 0.36540D0, 0.38132D0 /)
       
!       data  
       B_12=(/ 0.31920D0, 0.32430D0, 0.33060D0, 0.34330D0, &
       0.35370D0, 0.37050D0, 0.40310D0, 0.42620D0, &                                               
       0.44560D0, 0.46210D0, 0.47630D0, 0.48890D0, &                                                    
       0.49990D0, 0.51000D0, 0.51910D0, 0.53470D0, &                                                   
       0.54770D0, 0.55880D0, 0.56820D0, 0.57640D0, &                                                   
       0.60420D0, 0.62060D0, 0.63240D0, 0.64250D0, &                                                  
       0.66172D0, 0.68144D0, 0.70086D0, 0.72144D0, &                                                 
       0.74178D0, 0.76190D0, 0.78155D0 /)
!       write(*,*)'chloc(1)','chloc(31)',chloc(1),chloc(31)
!       write(*,*)'A_11(1)','A_11(31)',A_11(1),A_11(31)
!       write(*,*)'A_12(1)','A_12(31)',A_12(1),A_12(31)
!       write(*,*)'B_11(1)','B_11(31)',B_11(1),B_11(31)
!       write(*,*)'B_12(1)','B_12(31)',B_12(1),B_12(31)	



!tabel const have defined 
!define the depth thick
            arg_max=35.0D0 
!             percm=0.01D0   
              c0=0.0D0       
              c1=1.0D0  ! 
              
         ztr(1)=dzp(1)  !*100
      do k=2,km
           ztr(k)=ztr(k-1)+dzp(k)  !*100 !
      enddo

      chlmin  = chloc(1)
      chlmax  = chloc(M_chl) 
      dlogchl = (log10(chlmax)-log10(chlmin))/real(nsub)

      logchl = log10(chlmin) - dlogchl
       
      do n=0,nsub
        logchl = logchl + dlogchl
        chlamnt = 10.D0**(logchl)  
        do m=1,M_chl-1
          if( chloc(m) .le. chlamnt .and.chlamnt .le. chloc(m+1) ) then
            mc = m
            goto 307
          endif
        end do


307     continue

        w12 = (chlamnt-chloc(mc))/(chloc(mc+1)-chloc(mc))
        w11 = 1.0D0 - w12
        A_1 = A_11(mc)*w11 + A_11(mc+1)*w12
        A_2 = A_12(mc)*w11 + A_12(mc+1)*w12
        B_1 = B_11(mc)*w11 + B_11(mc+1)*w12
        B_2 = B_12(mc)*w11 + B_12(mc+1)*w12
        

         Tr(0,n) = 1.0D0
          do k=1,km-1
#ifdef D_PRECISION
          arg = dmin1(B_1*ztr(k),arg_max) !*percm
#else
          arg = min(B_1*ztr(k),arg_max) !*percm
#endif
          Tr(k,n) = A_1*exp(-arg)
#ifdef D_PRECISION
          arg = dmin1(B_2*ztr(k),arg_max)  !*percm
#else
          arg = min(B_2*ztr(k),arg_max)  !*percm
#endif
          Tr(k,n) = Tr(k,n) + A_2*exp(-arg)
        end do
       end do     
      
!        for every month limited data             
	  
         DO j=JSM,JEM  !1,jmt
          DO i=1,IMM  !1,imt
             DO k=1,km-1
        chl_temp = chloro(i,j)  
!        chl_temp = max(chl_temp,chlmin)               ! set lowest allowed limit
        if(chl_temp.lt.chlmin) chl_temp=chlmin
!        chl_temp = min(chl_temp,chlmax)               ! set highest allowed limit(i,j)
         if(chl_temp.gt.chlmax) chl_temp=chlmax
!         if(chl_temp.lt.0.2D0)D0 chl_temp=0.2D0 D0 for  experiment 3 different
        chlindx = log10(chl_temp/chlmin)/dlogchl ! compute chlorphyll index(i,j)
        chlindx = max(chlindx,0)            ! minimum limit for chl index (i,j)
        chlindx = min(chlindx,nsub)         ! maximum limit for chl index (i,j)
	pen_chl(i,j,k)=Tr(k,chlindx)* OD0CP        !(i,j)
!         if(k.eq.1.and.chloro(i,j).gt.0.1.and.i.eq.10) &
!         write(*,*)'pen_chlsw',pen_chl(i,j,k)/OD0CP,&
!         pen_chl(i,j,k),chloro(i,j)
           END DO             
         END DO
        END DO 
        
#endif
                        
!end by lpf
          RETURN
        END SUBROUTINE sw_absor
