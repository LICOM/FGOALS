       subroutine thermal(tt,ss,pp,aa,bb,mask)
!
! calculate the thermal expansion coefficient and
! the saline contraction coefficient using McDougall (1987) equations
!
! input: tt  potential temperature in C
!        ss  salinity in (s-35)/1000
!        p  pressure in Pa
! output:alpha thermal expansion in 1/C
!        beta saline contraction in 1/psu
! 
       use precision_mod
       use param_mod
       use pconst_mod
      IMPLICIT NONE
       real(r8),dimension(imt,jmt,km) :: tt,ss,pp,tmp,aa,bb,mask
!
!  transform the unit
!  ss from (s-35.0)/1000 to s-35.0
!  pp from Pa to db
!
!       ss=ss*1000.0
!       p=p/10000.0
!
      DO K = 1,KM
         DO J = JST,JET
            DO I = 2,IMM
            tmp(i,j,k)=-pp(i,j,k)/OD0/10000.D0*mask(i,j,k)
            END DO
        if (nx_proc==1) then
        tmp(1,j,k)=tmp(imm,j,k)
        tmp(imt,j,k)=tmp(2,j,k)
        end if
         END DO   
      END DO
#ifdef SPMD                                                    
      call exch_boundary(tmp,km)
#endif
!
!      if(mytid.eq.0) write(*,'(6f20.14)')tmp
!	stop
!
      DO K = 1,KM
         DO J = JST,JET
            DO I = 2,IMM
       bb(i,j,k)=(0.785567d-3-0.301985d-5*tt(i,j,k)+0.555579d-7*tt(i,j,k)**2 &
           -0.415613d-9*tt(i,j,k)**3+(ss(i,j,k)*1000.0d0)*(-0.356603d-6+0.788212d-8*tt(i,j,k)&
           +0.408195d-10*tmp(i,j,k)-0.602281d-15*tmp(i,j,k)**2)+(ss(i,j,k)*1000.0d0)**2*(0.515032d-8)&
           +tmp(i,j,k)*(-0.121555d-7+0.192867d-9*tt(i,j,k)-0.2131127d-11*tt(i,j,k)**2)&
           +tmp(i,j,k)**2*(0.176621d-12-0.175379d-14*tt(i,j,k))+tmp(i,j,k)**3*(0.12155d-17))*mask(i,j,k)
!
!
       aa(i,j,k)=(0.665157d-1+0.170907d-1*tt(i,j,k)-0.203814d-3*tt(i,j,k)**2&
            +0.298357d-5*tt(i,j,k)**3-0.255019d-7*tt(i,j,k)**4+(ss(i,j,k)*1000.0d0)*(0.378110d-2&
            -0.846960d-4*tt(i,j,k)-0.164759d-6*tmp(i,j,k)-0.251520d-11*tmp(i,j,k)**2)&
            +(ss(i,j,k)*1000.0d0)**2*(-0.678662d-5)+tmp(i,j,k)*(0.380374d-4-0.933746d-6*tt(i,j,k)&
            +0.791325d-8*tt(i,j,k)**2)+0.512857d-12*tmp(i,j,k)**2*tt(i,j,k)**2&
            -0.302285d-13*tmp(i,j,k)**3)*bb(i,j,k)*mask(i,j,k)
!
!	print*,bb(i,j,k),aa(i,j,k)
            END DO
        if (nx_proc==1) then
        aa(1,j,k)=aa(imm,j,k)
        aa(imt,j,k)=aa(2,j,k)
        bb(1,j,k)=bb(imm,j,k)
        bb(imt,j,k)=bb(2,j,k)
        end if
         END DO   
      END DO
#ifdef SPMD
      call exch_boundary(aa,km)
      call exch_boundary(bb,km)
#endif

!
       return
       end subroutine thermal
