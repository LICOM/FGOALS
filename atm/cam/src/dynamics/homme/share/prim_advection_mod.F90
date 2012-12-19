#ifdef CAM
#else
! sometimes used for debugging REMAP
#undef ZEROVERT
#undef ZEROHORZ
#endif

#if 0
SUBROUTINES:
   prim_advec_tracers()
      Full Euler + hypervis
      oscillatory, non-conservative QNEG column fixer
   prim_advec_tracers_remap_rk2()
      SEM 2D RK2 + monotone remap + hyper viscosity
      SEM 2D RK2 can use sign-preserving or monotone reconstruction

Notes on Lagrange+REMAP advection

From dynamics, we have the velocity on the reference levels which
have density:  dp(t-1)    dp(t)      dp(t+1)

Note that in a hydrostatic model, the density still satisifies the
continuity equation exactly:  
(1)      dp(t+1)-dp(t-1)  + 2dt div(U dp(t)) + 2dt d( eta_dot_dpdn(t) ) = 0 

We introduce the vertically lagrangian levels, which have density
              dp_star(t-1)   dp(t)   dp_star(t+1)
We want dp_star(t) to be the density field which satisfies the 2D 
continuity equation: 
(2)       dp_star(t+1) - dp_star(t-1) + 2dt div(U dp(t) ) = 0

Combining (1) and (2), we have:
  dp_star(t+1) - dp_star(t-1) =  dp(t+1)-dp(t-1)  + 2dt d( eta_dot_dpdn(t) ) 

Thus it is natural to define:
(3)       dp_star(t+1) = dp(t+1) + dt d( eta_dot_dpdn(t) ) 
(4)       dp_star(t-1) = dp(t-1) - dt d( eta_dot_dpdn(t) ) 

If we use a forward-in-time advection scheme:
    dp(t+1)-dp(t)   = dt div(Udp1) + dt d(eta_dot_dpdn1) = 0 
    dp(t+2)-dp(t+1) = dt div(Udp2) + dt d(eta_dot_dpdn2) = 0 
    dp(t+3)-dp(t+2) = dt div(Udp3) + dt d(eta_dot_dpdn3) = 0 
    ---------------
    dp(t+3)-dp(t)   = 3dt div(Udp_sum/3) + 3dt d(eta_dot_dpdn_sum/3) = 0 

so we can define
    dp_star(t+1) = dp(t+1) + 3dt d( eta_dot_dpdn_ave(t) ) 

For RK2 advection of Q:
For consistency, if Q=1
  dp1  = dp(t)- dt div[ U1 dp(t)]        U1 = U(t)
  dp2  = dp1  - dt div[ U2 dp1  ]        U2 = U(t+1) on dp_star(t+1) surface
  dp*  = (dp(t) + dp2 )/2
       =  dp(t) - dt/2   div[ U1 dp(t) + U2 dp1 ]  
and:
  Qdp1 = Qdp(t) - dt div[ U1 Qdp(t)]   limit with dp1 computed above 
  Qdp2 = Qdp1  - dt div[ U2 Qdp1 ]    limit with dp2 computed above
  Qdp* = (Qdp(t) + Qdp2 )/2
OR, if limiter is not active:
(5)  Qdp* = Qdp(t) - .5 dt div[ U1 Qdp(t) +  U2 Qdp1 ]

last step:
  remap Qdp* to Qdp(t+1)   [ dp_star(t+1) -> dp(t+1) ]

Consistency will only work if we can define he mean flux above, and
then pick U1 and U2 so that Udp_ave = [ U1 dp(t) + U2 dp1 ]/2
best way to do this:
  U1 = Udp_ave/dp(t)
  U2 = Udp_ave/dp1
so that
  Qdp(t+1) = Qdp(t) + dt/2 div[  Udp_ave/dp(t) Qdp(t)  + Udp_ave/dp1 Qdp1 ]

Until we switch to forward-in-time scheme for dynamics, what should
we use for U2?  

(A) remapping U(t+1) from dp(t+1) to dp_star(t+1) is expensive.  

(B) Could use (similar Eq 3 above) U2 = U(t+1) + dt*v_vadv(t)
v_vadv can be computed from preq_vadv() subroutine with the
DSSd eta_dot_dpdn.  No need to DSS v_vadv.  

#endif





MODULE pcmpsm_no_pbc
  !
  !**************************************************************************************
  !
  !  Purpose:
  !  Construct sub-grid-scale polynomials using piecewise cubic method (PCM) or 
  !  piecewise spline method with optional monotone filters.
  !
  !  References: PCM - Zerroukat et al., Q.J.R. Meteorol. Soc., 2005. (ZWS2005QJR)
  !              PSM - Zerroukat et al., Int. J. Numer. Meth. Fluids, 2005. (ZWS2005IJMF)
  !
  !
  !      Date       Programmer       Affiliation          Description of change
  !      ====       ==========       ===========          =====================
  !    3/13/06      P.H.Lauritzen    CMS,NCAR             Original code
  !
  !**************************************************************************************
  !
  USE kinds, only : real_kind, int_kind
!  USE constants        ! common constants
  IMPLICIT NONE

  real(kind=real_kind) :: zero=0
  real(kind=real_kind) :: one=1
  real(kind=real_kind) :: two=2
  real(kind=real_kind) :: three=3
  real(kind=real_kind) :: four=4
  real(kind=real_kind) :: six=6
  real(kind=real_kind) :: nine=9
  real(kind=real_kind) :: twelfe=12
  real(kind=real_kind) :: half = .5
  real(kind=real_kind) :: third = 1/3d0
  REAL(KIND=real_kind) :: tiny = 1e-12

CONTAINS
  SUBROUTINE CUBIC_PARAMETERS(y_left_cv,mass,dy,a0,a1,a2,a3,piecewise,ns,nf)
    !
    ! CUBIC_PARAMETERS computes left and right boundary values of cv 
    ! and the slope at the centre of cv
    ! y_left_cv   = left position of CV 
    ! mass        = mass of CV
    ! dy          = size (1-DIMENSIONal lenght) of CV
    ! a0,a1,a2,a3 = coefficents of cubic 
    ! piecewise   =piecewise option (0(constant, 1 linear, 2 Parabolic and 3 Cubic)
    !
    IMPLICIT NONE
    INTEGER(KIND=Int_kind) ::  ns, nf, piecewise, j
    REAL(KIND=real_kind), DIMENSION(ns:nf+1) :: y_left_cv
    REAL(KIND=real_kind), DIMENSION(ns:nf):: mass,dy,a0,a1,a2,a3
    
    !  locals
    REAL(KIND=real_kind), DIMENSION(ns:nf):: rho_left_cv,rho_right_cv,slope_rho
    REAL(KIND=real_kind), DIMENSION(ns:nf):: y_centre_cv, rho_bar
    INTEGER(KIND=int_kind), PARAMETER :: standard_number_of_cvs = 4
    REAL(KIND=real_kind), DIMENSION(ns:nf+1) :: rho_left_cv_all,slope_im,slope_ip
    INTEGER(KIND=int_kind) :: j0, j1, cv_start
    REAL(KIND=real_kind) :: y0, y1, y2, y3, y4, y, m1, m2, m3, m4
    
    a0=zero; a1=zero;a2=zero;a3=zero
    
    j0 = standard_number_of_cvs/2
    j1 = standard_number_of_cvs - 1
    rho_left_cv_all=zero ; slope_im=zero; slope_ip=zero
    y_centre_cv = half * ( y_left_cv(ns:nf) + y_left_cv(ns+1:nf+1) ) 
    
    IF ((nf - ns + 1) < 4 ) THEN  ! use piecewise constant for less than 4 CVs
       a0 = mass/dy
       RETURN
    ENDIF
    
    DO j = ns , nf + 1
       !
       ! Use extrapolation if neccesary
       !
       cv_start = MIN( MAX( j - j0, ns ) + j1, nf ) - j1
       !
       y0 = y_left_cv(cv_start)
       y1 = y_left_cv(cv_start+1) - y0  !change of variable (y) to have y0=0
       y2 = y_left_cv(cv_start+2) - y0
       y3 = y_left_cv(cv_start+3) - y0
       y4 = y_left_cv(cv_start+4) - y0
       m1 =   mass(cv_start)
       m2 =   mass(cv_start+1) + m1
       m3 =   mass(cv_start+2) + m2
       m4 =   mass(cv_start+3) + m3 
       m1 = m1 / (y1*(y1-y2)*(y1-y3)*(y1-y4))
       m2 = m2 / (y2*(y2-y1)*(y2-y3)*(y2-y4))
       m3 = m3 / (y3*(y3-y1)*(y3-y2)*(y3-y4))      
       m4 = m4 / (y4*(y4-y1)*(y4-y2)*(y4-y3))
       
       y  = y_left_cv(j) - y0
       rho_left_cv_all(j) =                                                   &
            (y*(y*(four*y-three*(y2+y3+y4))+two*(y3*y2+y2*y4+y3*y4))-y3*y2*y4)*m1+   &
            (y*(y*(four*y-three*(y1+y3+y4))+two*(y1*y3+y1*y4+y3*y4))-y1*y3*y4)*m2+   &   
            (y*(y*(four*y-three*(y2+y4+y1))+two*(y2*y4+y1*y2+y1*y4))-y1*y2*y4)*m3+   &
            (y*(y*(four*y-three*(y2+y3+y1))+two*(y3*y2+y1*y2+y1*y3))-y1*y2*y3)*m4   
       y  = y_centre_cv( MAX( ns, j - 1 ) ) - y0
       slope_im(j) =                                                &       
            (y*(twelfe*y-six*(y2+y3+y4))+two*(y3*y2+y2*y4+y3*y4))*m1+    &
            (y*(twelfe*y-six*(y1+y3+y4))+two*(y1*y3+y1*y4+y3*y4))*m2+    &
            (y*(twelfe*y-six*(y2+y4+y1))+two*(y2*y4+y1*y2+y1*y4))*m3+    &
            (y*(twelfe*y-six*(y2+y3+y1))+two*(y3*y2+y1*y2+y1*y3))*m4
       y  = y_centre_cv( MIN( j, nf    ) ) - y0
       slope_ip(j) =                                                &
            (y*(twelfe*y-six*(y2+y3+y4))+two*(y3*y2+y2*y4+y3*y4))*m1+    &
            (y*(twelfe*y-six*(y1+y3+y4))+two*(y1*y3+y1*y4+y3*y4))*m2+    &
            (y*(twelfe*y-six*(y2+y4+y1))+two*(y2*y4+y1*y2+y1*y4))*m3+    &
            (y*(twelfe*y-six*(y2+y3+y1))+two*(y3*y2+y1*y2+y1*y3))*m4   
    ENDDO
    
    rho_bar = mass / dy 
    
    rho_left_cv_all(ns  ) = rho_bar(ns)
    rho_left_cv_all(ns+1) = half*(rho_bar(ns) + rho_bar(ns+1))
    rho_left_cv_all(nf+1) = rho_bar(nf)
    rho_left_cv_all(nf  ) = half*(rho_bar(nf) + rho_bar(nf-1))
    
    slope_rho(ns:nf) = half * ( slope_im(ns+1:nf+1) + slope_ip(ns:nf) )
    rho_left_cv(ns:nf) = rho_left_cv_all(ns:nf)
    rho_right_cv(ns:nf)= rho_left_cv_all(ns+1:nf+1)
    
    slope_rho = slope_rho * dy
    
    
    IF (piecewise == 0 ) THEN
       a0 = rho_bar
    ELSEIF (piecewise == 1 ) THEN
       WHERE ( (rho_bar-rho_left_cv) > (rho_bar-rho_right_cv) )
          a0=rho_left_cv; a1=two*(rho_bar-rho_left_cv)
       ELSEWHERE
          a0=two*rho_bar-rho_right_cv; a1=two*(rho_right_cv-rho_bar) 
       END WHERE
    ELSEIF (piecewise == 2 ) THEN
       a0 = rho_left_cv 
       a1 = -four*rho_left_cv - two*rho_right_cv + six*rho_bar   
       a2 =  three*rho_left_cv + three*rho_right_cv - six*rho_bar
    ELSEIF (piecewise == 3 ) THEN
       a0 = rho_left_cv 
       a1 = -six*rho_left_cv + six*rho_bar   - two*slope_rho
       a2 =  nine*rho_left_cv - three*rho_right_cv - six*rho_bar + six*slope_rho
       a3 = -four*rho_left_cv + four*rho_right_cv - four*slope_rho
    ENDIF
  END SUBROUTINE CUBIC_PARAMETERS


  !----------------------------------------------------------------------------------------
  ! This routine to modIFy the PARAMETERs of a piecewise cubic to satisfy monotonocity
  !-----------------------------------------------------------------------------------------
  SUBROUTINE  CUBIC_MONOTONE(a0,a1,a2,a3,rho_bar,                    &
       ns, nf, rho_min, rho_max, check_global, piecewise)
    
    ! Given a set of PARAMETERs {rho_left, rho_right, rho_bar, slope_rho} for
    ! each peacewise cubic for each CV(i),i=ns:nf
    ! this routine modIFy these PARAMETERs so that the piecewise cubic is
    ! monotone by:
    ! 1- detect which interval i, WHERE these PARAMETERs will be modIFied.
    !    This is DOne using one of the 3 dIFferent detection filters
    !    MONO_FILTER(1,2,3)
    ! 2- Once a CV is detected to be modIFied, the algorithm decide
    !    whether to reduce PCM to PPM or PLM or PCoM accordingly
    ! 3- Compute the new piecewise cubic coefficents  {a0,a1,a2,a3} after
    !    inforcing the monotonic constraint
    
    IMPLICIT NONE
    INTEGER(KIND=Int_kind) ::  ns, nf, check_global, piecewise
    REAL(KIND=real_kind), DIMENSION(ns:nf) :: rho_left, rho_right, rho_bar, slope_rho 
    REAL(KIND=real_kind) :: rho_min, rho_max  
    !  locals ------------------------------------------------            
    INTEGER(KIND=Int_kind), DIMENSION(ns:nf) ::  filter_code, filter_code1,peaks                                         
    REAL(KIND=real_kind), DIMENSION(ns:nf) ::  rho_left_r, rho_right_r, slope_r  &
         ,  a0_r, a1_r, a2_r, a3_r            &
         ,  a0, a1, a2, a3       
    REAL(KIND=real_kind), DIMENSION(ns:nf+1)  ::  rho_left_all
    REAL(KIND=real_kind), DIMENSION(2,ns:nf) :: peaks_val
    
    INTEGER(KIND=Int_kind) :: i, mono_degree, im, ip
    REAL(KIND=real_kind) :: temp1, temp2
    
    ! section 1: detect which CV WHERE the cubic-PARAMETERs are to be modIFied
    !            IF    filter_code(i)=0 : cubic-PARAMETERs are not modIFied
    !            ELSE  filter_code(i)=1 : cubic-PARAMETERs are to be modIFied 
    
    
    
    filter_code1 = 0; filter_code=0
    rho_left  = a0; rho_right(ns:nf-1)=rho_left(ns+1:nf)
    rho_right(nf)=a0(nf)+a1(nf)+a2(nf)+a3(nf)
    slope_rho = a1 + a2 + (three/four)* a3 
    
    rho_left_all(ns:nf)= rho_left
    rho_left_all(nf+1)=rho_right(nf)
    
    CALL MONO_FILTER4(rho_left_all, rho_bar, ns, nf,      &
         rho_min, rho_max, check_global, filter_code1 )
    
    
    rho_left  = rho_left_all(ns:nf)
    rho_right = rho_left_all(ns+1:nf+1)
    
    CALL CUBIC_COEFFS(rho_left, rho_right, rho_bar, slope_rho, &
         ns, nf, a0, a1, a2, a3, piecewise)
    CALL PROPERTIES_CUBIC(a0, a1, a2, a3, piecewise, ns, nf,  &
         peaks, peaks_val)
    
    CALL MONO_FILTER2(rho_left_all, ns, nf, peaks, peaks_val,   & 
         rho_min, rho_max, check_global, filter_code)
    
    
    CALL CUBIC_REDUCED(rho_left, rho_right, rho_bar, slope_rho,  &            
         2, ns, nf, rho_left_r, rho_right_r, slope_r)
    
    CALL CUBIC_COEFFS(rho_left_r, rho_right_r, rho_bar, slope_r,  &
         ns, nf, a0_r, a1_r, a2_r, a3_r, piecewise)
    
    DO i=ns,nf                                                        
       IF ( filter_code(i) > 0 ) THEN                                                  
          IF ( abs(a2_r(i)) > zero ) THEN
             temp1=-a1_r(i)/(two*a2_r(i))  ! temp1= position of MAXima
          ELSE
             temp1 = two
          ENDIF
          IF ( temp1 <= zero .OR. temp1 >= one ) THEN  ! MAXima is outside [0,1]
             a0(i)=rho_left(i); a1(i)=-four*rho_left(i)-two*rho_right(i)+six*rho_bar(i)
             a2(i)=three*rho_left(i)+three*rho_right(i)-six*rho_bar(i); a3(i)=zero
          ELSE      ! rho_bar is inside [rho_left,rho_right]             
             CALL QUADRATIC_REDUCE2(a0(i),a1(i),a2(i),   &
                  rho_left(i),rho_right(i),rho_bar(i))
             
             a3(i)=zero                      
          ENDIF
       ENDIF    ! END IF ( filter_code(i) > 0 );  
    ENDDO      ! END  for i=1:n;  
  END SUBROUTINE CUBIC_MONOTONE
  !------------------------------------------------------------------------------------  
  SUBROUTINE MONO_FILTER2(y,  ns, nf, peaks, peaks_val,        & 
       ymin, ymax, check_global, yf )
    !  this type of filter check for details of local cubic and neighbouring
    !  discrete left-values
    IMPLICIT NONE
    INTEGER(KIND=Int_kind) :: i, ns, nf, check_global
    REAL(KIND=real_kind),    DIMENSION(ns:nf+1):: y
    REAL(KIND=real_kind),    DIMENSION(ns:nf)  :: dy
    INTEGER(KIND=Int_kind), DIMENSION(ns:nf)  :: peaks, yf
    REAL(KIND=real_kind), DIMENSION(2,ns:nf)   :: peaks_val
    REAL(KIND=real_kind) :: ymin, ymax
    ! locals
    INTEGER(KIND=Int_kind) :: im1, im2, ip1, ip2 
    
    dy(ns:nf) = y(ns+1:nf+1)-y(ns:nf)                          
    WHERE ( ABS(dy) < tiny ) dy=zero   
    
    DO i=ns,nf
       im1=MAX(ns,i-1) ; im2=MAX(ns,i-2)  
       ip1=MIN(nf,i+1); ip2=MIN(nf,i+2);                              
       IF ( ( peaks(i) /= 0 )              .AND.         &
            ( (dy(im2)*dy(im1) <= tiny) .OR.          &
            (dy(ip1)*dy(ip2) <= tiny) .OR.          &
            (dy(im1)*dy(ip1) >= tiny) .OR.          &
            (dy(im1)*float(peaks(i)) <= tiny)     ) ) THEN                                                
          yf(i)=1                                                       
       ENDIF
       IF ( peaks(i) > 1 ) THEN
          yf(i)=1
       ENDIF
       IF ( (check_global == 1) .AND.      &
            ( (y(i) >= ymax)     .OR.      &
            (y(i) <= ymin)      )        ) THEN                                          
          yf(i) = 2                                                                          
       ENDIF
       IF ( (check_global == 1)                .AND.   & 
            ( (peaks_val(2,i) > ymax-tiny)  .OR.    &
            (peaks_val(1,i) < ymin+tiny)   )      ) THEN                                          
          yf(i) = 2                                                                         
       ENDIF
    ENDDO
  END SUBROUTINE MONO_FILTER2
      
 !-----------------------------------------------------------------------------------------     
  SUBROUTINE MONO_FILTER4(yl, ybar, ns, nf, ymin, ymax, check_global, yf)          
    IMPLICIT NONE
    INTEGER(KIND=Int_kind) ::  i, ns, nf, check_global
    REAL(KIND=real_kind),    DIMENSION(ns:nf+1):: yl        ! left position
    REAL(KIND=real_kind),    DIMENSION(ns:nf)  :: ybar, dy  ! average values & size of CVs
    INTEGER(KIND=Int_kind), DIMENSION(ns:nf)  :: yf        ! yf=0 (piecewise is already monotone
    ! yf=1 (piecewise is not monotone
    REAL(KIND=real_kind) :: ymin, ymax                      ! MIN & MAX values for ybar
    
    
    ! local 
    INTEGER(KIND=Int_kind) :: im1, im2, im3, ip1
    REAL(KIND=real_kind) ::  r1, r2
    
    
    dy(ns:nf-1) = ybar(ns+1:nf)-ybar(ns:nf-1); dy(nf)=dy(nf-1)                          
    WHERE ( ABS(dy) < tiny ) dy=zero
    
    DO i=ns,nf
       im1=MAX(ns,i-1); im2=MAX(ns,i-2); im3=MAX(ns,i-3)
       ip1=MIN(nf,i+1)
       
       IF ( (ybar(i)-yl(i))*(yl(i)-ybar(im1)) >= zero  )  THEN
          yf(i)=0            
       ELSEIF (      dy(im2)*(yl(i)-ybar(im1)) > zero     &
            .AND. dy(im2)*dy(im3)           > zero     &
            .AND. dy(i)*dy(ip1)             > zero     &
            .AND. dy(im2)*dy(i)             < zero     ) THEN                                              
          yf(i)=0
       ELSE
          yf(i)=1; yf(im1)=1
          r1=ABS(yl(i)-ybar(im1)); r2=ABS(yl(i)-ybar(i))                  
          IF (r1 >= r2 ) THEN
             yl(i)=ybar(i)
          ELSE   
             yl(i)=ybar(im1)
          ENDIF
       ENDIF
    ENDDO
    IF( check_global == 1) THEN
       WHERE (yl > ymax) 
          yl = ymax
       ELSEWHERE (yl < ymin) 
          yl = ymin              
       END WHERE
    ENDIF
  END SUBROUTINE MONO_FILTER4
      
  !--------------------------------------------------------------------------------      
  
  SUBROUTINE CUBIC_COEFFS(rho_left, rho_right, rho_bar, slope_rho, ns, nf, &
       a0, a1, a2, a3, piecewise)      
    ! this routine simply computes the 4 cubic coeffients {a0, a1, a2, a3}
    ! knowing the cubic PARAMETERs {rho_left,rho_right,rho_bar,slope_rho}
    ! Zerroukat et.al., Q.J.R. Meteorol. Soc., Vol. 128, pp. 2801-2820 (2002).
    
    IMPLICIT NONE
    INTEGER(KIND=Int_kind) ::  ns, nf,piecewise
    REAL(KIND=real_kind), DIMENSION(ns:nf)  ::  rho_left, rho_right, slope_rho, rho_bar  &
         ,  a0, a1, a2, a3
         
    IF (piecewise == 0 ) THEN
       a0 = rho_bar; a1=zero; a2=zero; a3=zero
    ELSEIF (piecewise == 1 ) THEN
       WHERE ( (rho_bar-rho_left) > (rho_bar-rho_right) )
          a0=rho_left; a1=two*(rho_bar-rho_left)
          a2=zero; a3=zero
       ELSEWHERE
          a0=two*rho_bar-rho_right; a1=two*(rho_right-rho_bar)
          a2=zero; a3=zero
       END WHERE
    ELSEIF (piecewise == 2 ) THEN
       a0 = rho_left 
       a1 = -four*rho_left - two*rho_right + six*rho_bar   
       a2 =  three*rho_left + three*rho_right - six*rho_bar
       a3 = zero
    ELSEIF (piecewise == 3 ) THEN
       a0 = rho_left 
       a1 = -six*rho_left + six*rho_bar   - two*slope_rho
       a2 =  nine*rho_left - three*rho_right - six*rho_bar + six*slope_rho
       a3 = -four*rho_left + four*rho_right - four*slope_rho
    ENDIF
  END SUBROUTINE CUBIC_COEFFS
  !---------------------------------------------------------------------------------
  SUBROUTINE CUBIC_REDUCED(rho_left, rho_right, rho_bar, slope_rho,          &     
       piecewise, ns, nf, rho_left_r, rho_right_r, slope_r) 
    
    ! Given (rho_left, rho_right, slope_rho) of a full piecewise cubic
    ! this routine RETURNs the modIFied/reduced (rho_left_r, rho_right_r,slope_rho_r) 
    ! IF the cubic is reduced to 
    ! parabola (when  piecewise == 2)
    ! Linear   (      piecewise == 1)
    ! constand (      piecewise == 0)
    INTEGER(KIND=Int_kind) ::  ns, nf, piecewise
    REAL(KIND=real_kind), DIMENSION(ns:nf)  ::  rho_left, rho_right, rho_bar, slope_rho   &
         ,  rho_left_r, rho_right_r, slope_r
    
    IF (piecewise == 2) THEN 
       rho_left_r=rho_left
       rho_right_r=rho_right
       slope_r=rho_right-rho_left 
    ELSEIF (piecewise == 1) THEN
       slope_r(ns+1:nf-1)=half*(rho_bar(ns:nf-2)-rho_bar(ns+2:nf))
       slope_r(ns)=rho_bar(ns+1)-rho_bar(ns)
       slope_r(nf)=rho_bar(nf)-rho_bar(nf-1)
       rho_left_r  =rho_bar-half*slope_r
       rho_right_r =rho_bar+half*slope_r 
    ELSEIF (piecewise == 0) THEN
       rho_left_r = rho_bar
       rho_right_r =rho_bar
       slope_r=zero     
    ENDIF  
  END SUBROUTINE CUBIC_REDUCED
!----------------------------------------------------------------------------
  SUBROUTINE PROPERTIES_CUBIC(a0, a1, a2, a3, piecewise,    &
       ns, nf, peaks, peaks_val  )
    
    ! subroutine to RETURN properties of a cubic f(x):
    ! f(x) = a0 + a1 * x + a2 * x^2 + a3 * x^3. 
    ! df/dx = a1 + 2 a2 x + 3 a3 x^2 = a x^2 + b x + c 
    ! d2f/dx2 = 2 a x + b   
    
    IMPLICIT NONE
    INTEGER(KIND=Int_kind) ::  i, ns, nf, piecewise
    REAL(KIND=real_kind), DIMENSION(ns:nf)  ::  rho_left, rho_right, slope_rho  &
         ,  a0, a1, a2, a3
    REAL(KIND=real_kind), DIMENSION(2, ns:nf)  :: peaks_val
    INTEGER(KIND=Int_kind), DIMENSION(ns:nf)  :: peaks
    ! locals 
    REAL(KIND=real_kind)    :: xm1, xm2, f_xm1, f_xm2, a, b, c, delta, ddf_xm1,ddf_xm2 
    
    peaks=0 
    peaks_val(1,ns:nf)=MIN(a0,a0+a1+a2+a3)
    peaks_val(2,ns:nf)=MAX(a0,a0+a1+a2+a3)  
    
    SELECT CASE(piecewise)
    CASE (:1)  ! piecewise constant or linear
       peaks=0;  
    CASE (2)   ! piecewise parabolic
       
       DO i=ns,nf
          
          xm1=-a1(i)/(two*a2(i))
          f_xm1=a0(i) + a1(i)*xm1 + a2(i)*xm1**two
          
          IF (xm1 <= zero.OR. xm1 >= one) THEN
             peaks(i)=0
          ELSE 
             IF (a2(i) > zero) THEN         ! d2f/dx2 > 0 = MINima
                peaks(i)=-1;  peaks_val(1,i)= f_xm1
             ELSE        ! d2f/dx2 < 0 = MAXima
                peaks(i)=+1; peaks_val(2,i)= f_xm1
             ENDIF
          ENDIF
          
       ENDDO
    CASE (3) ! piecewise cubic
       DO i=ns,nf        
          a=three*a3(i); b=two*a2(i); c=a1(i); delta= b**two-four*a*c
          IF (delta <= 0) THEN
             peaks(i)=0
          ELSE         
             xm1 = (-b-sqrt(delta)) / (two*a)     ! location of extrema 1
             xm2 = (-b+sqrt(delta)) / (two*a)     ! location of extrema 2
             f_xm1 = a0(i) + a1(i)*xm1 + a2(i)*xm1**two + a3(i)*xm1**three
             f_xm2 = a0(i) + a1(i)*xm2 + a2(i)*xm2**two + a3(i)*xm2**three
             ddf_xm1 = two*a*xm1 + b              ! value of second derivative @ extrema 1
             ddf_xm2 = two*a*xm2 + b              ! value of second derivative @ extrema 2
             IF ((xm1 <= zero.OR. xm1 >= one) .AND.    &
                  (xm2 <= zero.OR. xm2 >= one)          ) THEN ! both outside intervals
                peaks(i)=0 
             ELSEIF (xm1 > zero.AND. xm1 < one .AND.  &
                  xm2 > zero.AND. xm2 < one        )      THEN  ! both inside
                peaks(i)=2 
                peaks_val(1,i)=MIN(f_xm1, f_xm2)
                peaks_val(2,i)=MAX(f_xm1, f_xm2) 
             ELSEIF ((xm1 > zero) .AND. (xm1 < one)) THEN   ! only xm1 inside                                         
                IF (ddf_xm1 > zero) THEN                     ! xm1 is a MINima 
                   peaks(i)=-1; peaks_val(1,i)= f_xm1                                        
                ELSE                                    ! xm1 is a MAXima 
                   peaks(i)=+1; peaks_val(2,i)= f_xm1
                ENDIF
             ELSEIF (xm2 > zero.AND. xm2 < one ) THEN   ! only xm2 inside                                            
                IF (ddf_xm2 > zero) THEN                     ! xm2 is a MINima 
                   peaks(i)=-1 ; peaks_val(1,i)= f_xm2                                        
                ELSE                                     ! xm2 is a MAXima
                   peaks(i)=+1 ; peaks_val(2,i)= f_xm2
                ENDIF
             ENDIF
          ENDIF  ! IF delta > 0;                   
       ENDDO    ! END DO i=one..
    END SELECT
  END SUBROUTINE PROPERTIES_CUBIC


  !===================================================================================================
!  SUBROUTINE QUADRATIC_SPLINE(dx,m,rho_left,rho_right,a0,a1,a2,n,boundary_condition)
  SUBROUTINE QUADRATIC_SPLINE(dx,m,a0,a1,a2,n)
    !
    ! this routine computes the coefficents (a0,a1,a2) & nodes-values (rho_left, rho_right)
    ! of a piecewise quadratic-spline  (similar to PPM)
    ! dx = mesh size of CV (1-DIMENSIONal)
    ! m = mass of CVs
    !     
    IMPLICIT NONE
    INTEGER(KIND=Int_kind) :: n
    REAL(KIND=real_kind), DIMENSION(n) :: dx,m,rho_left,rho_right,a0,a1,a2
    ! locals
    REAL(KIND=real_kind), DIMENSION(n+1) :: rho_left_all
    REAL(KIND=real_kind), DIMENSION(n)   :: h, rho_bar
    REAL(KIND=real_kind), DIMENSION(n+1) :: rhs,lower_diag,diag,upper_diag
    INTEGER(KIND=Int_kind) :: i, boundary_condition
    REAL(KIND=real_kind) :: upper_corner,lower_corner 
    
    ! boundary_condition = 1  ! natural spline (f'=0) for non-cyclic boundaries
    ! boundary_condition = 2  ! cyclic boundary condition
    
    h = one/dx ; rho_bar = m * h          
    rhs = zero; lower_diag = zero; diag = zero; upper_diag = zero
    
    rhs(2:n) = three*(rho_bar(2:n)*h(2:n) + rho_bar(1:n-1)*h(1:n-1)) 
    lower_diag(2:n) = h(1:n-1)
    diag(2:n) = two*(h(2:n) + h(1:n-1))
    upper_diag(2:n) = h(2:n)

    boundary_condition = 1 !ADD PHL
    IF ( boundary_condition == 1 ) THEN      
       rhs(1)=three*rho_bar(1); lower_diag(1)=one; diag(1)=two; upper_diag(1)=one
       rhs(n+1)=three*rho_bar(n); lower_diag(n+1)=one; diag(n+1)=two; upper_diag(n+1)=one
!PHL       CALL TRIDIAG_SYSTEM(lower_diag,diag,upper_diag,rhs,rho_left_all,n+1)
       CALL tridiag(lower_diag,diag,upper_diag,rhs,n+1) !ADD PHL
       rho_left_all = rhs                               !ADD PHL
    ELSEIF ( boundary_condition == 2 ) THEN 
       rhs(1) = three*(rho_bar(1)*h(1) + rho_bar(n)*h(n)) 
       lower_diag(1)= h(n); diag(1) = two*(h(1)+h(n)); upper_diag(1)=h(1)             
       upper_corner = h(n); lower_corner = h(n)              
!PHL       CALL TRIDIAG_CYCLIC_SYSTEM_V2(lower_diag(1:n),diag(1:n),upper_diag(1:n), &
!PHL            rhs(1:n),rho_left_all(1:n),n               )
       CALL tridiag_per(lower_diag(1:n),diag(1:n),upper_diag(1:n),rhs(1:n),n) !ADD PHL
       rho_left_all(1:n) = rhs                                                !ADD PHL
       rho_left_all(n+1) = rho_left_all(1)
    ENDIF
    rho_left(1:n)= rho_left_all(1:n); rho_right(1:n)= rho_left_all(2:n+1)
    
    a0 = rho_left
    a1 = -four* rho_left - two* rho_right + six*rho_bar
    a2 = +three* rho_left + three* rho_right - six*rho_bar         
  END SUBROUTINE QUADRATIC_SPLINE
  
  
  !-----------------------------------------------------------------------------------------
  SUBROUTINE  QUADRATIC_SPLINE_MONOTONE(a0,a1,a2,rho_bar,dx,        &
       ns, nf, rho_min, rho_max, check_global)
    !PHL       ns, nf, rho_min, rho_max, check_global, piecewise)
    !-----------------------------------------------------------------------------------
    ! Given a set of PARAMETERs {a0,a1,a2,a3,a4} for each peacewise quartic 
    ! for each CV(i),i=ns:nf, this routine modIFy these PARAMETERs so that the piecewise
    ! is monotone.            
    !-----------------------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER(KIND=Int_kind) ::  ns, nf, check_global, piecewise, mz_debug,  &
        filter_option
    
    REAL(KIND=real_kind), DIMENSION(ns:nf) :: a0,a1,a2,rho_bar,dx
    REAL(KIND=real_kind) :: rho_min, rho_max  
    !  locals ------------------------------------------------ 
    REAL(KIND=real_kind), DIMENSION(ns:nf) :: rho_left, rho_right           
    INTEGER(KIND=Int_kind), DIMENSION(ns:nf) ::  filter_code,peaks 
    REAL(KIND=real_kind), DIMENSION(ns:nf+1)  ::  rho_left_all
    REAL(KIND=real_kind), DIMENSION(2,ns:nf)  ::  peaks_val
    
    INTEGER(KIND=Int_kind) :: i, mono_degree
    REAL(KIND=real_kind)    :: temp,alfa,pi
    
    !----------------------------------------------------------------------------------
    ! section 1: detect which CV WHERE the quadratic spline is to be modIFied
    !            IF    filter_code(i)=0 : piecewise-PARAMETERs are not modIFied
    !            ELSE  filter_code(i)=1 : piecewise-PARAMETERs are to be modIFied 
    !----------------------------------------------------------------------------------
    
    pi=acos(-one)                       
    filter_code = 0
    rho_left_all(ns:nf)= a0
    rho_left_all(nf+1) = a0(nf)+a1(nf)+a2(nf)
    
    CALL MONO_FILTER4(rho_left_all, rho_bar, ns, nf,  &
         rho_min, rho_max, check_global, filter_code  )
    
    rho_left  = rho_left_all(ns:nf)
    rho_right = rho_left_all(ns+1:nf+1)
    
    CALL QUADRATIC_SPLINE_COEFFS(rho_left,rho_right,rho_bar,a0,a1,a2,ns,nf)
    
    CALL PROPERTIES_QUADRATIC_SPLINE(a0,a1,a2,peaks,peaks_val,ns,nf)
    CALL MONO_FILTER2(rho_left_all, ns, nf, peaks, peaks_val,        &
         rho_min, rho_max, check_global, filter_code )
    
    
    !----------------------------------------------------------------------------------             
    ! section 2:   For those CVs WHERE the PARAMETERs are to be modIFied decide wether 
    !              it will be reduced to piecewise linear or constant accordingly
    !----------------------------------------------------------------------------------
    
    
    DO i=ns,nf                                                          
       IF (filter_code(i) > 0) THEN
          CALL QUADRATIC_REDUCE2(a0(i),a1(i),a2(i),  &
               rho_left(i),rho_right(i),rho_bar(i)  )
       ENDIF
    ENDDO
  END SUBROUTINE QUADRATIC_SPLINE_MONOTONE

!--------------------------------------------------------------------------------------
  SUBROUTINE QUADRATIC_SPLINE_COEFFS(rho_left,rho_right,rho_bar,a0,a1,a2,ns,nf) 
    !-----------------------------------------------------------------------------     
    ! this routine simply computes the 3 quadratic spline coeffients {a0, a1, a2}
    ! knowing the quadratic spline parameters {rho_left,rho_right,rho_bar}
    ! Zerroukat et.al., Q.J.R. Meteorol. Soc., Vol. 128, pp. 2801-2820 (2002).
    !----------------------------------------------------------------------------- 
    
    IMPLICIT NONE
    INTEGER(KIND=Int_kind) ::  ns, nf
    REAL(KIND=real_kind), DIMENSION(ns:nf) :: rho_left,rho_right,rho_bar,a0,a1,a2 
    
    a0 = rho_left 
    a1 = -four*rho_left - two*rho_right + six*rho_bar  
    a2 =  three*rho_left + three*rho_right - six*rho_bar 
  END SUBROUTINE QUADRATIC_SPLINE_COEFFS

!---------------------------------------------------------------------------------
  SUBROUTINE QUADRATIC_REDUCE2(a0,a1,a2,rho_left,rho_right,rho_bar)            
    REAL(KIND=real_kind) ::  a0,a1,a2,rho_left,rho_right,rho_bar
    !locals
    REAL(KIND=real_kind) :: level1, level2, level3, level4, level5
    
    level1 = rho_left
    level2 = two*third*rho_left+third*rho_right
    level3 = half*rho_left+half*rho_right 
    level4 = third*rho_left+two*third*rho_right
    level5 = rho_right
    
    IF (rho_right >= rho_left ) THEN 
       IF ( rho_bar <= level1 .OR. rho_bar >= level5) THEN       !out range => constant
          a0 = rho_bar ; a1 = zero; a2=zero                                           
       ELSEIF( rho_bar > level1 .AND. rho_bar < level2 ) THEN    !parabola zero slope left         
          a0=rho_left; a1=zero; a2=three*(rho_bar-rho_left)
       ELSEIF (  rho_bar >  level4 .AND. rho_bar < level5 ) THEN !parabola zero slope right
          a0=-two*rho_right+three*rho_bar 
          a1=+six*rho_right-six*rho_bar
          a2=-three*rho_right+three*rho_bar                                           
       ENDIF
    ELSE                                           
       IF ( rho_bar >= level1 .OR. rho_bar <= level5) THEN       !out range => constant
          a0 = rho_bar ; a1 = zero; a2=zero                                           
       ELSEIF( rho_bar < level1 .AND. rho_bar > level2 ) THEN    !parabola zero slope left         
          a0=rho_left; a1=zero; a2=three*(rho_bar-rho_left)
       ELSEIF (  rho_bar <  level4 .AND. rho_bar > level5 ) THEN !parabola zero slope right
          a0=-two*rho_right+three*rho_bar 
          a1=+six*rho_right-six*rho_bar
          a2=-three*rho_right+three*rho_bar                                           
       ENDIF
    ENDIF
  END SUBROUTINE QUADRATIC_REDUCE2
  
  !----------------------------------------------------------------------------
  SUBROUTINE PROPERTIES_QUADRATIC_SPLINE(a0,a1,a2,peaks,peaks_val,ns,nf)
    ! subroutine to RETURN properties (peaks,peaks_val) of function f(x)
    ! f(x) = a0 + a1 * x + a2 * x^2
    
    IMPLICIT NONE
    INTEGER(KIND=Int_kind) ::  i, ns, nf
    REAL(KIND=real_kind), DIMENSION(ns:nf)     :: a0, a1, a2
    REAL(KIND=real_kind), DIMENSION(2, ns:nf)  :: peaks_val  
    INTEGER(KIND=Int_kind), DIMENSION(ns:nf)  :: peaks   
    ! locals 
    REAL(KIND=real_kind)    :: xm, f_xm
    
    ! peaks = 0 IF f has no value > or < than the values of the boudaries
    ! peaks = 1 IF f has MAXima (concave curve)
    ! peaks =-1 IF f has MINima (convex curve)
    ! peaks_val(1,:) MINimum value of f within the interval [0,1]
    ! peaks_val(2,:) MAXimum value of f within the interval [0,1]
    
    peaks=0 
    peaks_val(1,ns:nf)=MIN(a0,a0+a1+a2)
    peaks_val(2,ns:nf)=MAX(a0,a0+a1+a2)  
    
    DO i=ns,nf
       IF ( ABS(a2(i)) > tiny ) THEN
          xm   = -a1(i)/(two*a2(i))            
          f_xm = a0(i) + a1(i)*xm + a2(i)*xm**two
          IF ( xm <= zero .OR. xm >= one ) THEN
             peaks(i)=0 
          ELSEIF (a2(i) > zero) THEN         ! d2f/dx2 > 0 = MINima
             peaks(i)=-1; peaks_val(1,i)= f_xm
          ELSEIF (a2(i) < zero) THEN         ! d2f/dx2 < 0 = MAXima
             peaks(i)=+1; peaks_val(2,i)= f_xm
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE PROPERTIES_QUADRATIC_SPLINE

    !
    !**********************************************
    !
    ! Solves a periodic tridiagonal system
    !
    !**********************************************
    !
    SUBROUTINE tridiag_per(a,b,c,f,jmx)
      !
      ! jmx = dimension of all arrays
      ! a   = sub (lower) diagonal
      ! b   = center diagonal
      ! c   = super (upper) diagonal
      ! f   = right-hand side
      !
      IMPLICIT NONE
      REAL(KIND=real_kind), DIMENSION(jmx), INTENT(IN)   :: a,b,c
      REAL(KIND=real_kind), DIMENSION(jmx), INTENT(INOUT):: f
      INTEGER(KIND=int_kind), INTENT(IN)                :: jmx
      !
      ! Local workspace
      !
      REAL(KIND=real_kind), DIMENSION(jmx)     :: q,s
      REAL(KIND=real_kind)                     :: p,fmx
      INTEGER(KIND=int_kind)                  :: j
      fmx=f(jmx)
      ! forward elimination
      q(1)=-c(1)/b(1)
      f(1)= f(1)/b(1)
      s(1)=-a(1)/b(1)
      DO j=2,jmx
         p    =  one/(b(j)+a(j)*q(j-1))
         q(j) = -c(j)*p
         f(j) =  (f(j)-a(j)*f(j-1))*p
         s(j) = -a(j)*s(j-1)*p
      ENDDO
      ! backward pass
      q(jmx)=zero
      s(jmx)=one
      DO j=jmx-1,1,-1
         s(j) = s(j)+q(j)*s(j+1)
         q(j) = f(j)+q(j)*q(j+1)
      ENDDO
      ! final pass
      f(jmx) = (fmx-c(jmx)*q(1)-a(jmx)*q(jmx-1))/&
               (c(jmx)*s(1)+a(jmx)*s(jmx-1)+b(jmx))
      DO j=1,jmx-1
         f(j)=f(jmx)*s(j)+q(j)
      ENDDO
      RETURN
    END SUBROUTINE tridiag_per
    !
    !**********************************************
    !
    ! Solves a tridiagonal system
    !
    !**********************************************
    !
    SUBROUTINE tridiag(a,b,c,f,jmx)
      !
      ! jmx = dimension of all arrays
      ! a   = sub (lower) diagonal
      ! b   = center diagonal
      ! c   = super (upper) diagonal
      ! f   = right-hand side
      !
      IMPLICIT NONE
      REAL(KIND=real_kind), DIMENSION(jmx), INTENT(IN)   :: a,b
      REAL(KIND=real_kind), DIMENSION(jmx), INTENT(INOUT):: c
      REAL(KIND=real_kind), DIMENSION(jmx), INTENT(INOUT):: f
      INTEGER(KIND=int_kind), INTENT(IN)                :: jmx
      !
      ! Local workspace
      !
      REAL(KIND=real_kind), DIMENSION(jmx)     :: q,s
      REAL(KIND=real_kind)                     :: p,fmx
      INTEGER(KIND=int_kind)                  :: j
      c(jmx)=zero
      ! forward elimination
      q(1)=-c(1)/b(1)
      f(1)= f(1)/b(1)
      DO j=2,jmx
         p    =  one/(b(j)+a(j)*q(j-1))
         q(j) = -c(j)*p
         f(j) =  (f(j)-a(j)*f(j-1))*p
      ENDDO
      ! backward pass
      DO j=jmx-1,1,-1
         f(j)=f(j)+q(j)*f(j+1)
      ENDDO
      RETURN
    END SUBROUTINE tridiag
end module





module remap_lauritzen
use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
contains
SUBROUTINE verremap2(plev ,plevmodel, parg,klev,pres , ireconstruct,pres_min,pres_max,check_global)

!  USE pcmpsm_no_pbc
!  USE constants
  use kinds, only              : real_kind, int_kind
  use pcmpsm_no_pbc, only : quadratic_spline,quadratic_spline_monotone,&   ! _EXTERNAL  
         cubic_parameters,cubic_monotone
  use parallel_mod, only : abortmp
  IMPLICIT NONE
  
  ! ireconstruct = 0: piecewise cubic method
  ! ireconstruct = 1: piecewise cubic method with UK Met Office monotonoicity constraints
  ! ireconstruct = 2: quadratic splines 
  ! ireconstruct = 3: quadratic splines with UK Met Office monotonoicity constraints
  !
  !
  !
  INTEGER(KIND=Int_kind), INTENT(IN) :: klev, ireconstruct
  REAL(KIND=real_kind)& 
       plevmodel(klev+1),&         !location of pressure levels in terms of ps and
                                      !and hybrid coefficients
          plev(klev+1),  &              !location of pressure levels implied by advection scheme
          parg(klev),    &              !variable to be vertiCALLy remapped
          pres(klev),    &              !remaped field
          pres_min,      &
          pres_max 
  !
  ! LOCAL WORKSPACE
  !
  INTEGER(KIND=Int_kind) :: jk, ilev, jl, piecewise, check_global,&
  !
  !    IN CASE OF SURFACE PRESSURE MISMATCH WHERE MODEL LEVELS ARE .GE. THE
  !     SURFACE PRESSURE IMPLIED BY TRANSPORT SCHEME JLMS REFERS TO
  !     THE MINIMUM INDEX FOR WHICH THIS IS THE CASE
  !
  jlms       ,&
  zkr(klev+1),&          ! in which eulerian cell is floor of cell jk located
  itop, ibot ,&          ! index of cell in which top/bottom wall is located
  jsubz
  REAL(KIND=real_kind) :: & 
        zgam(klev+1)&         ! dimensionless distance from floor of cell to eulerian cell floor
       ,zdp          &        ! dummy for pressure level thickness
       ,za0(klev)    &       ! coefficients for parabolaes
       ,za1(klev)   &
       ,za2(klev)  & 
       ,za3(klev)  & 
       ,zhdp(klev)     &      !pressure level thicknesses
       ,zaccintegerb    &     !entire cell mass
       ,zacctop,zaccbot  &    !mass accumulated towards top of cell jk
       ,zpsmodel,zpscisl  &   !surface pressure for cisl scheme and model
       ,ztmp, zmasstoadd&
       ,zeps&
       !         ,zmass0,zmass1&
       ,zarg(klev)&
       ,zmass0,zmass1

  real(kind=real_kind) :: zero=0

  LOGICAL LMS               !LOWER LEVEL PRESSURE MISMATCH
  IF (ABS(plevMODEL(KLEV+1)-plev(KLEV+1)).GE.0.000001) THEN
     WRITE(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
     WRITE(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
     WRITE(6,*) 'DATA FOR MODEL LEVELS'
     WRITE(6,*) 'PLEVMODEL=',plevMODEL(KLEV+1)
     WRITE(6,*) 'PLEV     =',plev     (KLEV+1)
     WRITE(6,*) 'DIFF     =',plevMODEL(KLEV+1)-plev(KLEV+1)
!     c         CALL ABORT
  ENDIF
!
!     INITIALIZE BEFORE INTEGRATION
  pres     = 0.0
  zpsmodel = plevmodel(klev+1)
  zpscisl  = plev(klev+1)
  
  zkr  = 99
  ilev = 2
  jlms = klev+1
  lms  = .false.
  !
  ! upper boundary is known
  !
  zkr(1)       = 1
  zgam(1)      = 0.0
  zkr(klev+1)  = klev
  zgam(klev+1) = 1.0
  
!
!  call t_startf('verremap-loopa')
  jlloop: DO jl = 2,klev
     !
     !  find nearest level jk, where plev(jk).ge.plevmodel(jl)
     !
     DO jk = ilev,klev+1
        IF (plev(jk).ge.plevmodel(jl)) THEN
           ilev      = jk
           zkr(jl)   = jk-1
           !
           !     dimensionless distance from celling of cell jk in which zplev is located
           !     cell jk is bounded by plev(jk) and plev(jk+1)
           !
           !           ======      ================ plev(1) = zplevmodel(1) ; zgam(1) = 0   ---
           !                                                                                 | zgam(2)
           !           cell 1      ---------------- zplevmodel(2)                           ---
           !                                                                    
           !           ======      ================ plev(2)                  
           !           cell 2
           !           ======      ================ plev(3)                   ---         --- 
           !                                                                   | zgam(2)   |
           !                       ---------------- zplevmodel(3)             ---          | zgam(3)
           !                       ---------------- zplevmodel(4)                         ---
           !                                                          
           !           ======      ================ plev(4)           
           !           cell 3              .
           !                               .
           !                               .
           !           ======      ================ plev(nlev) = zplevmodel(nlev)
           !
           !
           
           zdp        = plev(jk)-plev(jk-1)
           zgam(jl)   = (plevmodel(jl)-plev(jk-1))/zdp
           CYCLE jlloop
        ENDIF
     ENDDO
  ENDDO jlloop
!  call t_stopf('verremap-loopa')
  !
  !     ****************************************************
  !     actual remapping
  !     ****************************************************
  !
!  call t_startf('verremap-loopb')
  zmass0=0
  DO jk=1,klev
     !     compute levels thicknesses for levels implied by advection scheme
     zhdp(jk) = plev(jk+1)-plev(jk)
     zarg(jk) = parg(jk)/zhdp(jk)
     zmass0  = zmass0+zhdp(jk)*zarg(jk)
  ENDDO
  IF (ireconstruct <2) THEN
     !     calculate pcm-coefficients
     piecewise    = 3
     CALL cubic_parameters(plev,parg,zhdp,za0,za1,za2,za3,piecewise,1,klev)
     IF (ireconstruct==1)&
          CALL cubic_monotone(za0,za1,za2,za3,zarg, 1,klev,pres_min, pres_max, check_global, piecewise)
  ELSE IF (ireconstruct<4) THEN
     za3 = 0
     CALL quadratic_spline(zhdp,parg,za0,za1,za2,klev)
     IF (ireconstruct==3) &
          CALL quadratic_spline_monotone(za0,za1,za2,zarg,zhdp,1,klev,pres_min,pres_max,check_global)
  ELSE
     call abortmp("ireconstruct out of range")
  ENDIF
!  call t_stopf('verremap-loopb')
      
  !-----start iteration from top to bottom of atmosphere------------------------------            
  zaccintegerb = 0
  itop = 1
  zacctop = 0.0

!  call t_startf('verremap-loopc')
  DO jk = 1,jlms-1
     ibot = zkr(jk+1)
     DO jsubz=itop,ibot-1,1
        zaccintegerb = zaccintegerb + parg(jsubz)
     ENDDO
     CALL int_cubic(zero,zgam(jk+1),za0(ibot),za1(ibot),za2(ibot),za3(ibot),ztmp)
     zaccbot = zaccintegerb + ztmp*zhdp(ibot)
     pres(jk) = (zaccbot-zacctop)
     !
     !	convert to mixing ratio
     !
     !     $	            /(plevmodel(jk+1)-plevmodel(jk))
     !     prepare for next iteration
     zacctop        = zaccbot
     itop           = ibot
  ENDDO
!  call t_stopf('verremap-loopc')

  RETURN

  !
  !
  !	convert to mixing ratio
  !
  zmass1 = 0
  DO jk=1,klev
     zdp = plevmodel(jk+1)-plevmodel(jk)
     pres(jk) = pres(jk)/zdp
     zmass1 = zmass1+pres(jk)*zdp
  ENDDO
!  WRITE(*,*) "mass change in verremap2",zmass1-zmass0

END SUBROUTINE verremap2


!
!***********************
! Intergral of cubic
!***********************
!
SUBROUTINE int_cubic(l,r,za0,za1,za2,za3,mass)
!  USE constants
  USE kinds, only : real_kind
  use parallel_mod, only : abortmp
  IMPLICIT NONE
  REAL(KIND=real_kind), INTENT(IN) :: l,r             !left/right bound \in [0,1]
  REAL(KIND=real_kind), INTENT(IN) :: za0,za1,za2,za3 !coefficients
  REAL(KIND=real_kind), INTENT(OUT):: mass  !RETURN: INTEGREAL    
  !
  ! local workspace
  !
  REAL(KIND=real_kind) zr, zl
  REAL(KIND=real_kind) :: tiny = 1e-12
  !
  zr = r
  zl = l
  IF (ABS(zl).LT.tiny)     zl = 0
  IF (ABS(zr).LT.tiny)     zr = 0
  IF (ABS(zr-1).LT.tiny) zr = 1
  IF (ABS(zl-1).LT.tiny) zl = 1
  
  IF (zl.LT.0.OR.zl.GT.1.OR.&            !DBG
       zr.LT.0.OR.zr.GT.1) THEN          !DBG
     WRITE(*,*) 'r or l not in [0:1]'         !DBG
     WRITE(*,*) 'r,l=',zr,zl                  !DBG
     call abortmp(' ')                                     !DBG
  ENDIF                                       !DBG
  IF (zl.GT.zr) THEN                          !DBG
     call abortmp('r<l')                         !DBG
  ENDIF                                       !DBG
  IF (ABS(zr - zl) .LT. tiny) THEN
     mass = 0
  ELSE
     mass    =                                            &
          za0*(zr-zl)+(za1/2)*(zr**2-zl**2)+             &
          (za2/3)*(zr**3-zl**3)+(za3/4)*(zr**4-zl**4)
  ENDIF
  RETURN
END SUBROUTINE int_cubic
end module









module prim_advection_mod
!
! two formulations.  both are conservative
! u grad Q formulation:
!
!    d/dt[ Q] +  U grad Q  +  eta_dot dp/dn dQ/dp  = 0
!                            ( eta_dot dQ/dn )
!
!    d/dt[ dp/dn ] = div( dp/dn U ) + d/dn ( eta_dot dp/dn )
!
! total divergence formulation:
!    d/dt[dp/dn Q] +  div( U dp/dn Q ) + d/dn ( eta_dot dp/dn Q ) = 0
!
! for convience, rewrite this as dp Q:  (since dn does not depend on time or the horizonal): 
! equation is now:
!    d/dt[dp Q] +  div( U dp Q ) + d( eta_dot_dpdn Q ) = 0
!
!  
  use kinds, only              : real_kind
  use dimensions_mod, only     : nlev, nlevp, nv, qsize
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, cp, cpwater_vapor
  use derivative_mod, only     : gradient, vorticity, gradient_wk, derivative_t, divergence, &
                                 laplace_sphere_wk, gradient_sphere, divergence_sphere
  use element_mod, only        : element_t
  use filter_mod, only         : filter_t, filter_V
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, smooth
  use prim_si_mod, only        : preq_pressure
  use diffusion_mod, only      : scalar_diffusion, diffusion_init
  use control_mod, only        : integration, test_case, filter_freq_advection,  hypervis_order, &
        statefreq, moisture, TRACERADV_TOTAL_DIVERGENCE, TRACERADV_UGRADQ, &
        tracer_advection_formulation, prescribed_wind, rk_stage_user, nu_q, &
        compute_mean_flux, limiter_option
  use edge_mod, only           : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, initedgebuffer, edgevunpackmin
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use viscosity_mod, only      : biharmonic_wk_scalar, biharmonic_wk_scalar_minmax


  implicit none
  
  private  

  public :: Prim_Advec_Init, Prim_Advec_Tracers_remap_rk2, Prim_Advec_Tracers_lf
  type (EdgeBuffer_t) :: edgeAdv, edgeAdvQ3, edgeAdv_p1

  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  integer,parameter :: USEconsistent = -1
  integer,parameter :: USEv = 0
  integer,parameter :: USEvn0 = 1
  integer,parameter :: USEave = 2
  integer,parameter :: USEvstar = 3


contains
  subroutine Prim_Advec_Init()
    use dimensions_mod, only : nlev, qsize

    call initEdgeBuffer(edgeAdv,qsize*nlev)
    call initEdgeBuffer(edgeAdv_p1,qsize*nlev + nlev) 
    call initEdgeBuffer(edgeAdvQ3,qsize*nlev*3)  ! Qtens,Qmin, Qmax

  end subroutine Prim_Advec_Init



!
! forward-in-time 2 level vertically lagrangian step
!  this code takes a lagrangian step in the horizontal 
! (complete with DSS), and then applies a vertical remap
!
! This routine may use dynamics fields at timelevel np1
! In addition, other fields are required, which have to be 
! explicitly saved by the dynamics:  (in elem(ie)%derived struct)
!
! Fields required from dynamics: (in 
!    omega_p   it will be DSS'd here, for later use by CAM physics
!              we DSS omega here because it can be done for "free"
!    Consistent mass/tracer-mass advection (used if subcycling turned on)
!       dp()   dp at timelevel n0
!       vn0()  mean flux  < U dp > going from n0 to np1
!    Non-consistent scheme used with leapfrog dynamics, no subcycling
!       vn0()           U at timelevel n0 
!       eta_dot_dpdn()  mean vertical velocity
! 
!
! two stage:  
!    Euler step from t -> t+1
!    Euler step from t+1 -> t+2
!    u(t) = (u(t)+u(t+2))/2
!
! 3 stage
!    Euler step from t     -> t+.5
!    Euler step from t+.5  -> t+1.0
!    Euler step from t+1.0 -> t+1.5
!    u(t) = u(t)/3 + u(t+2)*2/3
!

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine Prim_Advec_Tracers_remap_rk2(elem, deriv,hvcoord,flt,hybrid,&
        dt,tl,nets,nete, compute_diagnostics)
    use perf_mod, only : t_startf, t_stopf            ! _EXTERNAL
    use derivative_mod, only : divergence_sphere

    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (filter_t)                   :: flt
    type (hybrid_t),     intent(in):: hybrid
    type (TimeLevel_t)                :: tl

    real(kind=real_kind) , intent(in) :: dt
    integer,intent(in)                :: nets,nete

    logical,intent(in)                :: compute_diagnostics
    real (kind=real_kind), dimension(nv,nv,2)    :: gradQ


    real (kind=real_kind), dimension(nv,nv,nlev) :: Q_vadv 
    real (kind=real_kind), dimension(nv,nv,nlev)    :: dp_star
    real (kind=real_kind), dimension(nv,nv,nlev)    :: dp_np1
   
    integer :: i,j,k,l,ie,q,nmin
    integer :: n0,np1,nfilt,rkstage,rhs_multiplier

    real(kind=real_kind) :: notreliable



    call t_startf('prim_advec_tracers_remap_rk2')
    n0  = tl%n0
    np1 = tl%np1
    rkstage=2

    if (rk_stage_user==3) rkstage=3

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RK2 2D advection step
    ! note: stage 1 we take the oppertunity to DSS eta_dot_dpdn
    ! note: stage 2 we take the oppertunity to DSS omega
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(compute_mean_flux==1)then
      ! mean flux and dp(n0) was computed by dynamics
      ! use these for consistent advection (preserve Q=1)
      ! derived%vn0             =  mean horiz. flux:   U*dp
      ! derived%dp              =  dp at n0.  needed for remap.  
      ! derived%omega_p         =  advection code will DSS this for the physics, but otherwise 
      !                            it is not needed 
      ! Also: save a copy of div(U dp) in derived%div(:,:,:,1), which will be DSS'd 
      !       and a DSS'ed version stored in derived%div(:,:,:,2)
      do ie=nets,nete
	do k=1,nlev
	  ! div( U dp Q), 
	    gradQ(:,:,1)=elem(ie)%derived%vn0(:,:,1,k)
	    gradQ(:,:,2)=elem(ie)%derived%vn0(:,:,2,k)
	    elem(ie)%derived%divdp(:,:,k) = divergence_sphere(gradQ,deriv,elem(ie))
	enddo
	elem(ie)%derived%divdp_proj(:,:,:) = elem(ie)%derived%divdp(:,:,:)
      enddo

      if (rkstage .eq. 2) then
      !rhs_multiplier is for obtaining dp_tracers at each stage:
      !dp_tracers(stage) = dp - rhs_multiplier*dt*divdp_proj
	  rhs_multiplier = 0
	  call euler_step(np1,n0,dt,elem,hvcoord,hybrid,deriv,nets,nete,&
		compute_diagnostics,USEconsistent,DSSdiv_vdp_ave,rhs_multiplier)

	  rhs_multiplier = 1
	  call euler_step(np1,np1,dt,elem,hvcoord,hybrid,deriv,nets,nete,&
		.false.,USEconsistent,DSSomega,rhs_multiplier)
      else
      !rhs_multiplier is for obtaining dp_tracers at each stage:
      !dp_tracers(stage) = dp - rhs_multiplier*dt*divdp_proj
	  rhs_multiplier = 0
          call euler_step(np1,n0,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,&
		compute_diagnostics,USEconsistent,DSSdiv_vdp_ave,rhs_multiplier)

	  rhs_multiplier = 1
	  call euler_step(np1,np1,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,&
		.false.,USEconsistent,DSSomega,rhs_multiplier)

	  rhs_multiplier = 2
	  call euler_step(np1,np1,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,&
		.false.,USEconsistent,DSSno_var,rhs_multiplier)
      endif

    else
      ! non-consistent code uses (computed in the dynamics)
      ! derived%vn0             =  horiz velocity at n0
      ! state%v                 =  horiz velocity at np1
      ! derived%eta_dot_dpdn    =  mean vertical velocity (used for remap)
      ! derived%omega_p         =  advection code will DSS this for the physics, but otherwise 
      !                            it is not needed 

      if (rkstage .eq. 2) then

	! STAGE 1 uses U(t) which was saved by dynamics in elem%derived%vn0():
	call euler_step(np1,n0,dt,elem,hvcoord,hybrid,deriv,nets,nete,&
	      compute_diagnostics,USEvn0,DSSeta,0)

	! STAGE 2 needs U(t+1) on dp_star(t=1) surface
	call euler_step(np1,np1,dt,elem,hvcoord,hybrid,deriv,nets,nete,&
	      .false.,USEv,DSSomega,0)

      else
	! STAGE 1 needs U(t) which was saved in elem%state%vn0():
	call euler_step(np1,n0,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,&
	    compute_diagnostics,USEvn0,DSSeta,0)
	
	! compute U* = U(t+1) on dp_star levels, store in elem%state%vstar:
	call remap_velocity(np1,dt,elem,hvcoord,hybrid,deriv,nets,nete)

	! STAGE 2: needs U(t+.5) to second order.  Use average of vn0 and vstar:
	call euler_step(np1,np1,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,&
	    .false.,USEave,DSSomega,0)

	! STAGE 3: needs U(t+1) on dp_star(t).  Use vstar:
	call euler_step(np1,np1,dt/2,elem,hvcoord,hybrid,deriv,nets,nete,&
	    .false.,USEvstar,DSSno_var,0)
      endif
    endif



    do ie=nets,nete
      ! finish 2D advection step:
#ifdef ZEROHORZ
        elem(ie)%state%Qdp(:,:,:,:,np1)=elem(ie)%state%Qdp(:,:,:,:,n0) ! ignore 2D step
#else
       !  take average of t and t+2 results:

        elem(ie)%state%Qdp(:,:,:,:,np1) = (elem(ie)%state%Qdp(:,:,:,:,n0)+&
            (rkstage-1)*elem(ie)%state%Qdp(:,:,:,:,np1) )/rkstage
#endif

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !  VERTICAL and FORCING
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! compute dp_star and dp_np1
        if(compute_mean_flux==1 .and. prescribed_wind==0)then
	      ! consistent advection.  dp_star is horizontal advection of dp
	      do k=1,nlev
		  dp_np1(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
		      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
		  dp_star(:,:,k) =  elem(ie)%derived%dp(:,:,k)-dt*elem(ie)%derived%divdp_proj(:,:,k)
	      enddo
	else
	      ! otherwise, use mean eta_dot_dpdn from dynamics
	      ! also must be used if dp() and eta_dot_dpdn() prescribed
	      do k=1,nlev
		  dp_np1(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
		      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
		  ! compute dp_star from eta_dot_dpdn(t+.5)
		  dp_star(:,:,k) = dp_np1(:,:,k) + &
		      dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -  &
		      elem(ie)%derived%eta_dot_dpdn(:,:,k) ) 
	      enddo
	endif

#ifdef ZEROVERT
       dp_star=dp_np1  ! ignore the vertical motion
#endif

       do q=1,qsize
          ! remap Q.  also return Q_vadv for diagnostics
          call preq_vertadv_remap(elem(ie)%state%Qdp(:,:,:,q,np1),&
               dp_star,dp_np1,dt,hvcoord,Q_vadv,.true.)
#ifdef ENERGY_DIAGNOSTICS
          if (compute_diagnostics .and. q==1) then
             ! IEvert1_wet():  (Cpv-Cp) T Qdp_vadv  (Q equation)
             ! IEhorz1_wet():  (Cpv-Cp) T Qdp_hadv  (Q equation)
             do k=1,nlev
                elem(ie)%accum%IEvert1_wet(:,:) = elem(ie)%accum%IEvert1_wet(:,:) +&
                  (Cpwater_vapor-Cp)*elem(ie)%state%T(:,:,k,n0)*Q_vadv(:,:,k)
             enddo
          endif
#endif
       enddo
    enddo



#ifndef ZEROHORZ
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Dissipation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! limiter=8 applies dissipation in RHS, not here
    if ((limiter_option /= 8).and.(limiter_option /= 84)) then
       call advance_hypervis_scalar(edgeadv,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt)
    endif
#endif

    ! update Q from Qdp
    do ie=nets,nete
       do q=1,qsize       
       do k=1,nlev
	    dp_np1(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
            elem(ie)%state%Q(:,:,k,q,np1)=elem(ie)%state%Qdp(:,:,k,q,np1)/dp_np1(:,:,k) 
       enddo
       enddo
    enddo


    call t_stopf('prim_advec_tracers_remap_rk2')
  end subroutine prim_advec_tracers_remap_rk2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine Prim_Advec_Tracers_lf(elem, deriv,hvcoord,flt,hybrid,dt,tl,nets,nete, compute_diagnostics)
    use perf_mod, only : t_startf, t_stopf              ! _EXTERNAL
    implicit none
    type (element_t), intent(inout)   :: elem(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (filter_t)                   :: flt

    type (hybrid_t),     intent(in):: hybrid

    real(kind=real_kind) , intent(in) :: dt
    type (TimeLevel_t)                :: tl

    integer,intent(in)                :: nets,nete
    logical,intent(in)                :: compute_diagnostics

    ! local
    real(kind=real_kind) :: dp, pmid, dt2
    integer :: i,j,k,l,ie,q
    integer :: n0,nm1,np1,nfilt, nstep


    call t_startf('prim_advec_tracers_lf')

    n0  = tl%n0
    nm1 = tl%nm1
    np1 = tl%np1

    nfilt = n0
    nstep = tl%nstep
    dt2 = 2 * dt


    if (filter_freq_advection>0) then
    if(nstep>0 .and.  modulo(nstep,filter_freq_advection)==0) then
       do ie=nets,nete
          do q=1,qsize	
             do k=1,nlev
                call filter_V(elem(ie)%state%Q(:,:,k,q,nfilt),flt)
                elem(ie)%state%Q(:,:,k,q,nfilt) = elem(ie)%mv(:,:)*elem(ie)%state%Q(:,:,k,q,nfilt)
             end do
          end do
          call edgeVpack(edgeadv,elem(ie)%state%Q(:,:,:,:,nfilt),nlev*qsize,0,elem(ie)%desc)
       end do
       call bndry_exchangeV(hybrid,edgeadv)
       do ie=nets,nete
          call edgeVunpack(edgeadv,elem(ie)%state%Q(:,:,:,:,nfilt),nlev*qsize,0,elem(ie)%desc)
          do q=1,qsize
             do k=1,nlev
                elem(ie)%state%Q(:,:,k,q,nfilt) = elem(ie)%rmv(:,:)*elem(ie)%state%Q(:,:,k,q,nfilt)
             enddo
          end do
       end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
    end if
    end if


    if (nstep==0) then
       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
       call compute_and_apply_rhs(np1,n0,n0,dt/2,elem,hvcoord,hybrid,&
            deriv,nets,nete,.false.)
       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
       call compute_and_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics)
    else
       call compute_and_apply_rhs(np1,nm1,n0,dt2,elem,hvcoord,hybrid,&
            deriv,nets,nete,compute_diagnostics)
    endif


    if (hypervis_order==0) then
       call scalar_diffusion(elem, nets,nete,np1,deriv,dt2,hybrid)
    else
       ! hypervis_order==1:  weak form laplacian  (only 1 DSS)
       ! hypervis_order==2:  weak form biharmonic (2 DSS's) 
       call advance_hypervis_scalar_lf(edgeadv,elem,hvcoord,hybrid,deriv,np1,n0,nets,nete,dt2)
    endif

    if (nu_q>0) then
    ! if nu_q=0, we are running an inviscid test, skip fixer
    !
    ! apply negative Q fixer
    !
    do ie=nets,nete
       do q=1,qsize
          do k=1,nlev
             elem(ie)%state%Q(:,:,k,q,np1) = elem(ie)%spheremv(:,:)*elem(ie)%state%Q(:,:,k,q,np1)
          enddo

          ! limiter3d_noncon: no negative values, even if mass is added
          call limiter3d_noncon(elem(ie)%state%Q(:,:,:,q,np1),&
              elem(ie)%state%ps_v(:,:,n0),&
              hvcoord,elem(ie)%accum%mass_added(q))

       end do
       call edgeVpack(edgeadv,elem(ie)%state%Q(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
    end do
    call bndry_exchangeV(hybrid,edgeadv)
    do ie=nets,nete
       call edgeVunpack(edgeadv,elem(ie)%state%Q(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
       do q=1,qsize
          do k=1,nlev
             elem(ie)%state%Q(:,:,k,q,np1) = elem(ie)%rspheremv(:,:)*elem(ie)%state%Q(:,:,k,q,np1)
          enddo
       end do
    end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
    endif

    call t_stopf('prim_advec_tracers_lf')

    if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
       do ie=nets,nete
          do q=1,qsize
             do k=1,nlev
                do j=1,nv
                   do i=1,nv
                      ! timestep was done in Q.  copy over to Qdp:                                              
                      elem(ie)%state%Qdp(i,j,k,q,np1)=elem(ie)%state%Q(i,j,k,q,np1)
                      ! recompute Q from dpQ for consistency                                                    
                      dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                           ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)
                      elem(ie)%state%Q(i,j,k,q,np1) =elem(ie)%state%Qdp(i,j,k,q,np1)/dp
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
       

  end subroutine Prim_Advec_Tracers_lf

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine compute_and_apply_rhs(np1,nm1,n0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,compute_diagnostics)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine is normally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
  !
  ! Combining the RHS and DSS pack operation in one routine 
  ! allows us to fuse these two loops for more cache reuse
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : nv, np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use hybvcoord_mod, only : hvcoord_t

  implicit none
  integer :: np1,nm1,n0,nets,nete
  real (kind=real_kind) :: dt2
  logical  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv

  ! local
  real (kind=real_kind) ::  qtens(nv,nv,nlev)
  real (kind=real_kind) ::  Q_vadv(nv,nv,nlev)
  real (kind=real_kind) ::  rpdel(nv,nv,nlev)
  real (kind=real_kind) ::  gradQ(nv,nv,2)
  real (kind=real_kind) ::  divdp(nv,nv)
  real (kind=real_kind) ::  v1,v2
  integer :: i,j,k,ie,q

  logical ::  use_explicit_eta_dot=.true. ! reuse from dynamics or recompute?


  do ie=nets,nete
     do q=1,qsize
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   2D contribution
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
           do k=1,nlev
              ! div( U dp Q),                                                                                 
              gradQ(:,:,1)=elem(ie)%state%v(:,:,1,k,n0)*elem(ie)%state%Qdp(:,:,k,q,n0)
              gradQ(:,:,2)=elem(ie)%state%v(:,:,2,k,n0)*elem(ie)%state%Qdp(:,:,k,q,n0)
              divdp = divergence_sphere(gradQ,deriv,elem(ie))
              do j=1,nv
                 do i=1,nv
                    qtens(i,j,k)=-divdp(i,j)
                 enddo
              enddo
           enddo
        else
           !   UGRADQ formulation
           do k=1,nlev
              gradQ = gradient_sphere(elem(ie)%state%Q(:,:,k,q,n0),deriv,elem(ie))
              do j=1,nv	
                 do i=1,nv
                    v1    = elem(ie)%state%v(i,j,1,k,n0)
                    v2    = elem(ie)%state%v(i,j,2,k,n0)
                    Qtens(i,j,k) = -( v1*gradQ(i,j,1) + v2*gradQ(i,j,2)  )
                 enddo
              enddo
           enddo
        endif
        


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! vertical advection
        ! evaluate at np1 for time-split scheme
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
           call preq_vertadv_dpQ(elem(ie)%state%Q(:,:,:,q,n0),elem(ie)%derived%eta_dot_dpdn,Q_vadv)
           ! advance in time, into Q, apply mass matrix
           do k=1,nlev
              elem(ie)%state%Q(:,:,k,q,np1) = elem(ie)%spheremv(:,:)*&
                   ( elem(ie)%state%Qdp(:,:,k,q,nm1)  + &
                   dt2*(qtens(:,:,k)-Q_vadv(:,:,k)) )
           enddo
        else
           if (use_explicit_eta_dot) then
              ! vertical advection term, using eta_dot_dpdn from advection timestep
              do k=1,nlev
                 do j=1,nv
                    do i=1,nv
                       v1 = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(i,j,n0)
                       v2 = hvcoord%hyai(k+1)*hvcoord%ps0 + hvcoord%hybi(k+1)*elem(ie)%state%ps_v(i,j,n0)
                       rpdel(i,j,k) = 1.0D0/(v2-v1)
                    end do
                 end do
              enddo
              call preq_vertadvQ(elem(ie)%state%Q(:,:,:,q,n0),elem(ie)%derived%eta_dot_dpdn,rpdel,Q_vadv)
           else
              ! recompute eta_dot_dpdn using velocity at level n0
              ! and elem%derived%grad_lnps
              call preq_impsysQ(elem(ie),hvcoord,np1,n0,nm1,elem(ie)%state%Q(:,:,:,q,n0),Q_vadv)
           endif
           ! advance in time, apply mass matrix
           do k=1,nlev
              elem(ie)%state%Q(:,:,k,q,np1) = elem(ie)%spheremv(:,:)*&
                   ( elem(ie)%state%Q(:,:,k,q,nm1)  + &
                   dt2*(qtens(:,:,k)-Q_vadv(:,:,k)) )
           enddo
        endif
        
        if (nu_q>0) then
           ! if nu_q=0, we are running an inviscid test, skip fixer
           call limiter2d_zero(elem(ie)%state%Q(:,:,:,q,np1),&
             elem(ie)%state%ps_v(:,:,n0),hvcoord)
        endif
        
     end do
     call edgeVpack(edgeadv,elem(ie)%state%Q(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
  end do
  call bndry_exchangeV(hybrid,edgeadv)
  
  do ie=nets,nete
     call edgeVunpack(edgeadv,elem(ie)%state%Q(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
     
     do q=1,qsize
        do k=1,nlev
           elem(ie)%state%Q(:,:,k,q,np1) = elem(ie)%rspheremv(:,:)*elem(ie)%state%Q(:,:,k,q,np1)
        enddo
     end do
  end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
end subroutine compute_and_apply_rhs

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine remap_velocity(nl,dt,elem,hvcoord,hybrid,deriv,nets,nete)
  ! 
  ! compute u*(t+1): velocity at t+1 from a Lagrange step 
  ! store result in elem%derived%vstar
  !
  ! NOTE: routine uses elem%derived%eta_dot_dpdn, which is assumed
  ! already computed and DSS'd
  !
  ! Leapfrog: 
  !    u(t+1) = u(t-1) - 2dt*2Dterms(t) - 2dt*eta_dot_dpdn d/dn U  = 0
  ! Leapfrog Lagrangian would look like this:
  !    u*(t+1) = u*(t-1) - 2dt*2Dterms(t) = 0
  ! Which is satisfied if we take:
  !    u*(t+1) = u(t+1) + dt eta_dot_dpdn d/dn U
  !    u*(t-1) = u(t-1) - dt eta_dot_dpdn d/dn U
  !   
  ! Two methods to compute this:  
  ! REMAP:
  !   Reference surface:  dp 
  !   Lagrangian surface: dp_star
  !      dp_star = dp  + dt d/dn[eta_dot_dpdn] 
  !
  !   Remap U(nl) given on dp(nl) to Ustar given on dp_star
  !
  ! ADVECTION:
  !   advect U(nl) from dp to dp_star using DSS'd v_vadv(), which
  !   we recompute from the DSS'd eta_dot_dpdn()                                    
  !   use dynamics subroutine which computes v_adv() 
  !   not yet coded.
  !
  use kinds, only : real_kind
  use dimensions_mod, only : nv, np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use hybvcoord_mod, only : hvcoord_t

  implicit none
  integer :: nl,nets,nete
  real (kind=real_kind), intent(in)  :: dt

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv

  ! local
  real(kind=real_kind), dimension(nv,nv,nlev) :: work
  real(kind=real_kind), dimension(nv,nv,nlev) :: dp
  real(kind=real_kind), dimension(nv,nv,nlev) :: dp_star
  real(kind=real_kind), dimension(nv,nv,nlev) :: Ustar
  real(kind=real_kind), dimension(nv,nv,nlev) :: Vstar
  integer :: ie,i,j,k


  do ie=nets,nete
     ! remap U(nl) from dp to dp_star 
     do k=1,nlev
        dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,nl)
        ! compute dp_star from eta_dot_dpdn(t+.5)
        dp_star(:,:,k) = dp(:,:,k) + &
             dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -  &
             elem(ie)%derived%eta_dot_dpdn(:,:,k) ) 
     enddo
     Ustar(:,:,:) = elem(ie)%state%v(:,:,1,:,nl)*dp(:,:,:)
     Vstar(:,:,:) = elem(ie)%state%v(:,:,2,:,nl)*dp(:,:,:)
     call preq_vertadv_remap(Ustar,dp,dp_star,dt,hvcoord,work,.false.)
     call preq_vertadv_remap(Vstar,dp,dp_star,dt,hvcoord,work,.false.)
     elem(ie)%derived%vstar(:,:,1,:) = Ustar(:,:,:)/dp_star(:,:,:)
     elem(ie)%derived%vstar(:,:,2,:) = Vstar(:,:,:)/dp_star(:,:,:)
  enddo
  end subroutine remap_velocity

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine euler_step(np1,n0,dt,elem,hvcoord,hybrid,deriv,nets,nete,&
      compute_diagnostics,Uopt,DSSopt,rhs_multiplier)
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.  
  !
  ! Uopt = USEv      use elem%state%v(n0)  
  ! Uopt = USEvn0    use elem%derived%vn0   (precomputed by calling program)
  ! Uopt = USEvstar   use elem%derived%vstar (precomputed by calling program)
  ! Uopt = USEave    average of above vstar and vn0
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : nv, np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod, only : edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use hybvcoord_mod, only : hvcoord_t
  use control_mod, only : limiter_option

  implicit none
  integer :: np1,nm1,n0,nets,nete,Uopt,DSSopt,rhs_multiplier
  real (kind=real_kind), intent(in)  :: dt
  logical  :: compute_diagnostics

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv

  ! local
  real (kind=real_kind), dimension(nv,nv)    :: divdp
  real (kind=real_kind), dimension(nv,nv,2)    :: gradQ
  real(kind=real_kind), dimension(nv,nv,nlev) :: Qtens
  real(kind=real_kind), dimension(nv,nv,nlev) :: dp,dp_star
  real(kind=real_kind), dimension(nv,nv,2,nlev) :: Vstar
  real (kind=real_kind), pointer, dimension(:,:,:)   :: DSSvar
  real(kind=real_kind), dimension(nlev,qsize,nets:nete) :: qmin,qmax
  real(kind=real_kind), dimension(nv,nv,nlev,qsize,nets:nete) :: Qtens_biharmonic
  real(kind=real_kind) :: dp0
  integer :: ie,q,i,j,k
  real (kind=real_kind) :: notreliable

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute tendencies from biharmonic diffusion, 
  !   if this is to be done on the RHS instead of
  !   time-split
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if((limiter_option == 8).or.(limiter_option == 84))then
     do ie=nets,nete
        ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
        dp(:,:,:) = elem(ie)%derived%dp(:,:,:) - &
             rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,:) 
        do q=1,qsize
           Qtens_biharmonic(:,:,:,q,ie) = elem(ie)%state%Qdp(:,:,:,q,n0)/dp(:,:,:)
        enddo
     enddo
     call biharmonic_wk_scalar_minmax(elem,qtens_biharmonic,deriv,edgeAdvQ3,hybrid,&
          nets,nete,qmin,qmax)

     ! output has mass matrix already applied.  un-apply since we want to apply
     ! the mass matrix later
     do ie=nets,nete
        do q=1,qsize
           do k=1,nlev
              dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                   ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
              qtens_biharmonic(:,:,k,q,ie) =  -dt*nu_q*dp0*Qtens_biharmonic(:,:,k,q,ie)&
                   / elem(ie)%spheremv(:,:)
           enddo
        enddo
     enddo
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do ie=nets,nete

     ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
     ! all zero so we only have to DSS 1:nlev
     if ( DSSopt == DSSeta) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
     if ( DSSopt == DSSomega) DSSvar => elem(ie)%derived%omega_p(:,:,:)
     if ( DSSopt == DSSdiv_vdp_ave) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)


     ! dp for RK stage initial time 
     dp(:,:,:) = elem(ie)%derived%dp(:,:,:) - &
          rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,:) 

     ! un-DSS'ed dp for RK stage end time
     ! dp_star = dp - dt*elem(ie)%derived%divdp(:,:,:)   computed below, if needed 

     ! Compute velocity used to advance Qdp 
     if (Uopt==USEv) then
       Vstar = elem(ie)%state%v(:,:,:,:,n0)
     else if (Uopt==USEvn0) then
       Vstar = elem(ie)%derived%vn0(:,:,:,:)
     else if (Uopt==USEvstar) then
       Vstar = elem(ie)%derived%vstar(:,:,:,:)
     else if (Uopt==USEave) then
       Vstar(:,:,:,:) = (elem(ie)%derived%vstar(:,:,:,:) + elem(ie)%derived%vn0(:,:,:,:))/2
     else if (Uopt==USEconsistent) then

!if consistent subcycling then vn0 = UR, mean flux
!UR is defined in prim_advance_exp so that eqn for density is
!dp(n+1)=dp(n)-qsplit*dt_dynamics*DIV(UR), n is dynamics timestep
!for example, for qsplit=4 UR=((udp)_2+(udp)_4), _2 and _4 are 2nd and
!4th stages in 1 dynamics timestep

!also, Vstar is velocity for physics defined so that consistency holds.
!it turns out that Vstar at each physics stage should be then UR/(dp_tracers)
!therefore, dp_tracers(stage+1)=dp(stage)-dt_physics * DIV(UR)

!note that UR comes unprojected ie it requires a DSS, so it is done
!during the first stage in physics, because dp_tracers(stage=1)=dp(n)
!and rhs_multiplier=0

       Vstar(:,:,1,:) = elem(ie)%derived%vn0(:,:,1,:)/dp(:,:,:)
       Vstar(:,:,2,:) = elem(ie)%derived%vn0(:,:,2,:)/dp(:,:,:)
     else
       stop 'ERROR:  bad Uopt'
     endif


     ! advance Qdp
     do q=1,qsize


        do k=1,nlev
           ! div( U dp Q), 
           gradQ(:,:,1)=Vstar(:,:,1,k)*elem(ie)%state%Qdp(:,:,k,q,n0)
           gradQ(:,:,2)=Vstar(:,:,2,k)*elem(ie)%state%Qdp(:,:,k,q,n0)
           divdp = divergence_sphere(gradQ,deriv,elem(ie))
           Qtens(:,:,k)=elem(ie)%state%Qdp(:,:,k,q,n0) - dt*divdp(:,:)
           
#ifdef ENERGY_DIAGNOSTICS
           if (compute_diagnostics .and. q==1) then
              ! IEvert1_wet():  (Cpv-Cp) T Qdp_vadv  (Q equation)
              ! IEhorz1_wet():  (Cpv-Cp) T Qdp_hadv  (Q equation)
              elem(ie)%accum%IEhorz1_wet(:,:) = elem(ie)%accum%IEhorz1_wet(:,:) +&
                   (Cpwater_vapor-Cp)*elem(ie)%state%T(:,:,k,n0)*divdp(:,:)
           endif
#endif
        enddo

        if(limiter_option == 8)then
           ! add in hyperviscosity computed above:
           Qtens(:,:,:) = Qtens(:,:,:) + Qtens_biharmonic(:,:,:,q,ie)
           ! UN-DSS'ed dp at timelevel n0+1:  
           dp_star(:,:,:) = dp(:,:,:) - dt*elem(ie)%derived%divdp(:,:,:)  
           ! apply limiter to Q = Qtens / dp_star 
           ! todo: replace this with optimal limiter
           !       check that qmin/qmax are correct
           !call limiter2d_minmax(Qtens,dp_star,hvcoord,elem(ie)%spheremv,&
           !     qmin(:,q,ie),qmax(:,q,ie))

	   call limiter_optim_iter_full(Qtens(:,:,:),elem(ie)%spheremv(:,:),&
	      qmin(:,q,ie),qmax(:,q,ie),.true.,.true.,&
	      notreliable,dp_star(:,:,:),.false.)

	endif

        if(limiter_option == 84)then
           ! add in hyperviscosity computed above:
           Qtens(:,:,:) = Qtens(:,:,:) + Qtens_biharmonic(:,:,:,q,ie)
           ! UN-DSS'ed dp at timelevel n0+1:  
           dp_star(:,:,:) = dp(:,:,:) - dt*elem(ie)%derived%divdp(:,:,:)  
           ! apply limiter to Q = Qtens / dp_star 
           ! todo: replace this with optimal limiter
           !       check that qmin/qmax are correct
	   qmin(:,q,ie)=0.0d0
	   call limiter_optim_iter_full(Qtens(:,:,:),elem(ie)%spheremv(:,:),&
	      qmin(:,q,ie),qmax(:,q,ie),.true.,.false.,&
	      notreliable,dp_star(:,:,:),.false.)

	endif
        ! apply mass matrix, overwrite np1 with solution:
        ! dont do this earlier, since we allow np1 to be the same as n0
        ! and we dont want to overwrite n0 until we are done using it

        do k=1,nlev
           elem(ie)%state%Qdp(:,:,k,q,np1) = elem(ie)%spheremv(:,:)*Qtens(:,:,k) 
        enddo

        if(limiter_option == 4)then
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	  ! sign-preserving limiter, applied after mass matrix
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	  call limiter2d_zero(elem(ie)%state%Qdp(:,:,:,q,np1),&
             elem(ie)%state%ps_v(:,:,np1),hvcoord) ! ps_v argument not used
        endif


     end do

     if(DSSopt==DSSno_var)then
	call edgeVpack(edgeAdv,elem(ie)%state%Qdp(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
     else
	call edgeVpack(edgeAdv_p1,elem(ie)%state%Qdp(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
	! also DSS extra field
	do k=1,nlev
	    DSSvar(:,:,k) = elem(ie)%spheremv(:,:)*DSSvar(:,:,k) 
	enddo
	call edgeVpack(edgeAdv_p1,DSSvar(:,:,1:nlev),nlev,nlev*qsize,elem(ie)%desc)
     endif

  end do

  if(DSSopt==DSSno_var)then
     call bndry_exchangeV(hybrid,edgeAdv)
  else
     call bndry_exchangeV(hybrid,edgeAdv_p1)
  endif

  do ie=nets,nete

     if ( DSSopt == DSSeta) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
     if ( DSSopt == DSSomega) DSSvar => elem(ie)%derived%omega_p(:,:,:)
     if ( DSSopt == DSSdiv_vdp_ave) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

     if(DSSopt==DSSno_var)then
	call edgeVunpack(edgeAdv,elem(ie)%state%Qdp(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
	do q=1,qsize
	    do k=1,nlev
	      elem(ie)%state%Qdp(:,:,k,q,np1) = elem(ie)%rspheremv(:,:)*elem(ie)%state%Qdp(:,:,k,q,np1)
	    enddo
	end do
     else
	call edgeVunpack(edgeAdv_p1,elem(ie)%state%Qdp(:,:,:,:,np1),nlev*qsize,0,elem(ie)%desc)
	do q=1,qsize
	    do k=1,nlev
	      elem(ie)%state%Qdp(:,:,k,q,np1) = elem(ie)%rspheremv(:,:)*elem(ie)%state%Qdp(:,:,k,q,np1)
	    enddo
	end do
	call edgeVunpack(edgeAdv_p1,DSSvar(:,:,1:nlev),nlev,qsize*nlev,elem(ie)%desc)

	do k=1,nlev
	  DSSvar(:,:,k)=DSSvar(:,:,k)*elem(ie)%rspheremv(:,:)
	enddo

     endif
  end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
  
  end subroutine euler_step

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter_optim_iter_full(ptens,sphweights,minp,maxp,checkmin,checkmax, notreliable,dpmass,output)
!The idea here is the following: We need to find a grid field which is closest
!to the initial field (in terms of weighted sum), but satisfies the constraints.
!So, first we find values which do not satisfy constraints and bring these values
!to a closest constraint. This way we introduce some mass change (addmass),
!so, we redistribute addmass in the way that l2 error is smallest. 
!This redistribution might 
!violate constraints (though I think the solution is given by one iteration only if the 
!problem is well-posed) due to round off, for example; thus, we do a few iterations. 

!!!!!3d lim8 uses nv x nv 
!!! 2d uses np x np ????

    use kinds, only : real_kind
    use dimensions_mod, only : nv, np, nlev
    use control_mod, only : tol_limiter

    logical, intent(in) :: checkmin,checkmax
    real (kind=real_kind), dimension(nlev), intent(in)   :: minp, maxp
    real (kind=real_kind), dimension(nv,nv,nlev), intent(inout)   :: ptens
    real (kind=real_kind), dimension(nv,nv,nlev), intent(in), optional   :: dpmass
    real (kind=real_kind), dimension(nv,nv), intent(in)   :: sphweights
    real (kind=real_kind),  intent(out) :: notreliable
 
    real (kind=real_kind), dimension(nv,nv,nlev) :: weights
    real (kind=real_kind), dimension(nv,nv) :: ptens_mass
    integer  k1, k, i, j, iter
    logical :: whois_neg(np*np), whois_pos(nv*nv)

    logical, intent(in ) :: output

    real (kind=real_kind) :: addmass, weightssum, mass
    real (kind=real_kind) :: x(nv*nv),c(nv*nv)
    integer, parameter :: maxiter = 10


    do k=1,nlev
      weights(:,:,k)=sphweights(:,:)
    enddo
    notreliable=0

!if dpmass is given, we first get (\rho Q)/(\rho) fields 
    if(present(dpmass))then   
	do k=1,nlev
	  weights(:,:,k)=weights(:,:,k)*dpmass(:,:,k)
	  ptens(:,:,k)=ptens(:,:,k)/dpmass(:,:,k)
	enddo
    endif

    do k=1,nlev

      k1=1
      do i=1,nv
	do j=1,nv
	  c(k1)=weights(i,j,k)
	  x(k1)=ptens(i,j,k)
	  k1=k1+1
	enddo
      enddo

      mass=sum(c*x)


!first, check if the problem has a solution
!if it seems that it does not, we initialize a flag, 
!but the limiter still might do a good job
      if((mass-minp(k)*sum(c)<-tol_limiter).and.checkmin)then
	  notreliable=1
      endif
      if((mass-maxp(k)*sum(c)>tol_limiter).and.checkmax)then
	  notreliable=1
      endif

      iter=0
      do 
	iter=iter+1

	if(iter>maxiter)then
	  notreliable=1.0
	  exit
	endif

	addmass=0.0d0

	do k1=1,nv*nv
	  whois_neg(k1)=.true.
	  whois_pos(k1)=.true.
	  if((x(k1)>=maxp(k)).AND.(checkmax))then
	    addmass=addmass+(x(k1)-maxp(k))*c(k1)
	    x(k1)=maxp(k)
	    whois_pos(k1)=.false.
	  endif
	  if((x(k1)<=minp(k)).AND.(checkmin))then
	    addmass=addmass-(minp(k)-x(k1))*c(k1)
	    x(k1)=minp(k)
	    whois_neg(k1)=.false.
	  endif
	enddo

	if(abs(addmass)>tol_limiter*abs(mass))then

	  weightssum=0.0d0
	  if(addmass>0)then
	    do k1=1,nv*nv
	      if(whois_pos(k1))then
		weightssum=weightssum+c(k1)
	      endif
	    enddo
	    if(weightssum>0.0)then
	      do k1=1,nv*nv
		if(whois_pos(k1))then
		  x(k1)=x(k1)+addmass/weightssum
		endif
	      enddo
	    else
	      x=x+addmass/sum(c)
	      notreliable=1
	      exit
	    endif
	  else
	    do k1=1,nv*nv
	      if(whois_neg(k1))then
		weightssum=weightssum+c(k1)
	      endif
	    enddo
	    if(weightssum>0.0)then
	      do k1=1,nv*nv
		if(whois_neg(k1))then
		  x(k1)=x(k1)+addmass/weightssum
		endif
	      enddo
	    else
	      x=x+addmass/sum(c)
	      notreliable=1
	      exit
	    endif
	  endif

	else
	  x=x+addmass/sum(c)
	  exit
	endif

      enddo!end of iteration

      k1=1
      do i=1,nv
	do j=1,nv
	  ptens(i,j,k)=x(k1)
	  k1=k1+1
	enddo
      enddo

    enddo

    if(present(dpmass))then   
	do k=1,nlev
	  ptens(:,:,k)=ptens(:,:,k)*dpmass(:,:,k)
	enddo
    endif


  end subroutine limiter_optim_iter_full

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter2d_minmax(Q,dp,hvcoord,spheremv,qmin,qmax)
!
! mass conserving limiter (2D only).  to be called just before DSS
!
! in pure 2D advection, the element mass will not be negative before DSS
! this routine will redistribute to remove negative values (conservative)
!
! if used in 3D, should be applied with 2D/vertical split advection
! 
! call with Qdp and assocated dp
!
!
  implicit none
  real (kind=real_kind), intent(inout)  :: Q(nv,nv,nlev)
  real (kind=real_kind), intent(in)  :: spheremv(nv,nv)
  real (kind=real_kind), intent(in) ::  dp(nv,nv,nlev)
  type (hvcoord_t)                 :: hvcoord

  ! local
  integer i,j,k
  real (kind=real_kind) :: mass,mass_new,area,qmin(nlev),qmax(nlev),mass2


  do k=1,nlev
     mass=sum( Q(:,:,k)*spheremv(:,:) )
     area=sum( dp(:,:,k)*spheremv(:,:) )

     Q(:,:,k)=Q(:,:,k)/dp(:,:,k)  ! % convert to concentration


!     if (mass>0) print *,k,mass/area,qmin(k),qmax(k)
     ! max limiter
     if ( maxval(Q(:,:,k)) > qmax(k) ) then
        
        Q(:,:,k)=qmax(k)-Q(:,:,k)      ! some of these will be negative
        mass2 = area*qmax(k) - mass

        if (mass2 < 0) Q(:,:,k)=-Q(:,:,k) 
        mass_new=0
        do j=1,nv	
        do i=1,nv
           if (Q(i,j,k)<0) then
              Q(i,j,k)=0
           else
              mass_new = mass_new + Q(i,j,k)*dp(i,j,k)*spheremv(i,j)
           endif
        enddo
        enddo
     
        ! now scale the all positive values to restore mass
        if (mass_new>0) Q(:,:,k) = Q(:,:,k)*abs(mass2)/mass_new
        if (mass2 < 0) Q(:,:,k)=-Q(:,:,k) 
        
        Q(:,:,k)=qmax(k)-Q(:,:,k)
     endif


     ! min limiter
     if ( minval(Q(:,:,k)) < qmin(k) ) then
        Q(:,:,k)=Q(:,:,k)-qmin(k)
        mass2 = mass - area*qmin(k)
        ! negative mass.  so reduce all postive values to zero 
        ! then increase negative values as much as possible
        if (mass2 < 0) Q(:,:,k)=-Q(:,:,k) 
        mass_new=0
        do j=1,nv	
           do i=1,nv
              if (Q(i,j,k)<0) then
                 Q(i,j,k)=0
              else
                 mass_new = mass_new + Q(i,j,k)*dp(i,j,k)*spheremv(i,j)
              endif
           enddo
        enddo
        
        ! now scale the all positive values to restore mass
        if (mass_new>0) Q(:,:,k) = Q(:,:,k)*abs(mass2)/mass_new
        if (mass2 < 0) Q(:,:,k)=-Q(:,:,k) 

        Q(:,:,k)=Q(:,:,k)+qmin(k)
     endif

     Q(:,:,k)=Q(:,:,k)*dp(:,:,k)

  enddo

end subroutine 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter2d_zero(Q,ps,hvcoord)
!
! mass conserving zero limiter (2D only).  to be called just before DSS
!
! this routine is called inside a DSS loop, and so Q had already
! been multiplied by the mass matrix.  Thus dont include the mass
! matrix when computing the mass = integral of Q over the element
!
! ps is only used when advecting Q instead of Qdp
! so ps should be at one timelevel behind Q
!
  implicit none
  real (kind=real_kind), intent(inout)  :: Q(nv,nv,nlev)
  real (kind=real_kind), intent(in)  :: ps(nv,nv)
  type (hvcoord_t)                 :: hvcoord

  ! local
  real (kind=real_kind) ::  dp(nv,nv)
  integer i,j,k
  real (kind=real_kind) :: mass,mass_new,ml


  do k=nlev,1,-1

     mass=0
     do j=1,nv	
     do i=1,nv
        !ml = Q(i,j,k)*dp(i,j)*spheremv(i,j)  ! see above
        if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
           ml = Q(i,j,k)
        else
           dp(i,j) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*ps(i,j)
           ml = Q(i,j,k)*dp(i,j)
        endif
        mass = mass + ml
     enddo
     enddo

     ! negative mass.  so reduce all postive values to zero 
     ! then increase negative values as much as possible
     if (mass < 0) Q(:,:,k)=-Q(:,:,k) 
     mass_new=0
     do j=1,nv	
     do i=1,nv
        if (Q(i,j,k)<0) then
           Q(i,j,k)=0
        else
           if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
              ml = Q(i,j,k)
           else
              ml = Q(i,j,k)*dp(i,j)
           endif
           mass_new = mass_new + ml
        endif
     enddo
     enddo

     ! now scale the all positive values to restore mass
     if ( mass_new > 0) Q(:,:,k) = Q(:,:,k)*abs(mass)/mass_new
     if (mass < 0) Q(:,:,k)=-Q(:,:,k) 
  enddo

end subroutine 

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter3d_noncon(Q,ps,hvcoord,mass_added)
!
! local element negative mass "fixer"  
! Trucate negative values, then use CAM style mass fixer to restore
! original mass.  Mass fixer is local to the element.  Will break C0
! do DSS needs to be done afterwards.
! 
! this routine is called inside a DSS loop, and so Q had already
! been multiplied by the mass matrix.  Thus dont include the mass
! matrix when computing the mass = integral of Q over the element
!
! ps is only used when advecting Q instead of Qdp
! so ps should be at one timelevel behind Q
!
  implicit none
  real (kind=real_kind), intent(inout)  :: Q(nv,nv,nlev)
  real (kind=real_kind), intent(in)  :: ps(nv,nv)
  real (kind=real_kind), intent(out)  :: mass_added
  type (hvcoord_t)                 :: hvcoord

  ! local
  real (kind=real_kind) ::  dp(nv,nv)
  integer i,j,k
  real (kind=real_kind) :: mass,mass_new,ml

  mass_added=0
  mass=0
  mass_new=0


  if (tracer_advection_formulation==TRACERADV_TOTAL_DIVERGENCE) then
  do k=nlev,1,-1
     do j=1,nv	
        do i=1,nv
!           ml = Q(i,j,k)*dp(i,j)*spheremv(i,j)  ! see above
           mass = mass + Q(i,j,k)
           if (Q(i,j,k)<0) then
              Q(i,j,k)=0
           else
              mass_new = mass_new + Q(i,j,k)
           endif
        enddo
     enddo
  enddo
  else
  do k=nlev,1,-1
     dp(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*ps(:,:)
     do j=1,nv	
        do i=1,nv
!           ml = Q(i,j,k)*dp(i,j)*spheremv(i,j)  ! see above
           ml = Q(i,j,k)*dp(i,j)
           mass = mass + ml
           
           if (Q(i,j,k)<0) then
              Q(i,j,k)=0
           else
              mass_new = mass_new + ml
           endif
        enddo
     enddo
  enddo
  endif
 
  ! now scale the all positive values to restore mass
  if ( mass > 0) then
     ! rescale positive values to restore original mass
     Q(:,:,:) = Q(:,:,:)*mass/mass_new
  else
     ! mass was negative.  set all values to zero
     Q(:,:,:) = 0
     mass_added = mass_added -mass
  endif
  
  
  end subroutine 
  
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine preq_impsysQ(elemin,hvcoord,np1,n0,nm1,Q,Q_vadv)
    use perf_mod, only : t_startf, t_stopf                    ! _EXTERNAL
    implicit none

    type (element_t), intent(in)     :: elemin
    type (hvcoord_t)                 :: hvcoord
    integer, intent(in)              :: np1
    integer, intent(in)              :: n0
    integer, intent(in)              :: nm1

    real(kind=real_kind), intent(in),target, dimension(nv,nv,nlev)   :: Q         
    real(kind=real_kind), intent(inout), dimension(nv,nv,nlev)       :: Q_vadv    

    ! ==============================
    ! Local Variables
    ! ==============================
    
    real(kind=real_kind), dimension(nlevp)      :: hyai
    real(kind=real_kind), dimension(nlevp)      :: hybi
    real(kind=real_kind), dimension(nlev)       :: hyam
    real(kind=real_kind), dimension(nlev)       :: hybm
    real(kind=real_kind), dimension(nlev)       :: hybd 
    real(kind=real_kind), dimension(nv,nv,nlev) :: div

    real(kind=real_kind), dimension(nv,nv)      :: ps                ! surface pressure
    real(kind=real_kind), dimension(nv,nv)      :: rps               ! 1/ps
    real(kind=real_kind), dimension(nv,nv,nlev) :: rpdel             ! 1./pdel

    real(kind=real_kind), dimension(nv,nv,nlevp)  :: eta_dot_dp_deta       ! eta dot * dp/deta at time level n
    real(kind=real_kind), dimension(nv,nv,nlev)   :: vgrad_ps              ! ps*(v.grad(lnps))

    real(kind=real_kind) :: pint(nv,nv,nlevp)
    real(kind=real_kind) :: pdel(nv,nv,nlev)
    real(kind=real_kind) :: pmid(nv,nv,nlev)

    real(kind=real_kind) :: v1, v2, vcon1, vcon2

    integer :: i,j,k,l

    call t_startf('preq_impsysQ')




    hyai   = hvcoord%hyai
    hybi   = hvcoord%hybi
    hyam   = hvcoord%hyam
    hybm   = hvcoord%hybm
    hybd   = hvcoord%hybd

    ps(:,:)  = EXP(elemin%state%lnps(:,:,n0))
    rps(:,:) = 1.0_real_kind/ps(:,:)

    call preq_pressure(hvcoord%ps0,  ps,&
         hyai, hybi, hyam, hybm,&
         pint, pmid, pdel)

    rpdel = 1.0_real_kind/pdel

!   v should be contravariant for explicit time step

    do k=1,nlev
       do j=1,nv
          do i=1,nv
             v1 = elemin%state%v(i,j,1,k,n0)
             v2 = elemin%state%v(i,j,2,k,n0)

!            vcon1 = elemin%metinv(1,1,i,j)*v1 + elemin%metinv(1,2,i,j)*v2
!            vcon2 = elemin%metinv(2,1,i,j)*v1 + elemin%metinv(2,2,i,j)*v2

             vcon1 = elemin%Dinv(1,1,i,j)*v1 + elemin%Dinv(1,2,i,j)*v2
             vcon2 = elemin%Dinv(2,1,i,j)*v1 + elemin%Dinv(2,2,i,j)*v2

             vgrad_ps(i,j,k) = ps(i,j)*(vcon1*elemin%derived%grad_lnps(i,j,1) + vcon2*elemin%derived%grad_lnps(i,j,2))
          end do
       end do
    end do

    eta_dot_dp_deta(:,:,1) = 0.0_real_kind

    do k=1,nlev
       div(:,:,k) = elemin%derived%div(:,:,k,n0)
       do j=1,nv
          do i=1,nv
             eta_dot_dp_deta(i,j,k+1) = eta_dot_dp_deta(i,j,k) + vgrad_ps(i,j,k)*hybd(k) + div(i,j,k)*pdel(i,j,k)
          end do
       end do
    end do

    do k=1,nlev-1
       do j=1,nv
          do i=1,nv
             eta_dot_dp_deta(i,j,k+1) = hybi(k+1)*eta_dot_dp_deta(i,j,nlev+1) - eta_dot_dp_deta(i,j,k+1)
          end do
       end do
    end do

    eta_dot_dp_deta(:,:,nlev+1) = 0.0_real_kind

    call preq_vertadvQ(Q(:,:,:),eta_dot_dp_deta,rpdel,Q_vadv)
    call t_stopf('preq_impsysQ')

  end subroutine preq_impsysQ



  subroutine preq_vertadvQ(Q,eta_dot_dp_deta, rpdel, Q_vadv)
    use perf_mod, only : t_startf, t_stopf                   ! _EXTERNAL
    implicit none

    real (kind=real_kind), intent(in)  :: Q(nv,nv,nlev)
    real (kind=real_kind), intent(in)  :: eta_dot_dp_deta(nv,nv,nlevp)
    real (kind=real_kind), intent(in)  :: rpdel(nv,nv,nlev)

    real (kind=real_kind), intent(inout) :: Q_vadv(nv,nv,nlev)

    ! ========================
    ! Local Variables
    ! ========================

    integer                            :: i,j,k,l
    real (kind=real_kind)              :: facp, facm

    ! ===========================================================
    ! Compute vertical advection of Q from eq. (3.b.1)
    !
    ! k = 1 case:
    ! ===========================================================
    call t_startf('preq_vertadvQ')

    k=1
    do j=1,nv
       do i=1,nv 
          facp  = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)             
          Q_vadv(i,j,k) = facp*(Q(i,j,k+1)- Q(i,j,k))
       end do
    end do

    ! ===========================================================
    ! vertical advection
    !
    ! 1 < k < nlev case:
    ! ===========================================================

    do k=2,nlev-1
       do j=1,nv
          do i=1,nv
             facp = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k+1)
             facm = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)                
             Q_vadv(i,j,k) = facp*(Q(i,j,k+1)- Q(i,j,k)) + facm*(Q(i,j,k)- Q(i,j,k-1))
          enddo
       end do
    end do

    ! ===========================================================
    ! vertical advection
    !
    ! k = nlev case:
    ! ===========================================================

    k=nlev
    do j=1,nv
       do i=1,nv
          facm = (0.5_real_kind*rpdel(i,j,k))*eta_dot_dp_deta(i,j,k)             
          Q_vadv(i,j,k) = facm*(Q(i,j,k)- Q(i,j,k-1))
       enddo
    end do
    call t_stopf('preq_vertadvQ')

  end subroutine preq_vertadvQ




! compute the vertical advection term for this equation:
!    d/dt[dp Q] +  div( U dp Q ) + d( eta_dot_dpdn Q ) = 0
!
! qsplit=1 case:  
!   dp(t+1) = dp(t-1) + 2dt d[ eta_dot_dpdn(t) ]
!  
!  z1 = grid with intervales of size dp(t-1)    grid at time t-1
!  Lagrangian advection (mass preserving):  z1 grid unchanged, Q unchanged
!  z2 = grid with intervales of size dp(t+1)    grid at time t+1
!
!  map Q(t-1) on z1 grid to Q(t+1) on z2 grid
! 
!  Qtend = [ Q(t+1) dp(t+1) - Q(t-1) dp(t-1) ]  / 2 dt
!
!  Note: Q=1 is preserved by remap code.  Thus taking Q=1 will result in
!  [ dp(t+1) - dp(t-1) ] / 2 dt  = d[eta_dot_dpdn(t) ] which is the same 
!  term in the Primative Equations mass continutity equation (implicit
!  in the definition of eta_dot_dpdn)
!  (wont be exact if qsplit>1)
!
  subroutine preq_vertadv_remap(Q,dp1,dp2,dt2,hvcoord,Q_vadv,use_limiter)
    use remap_lauritzen, only :  verremap2                 ! _EXTERNAL 
    use perf_mod, only : t_startf, t_stopf                 ! _EXTERNAL

    implicit none
    real (kind=real_kind), intent(inout)  :: Q(nv,nv,nlev)
    real (kind=real_kind), intent(out)    :: Q_vadv(nv,nv,nlev)
    real (kind=real_kind), intent(in)     :: dp1(nv,nv,nlev),dp2(nv,nv,nlev)
    type (hvcoord_t)                      :: hvcoord
    logical, intent(in)                   :: use_limiter

    ! ========================
    ! Local Variables
    ! ========================
    integer :: check_global=1
    real (kind=real_kind)  :: qmin=0
    real (kind=real_kind)  :: qmax=1d50

    integer                            :: i,j,k,l
    real (kind=real_kind)  :: dt2
    real (kind=real_kind)  :: Qold(nlev,nv,nv),Qnew(nlev,nv,nv)
    real (kind=real_kind)  :: z2c(nlevp,nv,nv),z1c(nlevp,nv,nv)
    integer :: ipr(nlev),ik,ippm


    call t_startf('preq_vertadv_remap')
    if (use_limiter) then
       check_global=1
!       ippm=1   ! about 2x slower
       ippm=3
    else
       check_global=0
       ippm=2
    endif

    z1c(1,:,:)=0
    z2c(1,:,:)=0
    do k=1,nlev
       Qold(k,:,:)=Q(:,:,k)
       z1c(k+1,:,:) = z1c(k,:,:)+dp1(:,:,k)
       z2c(k+1,:,:) = z2c(k,:,:)+dp2(:,:,k)
    enddo
    do i=1,nv
    do j=1,nv

       ! 0 piecewise cubic
       ! 1 piecewise cubic with UK met office monotonicity constraints
       ! 2 quadratic splies
       ! 3 quadratic splies with UK met office monotonicity constraints
       call VERREMAP2(z1c(:,i,j),z2c(:,i,j),Qold(:,i,j),nlev,Qnew(:,i,j),&
           ippm,qmin,qmax,check_global) 
    enddo
    enddo
    do k=1,nlev
       ! for diagnostics:  Qdp(t-1) - Qdp(t+1) / dt2
       Q_vadv(:,:,k) = (  Q(:,:,k) - Qnew(k,:,:)   ) / dt2
       Q(:,:,k) = Qnew(k,:,:)
    enddo
    call t_stopf('preq_vertadv_remap')
end subroutine preq_vertadv_remap
       


  subroutine advance_hypervis_scalar_lf(edgeAdv,elem,hvcoord,hybrid,deriv,nt,n0,nets,nete,dt2)
  !
  !  hyperviscosity operator used by leapfrog scheme
  !  take one timestep of:  
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use kinds, only : real_kind
  use dimensions_mod, only : nv, np, nlev
  use control_mod, only : hypervis_order, hypervis_subcycle_q
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use perf_mod, only: t_startf, t_stopf                  ! _EXTERNAL
  implicit none

  type (hvcoord_t), intent(in)      :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edgeAdv
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: dt2
  integer :: nets,nete,nt,n0

  
  ! local
  integer :: k,kptr,i,j,ie,ic,q
  real (kind=real_kind), dimension(nv,nv,nlev,qsize,nets:nete) :: Qtens

! NOTE: PGI compiler bug: when using spheremv, rspheremv and ps as pointers to elem(ie)% members,
!       data is incorrect (offset by a few numbers actually)
!       removed for now.  
!  real (kind=real_kind), dimension(:,:), pointer :: spheremv,rspheremv

  real (kind=real_kind), dimension(nv,nv) :: lap_p
  real (kind=real_kind) :: v1,v2,dt,nu_scale,dp,dp0
  integer :: density_scaling = 0
  
  if (nu_q == 0) return;
  call t_startf('advance_hypervis_scalar_lf')

  if (tracer_advection_formulation == TRACERADV_UGRADQ ) then
     ! conservative advection using non-conservative form of the equations
     ! in this case, we use a density scaled viscosity coefficieint: 
     density_scaling = 1
  endif
  
  
  dt=dt2/hypervis_subcycle_q
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     do ic=1,hypervis_subcycle_q
        do ie=nets,nete
           !spheremv     => elem(ie)%spheremv
           do q=1,qsize
              do k=1,nlev
                 lap_p=laplace_sphere_wk(elem(ie)%state%Q(:,:,k,q,nt),deriv,elem(ie))
                 ! advace in time.  (note: DSS commutes with time stepping, so we
                 ! can time advance and then DSS.
                 do j=1,nv
                 do i=1,nv             
                    elem(ie)%state%Q(i,j,k,q,nt)=elem(ie)%state%Q(i,j,k,q,nt)*elem(ie)%spheremv(i,j)  +  dt*nu_q*lap_p(i,j) 
                 enddo
                 enddo
              enddo
           enddo
           call edgeVpack(edgeAdv, elem(ie)%state%Q(:,:,:,:,nt),nlev*qsize,0,elem(ie)%desc)
        enddo
           
        call bndry_exchangeV(hybrid,edgeAdv)
        
        do ie=nets,nete
           call edgeVunpack(edgeAdv, elem(ie)%state%Q(:,:,:,:,nt),nlev*qsize,0,elem(ie)%desc)
           do q=1,qsize
              !rspheremv     => elem(ie)%rspheremv
              ! apply inverse mass matrix
              do k=1,nlev
              do j=1,nv
              do i=1,nv             
                 elem(ie)%state%Q(i,j,k,q,nt)=elem(ie)%rspheremv(i,j)*elem(ie)%state%Q(i,j,k,q,nt)
              enddo
              enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
     enddo  ! subcycle
  endif
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle_q
        do ie=nets,nete
           if (density_scaling==1) then
              ! state%Q really is Q !
              Qtens(:,:,:,:,ie)=elem(ie)%state%Q(:,:,:,:,nt)
           else
              ! state%Q is really Qdp.  but we only apply diffusion on Q
              do k=1,nlev
              do j=1,nv
              do i=1,nv
                 ! note: use ps(t+1) to get exact consistency
                 dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
                 Qtens(i,j,k,:,ie)=elem(ie)%state%Q(i,j,k,:,nt)/dp
              enddo
              enddo
              enddo
           endif
        enddo
        ! compute biharmonic operator. Qtens = input and output 
        call biharmonic_wk_scalar(elem,Qtens,deriv,edgeAdv,hybrid,nt,nets,nete)
        do ie=nets,nete
           !spheremv     => elem(ie)%spheremv
           do q=1,qsize
           do k=1,nlev
           do j=1,nv
           do i=1,nv

              dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                   ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0

              if (density_scaling==1) then
                 ! advection Q.  For conservation:     dp0/dp DIFF(Q)
                 ! scale velosity by 1/rho (normalized to be O(1))
                 ! dp = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                 ! NOTE: use ps(t) to get exact mass conservation
                 dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,n0)
                 nu_scale = dp0/dp
                 elem(ie)%state%Q(i,j,k,q,nt)  =  elem(ie)%state%Q(i,j,k,q,nt)*elem(ie)%spheremv(i,j) &
                              -dt*nu_q*nu_scale*Qtens(i,j,k,q,ie)
              else
                 ! advection Qdp.  For mass advection consistency:
                 ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )  
!                 elem(ie)%state%Q(i,j,k,q,nt)  =  elem(ie)%state%Q(i,j,k,q,nt)*elem(ie)%spheremv(i,j) &
!                        -dt*nu_q*Qtens(i,j,k,q,ie)
                 elem(ie)%state%Q(i,j,k,q,nt)  =  elem(ie)%state%Q(i,j,k,q,nt)*elem(ie)%spheremv(i,j) &
                        -dt*nu_q*dp0*Qtens(i,j,k,q,ie)
              endif
           enddo
           enddo
           enddo

           ! smooth some of the negativities introduced by diffusion:
           ! note: ps_v not used if advecting Qdp
           call limiter2d_zero(elem(ie)%state%Q(:,:,:,q,nt),&
                elem(ie)%state%ps_v(:,:,n0),hvcoord)

           enddo
           call edgeVpack(edgeAdv,elem(ie)%state%Q(:,:,:,:,nt),qsize*nlev,0,elem(ie)%desc)
        enddo

        call bndry_exchangeV(hybrid,edgeAdv)

        do ie=nets,nete
        call edgeVunpack(edgeAdv, elem(ie)%state%Q(:,:,:,:,nt), qsize*nlev, 0, elem(ie)%desc)
        !rspheremv     => elem(ie)%rspheremv
        do q=1,qsize    
           ! apply inverse mass matrix
           do k=1,nlev
           do j=1,nv
           do i=1,nv
              elem(ie)%state%Q(i,j,k,q,nt)=elem(ie)%rspheremv(i,j)*elem(ie)%state%Q(i,j,k,q,nt)
           enddo
           enddo
           enddo
        enddo
        enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
     enddo
  endif

  call t_stopf('advance_hypervis_scalar_lf')
  
  end subroutine advance_hypervis_scalar_lf



  subroutine advance_hypervis_scalar(edgeAdv,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2)
  !
  !  hyperviscsoity operator for foward-in-time scheme
  !  take one timestep of:  
  !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use kinds, only : real_kind
  use dimensions_mod, only : nv, np, nlev
  use control_mod, only : hypervis_order, hypervis_subcycle_q
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  use perf_mod, only: t_startf, t_stopf                          ! _EXTERNAL
  implicit none

  type (hvcoord_t), intent(in)      :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edgeAdv
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: dt2
  integer :: nets,nete,nt

  
  ! local
  integer :: k,kptr,i,j,ie,ic,q
  real (kind=real_kind), dimension(nv,nv,nlev,qsize,nets:nete) :: Qtens
  real (kind=real_kind), dimension(nv,nv,nlev) :: dp

  real (kind=real_kind), dimension(nlev,qsize,nets:nete) :: min_neigh
  real (kind=real_kind), dimension(nlev,qsize,nets:nete) :: max_neigh

! NOTE: PGI compiler bug: when using spheremv, rspheremv and ps as pointers to elem(ie)% members,
!       data is incorrect (offset by a few numbers actually)
!       removed for now.  
!  real (kind=real_kind), dimension(:,:), pointer :: spheremv,rspheremv

  real (kind=real_kind), dimension(nv,nv) :: lap_p
  real (kind=real_kind) :: v1,v2,dt,dp0
  integer :: density_scaling = 0
  
  if (nu_q == 0) return;
  if (hypervis_order /= 2) return
  call t_startf('advance_hypervis_scalar')
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dt=dt2/hypervis_subcycle_q


  do ic=1,hypervis_subcycle_q
     do ie=nets,nete
        ! Qtens = Q/Qdp
        do k=1,nlev
           dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,nt)
           do q=1,qsize
              Qtens(:,:,k,q,ie)=elem(ie)%state%Qdp(:,:,k,q,nt)/dp(:,:,k)
           enddo
        enddo
     enddo
     ! compute biharmonic operator. Qtens = input and output 
     call biharmonic_wk_scalar(elem,Qtens,deriv,edgeAdv,hybrid,nt,nets,nete)
     do ie=nets,nete
        !spheremv     => elem(ie)%spheremv
        do q=1,qsize
           do k=1,nlev
           do j=1,nv
           do i=1,nv
              dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                   ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
              ! advection Qdp.  For mass advection consistency:
              ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )  
              elem(ie)%state%Qdp(i,j,k,q,nt)  =  elem(ie)%state%Qdp(i,j,k,q,nt)*elem(ie)%spheremv(i,j) &
                   -dt*nu_q*dp0*Qtens(i,j,k,q,ie)
           enddo
           enddo
           enddo

           ! smooth some of the negativities introduced by diffusion:
           call limiter2d_zero(elem(ie)%state%Qdp(:,:,:,q,nt),&
                elem(ie)%state%ps_v(:,:,nt),hvcoord)
        enddo
        call edgeVpack(edgeAdv,elem(ie)%state%Qdp(:,:,:,:,nt),qsize*nlev,0,elem(ie)%desc)
     enddo

     call bndry_exchangeV(hybrid,edgeAdv)
     
     do ie=nets,nete
        call edgeVunpack(edgeAdv, elem(ie)%state%Qdp(:,:,:,:,nt), qsize*nlev, 0, elem(ie)%desc)
        !rspheremv     => elem(ie)%rspheremv
        do q=1,qsize    
           ! apply inverse mass matrix
           do k=1,nlev
              elem(ie)%state%Qdp(:,:,k,q,nt)=elem(ie)%rspheremv(:,:)*elem(ie)%state%Qdp(:,:,k,q,nt)
           enddo
        enddo
     enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
  enddo
  call t_stopf('advance_hypervis_scalar')
  end subroutine advance_hypervis_scalar









  subroutine preq_vertadv_dpQ(Q,eta_dot_dp_deta, Q_vadv)
    implicit none
    
    real (kind=real_kind), intent(in)  :: Q(nv,nv,nlev)
    real (kind=real_kind), intent(in)  :: eta_dot_dp_deta(nv,nv,nlevp)
    real (kind=real_kind), intent(inout) :: Q_vadv(nv,nv,nlev)

    ! ========================
    ! Local Variables
    ! ========================

    integer                            :: i,j,k,l
    real (kind=real_kind)              :: qp,qm
    !
    ! compute  d(eta_dot_dp_deta Q)
    ! 
    !
    qp=0
    qm=0
    do k=1,nlev
       do i=1,nv
       do j=1,nv
       ! at k=1,nlev, eta_dot_dp_deta is zero, so value of q does not matter
#undef QADV_UPWIND
#ifdef QADV_UPWIND
       ! UPWIND, 1st order:
          if (k/=nlev) then
             if (eta_dot_dp_deta(i,j,k+1) > 0 ) then
                qp=Q(i,j,k)
             else
                qp=Q(i,j,k+1)
             endif
          endif
          if (k/=1) then
             if (eta_dot_dp_deta(i,j,k) > 0 ) then
                qm=Q(i,j,k-1)
             else
                qm=Q(i,j,k)
             endif
          endif
#else
          ! CENTERED, 2nd order:
          if (k/=nlev) then
             qp=.5*(Q(i,j,k)+Q(i,j,k+1))
          endif
          if (k/=1) then
             qm=.5*(Q(i,j,k)+Q(i,j,k-1))
          endif
#endif
          Q_vadv(i,j,k)= eta_dot_dp_deta(i,j,k+1)*qp - eta_dot_dp_deta(i,j,k)*qm
       enddo
      enddo
    enddo

end subroutine preq_vertadv_dpQ
end module prim_advection_mod












 
