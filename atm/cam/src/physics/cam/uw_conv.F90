
module uw_conv
  !
  ! $Id$
  !
  use cam_logfile, only: iulog
  implicit none
  private
  save
  !
  ! Public interfaces
  !
  public init_uw_conv      !  Initialization of data for moist convection
  public compute_uw_conv   !  UW shallow convection

  !
  !  Private data
  !
  integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
  real(r8) cpair          ! specific heat of dry air
  real(r8) gravit        ! gravitational constant       
  real(r8) rair        ! gas constant for dry air
  real(r8) latvap
  real(r8) latice
  real(r8) latsub
  real(r8) zvir
  real(r8) epsilo
  real(r8) cappa
  real(r8) tmelt
  real(r8) rgasv

contains

  subroutine init_uw_conv(kind, rair_in    ,cpair_in   ,gravit_in  ,latvap_in, &
       latice_in, zvir_in, epsilo_in, cappa_in, tmelt_in, rgasv_in  )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Initialize UW shallow convection scheme
    ! 
    !-----------------------------------------------------------------------
    !
    ! Input arguments
    !
    integer, intent(in)  :: kind      ! kind of reals being passed in
    real(r8), intent(in) :: rair_in              ! gas constant for dry air
    real(r8), intent(in) :: cpair_in             ! specific heat of dry air
    real(r8), intent(in) :: gravit_in            ! acceleration due to gravity
    real(r8), intent(in) :: latvap_in            ! latent heat of vaporization
    real(r8), intent(in) :: latice_in            ! latent heat of vaporization
    real(r8), intent(in) :: zvir_in
    real(r8), intent(in) :: epsilo_in
    real(r8), intent(in) :: cappa_in
    real(r8), intent(in) :: tmelt_in
    real(r8), intent(in) :: rgasv_in

    if ( kind .ne. r8 ) then
       write(iulog,*) 'KIND of reals passed to init_diffusvity -- exiting.'
       stop 'init_eddy_diff'
    endif

    !
    !-----------------------------------------------------------------------
    !
    ! Initialize physical constants for moist convective mass flux procedure
    !
    cpair     = cpair_in         ! specific heat of dry air
    gravit   = gravit_in        ! gravitational constant
    rair   = rair_in          ! gas constant for dry air
    latvap = latvap_in
    latice = latice_in
    latsub = latvap + latice
    zvir = zvir_in
    epsilo = epsilo_in
    cappa = cappa_in
    tmelt = tmelt_in
    rgasv = rgasv_in
    !
    return
  end subroutine init_uw_conv

  subroutine compute_uw_conv(pcols, pver, ncol, ncnst, &
       ztodt, pmid,  pdel, &
       zmid, pblht, tb, q,   &
       pint,  zint,  ub, vb,  tkeb, cush,      &
       sten, qvten,    &
       cmf, cmfdqr, cmfsl, cmflq,  prec_cmf,   &
       qc, toplev, botlev, cldqc, rliq,    &
       uten, vten, cldfrac, qsat)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !     SHALLOW CONVECTION SCHEME
    !     Described in McCaa, Bretherton, and Grenier:
    !     (submitted to MWR, December 2001)
    !     For info contact Jim McCaa: mccaa@u.washington.edu
    !
    !     Inputs: pressure, temperature, heights, vert. vel.,tke,
    !     specific humidity, cloud water mixing ratio, horizontal winds
    !
    !     Outputs: updated tendencies of specific humidity, temperature,
    !     cloud, ice, rain, and snow mixing ratios
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    
    use error_function, only: erfc

    implicit none

    integer, intent(in) :: pcols
    integer, intent(in) :: pver
    integer, intent(in) :: ncol                       ! number of atmospheric columns
    integer, intent(in) :: ncnst
    real(r8), intent(in) :: ztodt         ! 2 delta-t (seconds)
    real(r8), intent(in) :: pmid(pcols,pver)    ! Pressure at model mid-levels (pascals)
    real(r8), intent(in) :: pdel(pcols,pver)    ! delta-p 
    real(r8), intent(in) :: zmid(pcols,pver)    ! Height at model mid-levels (m)
    real(r8), intent(in) :: pblht(pcols) ! pbl height
    real(r8), intent(in) :: tb(pcols,pver) ! temperature profile (K)
    real(r8), intent(in) :: q(pcols,pver,ncnst) ! specific humidity (g/g))
    real(r8), intent(in) :: pint(pcols,pver+1)   ! Pressure at model interfaces (pascals)
    real(r8), intent(in) :: zint(pcols,pver+1)   ! height at model interfaces (m)
    real(r8), intent(in) :: ub(pcols,pver)
    real(r8), intent(in) :: vb(pcols,pver) ! wind profile (m/s)
    real(r8), intent(in) :: tkeb(pcols,pver+1) ! turubulence kinetic energy (m2/s2) at mid-levels
    real(r8), intent(inout) :: cush(pcols)    ! convective scale height (m) 
    real(r8), intent(out) :: sten(pcols,pver) ! tendency of dry static energy
    real(r8), intent(out) :: qvten(pcols,pver) ! tendency of specific humidity
    real(r8), intent(out) :: cmf(pcols,pver+1)! cloud mass flux at level above layer (kg/m2/s)
    real(r8), intent(out) :: cmfdqr(pcols,pver)  ! dq/dt due to moist convective rainout
    real(r8), intent(out) :: cmfsl(pcols,pver)   ! Moist convection lw stat energy flux
    real(r8), intent(out) :: cmflq(pcols,pver)   ! Moist convection total water flux
    real(r8), intent(out) :: prec_cmf(pcols) ! rain rate (kg/m2/s)
    real(r8), intent(out) :: qc(pcols,pver)      ! dq/dt due to rainout terms (same as cmfdqr)
    real(r8), intent(out) :: toplev(pcols)   ! top level of cloud
    real(r8), intent(out) :: botlev(pcols)   ! bottom level of cloud
    real(r8), intent(out) :: cldqc(pcols,pver) ! in-cloud liquid water mixing ratio
    real(r8), intent(out) :: rliq(pcols)
    real(r8), intent(out) :: uten(pcols,pver)! tendency of wind
    real(r8), intent(out) :: vten(pcols,pver)! tendency of wind
    real(r8), intent(out) :: cldfrac(pcols,pver) ! convective cloud fraction
    integer, external     :: qsat
    !
    !...DEFINE LOCAL VARIABLES...
    !
    !     VARIABLES WHICH DESCRIBE THE SOURCE LAYER
    real(r8) thlsrc,qtsrc,usrc,vsrc
    !     VARIABLES WHICH DESCRIBE THE UPDRAFT AT ITS LCL
    real(r8) zlcl,wrel
    !     VARIABLES WHICH DESCRIBE THE ENVIRONMENT AT THE UPDRAFT LCL
    real(r8)          &
         ee2,ud2,wtw, &
         cldhgt,dpsum,&
         xs,          &       ! fraction of environmental air in a neutrally buoyant mixture
         cin                  ! convective inhibition (m2/s2)
    integer           &
         K,           &        ! release level (half-sigma level just below lcl)
         KLCL,        &        ! first half-sigma level above lcl
         ktoppbl,     &        ! top of the
         i,           &
         let,         &        ! level where buoyancy becomes negative
         ltop,        &       ! highest level with positive vertical velocity
         nk,          &
         j,ibeg,iend,igauss,iteration
    integer iout,ik1,jk1,ik2,jk2,ik3,jk3 ! for diagnostic column output
    integer kinv ! inversion level
    integer km1,kthvmin
    integer krel ! level of release
    integer kp1,kl,klm

    !     And starring... (in the order of their appearance)
    real(r8) scaleh,rmfk1,rmfk2,nu, tkeavg, thvuinv,qct0bot,     &
         thj, qvj, qlj, qij, tscaleh,ssthc0b,ssqct0b,        &
         pinv, zinv, cinlcl, qj,rbuoy,rdrag,exnerlcl,exnerinv,             &
         cbmf, wexp, wcrit, sigmaw, ufrc,ssthc0a,ssqct0a,thc0bot,    &
         thvj, qrj, qsj, thcten, qctten, bigc,rle,rkm,rpen,thc0top,    &
         qct0top,thc0lcl,qct0lcl,thv0lcl,rho0lcl,thvubot,thvutop,    &
         rkfre,prel,thv0rel,thv0t1,porig,dporig,thv0j,rho0j,tj,rdqsdt,    &
         rbeta,ths0,aquad,bquad,cquad,xs1,xs2,excessu,excess0,xsat,    &
         bogbot,bogtop,delbog,expfac,rhos0j,ppen,ufrcbelow,    &
         rmaxfrac,rlwp,qconbelow,qtdef,drage,exnj,psorig
    real*4 erfarg,erfval
    real(r8) latvapocp ! latent heat of vaporization/cpair
    real(r8) leff
    integer :: kpbl2d(pcols)   ! output for uwpbl scheme
    !

    real(r8)               &     ! VARIABLES WHICH DESCRIBE THE ENVIRONMENT
         P0(pver),         &     ! environmental pressure
         PS0(0:pver),      &     ! environmental pressure at full sigma levels
         exner0(pver),     &     ! (p0/p00)**kapa
         exners0(0:pver),  &     ! (ps0/p00)**kapa
         DP(pver),         &     ! environmental layer pressure thickness
         z0(pver),         &     ! environmental height at half-levels
         ZS0(0:pver),      &     ! environmental height at sigma levels
         DZ(pver),         &     ! layer thickness, meters
         rh(pver),         &     ! relative humidity
         thv0bot(pver),    &     ! environmental virtual potential temperature,bottom of layer
         thv0top(pver),    &     ! environmental virtual potential temperature,top of layer
         ssthc0(pver),     &     ! slope of environmental liquid potential temperature
         Q0,               &     ! environmental water vapor mixing ratio
         qct0(pver),       &     ! environmental Total water mixing ratio
         thc0(pver),       &     ! environmental theta_c
         ssQCT0(pver),     &     ! slope of environmental total water mixing ratio
         u0(pver),         &     ! environmental zonal wind speed
         v0(pver),         &     ! environmental meridional wind speed
         u1(pver),         &     ! environmental zonal wind speed
         v1(pver),         &     ! environmental meridional wind speed
         dudp0(2:pver-1),  &     ! environmental zonal wind speed vertical gradient
         dvdp0(2:pver-1)         ! environmental meridional wind speed vertical gradient
    !
    real(r8) plcl(pcols)   ! pressure of lifting condensation level (Pa)
    real(r8) plfc(pcols)   ! pressure of level of free convection (Pa)
    real(r8) cino(pcols)   ! cin (m2/s2)
    real(r8)               &     ! VARIABLES WHICH DESCRIBE THE UPDRAFT
         WU(0:pver),       &     ! updraft vertical velocity at top of layer
         UMF(0:pver) ,     &     ! updraft mass flux at top of layer
         EMF(0:pver) ,     &     ! entrainment mass flux at top of layer
         THCU(0:pver),     &     ! updraft liquid potential temperature at top of layer
         uu(0:pver),       &     ! updraft zonal wind speed
         vu(0:pver),       &     ! updraft meridional wind speed
         qctu(0:pver),     &     ! updraft total water at top of layer
         thcflx(0:pver),   &     ! flux of thc due to convection
         qctflx(0:pver),   &     ! flux of qct due to convection
         umflx(0:pver),    &     ! flux of zonal momentum due to convection
         vmflx(0:pver),    &     ! flux of meridional momentum due to convection
         THVU(0:pver),     &     ! updraft virtual potential temperature at top of layer
         THVUT(0:pver),     &     ! updraft virtual potential temperature at top of layer
         REI(pver),        &     ! updraft mixing rate of with environment
         fer(pver),        &     ! updraft fracional entrainment rate
         fdr(pver),        &     ! updraft fracional detrainment rate
         PPTR(pver),       &     ! updraft production rate of rain (kg/m2/s)
         cwten(pver),      &     ! cloud water export tendency (kg/m2/        s)
         qcto(pver),       &     ! detraining total water
         thco(pver),       &     ! diagnostic thc
         qo(pver),         &     ! detraining water vapor
         qlo(pver),        &     ! detraining liquid water
         qio(pver)               ! detraining cloud ice
    real(r8) :: es(1)                     ! saturation vapor pressure
    real(r8) :: qs(1)                     ! saturation spec. humidity
    real(r8) :: gam(1)                    ! (L/cp)*dqs/dT
    integer  :: status
    !------pzhu-----
    real(r8) a_qvten,a_sten,dissip,a_uvten,a_uvten1,ktoppbl_1
    real(r8) temsrc,esrc,tdsrc,tlcl
    real(r8) thlfs,qtfs,thvfs,yy1,xbuo0,rmfk3 
    real(r8) sfbuoflx,wstar,cstar,cbmf0
    real(r8) tem0(pver),qc0(pver),qi0(pver)
    real(r8) plcl_ini
    real(r8) rain_ch	
    integer id_check, kfree, klimit	
    logical, parameter :: cin_lfc = .false. ! true when lfc is used to define CIN
    !------end of pzhu
    !
    !     For lateral entrainment

    !pzhu---if igauss=1 then
    parameter (cstar=0.385)  ! for estimating w'2 
    parameter (rle = 0.10)   ! for critical stopping distance for entrainment
    parameter (rkm = 16.0)    ! for fractional mixing rate
    parameter (rpen = 5.0)   ! for entrainment efficiency
    parameter (rkfre = 0.05)  ! vertical velocity variance as fraction of tke
    parameter (rmaxfrac = 0.05) ! maximum allowable updraft fraction
    parameter (rbuoy = 1.0)  ! for nonhydrostatic pressure effects on updraft
    parameter (rdrag = 1.0)
    parameter (rain_ch = 1.0e-3) ! threshold for precipitation

    !     For momentum transfer
    parameter (bigc = 0.7)
    !     Some options:
    !     Determine mass flux with cin/gaussian closure vs cin/tke closure
    !     1 = cin/gaussian, 0 = cin/tke
    !-----igauss=1: using tke to compute CIN. 
    !     igauss=2: using W* to compute CIN. 
    !     igauss=0: This is for the mapes-style closure
    !     These yield results similar to my new closure: (rmfk1=0.2821,rmfk2=1.333333)
    parameter(igauss = 1)
    parameter (rmfk1 = 0.3, rmfk2 = 5.0, rmfk3=3.0)
    !      parameter (thvlmj=0.01)    ! thvl minimum jump across inversion
    !
    real(r8) p00
    DATA p00/1.E5/

    latvapocp=latvap/cpair
    kl=pver
    klm=kl-1

    !
    !     Start the big i loop
    !
    DO 2112 I=1,ncol     ! start of big i loop
       cino(i)=0.
       toplev(i)=real(pver+1, r8)
       botlev(i)=0.
       plcl(i)=0.
       plfc(i) = 0.
       rliq(i) = 0.
       prec_cmf(i) = 0.
       do k=1,pver
          sten(i,k)=0.
          qvten(i,k)=0.
          uten(i,k)=0.
          vten(i,k)=0.
          cldfrac(i,k)=0.0
          cldqc(i,k)=0.0
          cmfsl(i,k) = 0
          cmflq(i,k) = 0
          cmfdqr(i,k) = 0
          qc(i,k) = 0
          cwten(k)= 0.0
          pptr(k)= 0.0
          thcflx(k)=0.
          qctflx(k)=0.
          umflx(k)=0.
          vmflx(k)=0.
          cmf(i,k) = 0.0
       enddo
       cmf(i,pver+1) = 0.0
       !        There are numerous exit points to the 2112 loop, generally to save time.
       !
       tscaleh = cush(i)
       cush(i)=-1.
       !
       !_____1. Get the environmental conditions
       !     Model layers are numbered from bottom up!
       !
       zs0(0) = 0.0
       ps0(0) = pint(i,pver+1)
       exners0(0)=(ps0(0)/p00) ** cappa
       DO K=1,pver
          nk=pver-k+1
          p0(k)=pmid(i,nk)
          ps0(k) = pint(i,nk)
          zs0(k) = zint(i,nk)
          dp(k)=PS0(K-1)-PS0(K)
          z0(k)=zmid(i,nk)
          dz(k) = zs0(k)-zs0(k-1)
          qc0(k) = q(i,nk,2)
          qi0(k) = q(i,nk,3)
          qct0(k) = q(i,nk,1)+q(i,nk,2)+q(i,nk,3)
          u0(k)=ub(i,nk)
          v0(k)=vb(i,nk)
          nu = max(min((268. - tb(i,nk) )/20.,1._r8),0._r8)
          leff = (1-nu)*latvap + nu*latsub
          exner0(k)=(p0(k)/p00) ** cappa
          exners0(k)=(ps0(k)/p00) ** cappa
          thc0(k)=(tb(i,nk) - leff*qc0(k)/cpair)/exner0(k)  ! theta_l for enviro.
          tem0(k)=tb(i,nk)
          status = qsat(tb(i,nk),ps0(k),es(1),qs(1),gam(1),1)
          rh(k)=qct0(k) / qs(1)
          umf(k)=0.
          emf(k)=0.
       END DO

       !--get index of PBL height!
       do k = kl-1,2,-1
          if ( (pblht(i)+1.-zs0(k))*(pblht(i)+1.-zs0(k+1)) .lt. 0. ) goto 106
       end do
106    ktoppbl=k
       
       ssthc0b = (thc0(2)-thc0(1))/(p0(2)-p0(1))
       ssqct0b = (qct0(2)-qct0(1))/(p0(2)-p0(1))

       DO K=2,KL
          ssthc0a = (thc0(k)-thc0(k-1))/(p0(k)-p0(k-1))
          if(ssthc0a.gt.0)then
             ssthc0(k-1) = max(0._r8,min(ssthc0a,ssthc0b))
          else
             ssthc0(k-1) = min(0._r8,max(ssthc0a,ssthc0b))
          endif
          ssthc0b = ssthc0a
          ssqct0a = (qct0(k)-qct0(k-1))/(p0(k)-p0(k-1))
          if(ssqct0a.gt.0)then
             ssqct0(k-1) = max(0._r8,min(ssqct0a,ssqct0b))
          else
             ssqct0(k-1) = min(0._r8,max(ssqct0a,ssqct0b))
          endif
          ssqct0b = ssqct0a
       enddo
       !     Wind shear
       do k = 2,kl-1
          dudp0(k) = (u0(k+1)-u0(k-1))/(p0(k+1)-p0(k-1))
          dvdp0(k) = (v0(k+1)-v0(k-1))/(p0(k+1)-p0(k-1))
       END DO
       ssthc0(kl)=ssthc0(kl-1)
       ssqct0(kl)=ssqct0(kl-1)
       do k = 1,kl
          thc0bot = thc0(k)+ssthc0(k)*(ps0(k-1)-p0(k))
          qct0bot  =  qct0(k)+ ssqct0(k)*(ps0(k-1)-p0(k))
          call conden(ps0(k-1),thc0bot,qct0bot,thj,qvj,qlj,qij,id_check,qsat)
          if(id_check.eq.1) go to 2112 
          thv0bot(k)= thj * (1.+zvir*qvj-qlj-qij)

          thc0top = thc0(k)+ssthc0(k)*(ps0(k)-p0(k))
          qct0top  =  qct0(k)+ ssqct0(k)*(ps0(k)-p0(k))
          call conden(ps0(k),thc0top,qct0top,thj,qvj,qlj,qij,id_check,qsat)
          if(id_check.eq.1) go to 2112
          thv0top(k)= thj * (1.+zvir*qvj-qlj-qij)
       enddo
       !

       !_____2. Determine precise height of the inversion. The inversion is considered
       !        to be the top of PBL that is defined at one interface level higher 
       !	 than the interface level with Ri<0. see BG PBL scheme for detail. 

       kinv = ktoppbl+1
       pinv=ps0(kinv-1)
       zinv=zs0(kinv-1)

       !------calculate the environmental virtual potential temperature at INV

       call conden(pinv,thc0(kinv)+ssthc0(kinv)*(ps0(kinv-1)-p0(kinv))  &
            ,qct0(kinv)+ssqct0(kinv)*(ps0(kinv-1)-p0(kinv)),thj,qvj,     &
            qlj,qij,id_check,qsat)
       if(id_check.eq.1) go to 2112
       thvuinv = thj * (1.+zvir*qvj-qlj-qij)

       !
       !_____3. Let's get some source air, and find its LCL
       !
       usrc = u0(ktoppbl)
       vsrc = v0(ktoppbl)

       !-------source air is simply set to be the first layer air, which is different
       !       from the Jim's code in MM5. Numerical tests indicate that such a 
       !	difinition will produce a much stable result, i.e., results are 
       !	insensitive to the model time step interval. 

       kthvmin = 1
       do k=1,ktoppbl
          if(thv0top(k).lt.thv0top(kthvmin))kthvmin = k
       enddo

       qtsrc=qct0(kthvmin)
       thlsrc=thc0(kthvmin)
       esrc=qtsrc*ps0(kthvmin)/100./(qtsrc+epsilo)  ! water vapor pressure
       tdsrc=tmelt/(1-tmelt*rgasv*log(esrc/6.11)/latvap) !dew-point of source air
       temsrc=thlsrc*(ps0(kthvmin)/p00)**cappa        ! temperature of source air
       zlcl=123.5*(temsrc-tdsrc)+zs0(kthvmin)         ! from sea-level
       tlcl=temsrc-0.0098*(zlcl-zs0(kthvmin))
       plcl(i)=ps0(kthvmin)*(tlcl/temsrc)**(1./cappa)
       plcl_ini=plcl(i)

!       if(plcl(i).ge.ps0(kinv)) then
!	  plcl(i)=ps0(kinv)
!       endif
       do k=1,kl
          klcl=k    !KLCL is the layer containing the lcl, i.e., ps0(klcl)<=plcl(i)
          if(ps0(K).le.plcl(i))GOTO 35
       END DO
       GOTO 2112
35     klcl=max(klcl,2)

       do k=klcl,kl
          call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
          thv0lcl = thj * (1.+zvir*qvj-qlj-qij)
          if(thv0lcl.ge.thv0bot(k)) then
             go to 115
          endif
       enddo
       go to 2112  ! no free convection level is found
115    kfree=k
       if(zs0(kfree).gt.3000.) go to 2112  ! convection is not related to PBL
       do k=kfree,kl
          call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
          thv0lcl = thj * (1.+zvir*qvj-qlj-qij)
          if(thv0lcl.lt.thv0bot(k)) then
             go to 116
          endif
       enddo
       go to 75
116    klimit=k
       if(z0(klimit).lt.1500.) go to 2112 ! convective layer is too shallow
       ! and must be associated with
       ! stratus clouds
75     continue


       !-----get environmental properties at KLCL

       thc0lcl = thc0(klcl)+ssthc0(klcl)*(plcl(i)-p0(klcl))
       qct0lcl = qct0(klcl)+ssqct0(klcl)*(plcl(i)-p0(klcl))
       call conden(plcl(i),thc0lcl,qct0lcl,thj,qvj,qlj,qij,id_check,qsat)
       if(id_check.eq.1) go to 2112
       thv0lcl = thj * (1.+zvir*qvj-qlj-qij)
       rho0lcl = plcl(i)/(rair*thv0lcl*(plcl(i)/p00)**cappa)

       !_____4. Determine the convective inhibition (CIN)

       !     Initialize CIN
       CIN = 0.
       cinlcl = 0.
       plfc(i) = 0.

       if( cin_lfc ) then    ! define CIN based on LFC  
          do k = kinv,kl-1
             if(k.eq.klcl-1) then !        Klcl-1 < layer < klcl
                call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
                if(id_check.eq.1) go to 2112
                thvubot=thj * (1.+zvir*qvj-qlj-qij) 
                call conden(plcl(i),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
                if(id_check.eq.1) go to 2112
                thvutop=thj * (1.+zvir*qvj-qlj-qij) 
                call getbuoy(ps0(k),thv0top(k),plcl(i),thv0lcl,      &
                     thvubot,thvutop,plfc(i),cin)
                cinlcl = cin
                thvubot=thvutop
                call conden(ps0(k+1),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
                if(id_check.eq.1) go to 2112
                thvutop= thj * (1.+zvir*qvj-qlj-qij)
                call getbuoy(plcl(i),thv0lcl,ps0(k+1),thv0top(k+1),      &
                     thvubot,thvutop,plfc(i),cin) 	
                if(plfc(i).gt.0.)goto 668 	
             else
                call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat) 
                if(id_check.eq.1) go to 2112
                thvubot=thj* (1.+zvir*qvj-qlj-qij)
                call conden(ps0(k+1),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
                if(id_check.eq.1) go to 2112
                thvutop=thj* (1.+zvir*qvj-qlj-qij)
                call getbuoy(ps0(k),thv0top(k),ps0(k+1),thv0top(k+1),    &
                     thvubot,thvutop,plfc(i),cin)
                if(plfc(i).gt.0.)goto 668 
             endif
          enddo
          !         write(iulog,*) 'No LFC for undilute parcel ascent.  Bailing.'
          !	   go to 2112
          cin=100.
       else
          do k=kinv,klcl-1
             if(k.eq.(klcl-1)) then !        Klcl-1 < layer < klcl
                call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
                if(id_check.eq.1) go to 2112
                thvubot=thj * (1.+zvir*qvj-qlj-qij) 
                call conden(plcl(i),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
                if(id_check.eq.1) go to 2112
                thvutop=thj * (1.+zvir*qvj-qlj-qij) 
                !            call getbuoy(ps0(k),thv0top(k),plcl(i),thv0lcl,      &
                call getbuoy(ps0(k),thv0bot(k+1),plcl(i),thv0lcl,      &
                     thvubot,thvutop,plfc(i),cin)
                cinlcl = cin
             else
                call conden(ps0(k),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
                if(id_check.eq.1) go to 2112 
                thvubot=thj* (1.+zvir*qvj-qlj-qij)
                call conden(ps0(k+1),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
                if(id_check.eq.1) go to 2112
                thvutop=thj* (1.+zvir*qvj-qlj-qij)
                !            call getbuoy(ps0(k),thv0top(k),ps0(k+1),thv0top(k+1),    &
                call getbuoy(ps0(k),thv0bot(k+1),ps0(k+1),thv0top(k+1),    &
                     thvubot,thvutop,plfc(i),cin)
             endif
          enddo

       endif  ! CIN has been estimated 

668    continue

       !_____5. Calculate updraft cloud base mass flux
       dpsum=0.
       tkeavg = 0.
       do K=1,ktoppbl
	  dpsum=dpsum+dp(k)
	  tkeavg = tkeavg + dp(k)*tkeb(i,pver-k+1)
       enddo
       tkeavg = tkeavg/dpsum
       tkeavg = max(0.1_r8,tkeavg)
       if(igauss.eq.0)then
          !           Use cin and pbl tke
          cbmf = rmfk1*rho0lcl*sqrt(tkeavg)*exp(-rmfk2*cin/tkeavg)
          !           Updraft vertical velocity at release height depends on tke
          wexp = rmfk3*sqrt(tkeavg)
       elseif(igauss.eq.1)then
          !           Use cin and gaussian distribution of w
          wcrit = sqrt(2 * cin * rbuoy)
          sigmaw = sqrt(rkfre*tkeavg)
          cbmf =  rho0lcl * sigmaw / 2.5066 * exp(-0.5*((wcrit/sigmaw)**2))
          !           Diagnose updraft fraction        ! sqrt(2.) = 1.4142
          erfarg=wcrit / (1.4142 * sigmaw)
          if(erfarg.lt.20.)then
             erfval=erfc(erfarg)
             ufrc = min(rmaxfrac,0.5_r8*erfval)
          else
             ufrc = 0.
          endif
          if(ufrc.gt.0.001)then
             !              Diagnose expected value of cloud base vertical velocity
             wexp =  cbmf / rho0lcl / ufrc
          else
             goto 2112
          endif
       endif
       !
       !_____6.  Determine release height and the vertical velocity at the LCL
       !         Must specify:  krel, prel, thv0rel, thv0t1, thvurel, wtw
       if(plcl(i).gt.pinv)then
          krel=kinv
          prel=pinv
          thv0rel = thvuinv
          thv0t1  = thv0top(kinv)
       else
          krel=klcl
          prel=plcl(i)
          thv0rel = thv0lcl
          thv0t1 = thv0top(krel)
       endif
       wexp=min(wexp,50._r8)
       wtw = wexp * wexp - 2 * cin * rbuoy
       !----------------------------------------------
       if(wtw.le.0.)then
          goto 2112
       endif
       wrel = sqrt(wtw)
       !     For these (krel-1) represents the bottom of the updraft
       psorig = ps0(krel-1)  
       ps0(krel-1) = prel
       thcu(krel-1)= thlsrc
       qctu(krel-1)= qtsrc
       !         thvu(krel-1)= thvuinv  !theta_v is not conserved. 
       call conden(ps0(krel-1),thlsrc,qtsrc,thj,qvj,qlj,qij,id_check,qsat)
       if(id_check.eq.1) go to 2112
       thvu(krel-1)=thj*(1.+zvir*qvj-qlj-qij)
       uu(krel-1)  = usrc
       vu(krel-1)  = vsrc
       umf(krel-1) = cbmf
       wu(krel-1)  = wrel

       !     And for these (krel) represents the first partial updraft layer
       !     The first ones are special because they must be restored later
       porig = p0(krel)
       p0(krel) = 0.5*(prel+ps0(krel))
       dporig = dp(krel)
       dp(krel) = prel - ps0(krel)
       !
       thv0bot(krel) = thv0rel
       thv0top(krel) = thv0t1
       if(krel.eq.kinv)then
          thc0(krel) = thc0(kinv)
          qct0(krel) = qct0(kinv)
       else
          thc0(krel) = thc0(krel)+ssthc0(krel)*(p0(krel)-porig)
          qct0(krel) = qct0(krel)+ssqct0(krel)*(p0(krel)-porig)
       endif
       !*******************************************************************
       !                                                                  *
       !_____7.     Compute updraft properties above the LCL
       !                                                                  *
       !*******************************************************************
       !
       let=krel
       scaleh = tscaleh
       if(tscaleh.lt.0.0)scaleh = 1000.
       DO 60 k=krel,klm
          km1=k-1
          !-----A.  Entrainment and Detrainment
          !     first, to determine fraction (xsat) of mixture that is to be detrained out 
          !     of clouds, i.e., the mixture with negative buoyancy. We consider a thin 
          !     layer between two interfaces, so using mid-point value to represent the 
          !     mean value of the layer. The properties of updraft at midpoint is assumed
          !     to be undiluted from the lower interface.  

          !-----calculate fraction of mixture that is just saturated

          status = qsat(thcu(km1)*exner0(k),p0(k),es(1),qs(1),gam(1),1)
          excessu = qctu(km1) - qs(1)
          excessu = max(excessu,0._r8)

          status = qsat(thc0(k)*exner0(k),p0(k),es(1),qs(1),gam(1),1)
          excess0 = qct0(k) - qs(1)

          if(excessu*excess0.le.0)then
             xsat = -excessu/(excess0-excessu)
          else
             xsat = 1.0
          endif
          thlfs=(1.-xsat)*thcu(km1)+xsat*thc0(k)
          qtfs=(1.-xsat)*qctu(km1)+xsat*qct0(k)
          thvfs=thlfs*(1+zvir*qtfs)


          call conden(p0(k),thcu(km1),qctu(km1),thj,qj,qlj,qij,id_check,qsat)
          if(id_check.eq.1) go to 2112
          thvj=thj*(1.+zvir*qj-qlj-qij)   ! theta_v of updraft
          call conden(p0(k),thc0(k),qct0(k),thj,qj,qlj,qij,id_check,qsat)
          if(id_check.eq.1) go to 2112
          thv0j=thj*(1.+zvir*qj-qlj-qij)  ! theta_v of environment
          rho0j = p0(k)/(rair*thv0j*(p0(k)/p00)**cappa)

          !-----calculate fraction of mixture with zero buoyancy
          if(thvfs.ge.thv0j) then
	     xbuo0=xsat
          elseif(thvj.le.thv0j) then
	     xbuo0=0.
          else
!	     xbuo0=xsat*(thv0j-thvfs)/(thvj-thvfs)
	     xbuo0=xsat*(thvj-thv0j)/(thvj-thvfs)
          endif

          !-----calculate fraction of mixture with negative buoyancy but can 
          !     penetrate a critical distance lc=rle*scaleh
          if(thvfs.ge.thv0j.or.xsat.le.0.05) then
	     xs=xsat ! mixture has to be saturated
          else
	     aquad = wu(km1)**2
	     bquad=-(2*wu(km1)**2+2*rbuoy*gravit*rle*scaleh*(thvj-thvfs)/thv0j/xsat)  
	     cquad=wu(km1)**2-2*rbuoy*gravit*rle*scaleh*(1-thvj/thv0j)
             call roots(aquad,bquad,cquad,xs1,xs2)
	     xs=min(xs1,xs2)
          endif
          xs=min(xs,xsat)
          xs=max(xbuo0,xs)
          xs=min(1._r8,xs)

          ee2 = xs**2
          ud2 = 1 - 2*xs + xs**2
          rei(k) = rkm/scaleh/gravit/rho0j
          fer(k)=rei(k)*ee2
          fdr(k)=rei(k)*ud2

          !-----B.  Calculate the mass flux
          umf(k)=umf(km1)*exp(dp(k)*(fer(k)-fdr(k)))
          emf(k)=0.0
          !	write(iulog,'(7e14.6)') xs,fer(k),fdr(k),umf(k),xsat,excess0,excessu

          !-----C. Now thermodynamics for the dilute plume
          thcu(k)=thc0(k)-(thc0(k)-thcu(km1))*exp(-fer(k)*dp(k))
          qctu(k)=qct0(k)-(qct0(k)-qctu(km1))*exp(-fer(k)*dp(k))
          if(fer(k)*dp(k).lt.1.e-4)then
             uu(k)=uu(km1) - dudp0(k)*dp(k)
             vu(k)=vu(km1) - dvdp0(k)*dp(k)
          else
             uu(k)=u0(k)-bigc*dudp0(k)/fer(k)-exp(-fer(k)*dp(k))*        &
                  (u0(k) - bigc*dudp0(k)/fer(k) - uu(km1))
             vu(k)=v0(k)-bigc*dvdp0(k)/fer(k)-exp(-fer(k)*dp(k))*        &
                  (v0(k) - bigc*dvdp0(k)/fer(k) - vu(km1))
          endif

          !     Precip at the flux level
          call conden(ps0(k),thcu(k),qctu(k),thj,qj,qlj,qij,id_check,qsat)
          if(id_check.eq.1) go to 2112
          thvu(k) = thj * (1.+zvir*qj-qlj-qij)
          rho0j = p0(k)/(rair*thv0j*(p0(k)/p00)**cappa)
          cwten(k) = min((qlj+qij),0.5*q(i,pver+1-k,1))/ztodt
          cldqc(i,pver-k+1)=cwten(k)*ztodt
          pptr(k) = min(rain_ch*umf(k)*(qlj+qij+qc0(k)+qi0(k))/rho0j,cwten(k))
          !
          !-----D.  Calculate vertical velocity
          !
          bogbot = (thvu(km1)/thv0bot(k) - 1.)
          bogbot =  bogbot*rbuoy
          bogtop = (thvu(k)/thv0top(k) - 1.)
          bogtop =  bogtop*rbuoy

          if(bogbot.gt.0.and.bogtop.gt.0)let = k
          !            if(bogbot.gt.0)let = k
          delbog = bogtop - bogbot
          drage = fer(k) * ( 1. + rdrag )
          expfac = exp(-2.*drage*dp(k)) ! dp(k) = - delta p
          if(drage*dp(k).gt.1.e-3)then
             wtw = wtw*expfac + (delbog + (1.-expfac) *               &
                  (bogbot+delbog/(-2.*drage*dp(k))))/(rho0j*drage)
          else
             wtw = wtw + dp(k) * (bogbot+bogtop)/rho0j
          endif
          if(wtw.le.0.) goto 65
          wu(k) = sqrt(wtw)
          if(wu(k).gt.100.)then
             write(iulog,*) 'big wu',bogbot,bogtop,expfac,fer(k)
             stop
          endif

          rhos0j = ps0(k) /(rair*0.5*(thv0bot(k+1)+thv0top(k))*exners0(k))
          ufrc = umf(k)/(rhos0j*wu(k))
          if(ufrc.gt.rmaxfrac)then
             ufrc = rmaxfrac
             fdr(k)= fer(k) -log(rmaxfrac*rhos0j*wu(k)/umf(km1))/dp(k)
             umf(k)= rmaxfrac * rhos0j * wu(k)
          endif
60     end do
       !     End of Updraft Loop
       !...  ltop is the first level with negative vertical velocity at top
65     ltop = k
       umf(ltop)=0.0
       cldhgt=z0(ltop)-zlcl

       !        convection too deep 
       if(cldhgt.ge.4.e3)then

          goto 2112 ! Cloud too deep
       endif

       !        Calculate convective scale height
       cush(i)=z0(ltop)

       !     Calculate penetrative entrainment
       emf(ltop)=0.0
       do k=ltop-1,let,-1
          qcto(k+1)=qct0(k)
          !            qo(k+1)=qcto(k+1)
          !            qlo(k+1)=0.0
          !            qio(k+1)=0.0
          rhos0j = ps0(k) /(rair*0.5*(thv0bot(k+1)+thv0top(k))*exners0(k))
          if(k.eq.ltop-1)then

             !         7. Calculate ppen

             bogbot = (thvu(k)/thv0bot(ltop) - 1.)
             bogbot =  bogbot*rbuoy
             bogtop = (thvu(ltop)/thv0top(ltop) - 1.)
             bogtop =  bogtop*rbuoy
             aquad = (bogtop - bogbot) / (ps0(ltop)-ps0(k))
             bquad = 2*bogbot
             cquad = -wu(k)*ps0(k)/(rair*thv0bot(ltop)*exners0(k))
             call roots(aquad,bquad,cquad,xs1,xs2)
             if(xs1.le.0..and.xs2.le.0.)then
                ppen = max(xs1,xs2)
             else
                ppen = min(xs1,xs2)
             endif
             ppen = min(0._r8,max(-dp(k+1),ppen))
             if(xs1.eq.-9.99e33.or.xs2.eq.-9.99e33)ppen=0.
             !_____8. Calculate returning mass flux
             emf(k)= max(umf(k)*ppen*rei(ltop)*rpen,-0.1*rhos0j)
             thcu(k)=thc0(ltop)+ssthc0(ltop)*(ps0(k)-p0(ltop))
             qctu(k)=qct0(ltop)+ssqct0(ltop)*(ps0(k)-p0(ltop))
          else
             emf(k)=max(emf(k+1)-umf(k)*dp(k+1)*rei(k+1)*rpen,-0.1*rhos0j)
             thcu(k)=(thcu(k+1)*emf(k+1)+thc0(k+1)*(emf(k)-emf(k+1)))/emf(k)
             qctu(k)=(qctu(k+1)*emf(k+1)+qct0(k+1)*(emf(k)-emf(k+1)))/emf(k)
          endif
          !            umf(k)=0.0
       enddo
       !     Restore special values
       ps0(krel-1) = psorig
       p0(krel) = porig
       dp(krel) = dporig
       !*****************************************************************
       !     Done describing the updraft
       !*****************************************************************
       !_____9.  Output to model
       !     Calculate Fluxes of heat, moisture, momentum
       do k=1,pver
          nk = pver+1-k
          cmf(i,nk) = umf(k)
          thcflx(k)=0.
          qctflx(k)=0.
          umflx(k)=0.
          vmflx(k)=0.
       enddo

       qctflx(0)=0.0
       thcflx(0)=0.0
       umflx(0)=0.0
       vmflx(0)=0.0
       dpsum = 0.0
       do k = 1, krel-1
          dpsum = dpsum + dp(k)
       enddo
       qtdef = max(0._r8,umf(krel)*(qctu(krel) - qct0(krel)))
       yy1 = min(0._r8,umf(krel)*(thcu(krel) - thc0(krel)))
       do k=2,krel-1
          !            thcflx(k)=0.0
          !           qctflx(k)=0.0
          thcflx(k)=thcflx(k-1) + yy1*dp(k)/dpsum
          qctflx(k)=qctflx(k-1) + qtdef*dp(k)/dpsum
          umflx(k)=0.0
          vmflx(k)=0.0
          !            pptr(k)=0.0
          !	    cwten(k)=0.0
       enddo
       !         do k = krel,ltop-1
       do k = krel,ltop
          kp1 = k+1
          thcflx(k)= umf(k) *(thcu(k)-(thc0(kp1)+ssthc0(kp1)*(ps0(k)-p0(kp1)))) +   &
               emf(k) * (thcu(k)-(thc0(k)+ssthc0(k)*(ps0(k)-p0(k))))
          qctflx(k)= umf(k) *(qctu(k)-(qct0(kp1)+ssqct0(kp1)*(ps0(k)-p0(kp1)))) +   &
               emf(k) * (qctu(k)-(qct0(k)+ssqct0(k)*(ps0(k)-p0(k))))
          umflx(k) =umf(k) * (  uu(k)-  u0(kp1)) + emf(k) * (  u0(kp1)-  u0(k))
          vmflx(k) =umf(k) * (  vu(k)-  v0(kp1)) + emf(k) * (  v0(kp1)-  v0(k))

          !----using top-interface
          !               thc0top=thc0(k)+ssthc0(k)*(ps0(k)-p0(k))
          !               qct0top=qct0(k)+ssqct0(k)*(ps0(k)-p0(k))
          !            thcflx(k)= umf(k) *(thcu(k)-thc0top)
          !            qctflx(k)= umf(k) *(qctu(k)-qct0top)

          !----using bot-interface
          !               thc0bot=thc0(kp1)+ssthc0(kp1)*(ps0(k)-p0(kp1))
          !               qct0bot=qct0(kp1)+ssqct0(kp1)*(ps0(k)-p0(kp1))
          !            thcflx(k)= umf(k) *(thcu(k)-thc0bot)
          !            qctflx(k)= umf(k) *(qctu(k)-qct0bot)

          !            umflx(k) =umf(k) * (uu(k)-u0(kp1))
          !            vmflx(k) =umf(k) * (vu(k)-v0(kp1))
       enddo
       thcflx(ltop+1) = 0.0
       qctflx(ltop+1) = 0.0
       umflx(ltop+1) = 0.0
       vmflx(ltop+1) = 0.0
       !
       !        Calculate model tendencies

       do k = ltop+1,1,-1
          km1 = k-1
          qctten=(qctflx(km1)-qctflx(k))*gravit/dp(k)
          nk = pver+1-k
          if((q(i,nk,1)+(qctten-cwten(k))*ztodt).lt.1.e-12) then
             go to 2112
          endif
       enddo

       a_qvten=0.
       a_sten=0.
       do k = pver,1,-1
          km1 = k-1
          nk = pver+1-k
          uten(i,nk) = (umflx(km1)-umflx(k))*gravit/dp(k)
          vten(i,nk) = (vmflx(km1)-vmflx(k))*gravit/dp(k)
          u1(k)=u0(k)+uten(i,nk)*ztodt
          v1(k)=v0(k)+vten(i,nk)*ztodt
       enddo
       !         do k=pver,2,-1
       do k=ltop+1,2,-1
          km1 = k-1
          nk = pver+1-k
          dissip=(umflx(km1)*(u1(km1)-u1(k)+u0(km1)-u0(k))+    &
               vmflx(km1)*(v1(km1)-v1(k)+v0(km1)-v0(k)))*0.5*ztodt

          !           dissip=(umflx(km1)-umflx(k))*(u1(k)+u0(k))+    &
          !                  (vmflx(km1)-vmflx(k))*(v1(k)+v0(k))
          !           dissip=dissip*0.5*ztodt

          thcten = ( thcflx(km1) - thcflx(k) ) *gravit/dp(k)
          qctten = ( qctflx(km1) - qctflx(k) ) *gravit/dp(k)
          qvten(i,nk) = qctten-cwten(k)
          qc(i,nk) =cwten(k)-pptr(k)
          cmfdqr(i,nk) = pptr(k)
          sten(i,nk) = cpair*thcten + latvap*cwten(k) + dissip*gravit/dp(k)/ztodt
          !
          cmfsl(i,nk) = cpair*thcflx(k)
          cmflq(i,nk) = qctflx(k)*latvap ! total water flux
          !     precipitation at bottom of layer
       enddo

       !     Diagnostic Outputs, and convective cloud fields for radiation scheme
       rhos0j = ps0(krel-1) / (rair*0.5*(thv0bot(krel)+thv0top(krel-1))*   &
            exners0(krel-1))
       ufrcbelow = cbmf/(rhos0j*wu(krel-1))
       ufrcbelow = cbmf
       qconbelow = 0.
       do k = krel, ltop-1
          nk = pver - k + 1
          rhos0j = ps0(k) / (rair*0.5*(thv0bot(k+1)+thv0top(k))*exners0(k))
          ufrc = umf(k)/(rhos0j*wu(k))
          call conden(ps0(k),thcu(k),qctu(k),thj,qj,qlj,qij,id_check,qsat)
          if(id_check.eq.1) go to 2112
          cldfrac(i,nk) = 0.5*(ufrcbelow + ufrc) ! 2 * 0.5
          ufrcbelow = ufrc
          qconbelow = qlj+qij
       enddo
       toplev(i) = pver-ltop+1
       botlev(i) = pver-krel+1

       do k = 1, pver
          rliq(i) = rliq(i) + qc(i,k)*pdel(i,k)/gravit
          prec_cmf(i) = prec_cmf(i) + cmfdqr(i,k)*pdel(i,k)/gravit  
       end do
       rliq(i) = rliq(i) /1000.
       prec_cmf(i) = prec_cmf(i) /1000.

2112 end do

    return 
  END subroutine compute_uw_conv

  subroutine conden(p,thc,qt,th,qv,ql,qi,id_check,qsat)
    implicit none
    integer, external     :: qsat
    real(r8), intent(in)  :: p
    real(r8), intent(in)  :: thc
    real(r8), intent(in)  :: qt
    real(r8), intent(out) :: th
    real(r8), intent(out) :: qv
    real(r8), intent(out) :: ql
    real(r8), intent(out) :: qi
    integer,  intent(out) :: id_check
    real(r8) :: p00
    real(r8) :: tc
    real(r8) :: exn
    real(r8) :: leff,nu,qc,temps,tc1
    integer  :: iteration
    real(r8) :: es(1)                     ! saturation vapor pressure
    real(r8) :: qs(1)                     ! saturation spec. humidity
    real(r8) :: gam(1)                    ! (L/cp)*dqs/dT
    integer  :: status
    DATA p00/1.E5/
    
    exn = (p/p00)**cappa
    tc = thc * exn
    nu = max(min((268-tc)/20.,1._r8),0._r8)
    leff = (1-nu)*latvap + nu*latsub

    temps = tc
    status = qsat(temps,p,es(1),qs(1),gam(1),1)
    
    if(qs(1).gt.qt) then
       id_check=0
    else
       do iteration = 1,20
          temps = temps + ((tc-temps)*cpair/leff + (qt -qs(1)))/        &
               (cpair/leff+epsilo*leff*qs(1)/rair/temps/temps)
          ! temps = temps + ((tc-temps)*cpair/leff + (qt -qs(1)))/        &
          !     (cpair/leff+1.5*leff*qs(1)/rgasv/temps/temps)
          status = qsat(temps,p,es(1),qs(1),gam(1),1)
       enddo
       
       tc1=temps-leff/cpair*(qt-qs(1))
       if(abs(tc1-tc).lt.1.0) then
          id_check=0
       else
          id_check=1
       endif
    endif
    qc = max(qt-qs(1),0._r8)
    qv = qt - qc
    ql = (1-nu)*qc
    qi = nu*qc
    th = temps/exn

    return
  end subroutine conden


  subroutine getbuoy(pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin)
    implicit none
    real(r8) pbot,thv0bot,ptop,thv0top,thvubot,thvutop,plfc,cin,frc
    real(r8) p00
    DATA p00/1.E5/

    if(thvubot.gt.thv0bot.and.thvutop.gt.thv0top)then
       plfc = pbot
       return
    elseif(thvubot.le.thv0bot.and.thvutop.le.thv0top)then ! got cin
       cin = cin - ((thvubot/thv0bot-1.)+(thvutop/thv0top-1.)) *(pbot-ptop)/     &
            (pbot/(rair*thv0bot*(pbot/p00)**cappa) + ptop/(rair*thv0top*(ptop/p00)**cappa))
    elseif(thvutop.le.thv0top.and.thvubot.gt.thv0bot)then ! fractional cin
       frc = (thvutop/thv0top-1.) /((thvutop/thv0top-1.)-(thvubot/thv0bot-1.))
       cin = cin - (thvutop/thv0top-1.) * ( (ptop + frc * (pbot-ptop)) - ptop)/   &
            (pbot/(rair*thv0bot*(pbot/p00)**cappa) + ptop/(rair*thv0top*(ptop/p00)**cappa))
    else                      ! last time through
       frc = (thvubot/thv0bot-1.) /((thvubot/thv0bot-1.) - (thvutop/thv0top-1.))
       plfc = pbot - frc * (pbot-ptop)
       cin = cin - (thvubot/thv0bot-1.) * (pbot-plfc)/                            &
            (pbot/(rair*thv0bot*(pbot/p00)**cappa) + ptop/(rair*thv0top*(ptop/p00)**cappa))
    endif
    if(cin.lt.0)stop
    return
  end subroutine getbuoy

  subroutine roots(a,b,c,r1,r2)
    implicit none
    real(r8) a,b,c,r1,r2,q

    if(a.eq.0)then            ! form b*x + c = 0
       if(b.eq.0)then         ! failure: c = 0
          r1 = -9.99e33
       else                   ! b*x + c = 0
          r1 = -c / b
       endif
       r2 = r1
    else
       if(b.eq.0.)then        ! form a*x**2 + c = 0
          if(a*c.gt.0.)then   ! failure: x**2 = -c/a < 0
             r1 =  -9.99e33
          else                ! x**2 = -c/a
             r1 = sqrt(-c/a)
          endif
          r2 = -r1
       else 
          if((b**2 - 4*a*c).lt.0.)then ! failure, no real(r8) roots
             r1 =  -9.99e33
             r2 = -r1
          else
             q = - 0.5 * ( b + sign(1._r8,b) * sqrt(b**2 - 4*a*c) )
             r1 = q/a
             r2 = c/q
          endif
       endif
    endif
    return
  end subroutine roots
end module uw_conv
