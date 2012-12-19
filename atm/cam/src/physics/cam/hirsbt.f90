
    module hirsbt

! ----- Modules -----
    use shr_kind_mod, only: r8 => shr_kind_r8
    use physconst, only: gravit, rair
    use ppgrid
    use hirsbtpar

    implicit none

! public subroutines
    public :: hirsrtm, hirsbt_init

    contains

! -------------------

    subroutine hirsrtm(lchnk, ncol, pp, tt, rmix, co2mmr, o3mix, ts, oro, &
                       tb_ir, britemp)

! mji hirsrtm
! Code includes modifications for F90 formatting and for application
! to NCAR CAM based on the original code in CCM3.
! M. J. Iacono, AER Inc., May 2004
! Further structural revisions for CAM3.5
! M. J. Iacono, AER Inc., April 2008

!      SUBROUTINE HIRSRTM(PP, TT, RMIX, O3MIX, TS, ORO,
!     $                   TB_IR,BRITEMP)
! -- ------------------------------------------------------------------2
!
!                  ***       VERSION 2.0        ***
!
! This subroutine calculates brightness temperatures at the top of the
! atmosphere for 7 TOVS/HIRS channels (2,4,6,8,10,11,12) and 4 TOVS/MSU 
! channels (1,2,3,4).
!
! On Input:
!    pp      -  level pressure (hPa) at which fluxes are evaluated
!               from top of the atmosphere to surface
!    tt      -  layer temperatures (K)
!    rmix    -  layer H2O mass mixing ratio in kg/kg
!    co2mix  -  layer CO2 mass mixing ratio kg/kg 
!    o3mix   -  layer ozone mass mixing ratio kg/kg 
!    ts      -  surface temperature (K)
!    oro     - land-sea flag (sea=0, land=1)
!
! On Ouput:
!    tb_ir   -  infrared brightness temperatures
!    britemp -  microwave brightness temperatures
!
!
!  A flag to include the 4 MSU channels can be switched on (msu_flag=1)
!  and off (msu_flag=0). To decrease the amount of computer time, the
!  microwave routine is changed to a lookup table with almost the same 
!  accuracy as the original routine.
!
! **  last revised 3/31/97 by Richard Engelen **
!
!   This version differs from original version:
!
!     1.  New NOAA10 coefficients
!     2.  Continuum added
!     3.  Any level exceeding 100% RH is changed to 100% RH
!     4.  New channels (2,4,6,8,10,11,12)
!
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: pp(pcols,pverp)      ! Pressure for each column and layer (hPa)
   real(r8), intent(in) :: tt(pcols,pver)       ! Temperature of each column and layer (K)
   real(r8), intent(in) :: ts(pcols)            ! Surface temperature (K)
   real(r8), intent(in) :: rmix(pcols,pver)     ! Water vapor mass mixing ratio (kg/kg)
   real(r8), intent(in) :: co2mmr(pcols)        ! CO2 mass mixing ratio (kg/kg)
   real(r8), intent(in) :: o3mix(pcols,pver)    ! Ozone mass mixing ratio (kg/kg)
   real(r8), intent(in) :: oro(pcols)           ! Land surface flag, sea=0, land=1
!
! Output arguments
!
   real(r8), intent(out) :: britemp(pcols,pnf_msu)   ! HIRS brightness temperatures
   real(r8), intent(out) :: tb_ir(pcols,pnb_hirs)    ! MSU brightness temperatures

!------------------------------Local variables--------------------------

!      real(r8), parameter :: rmwco2 = mwco2/mwdry   ! ratio of molecular weights of co2 to dry air

      real(r8)  uh2o(pver)
      real(r8)  uo3(pver)
      real(r8)  a_hirs(pnb_hirs)
      real(r8)  b_hirs(pnb_hirs)
      real(r8)  xh2o(pnb_hirs)
      real(r8)  yh2o(pnb_hirs)
      real(r8)  xo3(pnb_hirs)
      real(r8)  yo3(pnb_hirs)
      real(r8)  xco2(pnb_hirs)
      real(r8)  yco2(pnb_hirs)
      real(r8)  b_ir(pnb_hirs)
      real(r8)  otrans(pnb_hirs)
      real(r8)  tband(pnb_hirs)
      real(r8)  scoef10(pnb_hirs)
      real(r8)  fcoef10(pnb_hirs)
      real(r8)  cwn(pnb_hirs)
      real(r8)  dtrans(pnb_hirs)
      real(r8)  rad_lay(pnb_hirs)
      real(r8)  radir(pnb_hirs)
      real(r8)  rad(pnb_hirs)
      real(r8)  rad2(pnb_hirs)
      real(r8)  rad3(pnb_hirs)
      real(r8)  refl(pnb_hirs)
      real(r8)  otr_mw(pnf_msu)
      real(r8)  tau(pnf_msu)
      real(r8)  trans(pnf_msu)
      real(r8)  otr_mw2(pnf_msu)
      real(r8)  tau2(pnf_msu)
      real(r8)  trans2(pnf_msu)
!     real(r8)  britemp1(pcols,pnf_msu)
!     real(r8)  britemp2(pcols,pnf_msu)
!     real(r8)  britemp3(pcols,pnf_msu)
      real(r8)  freq(pnf_msu)
      real(r8)  upath_h2o
      real(r8)  upath_co2
      real(r8)  upath_o3
      real(r8)  ucont1
      real(r8)  ucont2
      real(r8)  ppath_h2o
      real(r8)  ppath_co2
      real(r8)  ppath_o3
      real(r8)  sfctemp                     ! surface temperature
      real(r8)  dp                          ! pressure depth of layer
      real(r8)  tlay                        ! layer temperature
      real(r8)  play                        ! pressure
      real(r8)  rlay                        ! water vapor mixing ratio
      real(r8)  o3lay                       ! ozone mixing ratio.
      real(r8)  rhoair                      ! layer density
      real(r8)  e                           ! saturation pressure
      real(r8)  delz                        ! height of layer
      real(r8)  dp2                         ! pressure depth of layer below
      real(r8)  tlay2
      real(r8)  play2
      real(r8)  rlay2
      real(r8)  rhoair2
      real(r8)  e2
      real(r8)  delz2
      real(r8)  b_mw
      real(r8)  b_mw2
      real(r8)  dtr_mw
      real(r8)  uco2(pver)
!      real(r8)  uco2
      real(r8)  pw
      real(r8)  psc_h2o                   ! partial pressure of h2o
      real(r8)  psc_co2                   ! partial pressure of co2
      real(r8)  psc_o3                    ! partial pressure of o3
      real(r8)  t_cont
      real(r8)  t_h2o
      real(r8)  t_co2
      real(r8)  t_o3
      real(r8)  abs(pnf_msu)              ! absorption
      real(r8)  abs2(pnf_msu)
      real(r8)  dtr_mw2
!
!      real(r8) t_malk, plnck, btemp
!      external t_malk, plnck, btemp
      integer  i, ib, if, icol    ! loop control

!------------------------------data statments---------------------------
      data scoef10/3.09,2.64,2.28,1.004,0.429,1.945,11.95/
      data fcoef10/.00434,.00301,.0018,4.33e-5,0.000393,0.0738,1.110/
      data a_hirs/.0183,-.00203,.0653,0.21797,0.29846,0.04612,0.06453/
      data b_hirs/.99992,.99994,.9998,.99957,.9996,.99963,1.0006/
      data cwn/680.23,704.33,733.13,899.5,1224.07,1363.32,1489.42/
      data freq/50.31,53.73,54.96,57.95/


      data xh2o/0.41,0.52,0.17,0.05,0.68,5.02,17.24/
      data yh2o/0.70,1.41,0.14,0.02,1.09,125.9,1194.5/

      data xco2/64.32,13.29,5.04,0.05,0.03,0.25,0.01/
      data yco2/1325.2,121.01,14.36,0.01,0.03,0.02,0.01/

      data xo3/32.83,28.9,27.44,0.01,0.68,0.37,0.01/
      data yo3/77.44,66.7,67.44,1.6,4.17,0.07,0.01/

      real(r8) grav, r
!      real(r8) grav, rco2, r
!      data grav,rco2,r/9.8,0.523e-3,287./
      real(r8) eps
      data eps/0.5/

! use values for constants consistent with radiation code.
      grav = gravit
      r = rair
!      rco2 = co2vmr*rmwco2

      do icol=1,ncol

        sfctemp=ts(icol)

        do ib=1,pnb_hirs
          radir(ib)=0.
        enddo

        upath_h2o=0.0
        upath_co2=0.0
        upath_o3=0.0
        ucont1=0.
        ucont2=0.
        ppath_h2o=0.0
        ppath_co2=0.0
        ppath_o3=0.
!       if(msu_flag.eq.1) then
          do if=1,pnf_msu
             tau(if)=0.
             trans(if)=0.
             rad(if)=0.
             otr_mw(if)=1.
             tau2(if)=0.
             trans2(if)=0.
             rad2(if)=0.
             otr_mw2(if)=1.
          end do
!       endif
        do ib=1,pnb_hirs
           tband(ib)=1.0
           otrans(ib)=1.
          end do


        do i=1,pver
      
          dp=pp(icol,i+1)-pp(icol,i)
          tlay=tt(icol,i)
          play=sqrt(pp(icol,i)*pp(icol,i+1))
          rlay=rmix(icol,i)
          o3lay=o3mix(icol,i)
          
          rhoair=play*100./(r*tt(icol,i))
          e = play*rlay/(rlay+0.6220)
          delz=dp*0.1/(rhoair*grav)
          
!         if(msu_flag.eq.1) then
            dp2=pp(icol,pver+2-i)-pp(icol,pver+1-i)
            tlay2=tt(icol,pver+1-i)
            play2=sqrt(pp(icol,pver+1-i)*pp(icol,pver+2-i))
            rlay2=rmix(icol,pver+1-i)

            rhoair2=play2*100./(r*tt(icol,pver+1-i))
            e2 = play2*rlay2/(rlay2+0.6220)
            delz2=dp2*0.1/(rhoair2*grav)
          
!
!         microwave transfer
!
            call lookup ( play-e, tlay, abs)
            call lookup (play2-e2, tlay2, abs2)
            do if=1,pnf_msu
              tau(if)=tau(if)+abs(if)*delz
              trans(if) = exp(-tau(if))
              b_mw = 1.47445e-23*freq(if)**3/(exp(0.047981*freq(if) &
                     /tlay)-1)
              dtr_mw = otr_mw(if)-trans(if)
              rad(if)=rad(if) + b_mw*dtr_mw
              otr_mw(if)=trans(if)
           
              tau2(if)=tau2(if)+abs2(if)*delz2
              trans2(if) = exp(-tau2(if))
              b_mw2 = 1.47445e-23*freq(if)**3/(exp(0.047981*freq(if) &
                   /tlay2)-1)
              dtr_mw2 = otr_mw2(if)-trans2(if)
              rad2(if)=rad2(if) + b_mw2*dtr_mw2
              otr_mw2(if)=trans2(if)
            end do
!         endif

!
!                ir transfer
!
          uh2o(i)=rlay*100.*dp/grav
          uo3(i)=o3lay*100.*dp/grav
          uco2(i)=co2mmr(icol)*100.*dp/grav
!          uco2=rco2*100.*dp/grav
          pw=play*rlay*28.97/18.
          ucont1=ucont1+uh2o(i)*((play-pw)/1013.)*(296./tlay)
          ucont2=ucont2+uh2o(i)*(pw/1013.)*(296./tlay)

          upath_h2o=upath_h2o+uh2o(i)
          ppath_h2o=ppath_h2o+uh2o(i)*play
          psc_h2o=ppath_h2o/upath_h2o

          upath_co2=upath_co2+uco2(i)
          ppath_co2=ppath_co2+uco2(i)*play
!          upath_co2=upath_co2+uco2
!          ppath_co2=ppath_co2+uco2*play
          psc_co2=ppath_co2/upath_co2

          upath_o3=upath_o3+uo3(i)
          ppath_o3=ppath_o3+uo3(i)*play
          psc_o3=ppath_o3/upath_o3

          do ib=1,pnb_hirs
            t_cont=exp(-scoef10(ib)*ucont2-fcoef10(ib)*ucont1)

            t_h2o=t_malk(upath_h2o,psc_h2o,xh2o(ib),yh2o(ib))
            t_co2=t_malk(upath_co2,psc_co2,xco2(ib),yco2(ib))
            t_o3=t_malk(upath_o3,psc_o3,xo3(ib),yo3(ib))
            tband(ib)=t_co2*t_o3*t_h2o*t_cont


            b_ir(ib)=plnck(cwn(ib),a_hirs(ib),b_hirs(ib),tlay)
            dtrans(ib)=otrans(ib)-tband(ib)
            rad_lay(ib)=  b_ir(ib)*dtrans(ib)
            radir(ib)= radir(ib)+ rad_lay(ib)
            otrans(ib) = tband(ib)

          end do
        end do     ! end of loop over vertical levels


!
!    add in the surface contribution to the radiance, including
!    reflection
!
        
!       if(msu_flag.eq.1) then
          if(oro(icol).eq.0)then
            eps=0.65      ! ocean
          else
            eps=0.9    ! land or sea-ice
          end if
      
          do if=1,pnf_msu
            b_mw = 1.47445e-23*freq(if)**3/(exp(0.047981*freq(if)/ &
               sfctemp)-1)
            rad(if) = rad(if) + eps*trans(if)*b_mw
            refl(if)=(1-eps)*rad2(if)*trans(if)
!           britemp1(icol,if) = 0.047981*freq(if)/log(1.0 
!    $                + 1.47445e-23*freq(if)**3/rad(if))
!           britemp2(icol,if) = 0.047981*freq(if)/log(1.0 
!    $                + 1.47445e-23*freq(if)**3/rad2(if))
!           britemp3(icol,if) = 0.047981*freq(if)/log(1.0 
!    $                + 1.47445e-23*freq(if)**3/refl(if))
     
            rad3(if) = rad(if) + refl(if)
            britemp(icol,if) = 0.047981*freq(if)/log(1.0  &
                      + 1.47445e-23*freq(if)**3/rad3(if))
          end do
!       endif
      
        do ib=1,pnb_hirs
          b_ir(ib)=plnck(cwn(ib),a_hirs(ib),b_hirs(ib),sfctemp)
          radir(ib)=radir(ib)+b_ir(ib)*tband(ib)
          tb_ir(icol,ib)=btemp(cwn(ib),a_hirs(ib),b_hirs(ib),radir(ib))
        end do

      end do   ! end of loop over columns

    end subroutine hirsrtm
!
!==================================================================
!
    function btemp(wvn,a,b,chnrad)
!
!*  calculates the brightness temperature given the channel radiance.
!   uses planck function tuned to tovs frequencies
!
      real(r8) btemp

!------------------------------arguments--------------------------------
      real(r8) wvn
      real(r8) chnrad
      real(r8) a
      real(r8) b
!------------------------------local variables--------------------------
      real(r8) c1
      real(r8) c2
!
!   planck function parameters
!   note: parameters a and b are temperature correction factors
!   which are dependent on the channel and satellite.  these
!   parameters were extracted from the rttovs model.
!
      parameter( c1=1.191066e-08 )
      parameter( c2=1.438833 )

      btemp=(c2*wvn/log(c1*wvn**3/chnrad+1.)-a)/b

    end function btemp
!
!==================================================================
!
    function plnck(wvn,a,b,t)

!
! planck function
!
      real(r8) plnck

      real(r8) wvn,a,b,t,c1,c2

      parameter( c1=1.191066e-08 )
      parameter( c2=1.438833 )

      plnck=c1*(wvn)**3/(exp(c2*wvn/(a+b*t))-1.)

    end function plnck
!
!===================================================================
!
    function t_malk(u,p,x,y)

      real(r8) t_malk

      real(r8) u, p, x, y
      real(r8) p0, dnu, pi, b, bp

      parameter( p0=1013. )
      parameter( dnu=10. )

      pi=acos(-1.)
      b=4.*x**2/(pi*y*dnu)
      bp=pi*b*p/p0
      t_malk = exp(-0.5*bp*(sqrt(1.+4.*y*u/(dnu*bp))-1.))

    end function t_malk
!
!===================================================================
!
    subroutine lookup(p,t,abs)

      real(r8) p, t, abs(pnf_msu)
!
      integer if, i, j, n, m
      parameter( n=17 )
      parameter( m=16 )
      real(r8) xx(n), yy(m)
      real(r8) zz1(n,m), zz2(n,m), zz3(n,m), zz4(n,m)
      real(r8) t1, t2
      data yy/5.0,7.5,10.,25.,50.,100.,200.,300., &
              400.,500.,600.,700.,800.,900.,1000,1050./
      data xx/160.,170.,180.,190.,200.,210.,220.,230.,240.,250., &
              260.,270.,280.,290.,300.,310.,320./
      data zz1/0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000, &
               0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000, &
               0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000, &
               0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000, &
               0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000, &
               0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000, &
               0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000, &
               0.0000,0.0000,0.0002,0.0002,0.0001,0.0001,0.0001, &
               0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0000, &
               0.0000,0.0000,0.0000,0.0000,0.0000,0.0008,0.0007, &
               0.0006,0.0005,0.0004,0.0004,0.0003,0.0003,0.0003, &
               0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0001, &
               0.0001,0.0032,0.0027,0.0023,0.0020,0.0017,0.0015, &
               0.0013,0.0011,0.0010,0.0009,0.0008,0.0008,0.0007, &
               0.0006,0.0006,0.0006,0.0005,0.0130,0.0108,0.0091, &
               0.0078,0.0068,0.0059,0.0052,0.0046,0.0041,0.0037, &
               0.0033,0.0030,0.0028,0.0025,0.0023,0.0022,0.0020, &
               0.0292,0.0244,0.0206,0.0176,0.0152,0.0133,0.0117, &
               0.0103,0.0092,0.0083,0.0074,0.0067,0.0061,0.0056, &
               0.0052,0.0048,0.0045,0.0521,0.0434,0.0367,0.0314, &
               0.0271,0.0237,0.0208,0.0184,0.0164,0.0147,0.0132, &
               0.0120,0.0109,0.0100,0.0092,0.0085,0.0079,0.0818, &
               0.0681,0.0576,0.0493,0.0426,0.0371,0.0326,0.0289, &
               0.0257,0.0230,0.0207,0.0188,0.0170,0.0156,0.0143, &
               0.0132,0.0122,0.1182,0.0985,0.0833,0.0712,0.0616, &
               0.0537,0.0472,0.0418,0.0372,0.0333,0.0300,0.0271, &
               0.0246,0.0225,0.0206,0.0189,0.0175,0.1617,0.1348, &
               0.1139,0.0975,0.0842,0.0734,0.0645,0.0571,0.0508, &
               0.0455,0.0409,0.0370,0.0336,0.0306,0.0281,0.0258, &
               0.0238,0.2123,0.1770,0.1496,0.1280,0.1106,0.0965, &
               0.0848,0.0750,0.0667,0.0597,0.0537,0.0486,0.0441, &
               0.0402,0.0368,0.0338,0.0312,0.2701,0.2253,0.1905, &
               0.1630,0.1408,0.1228,0.1079,0.0954,0.0849,0.0760, &
               0.0684,0.0618,0.0560,0.0511,0.0467,0.0429,0.0396, &
               0.3353,0.2798,0.2366,0.2024,0.1749,0.1525,0.1340, &
               0.1185,0.1055,0.0944,0.0849,0.0767,0.0695,0.0634, &
               0.0579,0.0532,0.0490,0.3707,0.3094,0.2617,0.2239, &
               0.1935,0.1687,0.1482,0.1311,0.1166,0.1043,0.0938, &
               0.0848,0.0769,0.0700,0.0640,0.0588,0.0541/ 
      data zz2/0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0001, &
               0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001, &
               0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001, &
               0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001, &
               0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0002, &
               0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002, &
               0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002, &
               0.0002,0.0002,0.0010,0.0010,0.0010,0.0011,0.0011, &
               0.0012,0.0012,0.0013,0.0013,0.0014,0.0014,0.0014, &
               0.0015,0.0015,0.0015,0.0015,0.0015,0.0040,0.0038, &
               0.0038,0.0039,0.0041,0.0042,0.0044,0.0046,0.0048, &
               0.0050,0.0051,0.0053,0.0054,0.0055,0.0056,0.0056, &
               0.0057,0.0144,0.0135,0.0131,0.0130,0.0132,0.0135, &
               0.0140,0.0145,0.0150,0.0155,0.0161,0.0166,0.0170, &
               0.0174,0.0177,0.0180,0.0182,0.0525,0.0473,0.0439, &
               0.0417,0.0405,0.0399,0.0399,0.0402,0.0407,0.0414, &
               0.0421,0.0429,0.0437,0.0445,0.0452,0.0458,0.0464, &
               0.1135,0.1005,0.0913,0.0848,0.0802,0.0772,0.0752, &
               0.0740,0.0734,0.0732,0.0733,0.0736,0.0739,0.0744, &
               0.0749,0.0753,0.0757,0.1971,0.1730,0.1553,0.1422, &
               0.1326,0.1255,0.1202,0.1164,0.1137,0.1118,0.1104, &
               0.1095,0.1088,0.1084,0.1081,0.1079,0.1077,0.3029, &
               0.2645,0.2358,0.2140,0.1975,0.1848,0.1750,0.1675, &
               0.1617,0.1572,0.1537,0.1509,0.1487,0.1468,0.1453, &
               0.1440,0.1429,0.4307,0.3749,0.3325,0.3000,0.2747, &
               0.2550,0.2395,0.2272,0.2174,0.2095,0.2030,0.1978, &
               0.1934,0.1897,0.1865,0.1837,0.1812,0.5801,0.5037, &
               0.4452,0.3998,0.3642,0.3360,0.3134,0.2953,0.2805, &
               0.2684,0.2584,0.2500,0.2429,0.2368,0.2316,0.2269, &
               0.2228,0.7503,0.6505,0.5734,0.5131,0.4655,0.4273, &
               0.3966,0.3715,0.3509,0.3339,0.3196,0.3075,0.2972, &
               0.2882,0.2805,0.2736,0.2674,0.9408,0.8148,0.7168, &
               0.6396,0.5782,0.5288,0.4887,0.4557,0.4284,0.4056, &
               0.3864,0.3700,0.3560,0.3437,0.3331,0.3236,0.3151, &
               1.1504,0.9956,0.8746,0.7787,0.7021,0.6400,0.5893, &
               0.5475,0.5126,0.4834,0.4586,0.4374,0.4191,0.4032, &
               0.3892,0.3768,0.3657,1.2620,1.0919,0.9586,0.8528, &
               0.7679,0.6991,0.6427,0.5961,0.5572,0.5245,0.4967, &
               0.4728,0.4523,0.4343,0.4186,0.4046,0.3921/
      data zz3/0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0001, &
               0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001, &
               0.0001,0.0001,0.0001,0.0004,0.0004,0.0003,0.0003, &
               0.0003,0.0003,0.0003,0.0003,0.0003,0.0003,0.0003, &
               0.0003,0.0003,0.0003,0.0003,0.0002,0.0002,0.0006, &
               0.0006,0.0006,0.0006,0.0006,0.0006,0.0006,0.0006, &
               0.0006,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005, &
               0.0004,0.0004,0.0040,0.0039,0.0038,0.0038,0.0037, &
               0.0037,0.0036,0.0035,0.0034,0.0033,0.0032,0.0031, &
               0.0030,0.0029,0.0028,0.0027,0.0026,0.0154,0.0150, &
               0.0147,0.0145,0.0143,0.0141,0.0139,0.0136,0.0133, &
               0.0130,0.0126,0.0122,0.0119,0.0115,0.0111,0.0107, &
               0.0103,0.0543,0.0528,0.0518,0.0511,0.0504,0.0498, &
               0.0491,0.0484,0.0475,0.0465,0.0454,0.0443,0.0431, &
               0.0419,0.0406,0.0393,0.0380,0.1702,0.1617,0.1560, &
               0.1521,0.1492,0.1470,0.1449,0.1430,0.1409,0.1387, &
               0.1363,0.1338,0.1310,0.1281,0.1251,0.1219,0.1187, &
               0.3229,0.3005,0.2846,0.2729,0.2642,0.2573,0.2516, &
               0.2467,0.2421,0.2377,0.2334,0.2290,0.2246,0.2200, &
               0.2153,0.2104,0.2055,0.5089,0.4673,0.4363,0.4128, &
               0.3946,0.3801,0.3681,0.3579,0.3489,0.3407,0.3330, &
               0.3257,0.3185,0.3115,0.3045,0.2975,0.2906,0.7263, &
               0.6605,0.6104,0.5716,0.5410,0.5162,0.4957,0.4783, &
               0.4631,0.4496,0.4373,0.4258,0.4149,0.4046,0.3945, &
               0.3848,0.3752,0.9737,0.8791,0.8060,0.7486,0.7028, &
               0.6654,0.6344,0.6080,0.5852,0.5651,0.5469,0.5303, &
               0.5149,0.5005,0.4867,0.4735,0.4609,1.2496,1.1217, &
               1.0219,0.9429,0.8792,0.8271,0.7836,0.7467,0.7149, &
               0.6870,0.6621,0.6395,0.6188,0.5995,0.5815,0.5645, &
               0.5482,1.5516,1.3866,1.2568,1.1532,1.0693,1.0004, &
               0.9428,0.8939,0.8518,0.8150,0.7824,0.7531,0.7263, &
               0.7018,0.6789,0.6575,0.6373,1.8769,1.6714,1.5088, &
               1.3782,1.2720,1.1844,1.1111,1.0489,0.9954,0.9488, &
               0.9076,0.8707,0.8374,0.8069,0.7788,0.7527,0.7282, &
               2.2222,1.9737,1.7758,1.6162,1.4859,1.3780,1.2876, &
               1.2109,1.1450,1.0876,1.0371,0.9921,0.9516,0.9147, &
               0.8809,0.8497,0.8206,2.4013,2.1306,1.9144,1.7396, &
               1.5966,1.4781,1.3787,1.2944,1.2219,1.1588,1.1034, &
               1.0541,1.0098,0.9696,0.9328,0.8989,0.8673/ 
      data zz4/0.0027,0.0023,0.0019,0.0016,0.0014,0.0012,0.0011, &
               0.0009,0.0008,0.0007,0.0006,0.0006,0.0005,0.0005, &
               0.0004,0.0004,0.0003,0.0059,0.0050,0.0043,0.0037, &
               0.0032,0.0027,0.0024,0.0021,0.0018,0.0016,0.0014, &
               0.0013,0.0011,0.0010,0.0009,0.0008,0.0007,0.0105, &
               0.0089,0.0076,0.0065,0.0056,0.0049,0.0042,0.0037, &
               0.0033,0.0029,0.0025,0.0023,0.0020,0.0018,0.0016, &
               0.0015,0.0013,0.0649,0.0550,0.0468,0.0402,0.0347, &
               0.0301,0.0262,0.0230,0.0202,0.0178,0.0158,0.0141, &
               0.0126,0.0113,0.0101,0.0091,0.0082,0.2465,0.2097, &
               0.1794,0.1544,0.1336,0.1162,0.1015,0.0891,0.0785, &
               0.0695,0.0617,0.0550,0.0491,0.0441,0.0396,0.0358, &
               0.0323,0.8274,0.7127,0.6168,0.5364,0.4684,0.4108, &
               0.3616,0.3196,0.2834,0.2522,0.2251,0.2015,0.1810, &
               0.1630,0.1471,0.1332,0.1208,2.1322,1.8745,1.6542, &
               1.4649,1.3016,1.1601,1.0371,0.9298,0.8358,0.7532, &
               0.6805,0.6163,0.5593,0.5088,0.4638,0.4236,0.3876, &
               3.2645,2.8943,2.5762,2.3012,2.0625,1.8542,1.6717, &
               1.5114,1.3698,1.2446,1.1334,1.0344,0.9460,0.8668, &
               0.7958,0.7319,0.6743,4.2737,3.8017,3.3956,3.0442, &
               2.7385,2.4714,2.2370,2.0305,1.8478,1.6858,1.5415, &
               1.4126,1.2972,1.1936,1.1004,1.0162,0.9400,5.2092, &
               4.6423,4.1542,3.7314,3.3633,3.0414,2.7585,2.5090, &
               2.2882,2.0920,1.9171,1.7607,1.6205,1.4945,1.3809, &
               1.2782,1.1851,6.0871,5.4326,4.8683,4.3790,3.9525, &
               3.5791,3.2507,2.9608,2.7039,2.4754,2.2716,2.0892, &
               1.9256,1.7782,1.6453,1.5251,1.4161,6.9131,6.1782, &
               5.5437,4.9928,4.5120,4.0906,3.7196,3.3917,3.1008, &
               2.8419,2.6106,2.4035,2.2175,2.0500,1.8987,1.7617, &
               1.6374,7.6903,6.8822,6.1833,5.5755,5.0445,4.5784, &
               4.1676,3.8041,3.4813,3.1937,2.9366,2.7061,2.4989, &
               2.3121,2.1432,1.9903,1.8514,8.4213,7.5468,6.7890, &
               6.1290,5.5515,5.0440,4.5961,4.1994,3.8467,3.5321, &
               3.2507,2.9981,2.7708,2.5657,2.3801,2.2119,2.0590, &
               9.1086,8.1740,7.3627,6.6548,6.0345,5.4886,5.0062, &
               4.5785,4.1978,3.8579,3.5535,3.2801,3.0338,2.8114, &
               2.6100,2.4272,2.2610,9.4365,8.4743,7.6380,6.9077, &
               6.2673,5.7033,5.2047,4.7622,4.3683,4.0163,3.7009, &
               3.4175,3.1621,2.9314,2.7224,2.5326,2.3600/ 
                    
      if(p.le.5) then
        do if = 1, pnf_msu
          abs(if)=0.0
        end do
        return
      endif
      
      call locate(xx,n,t,i)
      call locate(yy,m,p,j)
            
      t1=(t-xx(i))/(xx(i+1)-xx(i))
      t2=(p-yy(j))/(yy(j+1)-yy(j))
      abs(1)=(1-t1)*(1-t2)*zz1(i,j)+t1*(1-t2)*zz1(i+1,j)+ &
                      t1*t2*zz1(i+1,j+1)+(1-t1)*t2*zz1(i,j+1)
      abs(2)=(1-t1)*(1-t2)*zz2(i,j)+t1*(1-t2)*zz2(i+1,j)+ &
                      t1*t2*zz2(i+1,j+1)+(1-t1)*t2*zz2(i,j+1)
      abs(3)=(1-t1)*(1-t2)*zz3(i,j)+t1*(1-t2)*zz3(i+1,j)+ &
                      t1*t2*zz3(i+1,j+1)+(1-t1)*t2*zz3(i,j+1)
      abs(4)=(1-t1)*(1-t2)*zz4(i,j)+t1*(1-t2)*zz4(i+1,j)+ &
                      t1*t2*zz4(i+1,j+1)+(1-t1)*t2*zz4(i,j+1)
           
    end subroutine lookup
!
!===================================================================
!
    subroutine locate(xx,n,x,j)

      integer n, j
      real(r8) xx(n), x
!
      integer jl, ju, jm
!
      jl=0
      ju=n+1
      do while (ju-jl.gt.1)
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      end do
      j=jl

    end subroutine locate

    subroutine hirsbt_init ()

! Initialize values used for the HIRS brightness temperature calculation.
!
    use time_manager, only: get_step_size

!------------------------------Local variables--------------------------
    integer :: dtime      ! integer timestep size

!*******************************************************************************************
!     Set constants to values used in the model
!     Leave as they are for now.  Later if this gets committed back to the main 
!     development trunk we should use the values used in the rest of the radiation code.
! mji
! This step is done in hirsrtm.f90 to use values consistent with the radiation code.
!
!     GRAV = GRAVIT
!     R = RAIR
!     RCO2 = CO2VMR*RMWCO2
!
!*******************************************************************************************

    msuname(1)  = 'MSU_1   '
    msuname(2)  = 'MSU_2   '
    msuname(3)  = 'MSU_3   '
    msuname(4)  = 'MSU_4   '
    hirsname(1) = 'HIRS_2  '
    hirsname(2) = 'HIRS_4  '
    hirsname(3) = 'HIRS_6  '
    hirsname(4) = 'HIRS_8  '
    hirsname(5) = 'HIRS_10 '
    hirsname(6) = 'HIRS_11 '
    hirsname(7) = 'HIRS_12 '

    dtime  = get_step_size()

! These should be namelist variables; 
! Set flag to do HIRS brightness temperature calculation 
    dohirs = .true.
! Set frequency of HIRS calculation
! ihirsfq is in timesteps if positive, or hours if negative; 6 hours is recommended
    ihirsfq = -6

! Convert ihirsfq from hours to timesteps if necessary
    if (ihirsfq < 0) ihirsfq = nint((-ihirsfq*3600._r8)/dtime)

    end subroutine hirsbt_init

    end module hirsbt

