#include <def-undef.h>
!050208[Mardi Gras]	Note that in early march 2004 did a test with the zero shear 
!		parameterization and it seemed okay. 
!		Change a few single precision constants to double for consistency today.
!040422Zi1bj	While at MIT, futher correct write out of eplatidependmin(Found mis-spelling.)
!	        and explicitly zero out the undefined "al_back" in the ifepson2>0 case.
!040217Zi1b	Introduce a minimum foreground shear, "s2min". This was not necessary for
!	        the 3X3 25 layer NCOM with strong restoring in the arctic regions, but may
!		be needed for more realistic, higher resolution OGCMs. It is to avoid 
!		problems that may occur with the (l^2 S) turbulence dimensionalization  
!	        when the shear becomes extremely close to zero and (l^2 S) becomes an
!	        innaccurate dimensionalization. This is in lieu of using the zero shear
!	        parameterization, because "ifzeroshear.EQ.TRUE" has not been properly 
!		checked out yet. Patch of turb_2gi1a:
!031217-30Zi1a' Correct write out of floor on latitude dependent factor, eplatidependmin.
!030803Zi1a *CORRECT LATITUDE DEPENDENT MIXING ("ifdeeplat=1") OPTION.*
!	   *TEST FOR WHEN TO REVERT TO THE MINIMUM ON "eplatidepent"(\equiv \theta,N)"*
!	   *BECAUSE THE GREGG ET AL. FORMULA CEASED TO GIVE A REAL NUMBER*
!	   *HAD BEEN BOTCHED BY USING "N/f<1" WHEN IN FACT SHOULD BE "N/|f|<1"*
!	   *BECAUSE GREGG ET AL. USE "f" WHEN THEY MEAN THE ABSOLUTE VALUE OF*
!	   *THE CORIOLIS PARAMETER WHICH OTHERS CALL `f'.*
!	   [See NBp.030803-1,5 Vol.XIX .] Correction of turb_2gia:
!030717-22Z1a INTRODUCE OPTIONS "icondear=-1" AND "icondear=+1" FOR NO USE OF DEARDORFF
!	   LENGTHSCALE LIMITATION AND USE OF DEARDORFF LENGTHSCALE LIMITATION BUT
!	   WITHOUT ALTERING `\tau N' TO TRY TO DEAL WITH CONSISTENCY ISSUES BETWEEN
!	   PRODUCTION=DISSIPATION `\tau N' USED FOR S_{M,H,S} AND DEARDORFF LIMITED ONE.
!	   [See NBp.030714-2,3;15-2;17-6.] Modification of turb_2gi:
!030429-0504Z1  *INTRODUCE AN OPTION "ifdeeplat=1" FOR USE OF A LATITUDE AND STRATIFICATION*
!Z1	   *DEPENDENT FACTOR MULTIPLYING THE `epsilon/N^2' DEEP MIXING DIMENSIONALIZATION*
!Z1	   *BASED ON THE FORMULA (CITED AS FROM HENYEY ET. AL, JGR VOL.91 8487-8495,1986)*
!Z1	   *IN GREGG ET AL. NATURE VOL.422 513-515,2003, WHERE IT IS SHOWN CONFIRMED BY*
!Z1	   *OBSERVATIONS FOR LOWER LATITUDES EXCEPT FOR BEING LOW AT THE VERY EQUATOR.*
!Z1	   *I PLACE A MINIMUM ON THE FACTOR.*	
!0502Z1	   **Note that Gregg et. al.'s formula:**
!Z1	    "L(\theta,N) = (f cosh^{-1} (N/f))/(f_30^o cosh^{-1} (N_0/f_30^o)"**
!Z1	   **is only defined as a real number when N > f, since arccosh**
!Z1	   **[which is what they mean by cosh^{-1}] can only be defined as a real**
!Z1	   **for arguments of at least 1. At 1 it is zero.**
!Z1	   I DECIDE TO TAKE "L(\theta,N)" TO BE **ZERO** FOR (N/f < 1).
!Z1	   THIS CORRESPONDS TO SETTING A FLOOR OF 1 on (N/f). 
!Z1	   FOR FOREGROUND MIXING AT DEPTH DETACHED FROM THE MIXED-LAYER
!Z1	   I REVERT TO THE "DEEP" LENGTHSCALE, WHICH USES DENSITY GRADIENTS,
!Z1	   IN CASE (N/f)<1 TO TRY NOT TO MAKE DEEP ARCTIC&SUBARCTIC MIXING TOO SMALL.
!Z1	   Uses function eplatidepend_ dated 030430. Extension of turb_2g:
!030424Z   ****ADD A NEW INPUT ARGUMENT (PLACED AFTER SURFACE FORCING) FOR CORIOLIS****
!Z	   ****PARAMETER, USED FOR ROTATION'S EFFECT ON TURBULENCE[See NBp.030424-8.].****
!Z	   **INTRODUCE A CHOICE "ilomega" FOR ROTATIONAL EFFECT ON LENGTHSCALE**
!Z         **[See NBp.030424-13,25-2to4.]**
!Z	   Allow backward compatibility by using zero shear approximation *only when 
!Z	   the new logical parameter "ifzeroshear" is ".TRUE."* [See NBp.030424-10&12].
!Z	   Extension of turb_2f.
!030401-07Y   ****ADD A NEW INPUT ARGUMENT (PLACED AFTER THE BACKGROUND DIFFUSIVITIES)****
!Y	   ****FOR THE SQUARE OF THE BRUNT VAISALA FREQUENCY TO BE USED FOR THE CASE****
!Y         ****WHERE RICHARDSON NUMBER IS MORE NEGATIVE THAN THE TABLE MINIMUM.****
!Y	   ****THE ZERO SHEAR APPROXIMATION USED HERE ORIGINATED IN THE HYCOM VERSION.****
!03-04Y        [See NBp.030331-1to030401-4,030403-13 and 030404-01to2.]
!03Y	   Upgrade to incorporate improvements inspired by
!Y	   0 gradient problems found by Wallcraft&Halliwell in hycom turb_2cc1 of turb_2e.
!03Y	   ADDUCE "ttot"=`\tau_{\theta} \over \tau' TO THE TIMESCALES IN THE COMMON BLOCK
!Y   "/bb0/",REQUIRED FOR USE IN "quad_sal_pureconv" ROUTINE [See NBp030403-7to10,13to14].
!04Y	   Add N^2 to fort.92 and add an output of the zeroshear table.[See NBp.030404-8]
!07Y	   Correct error in version sent to Halliwell of recursive use of -mt0,-mt0+1 table
!Y	   values for zeroshear case before they were calculated by calculating out from
!Y	   zero in two pieces as was done for Ri and Ri_d. [See NBp.030407-11.]
!030124,28 oursal2_2 corrected.
!030123Y   ******INTRODUCE SUB-OPTION "ifast=1" IN "ifexpabs=1" CASE TO SPEED UP TABLE*****
!	   ******CALCULATIONS IN THE EXPONENTIAL ABSOLUTE VALUE REGIME BY ESTIMATING****** 
!	   ******THE INDICES WHICH CORRESPOND TO THE RICHARDSON NUMBERS [See NBp030122-6,7,10,11,23-2,4,5 .]******
!030116X1a George Halliwell, whom I'd sent my commented Natassia's version of the module, 
!	   based I think on my 000316-30 stage, reported an error of `isalback' for "ifsalback".
!	   It occurred in the line `IF(isalback.EQ.6) ibg=0', which I now change to
!	   "IF(ifsalback.EQ.6) ibg=0". Note that the ifsalback=6 option had never worked, 
!	   I don't know if this is why. Have not tested it in a long time. Still untested now.
!021210X1 Change the name of the submodule that calculates dimensionless turbulence 
!	  for the CC salinity model to "oursal2_2" which has an option for B1=y(0,0)^{3/4}.
!020925,26X Correct erroneous calculation of epsy (dissipation times (\tau S)^2). Remove bogus one half.
!X	 See NBp020925-2{extension on cover}. Corrects error on NBp020911-9 propogated NBp020912-12.
!020912-24X ******INTRODUCE SURFACE FLUXES AS ARGUMENTS TO THE TURBULENCE ROUTINE.******
!X	 ****RECEIVES FRICTION VELOCITY "ustar_" AND BUOYANCY FLUXES "buoytur" AND "buoysol" INTO OCEAN.****
!X	 **"isurfuse" CHOOSES ABOUT USING THEM. SEE "NBp020912-5to8" IN VOLUME XV.**
!X	CALCULATES THE COLUMN ARRAY "epsy" WHICH IS "\epsilon \times (\tau S)^2". 
!X	FOR "isurfuse=1" SET "K_X_foreground = (epsy/(2 S^2)) S_X". SEE "NBp020912-12" in VOLUME XV.
!020404 Correct writeout of epsilon.
!020213-22D  MODIFICATION TO INCLUDE AN OPTION 
!        FOR A SIMPLE ENHANCED BOTTOM MIXING AS A FUNCTION
!	 OF HEIGHT ABOVE THE BOTTOM OF OCEAN TURBULENCE II PAPER PROGRAM "turb_2cc1".
!000316-30  INTRODUCE AN OPTION TO USE FOR THE BACKGROUND MODEL AT POSITIVE "Ri"
!	   AT "theta_r" WHERE TURBULENCE EXISTS AT "ra_r=>infinity" 
!	   THE MINIMUM OF THE CALCULATED BACKGROUND "ra_r" 
!          AND THE FOREGROUND "ra_r".
!000309-15 INTRODUCE AN OPTION TO INTERPOLATE FROM A 1D "theta_r" TABLE FOR BG. 
!	   INTRODUCE A SEPARATE LENGTHSCALE, "L_DEEP", 
!          AND USE IT FOR "ifepson2=2" CASE WHEN "Ri<0".
!	   Correct errors in *INTERPOLATION* and table making of background.
!000229-0303  EXTEND WITH OPTION "ifepson2=2" to dimensionalize with "\epsilon/(N^2)"
!	 for the FOREGROUND CASE ALSO WHEN BENEATH A BACKGROUND-ONLY LEVEL.
!000215	 INTRODUCE OPTION "ifepson2=1" TO DIMENSIONALIZE 
!	 BACKGROUND DIFFUSIVITIES WITH "\epsilon/(N^2)" INSTEAD OF "l^2 S".
!000128-30,0201  Allow the timescale ratios to be adjustable
!	  by replacing oursal2 with oursal2_1 which leaves 
!	  the timescale ratio setting entirely to the smshsc$ program and
!	  by replacing smshsc2 with smshsc_a3 which
!	  CALCULATES THE "p's" FROM THE TIMESCALE RATIOS
!	  in setting the model constants.
!000107	  Make the tiny revision of adding a writeout of nmodel.
!991107-09  Version with option for constant Ri background turbulnce in the
!	  Temperature only model based on turb_2ca:
!990928   Version with option to use exponential regularization in Dubovikov's 
!	  salinity model based on turb_2c:
!990702   ************   New Module which incorporates option for Dubovikov's
!	              990621 2 point closure with salt.   ************
!         ************   BASED ON turb_2bi2o990513   ************
!990513-0629***CORRECT ERROR IN TURBULENCE MODEL CONSTANTS***
!990629   Update "Turbulence calculated by" model version output. Had forgot to.
!         Old routine "smshsc" had used an inconsistent set of model constants
!	  with p10 corrected from the erroneous value first given to me,
!	  but the a's and d's not recalculated using the new p10 as they should
!         have been. Since the p's appear both directly and indirectly 
!         through the a's and d's in the calculation of the S's, the old module
!         using "smshsc" had an unseemly mixture of old and new p10's.
!	  This version has the new routine "smshsc2" 
!	  which *CALCULATES* the a's and d's from the p's inside the routine
!	  to *ensure* consistency.
!990412-15 COVER SPECIAL CASES IMPROPERLY DEALT WITH IN bs 2 and bs>4.
!	When *Ri_T=0* use for bs2 and bs>4 an angle, \theta_r of \pi/2, 
!	if *Ri_C \ne 0*. When *Ri_T=Ri_C=0* set \theta_r = \pi/4 for bs2.
!	When Ri=0 set background diffusivity to zero like Ri<0 for bs>4.
!990303-23 ADD an Option for a Dubovikov-type background model in which
!	the internal wave Richardson number is a function of Ri_d, obtained by
!	taking ra_r such that S_M/(S l /q)[ra_r,\theta_r] = S_M/(S l/q)[ra_r=0].
!990226-0303 Correction to use double floats in integer multiplications of "dri"
!	AND ADD Option for a Dubovikov-type background model in which
!	the internal wave Richardson number is a function of Ri_d, 
!	obtained by taking ra_r = fraction * ra_r_critical(\theta_r),
!	where {ra_r}^2 \equiv ({Ri_T}^2 + {Ri_C}^2) 
!	Ri_T \equiv (cos(\theta_r)*ra_r), Ri_C \equiv (sin(\theta_r)*ra_r),
!	Ri \equiv (Ri_T + Ri_C), Ri_d \equiv (Ri_T - Ri_C) and
!	Ri_critical(Ri_d) \equiv (cos(\theta_r)*ra_r_critical(\theta_r)) 
!	with Ri_d = sin(\theta_r)*ra_r_critical(\theta_r).
!990204-08 Option to take background MOMENTUM, HEAT and SALT diffusivities
!	as being given by our usual turbulence model with the *inputs*
!	altered by having Shear squared become 
!	the Brunt Vaisala frequency squared divided by 
!	a constant internal wave Richardson number,
!	and the limiting lengthscale 
!	become a constant internal wave lengthscale.
!	This is to extend Mikhail Dubovikov's internal wave diffusivity
!	model, in which Ri is a constant in the deeper ocean, so that
!	S ~ N, to cover the upper ocean as well for continuity.
!	The constant internal wave Richardson number is an adjustable
!	parameter, but should be close to unity. 
!990201-04 Option to take background MOMENTUM, HEAT and SALT diffusivities
!	as being given by a small residual background plus
!	our usual turbulence model with the *inputs*
!	altered by having values of Shear squared and limiting lengthscale
!	taken as constants based on Canuto's reading of
!	"A Composite Spectrum of Vertical Shear in the Upper Ocean",
!	Gargett et. al. , JPO, Vol. 11, September 1981. 
!	This is to represent the nature of turbulence excited by shear
!	generated by internal waves presumably absent from the ocean model.
!	Brunt Vaisala frequency remains the literal ocean model value.
!981215-990119 Writeouts to results of some turbulence info added AND
!	Option to make background ratios of heat and salt diffusivities
!	when S_H and S_S have vanished the ratio that they *would have*
!	at the ocean point's Ri_C/Ri_T if Ri were just less than critical.
!981125 Option to make background ratios of heat and salt diffusivities
!	follow the ratios of the dimensionless turbulence functions S_H and S_S
!	using the ratio at the nearest point above when these are zero.
!981006-1110 FIRST ATTEMPT TO IMPLEMENT SALINITY MODEL WITH TABLE 
!980912	  Cheng model option to use the constants which Ye Cheng has been using
!	  in the atmosphere and which he believes give a better match to
!	  experimental data, although they also give a much smaller Ri_maximum.
!980821-24         Dubovikov model numerical solver error fixed following
!         Cheng's offline routine kirk:/data2/acamh/NCAR/Offline/ch_smshplot.f.
!980716-17 AH alteration to correct MLD calculation error, and
!	   replace Dubovikov mksmsh with the subroutine of the same name
!	   from /data1/acamh/TURB1/CHENG/OCode/Canuto/3D/Cheng/
!	   mike_12.f_mod980528 and change rimax to the value from
!	   mike_12.f_980528 in the same directory, and
!	   introduce an option for outputs to check turbulence model
!	   as proposed by Ye Cheng.
!980501-27 Modification by Armando Howard of Ye Cheng's turb_2 modular subroutine
!       for GISS turbulent vertical mixing models.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     beginning of turbulence models
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!030424Z****INTRODUCE THE CORIOLIS PARAMETER "2 \Omega sin \phi" FOR ROTATION MODELS****
!030401Y****INTRODUCE THE SQUARE OF THE BRUNT VAISALA FREQUENCY FOR ZERO SHEAR MODEL****
!020912X******INTRODUCE SURFACE FLUX INPUTS******
      subroutine turb_ini
      use precision_mod
      use param_mod, only: mytid
      use pconst_mod, only: dfricmx
      use canuto_mod

      implicit real(r8) (a-h,o-z)
!YU Jan. 20th, 2011
      real(r8)::  ttot,tcot,tctot,tptot,tpcot,tpvot
!X       and logical array lifepsy for points in column where use epsilon y dimensionalization.
      LOGICAL lifupper
!030424-25Z **Introduce an integer parameter for the effect of rotation on the lengthscale.**
!25Z        Introduce the complex variable zlomega for \sqrt{-B*/f^3} for diagnosis.
      COMPLEX*16 zlomega
!25Z        amldminlom is the minimum MLD for use of the lengthscale \sqrt{-B*/f^3}.
      DIMENSION aldeep(nbig)
!*****C
      ifirst = 0
!000107	Add nmodel to writeout. nmodel=1,2 for 1,2 pt. closure.
       if (mytid.eq.0) then
	WRITE(*,*) " "
	WRITE(*,*) " "
	WRITE(*,*) "************************************"
	WRITE(*,*) "************************************"
	WRITE(*,*) &
       "Turbulence calculated by turb_2gi1b 040217 version."
	WRITE(*,*) "************************************"
	WRITE(*,*) "nmodel=",nmodel
	WRITE(*,*) "************************************"
	WRITE(*,*) "************************************"
	WRITE(*,*) " "
	WRITE(*,*) " "
       endif
!*****C
!030401Y Take zero shear unstable case *foreground* calculation from my code for HYCOM
!Y	 in inigiss_fixed2.fs0 . [See NBp.030401-2to3.] 
!Y	 When the Richardson number is more negative than the most negative table value
!Y	 it is more accurate to use the zero shear approximation derived from Canuto's
!Y	 analytic pure convection formula [See NBp.030326-3to9.].
!030324-27AH Amend to make 1D table vs. (N_d^2/N^2) for zero shear unstable case.
!            Choose (N_d^2/N^2) table values to be the same as Ri_d table values.

!980501 Maximum diffusivity, fricmx, and surface minimum, wndmix, set externally
!       fricmx=80.
        visc_cbu_limit=dfricmx*1.0d4
        diff_cbt_limit=dfricmx*1.0d4
!       wndmix=10.
!990202 *SET THE VALUE OF THE KOLMOGOROV CONSTANT, K_0, HERE NAMED "ako".*
	ako = 1.6D0
!*****C
!991107 START OF SALINITY MODEL BACKGROUND LENGTHSCALE CALCULATION SECTION.
	IF(ifsali.EQ.1) THEN
!990202-26 Calculate constant lengthscale for the background for ifsalback=3,4,5
!  	\Delta_0 = {B_1 pi \over (3 Ko)^{3/2}} l_0
!	l_0 = {(3 Ko)^{3/2} \over B_1 pi} \Delta_0
!	"back_l_0" is the constant background l_0 in centimeters.
          IF(nmodel.EQ.1) THEN  
!981104	Need to pass back the value of B_1 from oursal2 for use here. 
!021210X1 Call version of submodule oursal2 which has an option to use B1={\tau S}^{3/2} (Ri=0,Rid=0).
	      CALL OURSAL2_2A(b1,0.D0,0.D0,slq2b_00, &
                     smb_00,shb_00,ssb_00, &
                     c_y0,c_y00,0,0,ii,jj,ttot,tcot,tctot,tptot,tpcot,tpvot)
	  ELSE IF(nmodel.EQ.2) THEN
!981104-990702	Need to pass back the value of B_1 from mikesal2 for use here. 
!990928 ***Replace mikesal2 with mikesal2a to include exponential regularization option.***
	    CALL MIKESAL2A(b1,0.D0,0.D0,slq2b_00, &
                     smb_00,shb_00,ssb_00, &
                     c_y0,c_y00,0,0,ii,jj)
	  END IF
      
	back_l_0 = (((3.D0*ako)**(3.D0/2.D0))/(b1*pi))*back_del_0
	IF(ifsalback.EQ.3) THEN
       if (mytid.eq.0) then
	  WRITE(*,*) " "
	  WRITE(*,*) "************************************"
	  WRITE(*,*) "Internal wave constants for background."
!990303 Add write-out of residual constant backgrounds.   
	  WRITE(*,*) "Residual Constant Background Diffusivities:"
	  WRITE(*,*) "K_M /(cm^2 sec^{-1})",v_back0
	  WRITE(*,*) "K_H /(cm^2 sec^{-1})",t_back0
	  WRITE(*,*) "K_S /(cm^2 sec^{-1})",s_back0
	  WRITE(*,*) "."
	  WRITE(*,*) " "
!*****C
	  WRITE(*,*) "Shear^2/(sec^{-2}) =",back_s2
	  WRITE(*,*) "Lengthscale, del_0/(cm) =",back_del_0
	  WRITE(*,*) "Lengthscale, l_0/(cm) =",back_l_0
	  WRITE(*,*) '"adjust_gargett="',adjust_gargett
	  WRITE(*,*) "************************************"
	  WRITE(*,*) " "
        endif
	ELSE IF(ifsalback.GE.4) THEN
       if (mytid.eq.0) then
	  WRITE(*,*) " "
	  WRITE(*,*) "************************************"
	  WRITE(*,*) "Dubovikov Internal wave constants for background."
       endif
          IF(ifsalback.EQ.4) THEN	
       if (mytid.eq.0) then
	    WRITE(*,*) "Internal wave Richardson number=",ri_internal
       endif
	  ELSE IF(ifsalback.EQ.5) THEN
!990301
       if (mytid.eq.0) then
	    WRITE(*,*) "Ratio of Background to Critical ra_r"// &
                 " [\\equiv ({Ri_T}^2 + {Ri_C}^2)^(1/2)]",backfrac
       endif
	  ELSE IF(ifsalback.EQ.6) THEN
!990303-04
       if (mytid.eq.0) then
	    WRITE(*,*) "Ratio of Background dimensionless K_M ="// &
                 " S_M/(S l/q) to its value at Ri_T=Ri_C=0", &
                 backfact
       endif
	  END IF
       if (mytid.eq.0) then
	  WRITE(*,*) "Lengthscale, del_0/(cm) =",back_del_0
	  WRITE(*,*) "Lengthscale, l_0/(cm) =",back_l_0
	  WRITE(*,*) "************************************"
	  WRITE(*,*) " "
       endif
	END IF
!*****C
	GO TO 10
!981015-990702 Calculate constants inside subroutine oursal2 
!              or mikesal2 in sali-temp model case.
      END IF
!*****C END OF SALINITY MODEL BACKGROUND LENGTHSCALE CALCULATION SECTION.
      if(nmodel.eq.1) then
!980912 Choose which set of model constants to use for Cheng model.
        IF(ifchengcon.EQ.0) THEN
          call ccoeff(b1,rimax,g_tur,d_tur,s_tur)
        ELSE IF(ifchengcon.EQ.1) THEN
          call ccoeff1(b1,rimax,g_tur,d_tur,s_tur)
        END IF
      else
          call mcoeff(b1,rimax,ttot,tcot,tctot,tptot,tpcot,tpvot)
      endif
!991107 Calculate constant lengthscale for the background for ifback > 2
!  	\Delta_0 = {B_1 pi \over (3 Ko)^{3/2}} l_0
!	l_0 = {(3 Ko)^{3/2} \over B_1 pi} \Delta_0
!	"back_l_0" is the constant background l_0 in centimeters.
      back_l_0 = (((3.D0*ako)**(3.D0/2.D0))/(b1*pi))*back_del_0
!*****C
!991107-08C Temperature=Salinity Diffusivity Models Writeouts for background
      IF(ifback.GE.4) THEN
       if (mytid.eq.0) then
        WRITE(*,*) " "
        WRITE(*,*) "************************************"
        WRITE(*,*) "Dubovikov Internal wave constants for background."
       endif
        IF(ifback.EQ.4) THEN	
       if (mytid.eq.0) then
          WRITE(*,*) "Internal wave Richardson number=",ri_internal
       endif
        ELSE IF(ifback.EQ.5) THEN
       if (mytid.eq.0) then
          WRITE(*,*) "Ratio of Background to Critical Ri = ",backfrac
       endif
        END IF
       if (mytid.eq.0) then
        WRITE(*,*) "************************************"
        WRITE(*,*) " "
       endif
      END IF
!******C
!020916X Switch for use of surface fluxes to dimensionalize turbulence model.
       if (mytid.eq.0) then
        WRITE(*,*) "************************************"
	WRITE(*,*) "isurfuse=",isurfuse
        WRITE(*,*) "************************************"
       endif
!******X
!       building the look-up tables of slq2,sm and sh vs. ri
        dri=(rimax-ri0)/float(ntbl-1)
        do k=1,ntbl
          ria(k)=ri0+DFLOAT(k-1)*dri
          if(k.eq.ntbl) ria(k)=rimax
          if(nmodel.eq.1) then
            call ourl2(b1,ria(k),slq2a(k),sma(k),sha(k),g_tur,d_tur,s_tur)
          else
            call mikel2(b1,ria(k),slq2a(k),sma(k),sha(k),ri)
          endif
        enddo
!       end of building look-up tables
       if (mytid.eq.0) then
        write(*,*) "nmodel=",nmodel," the table's ntbl=",ntbl
!991107
	WRITE(*,*) "Temperature=Salinity diffusivity model"
	WRITE(*,*) "rimax  =",rimax
	WRITE(*,*) "ifback =",ifback
!000215
!030721Z1a Write switch for Deardorff treatment.
	WRITE(*,*) " "
	WRITE(*,*) "icondear=",icondear
	IF(icondear.EQ.-1) THEN 
          WRITE(*,*) "Do *not* use Deardorff lengthscale modification." 
        ELSE IF(icondear.EQ.0) THEN
	  WRITE(*,*) "Ye Cheng's old Deardorff:"// &
                     "modify l and \tau N leaving S_X unmodified."
	ELSE IF(icondear.EQ.1) THEN
	  WRITE(*,*) "Ye Cheng's new Deardorff:"// &
                     "modify l but leave *both* \tau N and S_X unmodified."
	END IF
	WRITE(*,*) " "
!*****CZ1a
	WRITE(*,*) "ifepson2=",ifepson2
        IF(ifepson2.EQ.2) WRITE(*,*) & 
         "epsilon/(N^2) used even for strong mixing beneath weak mixing" 
!030429-0502Z1 Write switch for latitude dependence of background mixing.
	WRITE(*,*) "ifdeeplat=",ifdeeplat
	IF(ifdeeplat.GT.0) THEN
          WRITE(*,*)  "Use latitude dependence of interior ocean value of \\epsilon/N^2"
!Zi1a' Floor on latitude dependent factor in background mixing
	  WRITE(*,*) "eplatidependmin=",eplatidependmin
	END IF
!020404D,030429Z1 Must write out the parameter epson2__ for reference \epsilon/N^2 
!Z1	since epson2_ and epson2 are no longer constants.
	IF(ifepson2.GT.0) WRITE(*,*) "epson2__=",epson2__
	WRITE(*,*) " "
      endif
!*****C
	GO TO 17
   10   CONTINUE
!981019-22 Set step-size for *both* dimensions of 2D table here.
        IF(ifsali.EQ.1) dri = -ri0/DFLOAT(mt0)
        write(*,*) "dir",  dri, ri0, mt0
!*****C
!** BUILD SALINITY MODEL TABLES VS. "Ri = Ri_T + Ri_C" AND "Ri_d = Ri_T - Ri_C".
!981027 **Use separate loops for calculation of independent table variables.**
	DO iridsign=0,1
	iridstep=(-1)**iridsign 
	DO irid= 0,mt*iridstep,iridstep
!981019 Set Ri_d table values. (See NBP59,63=p#A27,30.)
           IF(ABS(irid).LE.mt0) THEN
             ridb(irid) = DFLOAT(irid)*dri
           ELSE
             mt0s = mt0*iridstep
             mtm1s = (mt0-1)*iridstep
!000320 1st day of spring introduction of exponential absolute val table option.
	     IF(ifexpabstable.EQ.0) THEN
               idifs = (ABS(irid)-mt0)*iridstep
               ridb(irid) = ridb(mt0s)*((ridb(mt0s)/ &
                      ridb(mtm1s))**(idifs**2))
	     ELSE IF(ifexpabstable.EQ.1) THEN
               idif = ABS(irid)-mt0
               ridb(irid) = ridb(mt0s)*((ridb(mt0s)/ &
                      ridb(mtm1s))**(idif))
	     END IF
           END IF
!*****C
	END DO 
	END DO
	DO irisign=0,1
	iristep=(-1)**irisign 
	DO iri= 0,mt*iristep,iristep
!981019 Set Ri table values. (See NBP59,63=p#A27,30.)
           IF(ABS(iri).LE.mt0) THEN
             rib(iri) = DFLOAT(iri)*dri
           ELSE
             mt0s = mt0*iristep
             mtm1s = (mt0-1)*iristep
!000320 1st day of spring introduction of exponential absolute val table option.
	     IF(ifexpabstable.EQ.0) THEN
               idifs = (ABS(iri)-mt0)*iristep
               rib(iri) = rib(mt0s)*((rib(mt0s)/ &
                      rib(mtm1s))**(idifs**2))
	     ELSE IF(ifexpabstable.EQ.1) THEN
               idif = ABS(iri)-mt0
               rib(iri) = rib(mt0s)*((rib(mt0s)/ &
                      rib(mtm1s))**(idif))
	     END IF
           END IF
!*****C
	END DO 
	END DO
!*****C
!011107yXI ***If using interp2d_expabs introduce ratio between adjacent Richardson numbers in nonlinear part of table.***
!030327 rnd2on2 is the ratio between adjacent N_d^2/N^2 in the nonlinear part of the zero shear 1D table.
        IF(ifastexpabs.EQ.1) THEN
          rri = rib(mt0)/rib(mt0-1)
!030424 Only calculate rnd2on2 when zero shear parameterization is enabled.
          IF(ifzeroshear) rnd2on2 = rri 	
        END IF
!yXI
	DO iridsign=0,1
	iridstep=(-1)**iridsign 
	DO irid= 0,mt*iridstep,iridstep
	   DO irisign=0,1
	   iristep=(-1)**irisign
	   DO iri= 0,mt*iristep,iristep
              IF(nmodel.EQ.1) THEN
!981104	Need to pass back the value of B_1 from oursal2 for use here. 
!021210X1 Call version of submodule oursal2 which has an option to use B1={\tau S}^{3/2} (0,0).
	        CALL OURSAL2_2A(b1,rib(iri),ridb(irid),slq2b(iri,irid), &
                       smb(iri,irid),shb(iri,irid),ssb(iri,irid), &
                       c_y0,c_y00,iri,irid,ii,jj,ttot,tcot,tctot,tptot,tpcot,tpvot)
              ELSE IF(nmodel.EQ.2) THEN
!981104-990702	Need to pass back the value of B_1 from mikesal2 for use here. 
!990928 *** Replace mikesal2 by mikesal2a ***
	        CALL MIKESAL2A(b1,rib(iri),ridb(irid),slq2b(iri,irid), &
                       smb(iri,irid),shb(iri,irid),ssb(iri,irid), &
                       c_y0,c_y00,iri,irid)
	      END IF
	      IF(slq2b(iri,irid).LT.0) THEN
	        irimax(irid) = iri - 1 
	        GO TO 15
	      END IF
	   END DO
   15      CONTINUE
	   END DO
	END DO	
	END DO
!**
!030421Z Only calculate the zero shear table when the zero shear option is enabled.
	IF(ifzeroshear) THEN
!030401-07Y Calculation of table for zero shear approximation from my HYCOM inigiss_fixed2.fs0 .
!030324-26AH Make 1 dimensional table of turbulence functions of N_d^2/N^2
!         to be used for the unstable case with negligible shear.
!         N_d^2 \equiv N_Heat^2 - N_Salt^2 . N_d^2/N^2 ranges from - to + infinity.
!         N_d^2 is analogous to Ri_d^2, so I try making its table values exactly the same.
!         oursal2_zeroshear assumes Shear^2=0 and N^2 < 0.
        dand2on2 = dri
! --- Set N_d^2/N^2 table values.
!07Y	Calculate moving out from zero index first postive indices and then negative ones.[See NBp.030407-08.]
        DO ind2on2sign=0,1
        ind2on2step=(-1)**ind2on2sign
        DO ind2on2 = 0,mt*ind2on2step,ind2on2step
           IF(ABS(ind2on2).LE.mt0) THEN
             and2on2a1(ind2on2) = DFLOAT(ind2on2)*dand2on2
           ELSE
             mt0s  = SIGN(mt0,ind2on2)
             mtm1s = SIGN(mt0-1,ind2on2)
! --- introduction of exponential absolute val table option.
             IF(ifexpabstable.EQ.0) THEN
               idifs = SIGN((ABS(ind2on2)-mt0),ind2on2)
               and2on2a1(ind2on2) = and2on2a1(mt0s)*((and2on2a1(mt0s)/ &
                      and2on2a1(mtm1s))**(idifs**2))
             ELSE IF(ifexpabstable.EQ.1) THEN
               idif = ABS(ind2on2)-mt0
               and2on2a1(ind2on2) = and2on2a1(mt0s)*((and2on2a1(mt0s)/ &
                      and2on2a1(mtm1s))**(idif))
             END IF
           END IF
        END DO
	END DO
          DO ind2on2 = -mt,mt
             CALL oursal2_zeroshear &
       (and2on2a1(ind2on2),amtaun2a1(ind2on2) &
        ,sma1(ind2on2),sha1(ind2on2),ssa1(ind2on2),ttot,tcot,tctot,tptot,tpcot,tpvot)
          END DO
!******Y
	END IF
!******Z

!981215-31 Add writes in salinity model case.
       if (mytid.eq.0) then
	WRITE(*,*) " "
	WRITE(*,*) " "
	WRITE(*,*) " "
	WRITE(*,*) "************************************************"
	WRITE(*,*) " "
	WRITE(*,*) "New Temperature-Salinity Model"
	WRITE(*,*) " "
!030722Z1a Add writeout of isurfuse here (See NBp.030722-2 Vol.XIX.)
!020916X Switch for use of surface fluxes to dimensionalize turbulence model.
        WRITE(*,*) "************************************"
        WRITE(*,*) "isurfuse=",isurfuse
        WRITE(*,*) "************************************"
!******X
!*****CZ1a
	WRITE(*,*) "ifsali=",ifsali
	WRITE(*,*) "ifsalback=",ifsalback
!040217Zi1b Write out switch for use of minimum shear and amount if used.
	IF(.NOT.ifzeroshear) THEN
	  WRITE(*,*) "ifshearmin=",ifshearmin
	  IF(ifshearmin) WRITE(*,*) "s2min=",s2min
	END IF
!030424Z Write switch for whether zero shear model is enabled.
	WRITE(*,*) "ifzeroshear=",ifzeroshear
!030424Z Write parameters for rotation's effect on lengthscale.
	WRITE(*,*) "ilomega=",ilomega
	WRITE(*,*) "amldminlom=",amldminlom
!000215
!030722Z1a Add writeout of icondear here (See NBp.030722-2 Vol.XIX.)
!030721Z1a Write switch for Deardorff treatment.
        WRITE(*,*) " "
        WRITE(*,*) "icondear=",icondear
        IF(icondear.EQ.-1) THEN
          WRITE(*,*) "Do *not* use Deardorff lengthscale modification."
        ELSE IF(icondear.EQ.0) THEN
          WRITE(*,*) "Ye Cheng's old Deardorff:"// &
     "modify l and \tau N leaving S_X unmodified."
        ELSE IF(icondear.EQ.1) THEN
          WRITE(*,*) "Ye Cheng's new Deardorff:"// &
     "modify l but leave *both* \tau N and S_X unmodified."
        END IF
        WRITE(*,*) " "
!*****CZ1a
!*****CZ1a
	WRITE(*,*) "ifepson2=",ifepson2
!020219D Bottom enhancemant writes included.
	IF(ifepson2.GT.0) THEN 
        IF(ifepson2.EQ.2) WRITE(*,*) &
       "epsilon/(N^2) used even for strong mixing beneath weak mixing" 
!030722Z1a Fix writeout of ideeplat here (See NBp.030722-2 Vol.XIX.)
!030429-0502Z1 Write switch for latitude dependence of background mixing.
        WRITE(*,*) "ifdeeplat=",ifdeeplat
        IF(ifdeeplat.GT.0) THEN
          WRITE(*,*) &
        "Use latitude dependence of interior ocean value of \\epsilon/N^2"
!Zi1a',040422Zi1bj Floor on latitude dependent factor in background mixing
	  WRITE(*,*) "eplatidependmin=",eplatidependmin
        END IF
!*****CZ1a
!020404D,030502Z1 Must write out the parameter epson2__ for reference \epsilon/N^2 
!Z1	since epson2_ and epson2 are no longer constants.
          WRITE(*,*) "epson2__=",epson2__
	  WRITE(*,*) "ifbotenhance=",ifbotenhance
	  IF(ifbotenhance.EQ.1) THEN
	    WRITE(*,*) "eps_bot0=",eps_bot0
	    WRITE(*,*) "scale_bot=",scale_bot
	  END IF
        END IF
!*****CD
!000317
	WRITE(*,*)"ifrafgmax=",ifrafgmax
	WRITE(*,*) " "
	WRITE(*,*)"ifbg_theta_interp=",ifbg_theta_interp
	WRITE(*,*) " "
!011108yXI Write out whether using exponential absolute value spacing for nonlinear part of table
!          and whether using new interpolation tailored to exponential absolute value in that part.
        WRITE(*,*) " "
        WRITE(*,*) "***********************************************"
        WRITE(*,*) "ifexpabstable=",ifexpabstable
        WRITE(*,*) "ifastexpabs=",ifastexpabs
        WRITE(*,*) "***********************************************"
        WRITE(*,*) " "
!yXI
        WRITE(*,*) &
        "    i      ", &
        "    rib(i)      ","    ridb(i)     ", &
        "irimax(i)  "
!       DO i= -mt,mt
!          WRITE(*,9050) i,rib(i),ridb(i),irimax(i)
!       END DO
!*****C
	WRITE(*,*) " "
	WRITE(*,*) "irid       Ri_d        Ri(irimax)  " &
	        // "S_M        S_H        S_S        " &
          // "S_M/S_H    S_S/S_H    "
!       DO irid= -mt,mt
!          WRITE(*,9100) irid,ridb(irid),rib(irimax(irid)), &
!          smb(irimax(irid),irid), &
!          shb(irimax(irid),irid), &
!          ssb(irimax(irid),irid), &
!          smb(irimax(irid),irid)/shb(irimax(irid),irid), &
!          ssb(irimax(irid),irid)/shb(irimax(irid),irid)
!       END DO
        endif
!000316 CALCULATE "R_r_Critical" USING CANUTO'S 000228 ANALYTIC FORMULA
!	FOR "R_rho_Critical". See NBp.000229-3 and 000316-4.
!	R_rho_Canuto \equiv -Ri_C/Ri_T \equiv -R_r .
!	In a sheet dated 000228 Canuto gave me:
!	"R_\rho^{cr} = {1 \over \Deta} [1 {+\over-} \sqrt{1 - \Delta^2}] 
!	\Delta \equiv {{\pi_2(1 + {15 \over 7} \pi_3)} \over
!	{\pi_3 - \pi_2 + (15 \over 14} \pi_3^2}} ".
!	Note that the + and - choices are reciprocals so this covers
!	both the Salt Fingering and Double Diffusive Critical R_\rho's.
!	From Ocean Turbulence III paper have: 
!	\pi_{1,2,3,4,5} = 
!	(\tau_pc,\tau_c\theta,\tau_c,\tau_p\theta,\tau_\theta)/\tau 
!	R_r_Crit = [-1 -/+ \sqrt{1 - \Delta^2}]/Delta
!	\Delta = {{{\tau_c\theta \over \tau} ( 1 + (15/7)*{\tau_c \over \tau})}
!	        \over {{\tau_c \over \tau} - {\tau_c\theta \over \tau} + 
!	                (15/14) {\tau_c \over \tau}^2}}
	deltanum = tctot*(1.D0 + ((15.D0/7.D0)*tcot))
	deltaden = tcot - tctot + ((15.D0/14.D0)*(tcot**2))
	delta = deltanum/deltaden
	rrcrn = (-1.D0 - SQRT(1.D0 - (delta**2)))/delta
	rrcrp = (-1.D0 + SQRT(1.D0 - (delta**2)))/delta
	theta_rcrn = ATAN(rrcrn)
	theta_rcrp = ATAN(rrcrp)
!990111-000316 Make sure the right choice of arctan(R_r)=[\theta_r] is made.
!	Arctan covers the range (-pi/2,pi/2) while 
!	\theta_r_Crit must be in the range (-pi/4,3pi/4) (The range of Ri>0.)
        IF(theta_rcrn.LT.-pi/4.D0) theta_rcrn = theta_rcrn + pi
        IF(theta_rcrp.LT.-pi/4.D0) theta_rcrp = theta_rcrp + pi
	theta_rcrn_deg = theta_rcrn*(180.D0/pi)
	theta_rcrp_deg = theta_rcrp*(180.D0/pi)
       if (mytid.eq.0) then
	WRITE(*,*) "   "
	WRITE(*,*) "   "
	WRITE(*,*) "   "
	WRITE(*,*) "   "
	WRITE(*,*) "R_r_Crit+ =",rrcrp
	WRITE(*,*) "R_r_Crit- =",rrcrn
	WRITE(*,*) "\\theta_r_Crit+ =",theta_rcrp
	WRITE(*,*) "\\theta_r_Crit- =",theta_rcrn
	WRITE(*,*) "\\theta_r_Crit+ in degrees =",theta_rcrp_deg
	WRITE(*,*) "\\theta_r_Crit- in degrees =",theta_rcrn_deg
	WRITE(*,*) "   "
	WRITE(*,*) "   "
!*****C
	WRITE(*,*) " "
	WRITE(*,*) " "
       endif
!	Increments in radial and angular coordinates in (Ri_T,Ri_C) plane.
	   delra_r = 1.D0/DFLOAT(mt_ra_r)
	deltheta_r = (pi/4.D0)/DFLOAT(n_theta_r_oct)
!	Calculate the ratio \sigma_sa_max \equiv S_S/S_H as a function 
!	of the angle \theta_r in Ri_T,Ri_C space,
!       \theta_r \equiv arctan(Ri_C/Ri_T). 
!981218 The range of angles where unrealizability occurs is 
!	a subset of theta_r = -pi/4 to 3pi/4.
       if (mytid.eq.0) then
	WRITE(*,*) "S_S/S_H at pre-maximum Ri as a function of" &
             // "\\theta_r \\equiv Arctan(Ri_C/Ri_T)" 
!981220    Absurd default on sisamax \equiv S_S/S_H.
	WRITE(*,*) "Arbitrarily show the absurd value -99.999"
	WRITE(*,*) "at angles where do not have "// &
	"a maximum Ri (or radius ra_r)."
	WRITE(*,*) " "
	WRITE(*,*) "  \\th_r ^o  ra_r      " &
          // "  Ri_T        Ri_C        Ri         Ri_d       " &
          // "  S_M       S_H       S_S      S_S/S_H  "
       endif
!*
!       For Ri_T and Ri_C positive find the realizability limits  
!	in polar coordinates in the (Ri_T,Ri_C) plane : (ra_r,theta_r).
	IF(ifpolartablewrite.EQ.1) then
       if (mytid.eq.0) then
        OPEN(UNIT=68,FILE="turb_ra_th",STATUS="NEW")
       endif
       endif
	DO itheta_r = -n_theta_r_oct,3*n_theta_r_oct
	   theta_r = DFLOAT(itheta_r)*deltheta_r
	   theta_r_deg = theta_r*(180.D0/pi)
!980119	Introduce jtheta_r, an angle index that begins at zero   
!	for the purposes of letting OURSAL2 know it starts at the origin.
	   jtheta_r = itheta_r + n_theta_r_oct
!981220-1 Initialize sisamax to the impossible negative value of -99.999 to 
!	let places where the realizability limit is not reached stand out.
	   sisamax(itheta_r) = -99.999
!981221  Initialize sm_r0,sh_r0,ss_r0 to the INCONSISTENT absurd value -9.999999.
	   sm_r0 = -9.999999
	   sh_r0 = -9.999999
	   ss_r0 = -9.999999
!990303 Flag ibg determines if the background value of ra_r has been calculated.
!030116X1a Option had never worked. George Halliwell points out I erroneously had "isalback" here before.
	   IF(ifsalback.EQ.6) ibg=0
!*****C
!990315 Flag ifunreal determines if realizability limit has been found.
	   ifunreal=0
!000315 Write to file "back_ra_r_values", the allowed values of background ra_r.
	   IF(itheta_r.EQ.-n_theta_r_oct) then
       if (mytid.eq.0) then
!lhl0711	   OPEN(UNIT=69,FILE="back_ra_r_values",STATUS="NEW")
	   OPEN(UNIT=69,FILE="back_ra_r_values")
       endif
       endif
!*****C
!000318 Make the ra_r max value not too large to try to avoid numerical trouble.
	   DO ira_r = 0,(mt_ra_r**2)/4
	      IF(ira_r.LE.mt_ra_r) THEN
                ra_r = DFLOAT(ira_r)*delra_r
	      ELSE
	        ra_r = ((1.D0+delra_r)**(ira_r - mt_ra_r)) &
                *(DFLOAT(mt_ra_r)*delra_r)
	      END IF
       if (mytid.eq.0) then
           IF(itheta_r.EQ.-n_theta_r_oct) WRITE(69,9200) ira_r,ra_r
       endif
!981216-990119	Convert radius and angle, (ra_r,theta_r), to rectangular coordinates.
	      rit = ra_r*COS(theta_r)
	      ric = ra_r*SIN(theta_r)
	      ri_r  = rit + ric
	      rid_r = rit - ric
!981216 Calculate turbulence functions at this radius and angle in (Ri_T,Ri_C).
              IF(nmodel.EQ.1) THEN
!021210X1 Call version of submodule oursal2 which has an option to use B1={\tau S}^{3/2} (0,0).
	        CALL OURSAL2_2A(b1,ri_r,rid_r,slq2_r, &
                       sm_r,sh_r,ss_r, &
                       c_y0,c_y00,ira_r,jtheta_r,ii,jj,ttot,tcot,tctot,tptot,tpcot,tpvot)
	      ELSE IF(nmodel.EQ.2) THEN
	        CALL MIKESAL2A(b1,ri_r,rid_r,slq2_r, &
                       sm_r,sh_r,ss_r, &
                       c_y0,c_y00,ira_r,jtheta_r)
	      END IF
!      if (mytid.eq.0) then
!             IF(ifpolartablewrite.EQ.1) WRITE(68,9001) &
!        itheta_r,theta_r_deg,ira_r,ra_r,slq2_r,sm_r,sh_r,ss_r
!      endif
!990303 Calculate S_M/(S l/q) and find where it's backfact of its origin value.
	      IF(ifsalback.EQ.6) THEN
	        smosloq_r = sm_r/SQRT(slq2_r)
	        IF(ira_r.EQ.0) smosloq_0r = smosloq_r
!	Use radius where dimensionless K_M falls below backfact*origin value.
	        IF((smosloq_r.LE.backfact*smosloq_0r).AND.(ibg.EQ.0)) &
           THEN
	          ra_r1            	= ra_r
	          rit1    	   	= rit
	      	  ric1  	 	= ric
	          ri_r1   		= ri_r
		  rid_r1  		= rid_r
	          slq2_r1(itheta_r)	= slq2_r
	          sm_r1(itheta_r)  	= sm_r
	          sh_r1(itheta_r)	= sh_r
	          ss_r1(itheta_r)	= ss_r
	          ibg=1
	        END IF
!*****C
	      END IF
	      IF(slq2_r.LE.0.D0) THEN 
!981216	Use value of last lattice point on this radius with "slq2" positive.
!	Calculate the ratio of the salt and heat diffusivities there.
	sisamax(itheta_r) = ss_r0/sh_r0 
!990301 Store in an array the maximum radius, ra_r, at this angle, theta_r,
!	in the polar (Ri_T,Ri_C) [that is the (theta_r,ra_r)] plane.
	ra_rmax(itheta_r) = ra_r0
!	Determine the background radius, ra_r, at this \theta_r.
	IF(ifsalback.EQ.5) THEN
!       Use a constant fraction of the maximum radius before model breakdown.
	  back_ra_r(itheta_r) = backfrac*ra_rmax(itheta_r)
!990303-15
	ELSE IF(ifsalback.EQ.6) THEN
	  back_ra_r(itheta_r) = ra_r1
	END IF
	ifunreal = 1 
!*****C
!981230 Skip straight to write out when last point reached.
	GO TO 16
	      END IF
	      ra_r0   = ra_r
	      rit0    = rit
	      ric0    = ric
	      ri_r0   = ri_r
	      rid_r0  = rid_r
	      slq2_r0 = slq2_r
	      sm_r0   = sm_r
	      sh_r0   = sh_r
	      ss_r0   = ss_r
!000319 Store c_y as c_y_0 for possible use as a  guess in background calc. .
	      c_y_r0(itheta_r) = c_y0
	   END DO
!000315 Close file with background ra_r values.
       if (mytid.eq.0) then
	IF(itheta_r.EQ.-n_theta_r_oct) CLOSE(69)
       endif
!981216-30 Write out stability functions, the S's and sisamax.
   16  continue
!      if (mytid.eq.0) then
!  	WRITE(*,9150) theta_r_deg,ra_r0,rit0,ric0,ri_r0,rid_r0, &
!               sm_r0,sh_r0,ss_r0,sisamax(itheta_r)
!      endif
!990315 Set background ra_r large at angles where unrealizability doesn't occur.
!000318 Make the ra_r max value not too large to try to avoid numerical trouble.
        IF(ifunreal.EQ.0) THEN
           ipenra_r = (mt_ra_r**2)/4-1
	   back_ra_r(itheta_r) = ((1.D0+delra_r)**(ipenra_r - mt_ra_r)) &
                         *(DFLOAT(mt_ra_r)*delra_r) 
	END IF
!*****C
!990315 For ifsalback=5 case get value for initialization of c_y calculation. 
	      IF(ifsalback.EQ.5) THEN
	        IF(jtheta_r.EQ.0) THEN
	          c_y001 = c_y0
	        END IF
	      END IF
!*****C
	END DO
       if (mytid.eq.0) then
	IF(ifpolartablewrite.EQ.1) CLOSE(68)
       endif
!990303-16 Write out stability functions at background ra_r .
	IF(ifsalback.GT.4) THEN
	  DO itheta_r = -n_theta_r_oct,3*n_theta_r_oct
	     theta_r = DFLOAT(itheta_r)*deltheta_r
	     theta_r_deg = theta_r*(180.D0/pi)
!981216-990119	Convert radius and angle, (ra_r,theta_r), to rectangular coordinates.
	     rit1 = back_ra_r(itheta_r)*COS(theta_r)
	     ric1 = back_ra_r(itheta_r)*SIN(theta_r)
	     ri_r1  = rit1 + ric1
	     rid_r1 = rit1 - ric1
!990315-16 Calculation of turbulence functions for ifsalback=5 case.
	     IF(ifsalback.EQ.5) THEN
!981216 Calculate turbulence functions at this radius and angle in (Ri_T,Ri_C).
	       jtheta_r = itheta_r + n_theta_r_oct
!990315 Set second table index to 1 to use last step's value except at start.
!000319 Transform that "last step" value from the most recent angle step to the
!	final realizable ra_r step at {\it this} angle in hope of more accuracy.
               IF(nmodel.EQ.1) THEN
!021210X1 Call version of submodule oursal2 which has an option to use B1={\tau S}^{3/2} (0,0).
	         CALL OURSAL2_2A(b1,ri_r1,rid_r1,slq2_r1(itheta_r), &
                 sm_r1(itheta_r),sh_r1(itheta_r),ss_r1(itheta_r), &
                        c_y_r0(itheta_r),c_y001,jtheta_r,1,ii,jj,ttot,tcot,tctot,tptot,tpcot,tpvot)
	       ELSE IF(nmodel.EQ.2) THEN
	         CALL MIKESAL2A(b1,ri_r1,rid_r1,slq2_r1(itheta_r), &
                 sm_r1(itheta_r),sh_r1(itheta_r),ss_r1(itheta_r), &
                        c_y_r0(itheta_r),c_y001,jtheta_r,1)
	       END IF
	     END IF
       if (mytid.eq.0) then
	     IF(itheta_r.EQ.-n_theta_r_oct) THEN
	       WRITE(*,*) " "
	       WRITE(*,*) &
          "Values at background ra_r=(Ri_T^2 + Ri_C^2)^(1/2)"
	       WRITE(*,*) "\\th_r ^o   ra_r       " &
           // "Ri_T       Ri_C       Ri         Ri_d       " &
           // "(Sl/q)^2   S_M       S_H       S_S       S_S/S_H  "
	       WRITE(*,*) " "
	     END IF
       endif
	     sisa1 = ss_r1(itheta_r)/sh_r1(itheta_r)
!      if (mytid.eq.0) then
!    	     WRITE(*,9160) theta_r_deg,back_ra_r(itheta_r), &
!                    rit1,ric1,ri_r1,rid_r1,slq2_r1(itheta_r), &
!                 sm_r1(itheta_r),sh_r1(itheta_r),ss_r1(itheta_r), &
!                 sisa1
!      endif
	     IF(slq2_r1(itheta_r).LT.0.D0) THEN
       if (mytid.eq.0) then
	     WRITE(*,*) &
        "Negative (Sl/q)^2 in table of Background vs. \\theta_r."
	     WRITE(*,*) "itheta_r=",itheta_r, &
	                 "   slq2_r1(itheta_r)=",slq2_r1(itheta_r)
	     WRITE(*,*) "Program is stopping in turb_2."
       endif
	     STOP
	     END IF
  	  END DO
	END IF
!*****C
!*
       if (mytid.eq.0) then
	WRITE(*,*) " "
	WRITE(*,*) "************************************************"
	WRITE(*,*) " "
!030404Y Write out to standard output and a file a table of dimensionless turbulence 
!Y	functions of N_d^2/N^2 for the zero shear unstable case.
!030424Z Only write in case zeroshear model enabled.
	IF(ifzeroshear) THEN 
	CLOSE(68)
	WRITE(*,*) "************************************************"
	WRITE(*,*)          "index (N_d^2)/(N^2)  -(\\tau N)^2    "// &
             "S_M            S_H            S_S            "
	WRITE(*,*) " "
!lhl0711	OPEN(UNIT=68,FILE="turb_nd2on2",STATUS="NEW")
	OPEN(UNIT=68,FILE="turb_nd2on2")
	WRITE(68,*) "************************************************"
	WRITE(68,*)          "index (N_d^2)/(N^2)  -(\\tau N)^2    "// &
             "S_M            S_H            S_S            "
	WRITE(68,*) " "
	DO ind2on2 = -mt,mt
!           WRITE(*,9268) ind2on2,and2on2a1(ind2on2), &
!     amtaun2a1(ind2on2),sma1(ind2on2),sha1(ind2on2),ssa1(ind2on2)
	   WRITE(68,9268) ind2on2,and2on2a1(ind2on2), &
      amtaun2a1(ind2on2),sma1(ind2on2),sha1(ind2on2),ssa1(ind2on2)
	END DO
	WRITE(*,*) " "
	CLOSE(68)
	WRITE(*,*) "************************************************"
!******Y
	END IF
!******Z
	WRITE(*,*) " "
	WRITE(*,*) " "
       endif
!******
   17   CONTINUE
        ifirst=1
!******
 9001 FORMAT(2(I8,'   ',1pe11.3),8(1pe11.3))
 9050 FORMAT(I8,'  ',2E16.4,I8,'  ')
 9100 FORMAT(I8,'  ',2E12.4,3F11.6,2F11.4)
 9150 FORMAT(F11.3,5E12.4,3F10.6,F9.3)
 9160 FORMAT(F11.3,6E10.4,3F10.6,F9.3)
 9200 FORMAT(I12,'    ',5E16.6)
 9268 FORMAT(I6,5F15.9)
      return
      end

	   SUBROUTINE oursal2_2a(b1_arg, &
                        ri,rid,slq2,sm,sh,sc,c_y0,c_y00,iri,irid,ii,jj, &
                        ttot,tcot,tctot,tptot,tpcot,tpvot)
!****************************************
!981016    Stripped and adapted from plot981007.f.
!**********
!981007-15 Program to generate contour and 1 variable plots vs. Ri,Ri_d based on
!981001-02 Program to generate contour plots vs. Ri_T and Ri_C based on
!980728-29 Program to generate plots vs. Ri_T at different Ri_C values based on
!980728	   CORRECTED PROGRAM WITH NEW VALUE OF "p10". 'p10 = tpt*tct/(tc**2)'
!980723-24 Program to generate K_X/((l^2) S) for Canuto based on plot980609.f:
!980609-16 Program to generate data for plots of turbulence functions including
!	S_{M,H,C} and Canuto's new y = (\tau_pv S)^2
!	and n,c as functions of stability parameters in the concentration theory
!	(structure is a 1 point closure like the generalized Mellor-Yamada, 
!        but the constants are derived based on Dubovikov's model according
!	 to Ye Cheng). The concentration theory dimensionless parameters
!	 associated with the squares of shear, temperature contribution to 
!        Brunt Vaisala frequency and concentration contribution to it,
!	 the new y,n,c are represented in this program by the variables
!        c_y,c_n,c_c. 
!	Adapted from Cheng's program mike_12.f_980528 for the Dubovikov model.

!      4/20/98
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     y=(tau*s)**2
!     tau=2*e/epsilon=b1*l/q
!     km=e*tau*sm=1/2*(b1*l)**2*s/y**(1/2)*sm
!     kh=e*tau*sh=1/2*(b1*l)**2*s/y**(1/2)*sh
!1016 ks=e*tau*ss=1/2*(b1*l)**2*s/y**(1/2)*ss
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!980723-1022  X = {M,H,C} . 
!	Cheng above gives K_X = (1/2)((B_1*l)^2) (S/(((\tau S)**2)^(1/2))) S_X
!	The "old" y used above is (\tau S)^2. 
!	The "new" y (c_y in the program) is (\tau_pv S)^2.
!	The program variable "slq2" is (S l/q)^2 = y (B_1)^(-2), 
!	since \tau=B_1 l/q. (S l/q)^2 = (\tau \over \tau_pv)^2 c_y (B_1)^(-2) .
!	c_y = (S l/q)^2 * [(B_1)^2 * (\tau_pv \over \tau)^2] .
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!000210 Take \tau_pv/\tau as being calculated in the smshsc routine instead.
!C980615	From the printed notes Canuto gave me on 980601 have:
!C	 \tau_pv = {2 \over 5} \tau			(B.1)
!         PARAMETER(tpvot = 0.4D0)
!C****** "tpv/tau" = 2/5
!980723 Ye Cheng says today that although part of the coefficients here come
!	from Dubovikov's model nevertheless B_1 is 16.6 in this model.
!021210A Rename the parameter from `b1' to "b1_0" to allow calculation of B1 by model.
	PARAMETER(b1_0=16.6D0)
!021210A Replace b1_0 by Cheng's value of  y(Ri_Tem=0,Ri_Sal=0)^{3/4} when ifchengb1 is .TRUE. .
	LOGICAL ifchengb1
	PARAMETER(ifchengb1=.FALSE.)
	INTEGER ib1set ! Flag positive if B1 has already been assigned its value.
	DATA ib1set / 0 /
!021211A ib1set at call of routine.
	INTEGER ib1set0
!******A
!980610-030403 Common block with ratios of timescales
!      real(r8)::  ttot,tcot,tctot,tptot,tpcot,tpvot
!980609 rit is the temperature's part of Ri and ric the concentration's.
      real(r8) ::rit,ric
      external fct_sal
!980609-1016  Need a guess for c_y for the solver for the neutral case, c_yst.
!980609  Take c_yst = 8.527882, the approximate value calculated at rit=ric=0.
      PARAMETER(c_yst0 = 8.527882D0)
!981016  A variable c_y00 is intended to hold the Ri=0 value of c_y from 
!	 the previous Ri_d row in a table the subroutine is being called to make
!	 and a variable c_y0 is intended to hold the previous Ri value
!	 from the current Ri_d row of that table.
!******
        integer :: ii,jj
!981104	Give B_1 an alias needed for its fortran role as an output argument.
!021211A Store incoming value of ib1set.
         ib1set = 0
	 ib1set0=ib1set
!021210A Calculate B1 later if use Cheng's formula instead of an a priori constant. 
	IF(.NOT.ifchengb1) THEN
!030124B Only set B1 when it has not been set before. [See NBp.030124-18 Vol.XVII .] 
	  IF(ib1set0.EQ.0) THEN
!030128B Set variable "b1" to parameter "b1_0". {See NBp.030128-4 .]
            b1 = b1_0
	    b1_arg=b1
	    ib1set=ib1set+1
	  END IF
	END IF
!*****C
!       Make dummy call with c_y=c_n=c_c=0 to get their values for initial use.
      call smshsc_a3(0.D0,0.D0,0.D0,sm,sh,sc,ttot,tcot,tctot,tptot,tpcot,tpvot)
!
      eps=1.D-6                                              
      iend=300                                              
!     rimax= ?
!-----rtwi finds the root of x=fct(x)                     
!980615  Need a guess at the root, c_yst. Use neighboring solution.
!981007 Initial guess for c_yst for this value of Ri_d.
      IF(iri.EQ.0.AND.irid.EQ.0) THEN
	c_yst = c_yst0
      ELSE IF(iri.EQ.0) THEN
	c_yst = c_y00
      ELSE 
	c_yst = c_y0
      END IF
!981007 Calculate Ri_T =(Ri + Ri_d)/2 and Ri_C =(Ri - Ri_d)/2.  
	 rit = (ri + rid)/2.D0
	 ric = (ri - rid)/2.D0
         call rtwi(c_y,val,fct_sal,c_yst,eps,iend,ier,ttot,tcot,tctot,tptot,tpcot,tpvot,ri,rit1,ric1,rit,ric)

         if(ier.ne.0) then
!981022-23 Make error message more specific.
       if (mytid.eq.0) then
	    WRITE(*,*) "In oursal2 subroutine"
	    WRITE(*,*) "c_y00=",c_y00,"	c_y0=",c_y0
	    WRITE(*,*) "ri=",ri,"	rid=",rid
	    WRITE(*,*) "rit=",rit,"	ric=",ric
	    WRITE(*,*) "Initial guess for rtwi c_yst=",c_yst
!*****C
            write(*,*) "rtwi call problem, ier=",ier
       endif
            stop
         endif
!981016 Store value of c_y for future guesses.
	 IF(c_y.GE.0) THEN
	    c_y0=c_y
	 ELSE 
!981022 Turbulence model becomes unphysical for c_y negative.
!981014-16  Realizability for negative Ri
           IF(ri.LT.0) THEN
       if (mytid.eq.0) then
	     WRITE(*,*) "c_y negative at negative Ri"
	     WRITE(*,*) "Ri=",ri," 	c_y=",c_y
	     WRITE(*,*) "Unstable realizability limit unexpected:" 
	     WRITE(*,*) "stopping in oursal2."
       endif
	     STOP
           END IF
	 END IF
!
!
	 IF(iri.EQ.0) c_y00=c_y
	 IF((iri.EQ.0).AND.(irid.EQ.0).AND.  &
      (ABS(c_y - c_yst0).GT.1.D-6)) THEN
       if (mytid.eq.0) then
        write(*,*) "ii,jj=",ii,jj
	   WRITE(*,*) "Inconsistency in neutral value of c_y"
	   WRITE(*,*) "Value used =",c_yst0
	   WRITE(*,*) "Value calculated =",c_y
	   WRITE(*,*) "Program stopping in oursal2"
       endif
	   STOP
	 END IF
!*****C	 
!021210A *Calculate Cheng's value of B1 \equiv y(Ri=0,Rid=0)^{3/4} .*
!A	 y \equiv {\tau S}^2 \equiv {\tau_pv S}^2 (\tau/\tau_pv)^2 \equiv c_y (\tau_pv/\tau)^{-2}
	IF(ifchengb1) THEN
	  IF(ib1set.EQ.0) THEN
	    taus2 = c_y/(tpvot**2)
	    taus3 = taus2*SQRT(taus2)
	    b1 = SQRT(taus3)
	    b1_arg = b1
	    ib1set=ib1set+1
	  END IF
	END IF
	IF(b1_arg.LE.0.D0) THEN
       if (mytid.eq.0) then
	  WRITE(*,*) "B1 <= ZERO."
	  WRITE(*,*) "B1=",b1_arg
	  WRITE(*,*) "Something must be wrong."
	  WRITE(*,*) "Program is stopping in oursal2."
       endif
	  STOP
	END IF
	IF(ib1set.NE.1) THEN
       if (mytid.eq.0) then
	  WRITE(*,*) "Problem in oursal2; B1 not properly set."
	  WRITE(*,*) "Number of times B1 set=",ib1set
	  WRITE(*,*) "b1_arg=",b1_arg
	  WRITE(*,*) "Program is stopping."
       endif
	  STOP
	END IF
!021211A Only do B1 writeout when first set. [See NBp021211-2.]
!    IF(ib1set0.EQ.0) THEN
!    if (mytid.eq.0) then
! WRITE(*,*) " "
! WRITE(*,*) "ifchengb1=",ifchengb1
! WRITE(*,*) "B_1=",b1_arg
! WRITE(*,*) " "
!      endif
!END IF
	ib1set0=ib1set
!******A
!021210A Move calculation of (S l/q)^2 until after the new B1 has been determined.
!981027 **Calculate (S l/q)^2[=program variable "slq2"] from c_y.**
!	 (S l/q)^2 = (\tau \over \tau_pv)^2 c_y (B_1)^(-2) .
!	 (S l/q)^2 = (\tau_pv \over \tau)^(-2) c_y (B_1)^(-2) .
	 slq2 = c_y/((b1*tpvot)**2)
!*****C
!980610 From last page (#5) of "980608 AH Concentration Work" handwritten sheetsC	have: 
!        n = -{{\tau_C \tau_{C\theta}} \over {\tau_{pv}}^2 } y Ri_T
!	 c = - {{\tau_C}^2 \over {\tau_{pv}}^2} y Ri_C
!000210 Decide to use the parameter "tpvot" instead of its value 2/5 \tau .
!        n = -{{(\tau_C/\tau) (\tau_{C\theta}/\tau)} \over {\tau_{pv}/\tau}^2 }
!           y Ri_T
!        c = - {{\tau_C/\tau}^2 \over {\tau_{pv}/\tau}^2} y Ri_C
!******
         c_n = -(tcot*tctot/(tpvot**2))*c_y*rit
         c_c = -((tcot**2)/(tpvot**2))*c_y*ric
         call smshsc_a3(c_y,c_n,c_c,sm,sh,sc,ttot,tcot,tctot,tptot,tpcot,tpvot)


 1003 format(12(I8))
 1004 format(12(1pe14.5))
      end



!030403Y Include "ttot'=`{\tau_\theta \over \tau}' in the common block with timescale 
!Y	 ratios, /bb0/. See NBp030403-8to11.
      function fct_sal(c_y,ttot,tcot,tctot,tptot,tpcot,tpvot,ri,rit1,ric1,rit,ric)
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
      real(r8)::  ttot,tcot,tctot,tptot,tpcot,tpvot
      real(r8)::  ri,rit1,ric1,rit,ric
!000210 Decide to use the parameter "tpvot" instead of its value 2/5 \tau .
      c_n = -((tcot*tctot)/(tpvot**2))*c_y*rit
      c_c = -((tcot**2)/(tpvot**2))*c_y*ric
      call smshsc_a3(c_y,c_n,c_c,sm,sh,sc,ttot,tcot,tctot,tptot,tpcot,tpvot)
!980609 y(S_\nu - Ri_T S_h - Ri_C S_c) = 8/25 . 8/25 = 0.32 . S_\nu = sm.
!	y = 0.32/(S_\nu - Ri_T S_h - Ri_C S_c). 
      fct_sal=(2.D0*(tpvot**2))/(sm-rit*sh-ric*sc)
      return                                          
      end                                            
                                                    


      subroutine smshsc_a3(y,n,c,sm,sh,sc,ttot,tcot,tctot,tptot,tpcot,tpvot)
!030403Y Include "ttot'=`{\tau_\theta \over \tau}' in the common block with timescale 
!Y       ratios, /bb0/. See NBp030403-8to11.
!000125-27 NEW SUBROUTINE WHICH calculates the "p's" from the timescale ratios.
!	BASED on "smshsc2":
!990513 SUBROUTINE WHICH CALCULATES "p's" from "sgmt". BASED ON "smshsc1":
!990513 NEW SUBROUTINE WHICH USES YE CHENG'S FORTRAN CODE TO CALCULATE CONSTANTS
!	FROM THE "p's" SENT TO ME BY HIM TODAY. BASED ON "smshsc0". 
!980728 **CORRECT THE VALUE OF "p10".**
!	p_10 = {\tau_{p \theta} \tau_{c \theta}} \over {\tau_c ^ 2}
!******
!980609-15 Replace Cheng's smsh with  smshsc, which includes concentration.
!980610 The y,n,c used here are Canuto's "y,n,c" called c_y,c_n,c_c 
!	elsewhere in this program.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
      REAL(r8) y,n,c
      REAL(r8) Nm,Nh,Nc
!      implicit real*8 (a-h,o-z)
!990513 Switch for whether(1) or not(0) to output p's a's and d's to a file.
      PARAMETER(ifmodelconstout=0)
!*****C
!000127	Add `\tau_pv \over \tau' to the common block with timescale ratios.
      real(r8)::  ttot,tcot,tctot,tptot,tpcot,tpvot
!980615	From the printed notes Canuto gave me on 980601 have:
!	 \tau_pv = {2 \over 5} \tau			(B.1)
         PARAMETER(tpvot0 = 0.4D0)
!****** "tpv/tau" = 2/5
!990513 Make "sgmt" a parameter. Standard value was 0.72.
 	PARAMETER(sgmt=0.72D0)
!000125	Determine {\tau_p\theta,tau_pc,tau_\theta,tau_c,tau_c\theta }/tau
	PARAMETER(tptot0=(1.D0/5.D0)*(1.D0/(1.D0+(1.D0/sgmt))))
	PARAMETER(tpcot0=tptot0)
	PARAMETER(ttot0=sgmt)
	PARAMETER(tcot0=ttot0)
	PARAMETER(tctot0=1.D0/3.D0)
!*****C

!****************************************
!000125 Set common block passable timescale ratios to parameter statement values
!000127 Include \tau_pv \over \tau .
	tpvot = tpvot0
	tptot = tptot0
	tpcot = tpcot0
	ttot  = ttot0
	tcot  = tcot0
	tctot = tctot0
!*****C
!990513 p calculation taken from today's Cheng email, cheng990513.result.2_1.
!results.2_1
!#  if sgmt=0.72 then
!
!      P1 = 0.832E0
!      P2 = 0.545E0
!      P3 = 0.2093023E0
!      P4 = 0.3229974E-1
!      P5 = 0.1550388E-1
!      P6 = 0.2422481E0
!      P7 = 0.48E0
!      P8 = 0.2093023E0
!      P9 = 0.872093E0
!      P10 = 0.1550388E-1
!      P11 = 0.1162791E0
!      P1m = 0.168E0
!      P2m = 0.455E0
!
!     P1 = 0.832E0
!     P2 = 0.545E0
!     P3 = 5.0/2.0/(5+5/sgmt)
!     P4 = 1/(5+5/sgmt)/sgmt**2/5
!     P5 = 2.0/15.0/(5+5/sgmt)/sgmt
!     P6 = 3.0/2.0/(5+5/sgmt)/sgmt**2
!     P7 = 2.0/3.0*sgmt
!     P8 = 5.0/2.0/(5+5/sgmt)
!     P9 = 15.0/2.0/(5+5/sgmt)/sgmt
!     P10 = 2.0/15.0/(5+5/sgmt)/sgmt
!     P11 = 1/(5+5/sgmt)/sgmt
!     P1m = 0.168E0
!     P2m = 0.455E0
 
!000125 Calculate the p's.
      p1  = 0.832D0
      p2  = 0.545D0
      p3  = (5.D0/2.D0)*tpcot
      p4  = (1.D0/5.D0)*tpcot*(tcot**(-2))
      p5  = tpcot*tctot*(tcot**(-2))
      p6  = (1.D0/5.D0)*(tcot**(-1))*(tctot**(-1))*tptot
      p7  = 5.D0*tctot
      p8  = (5.D0/2.D0)*tptot
      p9  = ttot*tptot*((tcot*tctot)**(-1))
      p10 = tctot*tptot*(tcot**(-2))
      p11 = tpcot*(tcot**(-1))
      p1m = 1.D0 - p1
      p2m = 1.D0 - p2
!*****C
!results.2_1
!990513 Values of a's and d's calculated from p's using Cheng's Fortran code
!	to do so, from today's email from him, cheng990513.results.2_1 .
!results.2_1
!##########################
!##  Fortran code:
!##########################
      A0 = 12
      A1 = p11*(12*p9+8*p6-30*p6*p8-5*p6*(p1m+3*p2m))
      A2 = 5*(2*p4*p6*p7-p4*p9-p6*p11)*(p1m+3*p2m)+8*p6*p11+8*p4*p9-16*p4&
       *p6*p7+12*p11*p9+12*p11*p10-12*p4*p7**2*p6-30*p6*p11*p8+30*p4*p6*&
        p7*p8+30*p6*p4*p7*p3-30*p4*p9*p3
      A3 = p10*(12*p11+8*p4-30*p3*p4-5*p4*(p1m+3*p2m))
      A4 = -p6*(8-30*p8-5*p1m-15*p2m)-12*p9-12*p11
      A5 = -p4*(8-30*p3-5*p1m-15*p2m)-12*p10-12*p11
      D0 = 24
      D1 = p11*((-p6-2*p9)*p1m**2+(p6+6*p9)*p2m**2+2*p6*p8*(p1m-3*p2m))
      D2 = (2*p4*p6*p7-p4*p9-p6*p11)*(p1m**2-p2m**2)+2*(-p11*p10-p11*p9+&
        p4*p7**2*p6)*(p1m**2-3*p2m**2)+2*(-p6*p4*p7*p3-p4*p6*p7*p8+p4*p9*p3&
        +p6*p11*p8)*(p1m-3*p2m)
      D3 = p10*((-p4-2*p11)*p1m**2+2*p4*p3*(p1m-3*p2m)+(6*p11+p4)*p2m**2)
      D4 = -4*p6*p11*(3*p9+2*p6)
      D5 = 4*p4*p6**2*p7*(4+3*p7)-4*p4*p9*(3*p11+2*p6)-4*p6*p11*(3*p9+3*&
          p10+2*p4+2*p6)
      D6 = 4*p4**2*p6*p7*(4+3*p7)-4*p4*p9*(2*p4+3*p11)-8*p4*p6*(p11+p10)&
          -12*p10*p11*(p4+p6)
      D7 = -4*p4*p10*(2*p4+3*p11)
      D8 = (2*p9+2*p11+p6)*p1m**2-2*p6*p8*(p1m-3*p2m)-(p6+6*p9+6*p11)*p2m**2
      D9 = (2*p10+p4+2*p11)*p1m**2-2*p4*p3*(p1m-3*p2m)-(p4+6*p10+6*p11)*p2m**2
      D10 = 8*p6**2+4*(7*p11+3*p9)*p6+24*p11*p9
      D11 = -8*(4+3*p7)*p4*p6*p7+4*p4*(4*p6+7*p9+3*p11)+4*p6*(3*p10+7*p11)&
             +24*p11*(p10+p9)
      D12 = 4*p10*(7*p4+6*p11)+4*p4*(2*p4+3*p11)
      D13 = 6*p2m**2-2*p1m**2
      D14 = -28*p6-24*p9-24*p11
      D15 = -24*p10-28*p4-24*p11
!results.2_1
!****************************************

!980728	Write out the p's.
!000125 Writeout the timescale ratios as well.
!lhl
	ifrecall=1
!lhl
	IF(ifrecall.EQ.0) THEN
       if (mytid.eq.0) then
	  WRITE(*,*) "tau_pv/tau     =",tpvot 
	  WRITE(*,*) "tau_ptheta/tau =",tptot
	  WRITE(*,*) "tau_pc/tau =",tpcot
	  WRITE(*,*) "tau_theta/tau  =",ttot
	  WRITE(*,*) "tau_c/tau  =",tcot
	  WRITE(*,*) "tau_ctheta/tau  =",tctot
	  WRITE(*,*) " "
	  WRITE(*,*) "p1 =",p1
	  WRITE(*,*) "p2 =",p2
	  WRITE(*,*) "p3 =",p3
	  WRITE(*,*) "p4 =",p4
	  WRITE(*,*) "p5 =",p5
	  WRITE(*,*) "p6 =",p6
	  WRITE(*,*) "p7 =",p7
	  WRITE(*,*) "p8 =",p8
	  WRITE(*,*) "p9 =",p9
	  WRITE(*,*) "p10=",p10
	  WRITE(*,*) "p11=",p11
!990513 Write out the a's and d's as well.
	  WRITE(*,*) "a0=",a0
	  WRITE(*,*) "a1=",a1
	  WRITE(*,*) "a2=",a2
	  WRITE(*,*) "a3=",a3
	  WRITE(*,*) "a4=",a4
	  WRITE(*,*) "a5=",a5
	  WRITE(*,*) "d0=",d0
	  WRITE(*,*) "d1=",d1
	  WRITE(*,*) "d2=",d2
	  WRITE(*,*) "d3=",d3
	  WRITE(*,*) "d4=",d4
	  WRITE(*,*) "d5=",d5
	  WRITE(*,*) "d6=",d6
	  WRITE(*,*) "d7=",d7
	  WRITE(*,*) "d8=",d8
	  WRITE(*,*) "d9=",d9
	  WRITE(*,*) "d10=",d10
	  WRITE(*,*) "d11=",d11
	  WRITE(*,*) "d12=",d12
	  WRITE(*,*) "d13=",d13
	  WRITE(*,*) "d14=",d14
	  WRITE(*,*) "d15=",d15
!990513 Output p#, a# and d# to the file model_constants if the switch is set.
!000125 Writeout the timescale ratios as well.
          IF(ifmodelconstout.EQ.1) THEN
            OPEN(UNIT=66,FILE='model_constants',STATUS='UNKNOWN')
	    WRITE(*,*) "tau_pv/tau     =",tpvot 
	    WRITE(*,*) "tau_ptheta/tau =",tptot
	    WRITE(*,*) "tau_pc/tau =",tpcot
	    WRITE(*,*) "tau_theta/tau  =",ttot
	    WRITE(*,*) "tau_c/tau  =",tcot
	    WRITE(*,*) "tau_ctheta/tau  =",tctot
	    WRITE(*,*) " "
	    WRITE(66,*) "p1 =",p1
	    WRITE(66,*) "p2 =",p2
	    WRITE(66,*) "p3 =",p3
	    WRITE(66,*) "p4 =",p4
	    WRITE(66,*) "p5 =",p5
	    WRITE(66,*) "p6 =",p6
	    WRITE(66,*) "p7 =",p7
	    WRITE(66,*) "p8 =",p8
	    WRITE(66,*) "p9 =",p9
	    WRITE(66,*) "p10=",p10
	    WRITE(66,*) "p11=",p11
	    WRITE(66,*) "a0 =",a0
	    WRITE(66,*) "a1 =",a1
	    WRITE(66,*) "a2 =",a2
	    WRITE(66,*) "a3 =",a3
	    WRITE(66,*) "a4 =",a4
	    WRITE(66,*) "a5 =",a5
	    WRITE(66,*) "d0 =",d0
	    WRITE(66,*) "d1 =",d1
	    WRITE(66,*) "d2 =",d2
	    WRITE(66,*) "d3 =",d3
	    WRITE(66,*) "d4 =",d4
	    WRITE(66,*) "d5 =",d5
	    WRITE(66,*) "d6 =",d6
	    WRITE(66,*) "d7 =",d7
	    WRITE(66,*) "d8 =",d8
	    WRITE(66,*) "d9 =",d9
	    WRITE(66,*) "d10=",d10
	    WRITE(66,*) "d11=",d11
	    WRITE(66,*) "d12=",d12
	    WRITE(66,*) "d13=",d13
	    WRITE(66,*) "d14=",d14
	    WRITE(66,*) "d15=",d15
            CLOSE(66)
	  END IF
!*****C
       endif
	END IF
	ifrecall = 1
!******

!980610 Modification of section of "sx" containing the den and nums of the "S"'s

!###############################################

         D = d0 + d1*y*n**2 + d2*y*n*c + d3*y*c**2 + d4*n**3 + d5*n**2*c &
      + d6*n*c**2 + d7*c**3 &
      + d8*y*n + d9*y*c + d10*n**2 + d11*n*c + d12*c**2 + d13*y &
      + d14*n + d15*c

!########################################################################
         Nm = a0 + a1*n**2 + a2*n*c + a3*c**2 + a4*n + a5*c

!###########################################################################


         Nh  = - (30.D0*n*p6 + 30.D0*c*p4 - 60.D0 &
         - ( 2.D0*p1m  + 15.D0*p2m**2 - 6.D0*p2m  - 5.D0*p1m**2 ) &
         * y ) &
         * (c*p4*p7 - c*p11 - n*p11 + 1.D0)


         Nc  =   (30.D0*n*p6 + 30.D0*c*p4 - 60.D0 &
         - ( 2.D0*p1m  + 15.D0*p2m**2 - 6.D0*p2m  - 5.D0*p1m**2 ) &
         * y ) &
         * (c*p10 - 1.D0 - n*p6*p7 + n*p9)


!980610-15 Modification of section of "sx" containing Sm, Sh, Sc

!*******************************************************************************
         Sm = (4.D0/15.D0) * tpvot * Nm/D

         Sh = (4.D0/15.D0) * tptot * Nh/D

         Sc = (4.D0/15.D0) * tpcot * Nc/D
!*******************************************************************************


      return
 1004 format(12(1pe14.5))
      end




                                                                        
      subroutine rtwi(x,val,fct,xst,eps,iend,ier,ttot,tcot,tctot,tptot,tpcot,tpvot,ri,rit1,ric1,rit,ric)
!     to solve general nonlinear equations of the form x=fct(x)       
!     by means of wegsteins iteration method                         
!     prepare iteration                                             
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
      real(r8) :: ttot,tcot,tctot,tptot,tpcot,tpvot
      real(r8) :: ri,rit1,ric1
!      implicit real*8 (a-h,o-z)
      ier=0                                                        
      tol=xst                                                     
      x=fct(tol,ttot,tcot,tctot,tptot,tpcot,tpvot,ri,rit1,ric1,rit,ric)
      a=x-xst                                                   
      b=-a                                                     
      tol=x                                                   
      val=x-fct(tol,ttot,tcot,tctot,tptot,tpcot,tpvot,ri,rit1,ric1,rit,ric)
!     start iteration loop                                 
      do 6 i=1,iend                                       
!981103 Crude fix to avoid mysterious problem which occurred with a close but not too close guess.
#if (defined D_PRECISION)
      IF(DABS(val).LT.1.D-12) val =0.D0
#else
      IF(ABS(val).LT.1.D-12) val =0.D0
#endif
      if(val) 1,7,1                                      
!     equation is not satisfied by x                    
 1    b=b/val-1.D0                                       
      if(b) 2,8,2                                     
!     iteration is possible                          
 2    a=a/b                                         
      x=x+a                                        
      b=val                                       
      tol=x                                      
      val=x-fct(tol,ttot,tcot,tctot,tptot,tpcot,tpvot,ri,rit1,ric1,rit,ric)
!     test on satisfactory accuracy            
      tol=eps                                 
      d=abs(x)                               
      if(d-1.d0) 4,4,3                        
 3    tol=tol*d                            
 4    if(abs(a)-tol) 5,5,6                
 5    if(abs(val)-10.D0*tol) 7,7,6         
 6    continue                          
!     end of iteration loop                                           
!     no convergence after iend iteration steps. error return.       
      ier=1                                          
 7    return                                        
!     error return in case of zero divisor         
 8    ier=2                                       
      return                                     
      end                                       


!030401-04Y Adaptation to the NCAR CSM Ocean Model of my calculation of the turbulence
!Y	functions of the one variable (N_d^2/N^2) for the zero shear unstable case
!Y      written for HYCOM [See NBp.030401-2to3]. I must make variables double precision
!Y	and remove the inclusion of the HYCOM file common_blocks_giss.h and instate the
!Y	common block bb0 which carries the timescale ratios [See NBp.030402-2].
!Y	[See NBp.030403-1.] Introduce a check on positivity of S_X in response to an 
!Y	error that I found in the calculation[See NBp.030404-5,6,8 .]
!030403Y Include "ttot'=`{\tau_\theta \over \tau}' in the common block with timescale 
!Y       ratios, /bb0/. See NBp030403-8to12.
!
!_______________________________________________________________________
