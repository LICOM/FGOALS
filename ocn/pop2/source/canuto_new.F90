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
      subroutine turb_2(z,t,s,rh,ri,rid,s2, &
                   fricmx,wndmix,v_back,t_back,s_back, &
                   an2,&					!030401 N^2
                	 ustar_,buoytur,buoysol,&		!020912 surface fluxes
                         Coriol,		&		!030424 f_coriolis
                         amld,akm,akh,aks,n,na,nmax,&
      			 isurfuse,			&	!020912 choice of use of surface fluxes
      		         ifextermld,ifoutput,ii,jj)
!******X
!**CyXI011107-09 ************INTRODUCE A ***NEW*** OPTION "ifastexpabs=1" TO GO WITH "ifexpabstable=1".************
!             ************THE CURRENT GENERIC "interp2d" STEPS ONE BY ONE THROUGH THE NONLINEAR************
!             ************PART OF MY LOOKUP TABLE IN "ri" and "rid". NOT READY NOW WITH A GENERIC FIX************
!             ************SO INSTEAD I TRY A QUICK FIX SPECIFIC TO A TABLE WITH AN EXPONENTIAL ABSOLUTE VALUE************
!             ************SPACING IN THE NONLINEAR PART. INVOKE A NEW ROUTINE "interp2d_expabs" TO DO THIS.************
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      use canuto_mod
      implicit real(r8) (a-h,o-z)
!lhl0606      implicit real*8 (a-h,o-z)
!-----------------------------------------------------------------------
!     level 2 turbulence models to be used in ogcm, 3-16-98
!     sm and sh depends only on the richardson number
!981015 In salinity model case level 2 means S_M,S_H,S_C depend only on Ri,Ri_d
!-----------------------------------------------------------------------
!     inputs:
!     
!     nmax       = total number of vertical grid points
!     n          = number of vertical grid points above sea bottom
!                  minus one
!020219D Need the level number of the bottom ocean level for enhanced bottom mixing schemes.
!     kbot       = number of vertical grid points above sea bottom (= n+1)
!*****CD
!980527 EXTEND TO CASES OF ONE AND ZERO OCEAN LEVEL. (n = 0,-1)
!     na         = Max(1,number of vertical grid points above sea bottom -1)
!     nb         = Max(0,number of vertical grid points above sea bottom -1)
!******
!020220D Increase number of vertical levels passed in by 1 because need to know 
!	 the bottom depth for purposes of enhanced bottom mixing option.
!     z(n+1)     = depth of vertical grid (cm, positive)
!*****CD
!     t(n)       = mean temperature (c)
!0501 s(n)       = mean salinity (parts per thousand minus 35)
!0501 rh(n)      = mean density (g cm**(-3))
!     ri(n)      = richardson number
!1015 rid(n)     = temperature component minus salinity component of Ri
!     s2(n)      = (dU/dz)**2+(dV/dz)**2, shear squared, (1/s)**2
!980501 Shear**2 units were mislabelled in turb_2, they should be frequency^2.
!     v_back(n)  = background turbulent viscosity (cm**2/sec)
!     t_back(n)  = background turbulent heat diffusivity (cm**2/sec)
!     s_back(n)  = background turbulent salt diffusivity (cm**2/sec)
!030401Y Add square of Brunt Vaisala frequency to the inputs.
!Y    an2(n)     = N^2 (1/s)**2
!981125-1216 s_back(n) becomes a quantity recalculated in the subroutine
!	when ifsalback>0.
!*****C
!991107 v_back,t_back,s_back become quantities calculated in the subroutine 
!       when ifback or ifsalback > 2
!020912-24,25X ******SURFACE FLUXES****** 
!     ustar_        = sfc friction velocity                    (cm/s)
!     buoytur       = sfc turbulent buoyancy forcing         (cm2/s3)
!     buoysol       = sfc radiative buoyancy forcing         (cm2/s3)
!X Switch for use of surface fluxes
!X    isurfuse  = 0 for don't use surface fluxes for dimensionalization of turbulence model. 
!X		    Rely on l_MLD^2 S, where l_MLD is the Blackadar lengthscale calculated using MLD.
!X              = 1 for use total surface buoyancy, buoytot=buoytur+buoysol, 
!X                  to calculate dissipation, epsilon. epsilon=-buoytot/(1-(1/2)yS_M). y=(\tau S)^2 .
!******X
!030424Z *****CORIOLIS PARAMETER*****
!     Coriol	 = 2 \Omega sin(latitude) 		     	(1/s)	
!******Z
!030424Z Choice of rotation influence on lengthscale.
!     ilomega    = 0 for l=l_Blackadar traditional Blackadar lengthscale without rotation
!
!                = 1 for l=(l_Blackadar^{-1} + l_\Omega^{-1})^{-1} ; when MLD>=500meters
!	                 l_\Omega \equiv \sqrt{-B*/f^3} for B*<0, and \infinity for B*>0 .
!******Z
!     
!
!     fricmx     = max viscosity (cm**2/sec)
!     wndmix     = min viscosity at bottom of 1st level to simulate
!                  missing high frequency windstress components
!                  (cm**2/sec)
!     visc_cbu_limit = largest "visc_cbu" in regions of gravitational
!                      instability (cm**2/sec)
!     diff_cbt_limit = largest "diff_cbt" in regions of gravitational
!                      instability (cm**2/sec)
!
!980501 For ifextermld = {0,1} calculate MLD {here, externally in calling unit).
!0501	ifextermld = MLD switch [0 or 1]
!980716 For ifoutput={0,1} don't do output fields to a file.
!******
!030717Z1a For icondear = 0 K_X = (1/2) B_1 {l_Deardorff}^2 S/(Sl/q)_Deardorff S_X_P=eps (traditional)
!	       icondear =-1 K_X = (1/2) B_1 {l_Blackadar}^2 S/{(Sl/q)_P=eps} S_X_P=eps (No Deardorff)
!	       icondear =+1 K_X = (1/2) B_1 {l_Deardorff}^2 S/{(Sl/q)_P=eps} S_X_P=eps (Deardorff length only) 
!*****CZ1a
!000215-29 For ifepson2=0  K_X = (1/2) B_1 l^2 S/(Sl/q) S_X
!  Background  ifepson2=1  K_X = (1/2) B_1^2 Ri (Sl/q)^2 (\epsilon/N^2) S_X
! Bg.&deep Fg. ifepson2=2  K_X = (1/2) B_1^2 Ri (Sl/q)^2 (\epsilon/N^2) S_X
!	Do *NOT* use Deardorff lengthscale limitation with (\epsilon/N^2) case.
!*****C
!030429Z1  For ifdeeplat=0	pelagic background (\epsilon/N^2) constant(for ifepson2>0)
!	       ifdeeplat=1      pelagic bg.(\epsilon/N^2)~L(latitude,N) fr. Gregg et al.03
!*****CZ1
!000316 FOR ifrafgmax = {0,1} {Don't,Do} limit background ra_r to be at maximum
!	foreground ra_r at \theta_r where turbulence exists as Ri=> +infinity.
!*****C
!990203 For ifsalback=3 case, v_back(n),t_back(n),s_back(n) become outputs.
!     outputs:
!
!     amld       = mixed layer depth (cm)
!     akm(nmax)  = turbulent viscosity (cm**2/sec)
!0501 akh(nmax)  = turbulent heat diffusivity (cm**2/sec)
!0501 aks(nmax)  = turbulent salt diffusivity (cm**2/sec)
!  
!     internal quantities:
!
!     nmodel     = 1: improved second order closure model
!                  2: parameter-free stochastic model
!     ntbl       = number of points in the look-up table
!     ria(ntbl)  = richardson number in the look-up table
!     slq2(ntbl) = shear number squared (s*l/q)**2 in the look-up table
!     sma(ntbl)  = sm in the look-up table
!     sha(ntbl)  = sh in the look-up table
!981015 Quantities for the temperature-salinity model. Enlarge the table.
!     rib(-mt:mt)      = Ri for the 2D look-up table with enlarged dimensions
!     ridb(-mt:mt)     = difference of Tem and Sal Ri's for the 2D look-up table
!     slq2b(-mt:mt,-mt:mt) = shear number squared (s*l/q)**2 in the 2D table
!     smb(-mt:mt,-mt:mt)   = sm in the look-up table
!     shb(-mt:mt,-mt:mt)   = sh in the look-up table
!     ssb(-mt:mt,-mt:mt)   = ss in the look-up table
!*****C
!     ri1        = Ri at the given level to be used in table interpolation
!     rid1       = Ri_d at the given level to be used in table interpolation
!     rimax      = maximum richardson number
!     visc_cbu_limit = largest "visc_cbu" in regions of gravitational
!                      instability (cm**2/sec)
!     diff_cbt_limit = largest "diff_cbt" in regions of gravitational
!                       instability (cm**2/sec)
!
!020912-24,26X ****Quantities for dimensionalization of foreground diffusivities using surface forcing.****
!X    buoytot	 = total buoyancy forcing into ocean (includes penetrative radiation) (cm**2/sec***3)
!X    epsy(na)       = dissipation, epsilon, times (\tau S)^2  (cm**2/sec***3)			
!X    lifupper    = *LOGICAL* variable true in upper ocean where can use l_MLD S^2 dimensionalization.
!X    lifepsy(na) = *LOGICAL* array for epsy dimensionaliation points (i.e. lifupper .TRUE. and epsy>=0).
!******X
!980501 Mixed Layer Depth definition
!0501 idefmld    = 0: MLD = 0.1 degrees centigrade off surface pot. temp. .
!	         = 1: MLD = deltemld deg. centigrade off surface pot. temp. .
!   	         = 2: MLD = delrhmld g cm^{-3} off surface pot. density .
!     deltemld   = pot. temperature difference criterion for idefmld = 1
!**** delrhmld   = pot. density difference criterion     for idefmld = 2
!980912 Switch for the 1 pt. closure model to allow use of the same set of model
!	constants which have been used in the atmosphere.
!	\gamma_{1 to 8} = [.96,.06,.16,.10,7.66,0.4,0,0.29] .
!	ifchengcon	= 0: old ocean constants,near-surface profile assumption
!	ifchengcon 	= 1: constants used in atmosphere,matched to experiments
!****
!*****
!991107 Switch for changing background diffusivities in old model with heat and
!	salt diffusivities equal
!	ifback=0   : Use input background values
!	ifback=4   : Background diffusivities from our model with S = N/sqrt(Ri)
!	             and fixed Ri = ri_internal
!	ifback=5   : Background diffusivities from our model with S = N/sqrt(Ri)
!		     and fixed Ri = backfrac*Ri_Critical
!*****C
!981006 SWITCH TO TURN ON SALINITY MODEL
!	ifsali	= 0: old model with heat and salt diffusivities equal
!	ifsali	= 1: Canuto's salinity model as worked out by the summer of '98
!*****
!040217Zi1b Minimum on shear to avoid singularity in ifzeroshear = .FALSE case.
!	ifshearmin = .FALSE.:leave out minimum (to allow my older runs to be reproduced).
!	ifshearmin = .TRUE. :minimum foreground shear of s2min to avoid shear=>0 problems.
!Zi1b
!030424Z SWITCH TO ENABLE ZERO SHEAR FORMULA (for Ri more negative than table)
!	ifzeroshear = .FALSE. : Old method: always use table of Ri and Ri_d.
!	ifzeroshear = .TRUE. : New method: use table of N_d^2/N^2 for Ri < Ri_tableminimum
!981125-0226 Switch for making salinity and heat backgrounds in ratio of S_S/S_H
!	using S_H and S_S from lowest level where they are nonzero for
!	levels where they are zero.(Leave them equal if S's zero at 1st level.)
!	ifsalback = 0: Use input salinity background values.
!	ifsalback = 1: Alter salinity background values to fit S_S/S_H ratio.
!	ifsalback = 2: As above but ratio taken at subcrit Ri at point's R_r.
!	ifsalback = 3: ALL backgrounds from our model using internal wave est.S.
!	ifsalback = 4: ALL bgs our model int. wave S=N/(Ri_i^(1/2)),Ri_i const..
!	ifsalback = 5: Like 4 but ra_r_i = constant factor * ra_r_crit.(theta_r)
!	ifsalback = 6: Like 4 but (S_M/(Sl/q))(ra_r_i)=constant*(S_M/(Sl/q))(0) 
!*****C
!990201-04 Constants used in ifsalback=3,4 cases.
!	back_l_0 = estimated lengthscale, l_0, for int. wave-generated turb.(cm)
!	back_k_0 = `minimum turbulent' wavenumber which yields back_l_0  (cm^-1)
!
!990204 Constants used in ifsalback=3 case.
!	back_s2  = estimated shear squared due to internal waves (sec^{-2})
!	back_sm2 = 1/back_s2 (sec^2)
!	v_back0 = residual background constant viscosity(may need for stability)
!	t_back0 = residual background constant heat diffusivity "             "
!	s_back0 = residual background constant salt diffusivity "             "
!*****C
!990202 Richardson numbers for the background turbulence in the ifsalback=3 case
!	back_ri1  = N^2/{S_{internal}}^2 = ({N_T}^2+{N_C}^2)/{S_{internal}^2
!	back_rid1 = {N_d}^2/(S_{internal})^2 = ({N_T}^2-{N_C}^2)/{S_{internal}^2
!	{S_internal}^2 \equiv back_s2
!*****C
!990204-0301 Ri and S^2 for the background turbulence in the ifsalback>=4 cases
!	ri_internal = Dubovikov model constant internal wave Richardson number 
!  	back_rid1 = {N_d}^2/(S_{internal})^2 = ({N_T}^2-{N_C}^2)/{S_{internal}^2
!	N^2 / {S_internal}^2  = ri_internal ; {S_internal}^2 = N^2 / ri_internal
!	back_rid1 = ri_internal {{N_T}^2 - {N_C}^2  \over N^2}, dividing by S^2,
!	back_rid1 = ri_internal {Ri_T - Ri_C \over Ri}
!	back_rid1 = ri_internal {Ri_d \over Ri}
!990205 (N^2 / S_internal^2) = ri_internal ;  S_internal^2 = (N^2 / ri_internal) . 
!	s2_back \equiv S_internal^2 (sec^(-2)).
!	s2_back   = N^2 / ri_internal
!990226-0303 ifsalback > 4 case :
!	back_ra_r = array for background ({Ri_T}^2 + {Ri_C}^2)^{1/2} at \theta_r
!000319 Add array for c_y at ra_r_max. [See NBp.000319-3,4 (Vol.IX)].
!	c_y_r0 = array for c_y at maximum ({Ri_T^2 + {Ri_C}^2)^{1/2} at \theta_r
!990315 Add arrays for dimensionless turbulence functions at background ra_r
!	sm_r1 = array for S_M at background ({Ri_T}^2 + {Ri_C}^2)^{1/2} at \theta_r
!	sh_r1 = array for S_H at background ({Ri_T}^2 + {Ri_C}^2)^{1/2} at \theta_r
!	ss_r1 = array for S_S at background ({Ri_T}^2 + {Ri_C}^2)^{1/2} at \theta_r
!000315 Add an array for (Sl/q)^2 at background ra_r also.
!	slq2_r1 = array for (Sl/q)^2 at background ({Ri_T}^2 + {Ri_C}^2)^{1/2} at \theta_r
!	Introduce  a switch for interpolation of 1D arrays versus theta_r 
!	for background dimensionless functions instead of 2D arrays vs. Ri,Ri_d.
!	ifbg_theta_interp=0: Interpolate 2D array with (Ri,Ri_d) indices.
!	ifbg_theta_interp=1: Interpolate 1D array with theta_r index.
!*****C
!*****C
!991107 ifback = 5 :
!	backfrac = constant fraction of Ri_critical used for background
!*****C
!	ifsalback =5 :
!	backfrac = constant fraction of ra_r_crit.(\theta_r) used for background
!	ifsalback =6 :
!	backfact = constant fraction of K_M(ra_r=0) which K_M background equals.
!*****C
!030429Z1 For case with ifdeeplat=1 interior value of (epson/N^2) is multiplied by
!Z1  a latitude and stratification dependent factor equal to unity when
!Z1  N =N_0= 5.24E-3 sec^{-1} and latitude = 30^o as per Gregg et al. Nature Vol.423.
!Z1  epson2__ is the constant =pelagic val. for ifdeeplat=0;30^o,N_0 val. for ifdeeplat=1. 
!Z1  epson2_ is the pelagic value = (\epsilon/N^2) for ifbotenhance=0.
!Z1  epson2 is the value of \epsilon/N^2 used. 
!020213-19D For case with both ifepson2>0 and ifbotenhance>0 background epsilon/(N^2) increased 
!	 from interior constant to higher values near the bottom.
!     epson2_    = non-enhanced backgound value for epson2 [epsilon/N^2)] (cm**2/sec).
!     eps_bot    = enhanced bottom dissipation of turbulent kinetic energy [epsilon] (cm**2 sec**(-3)).
!     epson2_bot = enhanced bottom epson2 [epsilon/N^2] (cm**2/sec). 
!*****CD
!*****CZ1
!000229 For case with "ifepson2=2", background AND foreground epsilon/(N^2)
!	taken as constant at levels below the highest in which foreground dies.
!000215	For case with "ifepson2=1", background epsilon/(N^2) taken as constant.
!     epson2     = (dissipation of turbulent kinetic energy)/(N^2), (cm**2/sec)
!****C
!     internal subroutines and functions:
!
!     formld,formld_te,formld_rh
!     for nmodel=1: ccoeff, ourl2 
!     for nmodel=2: mcoeff, mikel2, mksmsh, fct, rtwi
!
!     oursal2,mikesal2a,interp2d
!-----------------------------------------------------------------------
!******
      real(r8)::  ttot,tcot,tctot,tptot,tpcot,tpvot
!X       and logical array lifepsy for points in column where use epsilon y dimensionalization.
      LOGICAL lifupper,lifepsy(na)
!030424-25Z **Introduce an integer parameter for the effect of rotation on the lengthscale.**
!25Z        Introduce the complex variable zlomega for \sqrt{-B*/f^3} for diagnosis.
      COMPLEX*16 zlomega
!25Z        amldminlom is the minimum MLD for use of the lengthscale \sqrt{-B*/f^3}.
      DIMENSION aldeep(nbig)
!*****C
      real(r8) :: g_tur(8), d_tur(0:5), s_tur(0:9)

!030404Y Introduce array for an2.
      dimension z(na+1),t(na),s(na),rh(na),ri(na),rid(na),s2(na), &
          v_back(nmax),t_back(nmax),s_back(nmax) &
         ,an2(na)
      dimension akm(nmax),akh(nmax),aks(nmax)
!020923(First Day of Autumn)X Introduce an array for [\epsilon y], that is, dissipation times {\tau Shear}^2 .
      DIMENSION epsy(na)
!

!YU Jan. 20th, 2011
!000323 Raise flag and stop if number of levels exceeds nbig.
      IF(nmax.GT.nbig) THEN
       if (mytid.eq.0) then
	WRITE(*,*) " "
	WRITE(*,*) " "
	WRITE(*,*) "****************************"
	WRITE(*,*) "*PROBLEM IN TURB_2 ROUTINE.*"
	WRITE(*,*) "Number of model levels exceeds nbig."
	WRITE(*,*) "nbig=",nbig,"	nmax=",nmax
	WRITE(*,*) "If want to use this many levels"// &
             " increase nbig to be bigger than nmax."
	WRITE(*,*) "Program is stopping now."
       endif
        STOP
      END IF
!*****C
      nb=MAX(0,n)
!******
      kbot=n+1
!*****CD
!020912-17X ** Surface Buoyancy Flux (can be used for dimensionalization of turbulence model) ** 
!X	**Total Buoyancy Flux = Sum of "turbulent" and "solar" contributions**
	buoytot = buoytur + buoysol
!******X
!980501-0716 Choose the definition of MLD. 
!	     If ifextermld=1, keep the one that was defined externally.
      IF(ifextermld.EQ.0) THEN
!980527,050208 Use Mixed-Layer Routine only when there are at least two levels of sea.
        IF(n.GT.0) THEN
  	  IF(idefmld.EQ.0) THEN
            call formld(z,t,amld,n)
	  ELSE IF(idefmld.EQ.1) THEN
            call formld_te(z,t,deltemld,amld,n)
	  ELSE IF(idefmld.EQ.2) THEN
            call formld_rh(z,rh,delrhmld,amld,n)
	  END IF
        ELSE IF(n.EQ.0) THEN
  	  amld = z(1)
        ELSE 
	  amld = 0.D0
        END IF
      END IF
!
      al0=0.17*amld
!980717  Write internal turbulence quantities to fort.91 when writing enabled.
!	 Headers for each outputstep.
!981104 Add S_M,H,S to outputs.
!020924X Add (epsilon (\tau S}^2) to outputs in the case where it is calculated.
       if (mytid.eq.0) then
         IF(ifoutput.EQ.1) THEN
	   WRITE(91,*) " "
	   IF(isurfuse.EQ.0) THEN
	     WRITE(91,*) "z          al         slq2       "// &
                   "ri1        rid1       "// &
                   "sm         sh         ss         "// &
                   "v_back     t_back     s_back     "
	   ELSE IF(isurfuse.EQ.1) THEN
	     WRITE(91,*) "z          al         slq2       "// &
                   "ri1        rid1       "// &
                   "sm         sh         ss         "// &
                   "epsy       "// &
                   "v_back     t_back     s_back     "
	   END IF
	 END IF
       endif
!******
!  START OF FIRST LOOP THROUGH LEVELS
      IF(ifepson2.EQ.2) THEN
!000302 Initialize switch for sub(background-only) depth. 
        ifbelow=0      
      END IF
!020219D REFERENCE NOTE: START OF PRIMARY LOOP OVER DEPTH LEVELS.
      do 22 k=1,n
!030504-0803Zi1a Section for Gregg et al. parameterization case.
	  IF((ifepson2.EQ.2).AND.(ifdeeplat.GT.0)) THEN
!030502Z1 Initialize switch for reversion to deep lengthscale for (N/f)<1 for foreground
!Z1	  in ifepson2=2,ifdeeplat=1 case because Gregg et al. formula is then unusable.
	    ifnofsmall=0
!030504-0803Zi1a Calculate N from N^2 and set ifnofsmall flag when (N/f)<1
!Zi1a		 Note that when Gregg et al. use "f", it should really be interpreted as '|f|'.
!Zi1a	         See NBp.030803-3,5 .
	    IF(an2(k).GE.0.D0) an = SQRT(an2(k))
	    IF((an/ABS(Coriol)).LT.1.D0) THEN	!arccosh(N/f) is undefined can't use Gregg et al.
	      ifnofsmall=1
	    END IF
	  END IF
!*****CZi1a
          ri1=ri(k)
	  IF(ifsali.EQ.0) THEN
          if(ri1.ge.rimax) then
            ri(k)=rimax    ! ad hoc
!991107 Traditional received background case for Temperature model.
	    IF(ifback.EQ.0) THEN
              akm(k)=v_back(k)
              akh(k)=t_back(k)
!981102 Set background salt=heat diffusivity.
	      aks(k)=akh(k)
              goto 22
!991109C Set turbulence functions to values at rimax for ifback>0.
	    ELSE
             sm   = sma(ntbl)
	     sh   = sha(ntbl)
	     slq2 = slq2a(ntbl)	
	     ss = sh
!******C
	    END IF
!*****C
          elseif(ri1.le.ri0) then !asymptotically
            sm=sma(1)
            sh=sha(1)
!981102 Set background S_Salt = S_Heat
	    ss=sh
            slq2=slq2a(1)*ria(1)/ri1
          else ! linearly interpolate the look-up tables
            m=int((ri1-ri0)/dri)+1
            tmp=(ri1-ria(m))/(ria(m+1)-ria(m))
            sm=sma(m)+(sma(m+1)-sma(m))*tmp
            sh=sha(m)+(sha(m+1)-sha(m))*tmp
            slq2=slq2a(m)+(slq2a(m+1)-slq2a(m))*tmp
	    ss=sh
          endif
	  ELSE IF(ifsali.EQ.1) THEN
!981104 Use Ri_d = Ri_C - Ri_T in salinity-temperature turbulence model.
          rid1=rid(k)
!*****C

!030424Z ONLY ENABLE ZERO SHEAR APPROXIMATION WHEN "ifzeroshear" TRUE.[See NBp.030424-12.]
!030401Y Code for Ri => - \infinity approximation taken from my HYCOM turb_2_fixed2.fs0
!Y	 [See NBp.030401-1to4.]For HYCOM N^2 was reconstructed now it's an input argument.
!
!030324-28 For unstable case with shear too small for 2D table of (Ri,Ri_d) to handle,
!       treat as if Ri were minus infinity, the unstable zero shear case, and
!       interpolate 1D table of (N_d^2/N^2): N^2 = N_H^2 + N_S^2 , N_d^2 = N_H^2 - N_S^2 .
!       Estimate N^2 and N_d^2 from the Ri and Ri_d . Note that the minimum, "epsi",
!       on shear in the calling routine can cause Ri and Ri_d to be underestimated.
!       (N^d)^2/N^2 is a more accurate quantity.
        and2 = (rid(k)/ri(k))*an2(k)
        and2on2 = and2/an2(k)
!       Consider shear as being too small when N^2 < Ri_table_minimum * S^2 .
!       Note N^2/S^2, which can go to infinity, is the real Ri, as opposed to the
!       sometimes smaller N^2/(MAX(S^2,epsil)) used in the model which is always finite.
        IF((ifzeroshear).AND.(an2(k).LT.rib(-mt)*s2(k))) &
        THEN
          ifpureshear=1
        ELSE
          ifpureshear=0
        END IF
!******Z
        IF(ifpureshear.EQ.1) THEN
!030326 Introduce a modular 1D table interpolation derived
!       by stripping down the 2D table interpolation.
!011107yXI ***Option to call instead of generic interpolation one tailored to exponential absolute nonlinear part.***
!       imax is the maximum positive value of the table index
!       (equal to +mt when there is no realizability limit).
        imax = mt
            IF(ifastexpabs.EQ.0) THEN
              CALL INTERP1D(and2on2, &
                      and2on2a1,amtaun2a1,sma1,sha1,ssa1, &
                      amtaun2,sm,sh,ss, &
                      imax,mt,mt0,dand2on2)
            ELSE IF(ifastexpabs.EQ.1) THEN
              CALL INTERP1D_EXPABS(and2on2, &
                      and2on2a1,amtaun2a1,sma1,sha1,ssa1, &
                      amtaun2,sm,sh,ss, &
                      imax,mt,mt0,dand2on2,rnd2on2)
            END IF
!      write(6,9928) mt,mt0,and2on2,slq2,sm,sh,ss,dand2on2,rnd2on2
 9928  format('interp1dtest',2i5,1p,7e12.4)

!yXI

!       Reconstructed slq2 for the outputs - somewhat bogus.
          slq2 = (-amtaun2)/((b1**2)*ri(k))
          GO TO 5               !Skip 2D interpolation.
        END IF
!******
!******Y

!981015   Interpolate 2D table for salinity-temperature model case.
!011107yXI ***Option to call instead of generic interpolation one tailored to exponential absolute nonlinear part.***
            IF(ifastexpabs.EQ.0) THEN
              CALL INTERP2D(ri1,rid1, &
                      rib,ridb,slq2b,smb,shb,ssb, &
                      slq2,sm,sh,ss, &
                      irimax,mt,mt0,dri)
            ELSE IF(ifastexpabs.EQ.1) THEN
              CALL INTERP2D_EXPABS(ri1,rid1, slq2,sm,sh,ss)
            END IF
!yXI
	  END IF
!981216-990108 Check that "slq2" has been set to 0 where it might have been negative.
	IF(slq2.LT.0.D0) THEN
       if (mytid.eq.0) then
	  WRITE(*,*) "************************************************"
	  WRITE(*,*) "Error detected in turbulence module." 
	  WRITE(*,*) "'slq2' negative in turb_2 subroutine" &
               //" after interpolation."
	  WRITE(*,*) "k=",k,"     slq2=",slq2
	  WRITE(*,*) "sm=",sm,"   sh=",sh,"   ss=",ss
	  WRITE(*,*) "ri1=",ri1,"    rid1=",rid1
	  WRITE(*,*) "dri=",dri
	  WRITE(*,*) "Program will stop."
	  WRITE(*,*) "************************************************"
       endif
	  STOP
	END IF
!*****C
!000302 Assume region contiguous with surface where foreground model is
!	realizable has ended when get 0 "slq2".
	IF(slq2.EQ.0.D0) ifbelow = 1
!*****C
!030401 Skipped from 1D table interpolation to here for unstable zero shear approximation.
    5   CONTINUE

!020912-24,25X ******Calculate epsy(k) \equiv epsilon * (\tau S)^2 for surface forcing dimensionalization.******
!X	       Take epsilon = -buoytot/(1 - ((1/2)(\tau S)^2)S_M))
!X	       epsy = -buoytot/((b1^2 (Sl/q)^2)^{-1} - (1/2)S_M) [See NBp020912-8,12 
!X	       AND NBp020925-2{extension on cover}.(error of (1/4) for (1/2) had to be corrected.) ]
!020919X Keep track of total number of points to compare with number of problem points.
!020923X Define problem points as points where epsy<0 *AND* l^2 S dimensionalization is used 
!X       in isurfuse=0 case AND l is based on MLD. 
!020923X Define negative problem points: problem points with Ri<0 ,
!X	 Keep track of number of points where l^2 S dimensionalization is used AND l is based on MLD.
!020924,26X,030504Z1 Introduce lifupper, logical which is 1 only where l_MLD^2 S would be used for isurfuse=0 .
!26X-030803Zi1a      set lifupper false when unrealizable
	lifupper= (((ifepson2.LT.2).OR.(ifbelow.EQ.0)   &
              .OR.((ifdeeplat.GT.0).AND.(ifnofsmall.EQ.1)).OR.&	!030504 Revert to l^2 S for N/f<1 in Gregg et. al. case.  
            ((ri1.LT.0.D0).AND.(k.LE.2)))   &
            .AND.(slq2.GT.0.D0)) 
	IF(lifupper) THEN
          ilmldpoint=ilmldpoint+1
	  IF(ri1.LT.0.D0) ilmldpointneg=ilmldpointneg+1
	END IF
	ipoint=ipoint+1
        IF(k.EQ.1) icall=icall+1
!020924,26X Initialize epsy and lifepsy for this point.
	epsy(k)=0.D0
	lifepsy(k)=.FALSE.
	IF((isurfuse.EQ.1).AND.(n.GT.0)) THEN		
!020924X Although it seems to be being calculated to zero anyway in unrealizable Ri cases,
!X        set epsy to zero  when slq2 is zero for safety's sake.
	  IF(slq2.EQ.0.D0) THEN
	    epsy(k)=0.D0
	  ELSE
	    epsy(k) = -buoytot/((1.D0/((b1**2)*(slq2))) - 0.5D0*sm)
	  END IF
!020924,26X Introduce lifepsy(k), which is 1 only when epsy dimensionalization used.
	  lifepsy(k)= ((epsy(k).GE.0.D0).AND.lifupper)
!020924X Comment out write outs of epsy diagnosis.
!C020918-26X Write negative (epsilon y)'s to fort.67  and all values to fort.68
 	  IF((epsy(k).LT.0.D0).AND.lifupper) THEN
 	    iproblem=iproblem+1
	    IF(ri1.LT.0.D0) inegproblem=inegproblem+1
! 	    IF(k.EQ.1) THEN
!	      WRITE(67,*) "icall=",icall,
!    &         "  iproblem=",iproblem," inegproblem=",inegproblem
!              WRITE(67,*) "  ipoint=",ipoint,
!    &         " ilmldpoint=",ilmldpoint," ilmldpointneg=",ilmldpointneg
!	    END IF
! 	    WRITE(67,*) "epsy=",epsy(k)
! 	    WRITE(67,*) " k=",k
! 	    WRITE(67,*) "t=",t(k)," s=",s(k)
! 	    WRITE(67,*) "ri1=",ri1," rid1=",rid1
!	    WRITE(67,*) "slq2=",slq2
!	    WRITE(67,*) "sm=",sm
!	    WRITE(67,*) "buoytur=",buoytur," buoysol=",buoysol,
!    &                 "	buoytot=",buoytot
!	    WRITE(67,*) " "
	  END IF
! 	  IF(k.EQ.1) THEN
!	    WRITE(68,*) "icall=",icall,
!    &       "  iproblem=",iproblem," inegproblem=",inegproblem
!           WRITE(68,*) "  ipoint=",ipoint,
!    &       " ilmldpoint=",ilmldpoint," ilmldpointneg=",ilmldpointneg
!	  END IF
!	  WRITE(*,*) "icall=",icall,
!    &     "  iproblem=",iproblem," inegproblem=",inegproblem
!         WRITE(*,*) "  ipoint=",ipoint,
!    &     " ilmldpoint=",ilmldpoint," ilmldpointneg=",ilmldpointneg
! 	  WRITE(68,*) "epsy=",epsy(k)
! 	  WRITE(68,*) " k=",k
! 	  WRITE(68,*) "t=",t(k)," s=",s(k)
! 	  WRITE(68,*) "ri1=",ri1," rid1=",rid1
!	  WRITE(68,*) "slq2=",slq2
!	  WRITE(68,*) "sm=",sm	
!	  WRITE(68,*) "buoytur=",buoytur," buoysol=",buoysol,
!    &               "	buoytot=",buoytot
!	  WRITE(68,*) " "
	END IF
!******X
          akz=0.4*z(k)
          al=akz*al0/(al0+akz)
!030425Z **MODIFICATION OF LENGTHSCALE BY ROTATION OPTION**
!Z	 Use l=(l_Blackadar^{-1} + l_\Omega^{-1})^{-1} with l_\Omega \equiv \sqrt{-B*/f^3}
!Z	 *when* B*<0 and MLD extends to a set minimum [See NBp.030424-5to8&13&14,25-2&3.].
	  IF(ilomega.EQ.1) THEN
	    IF((buoytot.LT.0.D0).AND.(amld.GE.amldminlom)) THEN
	      rlomega = SQRT((Coriol**3)/(-buoytot))
	    ELSE
	      rlomega = 0.D0
	    END IF
	    rlblackadar = 1.D0/al  
	    rl = rlblackadar + rlomega
	    al = 1.D0/rl
	  END IF
!******Z
          al2=al*al
!000302 Do not use Deardorff limitation when use (\epsilon/N^2) dimensionalization.
!020924,26X FOR CONSISTENCY OF OUTPUTS WITH SURFACE BUOYANCY FORCING EPSILON DIMENSIONALIZATION, 
!X       ALTHOUGH "slq2" IS NOT USED, DO NOT APPLY DEARDORFF LIMITATION ON "(Sl/q)^2" .
	IF(.NOT.(((ifepson2.EQ.2).AND.(ifbelow.EQ.1)).OR.lifepsy(k))) &
   THEN
!         length scale reduction by buoyancy
!030717Z1a,050208 ***REMOVE DEARDORFF LIMITATION ON LENGTHSCALE IN "icondear=-1" CASE,***
!	   ***APPLY DEARDORFF LIMITATION ONLY TO LENGTHSCALE *NOT* TO {\tau N} IN "icondear=+1" CASE.***
!	   ***TRADITIONAL METHOD OF LIMITING {\tau N} AFTER S_X SET IS INCONSISTENT, "icondear=0" CASE.***
!	   DEARDORFF LIMITED {\tau N} BUT THIS IS NOT REALLY CONSISTENT WITH PRODUCTION=DISSIPATION.
          if(ri1.gt.0.D0) then
            anlq2=slq2*ri1
            if((anlq2.gt.0.281D0).AND.(icondear.GE.0)) then  !0.281=0.53**2
              al2=0.281D0/anlq2*al2
              IF(icondear.EQ.0) slq2=0.281D0/(ri1+1.D-20)
            endif
          endif
!*****CZ1a
        END IF
!*****C
!020219-20D Calculation of epsilon must be placed outside ifsali blocks (see NBp020219-1,2 Vol.XIII)
!030429-0504Z1 Calculate epsilon2_, the pelagic value at this latitude and stratification,
!Z1	  equal to epsilon2__ in the ifdeeplat=0 case, 
!Z1	  but in the ifdeeplat=1 case equal to epsilon2__ multiplied by L(\theta,N)
!Z1	  from equation (2) appearing in  "Reduced mixing from the breaking of 
!Z1	  internal waves in equatorial waters" by M.Gregg,T.Sanford & D. Winkel in
!Z1	  Nature /Vol.422, 3 April 2003 p.513-515, which cites as the source for
!Z1	  this formula "Energy and action flow through the internal wave field:
!Z1	  an eikonal approach." by Henyey,Wright&Flatte in JGR Vol.94 , 1989, 9686-9698.
!Z1	  L is one at latitude thirty degrees and Brunt Vaisala frequency N_0=5.24e-3.
!0501Z	  *NOTE IN GREGG ET AL.'S FORMULA "f" IS *THE ABSOLUTE VALUE* OF `f_Coriolis'.*
!Z1	  Because Gregg et al. who confirm the formula observationally, admit that there
!Z1	  is some wave-mixing at the equator where this formula would set it to zero,
!Z1	  I ARTIFICIALLY SET A MINIMUM ON L TO REPRESENT AS YET UNSTUDIED PROCESSES
!Z1	  WHICH LEAD TO EQUATORIAL BACKGROUND MIXING SEVERAL TIMES MOLECULAR.
!Z1	  (From Figure 1 of Gregg et al. I think the minimum must be greater than 0.02.)
!0502Z1   *SET L TO ZERO WHEN (N/f) BECOMES LESS THAN 1 BECAUSE IT IS THEN UNDEFINED.*
!0504Z1	   SET "ifnofsmall=1" AT START OF LOOP IN THIS CASE TO REVERT TO DEEP LENGTHSCALE.
	  IF(an2(k).LT.0.D0) THEN
	    epson2_ = epson2__	 	!Just to "hold the place". Value irrelevent here.
	  ELSE 
	    IF(ifdeeplat.EQ.0) THEN
	      epson2_ = epson2__
	    ELSE IF((ifdeeplat.EQ.1).OR.(ifdeeplat.EQ.2)) THEN !030803Zi1a (See NBp.030803-4,5.)
!0504Z1  Since (N/f)<1 => arccosh(N/f) undefined, can't use Gregg et al. formula.
	      IF(ifnofsmall.EQ.1) THEN	
	        eplatidepend = 0.D0
	      ELSE
	        eplatidepend = EPLATIDEPEND_(ABS(Coriol),an)
	      END IF
	        eplatidepend = & 
          MAX(eplatidepend,eplatidependmin)
	        epson2_ = epson2__*eplatidepend
	    END IF
	  END IF
!*****CZ1
	 IF(ifepson2.GE.1) THEN
!020214D Option for enhanced mixing near the bottom. 
!	  epson2_ is the "pelagic value" of epsilon/N^2 , that in the interior of the ocean.
!	  eps_bot is the enhanced bottom TKE dissipation, epsilon.
!	  epson2_bot = eps_bot / N^2 
	   IF(ifbotenhance.EQ.0) THEN 
	     epson2 = epson2_ 
	   ELSE IF(ifbotenhance.EQ.1) THEN
	     eps_bot = eps_bot0 * EXP((z(k) - z(kbot))/scale_bot)
	     epson2_bot = eps_bot/(ri(k)*s2(k))
	     epson2 = MAX(epson2_,epson2_bot)
	   END IF
         END IF
!*****CD
!*****CD
!************************************************************************
!************************************************************************
!991107 Change background diffusivities to value using ocean model N and
!	background internal wave model Ri for ifback >= 4 .
!************************************************************************
!	Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to 
!	diffusivities calculated using the turbulence model
!	with Ri and l_0 replaced by constants 'ri_internal' and 'back_l_0' 
!	and S^2 replaced by  (N^2 / Ri_internal) for N^2>=0 and 0 for N^2 <0
!	to represent a modified Dubovikov internal wave generated turbulence
!	with universal constant Richardson number for ifback=4 case.
	IF((ifback.GE.4).AND.(ifsali.EQ.0)) THEN
!990205 Use a constant background Ri estimate. 
	  IF(ifback.EQ.4) THEN
	    back_ri1  = ri_internal
	  ELSE
!************************************************************************
!      	Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to 
!	diffusivities calculated using the turbulence model
!	with l_0 replaced by a constant 'back_l_0' and 
!	Ri by a constant 
!	and S^2 replaced by  (N^2 / Ri_internal) for N^2>=0 and 0 for N^2 <0
!	to represent a modified Dubovikov internal wave generated turbulence
!	with constant Richardson number for ifback>4 case.
!
!  Calculate the background Ri .
	    back_ri1 = backfrac*rimax
!*****C	    
	  END IF
!       Use the calculated background Ri in the turbulence model.
!         linearly interpolate the look-up tables
          m=int((back_ri1-ri0)/dri)+1
          tmp=(back_ri1-ria(m))/(ria(m+1)-ria(m))
          sm_back=sma(m)+(sma(m+1)-sma(m))*tmp
          sh_back=sha(m)+(sha(m+1)-sha(m))*tmp
          slq2_back=slq2a(m)+(slq2a(m+1)-slq2a(m))*tmp
	  ss_back=sh_back
       if (mytid.eq.0) then
        write(*,*) 'i=',ii,'j=',jj,'m',m
        write(*,*)
        write(*,*) 'tmp',tmp
        write(*,*)
        write(*,*) 'sm_back', sm_back,'sma(m)',sma(m),'sma(m+1)',sma(m+1)
        write(*,*)
        write(*,*) 'sh_back', sh_back,'sha(m)',sha(m),'sha(m+1)',sha(m+1)
        write(*,*)
        write(*,*) 'slq2_back', slq2_back,'slq2a(m)',slq2a(m),'slq2a(m+1)',slq2a(m+1)
        write(*,*)
       endif
! Check that "slq2" has been set to 0 where it might have been negative.
          IF(slq2_back.LT.0) THEN
       if (mytid.eq.0) then
            WRITE(*,*) & 
      "************************************************"
            WRITE(*,*) "Error detected in turbulence module."
            WRITE(*,*) "'slq2_back' negative in turb_2 subroutine" &
                 //" after 1D interpolation of background."
            WRITE(*,*) "k=",k,"     slq2_back=",slq2_back
            WRITE(*,*) &
       "sm_back=",sm_back,"   sh_back=",sh_back
            WRITE(*,*) "back_ri1=",back_ri1
            WRITE(*,*) "dri=",dri
            WRITE(*,*) "Program will stop."
            WRITE(*,*) & 
      "************************************************"
       endif
            STOP
          END IF
!*****C
!990205 Calculate the square of the shear from the background Richardson number.
!	s2_back   = N^2 / ri_internal = (N^2 / S_ext^2) (S_ext^2 /ri_internal) 
!	          = (Ri_ext / ri_internal) S_ext^2
     	  s2_back = (ri1/back_ri1)*s2(k)
!990208-0301 Set square of shear to zero for unstable density stratification.
          IF(ri1.LE.0.D0) s2_back = 0.D0
!*****C
!990316 Set ill-defined S_M,H,S for unstable density stratification to zero.
	  IF(ri1.LT.0.D0) THEN
	    sm_back = 0.D0
	    sh_back = 0.D0
	    ss_back = 0.D0
	  END IF
!*****C
!000215 Skip background lengthscale calculation when using K_X/(epsilon/N^2) .
	IF(ifepson2.EQ.0) THEN
!990203,050208 Use the constant background l_0 lengthscale in the turbulence model.
          al0_back = back_l_0
          akz=0.4D0*z(k)
          al_back=akz*al0_back/(al0_back+akz)
          al2_back=al_back*al_back
!         length scale reduction by buoyancy
!030717Z1a ***REMOVE DEARDORFF LIMITATION ON LENGTHSCALE IN "icondear=-1" CASE,***
!	   ***APPLY DEARDORFF LIMITATION ONLY TO LENGTHSCALE *NOT* TO {\tau N} IN "icondear=+1" CASE.***
!	   ***TRADITIONAL METHOD OF LIMITING {\tau N} AFTER S_X SET IS INCONSISTENT, "icondear=0" CASE.***
!	   DEARDORFF LIMITED {\tau N} BUT THIS IS NOT REALLY CONSISTENT WITH PRODUCTION=DISSIPATION.
          if(back_ri1.gt.0.D0) then
            anlq2_back=slq2_back*back_ri1
            if((anlq2_back.gt.0.281D0).AND.(icondear.GE.0)) then  !0.281=0.53**2
              al2_back=0.281D0/anlq2_back*al2_back
              IF(icondear.EQ.0) slq2_back=0.281D0/(back_ri1+1.D-20)
            endif
          endif
!*****CZ1a
!990203-05,050208 Calculate the background diffusivities.
          tmp_back=0.5D0*b1*al2_back*sqrt(s2_back/(slq2_back+1.D-40))
!000215 Use K_X = K_X/(\epsilon/N^2) * (\epsilon/N^2)    
!		From NBp.000215-5, Volume IX : 
!        	 K_X/(\epsilon/N^2) = (1/2) B_1 Ri (S l/q)^2 S_X  .
!	    K_X = (((1/2) B_1^2 Ri (S l/q)^2)* (\epsilon/N^2)) * S_X 
	ELSE IF(ifepson2.GE.1) THEN
!040422Zi1bj-AH:MIT Dirty act as stopgap response to warning message produced by
!	  MIT ocean model. **Zero** al_back which has not been defined.
!	  [See NBp.040422-7&8.] Think should really calculated an equivalent
!	  lengthscale for the (l^2 S) background dimensionalization,
!	  but don't have time to work this out now. 
!	  al_back is *not* used to calculate the diffusivity in the 
!	  (\epsilon/N^2) dimensionalization case, but could be a diagnostic.
!	  Previously it had been printing out zero anyway because undefined.
	  al_back=0.D0
          tmp_back=0.5D0*b1**2*back_ri1*slq2_back*epson2
	END IF
          v_back(k)=tmp_back*sm_back
          t_back(k)=tmp_back*sh_back
          s_back(k)=tmp_back*ss_back
!************************************************************************
	END IF
!************************************************************************
!************************************************************************
!************************************************************************
!991108 BEGIN SECTION FOR SALINITY MODEL BACKGROUND DIFFUSIVITY CALCULATION.
	IF(ifsali.GT.0) THEN
!*****C
!************************************************************************
!************************************************************************
!981125 Change background salt diffusivity from input value to 
!	(S_S/S_H)*(background heat diffusivity) for ifsalback=1 case.
	IF(ifsalback.EQ.1) THEN
	  IF(slq2.NE.0.D0) THEN
	    sm_last = sm
	    sh_last = sh
	    ss_last = ss
	  ELSE IF(k.eq.1) THEN
!       Set salt background = heat background if turbulenceless 1st layer.
	    s_back(k) = t_back(k)
	    GO TO 20 
	  END IF
	  s_back(k) = (ss_last/sh_last)*t_back(k)
!************************************************************
!981216 Change background salt diffusivity to
!	(S_S/S_H)*(background heat diffusivity), 
!	but with S_S,H taken at Ri a little less than Ri_max
!	at the given Ri_T/Ri_C value for ifsalback=2 case.
!990412 When Ri_T = 0, 
!       correctly set the angle 'theta_r' in the (Ri_T,Ri_C) plane to 'pi'/2 , 
!	unless Ri_C = 0, in which case **arbitrarily** set 'theta_r' to 'pi'/4 .
	ELSE IF(ifsalback.EQ.2) THEN
	  IF(slq2.NE.0.D0) THEN
	    sisa = ss/sh
	  ELSE 
!  	Linearly interpolate sisamax array to the angle in (Ri_C,Ri_T) space .
!	Ri \equiv Ri_T + Ri_C 	; Ri_d \equiv Ri_T - Ri_C .
	    rit = (ri(k) + rid(k))/2.D0
	    ric = (ri(k) - rid(k))/2.D0
!990412-13 Find \theta_r for the Ri_T = 0 case. Treat "0/0 = 1".
            IF(rit.EQ.0.D0) THEN
	      IF(ric.EQ.0.D0) THEN
	        theta_r = ATAN(1.D0)
	      ELSE
	        theta_r = pi/2.D0	! Arctangent of infinity.
	      END IF
	    ELSE
	      theta_r = ATAN(ric/rit)
	    END IF
!*****C
!990111 Make sure the right choice of arctan(Ri_C/Ri_T) [\theta_r] is made.
!	Arctan only covers the range (-pi/2,pi/2) which theta_r may be outside.
!000323 Choose to have theta_r in range (-pi/4,3pi/4), stable case only.
	    IF(ABS(theta_r).GT.(pi/2.D0)) STOP
	    IF(theta_r.LT.(-pi)/4.D0) theta_r = theta_r + pi
!000309 MAKE 'jtheta' A NON-NEGATIVE INDEX - ZERO AT THETA = -PI/4 .
!	The fortran function "INT" rounds to the integer *NEAREST TO ZERO*
!	**I.E. ROUNDS **UP** FOR NEGATIVE NUMBERS**, DOWN ONLY FOR POSITIVES.
	    jtheta_r0 = INT((theta_r + (pi/4.D0))/deltheta_r)
!000309 INTRODUCE 'itheta' HERE FOR THE INDEX THAT IS ZERO AT THETA=0.
	    itheta_r0 = jtheta_r0 - n_theta_r_oct
	    itheta_r1 = itheta_r0+1
	    theta_r0 = itheta_r0*deltheta_r
	    theta_r1 = itheta_r1*deltheta_r
!000314 Angle in degrees.
	   theta_r_deg = theta_r*180.D0/pi
!*****C
!	Sound the alarm if have unrealizability outside expected range in angle.
	    IF((itheta_r1.GT.3*n_theta_r_oct).OR.  &
         (itheta_r0.LT.-n_theta_r_oct)) THEN
       if (mytid.eq.0) then
	      WRITE(*,*) &
         "************************************************"
	      WRITE(*,*) "Problem in turbulence module!"
	      WRITE(*,*) "Unrealizability outside Ri>0 region. "
	      WRITE(*,*) "slq2=",slq2,"    sm=",sm," sh=",sh," ss=",ss
	      WRITE(*,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
	      WRITE(*,*) "rit=",rit,"ric=",ric,"    theta_r=",theta_r
	      WRITE(*,*) "theta_r_deg=",theta_r_deg
	      WRITE(*,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
	      WRITE(*,*) "n_theta_r_oct=",n_theta_r_oct 
	      WRITE(*,*) " "
	      WRITE(*,*) "Program will stop."
	      WRITE(*,*) &
         "************************************************"
       endif
	      STOP
	    END IF
	    deltheta_r1 = theta_r - theta_r0
	    delsisa_r = sisamax(itheta_r1) - sisamax(itheta_r0)
	    dsisa_o_dtheta = delsisa_r/deltheta_r
	    sisa = sisamax(itheta_r0)+deltheta_r1*dsisa_o_dtheta
	  END IF
	  s_back(k) = sisa*t_back(k)
!************************************************************
!************************************************************************
!990202-03	Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to a
!	small constant plus diffusivities calculated using the turbulence model
!	with l_0 and S^2 replaced by constants intended to represent 
!	an internal-wave-generated turbulence for the ifsalback=3 case.
	ELSE IF(ifsalback.EQ.3) THEN
!990202 Use a constant background shear estimate to calculate background Ri,Ri_d
	  back_ri1  = ri1*s2(k)*back_sm2
	  back_rid1 = rid1*s2(k)*back_sm2
!990203 Use the calculated background Ri and Ri_d in the turbulence model.
!981015-990203   Interpolate 2D table for salinity-temperature model case.
!011107yXI ***Option to call instead of generic interpolation one tailored to exponential absolute nonlinear part.***
            IF(ifastexpabs.EQ.0) THEN
              CALL INTERP2D(back_ri1,back_rid1, &
                      rib,ridb,slq2b,smb,shb,ssb, &
                      slq2_back,sm_back,sh_back,ss_back, &
                      irimax,mt,mt0,dri)
            ELSE IF(ifastexpabs.EQ.1) THEN
              CALL INTERP2D_EXPABS(back_ri1,back_rid1, slq2_back,sm_back,sh_back,ss_back)
            END IF
!yXI
!981216-990108 Check that "slq2" has been set to 0 where it might have been negative.
        IF(slq2_back.LT.0) THEN
       if (mytid.eq.0) then
          WRITE(*,*) "************************************************"
          WRITE(*,*) "Error detected in turbulence module."
          WRITE(*,*) "'slq2_back' negative in turb_2 subroutine" &
               //" after interpolation of background."
          WRITE(*,*) "k=",k,"     slq2_back=",slq2_back
          WRITE(*,*) & 
          "sm_back=",sm_back,"   sh_back=",sh_back,"   ss_back=",ss_back
          WRITE(*,*) "back_ri1=",back_ri1,"   back_rid1=",back_rid1
          WRITE(*,*) "dri=",dri
          WRITE(*,*) "Program will stop."
          WRITE(*,*) "************************************************"
       endif
          STOP
        END IF
!*****C
!000215 Skip background lengthscale calculation when using K_X/(epsilon/N^2) .
	IF(ifepson2.EQ.0) THEN
!990203,050208 Use the constant background l_0 lengthscale in the turbulence model.
          al0_back = back_l_0
          akz=0.4D0*z(k)
          al_back=akz*al0_back/(al0_back+akz)
          al2_back=al_back*al_back
!         length scale reduction by buoyancy
!030717Z1a ***REMOVE DEARDORFF LIMITATION ON LENGTHSCALE IN "icondear=-1" CASE,***
!	   ***APPLY DEARDORFF LIMITATION ONLY TO LENGTHSCALE *NOT* TO {\tau N} IN "icondear=+1" CASE.***
!	   ***TRADITIONAL METHOD OF LIMITING {\tau N} AFTER S_X SET IS INCONSISTENT, "icondear=0" CASE.***
!	   DEARDORFF LIMITED {\tau N} BUT THIS IS NOT REALLY CONSISTENT WITH PRODUCTION=DISSIPATION.
          if(back_ri1.gt.0.D0) then
            anlq2_back=slq2_back*back_ri1
            if((anlq2_back.gt.0.281D0).AND.(icondear.GE.0)) then  !0.281=0.53**2
              al2_back=0.281D0/anlq2_back*al2_back
              IF(icondear.EQ.0) slq2_back=0.281D0/(back_ri1+1.D-20)
            endif
          endif
!*****CZ1a
!990203,050208 Calculate the background diffusivities.
          tmp_back=0.5D0*b1*al2_back*sqrt(back_s2/(slq2_back+1.D-40))
!000215 Use K_X = K_X/(\epsilon/N^2) * (\epsilon/N^2)    
!		From NBp.000215-5, Volume IX : 
!        	 K_X/(\epsilon/N^2) = (1/2) B_1 Ri (S l/q)^2 S_X  .
!	    K_X = (((1/2) B_1^2 Ri (S l/q)^2)* (\epsilon/N^2)) * S_X 
	ELSE IF(ifepson2.GT.0) THEN
          tmp_back=0.5D0*b1**2*back_ri1*slq2_back*epson2
	END IF
          v_back(k)=tmp_back*sm_back+v_back0
          t_back(k)=tmp_back*sh_back+t_back0
          s_back(k)=tmp_back*ss_back+s_back0
!************************************************************************
!************************************************************************
!990205-08	Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to 
!	diffusivities calculated using the turbulence model
!	with Ri and l_0 replaced by constants 'ri_internal' and 'back_l_0' 
!	and S^2 replaced by  (N^2 / Ri_internal) for N^2>=0 and 0 for N^2 <0
!	to represent a modified Dubovikov internal wave generated turbulence
!	with constant Richardson number for ifsalback=4 case.
	ELSE IF(ifsalback.GE.4) THEN
!990205 Use a constant background Ri estimate. 
	IF(ifsalback.EQ.4) THEN
	  back_ri1  = ri_internal
	  back_rid1 = (rid1/ri1)*ri_internal
	ELSE
!************************************************************************
!990301	Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to 
!	diffusivities calculated using the turbulence model
!	with l_0 replaced by a constant 'back_l_0' and 
!	Ri by a function of Ri_d
!	and S^2 replaced by  (N^2 / Ri_internal) for N^2>=0 and 0 for N^2 <0
!	to represent a modified Dubovikov internal wave generated turbulence
!	with stability-ratio dependent Richardson number for ifsalback>4 case.
!
!990412 When Ri_T = 0 and Ri_C \ne 0,
!       correctly set the angle 'theta_r' in the (Ri_T,Ri_C) plane to 'pi'/2 . 
!*****C
!990413	Skip background ra_r calculation in unstable *OR NEUTRAL* case.
!000323 Set background ra_r arbitrarily to zero in these cases.
	    IF(ri(k).LE.0.D0) THEN
	      back_ra_r1 = 0.D0
	      back_rit1  = 0.D0
	      back_ric1  = 0.D0
	      back_ri1   = 0.D0
	      back_rid1  = 0.D0
	      GO TO 19
	    END IF
!*****C
!  	Linearly interpolate back_ra_r array to this angle in (Ri_C,Ri_T) space.
!	Ri \equiv Ri_T + Ri_C 	; Ri_d \equiv Ri_T - Ri_C .
	    rit = (ri(k) + rid(k))/2.D0
	    ric = (ri(k) - rid(k))/2.D0
	    ra_r = SQRT((rit**2) + (ric**2))
!030403	Use same zero temp. gradient treatment as for isalback=2.[See NBp.030403-13.]
!990412-13 Find \theta_r for the Ri_T = 0 case. Treat "0/0 = 1".
            IF(rit.EQ.0.D0) THEN
              IF(ric.EQ.0.D0) THEN
                theta_r = ATAN(1.D0)
              ELSE
                theta_r = pi/2.D0       ! Arctangent of infinity.
              END IF
            ELSE
              theta_r = ATAN(ric/rit)
            END IF
!*****C
	    theta_r = ATAN(ric/rit)
!990111 Make sure the right choice of arctan(Ri_C/Ri_T) [\theta_r] is made.
!	Arctan only covers the range (-pi/2,pi/2) which theta_r may be outside.
!000323 Want to consider statically stable case only: Ri > 0.
	    IF(ABS(theta_r).GT.(pi/2.D0)) STOP
	    IF(theta_r.LT.(-pi)/4.D0) theta_r = theta_r + pi
!000309 MAKE 'jtheta' A NON-NEGATIVE INDEX - ZERO AT THETA = -PI/4 .
!	The fortran function "INT" rounds to the integer *NEAREST TO ZERO*
!	**I.E. ROUNDS **UP** FOR NEGATIVE NUMBERS**, DOWN ONLY FOR POSITIVES.
	    jtheta_r0 = INT((theta_r + (pi/4.D0))/deltheta_r)
	    jtheta_r1 = jtheta_r0+1
!000309 INTRODUCE 'itheta' HERE FOR THE INDEX THAT IS ZERO AT THETA=0.
	    itheta_r0 = jtheta_r0 - n_theta_r_oct
	    itheta_r1 = itheta_r0+1
!000330   ***WHEN THE ANGLE IS BETWEEN THE ANGLE FOR REALIZABILITY AT INFINITY***
!	  ***AND THE LAST TABLE ANGLE BEFORE THAT CRITICAL ANGLE, *** 
!	  ***SET IT TO THE LAST TABLE ANGLE BEFORE THE CRITICAL ANGLE.****
	    theta_r0 = itheta_r0*deltheta_r
	    theta_r1 = itheta_r1*deltheta_r
	    IF((theta_r0.LE.theta_rcrp).AND.(theta_r.GT.theta_rcrp)) THEN
	      theta_r = theta_r1
	      theta_r0 = theta_r1
	      itheta_r0 = itheta_r1 
	      itheta_r1 = itheta_r1+1
	      theta_r1 = theta_r1 + deltheta_r
	    ELSE IF((theta_r1.GE.theta_rcrn).AND.  &
              (theta_r.LT.theta_rcrn)) THEN
	      theta_r = theta_r0
	      theta_r1 = theta_r0
	      itheta_r1 = itheta_r0 
	      itheta_r0 = itheta_r0-1
	      theta_r0 = theta_r0 - deltheta_r
	    END IF
!*****C
!000314 Angle in degrees.
	   theta_r_deg = theta_r*180.D0/pi
!*****C
!	Sound the alarm if have unrealizability outside expected range in angle.
	    IF((itheta_r1.GT.3*n_theta_r_oct).OR.  &
         (itheta_r0.LT.-n_theta_r_oct)) THEN
       if (mytid.eq.0) then
	         WRITE(*,*) & 
        "************************************************"
	      WRITE(*,*) "Problem in turbulence module!"
	      WRITE(*,*) "Unrealizability outside Ri>0 region. "
	      WRITE(*,*) "slq2=",slq2,"    sm=",sm," sh=",sh," ss=",ss
	      WRITE(*,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
	      WRITE(*,*) "rit=",rit,"ric=",ric,"    theta_r=",theta_r
	      WRITE(*,*) "theta_r_deg =",theta_r_deg
	      WRITE(*,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
	      WRITE(*,*) "n_theta_r_oct=",n_theta_r_oct 
	      WRITE(*,*) " "
	      WRITE(*,*) "Program will stop."
	      WRITE(*,*) & 
        "************************************************"
       endif
	      STOP
	    END IF
	    deltheta_r1 = theta_r - theta_r0
	    delback_ra_r = back_ra_r(itheta_r1) - back_ra_r(itheta_r0)
	    dback_ra_r_o_dtheta = delback_ra_r/deltheta_r
	    back_ra_r1 = back_ra_r(itheta_r0) + & 
                  deltheta_r1*dback_ra_r_o_dtheta
!000316-17 In case choose ifrafgmax=1, ra_r is at maximum the ForeGround ra_r
!	at the "strong" double diffusive \theta_r's 
!       where have turbulence as Ri+> infinity. 
         ifrafglt=0
	 IF(ifrafgmax.EQ.1) THEN
	   IF((theta_r.LE.theta_rcrp).OR.(theta_r.GE.theta_rcrn)) THEN
	     IF(back_ra_r1.GT.ra_r) THEN
	       ifrafglt=1
	       back_ra_r1=ra_r
             END IF
	   END IF
	 END IF
   18    CONTINUE 
	 IF(back_ra_r1.LT.0.D0) THEN
       if (mytid.eq.0) then
	    WRITE(*,*) & 
      "************************************************"
	    WRITE(*,*) "Problem in turbulence module!"
	    WRITE(*,*) "Negative bg ra_r \\equiv (Ri_T^2+Ri_C^2)^(1/2)"
	    WRITE(*,*) "back_ra_r1 =", back_ra_r1
	    WRITE(*,*) "theta_r =", theta_r
	    WRITE(*,*) " "
	    WRITE(*,*) "slq2=",slq2,"    sm=",sm," sh=",sh," ss=",ss
	    WRITE(*,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
	    WRITE(*,*) "rit=",rit,"ric=",ric
	    WRITE(*,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
	    WRITE(*,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
	    WRITE(*,*) "theta_r_deg =",theta_r_deg
	    WRITE(*,*) "n_theta_r_oct=",n_theta_r_oct 
	    WRITE(*,*) " "
	    WRITE(*,*) " "
	    WRITE(*,*) "Program will stop."
	    WRITE(*,*) &
       "************************************************"
       endif
	    STOP
	 END IF 
!990301-0323 Calculate the background Ri and Ri_d .
	  back_rit1 = COS(theta_r)*back_ra_r1
	  back_ric1 = SIN(theta_r)*back_ra_r1
	  back_ri1  = back_rit1 + back_ric1
	  back_rid1 = back_rit1 - back_ric1
!*****C	    
	END IF
!000315 CALCULATE THE BACKGROUND DIMENSIONLESS TURBULENCE FUNCTIONS
!	USING TABLE OF VALUES FOR BACKGROUND "\theta_r"'S 
!	FOR "ifbg_theta_interp"=1.
!000317 Can only use theta_r table when do *not* reduce ra_r_BackGround
!	to a smaller ra_r_ForeGround.
	IF((ifbg_theta_interp.EQ.0).OR.(ifrafglt.EQ.1)) THEN
!990203 Use the calculated background Ri and Ri_d in the turbulence model.
!981015-990203   Interpolate 2D table for salinity-temperature model case.
!011107yXI ***Option to call instead of generic interpolation one tailored to exponential absolute nonlinear part.***
            IF(ifastexpabs.EQ.0) THEN
              CALL INTERP2D(back_ri1,back_rid1, &
                      rib,ridb,slq2b,smb,shb,ssb, &
                      slq2_back,sm_back,sh_back,ss_back, &
                      irimax,mt,mt0,dri)
            ELSE IF(ifastexpabs.EQ.1) THEN
              CALL INTERP2D_EXPABS(back_ri1,back_rid1, slq2_back,sm_back,sh_back,ss_back)
            END IF
!yXI
        ELSE IF(ifbg_theta_interp.EQ.1) THEN
!000315 Interpolate 1D table of background vs. theta_r instead.
	  deltheta_r1 = theta_r - itheta_r0*deltheta_r
	  delsm_back = sm_r1(itheta_r1) - sm_r1(itheta_r0)
	  dsm_back_o_dtheta = delsm_back/deltheta_r
	  sm_back = sm_r1(itheta_r0) + & 
                  deltheta_r1*dsm_back_o_dtheta
	  delsh_back = sh_r1(itheta_r1) - sh_r1(itheta_r0)
	  dsh_back_o_dtheta = delsh_back/deltheta_r

	  sh_back = sh_r1(itheta_r0) + & 
                  deltheta_r1*dsh_back_o_dtheta
	  delss_back = ss_r1(itheta_r1) - ss_r1(itheta_r0)
	  dss_back_o_dtheta = delss_back/deltheta_r
	  ss_back = ss_r1(itheta_r0) + & 
                  deltheta_r1*dss_back_o_dtheta
	  delslq2_back = slq2_r1(itheta_r1) - slq2_r1(itheta_r0)
	  dslq2_back_o_dtheta = delslq2_back/deltheta_r
	  slq2_back = slq2_r1(itheta_r0) + &
                   deltheta_r1*dslq2_back_o_dtheta

	ELSE
       if (mytid.eq.0) then
	  WRITE(*,*) "Problem with choice of background interpolation."
	  WRITE(*,*) "ifbg_theta_interp=",ifbg_theta_interp
	  WRITE(*,*) "ifrafglt=",ifrafglt
	  WRITE(*,*) "Program is stopping."
       endif
	  STOP
	END IF
!*****C
!981216-990108 Check that "slq2" has been set to 0 where it might have been negative.
        IF(slq2_back.LT.0) THEN
       if (mytid.eq.0) then
          WRITE(*,*) "************************************************"
          WRITE(*,*) "Error detected in turbulence module."
          WRITE(*,*) "'slq2_back' negative in turb_2 subroutine" &
               //" after interpolation of background."
          WRITE(*,*) "k=",k,"     slq2_back=",slq2_back
          WRITE(*,*) &
          "sm_back=",sm_back,"   sh_back=",sh_back,"   ss_back=",ss_back
          WRITE(*,*) "back_ri1=",back_ri1,"   back_rid1=",back_rid1
          WRITE(*,*) "dri=",dri
          WRITE(*,*) "Program will stop."
          WRITE(*,*) "************************************************"
       endif
          STOP
        END IF
!*****C
!990205 Calculate the square of the shear from the background Richardson number.
!	s2_back   = N^2 / ri_internal = (N^2 / S_ext^2) (S_ext^2 /ri_internal) 
!	          = (Ri_ext / ri_internal) S_ext^2
     	  s2_back = (ri1/back_ri1)*s2(k)
!990208-0301 Set square of shear to zero for unstable density stratification.
   19     IF(ri1.LE.0.D0) s2_back = 0.D0
!*****C
!990316,050208 Set ill-defined S_M,H,S for unstable density stratification to zero.
	  IF(ri1.LT.0.D0) THEN
	    sm_back = 0.D0
	    sh_back = 0.D0
	    ss_back = 0.D0
	  END IF
!*****C
	  IF((sm_back.LT.0.D0).OR.  &
       (sh_back.LT.0.D0).OR.  &
       (ss_back.LT.0.D0)) THEN
       if (mytid.eq.0) then
	         WRITE(*,*) & 
        "************************************************"
	      WRITE(*,*) "Problem in turbulence module!"
	      WRITE(*,*) "Negative Structure Function in Background."
        write(*,*) 'i=',ii,'j=',jj
        write(*,*) 'temperature'
        write(*,'(5d15.5)') t
        write(*,*)
        write(*,*) 'salinity'
        write(*,'(5d15.5)') s
        write(*,*)
        write(*,*) 'density'
        write(*,'(5d15.5)') rh
        write(*,*)
        write(*,*) 'richardson'
        write(*,'(5d15.5)') ri
        write(*,*)
        write(*,*) 'rid'
        write(*,'(5d15.5)') rid
        write(*,*)
        write(*,*) 'shear'
        write(*,'(5d15.5)') s2
        write(*,*)
        write(*,*) 'BV frequency'
        write(*,'(5d15.5)') an2
        write(*,*)
        write(*,*) 'ustart=',ustar_
        write(*,*)
        write(*,*) 'buoytur=',buoytur
        write(*,*)
        write(*,*) 'buoysol=',buoysol
	      WRITE(*,*) "slq2_back=",slq2_back
	      WRITE(*,*) "sm_back=",sm_back, &
                   " sh_back=",sh_back, &
                   " ss_back=",ss_back
	      WRITE(*,*) " "
	      WRITE(*,*) "back_ra_r1 =", back_ra_r1
	      WRITE(*,*) "theta_r =", theta_r
	      WRITE(*,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
	      WRITE(*,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
	      WRITE(*,*) " "
	      WRITE(*,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
	      WRITE(*,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
	      WRITE(*,*) "theta_r_deg=",theta_r_deg
	      WRITE(*,*) "n_theta_r_oct=",n_theta_r_oct 
	      WRITE(*,*) " "
	      WRITE(*,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
	      WRITE(*,*) "rit=",rit,"ric=",ric
	      WRITE(*,*) " "
	      WRITE(*,*) " "
	      WRITE(*,*) "Program will stop."
	      WRITE(*,*) & 
        "************************************************"
         endif
	      STOP
	    END IF
!000215 Skip background lengthscale calculation when using K_X/(epsilon/N^2) .
	IF(ifepson2.EQ.0) THEN
!990203,050208 Use the constant background l_0 lengthscale in the turbulence model.
          al0_back = back_l_0
          akz=0.4D0*z(k)
          al_back=akz*al0_back/(al0_back+akz)
          al2_back=al_back*al_back
!         length scale reduction by buoyancy
!030717Z1a ***REMOVE DEARDORFF LIMITATION ON LENGTHSCALE IN "icondear=-1" CASE,***
!	   ***APPLY DEARDORFF LIMITATION ONLY TO LENGTHSCALE *NOT* TO {\tau N} IN "icondear=+1" CASE.***
!	   ***TRADITIONAL METHOD OF LIMITING {\tau N} AFTER S_X SET IS INCONSISTENT, "icondear=0" CASE.***
!	   DEARDORFF LIMITED {\tau N} BUT THIS IS NOT REALLY CONSISTENT WITH PRODUCTION=DISSIPATION.
          if(back_ri1.gt.0.D0) then
            anlq2_back=slq2_back*back_ri1
            if((anlq2_back.gt.0.281D0).AND.(icondear.GE.0)) then  !0.281=0.53**2
              al2_back=0.281D0/anlq2_back*al2_back
              IF(icondear.EQ.0) slq2_back=0.281D0/(back_ri1+1.D-20)
            endif
          endif
!*****CZ1a
!990203-05 Calculate the background diffusivities.
          tmp_back=0.5D0*b1*al2_back*sqrt(s2_back/(slq2_back+1.D-40))
!000215 Use K_X = K_X/(\epsilon/N^2) * (\epsilon/N^2)    
!		From NBp.000215-5, Volume IX : 
!        	 K_X/(\epsilon/N^2) = (1/2) B_1 Ri (S l/q)^2 S_X  .
!	    K_X = (((1/2) B_1^2 Ri (S l/q)^2)* (\epsilon/N^2)) * S_X 
	ELSE IF(ifepson2.GT.0) THEN
          tmp_back=0.5D0*b1**2*back_ri1*slq2_back*epson2
	END IF
          v_back(k)=tmp_back*sm_back
          t_back(k)=tmp_back*sh_back
          s_back(k)=tmp_back*ss_back
!000309 Stop if background diffusivities are negative.
	  IF((v_back(k).LT.0.D0).OR.  &
       (t_back(k).LT.0.D0).OR.  &
       (s_back(k).LT.0.D0)) THEN
       if (mytid.eq.0) then
	         WRITE(*,*) &
         "************************************************"
	      WRITE(*,*) "Problem in turbulence module!"
	      WRITE(*,*) "Negative Background Diffusivity."
        write(*,*) 'i=',ii,'j=',jj
        write(*,*) 'temperature'
        write(*,'(5d15.5)') t
        write(*,*)
        write(*,*) 'salinity'
        write(*,'(5d15.5)') s
        write(*,*)
        write(*,*) 'density'
        write(*,'(5d15.5)') rh
        write(*,*)
        write(*,*) 'richardson'
        write(*,'(5d15.5)') ri
        write(*,*)
        write(*,*) 'rid'
        write(*,'(5d15.5)') rid
        write(*,*)
        write(*,*) 'shear'
        write(*,'(5d15.5)') s2
        write(*,*)
        write(*,*) 'BV frequency'
        write(*,'(5d15.5)') an2
        write(*,*)
        write(*,*) 'ustart=',ustar_
        write(*,*)
        write(*,*) 'buoytur=',buoytur
        write(*,*)
        write(*,*) 'buoysol=',buoysol
	      WRITE(*,*) "v_back=",v_back, &
                   " t_back=",t_back, &
                   " s_back=",s_back
	      WRITE(*,*) " "
	      WRITE(*,*) "slq2_back=",slq2_back
	      WRITE(*,*) "sm_back=",sm_back, &
                   " sh_back=",sh_back, &
                   " ss_back=",ss_back
	      WRITE(*,*) " "
	      WRITE(*,*) "back_ra_r1 =", back_ra_r1
	      WRITE(*,*) "theta_r =", theta_r, &
                   "   theta_r_deg=",theta_r_deg
	      WRITE(*,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
	      WRITE(*,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
	      WRITE(*,*) " "
	      WRITE(*,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
	      WRITE(*,*) "rit=",rit,"ric=",ric
	      WRITE(*,*) " "
	      WRITE(*,*) " "
	      WRITE(*,*) "Program will stop."
	      WRITE(*,*) &
         "************************************************"
        endif
	      STOP
	    END IF
!000314 Stop if background diffusivities are zero at positive Ri.
!                   v_back(k)=0.01D0
!                   t_back(k)=0.01D0
!                   s_back(k)=0.01D0
	  IF((ri(k).GT.0.D0).AND.((v_back(k).EQ.0.D0).OR.  &
       (t_back(k).EQ.0.D0).OR.  &
       (s_back(k).EQ.0.D0))) THEN
       if (mytid.eq.0) then
	         WRITE(*,*) & 
        "************************************************"
	      WRITE(*,*) "Problem in turbulence module!"
	      WRITE(*,*) "Zero Background Diffusivity in stable case."
        write(*,*) 'i=',ii,'j=',jj
        write(*,*) 'temperature'
        write(*,'(5d15.5)') t
        write(*,*)
        write(*,*) 'salinity'
        write(*,'(5d15.5)') s
        write(*,*)
        write(*,*) 'density'
        write(*,'(5d15.5)') rh
        write(*,*)
        write(*,*) 'richardson'
        write(*,'(5d15.5)') ri
        write(*,*)
        write(*,*) 'rid'
        write(*,'(5d15.5)') rid
        write(*,*)
        write(*,*) 'shear'
        write(*,'(5d15.5)') s2
        write(*,*)
        write(*,*) 'BV frequency'
        write(*,'(5d15.5)') an2
        write(*,*)
        write(*,*) 'ustart=',ustar_
        write(*,*)
        write(*,*) 'buoytur=',buoytur
        write(*,*)
        write(*,*) 'buoysol=',buoysol
	      WRITE(*,*) "v_back=",v_back(k), &
                   " t_back=",t_back(k), &
                   " s_back=",s_back(k)
	      WRITE(*,*) " "
	      WRITE(*,*) "slq2_back=",slq2_back
	      WRITE(*,*) "sm_back=",sm_back, &
                   " sh_back=",sh_back, &
                   " ss_back=",ss_back
	      WRITE(*,*) " "
	      WRITE(*,*) "slq2_r1(itheta_r0)=",slq2_r1(itheta_r0), &
                   " slq2_r1(itheta_r1)=",slq2_r1(itheta_r1)
	      WRITE(*,*) "sm_r1(itheta_r0)=",sm_r1(itheta_r0), &
                   " sm_r1(itheta_r1)=",sm_r1(itheta_r1)
	      WRITE(*,*) "sh_r1(itheta_r0)=",sh_r1(itheta_r0), &
                   " sh_r1(itheta_r1)=",sh_r1(itheta_r1)
	      WRITE(*,*) "ss_r1(itheta_r0)=",ss_r1(itheta_r0), &
                   " ss_r1(itheta_r1)=",ss_r1(itheta_r1)
	      WRITE(*,*) " "
	      WRITE(*,*) "back_ra_r1 =", back_ra_r1
	      WRITE(*,*) "theta_r =", theta_r
	      WRITE(*,*) "theta_r_deg =", theta_r_deg
	      WRITE(*,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
	      WRITE(*,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
	      WRITE(*,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
	      WRITE(*,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
	      WRITE(*,*) "n_theta_r_oct=",n_theta_r_oct 
	      WRITE(*,*) "deltheta_r=",deltheta_r
	      WRITE(*,*) " "
	      WRITE(*,*) " "
	      WRITE(*,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
	      WRITE(*,*) "rit=",rit,"ric=",ric
	      WRITE(*,*) " "
	      WRITE(*,*) " "
	      WRITE(*,*) "Program will stop."
	      WRITE(*,*) &
         "************************************************"
          endif
	      STOP
	    END IF
!************************************************************************
	END IF
 20     CONTINUE
!*****C
        END IF
!991108 END SECTION FOR SALINITY MODEL BACKGROUND DIFFUSIVITY CALCULATION.
!************************************************************************
!************************************************************************
!************************************************************************
!980717  Write internal turbulence quantities to fort.91 when writing active.
!981104 Add S_M,H,S to outputs.
!020924X,030404Y Add epsy= dissipation times (\tau Shear)^2 to outputs.
!X       [See NBp020912 Vol.V for epsy formula.]
         IF(ifoutput.EQ.1) THEN
       if (mytid.eq.0) then
	   IF(isurfuse.EQ.0) THEN
	     WRITE(91,9000) z(k),al,slq2,ri1,rid1,sm,sh,ss, &
                      v_back(k),t_back(k),s_back(k)
	   ELSE IF(isurfuse.EQ.1) THEN
	     WRITE(91,9000) z(k),al,slq2,ri1,rid1,sm,sh,ss, &
                      epsy(k), &
                      v_back(k),t_back(k),s_back(k)
	   END IF
       endif
!990208-0322
!991107	Include background turbulence function writeouts for Temp. model bg.
	   IF(ifback.EQ.4.OR.ifsalback.EQ.4) THEN
!*****C
       if (mytid.eq.0) then
	     IF(k.EQ.1) THEN
	       WRITE(94,*) " "
	       WRITE(94,*) &
         "z[cm]      l_back[cm] Ri-table   Ri_d-table Ri_d_back  " &
       //"s2_back    slq2_back  sm_back    sh_back    ss_back    "
	     END IF
	     WRITE(94,9000) z(k),al_back,ri1,rid1,back_rid1, &
                    s2_back,slq2_back,sm_back,sh_back,ss_back 
       endif
!991107	Include background turbulence function writeouts for Temp. model bg.
	   ELSE IF((ifback.GT.4).OR.(ifsalback.GT.4)) THEN
!*****C
       if (mytid.eq.0) then
	     IF(k.EQ.1) THEN
	       WRITE(94,*) " "
!020221[Mummy's 81st Birthday]-22D Introduce additional writeouts to fort.94:
!	                        epsilon/N^2 when used for dimensionalization,
!	                        eps_bot when enhanced bottom mixing is enabled
!	                        and eps_bot(z_bottom) only for special check case.
!22		                Also N^2 and S^2 in all cases.
!22				Note that (N^2)/(S^2) is the true Richardson number,
!22				while ri1, like rid1 can be modified by INTERP2D
!22			        to fit within the table. ri(k) and rid(k) should be pristine.

	       IF(ifepson2.EQ.0) THEN
	         WRITE(94,*) &
         "z[cm]      l_back[cm] Ri-table   Ri_d-table " &
       //"Ri_back    Ri_d_back  ra_r_back  " &
       //"s2_back    slq2_back  sm_back    sh_back    ss_back    " &
       //"N^2        S^2        "
	       ELSE
	         IF(ifbotenhance.EQ.0) THEN
	           WRITE(94,*) &
         "z[cm]      l_back[cm] Ri-table   Ri_d-table " &
       //"Ri_back    Ri_d_back  ra_r_back  " &
       //"s2_back    slq2_back  sm_back    sh_back    ss_back    " &
       //"N^2        S^2        " &
       //"epson2     "
	         ELSE
	           WRITE(94,*) & 
        "z[cm]      l_back[cm] Ri-table   Ri_d-table " &
       //"Ri_back    Ri_d_back  ra_r_back  " &
       //"s2_back    slq2_back  sm_back    sh_back    ss_back    " &
       //"N^2        S^2        " &
       //"epson2     eps_bot    "
	         END IF
	       END IF
	     END IF
        endif

	     IF(ifepson2.EQ.0) THEN
       if (mytid.eq.0) then
	       WRITE(94,9010) z(k),al_back,ri1,rid1, &
                      back_ri1,back_rid1,back_ra_r1, &
                    s2_back,slq2_back,sm_back,sh_back,ss_back, &
                    ri(k)*s2(k),s2(k)
        endif
	     ELSE
	       IF(ifbotenhance.EQ.0) THEN
       if (mytid.eq.0) then
	         WRITE(94,9020) z(k),al_back,ri1,rid1, &
                    back_ri1,back_rid1,back_ra_r1, &
                  s2_back,slq2_back,sm_back,sh_back,ss_back, &
                  ri(k)*s2(k),s2(k), &
	                epson2
       endif
	       ELSE
       if (mytid.eq.0) then
	         WRITE(94,9030) z(k),al_back,ri1,rid1, &
                back_ri1,back_rid1,back_ra_r1, &
              s2_back,slq2_back,sm_back,sh_back,ss_back, &
              ri(k)*s2(k),s2(k), &
	            epson2,eps_bot
       endif
!	Special check of bottom epsilon value.
		 IF((ifcheckbottomeps.EQ.1).AND.(k.EQ.n)) THEN
	           eps_bot__under = &
              eps_bot0 * EXP((z(k+1) - z(kbot))/scale_bot)
       if (mytid.eq.0) then
	           WRITE(94,9040) z(k+1),eps_bot__under
       endif
	         END IF
	       END IF
	     END IF
!*****CD	     
	   END IF
!*****C
	 END IF
!******
!030401 Apply zeroing of l_deep taken from HYCOM turb_2_fixed2.fs0
!030328,050208 To avoid confusion set l_deep to zero when it is not calculated.
        aldeep(k)=0.D0

!020924X ******WHEN HAVE SET "isurfuse=1" DO *NOT* DIMENSIONALIZE WITH THE "MLD" ******
!X       ******BASED "Blackadar" LENGTHSCALE. IN CASES WHERE "isurfuse=0" USES******
!X	 ******"l_MLD \equiv MLD based Blackadar lengthscale" USE INSTEAD******
!X       ******CANUTO'S SURFACE BUOYANCY BASED DIMENSIONALIZATION:******
!X       ******"K_X_foreground = ((epsilon*(\tau S)^2)/(2 S^2))S_X"(SEE NBp020912-12)******
!X       ******EXCEPT WHERE "epsilon*(\tau S)^2 \equiv epsy <0" IN WHICH CASE ******
!X       ******REVERT TO "al^2 S". LEAVE "l_deep^2 S" CASES ALONE.(SEE NBp020920-1,2)******
!000229-0302 In the case where the model is realizable at 
!       the Ri obtained from the external Shear,
!	*but* there is a level above where it is NOT thus realizable, 
!	USE THE "epsilon/(N^2)" DIMENSIONALIZATION FOR "ifepson=2".
!000311 EXCEPT IF  "Ri<0" DO *NOT* USE "epsilon/(N^2)" DIMENSIONALIZATION
!	BECAUSE IT PRODUCES NEGATIVE DIFFUSIVITIES IN THIS CASE.
!	*INSTEAD USE "l_deep^2 S", WHERE "l_deep" IS DERIVED FROM "rho" PROFILE.
!	"|{{d \rho / dz} \over {d2 \rho / dz^2}}| takes place of MLD in l_deep".
!000314 BUT *REVERT* TO "MLD" IN CASES OF FIRST TWO LEVELS !020924 Use epsy if isurfuse=1.
	  IF((ifepson2.EQ.2).AND.(ifbelow.EQ.1)) THEN
!030504Z1 **ALSO USE "l_deep^2 S" IN Gregg et al. LATITUDE DEPENDENCE CASE WHEN (N/f)<1.**
	    IF(ri1.GE.0.D0.OR.((ifdeeplat.EQ.2).AND.(ifnofsmall.EQ.1))) &
       THEN
              tmp=0.5D0*b1**2*ri1*slq2*epson2
	    ELSE IF(k.GT.2) THEN
              IF(k.EQ.n) THEN
                delz = z(k) - z(k-1)
                delrh = rh(k) - rh(k-1)
                del2rh = rh(k) - 2.D0*rh(k-1) + rh(k-2)
              ELSE
                delz = z(k+1) - z(k-1)
                delrh = rh(k+1) - rh(k-1)
                del2rh = rh(k+1) - 2.D0*rh(k) + rh(k-1)
              END IF
              dzrh = delrh/delz
              d2zrh = 4.D0*del2rh/(delz**2)
!000323 rdzlndzrh = *Reciprocal* of Dz_{ln(Dz_{rh})} = Dz_{rh}/Dz2_{rh} .
              rdzlndzrh = dzrh/d2zrh
              al0deep=0.17D0*ABS(rdzlndzrh)
              akz=0.4D0*z(k)
              aldeep(k)=akz*al0deep/(al0deep+akz)
              al2=aldeep(k)*aldeep(k)
!030401Y For case where use pure convection model taken from my HYCOM turb_2_fixed2.
!040217Zi1b When ifshearmin=.TRUE. introduce a minimum foreground shear to avoid
!	    the singularity that occurs in the (lengthscale^2 Shear) dimensionalization
!	    of the foreground turbulence model at zero shear, or else switch to
!	    (lengthscale^2 Brunt Vaisala frequency) dimensionalization. This deals with the
!	    problem that shear=0  implies (l^2 S)/(slq2+1.D-40) =0, which in turn makes the
!	    diffusivities zero, which is unphysical in the unstable case.
              IF(ifpureshear.EQ.1) THEN 
	        GO TO 21
	      ELSE IF(ifshearmin) THEN
		s2(k) = MAX(s2(k),s2min)
	      END IF
              tmp=0.5D0*b1*al2*sqrt(s2(k)/(slq2+1.D-40))
!X In isurfuse=1 case, when buoyancy forcing gives positive dissipation epsilon, 
!X and hence epsy, use epsilon instead of l^2 S as basis for dimensionalization.
	    ELSE
!030401Y For case where use pure convection model taken from my HYCOM turb_2_fixed2.
              IF(ifpureshear.EQ.1) THEN 
	        GO TO 21
	      ELSE IF(ifshearmin) THEN
		s2(k) = MAX(s2(k),s2min)
	      END IF
	      IF(lifepsy(k)) THEN
                tmp=0.5D0*epsy(k)/(s2(k)+1.D-40)
	      ELSE 
                tmp=0.5D0*b1*al2*sqrt(s2(k)/(slq2+1.D-40))
	      END IF
	    END IF
	  ELSE
!030401Y For case where use pure convection model taken from my HYCOM turb_2_fixed2.
              IF(ifpureshear.EQ.1) THEN 
	        GO TO 21
	      ELSE IF(ifshearmin) THEN
		s2(k) = MAX(s2(k),s2min)
	      END IF
!******Zi1b
	    IF(lifepsy(k)) THEN
              tmp=0.5D0*epsy(k)/(s2(k)+1.D-40)
	    ELSE 
              tmp=0.5D0*b1*al2*sqrt(s2(k)/(slq2+1.D-40))
	    END IF
	  END IF
!******X
!030401Y,050208 Calculate \epsilon \tau for the pure convection model from HYCOM turb_2_fixed2.
!030324-28 For very negative Ri use zero shear unstable case approximation.
!       K_X = e \tau S_X = (l^2 \sqrt(-N^2)) (B_1^2/2) (-(\tau N)^2)^{-1/2} S_X
   21    IF(ifpureshear.EQ.1) tmp=0.5D0*(b1**2)*al2*sqrt(-an2(k)/amtaun2)
!******Y
          akm(k)=min(tmp*sm+v_back(k),visc_cbu_limit)
          akh(k)=min(tmp*sh+t_back(k),diff_cbt_limit)
          aks(k)=min(tmp*ss+s_back(k),diff_cbt_limit)
 22   continue
!020219 REFERENCE NOTE: END OF PRIMARY LOOP OVER DEPTH LEVELS.
!980527 Zero diffusivity arrays at points involving rock or land.
      do 33 k=nb+1,nmax 
          akm(k)=0.D0
          akh(k)=0.D0
          aks(k)=0.D0
 33   continue
!980527 Introduce minimum windmixing only where there are at least two ocean levs
      IF(n.GT.0) THEN
        if(akm(1).lt.wndmix) akm(1)=wndmix
        if(akh(1).lt.wndmix) akh(1)=wndmix
        if(aks(1).lt.wndmix) aks(1)=wndmix
      END IF
!980716-19 Write turbulence inputs and outputs to fort.92 when writing active.
!	 Headers for each outputstep.
!981102    Add Ri_d \equiv Ri_T - Ri_C to the outputs.
!000311 Add foreground lengthscale to the outputs. 
!030404Y Add new turb input the square of the Brunt Vaisala frequency to the outputs.
!030425Z Add the Coriolis parameter and rotational lengthscale l_Omega to the outputs.
!25Z	 Also add the surface buoyancy forcing and friction velocity to the outputs.
!25Z	 Calculate the complex diagnostic l_Omega = \sqrt{-B* f_coriolis^{-3}} .
	 zlomega = SQRT(CMPLX(-buoytot/(Coriol**3)))
         IF(ifoutput.EQ.1) THEN
       if (mytid.eq.0) then
	   WRITE(92,*) " "
	   WRITE(92,*) "f_coriolis=",Coriol
	   WRITE(92,*) "MLD[cm] = ",amld
	   WRITE(92,*) "buoytot[cm^2/s^3] =",buoytot, &
                 "    ustar_[cm/s] =",ustar_
	   WRITE(92,*) "l_\\Omega[cm] =",zlomega
	   WRITE(92,*) "z[cm]      tem[C]     sal[ppt]   rho[g/cm3] "// &
                 "Ri         Ri_d	   S^2[/s2]   "// &
                 "N^2[/s2]   "// &
                 "K_M[cm2/s] K_H[cm2/s] K_S[cm2/s] "// &
                 "l_deep[cm] "
       endif
	DO k =1,n
       if (mytid.eq.0) then
	   WRITE(92,9000) z(k),t(k),s(k),rh(k), &
                    ri(k),rid(k),s2(k), &
                    an2(k), &
                    akm(k),akh(k),aks(k), &
                    aldeep(k)
       endif
!000309-12 STOP IF DIFFUSIVITY IS NEGATIVE.
	   IF((akm(k).LT.0.D0).OR.(akh(k).LT.0.D0).OR.(akm(k).LT.0.D0)) &
      THEN
       if (mytid.eq.0) then
	     WRITE(*,*) "Diffusivity is negative."
	     WRITE(*,*) "k=",k
	   WRITE(*,*) "z[cm]      tem[C]     sal[ppt]   rho[g/cm3] "// &
                 "Ri         Ri_d	   S^2[/s2]   "// &
                 "K_M[cm2/s] K_H[cm2/s] K_S[cm2/s] "
	   WRITE(*,9000) z(k),t(k),s(k),rh(k), &
                    ri(k),rid(k),s2(k), &
                    akm(k),akh(k),aks(k)
	     WRITE(*,*) " "
	     WRITE(*,*) "Program will stop."
        endif
	     STOP
	   END IF
	END DO
      END IF	
!******
 9000 FORMAT(12(1pe11.3))
!020221[Mummy's 81st Birthday]-22D
 9010 FORMAT(14(1pe11.3))
 9020 FORMAT(15(1pe11.3))
 9030 FORMAT(16(1pe11.3))
 9040 FORMAT(1pe11.3,154X,1pe11.3)
!*****D
 9001 FORMAT(2(I8,'   ',1pe11.3),8(1pe11.3))
 9050 FORMAT(I8,'  ',2E16.4,I8,'  ')
 9100 FORMAT(I8,'  ',2E12.4,3F11.6,2F11.4)
 9150 FORMAT(F11.3,5E12.4,3F10.6,F9.3)
 9160 FORMAT(F11.3,6E10.4,3F10.6,F9.3)
 9200 FORMAT(I12,'    ',5E16.6)
 9268 FORMAT(I6,5F15.9)
      return
      end



	FUNCTION eplatidepend_(f,an)
!030429-30 Function for use in turb_2 to calculate a latitude, and Brunt Vaisala frequency,
!	dependent factor by which to multiply the constant we had been using for 
!	(\epsilon/ N^2) . Use formula Arccosh = ln(x+\sqrt{x^2 - 1}) in place of "ACOSH".
!	Adapted from program gregglatifunc.f:
!030428 To calculate the latitude dependence of the paper of Gregg,Sanford & Winkel
!	in Nature Vol.422, "Reduced mixing from the breaking of internal waves 
!	in equatorial waters", equation (2).


!	f	= Coriolis parameter, 2 \Omega sin(latitude),  	[sec^{-1}]
!	an 	= Brunt Vaisala frequency, N,  			[sec^{-1}]
!
!	eplatidepend_ = L(\theta,N) from Gregg et. al Nature Vol.422,3 April 2003, eqn.(2).
!
!
!	an0	= Garerett and Munk reference stratification,	[sec^{-1}]
!	f_30	= Coriolis parameter at 30^o reference latitude,[sec^{-1}]
	
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
	REAL(r8) an,an0,f,f_30,eplatidepend_
	PARAMETER(an0=5.24D-3)		!See Gregg et. al in Methods section below eqn.(4).
	REAL(r8) anum,den
	REAL(r8) pi,omega
	
!030430 I use the statement function "ACOSH1", 
!	the functional equivalent of Arccosh I've worked out,
!	which I copy here from the test program gregglatifunc1.f created today,
!	because although there is a built-in f77 function "ACOSH" for kirk (IBM AIX),
!	on the SGI machines ra and halem "f77" does not have "ACOSH".
        ACOSH1(X) = LOG(X + SQRT((X**2) - 1.D0))
!030430 Also introduce a statement function "WAVELAT" for the function 
!	of the two variables `N' and `f' which appears in numerator and denominator 
!	of the Gregg et al. formula for the "latitude dependence" of the wave dissipation:
!	f cosh^{-1}(N/f) .
	WAVELAT(XF,YN) = XF*ACOSH1(YN/XF) 	


!	Calculate omega for consistency as it was done in ocean.F of the NCAR code.
!	The value of omega was used to set fcor, 
!	which is interpolated to Coriol in cctmix.F .
      pi     = 4.D0*atan(1.D0)
      omega  = pi/43082.0D0


!	f_30 = 2.*\Omega*sin(30^o), but sin(30^o)=sin(pi/6)=1/2, so
!	f_30 = 2.*\Omega*(1/2) = \Omega .
	f_30=omega

        den=WAVELAT(f_30,an0)

	anum = WAVELAT(f,an)
	eplatidepend_ = anum/den


	END
	
!-----------------------------------------------------------------------
!     finds mixed layer depth
!-----------------------------------------------------------------------
!980501 Choice of definitions.
      subroutine formld(z,t,amld,n)
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
!     z and amld are both positive
      dimension z(n),t(n)
      do k=1,n
          if (abs(t(k) - t(1)).gt.0.1) then
#if (defined D_PRECISION)
              tm = t(1) - sign(0.1D0,t(1) - t(k))
#else
              tm = t(1) - sign(0.1,t(1) - t(k))
#endif
              amld = z(k) + (z(k-1) - z(k))* &
        (tm - t(k))/(t(k-1) - t(k) + 1.e-20)
              goto 11
          endif
      enddo
      amld=z(n)
 11   continue
      return
      end

!0501 Temperature difference criterion used, but chose outside \Delta T = delte.
      subroutine formld_te(z,t,delte,amld,n)
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
!     z and amld are both positive
      dimension z(n),t(n)
      do k=1,n
          if (abs(t(k) - t(1)).gt.delte) then
              tm = t(1) - sign(delte,t(1) - t(k))
              amld = z(k) + (z(k-1) - z(k))* &
        (tm - t(k))/(t(k-1) - t(k) + 1.e-20)
              goto 11
          endif
      enddo
      amld=z(n)
 11   continue
      return
      end

!0501 Density difference criterion used, but chose outside \Delta \rho = delrh.
      subroutine formld_rh(z,t,delrh,amld,n)
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
!     z and amld are both positive
      dimension z(n),t(n)
      do k=1,n
          if (abs(t(k) - t(1)).gt.delrh) then
              tm = t(1) - sign(delrh,t(1) - t(k))
              amld = z(k) + (z(k-1) - z(k))* &
        (tm - t(k))/(t(k-1) - t(k) + 1.e-20)
              goto 11
          endif
      enddo
      amld=z(n)
 11   continue
      return
      end
!******


!-----------------------------------------------------------------------
!     beginning of improved turbulence model subroutines (nmodel=1)
!-----------------------------------------------------------------------
      subroutine ccoeff(b1,rimax,g_tur,d_tur,s_tur)
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
!
      real(r8) :: g_tur(8), d_tur(0:5), s_tur(0:9)
      real(r8) :: g1,g2,g3,g4,g5,g6,g7,g8
      real(r8) :: d0,d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9
!----------------------------------------------------------------------
!     defines constants and model coefficients
!     b1 == (q/ustar)**3,  qus == q/ustar
!     g2=(1.-3.*v2dq2)/2.
!     g3=(3.*u2dq2-1.-g2)/3.
!     w2dq2=(1+g2-3.*g3)/3.
!     shih & shabbir: g2=0.00336, g3=0.0906
!     apbl: g2=0.06, g3=0.16
!     beta5=0.6*(1-c5), if c5=.3 then beta5=.42
!     qty == us**2*theta2/(w*theta)**2
!     cthem1 == 1/(c sub theta)
!-----------------------------------------------------------------------
      prt0=1.0
      b1=16.6d0
      qus=b1**(1./3.)
      g2 = 0.00336
      g3 = 0.0906
      g6 = 0.4
      g7 = 0.
      beta5=0.42
      qty=3.1
      g1=4./3.*(3.*g3**2-g2**2)+4.*qus**(-4)
      g4=15./8.*beta5*g1
      temp1=1./6.*prt0*(1.+g2-3.*g3)*qus**4
      temp2=9.*(g6-g7)*(g6+g7+2.*prt0)/ &
      ( (1.+g2-3.*g3)*prt0 )**2 *qus**(-4)
      g5=temp1*(1.+sqrt(1.+temp2))
      cthem1=1./(prt0*qus**2)* qty  ! 0.952766
      c7=1./3.
      g8 = (1.-c7)*cthem1
!
      s0=3./2.*g1*g5*g5
      s1=-(g6+g7)*g4 + 2.*(g1-g3-1./3.*g2)*g4*g5 + 3./2.*g1*g5*g8
      s2=-3./8.*(g6*g6-g7*g7)*g1
      s3=3./2.*(g6 + g7)*g4 + (3.*g3 + g2)*g4*g5
      s4=2.*g5
      s5=2.*g4
      s6=-2./3.*g5*(g2*g2-3.*g3*g3)+1./2.*g1*g5*(g2-3.*g3)+3./4.*g1*(g6-g7)
      s7=3.*g5
      s8=3.*g4
      s9=g5*(3.*g3*g3-g2*g2)
!
      d0=3.*g5*g5
      d1=g5*(7.*g4 + 3.*g8)
      d2=g5*g5*(3.*g3*g3 - g2*g2) - 3./4.*(g6*g6 - g7*g7)
      d3=g4*(4.*g4 + 3.*g8)
      d4=g4*(g2*g6 - 3.*g3*g7 - g5*(g2*g2 - g3*g3)) + g5*g8*(3.*g3*g3 - g2*g2)
      d5=1./4.*(g2*g2 - 3.*g3*g3)*(g6*g6 - g7*g7)
!     find rimax:
      aa=(d3+g4)
      bb=(d4-s1/2.+s6/2.)
      cc=d5-s2/2.
      rimax=(-bb+sqrt(bb**2-4.*aa*cc))/(2*aa)
      rimax=rimax*0.999  !ad hoc
!
      g_tur(1)=g1
      g_tur(2)=g2
      g_tur(3)=g3
      g_tur(4)=g4
      g_tur(5)=g5
      g_tur(6)=g6
      g_tur(7)=g7
      g_tur(8)=g8
      s_tur(0)=s0
      s_tur(1)=s1
      s_tur(2)=s2
      s_tur(3)=s3
      s_tur(4)=s4
      s_tur(5)=s5
      s_tur(6)=s6
      s_tur(7)=s7
      s_tur(8)=s8
      s_tur(9)=s9
      d_tur(0)=d0
      d_tur(1)=d1
      d_tur(2)=d2
      d_tur(3)=d3
      d_tur(4)=d4
      d_tur(5)=d5
      return
      end

      subroutine ccoeff1(b1,rimax,g_tur,d_tur,s_tur)
!980912 Created to implement the option of using the same constants Ye Cheng
!	uses for his atmospheric code at this date: fit to experimental data.
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
      real(r8) :: g_tur(8), d_tur(0:5), s_tur(0:9)
      real(r8) :: g1,g2,g3,g4,g5,g6,g7,g8
      real(r8) :: d0,d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9
!----------------------------------------------------------------------
!     defines constants and model coefficients
!     b1 == (q/ustar)**3,  qus == q/ustar
!     g2=(1.-3.*v2dq2)/2.
!     g3=(3.*u2dq2-1.-g2)/3.
!     w2dq2=(1+g2-3.*g3)/3.
!     shih & shabbir: g2=0.00336, g3=0.0906
!     apbl: g2=0.06, g3=0.16
!     beta5=0.6*(1-c5), if c5=.3 then beta5=.42
!     qty == us**2*theta2/(w*theta)**2
!     cthem1 == 1/(c sub theta)
!-----------------------------------------------------------------------
      prt0=1.0
      b1=16.6d0
      qus=b1**(1./3.)
!980912 Values used in atmosphere given me by Ye Cheng on 980911.
      g1=0.192D0/2.D0
      g2=0.06D0
      g3=0.16D0
      g4=0.1D0
      g5=7.66D0
      g6=0.4D0
      g7=0.D0
      g8=0.29D0
!
      s0=3./2.*g1*g5*g5
      s1=-(g6+g7)*g4 + 2.*(g1-g3-1./3.*g2)*g4*g5 + 3./2.*g1*g5*g8
      s2=-3./8.*(g6*g6-g7*g7)*g1
      s3=3./2.*(g6 + g7)*g4 + (3.*g3 + g2)*g4*g5
      s4=2.*g5
      s5=2.*g4
      s6=-2./3.*g5*(g2*g2-3.*g3*g3)+1./2.*g1*g5*(g2-3.*g3) +3./4.*g1*(g6-g7)
      s7=3.*g5
      s8=3.*g4
      s9=g5*(3.*g3*g3-g2*g2)
!
      d0=3.*g5*g5
      d1=g5*(7.*g4 + 3.*g8)
      d2=g5*g5*(3.*g3*g3 - g2*g2) - 3./4.*(g6*g6 - g7*g7)
      d3=g4*(4.*g4 + 3.*g8)
      d4=g4*(g2*g6 - 3.*g3*g7 - g5*(g2*g2 - g3*g3)) + g5*g8*(3.*g3*g3 - g2*g2)
      d5=1./4.*(g2*g2 - 3.*g3*g3)*(g6*g6 - g7*g7)
!     find rimax:
      aa=(d3+g4)
      bb=(d4-s1/2.+s6/2.)
      cc=d5-s2/2.
      rimax=(-bb+sqrt(bb**2-4.*aa*cc))/(2*aa)
      rimax=rimax*0.999  !ad hoc
!
      g_tur(1)=g1
      g_tur(2)=g2
      g_tur(3)=g3
      g_tur(4)=g4
      g_tur(5)=g5
      g_tur(6)=g6
      g_tur(7)=g7
      g_tur(8)=g8
      s_tur(0)=s0
      s_tur(1)=s1
      s_tur(2)=s2
      s_tur(3)=s3
      s_tur(4)=s4
      s_tur(5)=s5
      s_tur(6)=s6
      s_tur(7)=s7
      s_tur(8)=s8
      s_tur(9)=s9
      d_tur(0)=d0
      d_tur(1)=d1
      d_tur(2)=d2
      d_tur(3)=d3
      d_tur(4)=d4
      d_tur(5)=d5
!
      return
      end

      subroutine ourl2(b1,ri,slq2,sm,sh,g_tur,d_tur,s_tur)
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
!     slq=s*l/q, slq2=slq**2, 1/2*q**2=e
      real(r8) :: g_tur(8), d_tur(0:5), s_tur(0:9)
      real(r8) :: g1,g2,g3,g4,g5,g6,g7,g8
      real(r8) :: d0,d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9
!
      g1=g_tur(1)
      g2=g_tur(2)
      g3=g_tur(3)
      g4=g_tur(4)
      g5=g_tur(5)
      g6=g_tur(6)
      g7=g_tur(7)
      g8=g_tur(8)
      s0=s_tur(0)
      s1=s_tur(1)
      s2=s_tur(2)
      s3=s_tur(3)
      s4=s_tur(4)
      s5=s_tur(5)
      s6=s_tur(6)
      s7=s_tur(7)
      s8=s_tur(8)
      s9=s_tur(9)
      d0=d_tur(0)
      d1=d_tur(1)
      d2=d_tur(2)
      d3=d_tur(3)
      d4=d_tur(4)
      d5=d_tur(5)
!
!
      aa = (d3+g4)*ri*ri+(d4-s1/2.+s6/2.)*ri+d5-s2/2.
      bb = (d1+g5)*ri+d2-s0/2.
      cc = d0
      temp = bb**2-4.*aa*cc
      if(temp.le.0.) then
         write(*,*) 'in turb_2, temp <= 0, stop, ri =',ri
         stop
      endif
#if (defined D_PRECISION)
      q = -1./2.*(bb+sign(1.D0,bb)*sqrt(temp))
#else
      q = -1./2.*(bb+sign(1.0,bb)*sqrt(temp))
#endif
      if(bb.lt.0.) then
          y=cc/q
      else
          y=q/aa
      endif
      if(y.lt.0.) then
       if (mytid.eq.0) then
          write(*,*) 'in turb_2, y < 0, stop; ri =', ri
       endif
          stop
      endif
      dd=d0+(d1*ri+d2)*y+(d3*ri*ri+d4*ri+d5)*y*y
      sm=(s0+(s1*ri+s2)*y)/dd
      sh=(s4+(s5*ri+s6)*y)/dd
      slq2=y/(b1*b1)
      return
      end

!---------------------------------------------------------------------
!     end of improved turbulence model subroutines (nmodel=1)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     beginning of parameter-free turbulence model subroutines(nmodel=2)
!-----------------------------------------------------------------------
      subroutine mcoeff(b1,rimax,ttot,tcot,tctot,tptot,tpcot,tpvot)
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
      real(r8)::  ttot,tcot,tctot,tptot,tpcot,tpvot
!      implicit real*8 (a-h,o-z)
!     defines constants and model coefficients
      b1=24.7D0
!980716 Use the value of rimax found in 
!	/data1/acamh/TURB1/CHENG/OCode/Canuto/3D/Cheng/mike_12.f_980528 .
      rimax=1.68D0
      return
      end
      subroutine mikel2(b1,ri,slq2,sm,sh,ri1)
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
      real(r8) :: ri1
      data ifirst/0/
      external fct
      ri1=ri
      eps=1.e-6                                              
      iend=10000                                             
!     rtwi finds the root of x=fct(x)                     
!     to make a guess at bst                     
      if(ifirst.eq.0) then
        yst=0.7   !at ri=-20
        ifirst=1
      else
        yst=y                                
      endif
      call rtwi(y,val,fct,yst,eps,iend,ier,ttot,tcot,tctot,tptot,tpcot,tpvot,ri1,rit1,ric1,rit,ric)
!      write(*,*) ri,ier,y,val
      if(ier.ne.0) then
       if (mytid.eq.0) then
         write(*,*) "rtwi call problem, ier=",ier
       endif
         stop
      endif
      x=ri*y
      call mksmsh(x,y,sm,sh)
      slq2=y/(b1*b1)
      return
      end
      function fct(y,ttot,tcot,tctot,tptot,tpcot,tpvot,ri,rit1,ric1,rit,ric)
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      implicit real*8 (a-h,o-z)
      real(r8) :: ri,rit1,ric1,rit,ric
      real(r8) :: ttot,tcot,tctot,tptot,tpcot,tpvot
      x=ri*y
      call mksmsh(x,y,sm,sh)
      fct=2./(sm-ri*sh)
      return                                          
      end                                            

!980716 Use as mksmsh:
!980716 The Dubovikov routine for S_{M,H} from 
!	/data1/acamh/TURB1/CHENG/OCode/Canuto/3D/Cheng/mike_12.f_mod980528.
      subroutine mksmsh(x,y,sm,sh)
!980717 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
!
      parameter(x1min=15.d0/8.d0)
      x1=0.5d0*sqrt(y)
      if(x1.lt.x1min) then
         fc=5.d0/9.d0
      else
         fc=1.d0-5.d0/(3.d0*x1)+25.d0/(16.d0*x1**2)
      endif
      y=(2.d0*x1)**2
      ri=x/y
      a=1.d0+1.d0/(1.+0.16d0*x1**2*ri)
      eta=0.05d0*x1**2*ri*a
      sm=1.d0/25.d0*fc/(1.d0+5.d0/(9.d0*fc)*eta)
      xi=x/60.d0
      sh=0.056d0/(1.d0+2.4d0*xi)
      return
      end
!981022 Comment out for use with 2D table since oursal2 routine has a copy.
!      subroutine rtwi(x,val,fct,xst,eps,iend,ier)                      
!C980501 Make double precision to conform to calling cctmix routine.
!      implicit real*8 (a-h,o-z)
!c     to solve general nonlinear equations of the form x=fct(x)       
!c     by means of wegsteins iteration method                         
!c     prepare iteration                                             
!      ier=0                                                        
!      tol=xst                                                     
!      x=fct(tol)                                                 
!      a=x-xst                                                   
!      b=-a                                                     
!      tol=x                                                   
!      val=x-fct(tol)                                         
!c     start iteration loop                                 
!      do 6 i=1,iend                                       
!      if(val) 1,7,1                                      
!c     equation is not satisfied by x                    
! 1    b=b/val-1.                                       
!      if(b) 2,8,2                                     
!c     iteration is possible                          
! 2    a=a/b                                         
!      x=x+a                                        
!      b=val                                       
!      tol=x                                      
!      val=x-fct(tol)                            
!c     test on satisfactory accuracy            
!      tol=eps                                 
!      d=abs(x)                               
!      if(d-1.) 4,4,3                        
! 3    tol=tol*d                            
! 4    if(abs(a)-tol) 5,5,6                
! 5    if(abs(val)-10.*tol) 7,7,6         
! 6    continue                          
!c     end of iteration loop                                           
!c     no convergence after iend iteration steps. error return.       
!      ier=1                                          
! 7    return                                        
!c     error return in case of zero divisor         
! 8    ier=2                                       
!      return                                     
!      end                                       
!*****C
!-----------------------------------------------------------------------
!     end of parameter-free turbulence model subroutines (nmodel=2)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     end of turbulence models
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!030403Y    Correct the omission of "ttot"=`{\tau_{\theta} \over \tau}' from
!Y	    common block /bb0/. This routine goes with turb_2f [See NBp.030403-8to10.].
!030128B    Correct error of omitting to set "b1" to "b1_0" [See NBp.030128-4 .]
!030124B    Correct error made in ifchengb1=.FALSE. option [See NBp.030124-18 Vol.XVII .]
!021210,11A    Introduce option to allow B1 to be calculated from Ye Cheng's formula: 
!A	    B1 = {y(Ri_Tem=0,Ri_Sal=0)}^{3/4}.
!000316     Make use of double precision more consistent.
!000210	    Replace the numerical value of 6.25 by 1/(tpvot**2) .
!000201     Version in which following OTsalche/plot000127 
!           the timescale ratios are calculated in the 'smshsc' routine 
!	    and passed back hrough the common block bb0/
!	    to simplify the process of adjustment of timescale ratios.
!981016-1110 Submodule to calculate turbulence functions (Sl/q)^2 and S_M,S_H,S_S
!          of Ri(=Ri_T+Ri_C) and Ri_d(=Ri_T-Ri_C) in our NCAR turbulence module.
!****************************************
!981104 Pass B_1 back to the calling routine as b1_arg.
!_______________________________________________________________________
!-----------------------------------------------------------------------
	   subroutine oursal2_zeroshear(and2on2,amtaun2,sm,sh,sc,ttot,tcot,tctot,tptot,tpcot,tpvot)
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
!
!030321-27AH SUBROUTINE TO CALCULATE THE TURBULENCE MODEL FOR THE ZERO SHEAR 
!AH	  *UNSTABLE CASE*. TO DIMENSIONALIZE DIFFUSIVITIES IN THAT CASE NEED
!AH	  THE POSITIVE QUANTITY -(\tau N)^2 WHICH I NAME "amtaun2".
! --- Submodule to calculate turbulence functions, -(\tau N)^2 and S_M,S_H,S_S ,
! --- of the difference of the ratios of the heat and salt contributions 
! --- to the square of the Brunt Vaisala frequency, N_d^2/N^2 , named "and2on2".
! --- Calls quad_sal_pureconv.
! --- Stripped down and altered from the winter 2003 hycom version of oursal2_1a.
! --- Version in which following OTsalche/plot000127 
! --- the timescale ratios are calculated in the 'smshsc' routine 
! --- and passed back through the common block bb0/
! --- to simplify the process of adjustment of timescale ratios.
! --- Stripped and adapted from plot981007.f.
! --- Program to generate contour and 1 variable plots vs. Ri,Ri_d based on
! --- Program to generate contour plots vs. Ri_T and Ri_C based on
! --- Program to generate plots vs. Ri_T at different Ri_C values based on
! --- CORRECTED PROGRAM WITH NEW VALUE OF "p10". 'p10 = tpt*tct/(tc**2)'
! --- Program to generate K_X/((l^2) S) for Canuto based on plot980609.f:
! --- Program to generate data for plots of turbulence functions including
! --- S_{M,H,C} and Canuto's new y = (\tau_pv S)^2
! --- and n,c as functions of stability parameters in the concentration theory
! --- (structure is a 1 point closure like the generalized Mellor-Yamada, 
! --- but the constants are derived based on Dubovikov's model according
! --- to Ye Cheng). The concentration theory dimensionless parameters
! --- associated with the squares of shear, temperature contribution to 
! --- Brunt Vaisala frequency and concentration contribution to it,
! --- the new y,n,c are represented in this program by the variables
! --- c_y,c_n,c_c. 
! --- Adapted from Cheng's program mike_12.f_980528 for the Dubovikov model.
!-----------------------------------------------------------------------
! --- Inputs:
! --- and2on2  	N_d^2/N^2   Difference of heat and salt fractions of N^2
! --- 
! --- Internal quantities:
! --- anh2on2 	N_H^2/N^2   Fraction of N^2 contributed by temperature gradient
! --- ans2on2 	N_S^2/N^2   Fraction of N^2 contributed by salinity gradient
! --- taunh2    (\tau N_H)^2 dimensionless variable for temperature gradient
! --- tauns2    (\tau N_S)^2 dimensionless variable for salinity gradient
! ---
! --- Outputs:
! --- amtaun2			-(\tau N)^2 
! --- sm			S_M
! --- sh			S_H
! --- ss			S_S
!-----------------------------------------------------------------------
!
! --- kx=e*tau*sx=1/2*(b1*l)**2*sqrt(-n^2)/(-(\tau n^2))**(1/2)*sx
!
! --- X = {M,H,C} . 
! --- The program variable "amtaun2" is -(\tau N^2). 
! --- Since \tau=B_1 l/q, (N l/q)^2 = (\tau N^2)/(B1^2) .
! --- Interested in the unstable case so (\tau N^2) is negative.
!
!     parameter(sgmt=0.72)    !Make "sgmt" a parameter.
!    .                          !Standard value was 0.72.
!     parameter(tptot0=(1./5.)*(1./(1.+(1./sgmt)))) ! \tau_p\theta over \tau
!     parameter(tpcot0=tptot0)                                !tau_pc over \tau
!     parameter( ttot0=sgmt)                                  !tau_\theta over \tau
!     parameter( tcot0=ttot0)                                 !tau_c over \tau
!     parameter(tctot0=1./3.)                             ! tau_c\theta } over \tau
!     parameter(tptot = tptot0)
!     parameter(tpcot = tpcot0)
!     parameter(ttot  = ttot0)
!     parameter(tcot  = tcot0)
!     parameter(tctot = tctot0)
!
!
!
! --- Commented excerpt from the file "sx"
!
!
! --- sgmt := 0.72;
!
! --- tpt := 1/(5*(1+1/sgmt))*tau;
! --- tpt  = .08372093019*tau
!
! --- tpc := 1/(5*(1+1/sgmt))*tau;
! --- tpc  = .08372093019*tau
!
! --- tt := sgmt*tau;
! --- tt   = .72*tau
!
! --- tc := sgmt*tau;
! --- tc   = .72*tau
!
! --- tct := 2/15*sgmt*tau;
! --- tct  = .09599999998*tau
!
! --- Calculate the timescale ratios in the 'smshsc' routine instead of here.
! --- Set \sigma_t0. sgmt = .72
!
! --- Calculate {\tau_C \over \tau} and {\tau_{C\theta} \over \tau}.
! --- tcot  = sgmt
! --- tctot = (2./15.)*sgmt
! --- "tpt/tau" and "tpc/tau" from the "sx" excerpt
! --- tptot = 1./(5.*(1+1/sgmt))
! --- tpcot = 1./(5.*(1+1/sgmt))
!
!980610-030403 Common block with ratios of timescales
      real(r8)::  ttot,tcot,tctot,tptot,tpcot,tpvot
! --- Timescale ratios are now calculated in the 'smshsc' subroutine.
! --- Make dummy call with c_y=c_n=c_c=0 to get their values for initial use.
#if (defined D_PRECISION)
      call smshsc_a3(0.D0,0.D0,0.D0,sm,sh,sc,ttot,tcot,tctot,tptot,tpcot,tpvot)
#else
      call smshsc_a3(0.,0.,0.,sm,sh,sc,ttot,tcot,tctot,tptot,tpcot,tpvot)
#endif
!
      eps=1.E-6                                              
      iend=300                                              

	
!
!030321-22AH  N^2 \equiv N_H^2 + N_S^2 ; N_d^2 \equiv N_H^2 - N_S^2 .  
! ---         N_H^2 = (N^2 + N_d^2)/2  ; N_S^2 = (N^2 - N_d^2)/2 .
! ---	      N_H^2/N^2 = (N^2 + N_d^2)/(2N^2) = (1 + N_d^2/N^2)/2 
! ---	      N_S^2/N^2 = (N^2 - N_d^2)/(2N^2) = (1 - N_d^2/N^2)/2 
	 anh2on2 = (1. + and2on2)/2.
	 ans2on2 = (1. - and2on2)/2.
!030321AH **SOLVE A QUADRATIC EQUATION TO FIND -(\tau N)^2 **
! ---	  **FOR PURE TURBULENT CONVECTION INCLUDING SALINITY (SHEAR=0).**
	 CALL QUAD_SAL_PURECONV(anh2on2,ans2on2,amtaun2,ier,ttot,tcot,tctot,tptot,tpcot,tpvot)
!
         if(ier.NE.0) then
       if (mytid.eq.0) then
	    WRITE(*,*) "In oursal2_zeroshear subroutine"
	    WRITE(*,*) "anh2on2",anh2on2,"	ans2on2=",ans2on2
            WRITE(*,*) "Error returned by quad_sal_pureconv ier=",ier
	    WRITE(*,*) "ier=1 means choice of root must be reconsidered" &
               //"for these model constants."
       endif
            STOP
         endif
!
!030322AH-0404Y CALCULATE TEMPERATURE AND SALINITY GRADIENT DIMENSIONLESS TURBULENCE FUNCTIONS.
	 taunh2 = -amtaun2*anh2on2
	 tauns2 = -amtaun2*ans2on2
!030322AH  c_y \equiv (\tau_pv/\tau)^2 x (1 - R_\rho) Ri^(-1) .
! ---	   Since the shear here is zero the reciprocal of Ri is zero   
! ---      and therefore c_y = 0 .
!
! --- From equations (12) and (15) of Ocean Turbulence II, JPO Vol.32 p.240 Jan.2002:
! --- "n = -(\tau_s\theta/\tau)(tau_s/\tau) x" , which in the variables used here becomes
! --- n = - {\tau_C \over \tau } {\tau_{C\theta} \over \tau } {\tau N_H}^2
! --- and
! --- "c = {\tau_s \over \tau}^2 R_rho x" , which becomes  
! --- c = - {\tau_C \over \tau}^2  {\tau N_S}^2
!
	 c_y = 0.
         c_n = -(tcot*tctot)*taunh2
         c_c = -(tcot**2)*tauns2
         call smshsc_a3(c_y,c_n,c_c,sm,sh,sc,ttot,tcot,tctot,tptot,tpcot,tpvot)
!030304Y Stop in case any of the S_X is negative.
	 IF((sm.LT.0.D0).OR.(sh.LT.0.D0).OR.(ss.LT.0.D0)) THEN
       if (mytid.eq.0) then
	   WRITE(*,*) " "
	   WRITE(*,*) &
      "In subroutine oursal2_zeroshear after call to smshsc_a3"
	   WRITE(*,*) "Error: negative structure function."
	   WRITE(*,*) "sm=",sm
	   WRITE(*,*) "sh=",sh
	   WRITE(*,*) "ss=",ss
	   WRITE(*,*) " "
	   WRITE(*,*) "c_y=",c_y
	   WRITE(*,*) "c_n=",c_n
	   WRITE(*,*) "c_c=",c_c
	   WRITE(*,*) " "
	   WRITE(*,*) "taunh2=",taunh2
	   WRITE(*,*) "tauns2=",tauns2
	   WRITE(*,*) "anh2on2=",anh2on2
	   WRITE(*,*) "ans2on2=",ans2on2
	   WRITE(*,*) "and2on2=",and2on2
	   WRITE(*,*) " "
	   WRITE(*,*) "Program is stopping."
       endif
	   STOP
	 END IF
!
!
 1003 format(12(I8))
 1004 format(12(1pe14.5))
      end



!_______________________________________________________________________
!-----------------------------------------------------------------------
!030403Y    I *MUST* add "ttot"=`{\tau_\theta \over \tau}' to the "/bb0/" common block with 
!Y	    timescales, because I use it in this routine. Trouble with this routine made me
!Y	    realize that the common block was missing it. [See NBp.030407-12.]
!030401-02Y Adaptation of the hycom routine to the NCOM must make explicity double precision.
!Y	    Fix error of hyphen for minus sign in amtaun2_negri_0rid test and
!Y	    remove common_blocks_giss.h inclusion since I don't have that file in current NCOM
!Y	    and instate the common block bb0 with the timescale ratios [See NBp.030402-2.]
      SUBROUTINE quad_sal_pureconv(anh2on2,ans2on2,amtaun2,ier,ttot,tcot,tctot,tptot,tpcot,tpvot)
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
!030322-27AH To calculate -(\tau N)^2 as a function of (N_H^2/N^2) and (N_S^2/N^2) .
!
!
!
!-----------------------------------------------------------------------
! --- Inputs:
! --- anh2on2   N_H^2/N^2    heat fraction of square of Brunt Vaisala frequency
! --- ans2on2   N_H^2/N^2    salt fraction of square of Brunt Vaisala frequency
! ---
! --- Internal quantities:
! --- rootplus   	   solution of quadratic for -(\tau N)^2 with + root
! --- rootminus  	   solution of quadratic for -(\tau N)^2 with - root
! ---
! --- Outputs:
! --- amtaun2  -(\tau N)^2  negative of dimensionless variable for density gradient 	
!-----------------------------------------------------------------------


!     parameter(sgmt=0.72)    !Make "sgmt" a parameter.
!    .                          !Standard value was 0.72.
!     parameter(tptot0=(1./5.)*(1./(1.+(1./sgmt)))) ! \tau_p\theta over \tau
!     parameter(tpcot0=tptot0)                                !tau_pc over \tau
!     parameter( ttot0=sgmt)                                  !tau_\theta over \tau
!     parameter( tcot0=ttot0)                                 !tau_c over \tau
!     parameter(tctot0=1./3.)                             ! tau_c\theta } over \tau
!     parameter(tptot = tptot0)
!     parameter(tpcot = tpcot0)
!     parameter(ttot  = ttot0)
!     parameter(tcot  = tcot0)
!     parameter(tctot = tctot0)
!980610-030403 Common block with ratios of timescales
      real(r8)::  ttot,tcot,tctot,tptot,tpcot,tpvot


!
!030325 Introduce a limit on the acceptable error in the quadratic.
      PARAMETER(errorbound=1.D-12)
!030326-0402 An approximate value for -(\tau N)^2 obtained from the model with shear 
!	at Ri_d=0 and Ri = -56.94668746578522 the minimal table value using 
!	ntbl=501,ri0=-4.,nextrtbl0=61,nextrtbl1=1000 and nposapprox=101, and B1=16.6.
      PARAMETER(amtaun2_negri_0rid= (16.6**2)* &
          (-(6.5489662174209907D-04)*(-56.94668746578522)))


	ier=0
	amtaun2=0. 	!Initialize the output variable.
	
!030323-26 Ocean Turbulence II eq. (15) includes x \equiv (\tau N)^2(1-R_\rho)^{-1}
!	this translates to x \equiv (\tau N)^2 (1 + Ri_Sal/Ri_Tem)^{-1}
!	x = (\tau N)^2 (Ri_Tem/Ri) = (\tau N_H)^2
!	eq.(12) pi_{1,2,3,4,5} = (\tau_{ps},\tau_{s\theta},\tau_s,\tau_{p\theta},\tau_\theta)\tau^{-1} 
!	Ocean Turbulence II gives the zero shear Production=Dissipation solution in eq.(18).
!	eq.(18b)
!	A + B x^{-1} - x^{-2} = 0
!	eq.(18c):
!	{15 \over 7} A = \pi_1(\mu - \pi_2 \pi_4)R_\rho - \pi_4(\eta + \pi_1 \pi_2 R_\rho)
!	                 -{a5 \over 7} (\eta \mu + \pi_1 \pi_2^2 \pi_4 R_\rho)
!	
!	{15 \over 7} B = \pi_1 R_\rho - \pi_4 - {15 \over 7}(\eta + \mu)
!
!	\eta = \pi_1(\pi_2 - \pi_3 R_\rho),
!	\mu  = \pi_4(\pi_5 - \pi_2 R_\rho)
!
!	030321-26 made the substitutions:
!	A''' = A''/N^4 = A' (S^4/N^4) = A (N_H^4/N^4)
!	B''' = B''/N^2 = B' (S^2/N^2) = B (N_H^2/N^2)
!	\eta''' = \eta''/N^2 = \eta' (S^2/N^2) = \eta (N_H^2/N^2)
!	\mu'''  = \mu'' /N^2 = \mu' (S^2/N^2)  = \mu  (N_H^2/N^2)
!
!       which used in (18b,c) yielded:
!	ax'''^2 + bx''' + c = 0 
!	a = A''' ; b = -B''' ; c = -1 .	x''' \equiv -(\tau N)^2 .
!	=> -(\tau N)^2 = {{-b +/- \sqrt{b^2 - 4ac}} \over 2a }
!030325 Numerical Recipes by Press,Flannery,Teukolsky and Vetterling c.1986, Section 5.5,
!	warns that solution of the quadratic  as "{-b +/- \sqrt{b^2 - 4ac}} \over 2a}"
!	is "asking for trouble". 
!	They recommend	instead using for the two roots "q/a" and "c/q" with
!	"q \equiv -{1 \over 2} [ b + sgn(b)\sqrt{b^2 - 4ac}]" .
!	************************
!	-(\tau N)^2 = q/a
!	or
!	-(\tau N)^2 = q/c

!
!	A''' = - (7/15) [ (\tau_{ps}/\tau)
!	                  (\mu''' - (\tau_{p\theta}/\tau)(\tau_{s\theta}/\tau)(N_H^2/N^2))
!                         (N_S^2/N^2) 
!	                 +(\tau_{p\theta}/\tau)
!		          (\eta''' - (\tau_{ps}/\tau)(\tau_{s\theta}/\tau)(N_S^2/N^2))
!		 	  (N_H^2/N^2) ]
!	               - ( \eta''' \mu''' 
!			  - (\tau_{p\theta}/\tau)(\tau_{s\theta}/\tau)
!		 	    (\tau_{s\theta}/\tau)^2 (N_H^2/N^2)(N_S^2/N^2) ) 
!
!	B''' = - (7/5) [ (\tau_{ps}/\tau)(N_S^2/N^2) + (\tau_{p\theta}/\tau)(N_H^2/N^2) ]
!		      - ( \eta''' + \mu''' ) 
!	
!
!	\eta''' = (\tau_{ps}/\tau) 
!	          ((\tau_{s\theta}/\tau)(N_H^2/N^2) + (\tau_s/\tau)(N_S^2/N^2))
!
!	\mu'''  = (\tau_{p\theta}/\tau) 
!	          ((\tau_{s\theta}/\tau)(N_S^2/N^2) + (\tau_\theta/\tau)(N_H^2/N^2))
!
!	"ttot"   = \tau_theta/\tau   	"tcot" = \tau_s/\tau  
!	"tctot" = \tau_{s\theta}/\tau     
!	"tptot" = \tau_{p\theta}/\tau   "tpcot" = \tau_{ps}/\tau
!	Note that our heat salt symmetry requires ttot=tcot and tptot=tpcot .
!	With this symmetry \eta''' changes into \mu''' under heat/salt exchange.


      eta3p = tpcot*(tctot*anh2on2 + tcot*ans2on2)
      amu3p = tptot*(tctot*ans2on2 + ttot*anh2on2)


      braasal = tpcot*(amu3p  - tptot*tctot*anh2on2)*ans2on2
      braatem = tptot*(eta3p - tpcot*tctot*ans2on2)*anh2on2
      braa = braasal+braatem
      para = (eta3p*amu3p - tptot*tpcot*(tctot**2)*anh2on2*ans2on2)
      a3p   = -(7./15.)*braa - para

      brabsal = tpcot*ans2on2
      brabtem = tptot*anh2on2
      brab = brabsal+brabtem
      parb = (eta3p + amu3p)
      b3p   = -(7./15.)*brab - parb

!-25     Solution of quadratic for -(\tau N)^2
	
      anum = a3p
      bnum = -b3p
      cnum = -1.
      radical = SQRT(bnum**2 - 4*anum*cnum)
      qnum = (-1./2.)*(bnum + SIGN(radical,bnum))
      rootplus  = qnum/anum
      rootminus = cnum/qnum
!030325-27 Check that calculated solutions actually satisfy the quadratic equation.
!030327 Make acceptable error a small fraction of the root size.
      errorplus  = a3p*(rootplus**2)-b3p*(rootplus)-1.
      errorminus = a3p*(rootminus**2)-b3p*(rootminus)-1.
      errorproplus  = errorplus/rootplus
      errorprominus = errorminus/rootminus
      IF(MAX(ABS(errorproplus),ABS(errorprominus)).GT.errorbound) THEN
       if (mytid.eq.0) then
	WRITE(*,*) "PROBLEM: Error too great in calculating amtaun2 ."
	WRITE(*,*) "In quad_sal_pureconv:"
	WRITE(*,*) "-(\\tau N)^2=q/a \\equiv first root:",rootplus
	WRITE(*,*) "A''' (- \\tau N)^2)^2 - B'''((-\\tau N)^2) - 1 =", &
             errorplus
	WRITE(*,*) "proportional error= ",errorproplus
	WRITE(*,*) "-(\\tau N)^2=c/q \\equiv second root:",rootminus
	WRITE(*,*) "A''' (- \\tau N)^2)^2 - B'''((-\\tau N)^2) - 1 =", &
             errorminus
	WRITE(*,*) "proportional error= ",errorprominus
	WRITE(*,*) "Theoretically both roots should give zero."
	WRITE(*,*) "Maximum proportional error we permit =",errorbound
	WRITE(*,*) "a=",anum," b=",bnum," c=",cnum
	WRITE(*,*) "q=",qnum
	WRITE(*,*) "Program is stopping."
       endif
	STOP
      END IF
!******
!030326-0402Y     I choose "rootminus" because it approximates the solution with shear
!     at very negative Ri and because it bends downward as (N_H^2/N_d^2)
!     departs from zero which is consistent with an additional driving of
!     mixing due to double-diffusivity even in the unstable case.
!     IF CHANGE MODEL CONSTANTS *MUST* REVISIT THE CHOICE OF ROOT.
      amtaun2 = rootminus
      IF(ABS(anh2on2-ans2on2).LT.errorbound) THEN
        IF(ABS(rootminus-amtaun2_negri_0rid).GT.1.) THEN 	
       if (mytid.eq.0) then
	  WRITE(*,*) "Check if have the right root."
	  WRITE(*,*) "rootminus=",rootminus
	  WRITE(*,*) "amtaun2_negri_0rid=",amtaun2_negri_0rid
       endif
	  ier=1
        END IF
      END IF
!
      return                                          
      end                                            
                                                    



!_______________________________________________________________________










!990928 Revision to use mksmshss1 with option for exponential regularization
!	of mikesal2:
!990701-02 SUBMODULE intended to calculate Mikhail Dubovikov's model with salt
!	for the OGCM TURBULENCE MODULE in analogy with oursal2 for Cheng's. 
!       Introduce arguments for last guess and last row's zero guess for c_y.  
!************************************************************C
      subroutine mikesal2a(b1_arg,ri,rid,slq2,sm,sh,ss,c_y0,c_y00,&
!******************************C
!990701 When calculating for (Ri,Ri_d) index0 should be Ri and index1 Ri_d.
!       When calculating for (ra_r,\th_r) index{0,1} should be {ra_r,\th_r}. 
                          index0,index1, ri1)
!************************************************************C
!990701-02 Main subroutine in module of the same name.
!       Adapted with borrowings from oursal2.fs2 981016-1110 version
!       from smshssplot1_990630's extract mikel2sal.fs2 : 
!990614-17 For Dubovikov's new salinity model. Based on mikel2.
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
!990701 Use c_y \equiv (\tau_pv S)^2 only to maintain analogy with oursal2.
!**C    Section taken from oursal2:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     y=(tau*s)**2
!     tau=2*e/epsilon=b1*l/q
!     km=e*tau*sm=1/2*(b1*l)**2*s/y**(1/2)*sm
!     kh=e*tau*sh=1/2*(b1*l)**2*s/y**(1/2)*sh
!1016 ks=e*tau*ss=1/2*(b1*l)**2*s/y**(1/2)*ss
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!980723-1022  X = {M,H,C} .
!       Cheng above gives K_X = (1/2)((B_1*l)^2) (S/(((\tau S)**2)^(1/2))) S_X
!       The "old" y used above is (\tau S)^2.
!       The "new" y (c_y in the program) is (\tau_pv S)^2.
!       The program variable "slq2" is (S l/q)^2 = y (B_1)^(-2),
!       since \tau=B_1 l/q. (S l/q)^2 = (\tau \over \tau_pv)^2 c_y (B_1)^(-2) .
!       c_y = (S l/q)^2 * [(B_1)^2 * (\tau_pv \over \tau)^2] .
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!980615 From the printed notes Canuto gave me on 980601 have:
!        \tau_pv = {2 \over 5} \tau                     (B.1)
         PARAMETER(tpvot = 0.4D0)
!****** "tpv/tau" = 2/5
!**C 
!*****C
!990701 Take the value of the parameter B_1 from mcoeff subroutine.
         PARAMETER(b1=24.7D0)
!*****C
      real(r8) :: ri1
      real(r8) :: rid1,ric1,rit1,rr1
      data ifirst/0/
!990701 Change name of fctsal to fctysal just to emphasize the fact that it
!	calculates y \equiv (\tau S)^2 instead of c_y like oursal2's fct_sal.
      external fctysal

!990701 Allow B_1 to be passed back in case its needed.
!**C    Section taken from oursal2:
!981104 Give B_1 an alias needed for its fortran role as an output argument.
        b1_arg=b1
!*****C
!**C
!*****C
      ri1=ri
      rid1=rid
      rit1=0.5D0*(ri+rid)
      ric1=0.5D0*(ri-rid)
      rr1 = ric1/rit1
!990702 Make the error tolerance for the solver a millionth of the last guess
!	to try to avoid failure to solve problems which occur *at some points*
!	near the realizability limit. Restore iend to 10000 as in mikel2.
      eps=(1.e-6)*c_y0                                             
!990924 Increase iend to avoid trouble with exponential regularization model.
      iend=50000
!*****C 
!     rtwi finds the root of x=fct(x)                     
!     to make a guess at bst                     
!990701 Alteration adapted from oursal2.
      if((ifirst.eq.0).OR.(index0.EQ.0.AND.index1.EQ.0)) then
!990616 Use Ri_T=0,Ri_C=0  analytic value of  y = ((25/36)(2+\sqrt{67})^2 .
	rty_00 = (5.D0/6.D0)*(2.D0 + SQRT(67.D0))	!positive solution of quadratic.
        yst=rty_00**2   !at rit=ric=0
!*****C
!990701 Calculate the alternate "square shear number" variable c_y.
!     y \equiv (\tau S)^2 =  c_y (\tau_pv \over \tau)^{-2}
!     c_y = y (\tau_pv \over \tau)^2
      c_yst0 = yst*(tpvot**2)
      c_yst = c_yst0
!*****C
        ifirst=1
      else if(index0.EQ.0) then
        c_yst = c_y00
!     y \equiv (\tau S)^2 = (\tau_pv {\tau \over \tau_pv} S)^2
!                         = (\tau_pv S)^2 (\tau \over \tau_pv)^2
!			  = c_y (\tau_pv \over \tau)^{-2}
        yst=c_yst/(tpvot**2)                                
      ELSE
        c_yst = c_y0
        yst=c_yst/(tpvot**2)                                
      endif
!*****C
      call rtwi(y,val,fctysal,yst,eps,iend,ier,ttot,tcot,tctot,tptot,tpcot,tpvot,ri,rit1,ric1,rit,ric)
!      write(*,*) ri,rid,ier,y,val
!     Failure when solver has not found a root.
      if(ier.ne.0) then
       if (mytid.eq.0) then
         write(*,*) "rtwi call problem, ier=",ier
!        stop
	 WRITE(*,*) "Could not solve turbulence model at:"
!990701-02 Elaborate error message when solver does not find root.
	    WRITE(*,*) "first index=",index0,"second index=",index1
!**C    Section taken from oursal2:
!981022-23-30 Make error message more specific.
            WRITE(*,*) "In mikesal2a subroutine"
            WRITE(*,*) "c_y00=",c_y00," c_y0=",c_y0
            WRITE(*,*) "ri=",ri,"       rid=",rid
            WRITE(*,*) "rit=",rit1,"     ric=",ric1
            WRITE(*,*) "Initial guess for rtwi c_yst=",c_yst
!*****C
!**C
	    WRITE(*,*) " "
!*****C
       endif
      endif
      slq2=y/(b1*b1)
!990701 Calculate the alternate "square shear number" variable c_y.
!     y \equiv (\tau S)^2 =  c_y (\tau_pv \over \tau)^{-2}
!     c_y = y (\tau_pv \over \tau)^2
      c_y = y*(tpvot**2)
!990701 Storage of c_y values for calling routine.
!	Adapted from oursal2.
!981016 Store value of c_y for future guesses.
         IF(c_y.GE.0) THEN
            c_y0=c_y
         ELSE
!981022 Turbulence model becomes unphysical for c_y negative. 
!981014-16  Realizability for negative Ri
           IF(ri.LT.0) THEN
       if (mytid.eq.0) then
             WRITE(*,*) "c_y negative at negative Ri"
             WRITE(*,*) "Ri=",ri,"      c_y=",c_y
	     WRITE(*,*) "first index=",index0,"second index=",index1
             WRITE(*,*) "Unstable realizability limit unexpected:"
             WRITE(*,*) "stopping in mikesal2a."
       endif
             STOP
           END IF
         END IF
         IF(index0.EQ.0) c_y00=c_y
         IF((index0.EQ.0).AND.(index1.EQ.0).AND.  &
      (ABS(c_y - c_yst0).GT.1.D-6)) THEN
       if (mytid.eq.0) then
           WRITE(*,*) "Inconsistency in neutral value of c_y"
           WRITE(*,*) "Value used =",c_yst0
           WRITE(*,*) "Value calculated =",c_y
           WRITE(*,*) "Program stopping in mikesal2a"
       endif
           STOP
         END IF
!*****C 
!*****C
!990614 Introduce dimensionless functions for \tau times the T and S parts of N^2.
      x=rit1*y
      z=ric1*y
      call mksmshss1(x,y,z,sm,sh,ss)
      return
      end



      function fctysal(y,ttot,tcot,tctot,tptot,tpcot,tpvot,ri,rit1,ric1,rit,ric)
!990614-15 Function for y \equiv (\tau S)^2 from the Production=Dissipation equation.
!980501 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
      real(r8) :: ri
      real(r8) :: rid1,ric1,rit1,rri,rit,ric
      real(r8) :: ttot,tcot,tctot,tptot,tpcot,tpvot
      x=rit1*y
      z=ric1*y
      call mksmshss1(x,y,z,sm,sh,ss)
!990614 Use the traditional Production=Dissipation equation:
!	y(S_M - Ri_T S_H - Ri_C S_C) = 2            
      fctysal=2.D0/(sm-rit1*sh-ric1*ss)
      return                                          
      end                                            


                                                                        
      subroutine mksmshss1(x,y,z,sm,sh,ss)
!991107 Write-out of regularization flag "ifnew" introduced.
!990922-24 Subroutine with option to calculate the new version of Dubovikov's
!	salinity model with exponential regularization,
!	using my 990923 correction of the Phi formula cleared on phone with mike
!       to drop the 1 term in the exponential, so it=1st order old expression,
!	altered from mksmshss0:
!990630-0702 Subroutine to calculate only the new 990621 version of
!	Dubovikov's salinity model. Stripped from mksmshss:
!990621(After Summer Solstice) Allow still another alternative 
!      for a new formulation Canuto and Dubovikov come up with today.
!990621 Allow for Canuto additional alternative where x_1**4 is used AND
!	add a term {F_C}^{-2} \eta_\theta \eta_sigma to {\sigma_m}^{-1}.
!990618 Allow alternate calculation of Dubovikov S_M and S_H
!	because this evening he admits that there may be an error:
!	instead of x_1**2 multiplying Ri_C and Ri_T, x_1**4 may be correct.
!990614-17 Calculate S_M,S_H,S_S for new Dubovikov model with salt using
!       Dubovikov's 990614 handwritten notes with F_i,F\twiddle_i = 1 for P=\ep.
!980717 Make double precision to conform to calling cctmix routine.
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!      implicit real*8 (a-h,o-z)
      REAL(r8) phi
!
!
!990924 Set ifnew = 0 for Summer Solstice '99 version with 1/1+x regularization.
!	    ifnew = 1 for 990923 corrected new version with e^-x regularization.
!
      PARAMETER(ifnew=0)
!
      parameter(x1min=15.d0/8.d0)
!990617 Artificially limit y to be at least a small positive number to avoid
!     problems due to the square root of y and the reciprocals.
      PARAMETER(calculation_epsilon=1.D-60)
!990621-30 Use the new formulation introducing the variable 
!       \phi = 1 + 0.08 {x_DB}^2 Ri
!990924 For option ifnew=1: \phi = e^{-0.08 {x_DB}^2 Ri}

!991107 INTRODUCE WRITEOUT TO INDICATE WHETHER EXPONENTIAL REGULARIZATION USED.
	IF(ifrecall.EQ.0) THEN
       if (mytid.eq.0) then
	  WRITE(*,*) " "
	  WRITE(*,*)"************************************************"
	  WRITE(*,*) " "
	  WRITE(*,*)"Regularization used for Dubovikov salinity model"
	  WRITE(*,*) "ifnew=",ifnew
	  IF (ifnew.EQ.0) THEN
	    WRITE(*,*) "1/1+x regularization"
	  ELSE IF(ifnew.EQ.1) THEN
	    WRITE(*,*) "e^-x regularization"
	  END IF  
	  WRITE(*,*) " "
	  WRITE(*,*)"************************************************"
	  WRITE(*,*) " "
        endif
	END IF
!*****************************************************************************C

      rit=x/y
      ric=z/y
      y = MAX(y,calculation_epsilon)
      x = rit*y
      z = ric*y
!*****C
      x1=0.5d0*sqrt(y)
      if(x1.lt.x1min) then
         fc_twid=5.d0/9.d0
      else
         fc_twid=1.d0-5.d0/(3.d0*x1)+25.d0/(16.d0*x1**2)
      endif
!990615 Total Richardson number, Ri = Ri_T + Ri_C .
      ri = rit + ric
!990614 From Equation (65b) of 990610 paper excerpt: F\twiddle_C = {5 \over 9} F_C  .
!	F_C = {9 \over 5} F\twiddle_C
      fc = 1.8D0*fc_twid
!990614 From Dubovikov's 990614 handwritten notes with F_i,F\twiddel_i = 1:
!	{\sigma_m}^{-1} = 1 + {F_C}^{-1} (\eta_\theta + \eta_s)
!	\eta_\s(Ri_T,Ri_C) = eta_\theta(Ri_C,Ri_T)       : Tem,Sal symmetry condition 
!990924 In ifnew =1 case use formula received 990913 : 
!	sigma_m = exp[-(\eta_\theta + \eta_\sigma) {F_C}^{-1}]
!990621 Dubovikov's newer equation for \psi:
!	\psi^{-1} = 
!         1 + 0.16 {x_1}^2 Ri + 0.013 {x_1}^4 Ri_T Ri_C \Phi
!       \Phi^{-1} = 1 + 0.08 x^2 Ri    
!990924 In ifnew=1 case use formula received 990913 ammended for \Phi 990923:
!	\psi_1 = exp[-0.16 {x_DB}^2 (1 - R_\rho) Ri_DB + 
!                0.013 {x_DB}^4 R_\rho Ri^2 \Phi]
!	Ri_DB \equiv Ri_T   	R_\rho \equiv -Ri_C/Ri_T
!	\psi_1 = e^{-0.16 {x_1}^2 (1 + Ri_C/Ri_T)Ri_T - 
!                    0.013 {x_1}^4 (Ri_C/Ri_T) {Ri_T}^2 \Phi}
!	\psi_1 = e^{-0.16 {x_1}^2 (Ri_T + Ri_C) - 
!                    0.013 {x_1}^4 Ri_C Ri_T \Phi}
!	\psi_1 = e^{-0.16 {x_1}^2 Ri - 0.013 {x_1}^4 Ri_C Ri_T \Phi}
!       \Phi = exp[ - 0.08 {x_DB}^2 (1 - R_\rho) Ri_DB ]    
!       \Phi = e^{ - 0.08 {x_1}^2 (1 + Ri_C/Ri_T) Ri_T }    
!       \Phi = e^{ - 0.08 {x_1}^2 (Ri_T + Ri_C)}    
!       \Phi = e^{ - 0.08 {x_1}^2 Ri}    
	IF(ifnew.EQ.0) THEN
	  phi = 1.D0/(1.D0 + 0.08*x1**2*ri)
        psi = 1.D0/ &
            (  1.D0 & 
            + 0.16D0*x1**2*ri &
             + 0.013D0*x1**4*rit*ric*phi )
	ELSE IF(ifnew.EQ.1) THEN 	!990923 version
	  phi =  EXP(-0.08*x1**2*ri)
          psi = EXP (-0.16D0*x1**2*ri  - 0.013D0*x1**4*rit*ric*phi)
	END IF
!*****C
!990621 Alternate equation 
      eta_t = 0.05D0*x1**2*rit*(1.D0 + (1.D0 +0.16D0*ric*x1**2)*psi) 
      eta_s = 0.05D0*x1**2*ric*(1.D0 + (1.D0 +0.16D0*rit*x1**2)*psi) 
      eta_sum = eta_t + eta_s
      IF(ifnew.EQ.0) THEN
	sigma_m = 1.D0/(1.D0 + eta_sum/fc)
      ELSE IF(ifnew.EQ.1) THEN		!990923 version
	sigma_m = EXP(-eta_sum/fc)
      END IF
!*****C
!990614 Equation (65a) of 990610 paper excerpt: 
!       S_M = {1 \over 25} F_1 F\twiddle_C \sigma_m
!	Set F_1 = 1 .
      sm=1.d0/25.d0*fc_twid*sigma_m
!990621 S_H = 0.056[1 + 0.08 x1**2 Ri_C - 0.006 x1**4 Ri_C Ri] \psi	
        sh=0.056d0*(1.D0 + 0.08D0*x1**2*ric*phi)*psi
        ss=0.056d0*(1.D0 + 0.08D0*x1**2*rit*phi)*psi
!*****C
!990924 Add Concentration Smagorinsky-like function for 990923 version.
!	'S_\tau = 0.056[1 + 0.08 {x_DB}^2 {Ri_DB} (1 - R_\rho) \Phi] \psi_1'
!       S_\tau = 0.056 (1 + 0.08 {x_1}^2 Ri_T (1 + Ri_C/Ri_T) \Phi] \psi_1
!       S_\tau = 0.056 (1 + 0.08 {x_1}^2 Ri \Phi] \psi_1
      IF(ifnew.EQ.1) THEN
        st=0.056d0*(1.D0 + 0.08D0*x1**2*ri*phi)*psi
      END IF
	

      return
      end



!-----------------------------------------------------------------------
!     end of parameter-free turbulence model subroutines (nmodel=2)
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
!     end of parameter-free turbulence model subroutines (nmodel=2)
!-----------------------------------------------------------------------

!030401-02Y Adaptation to the NCAR CSM Ocean Model of my 1D interpolation routine 
!	written for HYCOM [See NBp.030401-2to3]. I must make variables double precision.
	subroutine interp1d(x, &
                    x_1,slq2_1,sm_1,sh_1,ss_1, &
                    slq2,sm,sh,ss, &
                    ixmax,m,m0,delta)
      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!        IMPLICIT REAL*8 (a-h,o-z)

! --- 6-990112-C030327	Subroutine for a ONE VARIABLE modular interpolation calculation.
!	x is the independent variable in the turbulence model calculation.
!
! --- 1D input array with table spacing: 	x_1
! --- table limit value:			ixmax
! --- 1D input arrays with table values:	slq2_1,sm_1,sh_1,ss_1
! --- Output interpolated values:		slq2,sm,sh,ss

!
	INTEGER m,m0
	DIMENSION x_1(-m:m)
        integer ixmax
	DIMENSION slq2_1(-m:m), sm_1(-m:m),sh_1(-m:m),ss_1(-m:m)

!     print *,'start interp1d'
! --- 030326 Take values off the edge of the table back to the table limits.
	IF(x.GT.x_1(ixmax)) THEN
	    x = x_1(ixmax)
	ELSE IF(x.LT.x_1(-m)) THEN
	    x = x_1(-m)
	END IF
!
! --- Interpolate points within the table range.
! --- Table index ranges from -m to m with equal spacing for -m0 to m0.
!
	IF(ABS(x).LT.x_1(m0)) THEN
!
! --- Find Interpolation points in the equally spaced x part of the table.
#if (defined D_PRECISION)
          lx1 = INT(x/delta)+NINT(SIGN(DFLOAT(1),x))
#else
          lx1 = INT(x/delta)+NINT(SIGN(FLOAT(1),x))
#endif
	ELSE
!
! --- Find Interpolation points in unequally spaced x part of the table.
	  DO l=m0,m
!
! --- Keep stepping until pass input value.
#if (defined D_PRECISION)
	     IF(ABS(x).LT.ABS(x_1(NINT(SIGN(DFLOAT(l),x))))) THEN
#else
	     IF(ABS(x).LT.ABS(x_1(NINT(SIGN(FLOAT(l),x))))) THEN
#endif
!
! --- Cover both positive and negative index cases.
#if (defined D_PRECISION)
	       lx1 = NINT(SIGN(DFLOAT(l),x))
#else
	       lx1 = NINT(SIGN(FLOAT(l),x))
#endif
!
	       GO TO 250
	     END IF
	  END DO
!
! --- Special case where have a value which falls at the limit of the table.
#if (defined D_PRECISION)
	  lx0 = NINT(SIGN(DFLOAT(m),x))
#else
	  lx0 = NINT(SIGN(FLOAT(m),x))
#endif
	  lx1 = lx0
	  GO TO 252
!
	END IF
  250	CONTINUE
!
! --- Make lx0 one less or greater than lx1 according to sgn(x).
#if (defined D_PRECISION)
        lx0 = lx1 - NINT(SIGN(DFLOAT(1),x)) 
#else
        lx0 = lx1 - NINT(SIGN(FLOAT(1),x)) 
#endif
!
        IF(x.EQ.0.) lx1 = 1
  252   CONTINUE
!
! --- Check that the x value falls within the interpolation interval.
	IF((x.GT.0..AND.  &
      (x.LT.x_1(lx0).OR.x.GT.x_1(lx1))).OR.  &
	   (x.LT.0..AND.  &
      (x.GT.x_1(lx0).OR.x.LT.x_1(lx1)))) THEN
       if (mytid.eq.0) then
	   WRITE(*,*) "x is outside interpolation range in interp1d."
	   WRITE(*,*) "x=  ",x,"lx0= ",lx0,"lx1= ",lx1
	   WRITE(*,*) "x_1(lx0)=  ",x_1(lx0), &
                "   x_1(lx1)= ",x_1(lx1)
	   WRITE(*,*) "Program is stopping."
       endif
	   STOP
	END IF
!
!
!
!
!
  270	CONTINUE
!
!
!
! --- Interpolate turbulence fields.
! --- Introduce table spacing variables.
	deltaxta = x_1(lx1) - x_1(lx0)
	deltax = x - x_1(lx0)

! --- slq2
! --- Set delta field to zero in special cases falling at limit of the table. 
	IF(lx1.EQ.lx0) THEN
	  dslq2_x = 0.
	ELSE
	  dslq2_x = (slq2_1(lx1) - slq2_1(lx0))/ &
                deltaxta
	END IF
	slq2	 = slq2_1(lx0) + &
             dslq2_x*deltax
!
! --- sm
	IF(lx1.EQ.lx0) THEN
	  dsm_x   = 0.
	ELSE
	  dsm_x = (sm_1(lx1) - sm_1(lx0))/ &
              deltaxta
	END IF
	sm     = sm_1(lx0) + &
           dsm_x*deltax
!
! --- sh
	IF(lx1.EQ.lx0) THEN
	  dsh_x   = 0.
	ELSE
	  dsh_x = (sh_1(lx1) - sh_1(lx0))/ &
              deltaxta
	END IF
	sh     = sh_1(lx0) + &
           dsh_x*deltax
!
! --- ss
	IF(lx1.EQ.lx0) THEN
	  dss_x   = 0.
	ELSE
	  dss_x = (ss_1(lx1) - ss_1(lx0))/ &
              deltaxta
	END IF
	ss     = ss_1(lx0) + &
            dss_x*deltax
!
!
!     print *,'finish interp1d'
	RETURN
	END

!-----------------------------------------------------------
!030401 Adaptation to the NCAR CSM Ocean Model of my 1D interpolation routine
!       written for HYCOM [See NBp.030401-2to3]. 
	SUBROUTINE interp1d_expabs(x, &
                    x_1,slq2_1,sm_1,sh_1,ss_1, &
                    slq2,sm,sh,ss, &
                    ixmax,m,m0,delta,rat)
!030129-030328 Subroutine for a faster interpolation calculation in the ifexpabstable=1 case.
! --- 6-990112-C030326	Subroutine for a ONE VARIABLE modular interpolation calculation.
!	x is the independent variable in the turbulence model calculation.
!
! --- 1D input array with table spacing: 	x_1
! --- table limit value:			ixmax
! --- 1D input arrays with table values:	slq2_1,sm_1,sh_1,ss_1
! --- Output interpolated values:		slq2,sm,sh,ss
!01107yXI Input value of ratio between entries  rat

      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!        IMPLICIT REAL*8 (a-h,o-z)
	INTEGER ixmax
	INTEGER m,m0
	DIMENSION x_1(-m:m)
	DIMENSION slq2_1(-m:m),  sm_1(-m:m),sh_1(-m:m),ss_1(-m:m)

!030326 Take values off the edge of the table back to the table.
!     print *,'start interp1d_expabs'
	IF(x.GT.x_1(ixmax)) THEN
	    x = x_1(m)
	ELSE IF(x.LT.x_1(-m)) THEN
	    x = x_1(-m)
	END IF
!*****C
!981019 Interpolate points within the table range.
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
	IF(ABS(x).LT.x_1(m0)) THEN
!981019-27 Find Interpolation points in the equally spaced x part of the table.
#if (defined D_PRECISION)
	  lx1 = INT(x/delta)+NINT(SIGN(DFLOAT(1),x))
#else
	  lx1 = INT(x/delta)+NINT(SIGN(FLOAT(1),x))
#endif
!011107-08yXI Find Interpolation points in exponential absolute value spaced x part of the table.
	ELSE IF((ABS(x)).GE.(x_1(m))) &
   THEN
!981103 Special case where have a value which falls at the limit of the table.
#if (defined D_PRECISION)
	  lx0 = NINT(SIGN(DFLOAT(m),x))
#else
	  lx0 = NINT(SIGN(FLOAT(m),x))
#endif
	  lx1 = lx0
	  GO TO 252
!*****C
	ELSE
#if (defined D_PRECISION)
	  tabindx = SIGN( DFLOAT(m0) + ((LOG(ABS(x)) - LOG(x_1(m0)))/LOG(rat)), x)
#else
	  tabindx = SIGN( FLOAT(m0) + ((LOG(ABS(x)) - LOG(x_1(m0)))/LOG(rat)), x)
#endif
#if (defined D_PRECISION)
	  lx1 = INT(tabindx)+NINT(SIGN(DFLOAT(1),x))
#else
	  lx1 = INT(tabindx)+NINT(SIGN(FLOAT(1),x))
#endif
!yXI
	END IF
!011108yXI It is conceivable that rounding errors may in borderline cases 
!	   throw the calculated table indices for x off by one.
!	   Check and allow moving one to either side to take care of this.
 	IF((ABS(x_1(lx1))).LT.(ABS(x))) THEN
	  lx1 = lx1 + SIGN(1,lx1)
	ELSE IF((ABS(x_1(lx1-SIGN(1,lx1)))).GT.(ABS(x))) THEN
	  lx1 = lx1 - SIGN(1,lx1)
	END IF
!yXI
  250	CONTINUE
!981019-27 Make lx0 one less or greater than lx1 according to sgn(x).
#if (defined D_PRECISION)
        lx0 = lx1 - NINT(SIGN(DFLOAT(1),x)) 
#else
        lx0 = lx1 - NINT(SIGN(FLOAT(1),x)) 
#endif
!*****C
        IF(x.EQ.0.D0) lx1 = 1
  252   CONTINUE
!981019-28 Check that the x value falls within the interpolation interval.
	IF((x.GT.0.D0.AND.  &
      (x.LT.x_1(lx0).OR.x.GT.x_1(lx1))).OR.  &
	   (x.LT.0.D0.AND.  &
      (x.GT.x_1(lx0).OR.x.LT.x_1(lx1)))) THEN
       if (mytid.eq.0) then
	   WRITE(*,*) &
      "x is outside interpolation range in interp1d_expabs."
	   WRITE(*,*) "delta=",delta
	   WRITE(*,*) "m0=",m0," m=",m," rat=",rat
	   WRITE(*,*) "x=  ",x," lx0= ",lx0," lx1= ",lx1
	   WRITE(*,*) "x_1(lx0)=  ",x_1(lx0), &
                "   x_1(lx1)= ",x_1(lx1)
	   WRITE(*,*) "Program is stopping."
        endif
	   STOP
	END IF
!*****C
!981019-27 Interpolate turbulence fields.
!981027-990112 Introduce table spacing variables.
	deltaxta = x_1(lx1) - x_1(lx0)
	deltax = x - x_1(lx0)
!	slq2
!981103 Set delta field to zero in special cases falling at limit of the table. 
	IF(lx1.EQ.lx0) THEN
	  dslq2_x = 0.D0
	ELSE
	  dslq2_x = (slq2_1(lx1) - slq2_1(lx0))/ &
                deltaxta
	END IF
	slq2	 = slq2_1(lx0) + &
             dslq2_x*deltax
!	sm
	IF(lx1.EQ.lx0) THEN
	  dsm_x   = 0.D0
	ELSE
	  dsm_x = (sm_1(lx1) - sm_1(lx0))/ &
              deltaxta
	END IF
	sm     = sm_1(lx0) + &
           dsm_x*deltax
!	sh
	IF(lx1.EQ.lx0) THEN
	  dsh_x   = 0.D0
	ELSE
	  dsh_x = (sh_1(lx1) - sh_1(lx0))/ &
              deltaxta
	END IF
	sh     = sh_1(lx0) + &
           dsh_x*deltax
!	ss
	IF(lx1.EQ.lx0) THEN
	  dss_x   = 0.D0
	ELSE
	  dss_x = (ss_1(lx1) - ss_1(lx0))/ &
              deltaxta
	END IF
	ss     = ss_1(lx0) + &
           dss_x*deltax
!*****C


!     print *,'finish interp1d_expabs'
	RETURN
	END

!-----------------------------------------------------------------------------
	SUBROUTINE interp2d(ri,rid, &
                    ri_1,rid_1,slq2_2,sm_2,sh_2,ss_2, &
                    slq2,sm,sh,ss, &
                    irimax,m,m0,delta)

!981016-990112	Subroutine for a modular interpolation calculation.
!*
!	1D input arrays with table spacing: 	ri_1,rid_1
!	1D input array with table limits:	irimax
!	2D input arrays with table values:	slq2_2,sm_2,sh_2,ss_2
!	Output interpolated values:		slq2,sm,sh,ss

      use precision_mod
      use param_mod, only: mytid
      implicit real(r8) (a-h,o-z)
!	IMPLICIT REAL*8 (a-h,o-z)
	DIMENSION ri_1(-m:m),rid_1(-m:m),irimax(-m:m)
	DIMENSION slq2_2(-m:m,-m:m), &
            sm_2(-m:m,-m:m),sh_2(-m:m,-m:m),ss_2(-m:m,-m:m)

! 	Take values off the edge of the table back to the table on radial lines.
!981103 Must use ratio of Ri's taken *before* the cut-off has taken place.
	IF(ri.GT.ri_1(m)) THEN
	  IF(ABS(rid).LE.ri) THEN
	    rid = ri_1(m)*(rid/ri)
	    ri  = ri_1(m)
	  ELSE IF(rid.GT.ri) THEN
	    ri  = rid_1(m)*(ri/rid)
	    rid = rid_1(m)
	  ELSE IF(rid.LT.-ri) THEN
	    ri  = rid_1(-m)*(ri/rid)
	    rid = rid_1(-m)
	  END IF
	ELSE IF(ri.LT.ri_1(-m)) THEN
	  IF(ABS(rid).LE.-ri) THEN
	    rid = ri_1(-m)*(rid/ri)
	    ri  = ri_1(-m)
	  ELSE IF(rid.GT.-ri) THEN
	    ri  = rid_1(m)*(ri/rid)
	    rid = rid_1(m)
	  ELSE IF(rid.LT.ri) THEN
	    ri  = rid_1(-m)*(ri/rid)
	    rid = rid_1(-m)
	  END IF
	ELSE IF(rid.GT.rid_1(m)) THEN
	    ri  = rid_1(m)*(ri/rid)
	    rid = rid_1(m)
	ELSE IF(rid.LT.rid_1(-m)) THEN
	    ri  = rid_1(-m)*(ri/rid)
	    rid = rid_1(-m)
	END IF
!*****C
!981019 Interpolate points within the table range.
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
	IF(ABS(rid).LT.rid_1(m0)) THEN
!981019-27 Find Interpolation points in the equally spaced Ri_d part of the table.
#if (defined D_PRECISION)
	  lrid1 = INT(rid/delta)+NINT(SIGN(DFLOAT(1),rid))
#else
	  lrid1 = INT(rid/delta)+NINT(SIGN(FLOAT(1),rid))
#endif
	ELSE
!981019 Find Interpolation points in unequally spaced Ri_d part of the table.
	  DO l=m0,m
!981027 Keep stepping until pass input value.
#if (defined D_PRECISION)
	     IF(ABS(rid).LT.ABS(rid_1(NINT(SIGN(DFLOAT(l),rid))))) THEN
#else
	     IF(ABS(rid).LT.ABS(rid_1(NINT(SIGN(FLOAT(l),rid))))) THEN
#endif
!981027 Cover both positive and negative index cases.
#if (defined D_PRECISION)
	       lrid1 = NINT(SIGN(DFLOAT(l),rid))
#else
	       lrid1 = NINT(SIGN(FLOAT(l),rid))
#endif
!*****C
	       GO TO 250
	     END IF
	  END DO
!981103 Special case where have a value which falls at the limit of the table.
#if (defined D_PRECISION)
	  lrid0 = NINT(SIGN(DFLOAT(m),rid))
#else
	  lrid0 = NINT(SIGN(FLOAT(m),rid))
#endif
	  lrid1 = lrid0
	  GO TO 252
!*****C
	END IF
  250	CONTINUE
!981019-27 Make lrid0 one less or greater than lrid1 according to sgn(rid).
#if (defined D_PRECISION)
        lrid0 = lrid1 - NINT(SIGN(DFLOAT(1),rid)) 
#else
        lrid0 = lrid1 - NINT(SIGN(FLOAT(1),rid)) 
#endif
!*****C
        IF(rid.EQ.0.D0) lrid1 = 1
  252   CONTINUE
!981019-27 Check that the Ri_d value falls within the interpolation interval.
	IF((rid.GT.0.D0.AND.  &
      (rid.LT.rid_1(lrid0).OR.rid.GT.rid_1(lrid1))).OR.  &
	   (rid.LT.0.D0.AND.  &
      (rid.GT.rid_1(lrid0).OR.rid.LT.rid_1(lrid1)))) THEN
       if (mytid.eq.0) then
	   WRITE(*,*) "Ri_d is outside interpolation range in interp2d."
	   WRITE(*,*) "rid=  ",rid,"lrid0= ",lrid0,"lrid1= ",lrid1
	   WRITE(*,*) "rid_1(lrid0)=  ",rid_1(lrid0), &
                "   rid_1(lrid1)= ",rid_1(lrid1)
	   WRITE(*,*) "Program is stopping."
        endif
	   STOP
	END IF
!*****C
!C981022	Artificially reduce Ri if it threatens to surpass Ri_max(Ri_d).
!C	This is to conform to the 1D table's realizability limit treatment. 
!	IF(ri.GT.MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))) THEN
!	  ri = MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))
!	END IF
!C*****C
!981110 Set turbulence to zero if Ri threatens to surpass the realizability limit.
        IF(ri.GT.MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))) THEN
          slq2=0.D0
          sm = 0.D0
          sh = 0.D0
          ss = 0.D0
          RETURN
        END IF
!*****C
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
	IF(ABS(ri).LT.ri_1(m0)) THEN
!981019-27 Find Interpolation points in the equally spaced Ri part of the table.
#if (defined D_PRECISION)
	  lri1 = INT(ri/delta)+NINT(SIGN(DFLOAT(1),ri)) 
#else
	  lri1 = INT(ri/delta)+NINT(SIGN(FLOAT(1),ri)) 
#endif
	ELSE
!981019 Find Interpolation points in unequally spaced Ri part of the table.
	  DO l=m0,m
!981027 Keep stepping until pass input value.
#if (defined D_PRECISION)
	     IF(ABS(ri).LT.ABS(ri_1(NINT(SIGN(DFLOAT(l),ri))))) THEN
#else
	     IF(ABS(ri).LT.ABS(ri_1(NINT(SIGN(FLOAT(l),ri))))) THEN
#endif
!981027 Cover both positive and negative index cases.
#if (defined D_PRECISION)
	       lri1 = NINT(SIGN(DFLOAT(l),ri))
#else
	       lri1 = NINT(SIGN(FLOAT(l),ri))
#endif
!*****C
	       GO TO 270
	     END IF
	  END DO
!981103 Special case where have a value which falls at the limit of the table.
#if (defined D_PRECISION)
	  lri0 = NINT(SIGN(DFLOAT(m),ri))
#else
	  lri0 = NINT(SIGN(FLOAT(m),ri))
#endif
	  lri1 = lri0
	  GO TO 272
!*****C
  270	CONTINUE
	END IF
!981019-27 Make lri0 one less or greater than lri1 according to sgn(ri).
#if (defined D_PRECISION)
        lri0 = lri1 - NINT(SIGN(DFLOAT(1),ri)) 
#else
        lri0 = lri1 - NINT(SIGN(FLOAT(1),ri)) 
#endif
!*****C
        IF(ri.EQ.0.D0) lri1 = 1
  272	CONTINUE
!981019-27 Check that the Ri_d value falls within the interpolation interval.
	IF((ri.GT.0.D0.AND.(ri.LT.ri_1(lri0).OR.ri.GT.ri_1(lri1))).OR.  &
	   (ri.LT.0.D0.AND.(ri.GT.ri_1(lri0).OR.ri.LT.ri_1(lri1)))) THEN
       if (mytid.eq.0) then
	   WRITE(*,*) "Ri is outside interpolation range in interp2d."
	   WRITE(*,*) "ri=  ",ri,"lri0= ",lri0,"lri1= ",lri1
	   WRITE(*,*) "ri_1(lri0)=  ",ri_1(lri0), &
                "   ri_1(lri1)= ",ri_1(lri1)
	   WRITE(*,*) "Program is stopping."
       endif
	   STOP
	END IF
!*****C
!981019-27 Interpolate turbulence fields.
!981027-990112 Introduce table spacing variables.
	deltaridta = rid_1(lrid1) - rid_1(lrid0)
	deltarita  = ri_1(lri1)  - ri_1(lri0)
	deltarid = rid - rid_1(lrid0)
	deltari  = ri - ri_1(lri0)
!	slq2
!981103 Set delta field to zero in special cases falling at limit of the table. 
	IF(lrid1.EQ.lrid0) THEN
	  dslq2_rid = 0.D0
	ELSE
	  dslq2_rid = (slq2_2(lri0,lrid1) - slq2_2(lri0,lrid0))/ &
                deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dslq2_ri  = 0.D0
	ELSE
	  dslq2_ri = (slq2_2(lri1,lrid0) - slq2_2(lri0,lrid0))/ &
               deltarita
	END IF
	slq2	 = slq2_2(lri0,lrid0) + &
            dslq2_ri*deltari + dslq2_rid*deltarid
!	sm
	IF(lrid1.EQ.lrid0) THEN
	  dsm_rid   = 0.D0
	ELSE
	  dsm_rid = (sm_2(lri0,lrid1) - sm_2(lri0,lrid0))/ &
              deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dsm_ri    = 0.D0
	ELSE
	  dsm_ri = (sm_2(lri1,lrid0) - sm_2(lri0,lrid0))/ &
             deltarita
	END IF
	sm     = sm_2(lri0,lrid0) + &
            dsm_ri*deltari + dsm_rid*deltarid
!	sh
	IF(lrid1.EQ.lrid0) THEN
	  dsh_rid   = 0.D0
	ELSE
	  dsh_rid = (sh_2(lri0,lrid1) - sh_2(lri0,lrid0))/ &
              deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dsh_ri    = 0.D0
	ELSE
	  dsh_ri = (sh_2(lri1,lrid0) - sh_2(lri0,lrid0))/ &
              deltarita
	END IF
	sh     = sh_2(lri0,lrid0) + &
            dsh_ri*deltari + dsh_rid*deltarid
!	ss
	IF(lrid1.EQ.lrid0) THEN
	  dss_rid   = 0.D0
	ELSE
	  dss_rid = (ss_2(lri0,lrid1) - ss_2(lri0,lrid0))/ &
              deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dss_ri    = 0.D0
	ELSE
	  dss_ri = (ss_2(lri1,lrid0) - ss_2(lri0,lrid0))/ &
             deltarita
	END IF
!*****C
	ss     = ss_2(lri0,lrid0) + &
            dss_ri*deltari + dss_rid*deltarid
!*****C


	RETURN
	END


	SUBROUTINE interp2d_expabs(ri,rid,slq2,sm,sh,ss)

!030123z Remove the "sys_flush" from the model E interpolation routine that avoids
!	 in the exponential absolute value case the super-time-consuming stepping 
!	 through the nonlinear part of the table by estimating the table index near
!	 the value to be interpolated in order to use this more efficient routine
!	 in my traditional NCOM runs.
!011107-08yXI Subroutine for an interpolation calculation for a table whose nonlinear part 
!yXI          is exponential in the absolute value of the index. Based on:
!981016-990112	Subroutine for a modular interpolation calculation.
!*
!01107yXI Input value of ratio between entries  rat
!	1D input arrays with table spacing: 	ri_1,rid_1
!	1D input array with table limits:	irimax
!	2D input arrays with table values:	slq2_2,sm_2,sh_2,ss_2
!	Output interpolated values:		slq2,sm,sh,ss

      use precision_mod
      use param_mod, only: mytid
      use canuto_mod
      implicit real(r8) (a-h,o-z)
!	IMPLICIT REAL*8 (a-h,o-z)
! 	Take values off the edge of the table back to the table on radial lines.
!981103 Must use ratio of Ri's taken *before* the cut-off has taken place.
	IF(ri.GT.rib(mt)) THEN
	  IF(ABS(rid).LE.ri) THEN
	    rid = rib(mt)*(rid/ri)
	    ri  = rib(mt)
	  ELSE IF(rid.GT.ri) THEN
	    ri  = ridb(mt)*(ri/rid)
	    rid = ridb(mt)
	  ELSE IF(rid.LT.-ri) THEN
	    ri  = ridb(-mt)*(ri/rid)
	    rid = ridb(-mt)
	  END IF
	ELSE IF(ri.LT.rib(-mt)) THEN
	  IF(ABS(rid).LE.-ri) THEN
	    rid = rib(-mt)*(rid/ri)
	    ri  = rib(-mt)
	  ELSE IF(rid.GT.-ri) THEN
	    ri  = ridb(mt)*(ri/rid)
	    rid = ridb(mt)
	  ELSE IF(rid.LT.ri) THEN
	    ri  = ridb(-mt)*(ri/rid)
	    rid = ridb(-mt)
	  END IF
	ELSE IF(rid.GT.ridb(mt)) THEN
	    ri  = ridb(mt)*(ri/rid)
	    rid = ridb(mt)
	ELSE IF(rid.LT.ridb(-mt)) THEN
	    ri  = ridb(-mt)*(ri/rid)
	    rid = ridb(-mt)
	END IF
!*****C
!981019 Interpolate points within the table range.
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
	IF(ABS(rid).LT.ridb(mt0)) THEN
!981019-27 Find Interpolation points in the equally spaced Ri_d part of the table.
#if (defined D_PRECISION)
	  lrid1 = INT(rid/dri)+NINT(SIGN(DFLOAT(1),rid))
#else
	  lrid1 = INT(rid/dri)+NINT(SIGN(FLOAT(1),rid))
#endif
!011107-08yXI Find Interpolation points in exponential absolute value spaced Ri_d part of the table.
	ELSE IF((ABS(rid)).GE.(ridb(mt)))  THEN
!981103 Special case where have a value which falls at the limit of the table.
#if (defined D_PRECISION)
	  lrid0 = NINT(SIGN(DFLOAT(mt),rid))
#else
	  lrid0 = NINT(SIGN(FLOAT(mt),rid))
#endif
	  lrid1 = lrid0
	  GO TO 252
!*****C
	ELSE
#if (defined D_PRECISION)
	  tabindrid = SIGN(DFLOAT(mt0) + ((LOG(ABS(rid)) - LOG(ridb(mt0)))/LOG(rri)), rid)
#else
	  tabindrid = SIGN(FLOAT(mt0) + ((LOG(ABS(rid)) - LOG(ridb(mt0)))/LOG(rri)), rid)
#endif
#if (defined D_PRECISION)
	  lrid1 = INT(tabindrid)+NINT(SIGN(DFLOAT(1),rid))
#else
	  lrid1 = INT(tabindrid)+NINT(SIGN(FLOAT(1),rid))
#endif
!yXI
	END IF
!011108yXI It is conceivable that rounding errors may in borderline cases 
!	   throw the calculated table indices for Ri_d off by one.
!	   Check and allow moving one to either side to take care of this.
 	IF((ABS(ridb(lrid1))).LT.(ABS(rid))) THEN
	  lrid1 = lrid1 + SIGN(1,lrid1)
	ELSE IF((ABS(ridb(lrid1-SIGN(1,lrid1)))).GT.(ABS(rid))) THEN
	  lrid1 = lrid1 - SIGN(1,lrid1)
	END IF
!yXI
  250	CONTINUE
!981019-27 Make lrid0 one less or greater than lrid1 according to sgn(rid).
#if (defined D_PRECISION)
        lrid0 = lrid1 - NINT(SIGN(DFLOAT(1),rid)) 
#else
        lrid0 = lrid1 - NINT(SIGN(FLOAT(1),rid)) 
#endif
!*****C
        IF(rid.EQ.0.D0) lrid1 = 1
  252   CONTINUE
!981019-27 Check that the Ri_d value falls within the interpolation interval.
	IF((rid.GT.0.D0.AND.  &
      (rid.LT.ridb(lrid0).OR.rid.GT.ridb(lrid1))).OR.  &
	   (rid.LT.0.D0.AND.  &
      (rid.GT.ridb(lrid0).OR.rid.LT.ridb(lrid1)))) THEN
       if (mytid.eq.0) then
	   WRITE(*,*) "Ri_d is outside interpolation range in interp2d."
	   WRITE(*,*) "rid=  ",rid,"lrid0= ",lrid0,"lrid1= ",lrid1
	   WRITE(*,*) "rid_1(lrid0)=  ",ridb(lrid0), &
                "   rid_1(lrid1)= ",ridb(lrid1)
	   WRITE(*,*) "Program is stopping."
        endif
	   STOP
	END IF
!*****C
!C981022	Artificially reduce Ri if it threatens to surpass Ri_max(Ri_d).
!C	This is to conform to the 1D table's realizability limit treatment. 
!	IF(ri.GT.MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))) THEN
!	  ri = MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))
!	END IF
!C*****C
!981110 Set turbulence to zero if Ri threatens to surpass the realizability limit.
        IF(ri.GT.MIN(rib(irimax(lrid0)),rib(irimax(lrid1)))) THEN
          slq2=0.D0
          sm = 0.D0
          sh = 0.D0
          ss = 0.D0
          RETURN
        END IF
!*****C
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
	IF(ABS(ri).LT.rib(mt0)) THEN
!981019-27 Find Interpolation points in the equally spaced Ri part of the table.
#if (defined D_PRECISION)
	  lri1 = INT(ri/dri)+NINT(SIGN(DFLOAT(1),ri)) 
#else
	  lri1 = INT(ri/dri)+NINT(SIGN(FLOAT(1),ri)) 
#endif
!011107-08yXI Find Interpolation points in exponential absolute value spaced Ri part of the table.
	ELSE IF((ABS(ri)).GE.(rib(mt))) &
        THEN
!981103 Special case where have a value which falls at the limit of the table.
#if (defined D_PRECISION)
	  lri0 = NINT(SIGN(DFLOAT(mt),ri))
#else
	  lri0 = NINT(SIGN(FLOAT(mt),ri))
#endif
	  lri1 = lri0
	  GO TO 272
!*****C
	ELSE
#if (defined D_PRECISION)
	  tabindri = SIGN( DFLOAT(mt0) + ((LOG(ABS(ri)) - LOG(rib(mt0)))/LOG(rri)), ri)
#else
	  tabindri = SIGN( FLOAT(mt0) + ((LOG(ABS(ri)) - LOG(rib(mt0)))/LOG(rri)), ri)
#endif
#if (defined D_PRECISION)
	  lri1 = INT(tabindri)+NINT(SIGN(DFLOAT(1),ri))
#else
	  lri1 = INT(tabindri)+NINT(SIGN(FLOAT(1),ri))
#endif
!yXI
  270	CONTINUE
	END IF
!011108yXI It is conceivable that rounding errors will in borderline cases 
!	   throw the calculated table indices for Ri off by one.
!	   Check and allow moving one to either side to take care of this.
 	IF((ABS(rib(lri1))).LT.(ABS(ri))) THEN
	  lri1 = lri1 + SIGN(1,lri1)
	ELSE IF((ABS(rib(lri1-SIGN(1,lri1)))).GT.(ABS(ri))) THEN
	  lri1 = lri1 - SIGN(1,lri1)
	END IF
!yXI
!981019-27 Make lri0 one less or greater than lri1 according to sgn(ri).
#if (defined D_PRECISION)
        lri0 = lri1 - NINT(SIGN(DFLOAT(1),ri)) 
#else
        lri0 = lri1 - NINT(SIGN(FLOAT(1),ri)) 
#endif
!*****C
        IF(ri.EQ.0.D0) lri1 = 1
  272	CONTINUE
!981019-27 Check that the Ri_d value falls within the interpolation interval.
	IF((ri.GT.0.D0.AND.(ri.LT.rib(lri0).OR.ri.GT.rib(lri1))).OR.  &
	   (ri.LT.0.D0.AND.(ri.GT.rib(lri0).OR.ri.LT.rib(lri1)))) THEN
       if (mytid.eq.0) then
	   WRITE(*,*) "Ri is outside interpolation range in interp2d."
	   WRITE(*,*) "ri=  ",ri,"lri0= ",lri0,"lri1= ",lri1
	   WRITE(*,*) "ri_1(lri0)=  ",rib(lri0), &
                "   ri_1(lri1)= ",rib(lri1)
	   WRITE(*,*) "Program is stopping."
       endif
	   STOP
	END IF
!*****C
!981019-27 Interpolate turbulence fields.
!981027-990112 Introduce table spacing variables.
	deltaridta = ridb(lrid1) - ridb(lrid0)
	deltarita  = rib(lri1)  - rib(lri0)
	deltarid = rid - ridb(lrid0)
	deltari  = ri - rib(lri0)
!	slq2
!981103 Set delta field to zero in special cases falling at limit of the table. 
	IF(lrid1.EQ.lrid0) THEN
	  dslq2_rid = 0.D0
	ELSE
	  dslq2_rid = (slq2b(lri0,lrid1) - slq2b(lri0,lrid0))/ &
                deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dslq2_ri  = 0.D0
	ELSE
	  dslq2_ri = (slq2b(lri1,lrid0) - slq2b(lri0,lrid0))/ &
               deltarita
	END IF
	slq2	 = slq2b(lri0,lrid0) + &
            dslq2_ri*deltari + dslq2_rid*deltarid
!	sm
	IF(lrid1.EQ.lrid0) THEN
	  dsm_rid   = 0.D0
	ELSE
	  dsm_rid = (smb(lri0,lrid1) - smb(lri0,lrid0))/ &
              deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dsm_ri    = 0.D0
	ELSE
	  dsm_ri = (smb(lri1,lrid0) - smb(lri0,lrid0))/ &
             deltarita
	END IF
	sm     = smb(lri0,lrid0) + &
            dsm_ri*deltari + dsm_rid*deltarid
!	sh
	IF(lrid1.EQ.lrid0) THEN
	  dsh_rid   = 0.D0
	ELSE
	  dsh_rid = (shb(lri0,lrid1) - shb(lri0,lrid0))/ &
              deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dsh_ri    = 0.D0
	ELSE
	  dsh_ri = (shb(lri1,lrid0) - shb(lri0,lrid0))/ &
              deltarita
	END IF
	sh     = shb(lri0,lrid0) + &
            dsh_ri*deltari + dsh_rid*deltarid
!	ss
	IF(lrid1.EQ.lrid0) THEN
	  dss_rid   = 0.D0
	ELSE
	  dss_rid = (ssb(lri0,lrid1) - ssb(lri0,lrid0))/ &
              deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dss_ri    = 0.D0
	ELSE
	  dss_ri = (ssb(lri1,lrid0) - ssb(lri0,lrid0))/ &
             deltarita
	END IF
!*****C
	ss     = ssb(lri0,lrid0) + &
            dss_ri*deltari + dss_rid*deltarid
!*****C


	RETURN
	END

