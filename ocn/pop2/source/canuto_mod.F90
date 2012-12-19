module canuto_mod
#include <def-undef.h>
use precision_mod

!     parameter (nmodel=1,ntbl=2101)
      parameter (nmodel=1,ntbl= 501)
!981022 Minimum of linear range in lookup table.
!     PARAMETER(ri0=-20.D0)
      PARAMETER(ri0=- 4.D0)
!*
!020219D Use "D0" to try to ensure that assigned value is double precision.
!990226 e value
      PARAMETER(e=2.71828182845904509D0)
!981215 \pi value
      PARAMETER(pi=3.14159265358979312D0)
!*****CD
!991107 Temperature=Salt diffusivity model background model swith. 
      PARAMETER(ifback=5)
!*****C
!981006 Salinity model switch
      PARAMETER (ifsali=1)
!*
!040217Zi1b Minimum shear^2 and switch for its use in ifzeroshear=.FALSE. case.
      LOGICAL ifshearmin
      PARAMETER(ifshearmin=.TRUE.)
!lhl0607     PARAMETER(ifshearmin=.FALSE.)
      PARAMETER(s2min=1.D-14)
!******Zi1b
!030424Z Zero shear parameterization for stongly unstable case switch.
      LOGICAL ifzeroshear
      PARAMETER (ifzeroshear=.TRUE.)
!lhl0607      PARAMETER (ifzeroshear=.FALSE.)
!******Z
!	or (+1) to use a Deardorff lengthscale but *NOT* to use a Deardorff-modified
!	{\tau N}, which is consistent with P=\epsilon although not with Deardorff's paper.
!	[See NBp.030714-2,3 in Volume XIX.]
      INTEGER icondear
      PARAMETER(icondear=0)
!*****CZ1a
!000215 Background (epsilon/N^2) dimensionalization of diffusivities switch.
      PARAMETER(ifepson2=2)
!	Value of (epsilon/N^2)/(1 cm/sec^2) used. See NBp.000203-2, Vol.VIII .
!020213-19D Change the parameter's name from "epson2" to "epson2_" 
!	  because epson2 can vary with height above the bottom when ifbotenhance>0.
!030429Z1 Further change the parameter's name from "epson2_" to "epson2__"
!Z1	  because epson2_ can vary with latitude and stratification when ifdeeplat>0.
      REAL(r8) epson2__
      PARAMETER(epson2__=.288D0)
!	  Introduce parameters for bottom enhancement.
!	  ifbotenhance = 0 : no bottom enhancemant
!	               = 1 : epsilon exponentially decreasing with height above bottom to min.
!lhl0608      PARAMETER(ifbotenhance=1) 
      PARAMETER(ifbotenhance=0) 
!14-19 The value of epsilon at the bottom. St. Laurent et al. JPO2001 give 
! 	for slopes epsilon=(3--9)E-9 W/kg and decay scale = (150 +or- 50) meters
!	 and for crests and canyons (2--5)E-9 W/kg and decay scale = (500 +or- 100) meters .
!       1W/kg = 10^7 erg /(kg s) = 10^7 erg /(10^3 g s) = 10^{7-3} erg/(g s)=10^{4} cm^2/(s^3) .
      REAL(r8) eps_bot0
      PARAMETER(eps_bot0=2.D-5)
!14   The bottom enhancement epsilon value at the given height.[See NBp020214-1(Vol.XIII).] 
!19   The scale of decrease of bottom mixing with heigh in centimeters.
      REAL(r8) scale_bot
      PARAMETER(scale_bot=5.D4)
!******CD
!030429Z1 Parameters for a Coriolis-based latitude dependence of \epsilon/N^2.
!Z1   ifdeeplat = 0 : no Coriolis parameter dependence      	
!Z1             = 1 : depends on Coriolis parameter and stratification as in Gregg et al..
      PARAMETER(ifdeeplat=1)
!lhl0609      PARAMETER(ifdeeplat=1)
!Z1   Gregg et al. admit their formula (equation(2)) breaks down right at the equator where
!Z1   it would predict zero epsilon. Figure 1 of Gregg et al. suggests to me that the 
!Z1   value at the equator of (\epsilon/\epsilon_reference) is between 0.02 and 0.05 .
!Z1   In the v_3x3 NCOM the nearest tracer point to the equator is at 0.91^o and at 
!Z1   N_0=5.24e-3 sec^{-1} at this latitude their L(\theta,f) equals \approx 0.0538.
!Z1   Introduce a minimum on L(\theta,f_Coriolis), called eplatidepend in the program.
!Zi1a' 031217 Current guess is best value for eplatidependmin between 0.1 and 0.05 .
      PARAMETER(eplatidependmin=7.D-2)
!*****CZ1
!020221[Mummy's 81st Birthday]D Special check switch to output eps_bot at the bottom (should=eps_bot0).
      PARAMETER(ifcheckbottomeps=0)
      REAL(r8) eps_bot__under	! eps_bot at the level beneath (which is the bottom level for k=n).
!*****C
!000316 Switch for limiting BackGround ra_r to at most Foreground ra_r when Ri>0
!	for R_r in the [R_r_crit_DoubleDiffusion,R_r_crit_SaltFingers] regime.
	PARAMETER(ifrafgmax=1)
!	Need timescale ratios to calculate R_r_crit .
!980610-030403 Common block with ratios of timescales [See NBp.030403-8to10.]
!981125 Salinity background modification switch.
      PARAMETER(ifsalback=5)	
!*****C
!981015 Extend table to avoid use of analytic limiting behavior for sal model.
!000312 Introduce option to have the nonlinear part of the table increase
!	exponentially with the absolute value of the table index.
      PARAMETER(nextrtbl0= 62)
!	ifexpabstable = {0,1} for {don't,do} use ~e^|i| nonlinear table spacing.
      PARAMETER(ifexpabstable=1)
!011107CyXI ******Introduce option ifast=1, yielding ifastexpabs=1, to allow use of alternate interpolation scheme******
!           ******"interp2d_expabs" tailored to exponential absolute value case, should be faster than "interp2d".******
      PARAMETER(ifast=1)
      PARAMETER(ifastexpabs=ifast*ifexpabstable)
!yXI
!     Increase table size for e^|i| spacing to achieve sufficiently high values.
      PARAMETER(nextrtbl1=1000)
      PARAMETER(nextrtbl=nextrtbl0+ifexpabstable*nextrtbl1)
!*****C
!
!981019 mt0 sets the bounds for the part of the table with constant stepsize.
      PARAMETER(nposapprox=101)
      PARAMETER(mt0=ntbl-nposapprox)
      PARAMETER(mt=mt0+nextrtbl)
!*****C
!030401Y Take variables to be used in unstable case with zero shear approximation 
!Y       from my code for HYCOM in common_blocks_giss_fixed2 [See NBp.030401-2to3.].
!030324-28AH Add variables for use in Ri => -infinity case.
!030324-27 Add shearless table as a function of the one variable (N_d^2/N^2) .
       REAL(r8) and2on2a1(-mt:mt),amtaun2a1(-mt:mt),dand2on2 &
        ,sma1(-mt:mt),sha1(-mt:mt),ssa1(-mt:mt) ,rri,rnd2on2,dri, deltheta_r,b1 &
        ,visc_cbu_limit,diff_cbt_limit, theta_rcrp, theta_rcrn,ako,back_l_0 


!******Y
      PARAMETER (ifchengcon=0)
      PARAMETER (idefmld=0)
      PARAMETER (deltemld=0.5D0,delrhmld=0.125D-3)
!980527 For safety's sake make internal arrays' dimension at least 1 . 
!981029 Introduce height array for Ri_d.
!020220D Add z value for lowest ocean level to accomodate bottom enhancement.
!030404Y Introduce array for an2.
!030424-25Z **Introduce an integer parameter for the effect of rotation on the lengthscale.**
!25Z	    Introduce the complex variable zlomega for \sqrt{-B*/f^3} for diagnosis.
      INTEGER ilomega
      PARAMETER (ilomega=0)
!25Z	    amldminlom is the minimum MLD for use of the lengthscale \sqrt{-B*/f^3}.
      PARAMETER(amldminlom=5.D4)	! 50000 centimeters  =  500 meters
!******Z
!000323 Introduce array for a deep lengthscale used in ifepson2=2 case.
!	Do not want to make aldeep a parameter, but want it to be 
!       a big enough array to accomodate the model's number of levels.
      PARAMETER (nbig=100)
!*****C
!981015 Arrays for 2D table used in temperature-salinity model.
      DIMENSION rib(-mt:mt),ridb(-mt:mt), &
     slq2b(-mt:mt,-mt:mt), &
     smb(-mt:mt,-mt:mt),shb(-mt:mt,-mt:mt),ssb(-mt:mt,-mt:mt)
!981017 Table realizability limit array.
      DIMENSION	irimax(-mt:mt)
!*****C
!981215 K_S/K_H (= "sisamax") as a function of angle in Ri_C,Ri_T space 
!	at a radius just before the realizability limit.
      PARAMETER(mt_ra_r=nposapprox-1)
!981218 sisamax ranges from theta_r of -pi/4 to pi/4 to more than 
!	cover the unrealizable region.
      PARAMETER(n_theta_r_oct=INT(((pi/4.D0)*mt_ra_r)/15.D0)*15)
      DIMENSION	sisamax(-n_theta_r_oct:3*n_theta_r_oct)
!*****C
!990226-0301 Critical ra_r[({Ri_T}^2 + {Ri_C}^2)^(1/2)] array for ifsalback>4 .
      DIMENSION ra_rmax(-n_theta_r_oct:3*n_theta_r_oct)
!*****C
!000319 c_y(ra_rmax) array for ifsalback>4 to use as input guess for background.
      DIMENSION c_y_r0(-n_theta_r_oct:3*n_theta_r_oct)
!*****C
!991107	Background Ri constant for ifback >= 4
!	back_ri1
!	S_M,S_H constants -
!	sm_1,sh_1
!*****C
!990301 Background ra_r as a function of \theta_r array for ifsalback>4
      DIMENSION back_ra_r(-n_theta_r_oct:3*n_theta_r_oct)
!990315 Background S_M,S_H,S_S at ra_r as functions of \theta_r arrays for ifsalback>4
      DIMENSION sm_r1(-n_theta_r_oct:3*n_theta_r_oct)
      DIMENSION sh_r1(-n_theta_r_oct:3*n_theta_r_oct)
      DIMENSION ss_r1(-n_theta_r_oct:3*n_theta_r_oct)
!000315 Add (Sl/q)^2 as function of \theta_r array.
      DIMENSION slq2_r1(-n_theta_r_oct:3*n_theta_r_oct)
!000318 Switch to write out polar 2D turbulence table . 
!lhl0711      PARAMETER(ifpolartablewrite=1)
      PARAMETER(ifpolartablewrite=0)
!       Introduce flag for use of \theta_r arrays to interpolate background.
      PARAMETER(ifbg_theta_interp=1)
!*****C
!*****C
!990201-03 Parameters for ifsalback=3 case.
!     Gargett et. al. JPO Vol.11 p.1258-71 gives for "the deep record",
!	\phi_0 = 6 \times 10^{-5} s^{-2} cpm^{-1} . "cpm" is 'cycles per meter'.
!	\phi_0 = 6 \times 10^{-5} s^{-2} (2 pi/100)^{-1} cm
      PARAMETER(back_ph_0=(6.D-5)*(1.D2/(2.D0*pi)))	
!990202 Gargett et. al. favor the value, k_0 = 0.1 cpm. 
!	But k_0=0.05-0.2 cpm might be viable, see section 5 of their paper.
!	Take k_0 = 0.1 cpm * adjust_gargett, where adjust_gargett is adjustable.
! 	Convert to radians per cm: k_0 = 0.1 (2pi/100cm) * adjust_gargett.
!990202-04 used for ifsalback=4 case also, but set adjust_gargett=1 for ifsalback=4 
      PARAMETER(adjust_gargett=1.D0)
      PARAMETER(back_k_0=(0.1D0)*(2.D0)*pi*(1.D-2)*adjust_gargett) 
!       Introduce the lengthscale \Delta_0 \equiv pi/k_0 .
!	The units of \Delta_0 are centimeters, with k_0 in radians per cm.
      PARAMETER(back_del_0=pi/back_k_0)
!*****C
      PARAMETER(back_s2=back_ph_0*back_k_0)
      PARAMETER(back_sm2=1.D0/back_s2)
!990203-04
!	Residual constant background diffusivities to be added to model ones.	
      PARAMETER(v_back0=0.01)	
      PARAMETER(t_back0=0.01)	
      PARAMETER(s_back0=0.01)	
!*****C
!
!990201-03 Parameter for ifsalback=4 case.
      PARAMETER(ri_internal=1.D0)
!*****C
!990226-1107 Parameter for ifback or ifsalback=5 case.
      PARAMETER(backfrac = 85.D-2)
!*****C
!990303 Parameter for ifsalback=6 case.
      PARAMETER(backfact = e**(-1))
!*****C
      dimension ria(ntbl),slq2a(ntbl),sma(ntbl),sha(ntbl)

!
      integer :: ifirst
!020219D REFERENCE NOTE: START OF EXECUTABLE PROGRAM.
!       
!
!YU Jan. 20th, 2011
      end 
