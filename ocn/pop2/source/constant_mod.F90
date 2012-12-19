!=======================================================================
! CVS $Id: constant_mod.F90,v 1.2 2001/12/11 22:52:37 kauff Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/ocn/docn5/constant_mod.F90,v $
! CVS $Name: ccsm2_0_beta58 $
!=======================================================================

module constant_mod

  use shr_kind_mod
  use shr_const_mod

  implicit none

  !----- constants -----
  real(SHR_KIND_R8), parameter ::  latvap   = SHR_CONST_LATVAP ! latent heat of evap   ~ J/kg
  real(SHR_KIND_R8), parameter ::  Tzro     = SHR_CONST_TKFRZ  ! 0 degrees C                       ~ kelvin
  real(SHR_KIND_R8), parameter ::  Tfrz     = Tzro   - 1.8     ! temp of saltwater freezing ~ kelvin
  real(SHR_KIND_R8), parameter ::  pi       = SHR_CONST_PI     ! a famous math constant
  real(SHR_KIND_R8), parameter ::  omega    = SHR_CONST_OMEGA  ! earth's rotation  ~ rad/sec
  real(SHR_KIND_R8), parameter ::  g        = SHR_CONST_G      ! gravity ~ m/s^2
  real(SHR_KIND_R8), parameter ::  DEGtoRAD = PI/180.0         ! PI/180

end module constant_mod

