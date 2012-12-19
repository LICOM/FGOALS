!=======================================================================
! CVS $Id: buf_mod.F90,v 1.1 2001/05/08 16:25:38 mvertens Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/ocn/docn5/buf_mod.F90,v $
! CVS $Name: ccsm2_0_beta58 $
!=======================================================================

module buf_mod
  use shr_kind_mod
  implicit none

  !----- message send/receive parameters--------
  integer,parameter :: nibuff = 100  ! size of control buffers
  integer,parameter :: nsnd   =  7   ! number of 2d fields to send
  integer,parameter :: nrcv   = 14   ! number of 2d fields to receive

  !----- flags -----
  integer   ibuffs(nibuff) ! control flags sent to coupler
  integer   ibuffr(nibuff) ! control flags received from coupler

  !----- domain information -----
  integer               ::   nxg       ! size of x grid in the global domain
  integer               ::   nyg       ! size of y grid in the global domain
  integer               ::   nx        ! size of x grid
  integer               ::   ny        ! size of y grid
  real(SHR_KIND_R8), allocatable ::   xc(:,:)   ! x-grid center coordinates ~ deg east
  real(SHR_KIND_R8), allocatable ::   yc(:,:)   ! y-grid center coordinates ~ deg north
  ! lihuimin 2012.6.15
  integer               ::   nv        ! corners
  real(SHR_KIND_R8), allocatable ::   xv(:,:,:) ! x-grid vertex coordinates ~ deg east 
  real(SHR_KIND_R8), allocatable ::   yv(:,:,:) ! y-grid vertex coordinates ~ deg north 
  ! end modi
  real(SHR_KIND_R8), allocatable ::   area(:,:) ! grid cell area            ~ rad^2
  integer ,save, allocatable ::   mask(:,:) ! domain mask: 0 <=> cell NOT in domain

  !----- states sent to coupler -----
  real(SHR_KIND_R8), allocatable ::   T_cpl(:,:)   ! temperature               ~ Kelvin
  real(SHR_KIND_R8), allocatable ::   s_cpl(:,:)   ! state: salinity               ~ ppt
  real(SHR_KIND_R8), allocatable ::   u_cpl(:,:)   ! state: velocity, zonal        ~ m/s
  real(SHR_KIND_R8), allocatable ::   v_cpl(:,:)   ! state: velocity, meridional   ~ m/s
#ifdef USE_OCN_CARBON
  real(SHR_KIND_R8), allocatable ::   co2_cpl(:,:)   ! state: air-sea CO2 of ocean  ~ mol
#endif

  !----- fluxes sent to coupler -----
  real(SHR_KIND_R8), allocatable ::   dhdx (:,:)   ! dh/dx (surface slope)     ~ m/m
  real(SHR_KIND_R8), allocatable ::   dhdy (:,:)   ! dh/dy (surface slope)     ~ m/m
  real(SHR_KIND_R8), allocatable ::   Q    (:,:)   ! heat of fusion    (q > 0) ~ W/m^2

  !----- fluxes received from coupler -----
  real(SHR_KIND_R8), allocatable ::   taux (:,:)   ! wind stress, zonal        ~     N/m^2
  real(SHR_KIND_R8), allocatable ::   tauy (:,:)   ! wind stress, meridional   ~     N/m^2
  real(SHR_KIND_R8), allocatable ::   netsw(:,:)   ! heat flux: net short wave ~     W/m^2
  real(SHR_KIND_R8), allocatable ::   lat1(:,:)    ! heat flux: latent         ~     W/m^2
  real(SHR_KIND_R8), allocatable ::   sen  (:,:)   ! heat flux: sensible       ~     W/m^2
  real(SHR_KIND_R8), allocatable ::   lwup (:,:)   ! heat flux: long-wave up   ~     W/m^2
  real(SHR_KIND_R8), allocatable ::   lwdn (:,:)   ! heat flux: long-wave down ~     W/m^2
  real(SHR_KIND_R8), allocatable ::   melth(:,:)   ! heat flux: melt           ~     W/m^2
  real(SHR_KIND_R8), allocatable ::   salt (:,:)   ! salt flux                 ~ ppt/s/m^2
  real(SHR_KIND_R8), allocatable ::   prec (:,:)   ! water flux: precip        ~  kg/s/m^2
  real(SHR_KIND_R8), allocatable ::   evap (:,:)   ! water flux: evaporation   ~  kg/s/m^2
  real(SHR_KIND_R8), allocatable ::   meltw(:,:)   ! water flux: melt          ~  kg/s/m^2
  real(SHR_KIND_R8), allocatable ::   roff (:,:)   ! water flux: runoff        ~  kg/s/m^2
  real(SHR_KIND_R8), allocatable ::   duu10n(:,:)   ! 10m neutral wind speed squared ~m^2/s^2
#ifdef USE_OCN_CARBON
  real(SHR_KIND_R8), allocatable ::   pco2(:,:)   ! state: Partial CO2 Pressure  ~ mol
  real(SHR_KIND_R8), allocatable ::   pco2up(:,:)
#endif
  real(SHR_KIND_R8) pco2ups
  real(SHR_KIND_R8) pco2s

  !----- states received from coupler -----
  real(SHR_KIND_R8), allocatable ::   ifrac(:,:)   ! ice fraction              ~ in [0,1]
  real(SHR_KIND_R8), allocatable ::   patm (:,:)   ! atm sea level pressure    ~ Pa

  !----- send/receive buffers -------------
  real(SHR_KIND_R8), allocatable ::   buffs(:,:) ! contiguous array for sending data
  real(SHR_KIND_R8), allocatable ::   buffr(:,:) ! contiguous array for recving data

end module buf_mod
