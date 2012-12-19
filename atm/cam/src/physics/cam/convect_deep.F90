
module convect_deep
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to several deep convection interfaces. Currently includes:
!    Zhang-McFarlane (default)
!    Kerry Emanuel 
!
!
! Author: D.B. Coleman, Sep 2004
!
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use cam_logfile,  only: iulog

   implicit none

   save
   private                         ! Make default type private to the module

! Public methods

   public ::&
      convect_deep_register,           &! register fields in physics buffer
      convect_deep_init,               &! initialize donner_deep module
      convect_deep_tend,               &! return tendencies
      convect_deep_tend_2,             &! return tendencies
      deep_scheme_does_scav_trans             ! = .t. if scheme does scavenging and conv. transport
   
! Private module data
   character(len=16) :: deep_scheme    ! default set in phys_control.F90, use namelist to change
! Physics buffer indices 
   integer     ::  icwmrdp_idx      = 0 
   integer     ::  rprddp_idx       = 0 
   integer     ::  nevapr_dpcu_idx  = 0 
   integer     ::  cldtop_idx       = 0 
   integer     ::  cldbot_idx       = 0 
   integer     ::  cld_idx          = 0 
   integer     ::  fracis_idx       = 0 

!=========================================================================================
  contains 

!=========================================================================================
function deep_scheme_does_scav_trans()
!
! Function called by tphysbc to determine if it needs to do scavenging and convective transport
! or if those have been done by the deep convection scheme. Each scheme could have its own
! identical query function for a less-knowledgable interface but for now, we know that KE 
! does scavenging & transport, and ZM doesn't
!

  logical deep_scheme_does_scav_trans

  if ( deep_scheme .eq. 'ZM' ) then
     deep_scheme_does_scav_trans = .false.
  else
     deep_scheme_does_scav_trans = .true.
  endif

  return

end function deep_scheme_does_scav_trans

!=========================================================================================
subroutine convect_deep_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  use phys_buffer, only: pbuf_times, pbuf_add
  use zm_conv_intr, only: zm_conv_register
  use phys_control, only: phys_getopts

  implicit none

  integer idx

  ! get deep_scheme setting from phys_control
  call phys_getopts(deep_scheme_out = deep_scheme)

  select case ( deep_scheme )
  case('ZM') !    Zhang-McFarlane (default)
     call zm_conv_register
  end select

  call pbuf_add('ICWMRDP' , 'physpkg', 1,pver,      1,         icwmrdp_idx)
  call pbuf_add('RPRDDP' , 'physpkg', 1,pver,       1,          rprddp_idx)
  call pbuf_add('NEVAPR_DPCU' , 'physpkg', 1,pver,      1, nevapr_dpcu_idx)

end subroutine convect_deep_register

!=========================================================================================


subroutine convect_deep_init(hypi)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use pmgrid,        only: plevp
  use spmd_utils,    only: masterproc
  use zm_conv_intr,  only: zm_conv_init
  use abortutils,    only: endrun
  use phys_buffer,   only: pbuf_get_fld_idx


  implicit none

  real(r8),intent(in) :: hypi(plevp)        ! reference pressures at interfaces

  integer k

  select case ( deep_scheme )
  case('off') !     ==> no deep convection
     if (masterproc) write(iulog,*)'convect_deep: no deep convection selected'
  case('ZM') !    1 ==> Zhang-McFarlane (default)
     if (masterproc) write(iulog,*)'convect_deep initializing Zhang-McFarlane convection'
     call zm_conv_init(hypi)
  case default
     if (masterproc) write(iulog,*)'WARNING: convect_deep: no deep convection scheme. May fail.'
  end select

  cldtop_idx = pbuf_get_fld_idx('CLDTOP')
  cldbot_idx = pbuf_get_fld_idx('CLDBOT')
  cld_idx    = pbuf_get_fld_idx('CLD')
  fracis_idx = pbuf_get_fld_idx('FRACIS')

end subroutine convect_deep_init
!=========================================================================================
!subroutine convect_deep_tend(state, ptend, tdt, pbuf)

subroutine convect_deep_tend(prec    , &
     pblh    ,mcon    ,cme     ,          &
     tpert   ,dlf     ,pflx    ,zdu      , &
     rliq    , &
     ztodt   ,snow    ,&
     state   ,ptend   ,landfrac ,pbuf  )

   use physics_types, only: physics_state, physics_ptend, physics_tend, physics_ptend_init
   use phys_buffer,   only: pbuf_size_max, pbuf_fld 
   use constituents, only: pcnst
   use zm_conv_intr, only: zm_conv_tend
#if ( defined WACCM_PHYS )
   use gw_drag,         only: idx_zmdt
   use cam_history,     only: outfld
   use physconst,       only: cpair
#endif

! Arguments
   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
   type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer

   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblh(pcols)                 ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess 
   real(r8), intent(in) :: landfrac(pcols)                ! Land fraction
      

   real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux

   real(r8), intent(out) :: prec(pcols)   ! total precipitation
   real(r8), intent(out) :: snow(pcols)   ! snow from ZM convection 
   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals

   real(r8), pointer, dimension(:) :: jctop
   real(r8), pointer, dimension(:) :: jcbot
   real(r8), pointer, dimension(:,:,:) :: cld        
   real(r8), pointer, dimension(:,:) :: ql           ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

   real(r8), pointer, dimension(:,:) :: evapcdp      ! Evaporation of deep convective precipitation

  real(r8) zero(pcols, pver)

  integer i, k

#if ( defined WACCM_PHYS )
   real(r8), pointer, dimension(:,:) :: zmdt
   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
#endif

   jctop => pbuf(cldtop_idx)%fld_ptr(1,1:pcols,1,state%lchnk,1)
   jcbot => pbuf(cldbot_idx)%fld_ptr(1,1:pcols,1,state%lchnk,1)

  select case ( deep_scheme )
  case('off') !    0 ==> no deep convection
    zero = 0     
    mcon = 0
    dlf = 0
    pflx = 0
    cme = 0
    zdu = 0
    prec = 0
    snow = 0
    rliq = 0

!
! Associate pointers with physics buffer fields
!

   cld =>  pbuf(cld_idx)%fld_ptr(1,1:pcols,1:pver,state%lchnk,:) 
   ql =>   pbuf(icwmrdp_idx)%fld_ptr(1,1:pcols,1:pver,state%lchnk,1)
   rprd => pbuf(rprddp_idx)%fld_ptr(1,1:pcols,1:pver,state%lchnk,1)
   fracis  => pbuf(fracis_idx)%fld_ptr(1,1:pcols,1:pver,state%lchnk,1:pcnst)
   evapcdp => pbuf(nevapr_dpcu_idx)%fld_ptr(1,1:pcols,1:pver,state%lchnk,1)

    jctop = 0
    jcbot = 0
    cld = 0
    ql = 0
    rprd = 0
    fracis = 0
    evapcdp = 0



   call physics_ptend_init(ptend)
   ptend%name = "convect_deep"
  case('ZM') !    1 ==> Zhang-McFarlane (default)
     call zm_conv_tend(prec     , &
          pblh    ,mcon    ,cme     ,          &
          tpert   ,dlf     ,pflx    ,zdu      , &
          rliq    , &
          ztodt   ,snow    ,&
          jctop, jcbot , &
          state   ,ptend   ,landfrac ,pbuf  )

  end select

#if ( defined WACCM_PHYS )
   zmdt  => pbuf(idx_zmdt) %fld_ptr(1,:,:,state%lchnk,1)
   ftem(:state%ncol,:pver) = ptend%s(:state%ncol,:pver)/cpair
   zmdt(:state%ncol,:pver) = ftem(:state%ncol,:pver)
   call outfld('ZMDT    ',ftem           ,pcols   ,state%lchnk   )
   call outfld('ZMDQ    ',ptend%q(:,:,1) ,pcols   ,state%lchnk   )
#endif


end subroutine convect_deep_tend
!=========================================================================================


subroutine convect_deep_tend_2( state,  ptend,  ztodt, pbuf  )

   use physics_types, only: physics_state, physics_ptend
   use phys_buffer,   only: pbuf_size_max, pbuf_fld
   use constituents, only: pcnst
   use zm_conv_intr, only: zm_conv_tend_2

! Arguments
   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
   type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)

   if ( deep_scheme .eq. 'ZM' ) then  !    1 ==> Zhang-McFarlane (default)
      call zm_conv_tend_2( state,   ptend,  ztodt,  pbuf ) 
   end if


end subroutine convect_deep_tend_2


end module
