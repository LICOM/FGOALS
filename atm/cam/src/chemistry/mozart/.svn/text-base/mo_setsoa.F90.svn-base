module mo_setsoa

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog

  private
  public :: soa_inti, setsoa, has_soa

  save

  integer, target  :: spc_ndx(9)
  integer, pointer :: soa_ndx, oc1_ndx, oc2_ndx
  integer, pointer :: c10h16_ndx, o3_ndx, oh_ndx
  integer, pointer :: no3_ndx, bigalk_ndx, toluene_ndx
  integer :: rxn_soa(6),react_ndx(6,2)
  real(r8), dimension(6,2) :: alpha &                  ! mass-based stoichiometric coefficients
                             ,k_om                     ! equilibrium gas-particule partition
  real(r8), dimension(6)   :: bulk_yield               ! total yield of condensable compound (ug/m3/ppm)
  real(r8), dimension(6)   :: fraction                 ! fraction of VOC used in reaction

  real(r8), parameter :: avogadro = 6.023e23_r8        ! Avogadro number
  logical :: has_soa = .true.

contains

subroutine soa_inti

  use mo_chem_utls, only : get_spc_ndx, get_rxt_ndx
  use cam_history,  only : addfld, phys_decomp
  use ppgrid,       only : pver
  use spmd_utils,   only : masterproc

  implicit none

  soa_ndx     => spc_ndx(1)
  oc1_ndx     => spc_ndx(2)
  oc2_ndx     => spc_ndx(3)
  c10h16_ndx  => spc_ndx(4)
  o3_ndx      => spc_ndx(5)
  oh_ndx      => spc_ndx(6)
  no3_ndx     => spc_ndx(7)
  bigalk_ndx  => spc_ndx(8)
  toluene_ndx => spc_ndx(9)
!-----------------------------------------------------------------------      
! 	... set species index
!-----------------------------------------------------------------------      
  soa_ndx      = get_spc_ndx( 'SOA' )
  oc1_ndx      = get_spc_ndx( 'OC1' )
  oc2_ndx      = get_spc_ndx( 'OC2' )
  c10h16_ndx   = get_spc_ndx( 'C10H16')
  o3_ndx       = get_spc_ndx( 'OX' )
  if( o3_ndx < 1 ) then
     o3_ndx =  get_spc_ndx( 'O3' )
  end if
  oh_ndx       = get_spc_ndx( 'OH' )
  no3_ndx      = get_spc_ndx( 'NO3' )
  bigalk_ndx   = get_spc_ndx( 'BIGALK' )
  toluene_ndx  = get_spc_ndx( 'TOLUENE' )
  
  has_soa = all( spc_ndx(1:7) > 0 )
!-----------------------------------------------------------------------      
! 	... check if this is an aerosol simulation
!-----------------------------------------------------------------------      
  if( .not. has_soa ) then 
    return
  end if

!-----------------------------------------------------------------------      
! 	... set reaction indicies
!-----------------------------------------------------------------------      
  rxn_soa(1) = get_rxt_ndx( 'soa1' )
  rxn_soa(2) = get_rxt_ndx( 'soa2' )
  rxn_soa(3) = get_rxt_ndx( 'soa3' )
  rxn_soa(4) = get_rxt_ndx( 'soa4' )
  rxn_soa(5) = get_rxt_ndx( 'soa4' )
  rxn_soa(6) = get_rxt_ndx( 'soa5' )
  if( all( rxn_soa(:) < 1 ) ) then
     has_soa = .false.
     return
  else
     if (masterproc) then
        write(iulog,*) '-----------------------------------------'
        write(iulog,*) 'mozart will do soa aerosols'
        write(iulog,*) '-----------------------------------------'
     endif
  end if

!
! define reactants
!
  react_ndx(1,1) = c10h16_ndx
  react_ndx(1,2) = o3_ndx
  react_ndx(2,1) = c10h16_ndx
  react_ndx(2,2) = oh_ndx
  react_ndx(3,1) = c10h16_ndx
  react_ndx(3,2) = no3_ndx
  react_ndx(4,1) = toluene_ndx
  react_ndx(4,2) = oh_ndx
  react_ndx(5,1) = toluene_ndx
  react_ndx(5,2) = oh_ndx
  react_ndx(6,1) = bigalk_ndx
  react_ndx(6,2) = oh_ndx

  if ( masterproc ) then
     write(iulog,*)'soa_inti ',c10h16_ndx, o3_ndx, oh_ndx, no3_ndx, bigalk_ndx, toluene_ndx
     write(iulog,*)'soa_inti ',soa_ndx, oc1_ndx, oc2_ndx
     write(iulog,*)'soa_inti ',react_ndx
  endif
!
! define partitioning coefficients for each reaction
! bulk yields are from Seinfeld and Pandis (1998)
!
! c10h16 + o3 (from Chung and Seinfeld, JGR, 107, 2002)
!
  alpha(1,1)    = 0.067_r8
  alpha(1,2)    = 0.354_r8
  k_om (1,1)    = 0.184_r8
  k_om (1,2)    = 0.0043_r8
  fraction(1)   = 1._r8
  bulk_yield(1) = 762._r8
!
! c10h16 + oh (from Chung and Seinfeld, JGR, 107, 2002)
!
  alpha(2,1)    = 0.067_r8
  alpha(2,2)    = 0.354_r8
  k_om (2,1)    = 0.184_r8
  k_om (2,2)    = 0.0043_r8
  fraction(2)   = 1._r8
  bulk_yield(2) = 762._r8
! 
! c10h16 + no3 (from Chung and Seinfeld, JGR, 107, 2002)
!
  alpha(3,1)    = 1.000_r8
  alpha(3,2)    = 0.000_r8
  k_om (3,1)    = 0.0163_r8
  k_om (3,2)    = 0.0000_r8
  fraction(3)   = 1._r8
  bulk_yield(3) = 762._r8
!
! toluene + oh : toluene (from Odum et al., Environ. Sci. Technol., 1892, 1997)
!
  alpha(4,1)    = 0.038_r8
  alpha(4,2)    = 0.167_r8
  k_om (4,1)    = 0.042_r8
  k_om (4,2)    = 0.0014_r8
  fraction(4)   = 0.7_r8
  bulk_yield(4) = 424._r8
!
! toluene + oh : m-xylene (from Cocker et al., Atmos. Environ., 6079, 2001)
!
  alpha(5,1)    = 0.120_r8
  alpha(5,2)    = 0.019_r8
  k_om (5,1)    = 0.060_r8
  k_om (5,2)    = 0.010_r8
  fraction(5)   = 0.2_r8
  bulk_yield(5) = 419._r8
!
! bigalk + oh : only for alkanes >= heptane (assume low-yield aromatics as in Lack et al.)
!               (from Odum et al., Environ. Sci. Technol., 1892, 1997)
!
  alpha(6,1)    = 0.071_r8
  alpha(6,2)    = 0.138_r8
  k_om (6,1)    = 0.053_r8
  k_om (6,2)    = 0.0019_r8
  fraction(6)   = 0.1_r8
  bulk_yield(6) = 200._r8
!
  call addfld( 'SOA_PROD', 'kg/kg/s', pver, 'A', 'production of SOA', phys_decomp )

  return
end subroutine soa_inti
!
subroutine setsoa(dt,reaction_rates,tfld,vmr,xhnm,ncol,lchnk)
!
! secondary organic aerosol for mozart v2.5
!
! based on Lack et al., JGR, 109, D03203, 2004
!
! rewritten by Jean-Francois Lamarque for updated chemical
! mechanism (March 2004)
!
! adapted to CAM (May 2004)
!
  use ppgrid,       only : pcols, pver
  use chem_mods,    only : adv_mass, gas_pcnst, rxntot
  use mo_chem_utls, only : get_spc_ndx, get_rxt_ndx
  use cam_history,  only : outfld
  use abortutils,   only : endrun
!
  implicit none
!
!-----------------------------------------------------------------------
!      ... dummy arguments
!-----------------------------------------------------------------------
  integer, intent(in)  :: ncol                              ! number columns in chunkx  
  integer, intent(in)  :: lchnk                             ! chunk index
  real(r8), intent(in)     :: dt                                ! time step
  real(r8), intent(in)     :: reaction_rates(ncol,pver,rxntot)  ! reaction rates
  real(r8), intent(inout)  :: vmr(ncol,pver,gas_pcnst)          ! xported species ( vmr )
  real(r8), intent(in)     :: tfld(pcols,pver) &                ! temperature (K)
                             ,xhnm(ncol,pver)                   ! total atms density (mol/cm**3)

!-----------------------------------------------------------------------
!      ... local variables
!-----------------------------------------------------------------------
  integer :: i,k,n
  real(r8) :: m_0
  real(r8) :: mw_soa,yield,prod,soa_mass
  real(r8) :: soa_prod(ncol,pver)
!
! find molecular weight of SOA
!
  mw_soa = adv_mass(soa_ndx)
!
  do k=1,pver
    do i=1,ncol
!
! calculate initial mass of organic aerosols from OC1 and OC2 
! and convert to ug/m3
!
      m_0 = (vmr(i,k,oc1_ndx)+vmr(i,k,oc2_ndx)) * xhnm(i,k) * adv_mass(oc1_ndx)/avogadro * 1.e12_r8
!
! switch based on a minimum value of m_0.  The bulk approach is
! used to initiate the process
!
      if ( m_0 <= 0.2_r8 ) then
!
! bulk theory
!
        soa_mass = 0._r8
        do n=1,6
!
          if ( rxn_soa(n) <= 0 ) cycle
!
          yield = bulk_yield(n)
!
! define chemical production from gas-phase chemistry
!
          prod  = reaction_rates(i,k,rxn_soa(n)) * fraction(n) &
                * vmr(i,k,react_ndx(n,1)) * vmr(i,k,react_ndx(n,2)) * dt
!
! convert from mixing ratio to ppm
!
          prod = 1e6_r8 * prod
!
! collect into total SOA mass
!
          soa_mass = soa_mass + yield * prod
!
        end do
!
      else
!
! partitioning theory
!
        soa_mass = 0._r8
        do n=1,6
!
          if ( rxn_soa(n) <= 0 ) cycle
!
! define yield from available m_0
!
          yield = soa_yield(m_0,alpha(n,1:2),k_om(n,1:2))
!
! define chemical production from gas-phase chemistry
!
          prod  = reaction_rates(i,k,rxn_soa(n)) * fraction(n) &
                * vmr(i,k,react_ndx(n,1)) * vmr(i,k,react_ndx(n,2)) * dt
!
! convert from mixing ratio to mass (ug/m3)       
!
          prod = prod * xhnm(i,k) * mw_soa/avogadro * 1.e12_r8
!
! collect into total SOA mass
!
          soa_mass = soa_mass + yield * prod
!
        end do
!
      endif
!
! convert from ug/m3 to mixing ratio and update vmr
!
      vmr(i,k,soa_ndx) = vmr(i,k,soa_ndx) + soa_mass * 1.e-12_r8 * avogadro/(mw_soa*xhnm(i,k))
      if ( vmr(i,k,soa_ndx) > 1.e0_r8 ) then
        write(iulog,*)i,k,soa_mass,m_0
        call endrun('soa_yield: vmr(i,k,soa_ndx) > 1.e0_r8')
      endif
!
      soa_prod(i,k) = soa_mass*1.e-12*avogadro/(28.966_r8*xhnm(i,k)*dt)
    end do
  end do
!
  call outfld('SOA_PROD',soa_prod(:ncol,:),ncol, lchnk)
  return
end subroutine setsoa
!
real(r8) function soa_yield(m_0,alpha,k)
!
  implicit none
!
  real(r8) :: m_0
  real(r8), dimension(2) :: alpha,k
!
!!$  soa_yield = m_0 * ( ((alpha(1)*k(1))/(1._r8+(alpha(1)*k(1)*m_0))) &
!!$                    + ((alpha(2)*k(2))/(1._r8+(alpha(2)*k(2)*m_0))) )
  soa_yield = m_0 * ( ((alpha(1)*k(1))/(1._r8+k(1)*m_0)) &
                  +   ((alpha(2)*k(2))/(1._r8+k(2)*m_0)) ) 
!
  return
end function soa_yield
!
end module mo_setsoa
