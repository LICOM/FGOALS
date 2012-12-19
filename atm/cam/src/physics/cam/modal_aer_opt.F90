module modal_aer_opt

!     parameterizes aerosol coefficients using chebychev polynomial
!     parameterize aerosol radiative properties in terms of
!     surface mode wet radius and wet refractive index

!     Ghan and Zaveri, JGR 2007.

!     uses Wiscombe's (1979) mie scattering code


use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
use ppgrid,            only: pcols, pver, pverp
use constituents,      only: pcnst
use spmd_utils,        only: masterproc
use phys_control,      only: cam_chempkg_is
use physconst,         only: rhoh2o, rga, rair
use radconstants,      only: nswbands, nlwbands, idx_sw_diag
use rad_constituents,  only: rad_cnst_get_aer_idx, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props
use physics_types,     only: physics_state
use phys_buffer,       only: pbuf_fld, pbuf_get_fld_idx
use pio,               only: file_desc_t, var_desc_t, pio_inq_dimlen, pio_inq_dimid, pio_inq_varid, &
                             pio_get_var, pio_nowrite, pio_closefile
use cam_pio_utils,     only: cam_pio_openfile
use cam_history,       only: phys_decomp, addfld, add_default, outfld
use cam_history_support, only: fillvalue
use cam_logfile,       only: iulog
use perf_mod,          only: t_startf, t_stopf
use abortutils,        only: endrun

implicit none
private
save

public :: modal_aer_opt_readnl, modal_aer_opt_init, modal_aero_sw, modal_aero_lw


character(len=*), parameter :: unset_str = 'UNSET'

! Namelist variables:
character(shr_kind_cl)      :: modal_optics_file = unset_str   ! full pathname for modal optics dataset
character(shr_kind_cl)      :: water_refindex_file = unset_str ! full pathname for water refractive index dataset

! description of the modal aerosols
integer :: ntot_amode   ! number of modes
integer :: nspec_max    ! maximum number of species in a mode
integer, pointer :: nspec_amode(:) ! number of species in each mode
integer, parameter :: clen1 = 8   ! length of character strings containing species names
integer, parameter :: clen2 = 10  ! length of character strings containing species types
integer,              pointer :: spec_idx(:,:)       ! indices in the aerosol list for the species in each mode
character(len=clen1), pointer :: xname_massptr(:,:)  ! names of species in each mode
character(len=clen2), pointer :: xname_spectype(:,:) ! types of species in each mode

! coefficients for parameterizing aerosol radiative properties
! in terms of refractive index and wet radius
! Allocate in read_modal_optics after ntot_amode is read.
integer, parameter :: ncoef=5, prefr=7, prefi=10
real(r8), pointer :: extpsw(:,:,:,:,:) ! specific extinction
real(r8), pointer :: abspsw(:,:,:,:,:) ! specific absorption
real(r8), pointer :: asmpsw(:,:,:,:,:) ! asymmetry factor
real(r8), pointer :: absplw(:,:,:,:,:) ! specific absorption
real(r8), pointer :: sigma_logr_aer(:) ! geometric standard deviation of size distribution
real(r8), pointer :: alnsg_amode(:)    ! log(sigma_logr_aer)

! For BFB testing
real(r8) :: sigmag_amode(3) = (/ 1.800, 1.600, 1.800 /)


real(r8) :: xrmin, xrmax

! Table read in read_modal_optics
real(r8) :: refrtabsw(prefr,nswbands) ! table of real refractive indices for aerosols visible
real(r8) :: refitabsw(prefi,nswbands) ! table of imag refractive indices for aerosols visible
real(r8) :: refrtablw(prefr,nlwbands) ! table of real refractive indices for aerosols infrared
real(r8) :: refitablw(prefi,nlwbands) ! table of imag refractive indices for aerosols infrared

! refractive index for water read in read_water_refindex
complex :: crefwsw(nswbands) ! complex refractive index for water visible
complex :: crefwlw(nlwbands) ! complex refractive index for water infrared

! some types for constructing arrays of pointers
type c_ptr1d_t
   complex, pointer :: val(:)
end type c_ptr1d_t
type r_ptr2d_t
   real(r8), pointer :: val(:,:)
end type r_ptr2d_t

!===============================================================================
CONTAINS
!===============================================================================

subroutine modal_aer_opt_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'modal_aer_opt_readnl'

   namelist /modal_aer_opt_nl/ modal_optics_file, water_refindex_file
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'modal_aer_opt_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, modal_aer_opt_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(modal_optics_file,   len(modal_optics_file), mpichar, 0, mpicom)
   call mpibcast(water_refindex_file, len(water_refindex_file), mpichar, 0, mpicom)
#endif


end subroutine modal_aer_opt_readnl

!===============================================================================

subroutine modal_aer_opt_init()

   use ioFileMod,        only: getfil

   ! Local variables

   integer  :: i, m
   real(r8) :: rmmin, rmmax       ! min, max aerosol surface mode radius treated (m)

   character(len=256) :: locfile
   !----------------------------------------------------------------------------

   rmmin = 0.01e-6
   rmmax = 25.e-6
   xrmin = log(rmmin)
   xrmax = log(rmmax)

   call getfil(modal_optics_file, locfile)
   call read_modal_optics(locfile)
   if (masterproc) write(iulog,*) "modal_aer_opt_init: read modal optics file:", trim(locfile)

   allocate(alnsg_amode(ntot_amode))
   do m = 1, ntot_amode
      ! This version is for BFB testing
      !alnsg_amode(m) = log( sigmag_amode(m) )
      ! this version uses data read from the modal optics file
      alnsg_amode(m) = log( sigma_logr_aer(m) )
   end do

   ! lookup the names of the mode species in the aerosol list, and save the indices
   allocate(spec_idx(nspec_max,ntot_amode))
   do m = 1, ntot_amode
      do i = 1, nspec_amode(m)
         spec_idx(i,m) = rad_cnst_get_aer_idx(0, xname_massptr(i,m))
      end do
   end do

   call getfil(water_refindex_file, locfile)
   call read_water_refindex(locfile)
   if (masterproc) write(iulog,*) "modal_aer_opt_init: read water refractive index file:", trim(locfile)


   ! Add diagnostic fields to history output.
   call addfld ('EXTINCT','/m  ',pver,    'A','Aerosol extinction',phys_decomp, flag_xyfill=.true.)
   call addfld ('ABSORB','/m  ',pver,    'A','Aerosol absorption',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODVIS','  ',1,    'A','Aerosol optical depth 550 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODABS','  ',1,    'A','Aerosol absorption optical depth 550 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODMODE1','  ',1,    'A','Aerosol optical depth 550 nm mode 1',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODMODE2','  ',1,    'A','Aerosol optical depth 550 nm mode 2',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODMODE3','  ',1,    'A','Aerosol optical depth 550 nm mode 3',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODDUST1','  ',1,    'A','Aerosol optical depth 550 nm model 1 from dust',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODDUST2','  ',1,    'A','Aerosol optical depth 550 nm model 2 from dust',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODDUST3','  ',1,    'A','Aerosol optical depth 550 nm model 3 from dust',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDEN1','kg/m2',1,    'A','Aerosol burden mode 1',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDEN2','kg/m2',1,    'A','Aerosol burden mode 2',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDEN3','kg/m2',1,    'A','Aerosol burden mode 3',phys_decomp, flag_xyfill=.true.)
   call addfld ('SSAVIS','  ',1,    'A','Aerosol singel-scatter albedo',phys_decomp, flag_xyfill=.true.)

   call add_default ('EXTINCT', 1, ' ')
   call add_default ('ABSORB', 1, ' ')
   call add_default ('AODVIS', 1, ' ')
   call add_default ('AODABS', 1, ' ')
   call add_default ('AODMODE1', 1, ' ')
   call add_default ('AODMODE2', 1, ' ')
   call add_default ('AODMODE3', 1, ' ')
   call add_default ('AODDUST1', 1, ' ')
   call add_default ('AODDUST2', 1, ' ')
   call add_default ('AODDUST3', 1, ' ')
   call add_default ('BURDEN1', 1, ' ')
   call add_default ('BURDEN2', 1, ' ')
   call add_default ('BURDEN3', 1, ' ')
   call add_default ('SSAVIS', 1, ' ')

   if (cam_chempkg_is('trop_mam7')) then
      call addfld ('AODMODE4','  ',1,    'A','Aerosol optical depth 550 nm mode 4',phys_decomp, flag_xyfill=.true.)
      call addfld ('AODMODE5','  ',1,    'A','Aerosol optical depth 550 nm mode 5',phys_decomp, flag_xyfill=.true.)
      call addfld ('AODMODE6','  ',1,    'A','Aerosol optical depth 550 nm mode 6',phys_decomp, flag_xyfill=.true.)
      call addfld ('AODMODE7','  ',1,    'A','Aerosol optical depth 550 nm mode 7',phys_decomp, flag_xyfill=.true.)
      call addfld ('BURDEN4','kg/m2',1,    'A','Aerosol burden mode 4',phys_decomp, flag_xyfill=.true.)
      call addfld ('BURDEN5','kg/m2',1,    'A','Aerosol burden mode 5',phys_decomp, flag_xyfill=.true.)
      call addfld ('BURDEN6','kg/m2',1,    'A','Aerosol burden mode 6',phys_decomp, flag_xyfill=.true.)
      call addfld ('BURDEN7','kg/m2',1,    'A','Aerosol burden mode 7',phys_decomp, flag_xyfill=.true.)
      call add_default ('AODMODE4', 1, ' ')
      call add_default ('AODMODE5', 1, ' ')
      call add_default ('AODMODE6', 1, ' ')
      call add_default ('AODMODE7', 1, ' ')
      call add_default ('BURDEN4', 1, ' ')
      call add_default ('BURDEN5', 1, ' ')
      call add_default ('BURDEN6', 1, ' ')
      call add_default ('BURDEN7', 1, ' ')
   end if

end subroutine modal_aer_opt_init

!===============================================================================

subroutine modal_aero_sw(state, pbuf, nnite, idxnite, &
                         tauxar, wa, ga, fa)

   ! calculates aerosol sw radiative properties

   type(physics_state), intent(in) :: state          ! state variables
   type(pbuf_fld),      intent(in) :: pbuf(:)        ! physics buffer
   integer,             intent(in) :: nnite          ! number of night columns
   integer,             intent(in) :: idxnite(nnite) ! local column indices of night columns

   real(r8), intent(out) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth
   real(r8), intent(out) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo
   real(r8), intent(out) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor
   real(r8), intent(out) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction

   ! Local variables
   integer :: i, ifld, isw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk

   real(r8), pointer :: dgnumwet(:,:,:)     ! number mode diameter
   real(r8), pointer :: qaerwat(:,:,:)      ! aerosol water (g/g)

   real(r8) :: mass(pcols,pver)        ! layer mass
   real(r8) :: air_density(pcols,pver) ! (kg/m3)

   type(r_ptr2d_t), allocatable :: specmmr(:,:) ! species mass mixing ratio
   real(r8),        allocatable :: specdens(:,:) ! species density (kg/m3)
   type(c_ptr1d_t), allocatable :: specrefindex(:,:) ! species refractive index

   real(r8), pointer :: radsurf(:,:,:)    ! aerosol surface mode radius
   real(r8), pointer :: logradsurf(:,:,:) ! log(aerosol surface mode radius)
   real(r8), pointer :: cheb(:,:,:,:)

   complex  :: crefin(pcols)   ! complex refractive index
   real(r8) :: refr(pcols)     ! real part of refractive index
   real(r8) :: refi(pcols)     ! imaginary part of refractive index

   real(r8) :: vol(pcols)      ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)   ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: watervol(pcols) ! volume concentration of water in each mode (m3/kg)
   real(r8) :: wetvol(pcols)   ! volume concentration of wet mode (m3/kg)

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cext(pcols,ncoef), cabs(pcols,ncoef), casm(pcols,ncoef)
   real(r8) :: pext(pcols)     ! parameterized specific extinction (m2/kg)
   real(r8) :: specpext(pcols) ! specific extinction (m2/kg)
   real(r8) :: dopaer(pcols)   ! aerosol optical depth in layer
   real(r8) :: pabs(pcols)     ! parameterized specific absorption (m2/kg)
   real(r8) :: pasm(pcols)     ! parameterized asymmetry factor
   real(r8) :: palb(pcols)     ! parameterized single scattering albedo

   ! Diagnostics output for visible band only
   real(r8) :: extinct(pcols,pver)
   real(r8) :: absorb(pcols,pver)
   real(r8) :: aodvis(pcols)               ! extinction optical depth
   real(r8) :: aodabs(pcols)               ! absorption optical depth
   real(r8) :: ssavis(pcols)
   real(r8) :: dustvol(pcols)              ! volume concentration of dust in aerosol mode (m3/kg)
   real(r8), pointer :: aodmode(:,:)
   real(r8), pointer :: dustaodmode(:,:)   ! dust aod in aerosol mode
   real(r8), pointer :: burden(:,:)
   real(r8), pointer :: colext(:,:)

   logical :: savaervis ! true if visible wavelength (0.55 micron)

   ! debug output
   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf            ! volume fraction of insoluble aerosol
   character(len=*), parameter :: subname = 'modal_aero_sw'
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! allocate local storage
   allocate( &
      radsurf(pcols,pver,ntot_amode),    &
      logradsurf(pcols,pver,ntot_amode), &
      cheb(ncoef,ntot_amode,pcols,pver), &
      aodmode(pcols,ntot_amode),         &
      dustaodmode(pcols,ntot_amode),     &
      burden(pcols,ntot_amode),          &
      colext(pcols,ntot_amode))

   ! pointers to physics buffer
   ifld = pbuf_get_fld_idx('DGNUMWET')
   dgnumwet => pbuf(ifld)%fld_ptr(1,:,:,lchnk,:)

   ifld = pbuf_get_fld_idx('QAERWAT')
   if (associated(pbuf(ifld)%fld_ptr)) then
      qaerwat => pbuf(ifld)%fld_ptr(1,:,:,lchnk,:)
   else
      call endrun(subname//': pbuf for QAERWAT not allocated')
   end if

   ! calc size parameter for all columns
   call modal_size_parameters(ncol, dgnumwet, radsurf, logradsurf, cheb)

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8
   wa(:ncol,:,:)     = 0._r8
   ga(:ncol,:,:)     = 0._r8
   fa(:ncol,:,:)     = 0._r8

   ! zero'th layer does not contain aerosol
   tauxar(1:ncol,0,:)  = 0._r8
   wa(1:ncol,0,:)      = 0.925_r8
   ga(1:ncol,0,:)      = 0.850_r8
   fa(1:ncol,0,:)      = 0.7225_r8

   mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
   air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

   ! diagnostics for visible band
   extinct(1:ncol,:)     = 0.0_r8
   absorb(1:ncol,:)      = 0.0_r8
   aodvis(1:ncol)        = 0.0_r8
   aodabs(1:ncol)        = 0.0_r8
   ssavis(1:ncol)        = 0.0_r8
   aodmode(1:ncol,:)     = 0.0_r8
   dustaodmode(1:ncol,:) = 0.0_r8
   burden(:ncol,:ntot_amode) = 0._r8

   ! access the mixing ratio and properties of the modal species
   allocate( &
      specmmr(nspec_max,ntot_amode),     &
      specdens(nspec_max,ntot_amode),    &
      specrefindex(nspec_max,ntot_amode) )

   do m = 1, ntot_amode
      do l = 1, nspec_amode(m)
         call rad_cnst_get_aer_mmr(0, spec_idx(l,m),  state, pbuf, specmmr(l,m)%val)
         call rad_cnst_get_aer_props(0, spec_idx(l,m), density_aer=specdens(l,m), &
                                     refindex_aer_sw=specrefindex(l,m)%val)
      end do
   end do

   ! loop over all aerosol modes
   do m = 1, ntot_amode

      do isw = 1, nswbands
         savaervis = isw .eq. idx_sw_diag


         do k = 1, pver

            ! form bulk refractive index
            crefin(:ncol) = 0._r8
            dryvol(:ncol) = 0._r8
            dustvol(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec_amode(m)
               do i = 1, ncol
                  vol(i)      = specmmr(l,m)%val(i,k)/specdens(l,m)
                  dryvol(i)   = dryvol(i) + vol(i)
                  crefin(i)   = crefin(i) + vol(i)*specrefindex(l,m)%val(isw)
               end do

               ! compute some diagnostics for visible band only
               if (savaervis) then

                  do i = 1, ncol
                     burden(i,m) = burden(i,m) + specmmr(l,m)%val(i,k)*mass(i,k)
                  end do

                  if (trim(xname_spectype(l,m)) == 'dust') then
                     do i = 1, ncol
                        dustvol(i) = specmmr(l,m)%val(i,k)/specdens(l,m)
                     end do
                  end if

               end if

            end do ! species loop


            do i = 1, ncol
               watervol(i) = qaerwat(i,k,m)/rhoh2o
               wetvol(i) = watervol(i) + dryvol(i)
               if (watervol(i) < 0._r8) then
                  if (abs(watervol(i)) .gt. 1.e-1*wetvol(i)) then
                     write(iulog,'(a,4e10.2,a)') 'watervol,wetvol=', &
                        watervol(i), wetvol(i), ' in '//subname
                     !  call endrun()
                  end if
                  watervol(i) = 0._r8
                  wetvol(i) = dryvol(i)
               end if

               ! volume mixing
               crefin(i) = crefin(i) + watervol(i)*crefwsw(isw)
               crefin(i) = crefin(i)/max(wetvol(i),1.e-60_r8)
               refr(i)   = real(crefin(i))
               refi(i)   = abs(aimag(crefin(i)))
            end do

            ! call t_startf('binterp')

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(extpsw(:,:,:,m,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, cext)
            call binterp(abspsw(:,:,:,m,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, cabs)
            call binterp(asmpsw(:,:,:,m,isw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                         itab, jtab, ttab, utab, casm)

            ! call t_stopf('binterp')

            ! parameterized optical properties
            do i=1,ncol

               if (logradsurf(i,k,m) .le. xrmax) then
                  pext(i) = 0.5_r8*cext(i,1)
                  do nc = 2, ncoef
                     pext(i) = pext(i) + cheb(nc,m,i,k)*cext(i,nc)
                  enddo
                  pext(i) = exp(pext(i))
               else
                  pext(i) = 1.5_r8/(radsurf(i,k,m)*rhoh2o) ! geometric optics
               endif

               ! convert from m2/kg water to m2/kg aerosol
               specpext(i) = pext(i)
               pext(i) = pext(i)*wetvol(i)*rhoh2o
               pabs(i) = 0.5_r8*cabs(i,1)
               pasm(i) = 0.5_r8*casm(i,1)
               do nc = 2, ncoef
                  pabs(i) = pabs(i) + cheb(nc,m,i,k)*cabs(i,nc)
                  pasm(i) = pasm(i) + cheb(nc,m,i,k)*casm(i,nc)
               enddo
               pabs(i) = pabs(i)*wetvol(i)*rhoh2o
               pabs(i) = max(0._r8,pabs(i))
               pabs(i) = min(pext(i),pabs(i))

               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
               palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)

               dopaer(i) = pext(i)*mass(i,k)
            end do

            ! Save aerosol optical depth at longest visible wavelength
            ! sum over layers
            if (savaervis) then
               ! aerosol extinction (/m)
               do i = 1, ncol
                  extinct(i,k) = extinct(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                  absorb(i,k)  = absorb(i,k) + pabs(i)*air_density(i,k)
                  aodvis(i)    = aodvis(i) + dopaer(i)
                  aodabs(i)    = aodabs(i) + pabs(i)*mass(i,k)
                  aodmode(i,m) = aodmode(i,m) + dopaer(i)
                  if (wetvol(i) > 1.e-40_r8) then
                     dustaodmode(i,m) = dustaodmode(i,m) + dopaer(i)*dustvol(i)/wetvol(i)
                  endif
                  ssavis(i) = ssavis(i) + dopaer(i)*palb(i)
               end do
            endif

            do i = 1, ncol

               if ((dopaer(i) <= -1.e-10) .or. (dopaer(i) >= 30.)) then

                  write(iulog,*) 'dopaer(', i, ',', k, ',', m, ',', lchnk, ')=', dopaer(i)
                  ! write(iulog,*) 'itab,jtab,ttab,utab=',itab(i),jtab(i),ttab(i),utab(i)
                  write(iulog,*) 'k=', k, ' pext=', pext(i), ' specext=', specpext(i)
                  write(iulog,*) 'wetvol=', wetvol(i), ' dryvol=', dryvol(i), ' watervol=', watervol(i)
                  ! write(iulog,*) 'cext=',(cext(i,l),l=1,ncoef)
                  ! write(iulog,*) 'crefin=',crefin(i)
                  write(iulog,*) 'nspec_amode(m)=', nspec_amode(m)
                  ! write(iulog,*) 'cheb=', (cheb(nc,m,i,k),nc=2,ncoef)
                  do l = 1, nspec_amode(m)
                     volf = specmmr(l,m)%val(i,k)/specdens(l,m)
                     write(iulog,*) 'l=', l, 'vol(l)=', volf
                     write(iulog,*) 'isw=', isw, 'specrefindex(isw)=', specrefindex(l,m)%val(isw)
                     write(iulog,*) 'specdens=', specdens(l,m)
                  end do

                  nerr_dopaer = nerr_dopaer + 1
                  if (nerr_dopaer >= nerrmax_dopaer) then
                     ! write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                     ! call endrun('exit from '//subname)
                  end if

               end if
            end do

            do i=1,ncol
               tauxar(i,k,isw) = tauxar(i,k,isw) + dopaer(i)
               wa(i,k,isw)     = wa(i,k,isw)     + dopaer(i)*palb(i)
               ga(i,k,isw)     = ga(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)
               fa(i,k,isw)     = fa(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
            end do

         end do ! pver

      end do ! sw bands

   end do ! ntot_amode

   ! Output visible band diagnostics
   do i = 1, ncol
      if (aodvis(i) > 1.e-10) then
         ssavis(i) = ssavis(i)/aodvis(i)
      else
         ssavis(i) = 0.925_r8
      endif
      do m = 1, ntot_amode
         colext(i,m) = 0.001*aodmode(i,m)/burden(i,m)
         ! if (aodmode(i,m) .gt. 1.) write(iulog,*) 'm,burden,tau,ext=', m, &
         !    burden(i,m), aodmode(i,m), colext(i,m)
      end do
   end do

   do i = 1, nnite
      extinct(idxnite(i),:) = fillvalue
      absorb(idxnite(i),:)  = fillvalue
      aodvis(idxnite(i))    = fillvalue
      aodabs(idxnite(i))    = fillvalue
      aodmode(idxnite(i),:) = fillvalue
      dustaodmode(idxnite(i),:) = fillvalue
      burden(idxnite(i),:)  = fillvalue
      ssavis(idxnite(i))    = fillvalue
   end do

   call outfld('EXTINCT',  extinct,          pcols, lchnk)
   call outfld('ABSORB',   absorb,           pcols, lchnk)
   call outfld('AODVIS',   aodvis,           pcols, lchnk)
   call outfld('AODABS',   aodabs,           pcols, lchnk)
   call outfld('AODMODE1', aodmode(:,1),     pcols, lchnk)
   call outfld('AODMODE2', aodmode(:,2),     pcols, lchnk)
   call outfld('AODMODE3', aodmode(:,3),     pcols, lchnk)
   call outfld('AODDUST1', dustaodmode(:,1), pcols, lchnk)
   call outfld('AODDUST2', dustaodmode(:,2), pcols, lchnk)
   call outfld('AODDUST3', dustaodmode(:,3), pcols, lchnk)
   call outfld('BURDEN1',  burden(:,1),      pcols, lchnk)
   call outfld('BURDEN2',  burden(:,2),      pcols, lchnk)
   call outfld('BURDEN3',  burden(:,3),      pcols, lchnk)
   call outfld('SSAVIS',   ssavis,           pcols, lchnk)

   if (cam_chempkg_is('trop_mam7')) then
      call outfld('AODMODE4', aodmode(:,4), pcols, lchnk)
      call outfld('AODMODE5', aodmode(:,5), pcols, lchnk)
      call outfld('AODMODE6', aodmode(:,6), pcols, lchnk)
      call outfld('AODMODE7', aodmode(:,7), pcols, lchnk)
      call outfld('BURDEN4',  burden(:,4),  pcols, lchnk)
      call outfld('BURDEN5',  burden(:,5),  pcols, lchnk)
      call outfld('BURDEN6',  burden(:,6),  pcols, lchnk)
      call outfld('BURDEN7',  burden(:,7),  pcols, lchnk)
   end if

   ! deallocate local storage
   deallocate( &
      radsurf,     &
      logradsurf,  &
      cheb,        &
      aodmode,     &
      dustaodmode, &
      burden,      &
      colext)

   deallocate( &
      specmmr,    &
      specdens,   &
      specrefindex)

end subroutine modal_aero_sw

!===============================================================================

subroutine modal_aero_lw(state, pbuf, tauxar)

   ! calculates aerosol lw radiative properties


   type(physics_state), intent(in)  :: state    ! state variables
   type(pbuf_fld),      intent(in)  :: pbuf(:)  ! physics buffer

   real(r8), intent(out) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth

   ! Local variables
   integer :: i, ifld, ilw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk

   real(r8), pointer :: dgnumwet(:,:,:)  ! wet number mode diameter (m)
   real(r8), pointer :: qaerwat(:,:,:)  ! aerosol water (g/g)

   real(r8) :: xrad(pcols)
   real(r8), pointer :: cheby(:,:,:,:)  ! chebychef polynomials

   real(r8) :: mass(pcols,pver) ! layer mass

   type(r_ptr2d_t), allocatable :: specmmr(:,:) ! species mass mixing ratio
   real(r8),        allocatable :: specdens(:,:) ! species density (kg/m3)
   type(c_ptr1d_t), allocatable :: specrefindex(:,:) ! species refractive index

   real(r8) :: vol(pcols)       ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)    ! volume concentration of aerosol mode (m3/kg)
   real(r8) :: wetvol(pcols)    ! volume concentration of wet mode (m3/kg)
   real(r8) :: watervol(pcols)  ! volume concentration of water in each mode (m3/kg)
   complex  :: crefin(pcols)    ! complex refractive index
   real(r8) :: refr(pcols)      ! real part of refractive index
   real(r8) :: refi(pcols)      ! imaginary part of refractive index

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8) :: cabs(pcols,ncoef)
   real(r8) :: pabs(pcols)      ! parameterized specific absorption (m2/kg)
   real(r8) :: dopaer(pcols)    ! aerosol optical depth in layer

   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf             ! volume fraction of insoluble aerosol

   character(len=*), parameter :: subname = 'modal_aero_lw'

   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! allocate local storage
   allocate(cheby(ncoef,pcols,pver,ntot_amode))

   ! pointers to physics buffer
   ifld = pbuf_get_fld_idx('DGNUMWET')
   dgnumwet => pbuf(ifld)%fld_ptr(1,:,:,lchnk,:)

   ifld = pbuf_get_fld_idx('QAERWAT')
   if (associated(pbuf(ifld)%fld_ptr)) then
      qaerwat => pbuf(ifld)%fld_ptr(1,:,:,lchnk,:)
   else
      call endrun(subname//': pbuf for QAERWAT not allocated' )
   end if

   ! calc size parameter for all columns
   ! this is the same calculation that's done in modal_size_parameters, but there
   ! some intermediate results are saved and the chebyshev polynomials are stored
   ! in a array with different index order.  Could be unified.
   do m = 1, ntot_amode
      do k = 1, pver
         do i = 1, ncol
            ! convert from number diameter to surface area
            xrad(i) = log(0.5*dgnumwet(i,k,m)) + 2.0_r8*alnsg_amode(m)*alnsg_amode(m)
            ! normalize size parameter
            xrad(i) = max(xrad(i), xrmin)
            xrad(i) = min(xrad(i), xrmax)
            xrad(i) = (2*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
            ! chebyshev polynomials
            cheby(1,i,k,m) = 1.0_r8
            cheby(2,i,k,m) = xrad(i)
            do nc = 3, ncoef
               cheby(nc,i,k,m) = 2.0_r8*xrad(i)*cheby(nc-1,i,k,m)-cheby(nc-2,i,k,m)
            end do
         end do
      end do
   end do

   ! initialize output variables
   tauxar(:ncol,:,:) = 0._r8

   ! dry mass in each cell
   mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

   ! access the mixing ratio and properties of the modal species
   allocate( &
      specmmr(nspec_max,ntot_amode),     &
      specdens(nspec_max,ntot_amode),    &
      specrefindex(nspec_max,ntot_amode) )

   do m = 1, ntot_amode
      do l = 1, nspec_amode(m)
         call rad_cnst_get_aer_mmr(0, spec_idx(l,m),  state, pbuf, specmmr(l,m)%val)
         call rad_cnst_get_aer_props(0, spec_idx(l,m), density_aer=specdens(l,m), &
                                     refindex_aer_lw=specrefindex(l,m)%val)
      end do
   end do

   ! loop over all aerosol modes
   do m = 1, ntot_amode

      do ilw = 1, nlwbands

         do k = 1, pver

            ! form bulk refractive index. Use volume mixing for infrared
            crefin(:ncol) = 0._r8
            dryvol(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec_amode(m)
               do i = 1, ncol
                  vol(i)    = specmmr(l,m)%val(i,k)/specdens(l,m)
                  dryvol(i) = dryvol(i) + vol(i)
                  crefin(i) = crefin(i) + vol(i)*specrefindex(l,m)%val(ilw)
               end do
            end do

            do i = 1, ncol
               watervol(i) = qaerwat(i,k,m)/rhoh2o
               wetvol(i)   = watervol(i) + dryvol(i)
               if (watervol(i) < 0.0_r8) then
                  if (abs(watervol(i)) .gt. 1.e-1*wetvol(i)) then
                     write(iulog,*) 'watervol,wetvol,dryvol=',watervol(i),wetvol(i),dryvol(i),' in '//subname
                     !   call endrun()
                  end if
                  watervol(i) = 0._r8
                  wetvol(i)   = dryvol(i)
               end if

               crefin(i) = crefin(i) + watervol(i)*crefwlw(ilw)
               if (wetvol(i) > 1.e-40) crefin(i) = crefin(i)/wetvol(i)
               refr(i) = real(crefin(i))
               refi(i) = aimag(crefin(i))
            end do

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(absplw(:,:,:,m,ilw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtablw(:,ilw), refitablw(:,ilw), &
                         itab, jtab, ttab, utab, cabs)

            ! parameterized optical properties
            do i = 1, ncol
               pabs(i) = 0.5*cabs(i,1)
               do nc = 2, ncoef
                  pabs(i) = pabs(i) + cheby(nc,i,k,m)*cabs(i,nc)
               end do
               pabs(i)   = pabs(i)*wetvol(i)*rhoh2o
               pabs(i)   = max(0._r8,pabs(i))
               dopaer(i) = pabs(i)*mass(i,k)
            end do

            do i = 1, ncol

               if ((dopaer(i) <= -1.e-10) .or. (dopaer(i) >= 20.)) then

                  write(iulog,*) 'dopaer(',i,',',k,',',m,',',lchnk,')=', dopaer(i)
                  write(iulog,*) 'k=',k,' pabs=', pabs(i)
                  write(iulog,*) 'wetvol=',wetvol(i),' dryvol=',dryvol(i),     &
                     ' watervol=',watervol(i)
                  write(iulog,*) 'cabs=', (cabs(i,l),l=1,ncoef)
                  write(iulog,*) 'crefin=', crefin(i)
                  write(iulog,*) 'nspec_amode(m)=', nspec_amode(m)
                  do l = 1,nspec_amode(m)
                     volf = specmmr(l,m)%val(i,k)/specdens(l,m)
                     write(iulog,*) 'l=',l,'vol(l)=',volf
                     write(iulog,*) 'ilw=',ilw,' specrefindex(ilw)=',specrefindex(l,m)%val(ilw)
                     write(iulog,*) 'specdens=',specdens(l,m)
                  end do

                  nerr_dopaer = nerr_dopaer + 1
                  if (nerr_dopaer >= nerrmax_dopaer) then
                     write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                     call endrun()
                  end if

               end if
            end do

            do i = 1, ncol
               tauxar(i,k,ilw) = tauxar(i,k,ilw) + dopaer(i)
            end do

         end do ! pver

      end do  ! nlwbands

   end do ! ntot_amode

   ! deallocate local storage
   deallocate(cheby)

   deallocate( &
      specmmr,    &
      specdens,   &
      specrefindex)

end subroutine modal_aero_lw

!===============================================================================
! Private routines
!===============================================================================

subroutine read_modal_optics(infilename)

#ifdef MODAL_AERO
   use modal_aero_data, only: sigmag_amode
#endif

   ! read modal optics file and set module data

   character*(*), intent(in) :: infilename   ! modal optics filename

   ! Local variables

   integer            :: i, ierr, istat, l, m
   type(file_desc_t)  :: ncid              ! pio file handle
   integer            :: did               ! dimension ids
   integer            :: dimread           ! dimension lengths
   type(var_desc_t)   :: vid               ! variable ids
   !----------------------------------------------------------------------------

   ! open file
   call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

   ! inquire dimensions.  Check that file values match parameter values.

   ierr = pio_inq_dimid(ncid,'mode',did)
   ierr = pio_inq_dimlen(ncid, did, ntot_amode)

   ierr = pio_inq_dimid(ncid,'nspec_max',did)
   ierr = pio_inq_dimlen(ncid, did, nspec_max)

   ierr = pio_inq_dimid(ncid, 'clen1', did)
   ierr = pio_inq_dimlen(ncid, did, dimread)
   if (dimread .ne. clen1) then
      write(iulog,*) 'clen1 len=', dimread, ' from ', infilename, ' ne clen1=', clen1
      call endrun('read_modal_optics: bad clen1 value')
   endif

   ierr = pio_inq_dimid(ncid, 'clen2', did)
   ierr = pio_inq_dimlen(ncid, did, dimread)
   if (dimread .ne. clen2) then
      write(iulog,*) 'clen2 len=', dimread, ' from ', infilename, ' ne clen2=', clen2
      call endrun('read_modal_optics: bad clen2 value')
   endif

   ierr = pio_inq_dimid(ncid, 'lw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimread)
   if (dimread .ne. nlwbands) then
      write(iulog,*) 'lw_band len=', dimread, ' from ', infilename, ' ne nlwbands=', nlwbands
      call endrun('read_modal_optics: bad lw_band value')
   endif

   ierr = pio_inq_dimid(ncid, 'sw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimread)
   if (dimread .ne. nswbands) then
      write(iulog,*) 'sw_band len=', dimread, ' from ', infilename, ' ne nswbands=', nswbands
      call endrun('read_modal_optics: bad sw_band value')
   endif

   ierr = pio_inq_dimid(ncid, 'refindex_real', did)
   ierr = pio_inq_dimlen(ncid, did, dimread)
   if (dimread .ne. prefr) then
      write(iulog,*) 'refindex_real=', dimread, ' from ', infilename, ' ne prefr=', prefr
      call endrun('read_modal_optics: bad refindex_real value')
   endif

   ierr = pio_inq_dimid(ncid, 'refindex_im', did)
   ierr = pio_inq_dimlen(ncid, did, dimread)
   if (dimread .ne. prefi) then
      write(iulog,*) 'refindex_im=', dimread, ' from ', infilename, ' ne prefi=', prefi
      call endrun('read_modal_optics: bad refindex_im value')
   endif

   ierr = pio_inq_dimid(ncid, 'coef_number', did)
   ierr = pio_inq_dimlen(ncid, did, dimread)
   if (dimread .ne. ncoef) then
      write(iulog,*) 'coef_number=',dimread,' from ',infilename,' ne ncoef=',ncoef
      call endrun('read_modal_optics: bad coef_number value')
   endif

   ! Allocate storage for module variables
   allocate( &
      nspec_amode(ntot_amode),                       &
      xname_massptr(nspec_max,ntot_amode),           &
      xname_spectype(nspec_max,ntot_amode),          &
      extpsw(ncoef,prefr,prefi,ntot_amode,nswbands), &
      abspsw(ncoef,prefr,prefi,ntot_amode,nswbands), &
      asmpsw(ncoef,prefr,prefi,ntot_amode,nswbands), &
      absplw(ncoef,prefr,prefi,ntot_amode,nlwbands), &
      sigma_logr_aer(ntot_amode),                    &
      stat=istat)
   if (istat > 0) call endrun('read_modal_optics: allocate failed')

   ! read variables
   ierr = pio_inq_varid(ncid, 'nspec_amode', vid)
   ierr = pio_get_var(ncid, vid, nspec_amode)

   ierr = pio_inq_varid(ncid, 'extpsw', vid)
   ierr = pio_get_var(ncid, vid, extpsw)

   ierr = pio_inq_varid(ncid, 'abspsw', vid)
   ierr = pio_get_var(ncid, vid, abspsw)

   ierr = pio_inq_varid(ncid, 'asmpsw', vid)
   ierr = pio_get_var(ncid, vid, asmpsw)

   ierr = pio_inq_varid(ncid, 'absplw', vid)
   ierr = pio_get_var(ncid, vid, absplw)

   ierr = pio_inq_varid(ncid, 'refindex_real_sw', vid)
   ierr = pio_get_var(ncid, vid, refrtabsw)

   ierr = pio_inq_varid(ncid, 'refindex_im_sw', vid)
   ierr = pio_get_var(ncid, vid, refitabsw)

   ierr = pio_inq_varid(ncid, 'refindex_real_lw', vid)
   ierr = pio_get_var(ncid, vid, refrtablw)

   ierr = pio_inq_varid(ncid, 'refindex_im_lw', vid)
   ierr = pio_get_var(ncid, vid, refitablw)

   ierr = pio_inq_varid(ncid, 'xname_massptr', vid)
   ierr = pio_get_var(ncid, vid, xname_massptr)

   ierr = pio_inq_varid(ncid, 'xname_spectype', vid)
   ierr = pio_get_var(ncid, vid, xname_spectype)

   ierr = pio_inq_varid(ncid, 'sigma_logr_aer', vid)
   ierr = pio_get_var(ncid, vid, sigma_logr_aer)

#ifdef MODAL_AERO
   ! This check is done for the prognostic version of MAM since sigmag_amode is
   ! initialized independently in modal_aero_initialize_data
   do m = 1, ntot_amode
      if (sigma_logr_aer(m)/sigmag_amode(m) < 0.99 .or. &
         sigma_logr_aer(m)/sigmag_amode(m) > 1.01) then
         write(iulog,*) 'sigma_logr_aer,sigmag_amode=', sigma_logr_aer(m), sigmag_amode(m)
         call endrun("read_modal_optics: data check failed")
      end if
   end do
#endif

   call pio_closefile(ncid)

end subroutine read_modal_optics

!===============================================================================

subroutine read_water_refindex(infilename)

   ! read water refractive index file and set module data

   character*(*), intent(in) :: infilename   ! modal optics filename

   ! Local variables

   integer            :: i, ierr
   type(file_desc_t)  :: ncid              ! pio file handle
   integer            :: did               ! dimension ids
   integer            :: dimlen            ! dimension lengths
   type(var_desc_t)   :: vid               ! variable ids
   real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
   real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared
   !----------------------------------------------------------------------------

   ! open file
   call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

   ! inquire dimensions.  Check that file values match parameter values.

   ierr = pio_inq_dimid(ncid, 'lw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nlwbands) then
      write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
      call endrun('read_modal_optics: bad lw_band value')
   endif

   ierr = pio_inq_dimid(ncid, 'sw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nswbands) then
      write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
      call endrun('read_modal_optics: bad sw_band value')
   endif

   ! read variables
   ierr = pio_inq_varid(ncid, 'refindex_real_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refrwsw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refiwsw)

   ierr = pio_inq_varid(ncid, 'refindex_real_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refrwlw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refiwlw)

   ! set complex representation of refractive indices as module data
   do i = 1, nswbands
      crefwsw(i)  = cmplx(refrwsw(i), abs(refiwsw(i)))
   end do
   do i = 1, nlwbands
      crefwlw(i)  = cmplx(refrwlw(i), abs(refiwlw(i)))
   end do

   call pio_closefile(ncid)

end subroutine read_water_refindex

!===============================================================================

subroutine modal_size_parameters(ncol, dgnumwet, radsurf, logradsurf, cheb)

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: dgnumwet(:,:,:)   ! aerosol wet number mode diameter (m)
   real(r8), intent(out) :: radsurf(:,:,:)    ! aerosol surface mode radius
   real(r8), intent(out) :: logradsurf(:,:,:) ! log(aerosol surface mode radius)
   real(r8), intent(out) :: cheb(:,:,:,:)

   integer  :: m, i, k, nc
   real(r8) :: explnsigma
   real(r8) :: xrad(pcols) ! normalized aerosol radius

   do m = 1, ntot_amode

      explnsigma = exp(2.0_r8*alnsg_amode(m)*alnsg_amode(m))

      do k = 1, pver
         do i = 1, ncol
            ! convert from number mode diameter to surface area
            radsurf(i,k,m) = 0.5*dgnumwet(i,k,m)*explnsigma
            logradsurf(i,k,m) = log(radsurf(i,k,m))
            ! normalize size parameter
            xrad(i) = max(logradsurf(i,k,m),xrmin)
            xrad(i) = min(xrad(i),xrmax)
            xrad(i) = (2._r8*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
            ! chebyshev polynomials
            cheb(1,m,i,k) = 1._r8
            cheb(2,m,i,k) = xrad(i)
            do nc = 3, ncoef
               cheb(nc,m,i,k) = 2._r8*xrad(i)*cheb(nc-1,m,i,k)-cheb(nc-2,m,i,k)
            end do
         end do
      end do

   end do

end subroutine modal_size_parameters

!===============================================================================

      subroutine binterp(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)

!     bilinear interpolation of table
!
      implicit none
      integer im,jm,km,ncol
      real(r8) table(km,im,jm),xtab(im),ytab(jm),out(pcols,km)
      integer i,ix(pcols),ip1,j,jy(pcols),jp1,k,ic
      real(r8) x(pcols),dx,t(pcols),y(pcols),dy,u(pcols), &
             tu(pcols),tuc(pcols),tcu(pcols),tcuc(pcols)

      if(ix(1).gt.0)go to 30
      if(im.gt.1)then
        do ic=1,ncol
          do i=1,im
            if(x(ic).lt.xtab(i))go to 10
          enddo
   10     ix(ic)=max0(i-1,1)
          ip1=min(ix(ic)+1,im)
          dx=(xtab(ip1)-xtab(ix(ic)))
          if(abs(dx).gt.1.e-20_r8)then
             t(ic)=(x(ic)-xtab(ix(ic)))/dx
          else
             t(ic)=0._r8
          endif
	end do
      else
        ix(:ncol)=1
        t(:ncol)=0._r8
      endif
      if(jm.gt.1)then
        do ic=1,ncol
          do j=1,jm
            if(y(ic).lt.ytab(j))go to 20
          enddo
   20     jy(ic)=max0(j-1,1)
          jp1=min(jy(ic)+1,jm)
          dy=(ytab(jp1)-ytab(jy(ic)))
          if(abs(dy).gt.1.e-20_r8)then
             u(ic)=(y(ic)-ytab(jy(ic)))/dy
             if(u(ic).lt.0._r8.or.u(ic).gt.1._r8)then
                write(iulog,*) 'u,y,jy,ytab,dy=',u(ic),y(ic),jy(ic),ytab(jy(ic)),dy
             endif
          else
            u(ic)=0._r8
          endif
	end do
      else
        jy(:ncol)=1
        u(:ncol)=0._r8
      endif
   30 continue
      do ic=1,ncol
         tu(ic)=t(ic)*u(ic)
         tuc(ic)=t(ic)-tu(ic)
         tcuc(ic)=1._r8-tuc(ic)-u(ic)
         tcu(ic)=u(ic)-tu(ic)
         jp1=min(jy(ic)+1,jm)
         ip1=min(ix(ic)+1,im)
         do k=1,km
            out(ic,k)=tcuc(ic)*table(k,ix(ic),jy(ic))+tuc(ic)*table(k,ip1,jy(ic))   &
               +tu(ic)*table(k,ip1,jp1)+tcu(ic)*table(k,ix(ic),jp1)
	 end do
      enddo
      return
      end subroutine binterp



end module modal_aer_opt
