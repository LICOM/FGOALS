module phys_prop

! Properties of aerosols that are used by radiation and other parameterizations.

! *****N.B.*****
! This module is a utility used by the rad_constituents module.  The properties stored
! here are meant to be accessed via that module.  This module knows nothing about how
! this data is associated with the constituents that are radiatively active or those that
! are being used for diagnostic calculations.  That is the responsibility of the 
! rad_constituents module.

use shr_kind_mod,   only: r8 => shr_kind_r8
use ioFileMod,      only: getfil
use cam_pio_utils,  only: cam_pio_openfile
use pio,            only: file_desc_t, var_desc_t, pio_get_var, pio_inq_varid, &
                          pio_inq_dimlen, pio_inq_dimid , pio_nowrite, pio_closefile, &
                          pio_seterrorhandling, PIO_BCAST_ERROR, PIO_INTERNAL_ERROR, PIO_NOERR

use cam_logfile,    only: iulog
use abortutils,     only: endrun
use spmd_utils,     only: masterproc
use radconstants,   only: nrh, nlwbands, nswbands, idx_sw_diag

implicit none
private
save

integer, parameter, public :: ot_length = 32
public :: &
   physprop_accum_unique_files,  &! Make a list of the unique set of files that contain properties
                                  ! This is an initialization step that must be done before calling physprop_init
   physprop_init,                &! Initialization -- read the input datasets
   physprop_get_id,              &! Return ID used to access the property data from the input files
   physprop_get                   ! Return data for specified ID

! Data from one input dataset is stored in a structure of type(physprop_type).
type :: physprop_type
   character(len=256) :: sourcefile ! Absolute pathname of data file.
   character(len=ot_length)  :: opticsmethod ! one of {hygro,nonhygro}

   ! for hygroscopic species of externally mixed aerosols
   real(r8), pointer :: sw_hygro_ext(:,:)
   real(r8), pointer :: sw_hygro_ssa(:,:)
   real(r8), pointer :: sw_hygro_asm(:,:)
   real(r8), pointer :: lw_hygro_abs(:,:)

   ! for nonhygroscopic species of externally mixed aerosols
   real(r8), pointer :: sw_nonhygro_ext(:)
   real(r8), pointer :: sw_nonhygro_ssa(:)
   real(r8), pointer :: sw_nonhygro_asm(:)
   real(r8), pointer :: sw_nonhygro_scat(:)
   real(r8), pointer :: sw_nonhygro_ascat(:)
   real(r8), pointer :: lw_abs(:)

   ! complex refractive index
   complex,  pointer :: refindex_aer_sw(:)
   complex,  pointer :: refindex_aer_lw(:)

   ! for radius-dependent mass-specific quantities
   real(r8), pointer :: r_sw_ext(:,:)
   real(r8), pointer :: r_sw_scat(:,:)
   real(r8), pointer :: r_sw_ascat(:,:)
   real(r8), pointer :: r_lw_abs(:,:)
   real(r8), pointer :: mu(:)

   ! microphysics parameters.
   character(len=32) :: aername ! for output of number concentration
   real(r8) :: density_aer      ! density of aerosol (kg/m3)
   real(r8) :: hygro_aer        ! hygroscopicity of aerosol
   real(r8) :: dryrad_aer       ! number mode radius (m) of aerosol size distribution
   real(r8) :: dispersion_aer   ! geometric standard deviation of aerosol size distribution
   real(r8) :: num_to_mass_aer  ! ratio of number concentration to mass concentration (#/kg)
                                ! *** Is this actually (kg/#) ???
endtype physprop_type

! This module stores data in an array of physprop_type structures.  The way this data
! is accessed outside the module is via a physprop ID, which is an index into the array.
integer :: numphysprops = 0 ! an incremental total across ALL clim and diag constituents
type (physprop_type), pointer :: physprop(:)

! Temporary storage location for filenames in namelist, and construction of dynamic index
! to properties.  The unique filenames specified in the namelist are the identifiers of
! the properties.  Searching the uniquefilenames array provides the index into the physprop
! array.
integer, parameter :: maxuniquefiles = 50
character(len=256) :: uniquefilenames(maxuniquefiles) 
 
!================================================================================================
contains
!================================================================================================

subroutine physprop_accum_unique_files(radname, type, numaerosols)

   ! Count number of aerosols in input radname array.  Aerosols are identified
   ! as strings with a ".nc" suffix.
   ! Construct a cumulative list of unique filenames containing physical property data.

   character(len=*), intent(in)  :: radname(:)
   character(len=1), intent(in)  :: type(:)
   integer,          intent(out) :: numaerosols 

   integer :: ncnst, i
   character(len=*), parameter :: subname = 'physprop_accum_unique_files'
   !------------------------------------------------------------------------------------
  
   ncnst = ubound(radname, 1)

   numaerosols = 0
   do i = 1, ncnst

      ! check if radname is an aerosol
      if (type(i) == 'A') then
         numaerosols = numaerosols + 1

         ! check if this filename has been used by another aerosol.  If not
         ! then add it to the list of unique names.
         if (physprop_get_id(radname(i)) < 0) then
            numphysprops = numphysprops + 1
            if (numphysprops > maxuniquefiles) then
               write(iulog,*) subname//': request for more than ',maxuniquefiles, ' values'
               call endrun(subname//': need to increase maxuniquefiles value')
            end if
            uniquefilenames(numphysprops) = trim(radname(i))
         endif

      endif
   enddo

end subroutine physprop_accum_unique_files

!================================================================================================

subroutine physprop_init()

   ! Read properties from the aerosol data files.

   ! ***N.B.*** The calls to physprop_accum_unique_files must be made before calling
   !            this init routine.  physprop_accum_unique_files is responsible for building
   !            the list of files to be read here.

   ! Local variables
   integer            :: fileindex
   type(file_desc_t)  :: nc_id ! index to netcdf file
   character(len=256) :: locfn ! path to actual file used
   character(len=32)  :: aername_str ! string read from netCDF file -- may contain trailing
                                     ! nulls which aren't dealt with by trim()
   
   type(var_desc_T) :: aername_id, density_aer_id, dispersion_aer_id, dryrad_aer_id, hygro_aer_id, num_to_mass_aer_id

   integer :: ierr ! error codes from mpi
   logical :: debug = .true.
   !------------------------------------------------------------------------------------

   allocate(physprop(numphysprops))

   do fileindex = 1, numphysprops
      nullify(physprop(fileindex)%sw_hygro_ext)
      nullify(physprop(fileindex)%sw_hygro_ssa)
      nullify(physprop(fileindex)%sw_hygro_asm)
      nullify(physprop(fileindex)%lw_hygro_abs)

      nullify(physprop(fileindex)%sw_nonhygro_ext)
      nullify(physprop(fileindex)%sw_nonhygro_ssa)
      nullify(physprop(fileindex)%sw_nonhygro_asm)
      nullify(physprop(fileindex)%sw_nonhygro_scat)
      nullify(physprop(fileindex)%sw_nonhygro_ascat)
      nullify(physprop(fileindex)%lw_abs)

      nullify(physprop(fileindex)%refindex_aer_sw)
      nullify(physprop(fileindex)%refindex_aer_lw)

      nullify(physprop(fileindex)%r_sw_ext)
      nullify(physprop(fileindex)%r_sw_scat)
      nullify(physprop(fileindex)%r_sw_ascat)
      nullify(physprop(fileindex)%r_lw_abs)
      nullify(physprop(fileindex)%mu)

      call getfil(uniquefilenames(fileindex), locfn, 0)
      physprop(fileindex)%sourcefile = locfn

      ! Open the physprop file
      call cam_pio_openfile(nc_id, locfn, PIO_NOWRITE)

      call aerosol_optics_init(physprop(fileindex), nc_id)

      ! read microphys
      ierr = pio_inq_varid(nc_id, 'name', aername_id)
      ierr = pio_get_var(nc_id, aername_id, physprop(fileindex)%aername)

      ! use GLC function to remove trailing nulls and blanks.
      ! physprop(fileindex)%aername = aername_str(:GLC(aername_str))

      ierr = pio_inq_varid(nc_id, 'density', density_aer_id)
      ierr = pio_get_var(nc_id, density_aer_id, physprop(fileindex)%density_aer)

      ierr = pio_inq_varid(nc_id, 'sigma_logr', dispersion_aer_id)
      ierr = pio_get_var(nc_id, dispersion_aer_id, physprop(fileindex)%dispersion_aer)

      ierr = pio_inq_varid(nc_id, 'dryrad', dryrad_aer_id)
      ierr = pio_get_var(nc_id, dryrad_aer_id, physprop(fileindex)%dryrad_aer)
         
      ierr = pio_inq_varid(nc_id, 'hygroscopicity', hygro_aer_id)
      ierr = pio_get_var(nc_id, hygro_aer_id, physprop(fileindex)%hygro_aer)

      ierr = pio_inq_varid(nc_id, 'num_to_mass_ratio', num_to_mass_aer_id)
      ierr = pio_get_var(nc_id, num_to_mass_aer_id, physprop(fileindex)%num_to_mass_aer)
      
      ! Output select data to log file
      if (debug .and. masterproc) then
         if (trim(physprop(fileindex)%aername) == 'SULFATE') then
            write(iulog, '(2x, a)') '_______ hygroscopic growth in visible band _______'
            call aer_optics_log_rh('SO4', physprop(fileindex)%sw_hygro_ext(:,idx_sw_diag), &
               physprop(fileindex)%sw_hygro_ssa(:,idx_sw_diag), physprop(fileindex)%sw_hygro_asm(:,idx_sw_diag))
         end if
         write(iulog, *) ' physprop_init finished for ', trim(physprop(fileindex)%aername)
      end if

      ! Close the physprop file
      call pio_closefile(nc_id)

   enddo ! fileindex
end subroutine physprop_init

!================================================================================================

integer function physprop_get_id(filename)

   ! Look for filename in the global list of unique filenames (module data uniquefilenames).
   ! If found, return it's index in the list.  Otherwise return -1.

   character(len=*), intent(in) :: filename
   integer iphysprop

   physprop_get_id = -1
   do iphysprop = 1, numphysprops
     if(trim(uniquefilenames(iphysprop)) == trim(filename) ) then
       physprop_get_id = iphysprop
       return
     endif
   enddo

end function physprop_get_id

!================================================================================================

subroutine physprop_get(id, sourcefile, opticstype, sw_hygro_ext, sw_hygro_ssa, &
   sw_hygro_asm, sw_nonhygro_ext, sw_nonhygro_ssa, sw_nonhygro_asm,             &
   sw_nonhygro_scat, sw_nonhygro_ascat, lw_abs,     &
   aername, density_aer, hygro_aer, dryrad_aer, dispersion_aer, lw_hygro_abs,   &
   refindex_aer_sw, refindex_aer_lw, &
   r_sw_ext, r_sw_scat, r_sw_ascat, r_lw_abs, mu, &
   num_to_mass_aer)

   ! Return requested properties for specified ID.

   ! Arguments
   integer,                     intent(in)  :: id
   character(len=256),optional, intent(out) :: sourcefile ! Absolute pathname of data file.
   character(len=ot_length), optional, intent(out) :: opticstype
   real(r8),          optional, pointer     :: sw_hygro_ext(:,:)
   real(r8),          optional, pointer     :: sw_hygro_ssa(:,:) 
   real(r8),          optional, pointer     :: sw_hygro_asm(:,:) 
   real(r8),          optional, pointer     :: lw_hygro_abs(:,:)         
   real(r8),          optional, pointer     :: sw_nonhygro_ext(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ssa(:)
   real(r8),          optional, pointer     :: sw_nonhygro_asm(:)
   real(r8),          optional, pointer     :: sw_nonhygro_scat(:)
   real(r8),          optional, pointer     :: sw_nonhygro_ascat(:)
   real(r8),          optional, pointer     :: lw_abs(:)         
   complex,           optional, pointer     :: refindex_aer_sw(:)
   complex,           optional, pointer     :: refindex_aer_lw(:)
   character(len=20), optional, intent(out) :: aername           
   real(r8),          optional, intent(out) :: density_aer       
   real(r8),          optional, intent(out) :: hygro_aer         
   real(r8),          optional, intent(out) :: dryrad_aer        
   real(r8),          optional, intent(out) :: dispersion_aer
   real(r8),          optional, intent(out) :: num_to_mass_aer
   real(r8),          optional, pointer     :: r_sw_ext(:,:)
   real(r8),          optional, pointer     :: r_sw_scat(:,:)
   real(r8),          optional, pointer     :: r_sw_ascat(:,:)
   real(r8),          optional, pointer     :: r_lw_abs(:,:)
   real(r8),          optional, pointer     :: mu(:)

   ! Local variables
   character(len=*), parameter :: subname = 'physprop_get'
   !------------------------------------------------------------------------------------

   if (id <= 0 .or. id > numphysprops) then
      write(iulog,*) subname//': illegal ID value: ', id
      call endrun('physprop_get: ID out of range')
   end if

   if (present(sourcefile))        sourcefile        =  physprop(id)%sourcefile
   if (present(opticstype))        opticstype        =  physprop(id)%opticsmethod
   if (present(sw_hygro_ext))      sw_hygro_ext      => physprop(id)%sw_hygro_ext
   if (present(sw_hygro_ssa))      sw_hygro_ssa      => physprop(id)%sw_hygro_ssa
   if (present(sw_hygro_asm))      sw_hygro_asm      => physprop(id)%sw_hygro_asm
   if (present(lw_hygro_abs))      lw_hygro_abs      => physprop(id)%lw_hygro_abs
   if (present(sw_nonhygro_ext))   sw_nonhygro_ext   => physprop(id)%sw_nonhygro_ext
   if (present(sw_nonhygro_ssa))   sw_nonhygro_ssa   => physprop(id)%sw_nonhygro_ssa
   if (present(sw_nonhygro_asm))   sw_nonhygro_asm   => physprop(id)%sw_nonhygro_asm
   if (present(sw_nonhygro_scat))  sw_nonhygro_scat  => physprop(id)%sw_nonhygro_scat
   if (present(sw_nonhygro_ascat)) sw_nonhygro_ascat => physprop(id)%sw_nonhygro_ascat
   if (present(lw_abs))            lw_abs            => physprop(id)%lw_abs

   if (present(refindex_aer_sw))   refindex_aer_sw   => physprop(id)%refindex_aer_sw
   if (present(refindex_aer_lw))   refindex_aer_lw   => physprop(id)%refindex_aer_lw

   if (present(aername))         aername         =  physprop(id)%aername
   if (present(density_aer))     density_aer     =  physprop(id)%density_aer
   if (present(hygro_aer))       hygro_aer       =  physprop(id)%hygro_aer
   if (present(dryrad_aer))      dryrad_aer      =  physprop(id)%dryrad_aer
   if (present(dispersion_aer))  dispersion_aer  =  physprop(id)%dispersion_aer
   if (present(num_to_mass_aer)) num_to_mass_aer =  physprop(id)%num_to_mass_aer

   ! radius-dependent optics
   if (present(r_sw_ext))         r_sw_ext      => physprop(id)%r_sw_ext
   if (present(r_sw_scat))        r_sw_scat     => physprop(id)%r_sw_scat
   if (present(r_sw_ascat))       r_sw_ascat    => physprop(id)%r_sw_ascat
   if (present(r_lw_abs))         r_lw_abs      => physprop(id)%r_lw_abs
   if (present(mu))               mu            => physprop(id)%mu

end subroutine physprop_get

!================================================================================================
! Private methods
!================================================================================================

subroutine aerosol_optics_init(phys_prop, nc_id)

   ! Determine the opticstype, then call the 
   ! appropriate routine to read the data.

   type(physprop_type), intent(inout) :: phys_prop  ! data after interp onto cam rh mesh
   type(file_desc_t),   intent(inout) :: nc_id      ! indentifier for netcdf file

   integer :: opticslength_id, opticslength
   type(var_desc_t) :: op_type_id
   integer :: ierr ! mpi error codes
   character(len=ot_length)  :: opticstype_str ! string read from netCDF file -- may contain trailing
                                        ! nulls which aren't dealt with by trim()
   !------------------------------------------------------------------------------------

   ierr = pio_inq_dimid(nc_id, 'opticsmethod_len', opticslength_id)
   ierr = pio_inq_dimlen(nc_id, opticslength_id, opticslength)
   if ( opticslength .gt. ot_length ) then
      call endrun(" optics type length in "//phys_prop%sourcefile//" excedes maximum length of 32")
   endif
   ierr = pio_inq_varid(nc_id, 'opticsmethod', op_type_id)
   ierr = pio_get_var(nc_id, op_type_id,phys_prop%opticsmethod )

   select case (phys_prop%opticsmethod)
   case ('zero')
      call zero_optics_init(phys_prop, nc_id)

   case ('hygro')
      call hygro_optics_init(phys_prop, nc_id)

   case ('hygroscopic')
      call hygroscopic_optics_init(phys_prop, nc_id)

   case ('nonhygro')
      call nonhygro_optics_init(phys_prop, nc_id)
        
   case ('insoluble')
      call insoluble_optics_init(phys_prop, nc_id)
        
   case ('volcanic_radius')
      call volcanic_radius_optics_init(phys_prop, nc_id)

   case ('volcanic')
      call volcanic_optics_init(phys_prop, nc_id)
        
   ! other types of optics can be added here

   case default
      call endrun('aerosol_optics_init: unsupported optics type '//&
         phys_prop%opticsmethod//' in file '//phys_prop%sourcefile)
   end select

end subroutine aerosol_optics_init

!================================================================================================

subroutine hygro_optics_init(phys_prop, nc_id)

   ! Read optics data of type 'hygro' and interpolate it to CAM's rh mesh.

   type (physprop_type), intent(inout) :: phys_prop  ! data after interp onto cam rh mesh
   type (file_desc_t),   intent(inout) :: nc_id      ! indentifier for netcdf file

   ! Local variables
   integer :: ierr ! error flag

   integer :: rh_idx_id, lw_band_id, sw_band_id
   integer :: kbnd, krh
   integer :: rh_id, sw_ext_id, sw_ssa_id, sw_asm_id, lw_ext_id
   integer :: nbnd, swbands

   ! temp data from hygroscopic file before interpolation onto cam-rh-mesh
   integer  :: nfilerh ! number of rh values in file
   real(r8), allocatable, dimension(:) :: frh
   real(r8), allocatable, dimension(:,:)  :: fsw_ext
   real(r8), allocatable, dimension(:,:)  :: fsw_ssa
   real(r8), allocatable, dimension(:,:)  :: fsw_asm

   real(r8) :: rh ! real rh value on cam rh mesh (indexvalue)
   !------------------------------------------------------------------------------------

   allocate(phys_prop%sw_hygro_ext(nrh,nswbands))
   allocate(phys_prop%sw_hygro_ssa(nrh,nswbands))
   allocate(phys_prop%sw_hygro_asm(nrh,nswbands))
   allocate(phys_prop%lw_abs(nlwbands))

   ierr = pio_inq_dimid(nc_id, 'rh_idx', rh_idx_id)

   ierr = pio_inq_dimlen(nc_id, rh_idx_id, nfilerh)

   ierr = pio_inq_dimid(nc_id, 'lw_band', lw_band_id)

   ierr = pio_inq_dimid(nc_id, 'sw_band', sw_band_id)

   ierr = pio_inq_dimlen(nc_id, lw_band_id, nbnd)

   if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of lwbands')

   ierr = pio_inq_dimlen(nc_id, sw_band_id, swbands)

   if(swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
         ' has the wrong number of sw bands')

   ierr = pio_inq_varid(nc_id, 'rh', rh_id)

   ierr = pio_inq_varid(nc_id, 'ext_sw', sw_ext_id)

   ierr = pio_inq_varid(nc_id, 'ssa_sw', sw_ssa_id)

   ierr = pio_inq_varid(nc_id, 'asm_sw', sw_asm_id)

   ierr = pio_inq_varid(nc_id, 'abs_lw', lw_ext_id)

   ! specific optical properties on file's rh mesh
   allocate(fsw_ext(nfilerh,nswbands))
   allocate(fsw_asm(nfilerh,nswbands))
   allocate(fsw_ssa(nfilerh,nswbands))
   allocate(frh(nfilerh))

   ierr = pio_get_var(nc_id, rh_id, frh)

   ierr = pio_get_var(nc_id, sw_ext_id, fsw_ext)

   ierr = pio_get_var(nc_id, sw_ssa_id, fsw_ssa)

   ierr = pio_get_var(nc_id, sw_asm_id, fsw_asm)

   ierr = pio_get_var(nc_id, lw_ext_id, phys_prop%lw_abs)

   ! interpolate onto cam's rh mesh
   do kbnd = 1,nswbands
      do krh = 1, nrh
         rh = 1.0_r8 / nrh * (krh - 1)
         phys_prop%sw_hygro_ext(krh,kbnd) = &
            exp_interpol( frh, fsw_ext(:,kbnd) / fsw_ext(1,kbnd), rh ) &
            * fsw_ext(1, kbnd)
         phys_prop%sw_hygro_ssa(krh,kbnd) = &
            lin_interpol( frh, fsw_ssa(:,kbnd) / fsw_ssa(1,kbnd), rh ) &
            * fsw_ssa(1, kbnd)
         phys_prop%sw_hygro_asm(krh,kbnd) = &
            lin_interpol( frh, fsw_asm(:,kbnd) / fsw_asm(1,kbnd), rh ) &
            * fsw_asm(1, kbnd)
      enddo
   enddo

   deallocate (fsw_ext, fsw_asm, fsw_ssa, frh)

   ! read refractive index data if available
   call refindex_aer_init(phys_prop, nc_id)

end subroutine hygro_optics_init

!================================================================================================

subroutine zero_optics_init(phys_prop, nc_id)

   ! Read optics data of type 'nonhygro'

   type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
   type (file_desc_t),   intent(inout) :: nc_id      ! indentifier for netcdf file

   ! Local variables
   integer :: lw_band_id, sw_band_id
   integer :: sw_ext_id, sw_ssa_id, sw_asm_id, lw_ext_id
   integer :: swbands, nbnd
   integer :: ierr ! error flag
   !------------------------------------------------------------------------------------

   ! perhaps this doesn't even need allocated.
   allocate (phys_prop%sw_nonhygro_ext(nswbands))
   allocate (phys_prop%sw_nonhygro_ssa(nswbands))
   allocate (phys_prop%sw_nonhygro_asm(nswbands))
   allocate (phys_prop%lw_abs(nlwbands))

   phys_prop%sw_nonhygro_ext = 0._r8
   phys_prop%sw_nonhygro_ssa = 0._r8
   phys_prop%sw_nonhygro_asm = 0._r8
   phys_prop%lw_abs = 0._r8

end subroutine zero_optics_init

!================================================================================================

subroutine insoluble_optics_init(phys_prop, nc_id)

   ! Read optics data of type 'nonhygro'

   type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
   type (file_desc_t),   intent(inout) :: nc_id      ! indentifier for netcdf file

   ! Local variables
   integer :: lw_band_id, sw_band_id
   integer :: sw_ext_id, sw_ssa_id, sw_asm_id, lw_ext_id
   integer :: swbands, nbnd
   integer :: ierr ! error flag
   integer :: start(2), count(2)
   !------------------------------------------------------------------------------------

   allocate (phys_prop%sw_nonhygro_ext(nswbands))
   allocate (phys_prop%sw_nonhygro_ssa(nswbands))
   allocate (phys_prop%sw_nonhygro_asm(nswbands))
   allocate (phys_prop%lw_abs(nlwbands))

   ierr = pio_inq_dimid(nc_id, 'lw_band', lw_band_id)

   ierr = pio_inq_dimid(nc_id, 'sw_band', sw_band_id)

   ierr = pio_inq_dimlen(nc_id, lw_band_id, nbnd)

   if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of lwbands')

   ierr = pio_inq_dimlen(nc_id, sw_band_id, swbands)

   if (swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of sw bands')

   ! read file data
   ierr = pio_inq_varid(nc_id, 'ext_sw', sw_ext_id)
   ierr = pio_inq_varid(nc_id, 'ssa_sw', sw_ssa_id)
   ierr = pio_inq_varid(nc_id, 'asm_sw', sw_asm_id)
   ierr = pio_inq_varid(nc_id, 'abs_lw', lw_ext_id)

   start = 1
   count=(/1,swbands/)

   ierr = pio_get_var(nc_id, sw_ext_id, start, count, phys_prop%sw_nonhygro_ext)
   ierr = pio_get_var(nc_id, sw_ssa_id, start, count, phys_prop%sw_nonhygro_ssa)
   ierr = pio_get_var(nc_id, sw_asm_id, start, count, phys_prop%sw_nonhygro_asm)
   count = (/1,nbnd/)
   ierr = pio_get_var(nc_id, lw_ext_id, start, count, phys_prop%lw_abs)

   ! read refractive index data if available
   call refindex_aer_init(phys_prop, nc_id)

end subroutine insoluble_optics_init

subroutine volcanic_radius_optics_init(phys_prop, nc_id)

   ! Read optics data of type 'volcanic_radius'

   type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
   type (file_desc_t),   intent(inout) :: nc_id      ! indentifier for netcdf file

   ! Local variables
   integer :: lw_band_id, sw_band_id, mu_id, mu_did
   integer :: sw_ext_id, sw_scat_id, sw_ascat_id, lw_abs_id
   integer :: swbands, nbnd, n_mu_samples
   integer :: ierr ! error flag
   !------------------------------------------------------------------------------------

   ierr = pio_inq_dimid(nc_id, 'mu_samples', mu_did)
   ierr = pio_inq_dimlen(nc_id, mu_did, n_mu_samples)

   allocate (phys_prop%r_sw_ext(nswbands,n_mu_samples))
   allocate (phys_prop%r_sw_scat(nswbands,n_mu_samples))
   allocate (phys_prop%r_sw_ascat(nswbands,n_mu_samples))
   allocate (phys_prop%r_lw_abs(nlwbands,n_mu_samples))
   allocate (phys_prop%mu(n_mu_samples))

   ierr = pio_inq_dimid(nc_id, 'lw_band', lw_band_id)

   ierr = pio_inq_dimid(nc_id, 'sw_band', sw_band_id)

   ierr = pio_inq_dimlen(nc_id, lw_band_id, nbnd)

   if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of lwbands')

   ierr = pio_inq_dimlen(nc_id, sw_band_id, swbands)

   if (swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of sw bands')

   ! read file data
   ierr = pio_inq_varid(nc_id, 'bext_sw', sw_ext_id)
   ierr = pio_inq_varid(nc_id, 'bsca_sw', sw_scat_id)
   ierr = pio_inq_varid(nc_id, 'basc_sw', sw_ascat_id)
   ierr = pio_inq_varid(nc_id, 'babs_lw', lw_abs_id)
   ierr = pio_inq_varid(nc_id, 'mu_samples', mu_id)

   ierr = pio_get_var(nc_id, sw_ext_id, phys_prop%r_sw_ext)
   ierr = pio_get_var(nc_id, sw_scat_id, phys_prop%r_sw_scat)
   ierr = pio_get_var(nc_id, sw_ascat_id, phys_prop%r_sw_ascat)
   ierr = pio_get_var(nc_id, lw_abs_id, phys_prop%r_lw_abs)
   ierr = pio_get_var(nc_id, mu_id, phys_prop%mu)

end subroutine volcanic_radius_optics_init

subroutine volcanic_optics_init(phys_prop, nc_id)

   ! Read optics data of type 'volcanic'

   type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
   type (file_desc_t)  , intent(inout) :: nc_id      ! indentifier for netcdf file

   ! Local variables
   integer :: lw_band_id, sw_band_id
   integer :: sw_ext_id, sw_scat_id, sw_ascat_id, lw_abs_id
   integer :: swbands, nbnd
   integer :: ierr ! error flag
   !------------------------------------------------------------------------------------

   allocate (phys_prop%sw_nonhygro_ext(nswbands))
   allocate (phys_prop%sw_nonhygro_scat(nswbands))
   allocate (phys_prop%sw_nonhygro_ascat(nswbands))
   allocate (phys_prop%lw_abs(nlwbands))

   ierr = pio_inq_dimid(nc_id, 'lw_band', lw_band_id)
   ierr = pio_inq_dimid(nc_id, 'sw_band', sw_band_id)

   ierr = pio_inq_dimlen(nc_id, lw_band_id, nbnd)

   if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of lwbands')

   ierr = pio_inq_dimlen(nc_id, sw_band_id, swbands)
   if(masterproc) write(iulog,*) 'swbands',swbands

   if (swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of sw bands')

   ! read file data
   ierr = pio_inq_varid(nc_id, 'bext_sw', sw_ext_id)
   ierr = pio_inq_varid(nc_id, 'bsca_sw', sw_scat_id)
   ierr = pio_inq_varid(nc_id, 'basc_sw', sw_ascat_id)
   ierr = pio_inq_varid(nc_id, 'babs_lw', lw_abs_id)

   ierr = pio_get_var(nc_id, sw_ext_id, phys_prop%sw_nonhygro_ext)
   ierr = pio_get_var(nc_id, sw_scat_id, phys_prop%sw_nonhygro_scat)
   ierr = pio_get_var(nc_id, sw_ascat_id, phys_prop%sw_nonhygro_ascat)
   ierr = pio_get_var(nc_id, lw_abs_id, phys_prop%lw_abs)

end subroutine volcanic_optics_init

subroutine hygroscopic_optics_init(phys_prop, nc_id)

   ! Read optics data of type 'hygroscopic' and interpolate it to CAM's rh mesh.

   type (physprop_type), intent(inout) :: phys_prop  ! data after interp onto cam rh mesh
   type (file_desc_T),   intent(inout) :: nc_id      ! indentifier for netcdf file

   ! Local variables
   integer :: ierr ! error flag

   integer :: rh_idx_id, lw_band_id, sw_band_id
   integer :: kbnd, krh
   integer :: rh_id, sw_ext_id, sw_ssa_id, sw_asm_id, lw_ext_id
   integer :: nbnd, swbands

   ! temp data from hygroscopic file before interpolation onto cam-rh-mesh
   integer  :: nfilerh ! number of rh values in file
   real(r8), allocatable, dimension(:) :: frh
   real(r8), allocatable, dimension(:,:)  :: fsw_ext
   real(r8), allocatable, dimension(:,:)  :: fsw_ssa
   real(r8), allocatable, dimension(:,:)  :: fsw_asm
   real(r8), allocatable, dimension(:,:)  :: flw_abs

   real(r8) :: rh ! real rh value on cam rh mesh (indexvalue)
   character(len=*), parameter :: sub = 'hygroscopic_optics_init'
   !------------------------------------------------------------------------------------

   allocate(phys_prop%sw_hygro_ext(nrh,nswbands))
   allocate(phys_prop%sw_hygro_ssa(nrh,nswbands))
   allocate(phys_prop%sw_hygro_asm(nrh,nswbands))
   allocate(phys_prop%lw_hygro_abs(nrh,nlwbands))

   ierr = pio_inq_dimid(nc_id, 'rh_idx', rh_idx_id)
   ierr = pio_inq_dimlen(nc_id, rh_idx_id, nfilerh)
   ierr = pio_inq_dimid(nc_id, 'lw_band', lw_band_id)
   ierr = pio_inq_dimid(nc_id, 'sw_band', sw_band_id)
   ierr = pio_inq_dimlen(nc_id, lw_band_id, nbnd)

   if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of lwbands')

   ierr = pio_inq_dimlen(nc_id, sw_band_id, swbands)
   if(swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of sw bands')

   ierr = pio_inq_varid(nc_id, 'rh', rh_id)
   ierr = pio_inq_varid(nc_id, 'ext_sw', sw_ext_id)
   ierr = pio_inq_varid(nc_id, 'ssa_sw', sw_ssa_id)
   ierr = pio_inq_varid(nc_id, 'asm_sw', sw_asm_id)
   ierr = pio_inq_varid(nc_id, 'abs_lw', lw_ext_id)

   ! specific optical properties on file's rh mesh
   allocate(fsw_ext(nfilerh,nswbands))
   allocate(fsw_asm(nfilerh,nswbands))
   allocate(fsw_ssa(nfilerh,nswbands))
   allocate(flw_abs(nfilerh,nlwbands))
   allocate(frh(nfilerh))

   ierr = pio_get_var(nc_id, rh_id, frh)
   ierr = pio_get_var(nc_id, sw_ext_id, fsw_ext)
   ierr = pio_get_var(nc_id, sw_ssa_id, fsw_ssa)
   ierr = pio_get_var(nc_id, sw_asm_id, fsw_asm)
   ierr = pio_get_var(nc_id, lw_ext_id, flw_abs)

   ! interpolate onto cam's rh mesh
   do kbnd = 1,nswbands
      do krh = 1, nrh
         rh = 1.0_r8 / nrh * (krh - 1)
         phys_prop%sw_hygro_ext(krh,kbnd) = &
            exp_interpol( frh, fsw_ext(:,kbnd) / fsw_ext(1,kbnd), rh ) &
            * fsw_ext(1, kbnd)
         phys_prop%sw_hygro_ssa(krh,kbnd) = &
            lin_interpol( frh, fsw_ssa(:,kbnd) / fsw_ssa(1,kbnd), rh ) &
            * fsw_ssa(1, kbnd)
         phys_prop%sw_hygro_asm(krh,kbnd) = &
            lin_interpol( frh, fsw_asm(:,kbnd) / fsw_asm(1,kbnd), rh ) &
            * fsw_asm(1, kbnd)
      enddo
   enddo
   do kbnd = 1,nlwbands
      do krh = 1, nrh
         rh = 1.0_r8 / nrh * (krh - 1)
         phys_prop%lw_hygro_abs(krh,kbnd) = &
            exp_interpol( frh, flw_abs(:,kbnd) / flw_abs(1,kbnd), rh ) &
            * flw_abs(1, kbnd)
      enddo
   enddo

   deallocate (fsw_ext, fsw_asm, fsw_ssa, flw_abs, frh)

   ! read refractive index data if available
   call refindex_aer_init(phys_prop, nc_id)

end subroutine hygroscopic_optics_init

!================================================================================================

subroutine nonhygro_optics_init(phys_prop, nc_id)

   ! Read optics data of type 'nonhygro'

   type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
   type (file_desc_t)  , intent(inout) :: nc_id      ! indentifier for netcdf file

   ! Local variables
   integer :: lw_band_id, sw_band_id
   integer :: sw_ext_id, sw_ssa_id, sw_asm_id, lw_ext_id
   integer :: swbands, nbnd
   integer :: ierr ! error flag
   !------------------------------------------------------------------------------------

   allocate (phys_prop%sw_nonhygro_ext(nswbands))
   allocate (phys_prop%sw_nonhygro_ssa(nswbands))
   allocate (phys_prop%sw_nonhygro_asm(nswbands))
   allocate (phys_prop%lw_abs(nlwbands))

   ierr = pio_inq_dimid(nc_id, 'lw_band', lw_band_id)
   ierr = pio_inq_dimid(nc_id, 'sw_band', sw_band_id)

   ierr = pio_inq_dimlen(nc_id, lw_band_id, nbnd)

   if (nbnd .ne. nlwbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of lwbands')

   ierr = pio_inq_dimlen(nc_id, sw_band_id, swbands)

   if (swbands .ne. nswbands) call endrun(phys_prop%sourcefile// &
        ' has the wrong number of sw bands')

   ! read file data
   ierr = pio_inq_varid(nc_id, 'ext_sw', sw_ext_id)
   ierr = pio_inq_varid(nc_id, 'ssa_sw', sw_ssa_id)
   ierr = pio_inq_varid(nc_id, 'asm_sw', sw_asm_id)
   ierr = pio_inq_varid(nc_id, 'abs_lw', lw_ext_id)

   ierr = pio_get_var(nc_id, sw_ext_id, phys_prop%sw_nonhygro_ext)
   ierr = pio_get_var(nc_id, sw_ssa_id, phys_prop%sw_nonhygro_ssa)
   ierr = pio_get_var(nc_id, sw_asm_id, phys_prop%sw_nonhygro_asm)
   ierr = pio_get_var(nc_id, lw_ext_id, phys_prop%lw_abs)

   ! read refractive index data if available
   call refindex_aer_init(phys_prop, nc_id)

end subroutine nonhygro_optics_init

subroutine refindex_aer_init(phys_prop, nc_id)

!  Read refractive indices of aerosol

   type (physprop_type), intent(inout) :: phys_prop  ! storage for file data
   type (file_desc_T),   intent(inout) :: nc_id      ! indentifier for netcdf file

   ! Local variables
   integer :: i
   integer :: istat1, istat2, istat3     ! status flags
   integer :: vid_real, vid_im           ! variable ids
   real(r8), pointer :: ref_real(:), ref_im(:)  ! tmp storage for components of complex index
   character(len=*), parameter :: subname = 'refindex_aer_init'
   !------------------------------------------------------------------------------------

   ! assume that the dimensions lw_band and sw_band have already been checked
   ! by the calling subroutine

   ! Check that the variables are present before allocating storage and reading.
   ! Since we're setting complex data values, both the real and imaginary parts must
   ! be present or neither will be read.

   ! set PIO to return control to the caller when variable not found
   call pio_seterrorhandling(nc_id, PIO_BCAST_ERROR)

   istat1 = pio_inq_varid(nc_id, 'refindex_real_aer_sw', vid_real)
   istat2 = pio_inq_varid(nc_id, 'refindex_im_aer_sw',   vid_im)

   if (istat1 == PIO_NOERR  .and. istat2 == PIO_NOERR) then

      allocate(ref_real(nswbands), ref_im(nswbands))

      istat3 = pio_get_var(nc_id, vid_real, ref_real)
      if (istat3 /= PIO_NOERR) then
         call endrun(subname//': ERROR reading refindex_real_aer_sw')
      end if

      istat3 = pio_get_var(nc_id, vid_im, ref_im)
      if (istat3 /= PIO_NOERR) then
         call endrun(subname//': ERROR reading refindex_im_aer_sw')
      end if

      ! successfully read refindex data -- set complex values in physprop object
      allocate(phys_prop%refindex_aer_sw(nswbands))
      do i = 1, nswbands
         phys_prop%refindex_aer_sw(i) = cmplx(ref_real(i), abs(ref_im(i)))
      end do

      deallocate(ref_real, ref_im)

   end if

   istat1 = pio_inq_varid(nc_id, 'refindex_real_aer_lw', vid_real)
   istat2 = pio_inq_varid(nc_id, 'refindex_im_aer_lw',   vid_im)

   if (istat1 == PIO_NOERR  .and. istat2 == PIO_NOERR) then

      allocate(ref_real(nlwbands), ref_im(nlwbands))

      istat3 = pio_get_var(nc_id, vid_real, ref_real)
      if (istat3 /= PIO_NOERR) then
         call endrun(subname//': ERROR reading refindex_real_aer_lw')
      end if

      istat3 = pio_get_var(nc_id, vid_im, ref_im)
      if (istat3 /= PIO_NOERR) then
         call endrun(subname//': ERROR reading refindex_im_aer_lw')
      end if

      ! successfully read refindex data -- set complex value in physprop object
      allocate(phys_prop%refindex_aer_lw(nlwbands))
      do i = 1, nlwbands
         phys_prop%refindex_aer_lw(i) = cmplx(ref_real(i), abs(ref_im(i)))
      end do

      deallocate(ref_real, ref_im)

   end if

   ! reset PIO to handle errors internally
   call pio_seterrorhandling(nc_id, PIO_INTERNAL_ERROR)

end subroutine refindex_aer_init

!================================================================================================

function exp_interpol(x, f, y) result(g)
! Purpose:
!   interpolates f(x) to point y
!   assuming f(x) = f(x0) exp a(x - x0)
!   where a = ( ln f(x1) - ln f(x0) ) / (x1 - x0)
!   x0 <= x <= x1
!   assumes x is monotonically increasing
! Author: D. Fillmore

   implicit none

   real(r8), intent(in), dimension(:) :: x  ! grid points
   real(r8), intent(in), dimension(:) :: f  ! grid function values
   real(r8), intent(in) :: y                ! interpolation point
   real(r8) :: g                            ! interpolated function value

   integer :: k  ! interpolation point index
   integer :: n  ! length of x
   real(r8) :: a

   n = size(x)

   ! find k such that x(k) < y =< x(k+1)
   ! set k = 1 if y <= x(1)  and  k = n-1 if y > x(n)

   if (y <= x(1)) then
     k = 1
   else if (y >= x(n)) then
     k = n - 1
   else
     k = 1
     do while (y > x(k+1) .and. k < n)
       k = k + 1
     end do
   end if

   ! interpolate
   a = (  log( f(k+1) / f(k) )  ) / ( x(k+1) - x(k) )
   g = f(k) * exp( a * (y - x(k)) )
   return
end function exp_interpol

!================================================================================================

function lin_interpol(x, f, y) result(g)
! Purpose:
!   interpolates f(x) to point y
!   assuming f(x) = f(x0) + a * (x - x0)
!   where a = ( f(x1) - f(x0) ) / (x1 - x0)
!   x0 <= x <= x1
!   assumes x is monotonically increasing
! Author: D. Fillmore

   implicit none

   real(r8), intent(in), dimension(:) :: x  ! grid points
   real(r8), intent(in), dimension(:) :: f  ! grid function values
   real(r8), intent(in) :: y                ! interpolation point
   real(r8) :: g                            ! interpolated function value

   integer :: k  ! interpolation point index
   integer :: n  ! length of x
   real(r8) :: a

   n = size(x)

   ! find k such that x(k) < y =< x(k+1)
   ! set k = 1 if y <= x(1)  and  k = n-1 if y > x(n)

   if (y <= x(1)) then
     k = 1
   else if (y >= x(n)) then
     k = n - 1
   else
     k = 1
     do while (y > x(k+1) .and. k < n)
       k = k + 1
     end do
   end if

   ! interpolate
   a = (  f(k+1) - f(k) ) / ( x(k+1) - x(k) )
   g = f(k) + a * (y - x(k))
   return
end function lin_interpol

!================================================================================================

subroutine aer_optics_log(name, ext, ssa, asm)

   ! Purpose:
   !   write aerosol optical constants to log file

   ! Author: D. Fillmore

   character(len=*), intent(in) :: name
   real(r8), intent(in) :: ext(:)
   real(r8), intent(in) :: ssa(:)
   real(r8), intent(in) :: asm(:)

   integer :: kbnd, nbnd
   !------------------------------------------------------------------------------------

   nbnd = ubound(ext, 1)

   write(iulog, '(2x, a)') name
   write(iulog, '(2x, a, 4x, a, 4x, a, 4x, a)') 'SW band', 'ext (m^2 kg^-1)', ' ssa', ' asm'
   do kbnd = 1, nbnd
      write(iulog, '(2x, i7, 4x, f13.2, 4x, f4.2, 4x, f4.2)') kbnd, ext(kbnd), ssa(kbnd), asm(kbnd)
   end do

end subroutine aer_optics_log

!================================================================================================


subroutine aer_optics_log_rh(name, ext, ssa, asm)

   ! Purpose:
   !   write out aerosol optical properties
   !   for a set of test rh values
   !   to test hygroscopic growth interpolation

   ! Author: D. Fillmore

   character(len=*), intent(in) :: name
   real(r8), intent(in) :: ext(nrh)
   real(r8), intent(in) :: ssa(nrh)
   real(r8), intent(in) :: asm(nrh)

   integer :: krh_test
   integer, parameter :: nrh_test = 36
   integer :: krh
   real(r8) :: rh
   real(r8) :: rh_test(nrh_test)
   real(r8) :: exti
   real(r8) :: ssai
   real(r8) :: asmi
   real(r8) :: wrh
   !------------------------------------------------------------------------------------

   do krh_test = 1, nrh_test
      rh_test(krh_test) = sqrt(sqrt(sqrt(sqrt(((krh_test - 1.0_r8) / (nrh_test - 1))))))
   enddo
   write(iulog, '(2x, a)') name
   write(iulog, '(2x, a, 4x, a, 4x, a, 4x, a)') '   rh', 'ext (m^2 kg^-1)', '  ssa', '  asm'

   ! loop through test rh values
   do krh_test = 1, nrh_test
      ! find corresponding rh index
      rh = rh_test(krh_test)
      krh = min(floor( (rh) * nrh ) + 1, nrh - 1)
      wrh = (rh) *nrh - krh
      exti = ext(krh + 1) * (wrh + 1) - ext(krh) * wrh
      ssai = ssa(krh + 1) * (wrh + 1) - ssa(krh) * wrh
      asmi = asm(krh + 1) * (wrh + 1) - asm(krh) * wrh
      write(iulog, '(2x, f5.3, 4x, f13.3, 4x, f5.3, 4x, f5.3)') rh_test(krh_test), exti, ssai, asmi
   end do

end subroutine aer_optics_log_rh


!================================================================================================

end module phys_prop
