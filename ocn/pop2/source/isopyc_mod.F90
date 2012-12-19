!  CVS: $Id: isopyc_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module isopyc_mod
#include <def-undef.h>
use precision_mod
use param_mod
!     isopycnal diffusion variables:
!
!     ahisop = isopycnal tracer diffusivity (cm**2/sec) 
!     athkdf = isopycnal thickness diffusivity (cm**2/sec)
!     dptlim = depth limits for the reference pressure levels (in cm). 
!              the mid point of the two neighboring "dptlim" elements is
!              used as the reference pressure (depth) for that interval,
!              i.e.,
!
!             reference pressure level        reference depth
!              ----------------------- -----------------------------
!                          1           0.5*(dptlim(1)+dptlim(2))
!                          2           0.5*(dptlim(2)+dptlim(3))
!                          .           .
!                        nrpl          0.5*(dptlim(nrpl)+dptlim(nrpl+1))
!
!              REMARK: the first and the last elements of "dptlim" must
!                      be the depth at the top (0m) and the maximum bottom
!                      depth, respectively. also, the elements of "dptlim"
!                      must be in increasing order.
!     e      = scratch array 
!     K1     = (1,3) component of the isopycnal mixing tensor computed
!              at the center of the eastern face of the "t" grid cell
!     K2     = (2,3) component of the isopycnal mixing tensor computed
!              at the center of the northern face of the "t" grid cell
!     K3     = (3,.) components of the isopycnal mixing tensor 
!              computed at the top face of the "t" grid cell
!               (,,1) --> (3,1) component
!               (,,2) --> (3,2) component
!               (,,3) --> (3,3) component 
!              REMARK: (1,1) and (2,2) components of the tensor are
!                      assumed to be unity. also, (1,2) and (2,1)
!                      components are assumed to be zero.
!     fzisop = function containing the vertical variation of the isopycnal
!              diffusion coefficient. "fzisop" multiplies "ahisop".
!     kisrpl = isopycnal reference pressure levels for the "t" grid 
!              point levels computed based on the depth (pressure) at the
!              "t" grid points
!     krplin = indices indicating the location of the reference pressure 
!              depth in the 20-level table of polynomial expansion
!              variables
!     slmxr  = reciprocal of maximum slope of isopycnals allowed in mixing
!              scheme to prevent excessively large vertical mixing that
!              could create numerical instabilities. furthermore, the
!              form of the isopycnal diffusion tensor incorporates the
!              assumption that horizontal density gradients are much
!              smaller than vertical gradients. "slmxr" is also used to
!              satisfy this assumption. a value of 100 for "slmxr"
!              translates to a slope of 1:100.
!
!     rhoi    = potential density at "t" cell centers
!
!# ifdef isopycmixspatialvar
!     dciso1 = isopycnal tracer diffusivity coefficients modified based
!              on the slopes of the isopycnal surfaces on the east face
!              of "T" cells.
!     dciso2 = isopycnal tracer diffusivity coefficients modified based
!              on the slopes of the isopycnal surfaces on the north face
!              of "T" cells.
!     dslope = half length of the interval in which "ahisop" changes
!              with a steep slope from about 0.9*"ahisop" to about
!              0.1*"ahisop"
!     slopec = slope at which "ahisop" is equal to half of its original
!              value
!
!              REMARK: 0 <= "slopec", "dslope" <= 1/"slmxr".
!              REMARK: because the vertical gradient of density must be
!                      less than zero for isopycnal mixing, 1/"slmxr" is
!                      actually a negative maximum slope. this fact is
!                      taken into account in "isopyc.F". consequently,
!                      "slmxr" must be a positive number in "blkdta.F".
!# endif
!# ifdef gent_mcwilliams
!     adv_vetiso = zonal isopycnal mixing velocity computed at the 
!                  center of the eastern face of the "t" cells
!     adv_vntiso = meridional isopycnal mixing velocity computed at
!                  the center of the northern face of the "t" cells
!                  (Note: this includes the cosine as in "adv_vnt")
!     adv_vbtiso = vertical isopycnal mixing velocity computed at the
!                  center of the top face of the "t" cells
!     adv_fbiso  = "adv_vbtiso" * (tracer) evaluated at the center of
!                  the bottom face of the "t" cells
!# endif
!
!     nrpl = number of reference presure levels used in isopycnal
!
#if (defined ISO)
      integer,parameter:: nrpl=5
!
!YU
!YU
!

      real(r8):: ahisop, athkdf, dptlim(nrpl+1), fzisop(km), slmxr
      real(r8),dimension(:,:,:,:),allocatable:: e,rhoi
      real(r8),dimension(:,:,:,:),allocatable:: K1,K2,K3
      real(r8),dimension(:,:,:),allocatable:: adv_vetiso,adv_vntiso,adv_vbtiso
#ifdef isopycmixspatialvar
      real(r4) :: dciso1(imt,km,jmt), dciso2(imt,km,jmt) &
              , dslope, slopec
#endif
!
      real(r8):: kisrpl(km) 
      integer:: krplin(nrpl)
      real(r8)::  zt(km),dzw(0:km),dzwr(0:km),dzr(km),dytr(jmt),&
                      cstrdytr(jmt),dyur(jmt),tmask(imt,km,jmt)
      REAL(r8) , dimension(JMT) :: F3
#endif
end module isopyc_mod
