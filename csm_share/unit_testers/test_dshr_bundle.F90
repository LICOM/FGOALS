!
! test dshr_bundle
!
program test_dshr_bundle
   use shr_sys_mod
   use shr_kind_mod
   use dshr_kind
   use test_mod
   use dshr_bundle
   use dshr_domain
   use bundle_expected
   use shr_string_mod
   use shr_mpi_mod, only : shr_mpi_init, shr_mpi_finalize, shr_mpi_commRank
   use dshr_rest
   use shr_log_mod
   use dshr_iocdf, only : dshr_iocdf_setDebug, dshr_iocdf_create, dshr_iocdf_append, dshr_iocdf_appendAtt
   use shr_date_mod

implicit none

   logical :: status
   type(dshr_bundle_bundleType) :: bun, Gbun
   type(dshr_bundle_bundleType) :: expected_bun  ! expected bundle
   type(dshr_domain_domainType) :: domain    ! bundle's domain
   type(dshr_domain_domainType) :: Gdomain   ! bundle's global domain
   character(SHR_KIND_CL)       :: name      ! bundle's name
   character(SHR_KIND_CL)       :: dom_name  ! domain name
   character(SHR_KIND_CL)       :: field
   character(SHR_KIND_CX)       :: fieldList, fieldList_get ! list of data fields
   integer :: ni, nj, date, nf, sec, ni_get, nj_get, date_get, nf_get, sec_get
   integer :: i, j, f, test, rcode
   integer :: iarr(3), iarr_get(3)
   real(r8), pointer :: data(:,:), data3D_ptr(:,:,:)
   integer, parameter :: num_tests = 4
   integer :: mpicom, comPID, fid
   integer :: ngi,ngj, ndi,ndj, i0,j0
   character(*), parameter :: filename = "test_bundle_file.nc"
   type(shr_date) :: sdate
#include <mpif.h>

   call shr_mpi_init( )
   mpicom = MPI_COMM_WORLD
   call test_init( num_tests*19 )
   call shr_mpi_commRank(mpicom, comPID )
 
   shr_log_Level = 2
   call dshr_iocdf_setDebug( 4 )
   call dshr_bundle_setDebug( 4 )

   do test = 1, num_tests
      if (      test == 1 )then
         fieldList = "U:V:T:Q:PREC:RAD:O:Z:TS:VAL"
         ni = 4
         nj = 5
      else if ( test == 2 )then
         fieldList = "U:V:T"
         ni = 40
         nj = 50
      else if ( test == 3 )then
         fieldList = "U:T"
         ni = 100
         nj = 300
      else if ( test == 4 )then
         fieldList = "TBOT:WIND:QBOT:PRECTmms:FSDS:PSRF"
         ni = 48
         nj = 96
      else
         call shr_sys_abort( "test number out of bounds" )
      end if
      dom_name = "global_domain_test"
      call dshr_domain_fill(Gdomain,dom_name,ni,nj)
      call dshr_domain_getDims      (Gdomain,ngi,ngj)               ! get global dims
      call dshr_domain_decomp2d(mpicom,comPID, ngi,ngj, ndi,ndj, i0,j0, ni,nj) ! compute 2d decomp
      dom_name = "local_domain_test"
      call dshr_domain_create       (domain,dom_name,ni,nj)         ! create/alloc local domain
      call dshr_domain_putDecompInfo(domain,i0,j0,ngi,ngj,ndi,ndj) ! put decomp info into domain
      call dshr_domain_extractLocal (Gdomain,domain)               ! fill local domain with data


      nf = shr_string_listGetNum( fieldList )
      name = "bundle_test"
      call test_is( .not. dshr_bundle_checkInit(bun),          "Test that bundle is NOT initialized before create" )
      call test_is( .not. dshr_bundle_checkInit(expected_bun), "Test that expected bundle is NOT initialized before create" )
      call dshr_bundle_create(bun,domain,name,fieldList)
      call dshr_bundle_create(expected_bun,domain,name="expected_bundle",fieldList=fieldList)
      call test_is( .not. dshr_bundle_checkDate(bun),          "Test that bundle date is NOT set after create" )
      call test_is( .not. dshr_bundle_checkDate(expected_bun), "Test that expected bundle date is NOT set after create" )
      date = 20080527
      sec  = 3600
      call dshr_bundle_putDate( bun, date, sec )
      call dshr_bundle_putDate( expected_bun, date, sec )
      if ( test == 1 )then
         allocate( data(ni,nj) )
         do f = 1, nf
            call shr_string_listGetName( fieldList, f, field )
            do j = 1, nj
            do i = 1, ni
               data(i,j) = real(i,r8)/1000.0 + real(j,r8) + real(f,r8)*100.0
            end do
            end do
            call dshr_bundle_putField( bun, data, field )
         end do
         deallocate( data )
      else if ( test == 2 )then
          call dshr_bundle_fill( bun, type="index" )
      else if ( test == 3 )then
          call dshr_bundle_fill( bun, type="sincos" )
      else if ( test == 4 )then
          call dshr_bundle_fill( bun, type="nfield" )
      end if
      if ( test /= 4 )then
         call dshr_bundle_copyFields( bun, expected_bun )
      else
         allocate( data(ni,nj) )
         do f = 1, nf
            call shr_string_listGetName( fieldList, f, field )
            do j = 1, nj
            do i = 1, ni
               data(i,j) = real(f,r8)
            end do
            end do
            call dshr_bundle_putField( expected_bun, data, field )
         end do
         deallocate( data )
      end if
      call dshr_bundle_print( bun )
      call dshr_bundle_print( expected_bun )
      status = bundle_is_expected( bun, expected_bun )
      call test_is( status, "Test if bundle is same as expected" )
      status = bundle_closeto_expected( bun, expected_bun, eps=0.0_r8 )
      call test_is( status, "Test if bundle is close to expected within 0.0" )
      call test_is( dshr_bundle_checkInit(bun),          "Test if bundle is initialized" )
      call test_is( dshr_bundle_checkInit(expected_bun), "Test if bundle is initialized" )
      call test_is( dshr_bundle_checkDate(bun),          "Test if bundle date is set" )
      call test_is( dshr_bundle_checkDate(expected_bun), "Test if bundle expected date is set" )
      call dshr_bundle_getNumFields( bun, nf_get )
      call test_is( nf_get, nf, "Test that bundles have expected number of fields" )
      call dshr_bundle_getFieldList( bun, fieldList_get)
      call test_is( nf_get, nf, "Test that bundles get expected field list" )
      call test_is( fieldList_Get, fieldList, "Test that bundles have expected field list" )
      call dshr_bundle_getDims( bun, ni_get, nj_get, nf_get )
      iarr_get(1) = ni_get
      iarr_get(2) = nj_get
      iarr_get(3) = nf_get
      iarr(1)     = ni
      iarr(2)     = nj
      iarr(3)     = nf
      call test_is( iarr_get, iarr, "Test that bundles have expected dimensions" )
   
      call dshr_bundle_putDate( expected_bun, cdate=19640527, sec=0 )
      status = bundle_is_expected( bun, expected_bun )
      call test_is( .not. status, "Test that bundles are NOT same if date changes" )
      call dshr_bundle_putDate( expected_bun, cdate=date, sec=7200 )
      status = bundle_is_expected( bun, expected_bun )
      call test_is( .not. status, "Test that bundles are NOT same if sec changes" )
      call dshr_bundle_assignPtr( expected_bun, data3D_ptr )
      data3D_ptr = data3D_ptr(:,:,:) * 3.0_SHR_KIND_R4
      status = bundle_is_expected( bun, expected_bun )
      call dshr_bundle_putDate( expected_bun, cdate=date, sec=sec )
      call test_is( .not. status, "Test that bundles are NOT same if data changes" )
      data3D_ptr = data3D_ptr(:,:,:) / 3.00000001_SHR_KIND_R4
      status = bundle_closeto_expected( bun, expected_bun, eps=1.e-3_r8 )
      call test_is( status, "Test that bundles are close to each other if data off by a bit" )
      call dshr_bundle_print( bun )
      call dshr_bundle_print( expected_bun )

      call dshr_bundle_copyFields( bun, expected_bun )

      !
      ! Write restart file out, clean expected bundle and read into it, then make sure bundles are the same
      !
      write(*,*) "Gather local bundle to global one, then write it out"
      call dshr_bundle_create(Gbun,domain,name="global_bundle",fieldList=fieldList)
      call dshr_bundle_gather(Gbun,bun,mpicom,"test_bundle")
      call dshr_bundle_print( Gbun )
      if (comPID==0) then
         write(*,*) "Write file out on masterproc"
         sdate = shr_date_initCDate( cdate=date, ns=24 )
         call shr_sys_system( "/bin/rm -f "//trim(filename), rcode )
         call dshr_iocdf_create(filename)
         call dshr_iocdf_appendAtt(fid,"test_bundle_restart_fieldList",fieldList,filename)
         call dshr_iocdf_append(fid,sdate,Gbun,filename)
      end if
      call dshr_bundle_clean( expected_bun )
      call dshr_bundle_clean( Gbun )
      call dshr_bundle_create(expected_bun,domain,name="expected_bundle",fieldList=fieldList)
      call dshr_bundle_create(Gbun,domain,name="global_bundle",fieldList=fieldList)
      call dshr_rest_readBundle(filename,Gbun)
      call dshr_bundle_putDate( Gbun, cdate=date, sec=sec )
      call dshr_bundle_putDate( expected_bun, cdate=date, sec=sec )
      call dshr_bundle_scatter(Gbun,expected_bun, mpicom)
      status = bundle_is_expected( bun, expected_bun )
      call test_is( status, "Test that bundles written out and then read in are identical" )
   
      call dshr_bundle_clean( bun )
      call dshr_bundle_clean( expected_bun )
      call dshr_domain_clean(domain)

   end do

   call shr_mpi_finalize( )
   call test_final( status )
   if ( status ) call shr_sys_system( "/bin/rm -f "//trim(filename), rcode )

end program test_dshr_bundle
