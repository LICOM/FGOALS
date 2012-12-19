package Build::ChemPreprocess;
#-------------------------------------------------------------------------------------
# ($chem_nadv,$chem_cppdefs) = chem_preprocess( $cfg_ref, $prnt_lvl )
#
# This routine does the following:
# - Invokes the chemistry preprocessor
# - Checks consistancy of configure options
# - Determines the number of transported chemstical tracers
# - Sets the chemstry CPP definitions
#
# Date         Contributor      Modification
#-------------------------------------------------------------------------------------
# 19 Sep 2008  Francis Vitt     Created 
# 23 Oct 2009  Francis Vitt     moved to perl5lib/Build directory
#                               renamed to ChemPreprocess.pm and implemented as a module
#-------------------------------------------------------------------------------------

use strict;
use Exporter;
use File::Copy;

our @ISA = qw(Exporter);
our @EXPORT = qw(chem_preprocess get_species_list);
our $VERSION = 1.00;

my $print ;

sub chem_preprocess 
{
    my ($cfg_ref,$compiler,$prnt_lvl) = @_;

    $print = $prnt_lvl;

    my $chem_nadv = 0;
    my $chem_cppdefs = '';

    my $edit_chem_mech = $cfg_ref->get('edit_chem_mech');
    my $usr_mech_infile = $cfg_ref->get('usr_mech_infile');
    my $prog_species = $cfg_ref->get('prog_species');
    my $chem_pkg = $cfg_ref->get('chem');
    my $cam_root = $cfg_ref->get('cam_root');
    my $cam_bld = $cfg_ref->get('cam_bld');

    my $chem_proc_bld = "$cam_bld/chem_proc";
    my $chem_src_dir  = "$chem_proc_bld/source";
    $cfg_ref->set('chem_src_dir', $chem_src_dir);

    my $chem_preprocessor = "$cam_root/models/atm/cam/chem_proc";

    my $chem_mech_infile;

    if ($print>=2){ print "chem_preprocess:    prog_species = $prog_species \n"; }
    if ($print>=2){ print "chem_preprocess:        chem_pkg = $chem_pkg \n"; }
    if ($print>=2){ print "chem_preprocess: usr_mech_infile = $usr_mech_infile \n"; }
    if ($print>=2){ print "chem_preprocess:  edit_chem_mech = $edit_chem_mech \n"; }
    if ($print>=2){ print "chem_preprocess:        compiler = $compiler \n"; }

    if ($prog_species) {
	if ($chem_pkg =~ /"mozart"/) {
	    die "ERROR: -prog_species $prog_species is NOT allowed with -chem $chem_pkg \n";
	}
	if ($usr_mech_infile) {
	    die "ERROR: -prog_species $prog_species is NOT allowed with -usr_mech_infile $usr_mech_infile \n";
	}
    }
    if (!$chem_pkg) {
	if ($usr_mech_infile) {
	    die "ERROR: -usr_mech_infile $usr_mech_infile is NOT allowed without -chem option. \n";
	}
    }

    # create chem proc directory tree
    my $cmd = "mkdir -p $chem_proc_bld/tmp";
    run_shell_command($cmd);

    if (!$usr_mech_infile) {
	if ($prog_species) {
	    if ($usr_mech_infile) {
		die "ERROR: *** Cannot specify usr_mech_infile with prog_species  *** \n" ;
	    }
	    $chem_mech_infile = "$chem_proc_bld/chem_mech.inp";
	    write_chem_preproc($chem_mech_infile, $cfg_ref, $chem_preprocessor , $chem_proc_bld);
	} else {
	    $usr_mech_infile = "$cam_root/models/atm/cam/src/chemistry/pp_${chem_pkg}/chem_mech.in";
	}
    }

    if ($usr_mech_infile) {
	if ($prog_species) {
	    die "ERROR: *** Cannot specify usr_mech_infile with prog_species *** \n" ;
	}
	$chem_mech_infile = "$chem_proc_bld/chem_mech.inp";
	write_chem_mech($usr_mech_infile, $chem_mech_infile, $chem_preprocessor , $chem_proc_bld );
    }

    if ($chem_mech_infile) {

	if ($edit_chem_mech) {
	    edit_chem_preproc($chem_mech_infile);
	}

	$cfg_ref->set('chem_proc_bld', $chem_proc_bld);

	my $chem_proc_src = "$cam_bld/chem_proc/source";
	$cfg_ref->set('chem_proc_src',  $chem_proc_src  );

	my $chem_proc_exe = "campp";

	if (! -e "$chem_proc_bld/$chem_proc_exe") {
	    my $gmake = 'gmake';
	    build_chem_preproc($gmake,$chem_preprocessor ,$chem_proc_bld,$chem_proc_exe,$compiler);
	}


	run_chem_preproc($chem_proc_bld,$chem_proc_exe,$chem_mech_infile,$chem_proc_src,$cam_bld);

	if ( $chem_pkg =~ /trop/ or $prog_species ) {
	    $chem_cppdefs .= ' -DTROPCHEM';
	}
	if ( $prog_species =~ /DST/ || $prog_species =~ /SSLT/ ) {
	    $chem_cppdefs .= ' -DDUST';
	}
	if ( $prog_species =~ /SSLT/ ) {
	    $chem_cppdefs .= ' -DPROGSSLT';
	}

	# determine the number of transported chemical tracers
	open INPUT, "$chem_proc_src/chem_mods.F90";
	while ( my $line = <INPUT> ) {
	    if ( $line =~ /gas_pcnst\s*=/ ) {
		if($line =~ m/(\d+)/) { # extract the number of chem species
		    $chem_nadv += $1;
		    if ($print>=2) { print "total number of chemical species = $chem_nadv  \n"; }
		} else { 
		    die "**** Not able to determine total number of chemical species ****\n";
		}
	    }
	    if ( $line =~ /nslvd\s*=/ ) {
		if($line =~ m/(\d+)/) { # extract the number of chem species
		    $chem_nadv = $chem_nadv - $1;
		    if ($print>=2) { print "number of short-lived chemical species = $1  \n"; }
		    if ($print>=2) { print "number of transported chemical species = $chem_nadv  \n"; }
		} else { 
		    die "**** Not able to determine number of short-lived species ****\n";
		}
	    }
	}
	close INPUT;

        my $dst_cppdefs = '';
        my $sslt_cppdefs = '';
	my $chem_nwat = 0;


        my @species = get_species_list($chem_src_dir);
        foreach my $tracer (@species) {
	  if ( $tracer eq 'H2O' ) {
	      $chem_nwat = 1;
	  }
	  if ( $tracer eq 'DST01' ) {
	      $dst_cppdefs = ' -DDUST';
	  }
	  if ( $tracer eq 'SSLT01' ) {
	      $sslt_cppdefs = ' -DDUST -DPROGSSLT';
	  }
        }

	if ($print>=2) { print "number of water vapor species = $chem_nwat  \n"; }
	$chem_nadv -= $chem_nwat ;
        $chem_cppdefs .= $dst_cppdefs;
        $chem_cppdefs .= $sslt_cppdefs;

    }

    if ($print>=2) { print "Number of chem adv tracers: $chem_nadv \n"; }
    if ($print>=2) { print "Chemistry CCP defs: $chem_cppdefs \n"; }

    return ($chem_nadv,$chem_cppdefs);
}

#-----------------------------------------------------------------------------------------------
# Utility routines
#-----------------------------------------------------------------------------------------------

sub get_species_list
{
    my ($chem_src_dir) = @_;

    my @species_list ;

    my $end_of_rec = $/;
    $/ = "/)";

    open INPUT, "$chem_src_dir/mo_sim_dat.F90";

    while ( my $data = <INPUT> ) {
	if ( $data =~ /\s*solsym\(:/  ) {
	    chomp $data ;
 
	    my @list = split( /\//, $data );
	    my @spec_list = split( /\W+/, @list[ $#list ] );
	    foreach my $item (@spec_list) {
		if ( length($item) > 0 ){
		    push ( @species_list, $item );
		}
	    }

	}
    }
    close INPUT;
    $/ = $end_of_rec;

    return ( @species_list );
}

sub run_shell_command {

    my ($cmd) = @_;

    if ($print>=2) { print "cmd = $cmd\n";}

    my @out = `$cmd`;
    my $cmd_error = $? ; #CHILD_ERROR;
    foreach my $i (@out) {
        if ($print>=2) { print "$i";}
        if ($cmd_error || $i =~ /abort/ || $i =~ /Failed/ ) {
            die "**** FAILED ****\n$i\n";
        }
    }

}

#-----------------------------------------------------------------------------------------------

sub edit_chem_preproc
{
    my ($chem_proc_inp) = @_;
    my $cam_chem_editor = 'vi';
    
    if ($print>=2) { print "edit chemistry mechanism file.... \n";}

    if (defined $ENV{CAMCHEM_EDITOR}) {
	$cam_chem_editor = $ENV{CAMCHEM_EDITOR};
    }

    my $command = "$cam_chem_editor $chem_proc_inp";

    my $status = system("$command");
    if (($status >>=8) != 0) {
        die "Failed to run $command";
    }
    
    if ($print>=2) { print "edit chemistry mechanism file complete. \n";}
}

#-----------------------------------------------------------------------------------------------

sub run_chem_preproc
{
    my ($chem_proc_bld,$chem_proc_exe,$chem_proc_inp,$src_dir,$cam_bld) = @_;

    if ($print>=2) { print "run_chem_preproc.... \n";}

    # clean out old version
    my $cmd = "rm -rf $src_dir $chem_proc_bld/cam.subs.tar";
    run_shell_command($cmd);

    # run chem preprocessor
    my $cmd = "$chem_proc_bld/$chem_proc_exe $chem_proc_inp 2>&1";
    run_shell_command($cmd);
    
    if ($print) { print "creating $src_dir\n"; }

    # create dir to for new code
    my $cmd = "mkdir -p $src_dir";
    run_shell_command($cmd);

    # get new code from tar filex
    #my $cmd = "tar -xf $chem_proc_bld/cam.subs.tar -C $src_dir";
    my $cmd = "cd $src_dir && tar -xf $chem_proc_bld/cam.subs.tar";
    run_shell_command($cmd);

    copy( $chem_proc_inp        ,"$cam_bld/chem_mech.inp") or die "copy failed $! \n";
    copy( "$chem_proc_bld/chem_mech.doc" ,$cam_bld) or die "copy failed $! \n";

    if ($print>=2) { print "run_chem_preproc complete\n"; }
}

#-----------------------------------------------------------------------------------------------

sub build_chem_preproc
{
    my ($gmake,$chem_proc_src,$chem_proc_bld,$chem_proc_exe,$compiler) = @_;

    if ($print>=2) { print "build_chem_preproc.... \n"; }
    if ($print) { print "creating $chem_proc_bld/$chem_proc_exe\n"; }

    $ENV{'MODEL_EXEDIR'} = "$chem_proc_bld";
    $ENV{'EXENAME'} = "$chem_proc_exe";
    $ENV{'SRCLIST'} = "$chem_proc_src/src/Base_Srclist_f";
    $ENV{'SRCDIRS'} = "$chem_proc_src/src/cam_chempp";
    $ENV{'OBJ_DIR'} = "$chem_proc_bld";
    if ( $compiler =~ /lf95/ ) {
      $ENV{'COMPILER'} = "Lahey";
    } else {
      $ENV{'COMPILER'} = $compiler;
    }

    my $log_file = "$chem_proc_bld/MAKE.out";
    my $makefile = "$chem_proc_src/src/Makefile";
    my $cmd = "$gmake -f $makefile 2>&1";
    run_shell_command($cmd);

    my $cmd = "rm -f $chem_proc_bld/*.o $chem_proc_bld/*.mod";
    run_shell_command($cmd);

    if ($print>=2) { print "build_chem_preproc complete\n"; }
}

#-----------------------------------------------------------------------------------------------

sub write_chem_mech
{

    my ($file_in, $chem_proc_file, $proc_src, $proc_bld) = @_;

    my $fh_in = new IO::File;
    my $fh_out = new IO::File;

    if ($print) { print "creating $chem_proc_file\n"; }

    $fh_out->open(">$chem_proc_file") or die "** can't open chem preprocessor input file: $chem_proc_file\n";

print $fh_out <<"EOF";

BEGSIM
output_unit_number = 7
output_file        = chem_mech.doc
temp_path          = $proc_bld/tmp/
procout_path       = $proc_bld/
output_path        = $proc_bld/
src_path           = $proc_src/bkend/
procfiles_path     = $proc_src/procfiles/cam/
sim_dat_path       = $proc_bld/
sim_dat_filename   = chem_mech.dat

EOF

    # Copy the chemistry mechanism.
    $fh_in->open("<$file_in") or die "** can't open file: $file_in\n";
    while (<$fh_in>) {
	print $fh_out $_;
    }
    $fh_in->close;

print $fh_out <<"EOF";

ENDSIM

EOF

    $fh_out->close;
}

#-----------------------------------------------------------------------------------------------

sub write_chem_preproc
{

    my ($chem_proc_file, $cfg_ref, $proc_src, $proc_bld) = @_;

    my $prog_species = $cfg_ref->get('prog_species');
    my $fh = new IO::File;

    if ($print) { print "creating $chem_proc_file\n"; }

    $fh->open(">$chem_proc_file") or die "** can't open chem preprocessor input file: $chem_proc_file\n";

print $fh <<"EOF";

BEGSIM
output_unit_number = 7
output_file        = chem_mech.doc
temp_path          = $proc_bld/tmp/
procout_path       = $proc_bld/
output_path        = $proc_bld/
src_path           = $proc_src/bkend/
procfiles_path     = $proc_src/procfiles/cam/
sim_dat_path       = $proc_bld/
sim_dat_filename   = chem_mech.dat

Comments
     "This is a CAM simulation with :"
End Comments

      SPECIES

      Solution
EOF

    if ( $prog_species =~ /SO4/ ) {
	print $fh "    H2O2, SO2, SO4, DMS -> CH3SCH3\n";
    }
    if ( $prog_species =~ /OC/ ) {
	print $fh "    OC1 -> C, OC2 -> C\n";
    }
    if ( $prog_species =~ /BC/ ) {
        print $fh "    CB1 -> C, CB2 -> C\n";
    }
    if ( $prog_species =~ /GHG/ ) {
        print $fh "    CH4, N2O, CFC11 -> CFCl3, CFC12 -> CF2Cl2, H2O\n";
    }
    if ( $prog_species =~ /SSLT/ ) {
        print $fh "    SSLT01 -> NaCl, SSLT02 -> NaCl, SSLT03 -> NaCl, SSLT04 -> NaCl\n";
    }
    if ( $prog_species =~ /DST/ ) {
        print $fh "    DST01 -> AlSiO5, DST02 -> AlSiO5, DST03 -> AlSiO5, DST04 -> AlSiO5\n";
    }
    if ( $prog_species =~ /CARBON16/ ) {
	print $fh "    OFPHO -> C, BFPHO -> C, OBPHO -> C, BBPHO -> C \n";
	print $fh "    OOPHO -> C, BOPHO -> C, NOPHO -> C, MMPHO -> C \n";
	print $fh "    OFPHI -> C, BFPHI -> C, OBPHI -> C, BBPHI -> C \n";
	print $fh "    OOPHI -> C, BOPHI -> C, NOPHI -> C, MMPHI -> C \n";
    }

print $fh  <<"EOF";
      End Solution

      Fixed
    M, N2, O2
EOF
    if ( $prog_species =~ /SO4/ ) {
	print $fh "    O3, OH, NO3, HO2\n"; 
    } 
print $fh <<"EOF";
      End Fixed

      Col-int
 O3 = 0.
 O2 = 0.
      End Col-int

   End SPECIES

   Solution Classes
      Explicit
      End Explicit
      Implicit
EOF
    if ( $prog_species =~ /SO4/ ) {
	print $fh "    H2O2, SO2, SO4, DMS\n"; 
    }
    if ( $prog_species =~ /CARBON16/ ) {
	print $fh "    OFPHO,BFPHO,OBPHO,BBPHO,OOPHO,BOPHO,NOPHO,MMPHO \n";
	print $fh "    OFPHI,BFPHI,OBPHI,BBPHI,OOPHI,BOPHI,NOPHI,MMPHI \n";
    }
    if ( $prog_species =~ /BC/ ) {
	print $fh "    CB1, CB2\n";
    }
    if ( $prog_species =~ /OC/ ) {
        print $fh "    OC1, OC2\n";
    }
    if ( $prog_species =~ /GHG/ ) {
        print $fh "    CH4, N2O, CFC11, CFC12, H2O\n";
    }
    if ( $prog_species =~ /SSLT/ ) {
        print $fh "    SSLT01, SSLT02, SSLT03, SSLT04\n";
    }
    if ( $prog_species =~ /DST/ ) {
        print $fh "    DST01, DST02, DST03, DST04\n";
    }
print $fh <<"EOF";
      End Implicit
   End Solution Classes

 CHEMISTRY
      Photolysis
      End Photolysis

      Reactions
EOF
    if ( $prog_species =~ /CARBON16/ ) {
	print $fh "    OFPHO -> OFPHI    ; 1.006e-05 \n";
	print $fh "    BFPHO -> BFPHI    ; 1.006e-05 \n";
	print $fh "    OBPHO -> OBPHI    ; 1.006e-05 \n";
	print $fh "    BBPHO -> BBPHI    ; 1.006e-05 \n";
	print $fh "    OOPHO -> OOPHI    ; 1.006e-05 \n";
	print $fh "    BOPHO -> BOPHI    ; 1.006e-05 \n";
	print $fh "    NOPHO -> NOPHI    ; 1.006e-05 \n";
	print $fh "    MMPHO -> MMPHI    ; 1.006e-05 \n";
    }
    if ( $prog_species =~ /BC/ ) {
	print $fh "    CB1 -> CB2    ; 1.006e-05 \n";
    }
    if ( $prog_species =~ /OC/ ) {
	print $fh "    OC1 -> OC2    ; 1.006e-05 \n";
    }
    if ( $prog_species =~ /GHG/ ) {
	print $fh "    [ch4_loss]    CH4   -> 2.* H2O\n";
	print $fh "    [n2o_loss]    N2O   -> \n";
	print $fh "    [cfc11_loss]  CFC11 -> \n";
	print $fh "    [cfc12_loss]  CFC12 -> \n";
	print $fh "    [lyman_alpha] H2O   -> \n";
    }
print $fh <<"EOF";
      End Reactions

      Heterogeneous
EOF
    if ( $prog_species =~ /SO4/ ) {
        print $fh "    H2O2, SO2 \n";
    }
print $fh <<"EOF";
      End Heterogeneous

      Ext Forcing
EOF
    if ( $prog_species =~ /SO4/ ) {
        print $fh "    SO2 <- dataset\n";
        print $fh "    SO4 <- dataset\n";
    }
print $fh <<"EOF";
      End Ext Forcing

   END CHEMISTRY

   SIMULATION PARAMETERS

     Version Options
        model   = cam
        machine = intel
        architecture = hybrid
        vec_ftns  = on
        multitask = on
        namemod = on
        modules = on
     End Version Options

   END SIMULATION PARAMETERS

ENDSIM

EOF

    $fh->close;
}


1; # to appease require 
