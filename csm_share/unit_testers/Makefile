#-----------------------------------------------------------------------
# This Makefile is for doing csm_share unit testing
#------------------------------------------------------------------------
cpp_dirs := . ../shr ../dshr ../../utils/esmf_wrf_timemgr
ifneq ($(SPMD),TRUE)
cpp_dirs += ../../utils/mct/mpi-serial
endif
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Platform specific macros
include make.Macros

space := $(null) $(null)

ifneq ($(ESMF_BLD),$(null))
cpp_dirs += $(ESMF_LIB)
endif

# Expand any tildes in directory names.
VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
# Change spaces to colons.
VPATH    := $(subst $(space),:,$(VPATH))

all: test_shr_tInterp


ifneq ($(SPMD),TRUE)
   OBJS_NOMPI   := fort.o group.o collective.o comm.o list.o handles.o mpi.o recv.o req.o \
	      		send.o time.o
else
   OBJS_NOMPI   := $(null)
endif

OBJS       := test_shr_sys.o shr_sys_mod.o shr_kind_mod.o shr_mpi_mod.o shr_const_mod.o shr_log_mod.o \
	      $(OBJS_NOMPI)
OBJS_FILE  := test_shr_file.o shr_sys_mod.o shr_kind_mod.o shr_file_mod.o shr_mpi_mod.o shr_log_mod.o \
	      $(OBJS_NOMPI)
OBJS_IN    := test_shr_inputinfo.o shr_sys_mod.o shr_kind_mod.o shr_mpi_mod.o shr_const_mod.o shr_log_mod.o \
		shr_cal_mod.o shr_file_mod.o shr_string_mod.o shr_timer_mod.o \
		shr_inputinfo_mod.o shr_ncio_mod.o $(OBJS_NOMPI)
OBJS_ORB := test_shr_orb.o shr_sys_mod.o shr_kind_mod.o shr_orb_mod.o shr_mpi_mod.o shr_log_mod.o \
	    shr_const_mod.o $(OBJS_NOMPI)
OBJS_DATE := test_shr_date.o shr_date_mod.o shr_cal_mod.o shr_kind_mod.o shr_sys_mod.o shr_log_mod.o \
               shr_const_mod.o shr_mpi_mod.o test_mod.o  $(OBJS_NOMPI)
OBJS_STRMS  := test_shr_streams.o shr_kind_mod.o shr_stream_mod.o shr_sys_mod.o \
               shr_file_mod.o shr_string_mod.o shr_timer_mod.o shr_mpi_mod.o    \
               shr_cal_mod.o shr_date_mod.o shr_ncread_mod.o shr_const_mod.o    \
               shr_log_mod.o test_mod.o $(OBJS_NOMPI)
OBJS_DSBUNT := test_dshr_tInterp.o shr_kind_mod.o shr_stream_mod.o shr_sys_mod.o \
               shr_file_mod.o shr_string_mod.o shr_timer_mod.o shr_mpi_mod.o    \
               shr_cal_mod.o shr_date_mod.o shr_ncread_mod.o shr_const_mod.o    \
               shr_log_mod.o test_mod.o dshr_tInterp.o dshr_kind.o dshr_domain.o \
               dshr_const.o bundle_expected.o dshr_bundle.o shr_tInterp_mod.o    \
               shr_orb_mod.o $(OBJS_NOMPI)
OBJS_DSBUN  := test_dshr_bundle.o shr_kind_mod.o shr_stream_mod.o shr_sys_mod.o \
               shr_file_mod.o shr_string_mod.o shr_timer_mod.o shr_mpi_mod.o    \
               shr_cal_mod.o shr_date_mod.o shr_ncread_mod.o shr_const_mod.o    \
               shr_log_mod.o test_mod.o dshr_kind.o dshr_domain.o \
               dshr_const.o bundle_expected.o dshr_bundle.o dshr_rest.o dshr_iocdf.o $(OBJS_NOMPI)
OBJS_STIN   := test_shr_tInterp.o shr_kind_mod.o shr_const_mod.o shr_sys_mod.o \
               shr_cal_mod.o shr_date_mod.o shr_log_mod.o test_mod.o \
               shr_tInterp_mod.o shr_mpi_mod.o $(OBJS_NOMPI)
OBJS_MPI := test_shr_mpi.o shr_mpi_mod.o shr_kind_mod.o shr_sys_mod.o shr_const_mod.o shr_log_mod.o $(OBJS_NOMPI)

OBJS_NAN := test_shr_infnan.o shr_kind_mod.o shr_infnan_mod.o test_mod.o shr_sys_mod.o \
                shr_mpi_mod.o shr_log_mod.o shr_isnan.o $(OBJS_NOMPI)

WRFESMF_OBJS := ESMF_AlarmMod.o ESMF_BaseTimeMod.o ESMF_ClockMod.o ESMF_TimeMod.o \
                Meat.o ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_FractionMod.o   \
                ESMF_TimeIntervalMod.o ESMF_Stubs.o ESMF_Mod.o \
                ESMF_AlarmClockMod.o wrf_error_fatal.o wrf_message.o

ifeq ($(ESMF_BLD),$(null))
   OBJS_TIME   += $(WRFESMF_OBJS)
   OBJS_ESTATE += $(WRFESMF_OBJS)
endif

#
# Executables:
#
debug:
	@echo "VPATH:      " $(VPATH)
	@echo "ESMF_MOD:   " $(ESMF_MOD)
	@echo "ESMF_ARCH:  " $(ESMF_ARCH)
	@echo "FC:         " $(FC)
	@echo "INC_NETCDF: " $(INC_NETCDF)
	@echo "LIB_MPI:    " $(LIB_MPI)
test_shr_sys: $(OBJS)
	$(LD) -o test_shr_sys     $(OBJS) $(LDFLAGS)
test_shr_infnan: $(OBJS_NAN)
	$(LD) -o test_shr_infnan  $(OBJS_NAN) $(LDFLAGS)
test_shr_file: $(OBJS_FILE)
	$(LD) -o test_shr_file    $(OBJS_FILE) $(LDFLAGS)
test_shr_orb: $(OBJS_ORB)
	$(LD) -o test_shr_orb     $(OBJS_ORB) $(LDFLAGS)
test_shr_inputinfo: $(OBJS_IN)
	$(LD) -o test_shr_inputinfo   $(OBJS_IN) $(LDFLAGS)
test_shr_streams: $(OBJS_STRMS)
	$(LD) -o test_shr_streams $(OBJS_STRMS) $(LDFLAGS)
test_shr_date: $(OBJS_DATE)
	$(LD) -o test_shr_date    $(OBJS_DATE)  $(LDFLAGS)
test_shr_tInterp: $(OBJS_STIN)
	$(LD) -o test_shr_tInterp $(OBJS_STIN) $(LDFLAGS)
test_dshr_tInterp: $(OBJS_DSBUNT)
	$(LD) -o test_dshr_tInterp $(OBJS_DSBUNT) $(LDFLAGS)
test_dshr_bundle: $(OBJS_DSBUN)
	$(LD) -o test_dshr_bundle $(OBJS_DSBUN) $(LDFLAGS)
test_shr_mpi: $(OBJS_MPI)
	$(LD) -o test_shr_mpi $(OBJS_MPI) $(LDFLAGS)

clean:
	$(RM) *.mod *.o *.f *.f90 test_shr_inputinfo test_shr_sys test_shr_orb \
		libesmf.a

#
# Dependencies
#

# Test programs
test_shr_file.o: test_shr_file.F90 shr_file_mod.o shr_sys_mod.o shr_kind_mod.o
test_shr_date.o: test_shr_date.F90 shr_date_mod.o shr_cal_mod.o shr_kind_mod.o shr_sys_mod.o shr_log_mod.o \
               shr_const_mod.o shr_mpi_mod.o test_mod.o
test_shr_inputinfo.o: shr_inputinfo_mod.o test_shr_inputinfo.F90
test_shr_sys.o: shr_sys_mod.o shr_sys_mod.o test_shr_sys.F90
test_shr_infnan.o: shr_kind_mod.o shr_infnan_mod.o test_shr_infnan.F90 test_mod.o shr_sys_mod.o
test_shr_orb.o: shr_orb_mod.o test_shr_orb.F90
test_shr_mpi.o: test_shr_mpi.F90 shr_mpi_mod.o shr_sys_mod.o shr_kind_mod.o shr_const_mod.o
test_shr_streams.o: test_shr_streams.F90 shr_kind_mod.o shr_stream_mod.o shr_sys_mod.o \
               shr_file_mod.o shr_string_mod.o shr_timer_mod.o shr_log_mod.o           \
               shr_cal_mod.o shr_date_mod.o shr_mpi_mod.o shr_ncread_mod.o             \
               shr_const_mod.o test_mod.o

test_shr_tInterp.o: test_shr_tInterp.F90 shr_const_mod.o shr_sys_mod.o test_mod.o \
                    shr_kind_mod.o shr_cal_mod.o shr_date_mod.o shr_tInterp_mod.o
test_dshr_bundle.o: test_dshr_bundle.F90 shr_const_mod.o shr_sys_mod.o test_mod.o \
                    shr_kind_mod.o shr_cal_mod.o shr_date_mod.o dshr_rest.o \
                    shr_string_mod.o shr_mpi_mod.o shr_log_mod.o shr_ncread_mod.o \
                    shr_file_mod.o dshr_domain.o dshr_const.o    dshr_iocdf.o   \
                    dshr_kind.o dshr_bundle.o bundle_expected.o
test_dshr_tInterp.o: test_dshr_tInterp.F90 shr_const_mod.o shr_sys_mod.o test_mod.o \
                    shr_kind_mod.o shr_cal_mod.o shr_date_mod.o shr_tInterp_mod.o \
                    shr_string_mod.o shr_mpi_mod.o shr_log_mod.o shr_ncread_mod.o \
                    dshr_tInterp.o shr_file_mod.o dshr_domain.o dshr_const.o      \
                    dshr_kind.o dshr_bundle.o bundle_expected.o shr_orb_mod.o     \
                    shr_date_mod.o

# Test modules

# share codes
shr_alarm_mod.o : shr_alarm_mod.F90 shr_cal_mod.o shr_date_mod.o shr_sys_mod.o shr_kind_mod.o shr_log_mod.o
shr_cal_mod.o : shr_cal_mod.F90 shr_const_mod.o shr_sys_mod.o shr_kind_mod.o shr_log_mod.o
shr_const_mod.o : shr_const_mod.F90 shr_kind_mod.o 
shr_date_mod.o : shr_date_mod.F90 shr_cal_mod.o shr_sys_mod.o shr_kind_mod.o \
                 shr_log_mod.o shr_const_mod.o
shr_file_mod.o : shr_file_mod.F90 shr_sys_mod.o shr_kind_mod.o  shr_log_mod.o
shr_inputinfo_mod.o : shr_inputinfo_mod.F90 shr_file_mod.o shr_string_mod.o shr_sys_mod.o shr_ncio_mod.o \
                      shr_kind_mod.o shr_log_mod.o
shr_kind_mod.o : shr_kind_mod.F90  
shr_map_mod.o : shr_map_mod.F90 shr_const_mod.o shr_sys_mod.o shr_kind_mod.o shr_log_mod.o
shr_mpi_mod.o : shr_mpi_mod.F90 shr_kind_mod.o shr_log_mod.o
shr_msg_mod.o : shr_msg_mod.F90 shr_file_mod.o shr_log_mod.o
shr_ncio_mod.o : shr_ncio_mod.F90 shr_file_mod.o shr_string_mod.o shr_sys_mod.o shr_kind_mod.o shr_log_mod.o
shr_ncread_mod.o : shr_string_mod.o shr_kind_mod.o shr_sys_mod.o shr_file_mod.o \
                   shr_cal_mod.o shr_log_mod.o
shr_stream_mod.o : shr_sys_mod.o shr_kind_mod.o shr_const_mod.o shr_string_mod.o \
                   shr_ncread_mod.o shr_mpi_mod.o shr_log_mod.o
shr_orb_mod.o : shr_orb_mod.F90 shr_const_mod.o shr_sys_mod.o shr_kind_mod.o shr_log_mod.o
shr_string_mod.o : shr_string_mod.F90 shr_cal_mod.o shr_sys_mod.o shr_kind_mod.o shr_timer_mod.o shr_log_mod.o
shr_sys_mod.o : shr_sys_mod.F90 shr_mpi_mod.o shr_kind_mod.o shr_log_mod.o
shr_tInterp_mod.o : shr_tInterp_mod.F90 shr_const_mod.o shr_cal_mod.o \
                    shr_sys_mod.o shr_kind_mod.o shr_log_mod.o
shr_timer_mod.o : shr_timer_mod.F90 shr_sys_mod.o shr_kind_mod.o shr_log_mod.o
shr_mpi_mod.o : shr_kind_mod.F90 shr_log_mod.o shr_log_mod.o
test_mod.o : test_mod.F90 shr_kind_mod.o shr_sys_mod.o
shr_isnan.o : shr_isnan.c shr_isnan.h
shr_infnan_mod.o : shr_infnan_mod.F90 shr_kind_mod.o shr_isnan.h

# dshr code

dshr_domain.o: dshr_domain.F90 shr_sys_mod.o shr_string_mod.o shr_ncread_mod.o shr_mpi_mod.o shr_log_mod.o \
                    shr_file_mod.o dshr_kind.o dshr_const.o
dshr_const.o:  dshr_const.F90 shr_const_mod.o shr_kind_mod.o
dshr_kind.o:   dshr_kind.F90 shr_kind_mod.o
dshr_bundle.o: dshr_bundle.F90 dshr_domain.o shr_sys_mod.o shr_string_mod.o shr_ncread_mod.o shr_mpi_mod.o shr_log_mod.o \
                    shr_file_mod.o dshr_kind.o dshr_const.o
dshr_tInterp.o: dshr_tInterp.F90 shr_sys_mod.o shr_cal_mod.o shr_kind_mod.o shr_date_mod.o shr_tInterp_mod.o \
                shr_string_mod.o shr_mpi_mod.o shr_log_mod.o dshr_domain.o dshr_const.o      \
                dshr_kind.o dshr_bundle.o shr_orb_mod.o
bundle_expected.o: bundle_expected.F90 dshr_bundle.o dshr_kind.o shr_date_mod.o shr_sys_mod.o shr_orb_mod.o
dshr_rest.o : dshr_rest.F90 shr_sys_mod.o shr_timer_mod.o shr_date_mod.o shr_string_mod.o shr_stream_mod.o shr_log_mod.o \
              shr_file_mod.o dshr_kind.o dshr_bundle.o dshr_iocdf.o
dshr_iocdf.o: dshr_iocdf.F90 dshr_bundle.o dshr_domain.o dshr_kind.o dshr_const.o shr_sys_mod.o shr_date_mod.o shr_string_mod.o \
              shr_ncread_mod.o shr_log_mod.o


# ESMF code...
ifeq ($(ESMF_BLD),$(null))

include ../../utils/esmf_wrf_timemgr/unittests/Depends
AR  := ar
CPP := cpp

libesmf.a : $(WRFESMF_OBJS)
	$(RM) -f libesmf.a
	$(AR) $(ARFLAGS) libesmf.a $(WRFESMF_OBJS)
	$(RANLIB) libesmf.a

endif
wrf_error_fatal.o : wrf_error_fatal.F90 shr_sys_mod.o
wrf_message.o : wrf_message.F90
