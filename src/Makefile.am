F77 = $(FC)
FCFLAGS += @ARCH@ @FORTRANFLAGS@
AM_FFLAGS = $(FCFLAGS)

#-- Here is the simple test file

check_PROGRAMS = test 
TESTS = test

test_SOURCES= test.f
test_LDADD = libqtpie.la

#-- Here is a standalone calculator

bin_PROGRAMS = onexyz 

onexyz_SOURCES = onexyz.f
onexyz_LDADD = libqtpie.la

lib_LTLIBRARIES = libqtpie.la
libqtpie_la_SOURCES = atomicunits.f sparse.f parameters.f io.f cg.f geometry.f qtpie.f properties.f

# If Tinker exists, compile tinker API function

if WITH_TINKER
libqtpie_la_SOURCES += api_tinker.f
endif

## Here, compile Slater orbitals only if requested
## But really, need to patch source to do this for real
## For now, just compile everything

#if USE_SLATER
libqtpie_la_SOURCES += factorial.f sto-int.f
#else
libqtpie_la_SOURCES += gto-int.f
#endif

CLEANFILES = *.mod

if HAVE_DOXYGEN
doc:
	doxygen
endif

if WITH_TINKER
tinker: libqtpie.a
	echo Linking QTPIE with TINKER
	make -C @TINKER_DIR@ 
	make -C @TINKER_DIR@ all 
	make -C @TINKER_DIR@ rename
	make -C @TINKER_DIR@ create_links
endif

