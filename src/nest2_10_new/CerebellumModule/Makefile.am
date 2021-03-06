# Automake file for external dynamic modules for NEST
#
# Hans Ekkehard Plesser, April 2008
# Automake file for the Developer Module
#
# libmymodule is built as a normal, installable library.
# It will be installed to $prefix/lib by make install.
#
# Headers from this directory are not to be installed upon
# make install. They are therefore included in _SOURCES.


# 1. Exchange "my" in "mymodule" with the name of your model below
#    (ten locations).
# 2. Add all .cpp and .h files from your code as *_SOURCES. Header files
#    are given only so that they will be included in the tarball if you
#    run "make dist" on your module.
# 3. The mymodule* stuff creates a module that can be loaded at runtime.
#    It is called mymodule.so.
# 4. The libmymodule* stuff creates a library against which NEST can be
#    linked.

libdir= @libdir@/nest

# We need to install the module header for static linking on BlueGene
include_HEADERS = cerebellummodule.h

# All other source files
source_files=  cerebellummodule.cpp \
               ExponentialTable.cpp ExponentialTable.h \
               TrigonometricTable.cpp TrigonometricTable.h \
               stdp_sin_connection.h \
               stdp_cos_connection.h \
               archiving_node_cs.cpp archiving_node_cs.h \
               archiving_node_cos.cpp archiving_node_cos.h \
               histentry_cs.cpp histentry_cs.h \
               histentry_cos.cpp histentry_cos.h \
               iaf_cond_exp_cs.cpp iaf_cond_exp_cs.h \
               iaf_cond_exp_cos.cpp iaf_cond_exp_cos.h \
               cd_poisson_generator.cpp cd_poisson_generator.h \
               rbf_poisson_generator.cpp rbf_poisson_generator.h
               
if BUILD_DYNAMIC_USER_MODULES
  lib_LTLIBRARIES= libcerebellummodule.la cerebellummodule.la
  cerebellummodule_la_CXXFLAGS= @AM_CXXFLAGS@
  cerebellummodule_la_SOURCES=  $(source_files)
  cerebellummodule_la_LDFLAGS=  -module
else
  lib_LTLIBRARIES= libcerebellummodule.la
endif

libcerebellummodule_la_CXXFLAGS= @AM_CXXFLAGS@ -DLINKED_MODULE
libcerebellummodule_la_SOURCES=  $(source_files)

MAKEFLAGS= @MAKE_FLAGS@

AM_CPPFLAGS= @NEST_CPPFLAGS@ \
             @INCLTDL@

.PHONY: install-slidoc

pkgdatadir=@datadir@/nest

nobase_pkgdata_DATA=\
	sli/cerebellummodule-init.sli

install-slidoc:
	NESTRCFILENAME=/dev/null $(DESTDIR)$(NEST_PREFIX)/bin/sli --userargs="@HELPDIRS@" $(NEST_PREFIX)/share/nest/sli/install-help.sli

install-data-hook: install-exec install-slidoc

EXTRA_DIST= sli
