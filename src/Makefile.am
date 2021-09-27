# src/Makefile.am for fsim

include $(top_srcdir)/common.mk

noinst_LTLIBRARIES = libplane.la
libplane_la_SOURCES = gsl_if.hpp gsl_if.cpp plane.cpp plane.hpp
libplane_la_CPPFLAGS = ${AM_CPPFLAGS}
libplane_la_CXXFLAGS = -I. -I${top_srcdir}
libplane_la_CXXFLAGS += ${AM_CXXFLAGS}
libplane_la_CXXFLAGS += ${AM_CRAZY_MATH}

bin_PROGRAMS = fsim
fsim_SOURCES = fsim.cpp
fsim_CPPFLAGS = ${AM_CPPFLAGS}
fsim_CXXFLAGS = -I. -I${top_srcdir}
fsim_CXXFLAGS += ${AM_CXXFLAGS} ${AM_CRAZY_MATH}
fsim_LDADD = libplane.la
fsim_LDFLAGS = $(gsl_LIBS) -lboost_program_options -lboost_system
