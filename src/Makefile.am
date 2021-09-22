include $(top_srcdir)/common.mk

lib_LTLIBRARIES = plane.la
plane_la_SOURCES = plane.cpp plane.hpp
plane_la_CPPFLAGS=-fPIC
plane_la_CPPFLAGS+=${AM_CPPFLAGS}
plane_la_CXXFLAGS=-fPIC -I. -I${top_srcdir}
plane_la_CXXFLAGS+=${AM_CXXFLAGS}
plane_la_CXXFLAGS+=${AM_CRAZY_MATH}
plane_la_LDFLAGS= -shared -module -avoid-version
