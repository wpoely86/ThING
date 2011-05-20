AM_CPPFLAGS = -DDATADIR='"$(datadir)"'

lib_LTLIBRARIES = libThING.la

libThING_la_SOURCES = Gauss.cpp  Gauss.h  HF.cpp  HF.h  input.cpp  input.h  lapack.h  MxElem.cpp  MxElemFiller.cpp  MxElemFiller.h  MxElem.h  R.cpp  R.h  Diag.cpp  Diag.h

LIBTHING_LIBRARY_VERSION=1:0:0

libThING_la_LDFLAGS = -version-info $(LIBTHING_LIBRARY_VERSION)

include_HEADERS = ThING.h

ThINGdir = $(includedir)/ThING

ThING_HEADERS = Gauss.h  HF.h  input.h  lapack.h  MxElemFiller.h  MxElem.h  R.h  Diag.h

pkgconfigdir = $(libdir)/pkgconfig

pkgconfig_DATA = libThING.pc
