AM_CXXFLAGS =  -Wall -Wextra -Wpedantic -Werror -Wpedantic -Wshadow
AM_CXXFLAGS += -Wfatal-errors
AM_CXXFLAGS += -fanalyzer # Only for gcc 10
AM_CPPFLAGS = -DDATADIR='"$(pkgdatadir)"'
AM_CRAZY_MATH = -O3 -march=native -mtune=native
AM_CRAZY_MATH += -funroll-loops
AM_CRAZY_MATH += @X86_FEATURE_CFLAGS@
