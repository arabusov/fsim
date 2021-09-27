SUBDIRS = src share
ACLOCAL_AMFLAGS = -I m4
dist_doc_DATA = README

clean-local:
	@rm -f */*.i */*.ii */*.s
	@rm -f Makefile Makefile.in
	@rm -f config.status config.log
	@rm -rf autom4te.cache/
	@rm -f aclocal.m4
	@rm -f configure~
	@rm -rf .deps
	@rm -rf */.deps */Makefile */Makefile.in
	@rm -f compile missing install-sh depcomp
	@rm -f install.sh~
	@rm -f configure
	@rm -f build-aux/compile
	@rm -rf config.guess config.sub
	@rm -rf build-aux/config.guess~ build-aux/config.sub~
	@rm -rf build-aux/config.guess build-aux/config.sub build-aux/depcomp
	@rm -rf build-aux/install-sh build-aux/missing build-aux/ltmain.sh
	@rm -rf libtool m4
