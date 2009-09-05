#  Copyright 1999-2009 ImageMagick Studio LLC, a non-profit organization
#  dedicated to making software imaging solutions freely available.
#
#  You may not use this file except in compliance with the License.  You may
#  obtain a copy of the License at
#
#    http://www.imagemagick.org/script/license.php
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# 
#  Copyright (C) 2003 - 2008 GraphicsMagick Group
#
#  Makefile for building PerlMagick.

# If source files missing, see if they can be obtained via VPATH
perl-sources:
	@if test -n "$(VPATH)" ; then \
	  echo "Linking PerlMagick Sources ..." ; \
	  imagemagick=`(cd $(VPATH) ; pwd)` && \
	  ( cd $(PERLMAGICK) && \
	    sh $$imagemagick/lndir.sh $$imagemagick/$(PERLMAGICK) ) \
	fi ; \
	touch perl-sources

if WITH_PERL

PERLMAGICK=PerlMagick
PERLMAKEMAKER=$(PERLMAGICK)/Makefile.PL
PERLMAKEFILE=$(PERLMAGICK)/Makefile

PERLMAGICK_ALL_LOCAL_TARGETS = all-perl
PERLMAGICK_INSTALL_EXEC_LOCAL_TARGETS = install-exec-perl
PERLMAGICK_INSTALL_DATA_LOCAL_TARGETS = 
PERLMAGICK_UNINSTALL_LOCAL_TARGETS = uninstall-exec-perl
PERLMAGICK_CLEAN_LOCAL_TARGETS = clean-perl
PERLMAGICK_DISTCLEAN_LOCAL_TARGETS = clean-perl
PERLMAGICK_MAINTAINER_CLEAN_LOCAL_TARGETS = distclean-local
PERLMAGICK_CHECK_LOCAL_TARGETS = check-perl

if WITH_PERL_DYNAMIC

$(PERLMAKEFILE): perl-sources $(PERLMAKEMAKER)
	cd $(PERLMAGICK) && @PERL@ Makefile.PL $(PERL_MAKE_OPTIONS)

install-exec-perl: $(MAGICKCORE_LIBS) $(MAGICKWAND_LIBS) $(PERLMAKEFILE)
	cd $(PERLMAGICK) && $(MAKE) CC='@CC@' install

# add MAGICKCORE_LIBS in deps to fix parallel make issue
all-perl: $(MAGICKCORE_LIBS) $(MAGICKWAND_LIBS) $(PERLMAKEFILE)
	cd $(PERLMAGICK) && $(MAKE) CC='@CC@'

uninstall-exec-perl: $(PERLMAKEFILE)
	echo "Uninstall not supported for PerlMagick"

check-perl: $(MAGICKCORE_LIBS) $(MAGICKWAND_LIBS) $(PERLMAKEFILE)
	cd $(PERLMAGICK) && $(MAKE) CC='@CC@' test

else
if WITH_PERL_STATIC

PERLSTATICNAME=PerlMagick

$(PERLMAKEFILE): perl-sources $(PERLMAKEMAKER)
	cd $(PERLMAGICK) && @PERL@ Makefile.PL MAP_TARGET=$(PERLSTATICNAME) $(PERL_MAKE_OPTIONS)

$(PERLMAGICK)/$(PERLSTATICNAME): $(PERLMAKEFILE) magick/libMagickCore.la
	( rm -f $(PERLMAGICK)/$(PERLSTATICNAME) ; cd $(PERLMAGICK) && $(MAKE) CC='@CC@' $(PERLSTATICNAME) )

all-perl: $(PERLMAGICK)/$(PERLSTATICNAME)

install-exec-perl: $(PERLMAGICK)/$(PERLSTATICNAME)
	rm -f "$(DESTDIR)$(BIN_DIR)/$(PERLSTATICNAME)"
	if test "x$(DESTDIR)" = "x" -o "$(PERL_SUPPORTS_DESTDIR)" = 'yes' ; then \
	  ( cd $(PERLMAGICK) && \
	    $(MAKE) -f Makefile.aperl CC='@CC@' inst_perl MAP_TARGET=$(PERLSTATICNAME) \
	    INSTALLBIN="$(BIN_DIR)" \
	  ) ; \
	else \
	  ( cd $(PERLMAGICK) && \
	    $(MAKE) -f Makefile.aperl CC='@CC@' inst_perl MAP_TARGET=$(PERLSTATICNAME) \
	    INSTALLBIN="$(DESTDIR)$(BIN_DIR)" PREFIX="$(DESTDIR)$(prefix)" \
	  ) ; \
	fi

uninstall-exec-perl:
	rm -f '$(DESTDIR)$(BIN_DIR)/$(PERLSTATICNAME)'

check-perl: $(PERLMAGICK)/$(PERLSTATICNAME)
	cd $(PERLMAGICK) && $(MAKE) -f Makefile.aperl CC='@CC@' test

endif # WITH_PERL_STATIC
endif # WTIH_PERL_DYNAMIC


clean-perl: $(PERLMAKEFILE)
	(cd $(PERLMAGICK) && \
	$(MAKE) CC='@CC@' clean && \
	rm -f Makefile.old PerlMagick && \
	rm -f t/output* t/jng/*_tmp.jng t/*/output* )
	rm -f perl-sources

distclean-perl: clean-perl

else
# Satisfy makefile requirements if not building PERL
all-perl:
install-exec-perl:
uninstall-exec-perl:
check-perl:
clean-perl:
distclean-perl:
endif # WITH_PERL
