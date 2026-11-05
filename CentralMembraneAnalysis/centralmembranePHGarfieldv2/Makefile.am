##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libcentralmembranephgarfield.la

libcentralmembranephgarfield_la_SOURCES = \
  TpcCentralMembranePHGarfield.cc


libcentralmembranephgarfield_la_LIBADD = \
  -lFROG \
  -lffaobjects \
  -lphool \
  -lodbc++ \
  -lSubsysReco \
  -lg4detectors_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -ltrack_reco \
  -ltpc_io \
  -lPHGarfield \
  -lcentralmembraneanalysisv2

pkginclude_HEADERS = \
  TpcCentralMembranePHGarfield.h

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

# Just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_centralmembranephgarfield

testexternals_centralmembranephgarfield_SOURCES = \
  testexternals.cc

testexternals_centralmembranephgarfield_LDADD = \
  libcentralmembranephgarfield.la \
  -lPHGarfield \
  -lcentralmembraneanalysisv2

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f $(BUILT_SOURCES)