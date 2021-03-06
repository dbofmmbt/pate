#-----------------------------------------------------------------------
# File    : seqwog.mak
# Contents: build seqwog program (on Windows systems)
# Author  : Christian Borgelt
# History : 2010.08.06 file created from relim makefile
#           2010.08.22 module escape added (for module tabread)
#           2016.04.20 creation of dependency files added
#-----------------------------------------------------------------------
THISDIR  = ..\..\seqwog\src
UTILDIR  = ..\..\util\src
TRACTDIR = ..\..\tract\src

CC       = cl.exe
DEFS     = /D WIN32 /D NDEBUG /D _CONSOLE /D _CRT_SECURE_NO_WARNINGS
CFLAGS   = /nologo /W3 /O2 /GS- $(DEFS) /c $(ADDFLAGS)
INCS     = /I $(UTILDIR) /I $(TRACTDIR)

LD       = link.exe
LDFLAGS  = /nologo /subsystem:console /incremental:no
LIBS     = 

HDRS     = $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h   \
           $(UTILDIR)\symtab.h   $(UTILDIR)\error.h    \
           $(UTILDIR)\tabread.h  $(UTILDIR)\tabwrite.h \
           $(TRACTDIR)\tract.h   $(TRACTDIR)\patspec.h \
           $(TRACTDIR)\report.h
OBJS     = $(UTILDIR)\arrays.obj   $(UTILDIR)\idmap.obj    \
           $(UTILDIR)\escape.obj   $(UTILDIR)\tabread.obj  \
           $(UTILDIR)\tabwrite.obj $(UTILDIR)\scform.obj   \
           $(TRACTDIR)\taread.obj  $(TRACTDIR)\patspec.obj \
           $(TRACTDIR)\report.obj  seqwog.obj
PRGS     = seqwog.exe

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:          $(PRGS)

seqwog.exe:   $(OBJS) seqwog.mak
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) /out:$@

#-----------------------------------------------------------------------
# Main Programs
#-----------------------------------------------------------------------
seqwog.obj:   $(HDRS) seqwog.mak
	$(CC) $(CFLAGS) $(INCS) /D SWOG_MAIN seqwog.c /Fo$@

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(UTILDIR)\arrays.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak arrays.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\memsys.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak memsys.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\idmap.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak idmap.obj    ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\escape.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak escape.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\tabread.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak tabread.obj  ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\tabwrite.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak tabwrite.obj ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\scform.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak scform.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(TRACTDIR)\taread.obj:
	cd $(TRACTDIR)
	$(MAKE) /f tract.mak taread.obj  ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(TRACTDIR)\patspec.obj:
	cd $(TRACTDIR)
	$(MAKE) /f tract.mak patspec.obj ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(TRACTDIR)\report.obj:
	cd $(TRACTDIR)
	$(MAKE) /f tract.mak report.obj  ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)

#-----------------------------------------------------------------------
# Install
#-----------------------------------------------------------------------
install:
	-@copy $(PRGS) ..\..\..\bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
localclean:
	-@erase /Q *~ *.obj *.idb *.pch $(PRGS)

clean:
	$(MAKE) /f seqwog.mak localclean
	cd $(TRACTDIR)
	$(MAKE) /f tract.mak localclean
	cd $(UTILDIR)
	$(MAKE) /f util.mak clean
	cd $(THISDIR)
