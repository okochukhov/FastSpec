# code version
VERSION = 3.2

# HDF-enabled gfortran compiler
FF = h5fc

# OS, hostname, architecture
OSNAME := $(shell uname)
HOSTNAME := $(shell hostname -s)
MACHINE := $(shell uname -m)

# standard compiler specific flags.
FLAG = -O3 -w -fopenmp -flto=auto
LIBS = 

# EOS is compiled with standard gfortran
FF_EOS = gfortran
FLAG_EOS = -O3 -std=legacy

# Self-contained ("make static") link options. These are macOS/Homebrew specific.
MINOS      = -mmacosx-version-min=11.0
FLAG_STATIC = -O3 -w -flto=auto
HDF5DIR    = $(shell brew --prefix hdf5)/lib
AECDIR     = $(shell brew --prefix libaec)/lib
GCCDIR     = $(shell dirname $(shell $(FF_EOS) -print-file-name=libgfortran.a))
STATICLIBS = -static-libgfortran -static-libgcc -static-libquadmath \
             $(HDF5DIR)/libhdf5_hl_fortran.a $(HDF5DIR)/libhdf5_hl_f90cstub.a \
             $(HDF5DIR)/libhdf5_fortran.a $(HDF5DIR)/libhdf5_f90cstub.a \
             $(HDF5DIR)/libhdf5_hl.a $(HDF5DIR)/libhdf5.a \
             $(AECDIR)/libsz.a $(AECDIR)/libaec.a \
             $(GCCDIR)/libgomp.a -lz

# Location of EOS sources and include files.
EOSDIR_DEV = $(HOME)/Codes/Eos/100726_SME_fixed
EOSDIR = $(if $(wildcard $(EOSDIR_DEV)/eos.f),$(EOSDIR_DEV),$(SOURCEDIR))

# source file directory
SOURCEDIR = .

# build directory
BDIR = .

# set of routines to compile for normal version of the code
SRC = shared utils eos eos_eqns eos_math ioutils contop voigt hydtab_data hydtab hydropac hestark solver opacity convol synth fastspec

# relevant object file file dependencies
OBJ = $(addsuffix .o,$(SRC))

.PHONY: all static clean arch

# default build is for normal version
all:
	@echo "------------------------------------------------"
	@echo "  Building FastSpec for" $(OSNAME) on $(HOSTNAME)
	@echo "------------------------------------------------"
	make $(BDIR)/fastspec

$(BDIR)/fastspec: $(OBJ)
	$(FF) $(FLAG) -o $@ $^ $(LIBS)

# self-contained build: same objects, statically linked runtime/HDF5/szip
static: $(OBJ)
	@echo "------------------------------------------------"
	@echo "  Building static FastSpec for" $(OSNAME) on $(HOSTNAME) "($(MACHINE))"
	@echo "------------------------------------------------"
	$(FF_EOS) $(FLAG_STATIC) $(MINOS) -o $(BDIR)/fastspec_$(OSNAME)_$(MACHINE) $(OBJ) $(STATICLIBS)

%.o: $(SOURCEDIR)/%.f90 $(SOURCEDIR)/Makefile
	$(FF) -c $(FLAG) $< -o $@

%.o: $(SOURCEDIR)/%.f $(SOURCEDIR)/Makefile
	$(FF) -c $(FLAG) $< -o $@

fastspec.o: $(SOURCEDIR)/fastspec.f90 $(SOURCEDIR)/Makefile
	echo "CHARACTER(*), PARAMETER :: Version = '$(VERSION)'" > VERSION
	echo "CHARACTER(*), PARAMETER :: VerDate = '`date +%Y-%m-%d`'" >> VERSION
	$(FF) -c $(FLAG) $< -o $@
	rm VERSION

eos.o: $(EOSDIR)/eos.f $(EOSDIR)/SIZES.EOS $(EOSDIR)/DEFAULT.EOS Makefile
	$(FF_EOS) -c $(FLAG_EOS) $< -o $@

eos_eqns.o: $(EOSDIR)/eos_eqns.f $(EOSDIR)/SIZES.EOS Makefile
	$(FF_EOS) -c $(FLAG_EOS) $< -o $@

eos_math.o: $(EOSDIR)/eos_math.f Makefile
	$(FF_EOS) -c $(FLAG_EOS) $< -o $@

utils.o:    shared.o
ioutils.o:  shared.o utils.o
contop.o:   shared.o
hydropac.o: shared.o voigt.o hydtab.o
hydtab.o:   hydtab_data.o
hestark.o:  shared.o utils.o
opacity.o:  shared.o utils.o voigt.o hestark.o
solver.o:   shared.o utils.o
convol.o:   shared.o
synth.o:    shared.o utils.o ioutils.o contop.o voigt.o hydropac.o hestark.o solver.o opacity.o convol.o
fastspec.o: shared.o utils.o ioutils.o opacity.o convol.o synth.o

clean:
	rm -f *.o *.mod fastspec

arch:
	COPYFILE_DISABLE=1 tar --no-xattrs -cvzhf fastspec.tgz \
        Makefile shared.f90 utils.f90 ioutils.f90 opacity.f90 solver.f90 \
        hydtab_data.f90 hydtab.f90 hydropac.f90 hestark.f90 contop.f90 \
        voigt.f90 convol.f90 synth.f90 fastspec.f90 \
	-C "$(EOSDIR)" eos.f SIZES.EOS DEFAULT.EOS eos_eqns.f eos_math.f
