### prefix for running on HPC (installed in folder):
prefix = /mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/antinu/GLoBES
### prefix for running on personal Linux computer (installed with root privilege):
# prefix = /usr/local
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)

BIN = Matter_Oscillations.exe PMNS_m2.exe Mat_Os_Consts.exe
OBJ = Matter_Oscillations.o PMNS_m2.o Mat_Os_Consts.o
TXT = results.txt

all: $(BIN)

Matter_Oscillations.exe: Matter_Oscillations.o
	g++ -std=c++11 Matter_Oscillations.o -o Matter_Oscillations.exe $(LDFLAGS) $(local_LDFLAGS)

PMNS_m2.exe: PMNS_m2.o cplx.hpp
	g++ -std=c++11 PMNS_m2.o cplx.hpp -o PMNS_m2.exe

Mat_Os_Consts.exe: Mat_Os_Consts.o
	g++ -std=c++11 Mat_Os_Consts.o -o Mat_Os_Consts.exe

%.o : %.cpp
	g++ -c $< $(INCFLAGS)

clean:
	rm -f $(BIN) $(OBJ) $(TXT)

