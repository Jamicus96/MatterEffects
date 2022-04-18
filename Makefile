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

BIN = Matter_Oscillations.exe
OBJ = Matter_Oscillations.o
TXT = results.txt

all: $(BIN)

Matter_Oscillations.exe: Matter_Oscillations.o
	g++ -std=c++11 Matter_Oscillations.o -o Matter_Oscillations.exe $(LDFLAGS) $(local_LDFLAGS)

%.o : %.cpp
	g++ -std=c++11 -c $< $(INCFLAGS)

clean:
	rm -f $(BIN) $(OBJ) $(TXT)

