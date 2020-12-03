MODE = zero
FC = mpif90

ifeq ($(MODE),zero)
#### Using the GNU compiler
FOPTS = -g -fPIC -D_GNU_SOURCE -O0 -Wall -Wstrict-aliasing  -fstrict-aliasing -ffixed-line-length-132
FOPTS += -fbounds-check
FOPTS += -cpp
#FOPTS += -fsanitize=address
endif

ifeq ($(MODE),one)
#### Using the Portland Group/nVidia (PGI) compiler
FOPTS = -Mextend
endif

FOPTS += -DNO_DEBUG_

all: test

test:
	$(FC) $(FOPTS) -c llist_mod.f90 -o llist_mod.o
	$(FC) $(FOPTS) main.f90 llist_mod.o

#test:
#	$(FC) $(FOPTS) -c main.f90 -o main.o
#	$(FC) $(FOPTS) -o a.out main.o


clean:
	rm -f *.o *.so *.a a.out *.mod
