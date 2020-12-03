MODE = zero
FC = mpif90

ifeq ($(MODE),zero)
#### Using the GNU compiler
FOPTS = -g -fPIC -D_GNU_SOURCE -O0 -Wall -Wstrict-aliasing  -fstrict-aliasing -ffixed-line-length-132
FOPTS += -fbounds-check
#FOPTS += -fsanitize=address
endif

ifeq ($(MODE),one)
#### Using the Portland Group/nVidia (PGI) compiler
FOPTS = -Mextend
endif

FOPTS += -DNO_DEBUG_

all: test

test:
	$(FC) $(FOPTS) main.f90

#test:
#	$(FC) $(FOPTS) -c main.f90 -o main.o
#	$(FC) $(FOPTS) -o a.out main.o


clean:
	rm -f *.o *.so *.a a.out *.mod
