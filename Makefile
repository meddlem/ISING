FC = gfortran
FFLAGS = -ffast-math -Wall -march=native -O3 -fopenmp -mno-avx #-fcheck=bounds -Warray-temporaries #compiler flags
LDFLAGS = -fopenmp #link flags

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

PROG = main #program name

#required objects: 
OBJS =
OBJS += constants.o
OBJS += plotroutines.o
OBJS += output_processing.o
OBJS += io.o
OBJS += initialize.o
OBJS += Wolff.o
OBJS += Swendsen_Wang.o
OBJS += markov.o
OBJS += main.o

all: $(PROG)

main: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f95
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) *.mod
	$(RM) *.png *.txt *.plt *.dat
