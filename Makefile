FC = gfortran
FFLAGS = -ffast-math -Wall -march=native -O3 -fopenmp -mno-avx  #-fbounds-check #compiler flags
LDFLAGS = -fopenmp #link flags

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

PROG = main  main_automatic  #program name

#required objects for main: 
OBJS =
OBJS += constants.o
OBJS += plotroutines.o
OBJS += main_routines.o
OBJS += io.o
OBJS += initialize.o
OBJS += markov.o
OBJS += main.o

#required objects for main_automatic: 
OBJSA =
OBJSA += constants.o
OBJSA += plotroutines.o
OBJSA += main_routines.o
OBJSA += io.o
OBJSA += initialize.o
OBJSA += markov.o
OBJSA += main_automatic.o

all: $(PROG)

main: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

main_automatic: $(OBJSA) 
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f95
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) *.mod
	$(RM) *.png *.txt *.plt *.dat
