## Your favourite compiler
CC=gcc
#CC=icc

OMPFLAGS=-fopenmp
#OMPFLAGS=-openmp

## Other options
CFLAGS=-lm -Wall
OFLAGS=-Ofast
ADD_HEADERS=-I./INCLUDE

################################################################################

DEBUGFLAG=

MODULES=MODULES
SOURCE=SOURCE

BOUNDARYFILES=$(MODULES)/boundary_conditions.o
MODULEFILES=$(MODULES)/background_functions.o $(MODULES)/neutrino_distribution_function.o $(MODULES)/RungeKutta_solver_3fluids.o $(MODULES)/boltzmann_solver.o $(MODULES)/read_ini_file.o $(MODULES)/general_purpose.o $(MODULES)/write_output.o $(MODULES)/boundary_conditions.o
SOURCEFILE=$(SOURCE)/reps.o

all: modules source reps

reps: $(SOURCEFILE) $(MODULEFILES)
	$(CC) -o reps $(SOURCEFILE) $(MODULEFILES) $(CFLAGS) $(OMPFLAGS) $(OFLAGS) $(DEBUGFLAG) $(ADD_HEADERS)

source: $(SOURCEFILE)

$(SOURCEFILE): %.o : %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(OMPFLAGS) $(OFLAGS) $(DEBUGFLAG) $(ADD_HEADERS)

modules: $(MODULEFILES)

$(MODULEFILES): %.o : %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(OMPFLAGS) $(OFLAGS) $(DEBUGFLAG) $(ADD_HEADERS)

debug: DEBUGFLAG=-DDEBUG
debug: all

clean:
	rm -rf $(SOURCE)/*.o
	rm -rf $(MODULES)/*.o
	rm -f reps
