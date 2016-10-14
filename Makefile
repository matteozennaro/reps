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
## Do not touch things below this line #########################################
################################################################################

DEBUGFLAG=

BOUNDARY=BOUNDARY_CONDITIONS_MODULE
MODULES=MODULES
SOURCE=SOURCE

BOUNDARYFILES=$(BOUNDARY)/BC.o
MODULEFILES=$(MODULES)/background_functions.o $(MODULES)/neutrino_distribution_function.o $(MODULES)/RungeKutta_solver_3fluids.o $(MODULES)/boltzmann_solver.o $(MODULES)/read_ini_file.o $(MODULES)/general_purpose.o $(MODULES)/write_output.o
SOURCEFILE=$(SOURCE)/reps.o

all: boundary modules source bc reps

reps: $(SOURCEFILE) $(MODULEFILES)
	$(CC) -o reps $(SOURCEFILE) $(MODULEFILES) $(CFLAGS) $(OMPFLAGS) $(OFLAGS) $(DEBUGFLAG) $(ADD_HEADERS)

bc: $(BOUNDARYFILES) $(MODULEFILES)
	$(CC) -o BC $(BOUNDARYFILES) $(MODULEFILES) $(CFLAGS) $(OMPFLAGS) $(OFLAGS) $(DEBUGFLAG) $(ADD_HEADERS)

source: $(SOURCEFILE)

$(SOURCEFILE): %.o : %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(OMPFLAGS) $(OFLAGS) $(DEBUGFLAG) $(ADD_HEADERS)

modules: $(MODULEFILES)

$(MODULEFILES): %.o : %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(OMPFLAGS) $(OFLAGS) $(DEBUGFLAG) $(ADD_HEADERS)

boundary: $(BOUNDARYFILES)

$(BOUNDARYFILES): %.o : %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(OMPFLAGS) $(OFLAGS) $(DEBUGFLAG) $(ADD_HEADERS)

debug: DEBUGFLAG=-DDEBUG
debug: all

clean:
	rm -rf $(SOURCE)/*.o
	rm -rf $(MODULES)/*.o
	rm -rf $(BOUNDARY)/*.o
	rm reps
	rm BC
