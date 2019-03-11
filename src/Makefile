PROG =	dmat

SRCS =	ks_matrix_module.f90 dummy_main.f90 ks_array_module.f90 \
	matrix.f90 matrix_mapping.f90 numbers.f90 proc_mapping.f90 \
	scalapack_interfaces.f90

OBJS =	ks_matrix_module.o dummy_main.o ks_array_module.o matrix.o \
	matrix_mapping.o numbers.o proc_mapping.o scalapack_interfaces.o

LIBS =-lblacs-openmpi -lblacsF77init-openmpi -lscalapack-openmpi -llapack -lblas

TESTPROG = test

CC = cc
CFLAGS = -O
FC = mpif90
FFLAGS = -std=f2003 -fcheck=all -O -Wall -Wextra -fbacktrace -g
F90 = mpif90
F90FLAGS = -std=f2003 -fcheck=all -O -Wall -Wextra -fbacktrace -g
LDFLAGS = -fbacktrace -g

.PHONY: clean remake all objects
.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

all: $(PROG)

objects: $(OBJS)

remake: clean all

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod



ks_matrix_module.o: matrix.o numbers.o
dummy_main.o: ks_matrix_module.o matrix.o matrix_mapping.o \
	proc_mapping.o ks_array_module.o
ks_array_module.o: numbers.o ks_matrix_module.o
matrix.o: matrix_mapping.o numbers.o scalapack_interfaces.o
matrix_mapping.o: numbers.o proc_mapping.o
proc_mapping.o:
scalapack_interfaces.o: numbers.o
