PROG =	dmat

SRCS =	distributed_k_matrices.f90 dummy_main.f90 k_point_matrix.f90 \
	matrix.f90 matrix_mapping.f90 numbers.f90 proc_mapping.f90

OBJS =	distributed_k_matrices.o dummy_main.o k_point_matrix.o matrix.o \
	matrix_mapping.o numbers.o proc_mapping.o

LIBS =-lblacs-openmpi -lblacsF77init-openmpi	

CC = cc
CFLAGS = -O
FC = mpif90
FFLAGS = -std=f2003 -fcheck=all -O -Wall -Wextra -fbacktrace -g
F90 = mpif90
F90FLAGS = -std=f2003 -fcheck=all -O -Wall -Wextra -fbacktrace -g
LDFLAGS = -fbacktrace -g 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

dummy_main.o: matrix_mapping.o proc_mapping.o
matrix.o: matrix_mapping.o numbers.o
matrix_mapping.o: numbers.o proc_mapping.o
proc_mapping.o: 
