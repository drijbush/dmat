TEST = scf_step_real

#SRCS =	distributed_k_matrices.f90 k_point_matrix.f90 \
#	matrix.f90 matrix_mapping.f90 numbers.f90 proc_mapping.f90 \
#	scalapack_interfaces.f90

OBJS =	../src/distributed_k_matrices.o ../src/k_point_matrix.o ../src/matrix.o \
	../src/matrix_mapping.o ../src/numbers.o ../src/proc_mapping.o ../src/scalapack_interfaces.o

LIBS =-lblacs-openmpi -lblacsF77init-openmpi -lscalapack-openmpi -llapack -lblas

SRCDIR = ../src/

CC = cc
CFLAGS = -O
FC = mpif90
FFLAGS = -std=f2003 -fcheck=all -O -Wall -Wextra -fbacktrace -g
F90 = mpif90
F90FLAGS = -std=f2003 -fcheck=all -O -Wall -Wextra -fbacktrace -g
LDFLAGS = -fbacktrace -g

.PHONY: clean remake all
.SUFFIXES: $(SUFFIXES) .f90

all: $(TEST)

remake: clean all

$(OBJS): 
	$(MAKE) -C $(SRCDIR) objects


$(TEST): $(TEST).f90 $(OBJS)
	$(F90) $(LDFLAGS) -o $@ -I$(SRCDIR) $(OBJS) $(TEST).f90 $(LIBS)

clean:
	$(MAKE) -C $(SRCDIR) clean
	rm -f $(TEST) *.o *.mod

