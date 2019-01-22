PROG =	dmat

SRCS =	dummy_main.f90 k_point_matrix.f90 matrix.f90 numbers.f90

OBJS =	dummy_main.o k_point_matrix.o matrix.o numbers.o

LIBS =	

CC = cc
CFLAGS = -O
FC = gfortran
FFLAGS = -std=f2003 -fcheck=all -O -Wall -Wextra
F90 = gfortran
F90FLAGS = -std=f2003 -fcheck=all -O -Wall -Wextra
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

dummy_main.o: k_point_matrix.o
k_point_matrix.o: matrix.o numbers.o
matrix.o: numbers.o
