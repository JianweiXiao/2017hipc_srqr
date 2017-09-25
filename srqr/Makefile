include make.inc

TARGETS = test_srqr

all: 	$(TARGETS)

# executable files
test_srqr: test_srqr.o partial_dgeqpf.o randperm.o srqr.o
	$(F90) $^ -o $@ $(LIBS)

# objective fortran files
srqr.o: srqr.f90
	$(F90) -c $(F90FLAGS) $^

test_srqr.o: test_srqr.f90
	$(F90) -c $(F90FLAGS) $^

partial_dgeqpf.o: partial_dgeqpf.f
	$(F90) -c $(F90FLAGS) $^

randperm.o: randperm.f90
	$(F90) -c $(F90FLAGS) $^

clean:
	rm -f *.o
	rm -f TARGETS
	rm -f *.out

