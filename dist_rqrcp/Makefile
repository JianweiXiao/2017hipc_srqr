include make.inc

TARGETS	= xtimer_qr

all: 	$(TARGETS)

# executable files
xtimer_qr: timer_qr.o check_qrf.o partial_pdgeqrf.o partial_pdgeqpf.o  partial_pdgeqp3.o pdgeqp3s.o pdqppiv.o check_qpf.o rqrcp.o partial_qr_swap.o
	$(F90) $^ -o $@ $(LIBS)

# objective files
partial_pdgeqpf.o: partial_pdgeqpf.f
	$(F90) -c $(F90FLAGS) $^

partial_pdgeqp3.o: partial_pdgeqp3.f
	$(F90) -c $(F90FLAGS) $^

partial_pdgeqrf.o: partial_pdgeqrf.f
	$(F90) -c $(F90FLAGS) $^

partial_qr_swap.o: partial_qr_swap.f90
	$(F90) -c $(F90FLAGS) $^

pdgeqp3s.o: pdgeqp3s.f
	$(F90) -c $(F90FLAGS) $^

timer_qr.o: timer_qr.c
	$(CC) -c $(F90FLAGS) $^

rqrcp.o: rqrcp.f90
	$(F90) -c $(F90FLAGS) $^

check_qrf.o: check_qrf.f90
	$(F90) -c $(F90FLAGS) $^

check_qpf.o: check_qpf.f90
	$(F90) -c $(F90FLAGS) $^

pdqppiv.o: pdqppiv.f
	$(F90) -c $(F90FLAGS) $^

clean:
	rm -f *.o
	rm -f $(TARGETS)
	rm -f *.out

