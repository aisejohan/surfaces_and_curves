all:
	gcc -Wall -O3 -march=nocona -c basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c
	gcc -lgmp -O3 -march=nocona -o tester basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o

clean:
	rm -f tester tijdelijk gmon.out
	rm -f basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o test_scalars.o

debug:
	gcc -g -DKIJKEN -Wall -c basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c
	gcc -lgmp -g -Wall -o tester basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o

profiler:
	gcc -g -pg -DPROFILER -Wall -c basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c
	gcc -g -pg -Wall -lgmp -o tester basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o

test:
	gcc -O3 -Wall -c scalar.c pol.c helper.c test_scalars.c
	gcc -lgmp -O3 -Wall -o tester test_scalars.o pol.o helper.o scalar.o

maaklijst:
	gcc -g -DINPUT_F -DOUTPUT_LIST -Wall -c basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c
	gcc -lgmp -g -Wall -o tester basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o
