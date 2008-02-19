all:
	gcc -Wall -O3 -march=pentium-m -c basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c
	gcc -lgmp -O3 -march=pentium-m -o tester basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o

clean:
	rm -f tester tijdelijk gmon.out
	rm -f basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o test_scalars.o

debug:
	gcc -g -DKIJKEN -Wall -pedantic -std=c99 -c basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c
	gcc -lgmp -g -Wall -o tester basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o

profiler:
	gcc -pg -DPROFILER -O2 -march=pentium-m -Wall -c basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c
	gcc -pg -Wall -lgmp -O2 -march=pentium-m -o tester basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o

test:
	gcc -g -DKIJKEN -Wall -c scalar.c pol.c helper.c test_scalars.c
	gcc -g -lgmp -Wall -o tester test_scalars.o pol.o helper.o scalar.o

make_list:
	gcc -DLIST_F -Wall -c make_list.c silent_compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c
	gcc -lgmp -Wall -o tester make_list.o silent_compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o

input_pol:
	gcc -DINPUT_F -Wall -O3 -march=pentium-m -c basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c
	gcc -lgmp -O3 -march=pentium-m -o tester basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o

output_pol:
	gcc -DINPUT_F -DOUTPUT_LIST -Wall -c basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c
	gcc -lgmp -Wall -o tester basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o
