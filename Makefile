all:
	gcc -Wall -O3 -march=native -c main.c  basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c  char_p_0.c
	gcc -lgmp -O3 -march=native -o tester main.o  basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o  char_p_0.o

clean:
	rm -f tester tijdelijk gmon.out
	rm -f main.o  char_p_0.o  basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o  test_scalars.o  make_list.o

debug:
	gcc -g -DKIJKEN -Wall -pedantic -std=c99 -c main.c  basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c  char_p_0.c
	gcc -lgmp -g -Wall -o tester main.o  basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o  char_p_0.o

profiler:
	gcc -pg -DPROFILER -O1 -march=native -Wall -c  main.c  basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c  char_p_0.c
	gcc -pg -Wall -lgmp -O1 -march=native -o tester  main.o  basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o  char_p_0.o

test:
	gcc -O3 -Wall -c scalar.c pol.c helper.c test_scalars.c
	gcc -lgmp -O3 -Wall -o tester test_scalars.o pol.o helper.o scalar.o

input_pol:
	gcc -DINPUT_F -Wall -O3 -march=native -c main.c  basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c  char_p_0.c
	gcc -lgmp -O3 -march=native -o tester main.o  basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o  char_p_0.o

output_pol:
	gcc -DOUTPUT_LIST -Wall -c main.c  basis.c  compute.c  delta.c  grobner.c  helper.c  pol.c  reduce.c  scalar.c  char_p_0.c
	gcc -lgmp -Wall -o tester main.o  basis.o  compute.o  delta.o  grobner.o  helper.o  pol.o  reduce.o  scalar.o  char_p_0.o
