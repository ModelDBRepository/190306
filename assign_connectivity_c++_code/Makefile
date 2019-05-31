#CC=/usr/local/gcc-4.4.1/bin/g++
#module load openmpi/1.8

# gcc compiler
CC=mpic++
LL=mpic++
INCLUDE=-I/gpfs/runtime/opt/openmpi/1.8/include
CPPFLAGS=-O3 -ffast-math -std=c++0x -pg $(INCLUDE) -fopenmp -I/usr/inc
#CPPSUFF=-pg
CPPSUFF=$(INCLUDE)
#CPPLFLAGS=-fopenmp
CPPLFLAGS=-fopenmp
CPPFTTWFLAGS=-lfftw3 -lm


# Intel compiler
#CC=icpc
#LL=icpc
#CPPFLAGS=-O3 -openmp
#CPPLFLAGS=-openmp -openmp-lib=compat

OBJS14=srng.o cortical_assignment.o cortical_assignment_main_flex.o



all:	cortical_assignment_main_flex


cortical_assignment_main_flex: $(OBJS14)
	$(CC) $(CPPLFLAGS) $(OBJS14) -pg -o cortical_assignment_main_flex

srng.o: srng.c
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $< 
cortical_assignment.o: cortical_assignment.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $<	



clean:
	rm -f *.o

