#CC=/usr/local/gcc-4.4.1/bin/g++
# gcc compiler
CC=mpic++
LL=mpic++
#INCLUDE=-I/gpfs/runtime/opt/openmpi/1.6.2/include
INCLUDE=-I/gpfs/runtime/opt/openmpi/1.8/include
#openmpi/1.8 centos-libs/6.5
#INCLUDE=-I/gpfs/runtime/opt/mvapich2/1.9a/include/
#CPPFLAGS=-O3 -ffast-math -std=c++0x -pg $(INCLUDE) -fopenmp -I/usr/inc
CPPFLAGS=-O3 -D__FIGTWO__ -D__HETERO__ -ffast-math -std=c++0x -pg $(INCLUDE) -fopenmp -I/usr/inc # use this when there is no gaba depletion
#CPPFLAGS=-O3 -D__FIGTWO__ -D__HETERO__ -D__DEPLETE__ -ffast-math -std=c++0x -pg $(INCLUDE) -fopenmp -I/usr/inc  #use this for gaba depletion
#CPPFLAGS=-O3 -D__FIGTWO__ -D__DEPLETE__ -ffast-math -std=c++0x -pg $(INCLUDE) -fopenmp -I/usr/inc

#CPPFLAGS=-O3 -D__FIGTWO__ -ffast-math -std=c++0x -pg $(INCLUDE) -fopenmp -I/usr/inc

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


OBJS5=mapFunctions.o currentEquations.o synapseEquations.o integration_routines.o srng.o gasdev.o Neuron.o cellVariables.o cellParameters.o Neuron_read_write.o cortical_assignment.o LFP.o LFP_tools.o Diffusion.o Diffusion_tools.o test_neurons_openmp_dev.o
OBJS8=srng.o cortical_assignment.o post_result_process.o
OBJS9=srng.o cortical_assignment.o prep_spike_timing.o
#OBJS12=srng.o cortical_assignment.o sort_syn.o
OBJS10=srng.o cortical_assignment.o cortical_assignment_main_custom.o



all: test_neurons_openmp_dev post_result_process prep_spike_timing cortical_assignment_main_custom 


test_neurons_openmp_dev: $(OBJS5)
	$(CC) $(CPPLFLAGS) $(OBJS5) -pg -o test_neurons_openmp_dev
post_result_process: $(OBJS8)
	$(CC) $(CPPLFLAGS) $(OBJS8) -pg -o post_result_process
prep_spike_timing: $(OBJS9)
	$(CC) $(CPPLFLAGS) $(OBJS9) -pg -o prep_spike_timing
cortical_assignment_main_custom: $(OBJS10)
	$(CC) $(CPPLFLAGS) $(OBJS10) -pg -o cortical_assignment_main_custom
#sort_syn: $(OBJS12)
#	$(CC) $(CPPLFLAGS) $(OBJS12) -pg -o sort_syn
#cortical_assignment_main_flex: $(OBJS14)
#	$(CC) $(CPPLFLAGS) $(OBJS14) -pg -o cortical_assignment_main_flex


currentEquations.o: currentEquations.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF) $<  
synapseEquations.o: synapseEquations.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF) $<
integration_routines.o: integration_routines.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $< 
srng.o: srng.c
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $< 
gasdev.o: gasdev.c
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF) $<  
Neuron.o: Neuron.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $< 
cellVariables.o:  cellVariables.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $< 
cellParameters.o: cellParameters.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $< 
Neuron_read_write.o: Neuron_read_write.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $<
mapFunctions.o: mapFunctions.cc
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $<	
cortical_assignment.o: cortical_assignment.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $<	
LFP.o: LFP.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $<	
LFP_tools.o: LFP_tools.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $<
Diffusion.o: Diffusion.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $<	
Diffusion_tools.o: Diffusion_tools.cpp
	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $<		
#sort_syn.o: sort_syn.cpp
#	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $<	

#cortical_assignment_may21.o: cortical_assignment_may21.cpp
#	$(CC) -c $(CPPFLAGS) -o $@  $(CPPSUFF)  $<	


clean:
	rm -f *.o

