#!/bin/bash

cd ~/jcneurosci_code/assign_connectivity_c++_code/ # cd to the executables
numproc=8; # number of MPI processes used in simulation
fileprefix="/gpfs/scratch/XXXXX/9x9_test_${numproc}/" # Location of the directory that stores the connectivity matrix files

connectionfile="anderson_connection.DAT";
synfile="syn_sorted.DAT";
sizefile="syn_size.DAT";

infofile="conn_info.dat";

mkdir -p ${fileprefix}

./cortical_assignment_main_flex $numproc ${fileprefix}${connectionfile} ${fileprefix}${synfile} ${fileprefix}${sizefile} >  ${fileprefix}${infofile}

perl rename_syn_files.pl ${fileprefix}

