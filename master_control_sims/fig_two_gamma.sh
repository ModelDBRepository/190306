#!/bin/bash

# THIS IS THE HEAD SCRIPT THAT RUNS THE COMPILED C++ CODE FOR NETWORK SIMULATION.  IT SUPPLIES MOST OF THE PARAMETERS (ADJUSTABLE) TO THE SIMULATION, CREATES THE SUITABLE SLURM SCRIPTS FOR JOB SUBMISSION (AND SUBMITS IT) AND PREPARES THE OUTPUT DIRECTORY. (may, 2016)

#numcorespersim=16 # number of cores per simulation
#numprocinnode=16 # number of cores in a node
numcorespersim=8 # number of cores per simulation
numprocinnode=8 # number of cores in a node

numprocinnode_m1=$(echo $numprocinnode - $numcorespersim | bc)

exedirname="/users/XXXXX/jcneurosci_code/network_c++_code/"  # FULL PATH OF THE DIRECTORY THAT HAS THE COMPILED C++ CODE
mscriptdirname="/users/XXXXX/jcneurosci_code/master_control_sims/" # FULL PATH OF THE DIRECTORY THAT HAS THE HEAD BATCH SCRIPT (THAT IS, THIS SCRIPT)
connfiledirname="/gpfs/scratch/XXXXX/9x9_test_8/" # FULL PATH OF THE DIRECTORY THAT CONTAINS THE CONNECTIVITY TEMPLATE FILES
slurmscriptname="script_fig_two_gamma_" # BUILDING THE SLURM SCRIPT NAME HERE
i=0
y=0
cd $mscriptdirname
rm *.slurm # REMOVES ALL PREVIOUS SLURM SCRIPTS

# EXTERNAL DC CURRENT PARAMETERS

i_ext=0.0
i_extd=$(echo $i_ext | sed s/[.]/d/g)

i_extfs=0.0
i_extfsd=$(echo $i_extfs | sed s/[.]/d/g) #FS = inhib

i_extib=0.0
i_extibd=$(echo $i_extib | sed s/[.]/d/g)

# GLIAL COMPONENT PARAMETERS

k_for=0.0008
k_back=0.0008
rise=-0.15 # path -0.83, -0.25
Ko_eq_pump=3.6 # from 5
Ko_eq_glia=7.5 # from 7 # canonical 15
#max_pump_current=2 # path 5, 13
max_pump_current=1.45 # from 2 # path 5, 13

k_for_in=0.0008 # INs that go in DB
k_back_in=0.0008
rise_in=-0.15 # canonical -1.15
Ko_eq_pump_in=3 # from 3
Ko_eq_glia_in=7.5 # from 7
#max_pump_current_in=2 # path 10
max_pump_current_in=1.9 # from 2

k_for_in_2=0.0008 # spiking interneurons
k_back_in_2=0.0008
rise_in_2=-0.15 # canonical -1.15
Ko_eq_pump_in_2=3
Ko_eq_glia_in_2=7.5
#max_pump_current_in_2=2 # path 10
max_pump_current_in_2=1.9


# STIMULATION DC CURRENT PARAMETERS 
stim_start=40000 # ms
stim_end=42500 # ms 
stim_strength=2.5 # uA/cm2


stim_start_in=40000
stim_end_in=42500
stim_strength_in=2.5 #


# SYNAPTIC PARAMETERS (ee_syn = g(e->e) etc)
ee_syn=0.0007  # default 0.0005
ei_syn=0.0007
ie_syn=0.025 # 0.004
ii_syn=0.025 #0.004


# BACKGROUND SYNAPTIC ENVIRONMENT
gi_zero=0.084 # MEAN
si_zero=0.020 # SIGMA
inhib_ramp_factor=1.0 # RAMPING UP FACTOR OF EFFECTIVE SYNAPTIC INHIBITION (S^{TERM}_{MAX} IN MANUSCRIPT) DURING DC STIMULATION (SET TO 1 WHEN THERE IS NO GABA DEPLETION) 
# LOOP OVER A RANGE OF GE_ZERO VALUES
for ge_zero in `seq -f %.5f 0.01008 0.00002 0.01010`

			do 
# PREPARING THE OUTPUT DIRECTORY
		condpattern="test_new_prng_${i_extd}_${i_extfsd}"
		mkdir -p "/users/XXXXX/scratch/${condpattern}"
# SIGMA FOR BACKGROUND EXCITATION
for se_zero in `seq -f %.5f 0.00250 0.005 0.00250`


										do
												ee_synd=$(echo $ee_syn | sed s/[.]/d/g)
												ei_synd=$(echo $ei_syn | sed s/[.]/d/g)
												ie_synd=$(echo $ie_syn | sed s/[.]/d/g)
												ii_synd=$(echo $ii_syn | sed s/[.]/d/g)
												ge_zerod=$(echo $ge_zero | sed s/[.]/d/g)
												gi_zerod=$(echo $gi_zero | sed s/[.]/d/g)

												k_ford=$(echo $k_for | sed s/[.]/d/g)
												k_backd=$(echo $k_back | sed s/[.]/d/g)
												rised=$(echo $rise | sed s/[.]/d/g)
												Ko_eq_pumpd=$(echo $Ko_eq_pump | sed s/[.]/d/g)
												Ko_eq_gliad=$(echo $Ko_eq_glia | sed s/[.]/g/g)
												max_pump_currentd=$(echo $max_pump_current | sed s/[.]/g/g)
	
	
												k_for_ind=$(echo $k_for_in | sed s/[.]/d/g)
												k_back_ind=$(echo $k_back_in | sed s/[.]/d/g)
												rise_ind=$(echo $rise_in | sed s/[.]/d/g)
												Ko_eq_pump_ind=$(echo $Ko_eq_pump_in | sed s/[.]/d/g)
												Ko_eq_glia_ind=$(echo $Ko_eq_glia_in | sed s/[.]/g/g)
												max_pump_current_ind=$(echo $max_pump_current_in | sed s/[.]/g/g)


												k_for_in_2d=$(echo $k_for_in_2 | sed s/[.]/d/g)
												k_back_in_2d=$(echo $k_back_in_2 | sed s/[.]/d/g)
												rise_in_2d=$(echo $rise_in_2 | sed s/[.]/d/g)
												Ko_eq_pump_in_2d=$(echo $Ko_eq_pump_in_2 | sed s/[.]/d/g)
												Ko_eq_glia_in_2d=$(echo $Ko_eq_glia_in_2 | sed s/[.]/g/g)
												max_pump_current_in_2d=$(echo $max_pump_current_in_2 | sed s/[.]/g/g)

												stim_startd=$(echo $stim_start | sed s/[.]/d/g)
												stim_endd=$(echo $stim_end | sed s/[.]/d/g)
												stim_strengthd=$(echo $stim_strength | sed s/[.]/d/g)

												stim_start_ind=$(echo $stim_start_in | sed s/[.]/d/g)
												stim_end_ind=$(echo $stim_end_in | sed s/[.]/d/g)
												stim_strength_ind=$(echo $stim_strength_in | sed s/[.]/d/g)
												inhib_ramp_factord=$(echo $inhib_ramp_factor | sed s/[.]/d/g)



												#ggliad=$(echo $gglia | sed s/[.]/d/g)
												#kzeroinfd=$(echo $kzeroinf | sed s/[.]/d/g)
												#gglia_ind=$(echo $gglia_in | sed s/[.]/d/g)
												#kzeroinf_ind=$(echo $kzeroinf_in | sed s/[.]/d/g)

												se_zerod=$(echo $se_zero | sed s/[.]/d/g)
												si_zerod=$(echo $si_zero | sed s/[.]/d/g)

												patternprefix="multidim_${ge_zerod}_${se_zerod}_${gi_zerod}_${si_zerod}_${k_ford}_${k_backd}_${rised}_${Ko_eq_pumpd}_${Ko_eq_gliad}_${max_pump_currentd}_${stim_startd}_${stim_endd}_${stim_strengthd}_${k_for_ind}_${k_back_ind}_${rise_ind}_${Ko_eq_pump_ind}_${Ko_eq_glia_ind}_${max_pump_current_ind}_${k_for_in_2d}_${k_back_in_2d}_${rise_in_2d}_${Ko_eq_pump_in_2d}_${Ko_eq_glia_in_2d}_${max_pump_current_in_2d}_${stim_start_ind}_${stim_end_ind}_${stim_strength_ind}_${inhib_ramp_factord}_syn_${ee_synd}_${ei_synd}_${ie_synd}_${ii_synd}"  # gglia and kzeroinf have to be changed manually
# Dec 18th, 2015 (data node fixed; move back to data)
												dirname="/users/XXXXX/scratch/${condpattern}/${patternprefix}/"
											#	dirname="/users/XXXXX/data/XXXXX/${condpattern}/${patternprefix}/"
# End
												echo $dirname
												mkdir -p $dirname
												cd $mscriptdirname
												cp fig_two_gamma_run_template.sh  $dirname
												# cp mpi_test.sh $dirname
												cd $exedirname
												cd $dirname
												sed -e s,FILEPREFIX,${dirname},g -e s,CONNPREFIX,${connfiledirname},g -e s,PATTERN,${patternprefix},g -e s,MSCRIPTPREFIX,${mscriptdirname},g -e s,IEXTRS,${i_ext},g -e s,IEXTIB,${i_extib},g -e s,IEXTFS,${i_extfs},g -e s,GGLIARISEIN2,${rise_in_2},g -e s,GGLIARISEIN,${rise_in},g -e s,GGLIARISE,${rise},g -e s,GGLIAFORIN2,${k_for_in_2},g -e s,GGLIAFORIN,${k_for_in},g -e s,GGLIAFOR,${k_for},g -e s,GGLIABACKIN2,${k_back_in_2},g -e s,GGLIABACKIN,${k_back_in},g -e s,GGLIABACK,${k_back},g -e s,GGLIAKEQPUMPIN2,${Ko_eq_pump_in_2},g -e s,GGLIAKEQPUMPIN,${Ko_eq_pump_in},g -e s,GGLIAKEQPUMP,${Ko_eq_pump},g -e s,GGLIAKEQGLIAIN2,${Ko_eq_glia_in_2},g -e s,GGLIAKEQGLIAIN,${Ko_eq_glia_in},g -e s,GGLIAKEQGLIA,${Ko_eq_glia},g -e s,GGLIAMAXPUMPIN2,${max_pump_current_in_2},g -e s,GGLIAMAXPUMPIN,${max_pump_current_in},g -e s,GGLIAMAXPUMP,${max_pump_current},g -e s,STIMSTARTIN,${stim_start_in},g -e s,STIMENDIN,${stim_end_in},g -e s,STIMSTRENGTHIN,${stim_strength_in},g -e s,STIMSTART,${stim_start},g -e s,STIMEND,${stim_end},g -e s,STIMSTRENGTH,${stim_strength},g -e s,EE,${ee_syn},g -e s,EI,${ei_syn},g -e s,IE,${ie_syn},g -e s,II,${ii_syn},g -e s,Gee,${ge_zero},g -e s,See,${se_zero},g -e s,Gie,${gi_zero},g -e s,Sie,${si_zero},g -e s,Gei,${ge_zero},g -e s,Sei,${se_zero},g -e s,Gii,${gi_zero},g -e s,Sii,${si_zero},g -e s,POSS_RATE,0,g -e s,POSS_INCR,0,g -e s,INHIBRAMPFACTOR,${inhib_ramp_factor},g fig_two_gamma_run_template.sh > fig_two_gamma_run.sh
												chmod u+x fig_two_gamma_run.sh
												cd $mscriptdirname
# WRITES THE REQUIRED SLURM SCRIPTS
												if  [ $(( $y % $numprocinnode)) -eq 0 ]; then
													cat > ${slurmscriptname}${i}.slurm <<EOF
#!/bin/bash

# Runtime request:
#SBATCH --qos=bibs-truccolo-condo
# SBATCH --time=00:30:00
#SBATCH --time=30:00:00
# SBATCH --time=35:00:00
# Takes approximately 12 hours for 16x16 simulation of 10s using 1 core													
# #########################
# Queue request
# SBATCH -p default-batch
# SBATCH --constraint="e5-2600" 
# Sandy bridge (sucks--so slow)
# SBATCH --constraint="e5000"
# Regular
# ###############################################
# Use 8 nodes (4 MPI per node); altogether 32 MPI jobs, and there are 2 open mp threads per MPI job (so 8 tasks per node).
# SBATCH --nodes=1
# SBATCH --tasks-per-node=8
#SBATCH --nodes=1-1
#SBATCH --tasks-per-node=$numprocinnode
# SBATCH --exclusive													
# ###############################################
# Specify a job name:
#SBATCH -J fig_two_gamma

# Specify an output file
#SBATCH -o ./fig_two_gamma-%j.out
#SBATCH -e ./fig_two_gamma-%j.out

cd  $mscriptdirname
${dirname}fig_two_gamma_run.sh  &
EOF
#echo $numprocinnode
#echo $numprocinnode_m1 
											else
												cat >>${slurmscriptname}${i}.slurm<<EOF
${dirname}fig_two_gamma_run.sh  & 
EOF
											fi

											if [ $(( $y % $numprocinnode)) -eq $numprocinnode_m1 ]; then
												i=$(($i+1))
											fi
									y=$(($y+$numcorespersim))
											done
										done
									
								
						
										
for scriptfile in $(ls *.slurm)
	do
	cat >>$scriptfile<<EOF
wait
EOF
# FINALLY SUMBITS THE SLURM SCRIPTS
sbatch ${scriptfile}
	done
