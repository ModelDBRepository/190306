#!/bin/bash

numneuron=120
ranseed=-35
#gexcstart=0.0805 # normal
#gexcstart=0.081 #c1
#gexcstart=0.080 #cm1
#gexcstart=0.0815 #c2
gexcstart=0.0795 #cm2
ranseedm=$(echo $ranseed | sed s/-/m/g)
a=1
ad=$(echo $a | sed s/[.]/d/g)
d=0.81
dd=$(echo $d | sed s/[.]/d/g)

#k=0
i=0
y=0
#dirprefix="/scratch/ernestho/abc_homo_${ad}_${dd}/"
#dirprefix="/scratch/ernestho/abc_homo_${ad}_${dd}_c1/"
#dirprefix="/scratch/ernestho/abc_homo_${ad}_${dd}_cm1_nomp/"
#dirprefix="/scratch/ernestho/abc_homo_${ad}_${dd}_c2_nomp/"
dirprefix="/scratch/ernestho/abc_homo_${ad}_${dd}_cm2_nomp/"


fileprefix="simdata"
us="_"
filesuffix="abc_vmddiscont_homo.dat"

#for k in `seq 0.000 0.002 0.006`
#  do
#for taun in `seq 1.5 1.5 1.5`
#  do
#for gin in `seq -f %.3f 0.010 0.001 0.0551`
#for gin in `seq -f %.3f 0.002 0.001 0.0351`
#for gin in `seq -f %.3f 0.036 0.001 0.0701`
for gin in `seq -f %.3f 0.002 0.001 0.0701`

  do
  #for noise in `seq 0.00027 0.00025 0.001771`
  #for noise in `seq 0.00202 0.00025 0.002521`
 #for noise in `seq 0.00277 0.00025 0.002771`
for noise in `seq 0.00027 0.00025 0.003021`
    do
    gind=$(echo $gin | sed s/[.]/d/g)
    gexc=$gexcstart
      #gexc=$(echo "scale=4; $gexcstart+$k" | bc)
      #let gexc = $gexc+$k
    gexcd=$(echo $gexc | sed s/[.]/d/g)
    taund=$(echo $taun | sed s/[.]/d/g)
    noised=$(echo $noise | sed s/[.]/d/g)
      #echo  $dirprefix$fileprefix$us$numneuron$us$ranseedm$us$gind$us$noised$us$gexcd$us$taund$us$filesuffix
    filename=$dirprefix$fileprefix$us$numneuron$us$ranseedm$us$gind$us$ad$us$dd$us$gexcd$us$noised$us$filesuffix
    if [ -f ${filename} ]; then
	continue
    fi
    if [ $y -ge 840 ]; then
	continue
    fi

    if [ $(( $y % 8)) -eq 0 ]; then
	cat > script${i}.pbs <<EOF
#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=4:30:00
#PBS -N job_simulation_abc_noise_${i}
cd SimuCPP_parallel
 ./simulation_abc_jumpsynapse_homo $numneuron $ranseed $gin $a $d $gexc $noise > $filename &

EOF
    elif [ $(( $y % 8)) -lt 7 ]; then
	cat >> script${i}.pbs <<EOF
 ./simulation_abc_jumpsynapse_homo $numneuron $ranseed $gin $a $d $gexc $noise > $filename &
EOF
    else
	cat >> script${i}.pbs <<EOF
 ./simulation_abc_jumpsynapse_homo $numneuron $ranseed $gin $a $d $gexc $noise > $filename &
wait
EOF
	
	qsub script${i}.pbs
	i=$(($i+1))
    fi
    y=$(($y+1))
  done
done
#done

if [ $(($y % 8)) -ne 0 ]; then
    cat >> script${i}.pbs <<EOF
wait
EOF
    qsub script${i}.pbs
#rm script${i}.pbs
fi
