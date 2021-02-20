#!/bin/bash
#PBS -q large 
#PBS -l nodes=1:ppn=12
#PBS -l mem=120gb
#PBS -l walltime=84:00:00
#PBS -N MSG.job1.scenic
#PBS -j oe

cd $PBS_O_WORKDIR
JOBID=`echo $PBS_JOBID|tr '[' '_'|tr ']' '_'`
{
	pwd

	#module load python2/2.7.14
	module load R

	#source activate py27
	Rscript --version

	NUM=1
	#$PBS_ARRAYID
	export PYTHONPATH=0
	date
	####################################################
	case $NUM in
		1)
			Rscript step1.run.prepare.input.R
			Rscript step2.run.Corr.and.GENIE3.R 
			Rscript step3.run.scenic.R 
			;;
		2)
			;;
		3)
			;;
		*)
			echo no job
			exit
			;;
	esac




	date
	####################################################
	source deactivate
	echo finished
} 2>> Message_scnc_"$JOBID".txt 1>> Message_scnc_"$JOBID".txt
