#!/bin/bash
#PBS -q large 
#PBS -l nodes=1:ppn=10
#PBS -l mem=120gb
#PBS -l walltime=84:00:00
#PBS -N MSG.cufflinks-G-gtf
#PBS -j oe

cd $PBS_O_WORKDIR
JOBID=`echo $PBS_JOBID|tr '[' 'J'|tr ']' '_'`
{
	pwd


	NUM=5


	case $NUM in
	2)
		choice="second"
		;;
	3)
		choice="third"
		;;
	*)
		echo "wrong choice, exittng"
		exit
		;;
	esac

	echo choice=[$choice]

	INDIR="output-aligned/$choice/"
	module load samtools

	i=0
	for fni in `find $INDIR -name "*[CD|WT].bam"`
	do
		let i=i+1
		echo "<---------------------------start of $i--------------------------------->"
		date

		outbam=`echo $fni|sed 's/bam/sorted.bam/g'`
		echo outbam=[$outbam]

		samtools sort -o $outbam $fni
		samtools index $outbam


		prefix=`basename $fni|sed -r 's/\.bam$//g'`
		echo "working on $fni (prefix=[$prefix])"


		echo "[now cufflink]"
		cufflinks --library-type fr-firststrand \
			--max-bundle-frags 100000000 \
			-G /mouse-gencode/gencode.vM11.annotation.gtf \
			-o "cufflinks-G-gtf/$choice/$prefix-grcm38-cufflinks" \
			-p 10 \
			$outbam

		echo "[finish cufflinks]"
		gzip -vf cufflinks-G-gtf/$choice/$prefix-grcm38-cufflinks/*
		echo "<---------------------------end of $i--------------------------------->"
	done



	date
	echo "finished"

} 2>> Message_cfl_"$JOBID".txt 1>> Message_cfl_"$JOBID".txt
