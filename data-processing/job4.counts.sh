#!/bin/bash
#PBS -q large
#PBS -l nodes=1:ppn=12
#PBS -l mem=120gb
#PBS -l walltime=84:00:00
#PBS -N Message.htseq.counts.ear
#PBS -j oe



cd $PBS_O_WORKDIR
JOBID=`echo $PBS_JOBID|tr '[' 'J'|tr ']' '_'`
{
	pwd

	module load python3

	NUM=3

	case $NUM in
		2)
			choice="second"
			;;
		3)
			choice="third"
			;;
		*)
			echo "no choice made, exitting"
			exit
			;;
	esac

	OUTDIR=htseq-counts/$choice
	echo "choice=[$choice]"
	echo "OUTDIR=[$OUTDIR]"

	i=0
	for fni in `find ./output-aligned/$choice -name "*sorted.bam"`
	do
		let i=i+1
		echo "<---------------------------start of $i--------------------------------->"
		date


		prefix=`basename $fni | sed -r 's/_1.fastq.gz.grcm38.sorted.bam$//g'`
		folder=`dirname $fni`
		OUTFILE=$OUTDIR/$prefix.$choice.gtf.counts
		echo "working on $fni (prefix=[$prefix])"
		echo "OUTFILE=[$OUTFILE]"


		echo "begin of htseq count"

	#	 python -c "import numpy; print(numpy.version.version)"
	#continue
		if [ -f $OUTDIR/$prefix.$folder.gtf.counts ]; then
			echo "[$OUTDIR/$prefix.$folder.gtf.counts] exists"
			echo "skipping"
			continue
		fi
		#/usr/bin/python2.7 -m HTSeq.scripts.count \
		python3 -m HTSeq.scripts.count \
			-f bam \
			-r pos \
			-m union \
			-s no \
			$fni \
			/PROJECT/gencode/mouse-gencode/gencode.vM11.annotation.gtf \
			> $OUTFILE

		#echo $cmd

		#$cmd
		
		echo "end of htseq count"

		echo "<---------------------------end of $i--------------------------------->"
#	break
	done

	date
	echo finished
} 2>> Message_hts_"$JOBID".txt 1>> Message_hts_"$JOBID".txt

