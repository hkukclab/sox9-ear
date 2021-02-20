#!/bin/bash
#PBS -q large 
#PBS -l nodes=1:ppn=12
#PBS -l mem=120g
#PBS -l walltime=84:00:00
#PBS -N Message.hisat2.grcm38.single.cell.ear
#PBS -j oe
#PBS -J 1-2

cd $PBS_O_WORKDIR
JOBID=`echo $PBS_JOBID|tr '[' 'J'|tr ']' '_'`
{
	pwd



	hisatRef="/USER/hisat2/grcm38/genome"
	hisatDir="/USER/hisat2/hisat2-2.0.4"

	index=$PBS_ARRAY_INDEX
	NUM=3

	case $NUM in
	2)
		choice="second"
		OUTDIR="/PROJECT/output-aligned/second"
		INDIR="/PROJECT/2nd-batch-2019/CheahK_RNASeqSO_SS-190118-01a/primary_seq"
		;;
	3)
		choice="third"
		OUTDIR="/PROJECT/output-aligned/third"
		INDIR="/PROJECT/third-batch-aug2019/CheahK_RNASeqSO_SS-190715-02a/primary_seq/"
		;;
	*)
		echo "no choice makde, exitting"
		exit
		;;
	esac

	case $index in
		1) PATT="P9" ;;
		2) PATT="P10" ;;
		*) echo "no choice makde, exitting"; exit; ;;
	esac

	echo "choice=["$choice"]"
	echo "INDIR=["$INDIR"]"
	echo "OUTDIR=["$OUTDIR"]"
	echo "PATT=["$PATT"]"

	#exit

	module load samtools
	######################################################################################
	cd $OUTDIR
	pwd


	i=0
	for fni in `ls $INDIR/*1.fastq.gz |grep -e "$PATT"`
	do
		let i=i+1
		echo "<---------------------------start of $i--------------------------------->"
		date

		pref=`basename $fni|sed 's/_1.fastq.gz//g'`

		fq1=$fni
		fq2=`echo $fni|sed 's/_1.fastq/_2.fastq/g'`


		echo $fq1 
		echo $fq2 
		echo "prefix=[$pref]"

		countres=`find . -name "$pref*bam" |wc -l`

		if [ ! $countres -eq 0 ]; then
			echo "$pref.*.bam exists, skipping"
			continue
		fi
		
		echo "start aligning"

		$hisatDir/hisat2  \
			-x $hisatRef \
			-1 $fq1 \
			-2 $fq2 \
			-q \
			-p 12 \
			--no-unal \
			--ignore-quals \
			--novel-splicesite-outfile $pref.nov.spsite.txt \
			--add-chrname \
			--met-file $pref.met \
			--un-conc-gz /PROJECT/unaligned/second/$pref."$choice".un-conc.gz \
			-S $pref.sam

	echo "finished aligning"	

	samtools view -bS $pref.sam -o $pref.grcm38.bam
	rm  $pref.sam

	samtools sort -o $pref.grcm38.sorted.bam $pref.grcm38.bam
        samtools index $pref.grcm38.sorted.bam
	
	rm $pref.grcm38.bam

	echo "<---------------------------end of $i--------------------------------->"
	done


	date
	echo "finished"

} 2>> Message_align_"$JOBID".txt 1>> Message_align_"$JOBID".txt

