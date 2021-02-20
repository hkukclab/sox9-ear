#!/bin/bash
#PBS -q medium 
#PBS -l nodes=1:ppn=12
#PBS -l mem=50gb
#PBS -l walltime=24:00:00
#PBS -N Message.bam-to-wig-to-bigwig
#PBS -j oe
#PBS -J 111-126

cd $PBS_O_WORKDIR
JOBID=`echo $PBS_JOBID|tr '[' 'J'|tr ']' '_'`
{
	pwd

	index=$PBS_ARRAY_INDEX
	echo index=[$index]
	let from=1+index*10-10
	let to=10*index
	echo from=[$from]
	echo to=[$to]

	module load samtools

#	for bamfile in `ls output-aligned/third/*.sorted.bam -r`
	cnt=0
	for bamfile in `find  output-aligned/third/ -name "*.sorted.bam"`
	do
		let cnt=cnt+1
		date
		echo "-------------------------$cnt-----------------------------"
	#	if [ $cnt -lt $from ] || [ $cnt -gt $to ]; then
		if [ $cnt -ne "$index" ]; then
			echo $cnt not in range, skipping
			continue
		fi

		batch=`echo $bamfile|sed 's,output-aligned/,,g'|sed -r 's,/.*$,,g'`
		prefix=${batch}.`basename $bamfile |sed -r 's/.grcm38.sorted.bam|_1.fastq.gz.grcm38.sorted.bam//g'`
		wigfile=bam-to-wig/$prefix.wig
		bigwigfile=bam-to-wig/$prefix.bw

		echo "working on $bamfile (prefix=[$prefix])"
		echo "batch to [$batch]"
		echo "output to [$wigfile]"
		echo "output to [$bigwigfile]"

		if [ -f $bigwigfile ]; then
			echo $bigwigfile "exists, skipping"
			continue
		fi

		if [ -f ${wigfile}.gz ]; then
			echo ${wigfile}.gz "exists, directly wig-to-bigwig"
			gunzip ${wigfile}.gz

			/SOFTWARE/bam-to-wig/wigToBigWig \
				$wigfile \
				/USER/grcm38/chro.size.with.chr.extended.100bp \
				$bigwigfile

			gzip -vf ${wigfile}
			continue
		fi

		echo "start making wig"
		date
		java -jar /SOFTWARE/jvarkit/dist/bam2wig.jar \
			-w 3 -s 3   $bamfile  > $wigfile

		echo "start making big-wig"
		date
		/SOFTWARE/bam-to-wig/wigToBigWig \
			$wigfile \
			/USER/hisat2/grcm38/chro.size.with.chr.extended.100bp \
			$bigwigfile

		gzip -vf $wigfile


	done

	date
	echo finished
} 2>> Message_std_"$JOBID".txt 1>> Message_std_"$JOBID".txt
