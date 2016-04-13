# Create output folder
mkdir -p bam

# This will make the script stop on errors.
set -eu

IDX=annot/22
GTF=annot/22.gtf

for SAMPLE in UHR HBR;
do
	for REPLICATE in 1 2 3;
	do
		LABEL=${SAMPLE}_${REPLICATE}
		tophat2 -G $GTF -o tophat_${LABEL} $IDX data/${LABEL}_R1.fq  data/${LABEL}_R2.fq
		cp tophat_${LABEL}/accepted_hits.bam bam/${LABEL}.bam
		samtools index bam/${LABEL}.bam
	done
done