# Create output folder
mkdir -p bam

# This will make the script stop on errors.
set -eu

IDX=annot/22

for SAMPLE in UHR HBR;
do
	for REPLICATE in 1 2 3;
	do
		LABEL=${SAMPLE}_${REPLICATE}
		hisat2 ${IDX} -1 data/${LABEL}_R1.fq -2 data/${LABEL}_R2.fq > bam/${LABEL}.sam
		samtools view -b bam/${LABEL}.sam | samtools sort > bam/${LABEL}.bam
		samtools index bam/${LABEL}.bam
	done
done