Set up a few shortcuts to facilitate data exploration:
 	
 	# Alias to convert and sort SAM to BAM.
	alias sam2bam='samtools view -b - | samtools sort'

	# This is the reference genome.
	REF=annot/22.fa

	# This will be the name of the index.
	IDX=annot/22

	# These are the coordinates of the annotations.
	GTF=annot/22.gtf

	# The sequencing data that contains the first in pair.
	R1=data/UHR_1_R1.fq

	# The sequencing data that contains the second in pair.
	R2=data/UHR_1_R2.fq
	
	# Create the folder that will store the BAM files.
	mkdir -p bam