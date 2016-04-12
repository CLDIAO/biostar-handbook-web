Set up a few shortcuts to facilitate data exploration:
 	
 	# Convert SAM to BAM.
	alias sam2bam='samtools view -b - | samtools sort'

	# This is the reference genome.
	REF=annot/22.fa

	# This will be the name of the index
	IDX=annot/22

	# These are the coordinates of the genes.
	GTF=annot/22.gtf

	# The sequencing data that contain the first in pair
	R1=data/UHR_1_R1.fq

	# The sequencing data that contain the second in pair
	R2=data/UHR_1_R2.fq

	# We also need to build an index out of the genome
	hisat2-build $REF $IDX